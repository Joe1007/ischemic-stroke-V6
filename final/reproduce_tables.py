#!/usr/bin/env python3
"""V6 adaptation of manuscript comparisons; never label lesion T as original F.
All fits and selection use training patients only. No selection by outer scores.
LightGBM pred_contrib provides exact tree-path-dependent SHAP in raw log odds.
"""
import argparse, json, sys, re, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
INPUTS=HERE.parent/'inputs'
sys.path.insert(0,str(HERE/'deps'))
sys.path.insert(0,str(HERE.parent))
import model_v6 as v6
mb=v6.mb
FLAGS={'affected','firstorder_available','contrast_available','texture_available'}

def estimator(kind,seed,jobs):
    if kind=='LightGBM': return mb.lightgbm_model(seed,jobs)
    if kind=='RF': return RandomForestClassifier(n_estimators=100,max_depth=4,random_state=seed,n_jobs=jobs)
    if kind=='XGBoost':
        from xgboost import XGBClassifier
        return XGBClassifier(n_estimators=100,learning_rate=.015,random_state=seed,n_jobs=jobs,eval_metric='logloss')
    from catboost import CatBoostClassifier
    return CatBoostClassifier(loss_function='Logloss',eval_metric='AUC',learning_rate=.01,iterations=30,
        od_type='Iter',depth=9,random_seed=seed,thread_count=jobs,verbose=False,allow_writing_files=False)

def select(a,y,names,kind,selector,keep,seed,jobs):
    active=np.arange(a.shape[1]);trace=[]
    if selector=='backward': stages=[(keep,len(active))]
    else: stages=[(10000,40000),(950,2000),(50,150),(keep,1)]
    for target,step in stages:
        target=max(target,keep)
        while len(active)>target:
            m=estimator(kind,seed,jobs);m.fit(a[:,active],y)
            order=np.argsort(-np.asarray(m.feature_importances_),kind='stable')
            active=active[order[:max(target,len(active)-step)]];trace.append(len(active))
    return active,trace

def fit(a,b,y,names,kind,seed,jobs,explain):
    tx,ty=mb.smote_training_only(a,y,seed)
    m=estimator(kind,seed,jobs);m.fit(tx,ty)
    p=m.predict_proba(b)[:,1]
    sv=None
    if explain:
        sv=m.predict(b,pred_contrib=True)
        raw=m.predict(b,raw_score=True)
        np.testing.assert_allclose(sv.sum(axis=1),raw,atol=1e-8,rtol=1e-8)
    return p,sv

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=4)
    ap.add_argument('--tasks',nargs='+',default=list(mb.TASKS));ap.add_argument('--fold-limit',type=int)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--worker-task');ap.add_argument('--split-index',type=int)
    ap.add_argument('--output',type=Path,default=HERE/'tables_20260917');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    oldpath=INPUTS/'v5_features.npz';newpath=HERE.parent/'output/features.npz'
    old=np.load(oldpath);new=np.load(newpath);np.testing.assert_array_equal(old['ids'],new['ids'])
    import sklearn,lightgbm,xgboost,catboost
    cfg={'adaptation':'C=conservative clinical; R=V5 global lesion Original; S=AAL3 burden; T=V6 lesion-intersection plus flags; Fwhole=V5 whole-AAL3 Original unclip_p99. Neither T nor Fwhole is the original full-filter F.',
         'selection':'fixed comparisons, no outer performance selection; backward reproduces notebook single-fit top-k, not sequential elimination',
         'tasks':args.tasks,'folds':5,'repeats':3,'seed':20260909,'jobs':args.jobs,'fold_limit':args.fold_limit,
         'sources':{str(p.relative_to(ROOT)):mb.digest(p) for p in [Path(__file__),Path(v6.__file__),Path(mb.__file__),oldpath,newpath,ROOT/'data/20231124is0155-500NIHSSfinal.xlsx']},
         'versions':{m.__name__:m.__version__ for m in [np,pd,sklearn,lightgbm,xgboost,catboost]},
         'shap':'LightGBM held-out pred_contrib, tree-path-dependent, raw log odds; training tree paths include training-only SMOTE'}
    cp=out/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=cfg: raise ValueError('Resume configuration mismatch')
    if not cp.exists():mb.write_json(cp,cfg)
    if args.workers>1:
        assert args.worker_task is None and args.fold_limit is None
        def worker(item):
            task,split=item
            cmd=[sys.executable,str(Path(__file__)), '--jobs',str(args.jobs),'--tasks',*args.tasks,
                 '--output',str(out),'--worker-task',task,'--split-index',str(split)]
            with (out/f'{task}_{split:03d}.log').open('w') as log:
                subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
            print(task,split,'all methods complete',flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(worker,[(task,split) for task in args.tasks for split in range(15)]))
        mb.write_json(out/'COMPLETE.json',{'tasks':args.tasks})
        return
    clinical=mb.load_clinical();mb.valid_ids(clinical.ID.tolist())
    for task in args.tasks:
        if args.worker_task is not None and task!=args.worker_task:continue
        target=mb.TASKS[task].target_column
        frame=clinical.loc[clinical.ID.isin(new['ids'])&clinical[target].notna()].copy()
        ids=frame.ID.to_numpy();y=mb.make_target(frame[target],task)
        assert set(np.unique(y))=={0,1}
        lookup={s:i for i,s in enumerate(new['ids'])};rows=np.array([lookup[s] for s in ids])
        c,cn=mb.clinical_matrix(clinical,ids,'core',task)
        assert not any('mrs' in n.lower() or 'nihss' in n.lower() for n in cn)
        blocks={'C':(c,cn),'R':(old['R_unclip_p99'][rows],old['global_names']),
                'S':(old['burden'][rows],old['burden_names']),'T':(new['T'][rows],new['T_names']),
                'Fwhole':(old['F_unclip_p99'][rows],old['names'])}
        mask=np.array([n.rsplit('|',1)[-1] in FLAGS for n in new['T_names']])
        blocks['flags']=(new['T'][rows][:,mask],new['T_names'][mask])
        cohort={'task':task,'n':len(y),'positive':int(y.sum()),'negative':int((1-y).sum()),
                'missing_target_in_image_cohort':int((clinical.ID.isin(new['ids'])&clinical[target].isna()).sum())}
        if args.split_index is None or args.split_index==0:mb.write_json(out/f'{task}_cohort.json',cohort)
        desc=pd.DataFrame(c,columns=cn).describe().T;desc['missing']=np.isnan(c).sum(axis=0)
        if args.split_index is None or args.split_index==0:desc.to_csv(out/f'{task}_clinical_description.csv')
        for split,(tr,te) in enumerate(RepeatedStratifiedKFold(n_splits=5,n_repeats=3,random_state=20260909).split(c,y)):
            if args.fold_limit is not None and split>=args.fold_limit:break
            if args.split_index is not None and split!=args.split_index:continue
            assert not set(ids[tr])&set(ids[te])
            seed=20260909+split+1000*list(mb.TASKS).index(task)
            dest=out/f'{task}_{split:03d}';dest.mkdir(exist_ok=True)
            audit={'train_ids':ids[tr].tolist(),'test_ids':ids[te].tolist(),'y_test':y[te].tolist(),'split':split,'seed':seed}
            mb.write_json(dest/'split.json',audit)
            cache={}
            def signature(block,kind='LightGBM',selector='RFE'):
                key=(block,kind,selector)
                if key not in cache:
                    x,n=blocks[block]
                    a,b,nn,prep=mb.fit_train_preprocessor(x[tr],x[te],n)
                    a,b=mb.standardize(a,b)
                    chosen,trace=select(a,y[tr],nn,kind,selector,mb.TASKS[task].final_features,seed,args.jobs)
                    cache[key]=(a[:,chosen],b[:,chosen],np.array([block+'::'+str(k) for k in nn[chosen]]),{'preprocessing':prep,'trace':trace})
                return cache[key]
            def run(name,parts,kind='LightGBM',selector='RFE'):
                fp=dest/(name+'.json')
                if fp.exists():return
                chunks=[signature(part,kind,selector) for part in parts]
                a=np.hstack([s[0] for s in chunks]);b=np.hstack([s[1] for s in chunks]);nn=np.concatenate([s[2] for s in chunks])
                p,sh=fit(a,b,y[tr],nn,kind,seed,args.jobs,kind=='LightGBM')
                if sh is not None:
                    np.savez_compressed(dest/(name+'_shap.npz'),ids=ids[te],names=nn,values=b,shap=sh[:,:-1],expected=sh[:,-1])
                mb.write_json(fp,{'method':name,'features':nn.tolist(),'probability':p.tolist(),'preprocessing':[s[3] for s in chunks]})
                print(task,split,name,'complete',flush=True)
            for kind in ['LightGBM','RF','XGBoost','CatBoost']:
                for block in ['C','R','S','T','Fwhole']:run(f'{kind}_RFE_{block}',[block],kind)
            for block in ['C','R','S','T','Fwhole']:run('LightGBM_backward_'+block,[block],selector='backward')
            for regional in ['T','Fwhole']:
                for other in [('C',),('R',),('S',),('C','R'),('C','R','S')]:
                    parts=(regional,)+other;run('LightGBM_RFE_'+'_'.join(parts),parts)
            run('LightGBM_RFE_flags',['flags'])
        print(task,'finished',flush=True)
    if args.worker_task is None:
        mb.write_json(out/('SMOKE_COMPLETE.json' if args.fold_limit is not None else 'COMPLETE.json'),{'tasks':args.tasks})
if __name__=='__main__':main()
