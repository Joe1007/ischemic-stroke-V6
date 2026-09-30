#!/usr/bin/env python3
"""reproduce_tables.py (36 methods x 4 tasks) with S/T/Fwhole/flags taken from another atlas.

Identical splits, seeds, preprocessing, selection, classifiers and SHAP as the AAL3 run in
../tables_20260917. C and R do not depend on the atlas, so their methods must reproduce the
AAL3 predictions exactly (checked by report_atlas.py).
"""
import argparse, json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
HERE=Path(__file__).resolve().parent;FINAL=HERE.parent;ROOT=FINAL.parents[1];INPUTS=FINAL.parent/'inputs'
sys.path.insert(0,str(FINAL))
import reproduce_tables as rt
mb=rt.mb;FLAGS=rt.FLAGS

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--atlas',required=True)
    ap.add_argument('--jobs',type=int,default=2);ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--tasks',nargs='+',default=list(mb.TASKS))
    ap.add_argument('--worker-task');ap.add_argument('--split-index',type=int);args=ap.parse_args()
    out=HERE/f'tables_{args.atlas}';out.mkdir(parents=True,exist_ok=True)
    oldpath=INPUTS/'v5_features.npz';newpath=HERE/f'features_{args.atlas}.npz'
    old=np.load(oldpath);new=np.load(newpath);np.testing.assert_array_equal(old['ids'],new['ids'])
    cfg={'atlas':args.atlas,'base_script':'final/reproduce_tables.py','adaptation':'C, R unchanged (V5); S, T, Fwhole, flags from '+args.atlas,
         'tasks':args.tasks,'folds':5,'repeats':3,'seed':20260909,'jobs':args.jobs,'fold_limit':None,
         'sources':{str(p.relative_to(ROOT)):mb.digest(p) for p in [Path(__file__),Path(rt.__file__),Path(rt.v6.__file__),Path(mb.__file__),oldpath,newpath,ROOT/'data/20231124is0155-500NIHSSfinal.xlsx']}}
    cp=out/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=cfg:raise ValueError('Resume configuration mismatch')
    if not cp.exists():mb.write_json(cp,cfg)
    if args.workers>1:
        def worker(item):
            task,split=item
            cmd=[sys.executable,str(Path(__file__)),'--atlas',args.atlas,'--jobs',str(args.jobs),'--tasks',*args.tasks,
                 '--worker-task',task,'--split-index',str(split)]
            with (out/f'{task}_{split:03d}.log').open('w') as log:
                subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
            print(args.atlas,task,split,'all methods complete',flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(worker,[(t,s) for t in args.tasks for s in range(15)]))
        mb.write_json(out/'COMPLETE.json',{'tasks':args.tasks});return
    clinical=mb.load_clinical();mb.valid_ids(clinical.ID.tolist())
    for task in args.tasks:
        if args.worker_task is not None and task!=args.worker_task:continue
        target=mb.TASKS[task].target_column
        frame=clinical.loc[clinical.ID.isin(new['ids'])&clinical[target].notna()].copy()
        ids=frame.ID.to_numpy();y=mb.make_target(frame[target],task)
        lookup={s:i for i,s in enumerate(new['ids'])};rows=np.array([lookup[s] for s in ids])
        c,cn=mb.clinical_matrix(clinical,ids,'core',task)
        blocks={'C':(c,cn),'R':(old['R_unclip_p99'][rows],old['global_names']),
                'S':(new['burden'][rows],new['burden_names']),'T':(new['T'][rows],new['T_names']),
                'Fwhole':(new['F_unclip_p99'][rows],new['names'])}
        mask=np.array([n.rsplit('|',1)[-1] in FLAGS for n in new['T_names']])
        blocks['flags']=(new['T'][rows][:,mask],new['T_names'][mask])
        if args.split_index in (None,0):
            mb.write_json(out/f'{task}_cohort.json',{'task':task,'n':len(y),'positive':int(y.sum()),'negative':int((1-y).sum())})
        for split,(tr,te) in enumerate(RepeatedStratifiedKFold(n_splits=5,n_repeats=3,random_state=20260909).split(c,y)):
            if args.split_index is not None and split!=args.split_index:continue
            seed=20260909+split+1000*list(mb.TASKS).index(task)
            ref=json.loads((FINAL/'tables_20260917'/f'{task}_{split:03d}'/'split.json').read_text())
            assert ref['train_ids']==ids[tr].tolist() and ref['test_ids']==ids[te].tolist() and ref['seed']==seed
            dest=out/f'{task}_{split:03d}';dest.mkdir(exist_ok=True)
            mb.write_json(dest/'split.json',{'train_ids':ids[tr].tolist(),'test_ids':ids[te].tolist(),'y_test':y[te].tolist(),'split':split,'seed':seed})
            cache={}
            def signature(block,kind='LightGBM',selector='RFE'):
                key=(block,kind,selector)
                if key not in cache:
                    x,n=blocks[block]
                    a,b,nn,prep=mb.fit_train_preprocessor(x[tr],x[te],n)
                    a,b=mb.standardize(a,b)
                    chosen,trace=rt.select(a,y[tr],nn,kind,selector,mb.TASKS[task].final_features,seed,args.jobs)
                    cache[key]=(a[:,chosen],b[:,chosen],np.array([block+'::'+str(k) for k in nn[chosen]]),{'preprocessing':prep,'trace':trace})
                return cache[key]
            def run(name,parts,kind='LightGBM',selector='RFE'):
                fp=dest/(name+'.json')
                if fp.exists():return
                chunks=[signature(part,kind,selector) for part in parts]
                a=np.hstack([s[0] for s in chunks]);b=np.hstack([s[1] for s in chunks]);nn=np.concatenate([s[2] for s in chunks])
                p,sh=rt.fit(a,b,y[tr],nn,kind,seed,args.jobs,kind=='LightGBM')
                if sh is not None:
                    np.savez_compressed(dest/(name+'_shap.npz'),ids=ids[te],names=nn,values=b,shap=sh[:,:-1],expected=sh[:,-1])
                mb.write_json(fp,{'method':name,'features':nn.tolist(),'probability':p.tolist(),'preprocessing':[s[3] for s in chunks]})
                print(args.atlas,task,split,name,'complete',flush=True)
            for kind in ['LightGBM','RF','XGBoost','CatBoost']:
                for block in ['C','R','S','T','Fwhole']:run(f'{kind}_RFE_{block}',[block],kind)
            for block in ['C','R','S','T','Fwhole']:run('LightGBM_backward_'+block,[block],selector='backward')
            for regional in ['T','Fwhole']:
                for other in [('C',),('R',),('S',),('C','R'),('C','R','S')]:
                    parts=(regional,)+other;run('LightGBM_RFE_'+'_'.join(parts),parts)
            run('LightGBM_RFE_flags',['flags'])

if __name__=='__main__':main()
