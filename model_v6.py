#!/usr/bin/env python3
"""Fixed paired exploratory comparisons; training-only feature processing."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;INPUTS=HERE/'inputs'
import modeling_base as mb

def preprocess(x,names,tr,te):
    a,b,nn,audit=mb.fit_train_preprocessor(x[tr],x[te],names)
    a,b=mb.standardize(a,b)
    return a,b,nn,audit

def alternative(model,a,b,y,seed,jobs):
    if model=='ET':est=ExtraTreesClassifier(n_estimators=400,min_samples_leaf=5,max_features=.7,
                                           class_weight='balanced',random_state=seed,n_jobs=jobs)
    else:est=LogisticRegression(C=.1,class_weight='balanced',max_iter=5000,solver='liblinear',random_state=seed)
    est.fit(a,y)
    return est.predict_proba(b)[:,1]

def task_run(task,data,old,clinical,p,args):
    out=args.output;spec=p['modelling'];target=mb.TASKS[task].target_column
    usable=clinical.loc[clinical[target].notna()&clinical.ID.isin(data['ids'])].copy()
    ids=usable.ID.to_numpy();mb.valid_ids(ids);y=mb.make_target(usable[target],task)
    lookup={s:i for i,s in enumerate(data['ids'])};rows=np.array([lookup[s] for s in ids])
    cx,cn=mb.clinical_matrix(clinical,ids,'core',task)
    small=np.column_stack([cx,np.log1p(old['volume'][rows]),old['global_available'][rows]])
    sn=np.r_[cn,['log_lesion_volume_mm3','global_texture_available']]
    # Matched alternative baselines use ALL burden/global columns, without supervised filtering.
    fullbase=np.hstack([small,old['burden'][rows],old['R_unclip_p99'][rows]])
    bn=np.r_[sn,old['burden_names'],old['global_names']]
    detail=[];audit=[]
    cv=RepeatedStratifiedKFold(n_splits=spec['folds'],n_repeats=spec['repeats'],random_state=spec['seed'])
    for split,(tr,te) in enumerate(cv.split(cx,y)):
        seed=spec['seed']+split+1000*list(mb.TASKS).index(task)
        cp=out/f'{task}_split_{split:03d}.json'
        if cp.exists():
            saved=json.loads(cp.read_text())
            assert saved['audit']['train_ids']==ids[tr].tolist() and saved['audit']['test_ids']==ids[te].tolist()
        else:
            st,sv,_,_=preprocess(small,sn,tr,te)
            bt,bv,ba=mb.signature(old['burden'][rows][tr],old['burden'][rows][te],old['burden_names'],y[tr],20,seed,args.jobs)
            rt,rv,ra=mb.signature(old['R_unclip_p99'][rows][tr],old['R_unclip_p99'][rows][te],old['global_names'],y[tr],spec[task+'_features'],seed,args.jobs)
            bt,bv=np.hstack([st,bt,rt]),np.hstack([sv,bv,rv])
            at,av,_,aa=preprocess(fullbase,bn,tr,te)
            pred={'V1_base':mb.fit_predict(bt,bv,y[tr],seed,args.jobs)}
            for model in ['ET','LR']:pred[model+'_base']=alternative(model,at,av,y[tr],seed,args.jobs)
            fits={}
            for block in p['blocks']:
                x=data[block][rows];names=data[block+'_names']
                ft,fv,nn,fa=preprocess(x,names,tr,te)
                chosen,trace=mb.staged_rfe(ft,y[tr],nn,spec[task+'_features'],seed,args.jobs)
                fs,vs=ft[:,chosen],fv[:,chosen]
                pred['V1_'+block]=mb.fit_predict(fs,vs,y[tr],seed,args.jobs)
                pred['V1_base_'+block]=mb.fit_predict(np.hstack([bt,fs]),np.hstack([bv,vs]),y[tr],seed,args.jobs)
                pred['V1_blend_'+block]=.75*pred['V1_base']+.25*pred['V1_'+block]
                for model in ['ET','LR']:
                    pred[model+'_'+block]=alternative(model,ft,fv,y[tr],seed,args.jobs)
                    pred[model+'_base_'+block]=alternative(model,np.hstack([at,ft]),np.hstack([av,fv]),y[tr],seed,args.jobs)
                fits[block]={'preprocessing':fa,'selected_features':nn[chosen].tolist(),'rfe_remaining':trace}
            # Reproduce V5 predictions exactly, so comparisons cannot silently change the reference.
            prior=json.loads((INPUTS/'v5_reference_predictions'/f'{task}_split_{split:03d}.json').read_text())
            assert prior['audit']['test_ids']==ids[te].tolist() and prior['audit']['train_ids']==ids[tr].tolist()
            np.testing.assert_allclose(pred['V1_base'],prior['predictions']['B_R__unclip_p99'],atol=1e-12,rtol=0)
            saved={'audit':{'task':task,'split':split,'train_ids':ids[tr].tolist(),'test_ids':ids[te].tolist(),
                   'seed':seed,'burden':ba,'global':ra,'alternative_base_preprocessing':aa,'blocks':fits,
                   'v5_baseline_reproduced':True,'outer_selection':'none'},
                   'predictions':{k:v.tolist() for k,v in pred.items()}}
            mb.write_json(cp,saved)
        audit.append(saved['audit'])
        repeat,fold=divmod(split,spec['folds'])
        for method,probs in saved['predictions'].items():
            detail.extend({'ID':str(ids[i]),'y':int(y[i]),'repeat':repeat,'fold':fold,'method':method,'probability':float(prob)}
                          for i,prob in zip(te,probs))
        print(task,split+1,'/',spec['folds']*spec['repeats'],'complete',flush=True)
    detailed=pd.DataFrame(detail);detailed.to_csv(out/f'{task}_oof_by_repeat.csv',index=False)
    avg=detailed.groupby(['ID','method']).probability.mean().unstack().loc[ids]
    assert detailed.groupby(['ID','method']).size().eq(spec['repeats']).all()
    predictions={m:avg[m].to_numpy() for m in avg.columns}
    summary=mb.uncertainty(y,predictions,'V1_base',2000,spec['seed'])
    rm=[]
    for (repeat,method),g in detailed.groupby(['repeat','method']):
        assert len(g)==len(ids)==g.ID.nunique()
        rm.append(dict(repeat=int(repeat),method=method,**mb.metrics(g.y,g.probability)))
    rm=pd.DataFrame(rm)
    for row in summary:
        r=rm[rm.method==row['method']].auroc
        row.update(task=task,n=len(y),positives=int(y.sum()),repeat_auroc_mean=float(r.mean()),
                   repeat_auroc_min=float(r.min()),repeat_auroc_max=float(r.max()))
    pairs=[]
    for model in ['V1','ET','LR']:
        ref=model+'_base';methods=[m for m in predictions if m.startswith(model+'_base_') or m.startswith(model+'_blend_')]
        pairs.extend(mb.uncertainty(y,{ref:predictions[ref],**{m:predictions[m] for m in methods}},ref,2000,spec['seed']))
    pd.DataFrame({'ID':ids,'y':y,**predictions}).to_csv(out/f'{task}_oof_predictions.csv',index=False)
    pd.DataFrame(summary).to_csv(out/f'{task}_summary.csv',index=False)
    pd.DataFrame(pairs).to_csv(out/f'{task}_matched_comparisons.csv',index=False)
    rm.to_csv(out/f'{task}_repeat_metrics.csv',index=False);mb.write_json(out/f'{task}_audit.json',audit)
    print(pd.DataFrame(summary)[['method','auroc','delta_auc','repeat_auroc_mean']].to_string(index=False),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--jobs',type=int,default=6)
    parser.add_argument('--output',type=Path,default=HERE/'results/lesion_20260909');args=parser.parse_args()
    p=json.loads((HERE/'config/protocol.json').read_text());fm=json.loads((HERE/'output/features_manifest.json').read_text())
    assert mb.digest(HERE/'output/features.npz')==fm['archive_sha256']
    assert mb.digest(HERE/'config/protocol.json')==fm['protocol_sha256']
    oldpath=INPUTS/'v5_features.npz'
    assert mb.digest(oldpath)==fm['v5_archive_sha256']
    import sklearn,lightgbm
    config={'protocol':p,'features_manifest':fm,'jobs':args.jobs,
            'versions':{'sklearn':sklearn.__version__,'lightgbm':lightgbm.__version__,'numpy':np.__version__,'pandas':pd.__version__},
            'sources':{str(f.relative_to(ROOT)):mb.digest(f) for f in [Path(__file__),Path(mb.__file__),ROOT/'data/20231124is0155-500NIHSSfinal.xlsx']}}
    args.output.mkdir(parents=True,exist_ok=True);cp=args.output/'config.json'
    if cp.exists() and json.loads(cp.read_text())!=config:raise ValueError('Model resume mismatch')
    mb.write_json(cp,config)
    data=np.load(HERE/'output/features.npz',allow_pickle=False);old=np.load(oldpath,allow_pickle=False)
    np.testing.assert_array_equal(data['ids'],old['ids'])
    clinical=mb.load_clinical();mb.valid_ids(clinical.ID.tolist())
    for task in p['modelling']['tasks']:task_run(task,data,old,clinical,p,args)
    mb.write_json(args.output/'COMPLETE.json',{'tasks':p['modelling']['tasks']})
if __name__=='__main__':main()
