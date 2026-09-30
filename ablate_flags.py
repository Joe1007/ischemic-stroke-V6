#!/usr/bin/env python3
"""Post-hoc diagnostic: does V1_base_T gain come from non-radiomic flags?"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import model_v6 as m
mb=m.mb
HERE=m.HERE;ROOT=m.ROOT;INPUTS=m.INPUTS
OUT=HERE/'results/flags_ablation_20260910'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data=np.load(HERE/'output/features.npz');old=np.load(INPUTS/'v5_features.npz')
    p=json.loads((HERE/'config/protocol.json').read_text());clinical=mb.load_clinical()
    config={'purpose':'Post-hoc mechanistic diagnostic after mRS T selected 299 flags / 300 slots; not independent validation',
        'model':'V1 base plus separately training-RFE-selected flags, same keep count and splits',
        'sources':{str(f.relative_to(ROOT)):mb.digest(f) for f in [Path(__file__),Path(m.__file__),HERE/'output/features.npz']}}
    cp=OUT/'config.json'
    if cp.exists():assert json.loads(cp.read_text())==config
    mb.write_json(cp,config)
    # mRS only: targeted to the observed improvement, not a new model search.
    task='mRS';target=mb.TASKS[task].target_column
    usable=clinical.loc[clinical[target].notna()&clinical.ID.isin(data['ids'])].copy()
    ids=usable.ID.to_numpy();y=mb.make_target(usable[target],task)
    lookup={s:i for i,s in enumerate(data['ids'])};rows=np.array([lookup[s] for s in ids])
    cx,cn=mb.clinical_matrix(clinical,ids,'core',task)
    small=np.column_stack([cx,np.log1p(old['volume'][rows]),old['global_available'][rows]])
    sn=np.r_[cn,['log_lesion_volume_mm3','global_texture_available']]
    flags=data['T'][rows,-664:];fn=data['T_names'][-664:]
    cv=m.RepeatedStratifiedKFold(n_splits=5,n_repeats=3,random_state=p['modelling']['seed'])
    detail=[]
    for split,(tr,te) in enumerate(cv.split(cx,y)):
        seed=p['modelling']['seed']+split;path=OUT/f'mRS_split_{split:03d}.json'
        prior=json.loads((HERE/'results/lesion_20260909'/path.name).read_text())
        assert prior['audit']['train_ids']==ids[tr].tolist() and prior['audit']['test_ids']==ids[te].tolist()
        if path.exists():saved=json.loads(path.read_text())
        else:
            st,sv,_,_=m.preprocess(small,sn,tr,te)
            bt,bv,_=mb.signature(old['burden'][rows][tr],old['burden'][rows][te],old['burden_names'],y[tr],20,seed,8)
            rt,rv,_=mb.signature(old['R_unclip_p99'][rows][tr],old['R_unclip_p99'][rows][te],old['global_names'],y[tr],20,seed,8)
            ft,fv,fa=mb.signature(flags[tr],flags[te],fn,y[tr],20,seed,8)
            prob=mb.fit_predict(np.hstack([st,bt,rt,ft]),np.hstack([sv,bv,rv,fv]),y[tr],seed,8)
            saved={'train_ids':ids[tr].tolist(),'test_ids':ids[te].tolist(),'selection':fa,'probability':prob.tolist()}
            mb.write_json(path,saved)
        detail.extend({'ID':str(ids[i]),'y':int(y[i]),'repeat':split//5,'fold':split%5,'flags':float(a),'T':float(b),'base':float(c)}
            for i,a,b,c in zip(te,saved['probability'],prior['predictions']['V1_base_T'],prior['predictions']['V1_base']))
        print('flags mRS',split+1,'/15',flush=True)
    d=pd.DataFrame(detail);d.to_csv(OUT/'mRS_oof_by_repeat.csv',index=False)
    avg=d.groupby('ID')[['flags','T','base']].mean().loc[ids]
    summary=mb.uncertainty(y,{k:avg[k].to_numpy() for k in avg},'flags',2000,p['modelling']['seed'])
    pd.DataFrame(summary).to_csv(OUT/'mRS_summary.csv',index=False)
    mb.write_json(OUT/'COMPLETE.json',{'task':'mRS','post_hoc':True})
    print(pd.DataFrame(summary).to_string(index=False),flush=True)
if __name__=='__main__':main()
