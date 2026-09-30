#!/usr/bin/env python3
"""Audit and report completed V6-adapted manuscript tables and held-out SHAP."""
import argparse,json,re,sys,shutil,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
from sklearn.metrics import roc_auc_score,roc_curve
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent))
import model_v6 as v6
mb=v6.mb
FLAGS={'affected','firstorder_available','contrast_available','texture_available'}

def region(name):
    block,name=name.split('::',1)
    if block=='T' or block=='flags': return int(re.match(r'AAL3_(\d+)',name).group(1))
    if block=='Fwhole': return int(name.rsplit('_',1)[-1])
    return None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=HERE/'tables_20260917')
    ap.add_argument('--output',type=Path,default=HERE/'output')
    ap.add_argument('--bootstraps',type=int,default=2000);args=ap.parse_args();out=args.input
    assert (out/'COMPLETE.json').exists(),'Full run not complete; do not publish smoke output'
    cfg=json.loads((out/'config.json').read_text());assert cfg['fold_limit'] is None
    report=args.output;report.mkdir(parents=True,exist_ok=True)
    atlas=nib.load(str(ROOT/'data/MNI/AAL3_MNI.nii.gz'));labels=np.asarray(atlas.dataobj)
    anatomical_support=np.asarray(nib.as_closest_canonical(atlas).dataobj)>0
    labelnames={int(s.split(maxsplit=1)[0]):s.split(maxsplit=1)[1].strip() for s in (ROOT/'data/MNI/AAL3.txt').read_text().splitlines() if s.strip()}
    allsummary=[];audits={};cohorts=[]
    for task in cfg['tasks']:
        details=[];attributions={};values={};frequencies={};nfold=0
        cohort=json.loads((out/f'{task}_cohort.json').read_text());cohorts.append(cohort)
        for dest in sorted(out.glob(task+'_[0-9][0-9][0-9]')):
            a=json.loads((dest/'split.json').read_text());assert not set(a['train_ids'])&set(a['test_ids'])
            assert len(set(a['train_ids'])|set(a['test_ids']))==cohort['n'];nfold+=1
            methods=[]
            for fp in sorted(dest.glob('*.json')):
                if fp.name=='split.json':continue
                m=json.loads(fp.read_text());name=m['method'];methods.append(name)
                assert len(m['probability'])==len(a['test_ids'])
                assert np.isfinite(m['probability']).all() and np.all((np.array(m['probability'])>=0)&(np.array(m['probability'])<=1))
                details.extend({'ID':i,'y':y,'repeat':a['split']//5,'fold':a['split']%5,'method':name,'probability':p} for i,y,p in zip(a['test_ids'],a['y_test'],m['probability']))
                shp=dest/(name+'_shap.npz')
                if shp.exists():
                    # These archives were generated locally by reproduce_tables.py;
                    # pandas patient IDs were saved as object arrays in that run.
                    z=np.load(shp,allow_pickle=True);np.testing.assert_array_equal(z['ids'],a['test_ids'])
                    logits=z['shap'].sum(axis=1)+z['expected']
                    np.testing.assert_allclose(1/(1+np.exp(-logits)),m['probability'],atol=1e-8)
                    sf=pd.DataFrame(z['shap'],columns=z['names'],index=z['ids']);sf.index.name='ID'
                    vf=pd.DataFrame(z['values'],columns=z['names'],index=z['ids']);vf.index.name='ID'
                    attributions.setdefault(name,[]).append(sf);values.setdefault(name,[]).append(vf)
                    frequencies.setdefault(name,[]).extend(z['names'].tolist())
            assert len(methods)==36,(dest,len(methods))
        assert nfold==15
        detail=pd.DataFrame(details);assert detail.groupby(['ID','method']).size().eq(3).all()
        assert detail.groupby(['ID','method','repeat']).size().eq(1).all()
        assert detail.groupby('ID').y.nunique().eq(1).all()
        detail.to_csv(report/f'{task}_oof_by_repeat.csv',index=False)
        avg=detail.groupby(['ID','method']).probability.mean().unstack()
        y=detail.groupby('ID').y.first().loc[avg.index]
        avg.insert(0,'y',y);avg.to_csv(report/f'{task}_oof_predictions.csv')
        preds={m:avg[m].to_numpy() for m in avg.columns if m!='y'}
        summary=pd.DataFrame(mb.uncertainty(y.to_numpy(),preds,'LightGBM_RFE_C',args.bootstraps,cfg['seed']))
        summary['task']=task;summary['n']=len(avg)
        rm=pd.DataFrame([{'method':m,'repeat':r,**mb.metrics(g.y,g.probability)} for (m,r),g in detail.groupby(['method','repeat'])])
        summary['repeat_auroc_mean']=summary.method.map(rm.groupby('method').auroc.mean())
        rm.to_csv(report/f'{task}_repeat_metrics.csv',index=False);summary.to_csv(report/f'{task}_summary.csv',index=False)
        allsummary.append(summary)
        fig,ax=plt.subplots(figsize=(6,5))
        for name in ['LightGBM_RFE_C','LightGBM_RFE_R','LightGBM_RFE_T','LightGBM_RFE_Fwhole','LightGBM_RFE_T_C_R_S']:
            fpr,tpr,_=roc_curve(y,preds[name]);ax.plot(fpr,tpr,label=f'{name.replace("LightGBM_RFE_", "")} ({roc_auc_score(y,preds[name]):.3f})')
        ax.plot([0,1],[0,1],'k--',alpha=.4);ax.set(xlabel='False positive rate',ylabel='True positive rate',title=task+' averaged OOF ROC');ax.legend(fontsize=8)
        fig.tight_layout();fig.savefig(report/f'{task}_ROC.png',dpi=180);plt.close(fig)
        for name,parts in attributions.items():
            sh=pd.concat(parts).fillna(0);val=pd.concat(values[name])
            assert sh.index.value_counts().eq(3).all()
            # Zero contribution when a feature is not selected, denominated over all held-out visits.
            importance=sh.abs().mean().sort_values(ascending=False)
            freq=pd.Series(frequencies[name]).value_counts()/15
            tab=pd.DataFrame({'mean_abs_heldout_shap':importance,'selection_fraction':freq}).fillna(0).sort_values('mean_abs_heldout_shap',ascending=False)
            tab['region']=[region(n) for n in tab.index]
            tab['kind']=['availability' if n.rsplit('|',1)[-1] in FLAGS else 'radiomics' if n.startswith(('T::','Fwhole::','R::')) else 'clinical_or_burden' for n in tab.index]
            tab.to_csv(report/f'{task}_{name}_importance.csv',index_label='feature')
            top=importance.head(20).index.tolist()[::-1];fig,ax=plt.subplots(figsize=(12,8));rng=np.random.default_rng(0)
            for j,feat in enumerate(top):
                color=val[feat].to_numpy();finite=np.isfinite(color)
                lo,hi=np.nanpercentile(color,[5,95]) if finite.any() else (0,1)
                color=np.where(finite,np.clip((color-lo)/(hi-lo+1e-8),0,1),.5)
                ax.scatter(sh[feat],j+rng.uniform(-.25,.25,len(sh)),c=color,cmap='coolwarm',s=8,alpha=.55,vmin=0,vmax=1)
            ax.set_yticks(range(len(top)),top,fontsize=7);ax.set(xlabel='Held-out SHAP (raw log odds); red=high, blue=low fold-scaled value',title=task+' '+name)
            fig.tight_layout();fig.savefig(report/f'{task}_{name}_SHAP.png',dpi=160);plt.close(fig)
            for kind in ['radiomics','availability']:
                rg=tab[tab.kind.eq(kind)&tab.region.notna()].groupby('region').mean_abs_heldout_shap.sum()
                if rg.empty:continue
                rtab=rg.rename('mean_abs_shap_sum').to_frame();rtab['region_name']=[labelnames.get(int(k),'unknown') for k in rtab.index]
                rtab.sort_values('mean_abs_shap_sum',ascending=False).to_csv(report/f'{task}_{name}_{kind}_regions.csv')
                arr=np.zeros(labels.shape,dtype=np.float32)
                for k,w in rg.items():arr[labels==k]=w
                nif=nib.Nifti1Image(arr,atlas.affine,atlas.header.copy());nif.set_data_dtype(np.float32)
                nib.save(nif,str(report/f'{task}_{name}_{kind}_atlas.nii.gz'))
                canon=nib.as_closest_canonical(nif);data=np.asarray(canon.dataobj)
                zs=np.linspace(int(data.shape[2]*.25),int(data.shape[2]*.8),8).astype(int)
                fig,axs=plt.subplots(2,4,figsize=(12,6));maxval=max(float(data.max()),1e-12)
                for ax,z in zip(axs.flat,zs):
                    ax.imshow(anatomical_support[:,:,z].T,origin='lower',cmap='gray',vmin=0,vmax=1,alpha=.2)
                    im=ax.imshow(np.ma.masked_where(data[:,:,z].T==0,data[:,:,z].T),origin='lower',cmap='magma',vmin=0,vmax=maxval)
                    zmm=nib.affines.apply_affine(canon.affine,[0,0,z])[2]
                    ax.set_title(f'z={zmm:.0f} mm');ax.axis('off')
                    ax.text(0,.5,'L',transform=ax.transAxes);ax.text(.95,.5,'R',transform=ax.transAxes)
                fig.suptitle(task+' '+name+' '+kind+' (exploratory absolute SHAP)')
                fig.colorbar(im,ax=list(axs.flat),shrink=.6,label='Sum of mean absolute held-out SHAP')
                fig.savefig(report/f'{task}_{name}_{kind}_atlas.png',dpi=150);plt.close(fig)
        audits[task]={'n':len(avg),'folds':nfold,'methods':36,'oof_per_patient':3,'shap_probability_reconstruction':True,'disjoint_patient_ids':True}
    summary=pd.concat(allsummary);summary.to_csv(report/'all_metrics.csv',index=False)
    pd.DataFrame(cohorts).to_csv(report/'cohort_counts.csv',index=False)
    # Manuscript-style tables explicitly retain adapted feature names.
    for task in cfg['tasks']:
        sub=summary[summary.task.eq(task)].set_index('method')
        modeltable=pd.DataFrame({b:{m:sub.loc[f'{m}_RFE_{b}','auroc'] for m in ['RF','XGBoost','CatBoost','LightGBM']} for b in ['C','R','S','T','Fwhole']})
        modeltable.to_csv(report/f'{task}_table_4_2_to_4_5.csv')
    summary[summary.method.str.contains('LightGBM_(RFE|backward)_(C|R|S|T|Fwhole)$',regex=True)].to_csv(report/'table_4_6_selection.csv',index=False)
    summary[summary.method.str.match('LightGBM_RFE_(T|Fwhole)(_|$)')].to_csv(report/'table_4_7_combinations.csv',index=False)
    mb.write_json(report/'independent_audit.json',audits)
    lines=['# V6 核心的論文結果重現（方法有調整）','',
      '本輪為同 cohort 的探索性重跑。T 是 V6 lesion∩AAL3 加 availability；Fwhole 是 V5 完整腦區 Original 特徵對照。兩者都不能標示為原稿全濾波器 F 的逐字重現。C 使用保守臨床欄位；S 僅 AAL3 burden；R 為重抽 Original 病灶 radiomics。','',
      '表 4.2–4.5 對應各任務的 classifier CSV；表 4.6 為篩選法比較；表 4.7 為固定組合。backward 依原 notebook 實作為 single-fit top-k。所有比較固定列出，不以 outer 分數挑選後宣稱獨立驗證。','',
      'SHAP 為 held-out TreeSHAP、raw log odds，使用包含 training-only SMOTE 的訓練樹路徑背景。每人三次 held-out 解釋；未入選特徵貢獻設零。腦區圖把同區特徵的平均絕對 SHAP 加總，availability 與 radiomics 分開，不代表因果或經確認的生物機制。','',
      'AUC 使用每人三次 OOF 平均；另提供 mean-repeat AUC。bootstrap 為固定 OOF 的配對病人抽樣，不包含重訓或跨版本選擇不確定性。','',
      '| 任務 | n | T AUC | whole-region F AUC | C+R+S+T AUC |','|---|---:|---:|---:|---:|']
    for task in cfg['tasks']:
        s=summary[summary.task.eq(task)].set_index('method')
        lines.append(f"| {task} | {int(s.n.iloc[0])} | {s.loc['LightGBM_RFE_T','auroc']:.3f} | {s.loc['LightGBM_RFE_Fwhole','auroc']:.3f} | {s.loc['LightGBM_RFE_T_C_R_S','auroc']:.3f} |")
    lines+=['','## 已完成的 V6 核心精確重跑','']
    core=HERE/'output/core'
    assert (core/'reproduction_audit.json').exists(),'Run audit_core.py before packaging'
    (report/'core').mkdir(exist_ok=True)
    if core.resolve()!=(report/'core').resolve():
        for p in core.iterdir():
            if p.is_file() and p.suffix in {'.csv','.json','.md'}:shutil.copy2(p,report/'core'/p.name)
    lines+=['[核心結果與配對增量](core/RESULTS.md)｜[核心指標 CSV](core/core_results.csv)｜[精確重跑稽核](core/reproduction_audit.json)',
            '', '核心 V1_base_T 與此處 T+C+R+S 不是相同模型：核心另含體積、global availability，且各分支保留數不同，不可混用其數字或解釋圖。',
            '', '## 指標與下載','',
            '- [全部方法指標](all_metrics.csv)：AUROC、AUPRC（average precision）、accuracy、sensitivity、specificity、Brier、AUROC 95% CI 與 mean-repeat AUROC。',
            '- 分類指標固定使用 0.5 閾值；未使用測試結局挑閾值。ΔAUC 欄位以該表 reference 指定的模型比較，不能視為所有模型間的增量。',
            '- [樣本數與標籤分布](cohort_counts.csv)｜[表 4.6：篩選法](table_4_6_selection.csv)｜[表 4.7：特徵組合](table_4_7_combinations.csv)｜[稽核紀錄](independent_audit.json)',
            '', '## 各任務結果與圖表','']
    for task in cfg['tasks']:
        s=summary[summary.task.eq(task)]
        lines += [f'### {task}','',
                  f'[全部指標]({task}_summary.csv)｜[逐次 CV 指標]({task}_repeat_metrics.csv)｜[逐人 OOF]({task}_oof_predictions.csv)｜[分類器比較表]({task}_table_4_2_to_4_5.csv)','',
                  '| 方法 | AUC（conditional 95% CI） | AUPRC | Accuracy | Sensitivity | Specificity | Brier |',
                  '|---|---|---:|---:|---:|---:|---:|']
        for r in s.itertuples():
            lines.append(f'| {r.method} | {r.auroc:.3f} [{r.auc_low:.3f}, {r.auc_high:.3f}] | {r.auprc:.3f} | {r.accuracy:.3f} | {r.sensitivity:.3f} | {r.specificity:.3f} | {r.brier:.3f} |')
        lines += ['',f'![{task} ROC]({task}_ROC.png)','']
        for block in ['C','R','S','T','Fwhole']:
            stem=f'{task}_LightGBM_RFE_{block}'
            lines += [f'#### {block} 的 held-out SHAP','',f'![{task} {block} SHAP]({stem}_SHAP.png)','',f'[特徵重要性與入選頻率]({stem}_importance.csv)','']
            for kind in ['radiomics','availability']:
                f=report/f'{stem}_{kind}_atlas.png'
                if f.exists():
                    lines += [f'![{task} {block} {kind} 腦區圖]({f.name})','',
                              f'[腦區數值]({stem}_{kind}_regions.csv)｜[NIfTI 圖譜]({stem}_{kind}_atlas.nii.gz)','']
    lines+=['## 全部圖檔索引','', '下列亦包含 backward 與融合組合；完整數值檔與圖檔放在同一目錄。','']
    lines += [f'- [{p.name}]({p.name})' for p in sorted(report.glob('*.png'))]
    (report/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    for task in cfg['tasks']:
        for suffix in ['cohort.json','clinical_description.csv']:
            shutil.copy2(out/f'{task}_{suffix}',report/f'{task}_{suffix}')
    shutil.copy2(out/'config.json',report/'experiment_config.json')
    for name in ['NOVELTY_AND_MANUSCRIPT_DECISION.md','REPRODUCTION_SCOPE.md']:
        shutil.copy2(HERE/name,report/name)
    manifest={str(p.relative_to(report)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(report.rglob('*')) if p.is_file() and p.name not in {'MANIFEST.json','COMPLETE.json','STATUS.md'}}
    mb.write_json(report/'MANIFEST.json',manifest)
    mb.write_json(report/'COMPLETE.json',{'tasks':cfg['tasks'],'files_hashed':len(manifest),'report':'RESULTS.md'})
    print(json.dumps(audits,indent=2))
if __name__=='__main__':main()
