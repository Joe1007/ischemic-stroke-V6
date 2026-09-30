#!/usr/bin/env python3
"""Compare newly trained V6 predictions with archived predictions; no cached fold reuse."""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
HERE=Path(__file__).resolve().parent;OLD=HERE.parent/'results/lesion_20260909';NEW=HERE/'output/core'
assert (NEW/'COMPLETE.json').exists()
rows=[];audit={}
for task in ['mRS','NIHSS']:
    new=pd.read_csv(NEW/f'{task}_oof_predictions.csv').set_index('ID')
    old=pd.read_csv(OLD/f'{task}_oof_predictions.csv').set_index('ID').loc[new.index]
    np.testing.assert_array_equal(new.y,old.y)
    np.testing.assert_allclose(new.drop(columns='y'),old.drop(columns='y'),atol=1e-12,rtol=0)
    detail=pd.read_csv(NEW/f'{task}_oof_by_repeat.csv')
    assert detail.groupby(['ID','method','repeat']).size().eq(1).all()
    assert detail.groupby(['ID','method']).size().eq(3).all()
    for split in range(15):
        saved=json.loads((NEW/f'{task}_split_{split:03d}.json').read_text());a=saved['audit']
        assert not set(a['train_ids'])&set(a['test_ids'])
        assert set(a['train_ids'])|set(a['test_ids'])==set(new.index)
        assert a['v5_baseline_reproduced']
        g=detail[(detail.repeat==split//5)&(detail.fold==split%5)]
        for m,p in saved['predictions'].items():
            np.testing.assert_allclose(g[g.method.eq(m)].set_index('ID').loc[a['test_ids']].probability,p,atol=1e-12,rtol=0)
    summary=pd.read_csv(NEW/f'{task}_summary.csv').set_index('method')
    for m in summary.index:
        np.testing.assert_allclose(detail[detail.method.eq(m)].groupby('ID').probability.mean().loc[new.index],new[m],atol=1e-12,rtol=0)
        assert abs(roc_auc_score(new.y,new[m])-summary.loc[m,'auroc'])<1e-12
    audit[task]={'n':len(new),'folds':15,'methods':len(summary),'reproduces_all_archived_predictions_atol':1e-12,'patient_separation_verified':True}
    for m in ['V1_base','V1_base_L','V1_base_T','V1_base_W']:
        r=summary.loc[m];rows.append({'task':task,'method':m,'n':len(new),'auc':r.auroc,'low':r.auc_low,'high':r.auc_high,'delta':r.delta_auc,'delta_low':r.delta_low,'delta_high':r.delta_high,'mean_repeat_auc':r.repeat_auroc_mean})
(NEW/'reproduction_audit.json').write_text(json.dumps(audit,indent=2))
pd.DataFrame(rows).to_csv(NEW/'core_results.csv',index=False)
lines=['# V6 全模型重新訓練核對','',
'本次使用新的結果目錄，重新訓練全部 30 folds；不是重讀舊 checkpoint。所有 24 個方法的逐人 OOF 預測均在 1e-12 絕對誤差內重現原始 V6。特徵使用經 hash 驗證的既有抽取檔，本次沒有重新抽取 MRI。','',
'歷史最高分 V1_base_T 是事後指定核心，不取代原先 primary（base_L）。同 cohort 的數值重現不是新的獨立驗證。','',
'| 任務 | 方法 | n | AUC | ΔAUC 對 base（conditional 95% CI） |','|---|---|---:|---:|---|']
for r in rows:lines.append(f"| {r['task']} | {r['method']} | {r['n']} | {r['auc']:.4f} | {r['delta']:+.4f} [{r['delta_low']:+.4f}, {r['delta_high']:+.4f}] |")
lines+=['','V6 T 的 mRS 入選幾乎都是 flags；上述效能不能被歸因於 texture。AUC 是三次 OOF probability 平均後計分，mean-repeat AUC 另存 CSV。所有 CI 為固定 OOF 的條件 bootstrap，不涵蓋訓練或研究選擇不確定性。']
(NEW/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(audit,indent=2))
