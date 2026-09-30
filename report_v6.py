#!/usr/bin/env python3
"""Independent OOF integrity checks and complete non-cherry-picked report."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;INPUTS=HERE/'inputs';OUT=HERE/'results/lesion_20260909'

def main():
    assert (OUT/'COMPLETE.json').exists()
    fm=json.loads((HERE/'output/features_manifest.json').read_text());qc=pd.DataFrame(fm['qc'])
    qc.to_csv(HERE/'output/availability_qc.csv',index=False)
    tables={};matched={};checks={}
    for task in ['mRS','NIHSS']:
        oof=pd.read_csv(OUT/f'{task}_oof_predictions.csv')
        detail=pd.read_csv(OUT/f'{task}_oof_by_repeat.csv')
        summary=pd.read_csv(OUT/f'{task}_summary.csv').set_index('method')
        pairs=pd.read_csv(OUT/f'{task}_matched_comparisons.csv').set_index('method')
        ids=set(oof.ID);assert len(ids)==len(oof)
        assert detail.groupby(['ID','method']).size().eq(3).all()
        assert detail.probability.between(0,1).all()
        assert detail.groupby('ID').y.nunique().eq(1).all()
        old=pd.read_csv(INPUTS/'v5_reference_predictions'/f'{task}_oof_predictions.csv').set_index('ID').loc[oof.ID]
        np.testing.assert_allclose(oof.V1_base,old.B_R__unclip_p99,atol=1e-12,rtol=0)
        for split in range(15):
            saved=json.loads((OUT/f'{task}_split_{split:03d}.json').read_text())
            a=saved['audit'];tr=set(a['train_ids']);te=set(a['test_ids'])
            assert not tr&te and tr|te==ids
            assert a['v5_baseline_reproduced']
            g=detail[(detail.repeat==split//5)&(detail.fold==split%5)]
            assert set(g.ID)==te
            for method,probs in saved['predictions'].items():
                np.testing.assert_allclose(g[g.method==method].set_index('ID').loc[a['test_ids']].probability,probs,atol=1e-12,rtol=0)
        for method in summary.index:
            avg=detail[detail.method==method].groupby('ID').probability.mean().loc[oof.ID]
            np.testing.assert_allclose(avg,oof[method],atol=1e-12,rtol=0)
            auc=roc_auc_score(oof.y,oof[method]);assert abs(auc-summary.loc[method,'auroc'])<1e-12
            assert abs(auc-roc_auc_score(oof.y,oof.V1_base)-summary.loc[method,'delta_auc'])<1e-12
        tables[task]=summary;matched[task]=pairs
        checks[task]={'patients':len(ids),'outer_folds':15,'methods':len(summary),
                      'oof_per_patient':3,'disjoint_ids':True,'v5_baseline_reproduced':True,
                      'all_averages_aucs_and_reference_deltas_recomputed':True}
    (OUT/'independent_audit.json').write_text(json.dumps(checks,indent=2))
    lines=['# V6 結果：病灶相關 AAL3 radiomics','',
        '本輪完整執行 198 人特徵重抽、兩結局各 5 × 3 repeated CV。所有預先固定組合均列出。這是已反覆探索 cohort 的結果，並非外部或未碰觸測試集的效能。','',
        '## 特徵與核心 novelty','',
        '- L：逐區病灶一階統計與同源對側非病灶組織差異。',
        '- T：逐區病灶與 AAL3 交集的 93 個 Original radiomics。',
        '- W：完整 AAL3 腦區 radiomics 乘以該區病灶占比，最接近原稿 ROI。',
        '- base：臨床 core + 體積 + AAL3 burden + 全病灶 radiomics；V1_base 完全重現 V5 B_R__unclip_p99。',
        '- V1：原稿 LightGBM / staged RFE / SMOTE；ET：固定 ExtraTrees；LR：固定 L2 logistic。ET/LR 各自的 base 使用全部 training-preprocessed burden/global 特徵，必須以同模型比較增量。',
        '- blend：固定 75% V1_base + 25% 區域單模態預測，未依 outer 表現選權重。','',
        f"可用區域中位數：受影響 {qc.affected.median():.0f}、一階特徵 {qc.firstorder_available.median():.0f}、對側對照 {qc.contrast_available.median():.0f}、紋理 {qc.texture_available.median():.0f}。",
        f"沒有任何可用局部紋理的個案：{int((qc.texture_available==0).sum())}/198；沒有局部一階特徵：{int((qc.firstorder_available==0).sum())}/198。",'',
        '## 全部配對結果','',
        'AUC 是每人三次 OOF probability 平均後計分；mean-repeat AUC 為三個 repeat 分別計分後平均，兩者不可混稱。','',
        '| 方法 | mRS AUC | mRS mean-repeat | NIHSS AUC | NIHSS mean-repeat |',
        '|---|---:|---:|---:|---:|']
    for method in tables['mRS'].index:
        a=tables['mRS'].loc[method];b=tables['NIHSS'].loc[method]
        lines.append(f'| {method} | {a.auroc:.3f} | {a.repeat_auroc_mean:.3f} | {b.auroc:.3f} | {b.repeat_auroc_mean:.3f} |')
    lines+=['','## 預先指定主要比較','']
    a=matched['mRS'].loc['V1_base_L']
    lines.append(f"mRS V1_base_L 相對 V1_base：ΔAUC {a.delta_auc:+.3f}，conditional 95% CI [{a.delta_low:+.3f}, {a.delta_high:+.3f}]。")
    lines+=['','## 各模型融合特徵的配對增量','',
            '| 結局 | 方法 | 相同模型基準 | ΔAUC | conditional 95% CI |','|---|---|---|---:|---|']
    for task,pairs in matched.items():
        for method,a in pairs.iterrows():
            if method==a.reference:continue
            lines.append(f'| {task} | {method} | {a.reference} | {a.delta_auc:+.3f} | [{a.delta_low:+.3f}, {a.delta_high:+.3f}] |')
    lines+=['','## 入選特徵的種類與事後消融','']
    flag_names={'affected','firstorder_available','contrast_available','texture_available'}
    for task in ['mRS','NIHSS']:
        for block in ['L','T','W']:
            selected=[f for split in range(15) for f in json.loads((OUT/f'{task}_split_{split:03d}.json').read_text())['audit']['blocks'][block]['selected_features']]
            nflags=sum(f.rsplit('|',1)[-1] in flag_names for f in selected)
            lines.append(f'- {task} {block}：{len(selected)} 個跨 fold 入選名額中，{nflags} 個是病灶／availability 指標，{len(selected)-nflags} 個是 radiomics 衍生特徵。')
    ap=HERE/'results/flags_ablation_20260910/mRS_summary.csv'
    if ap.exists():
        af=pd.read_csv(ap).set_index('method');t=af.loc['T'];f=af.loc['flags']
        lines+=['',f"事後機制消融（並非新驗證）：mRS 相同 V1 baseline 加上僅有 flags 的 RFE 分支 AUC {f.auroc:.3f}；T 分支 {t.auroc:.3f}，相對 flags-only 的 ΔAUC {t.delta_auc:+.3f}，conditional 95% CI [{t.delta_low:+.3f}, {t.delta_high:+.3f}]。",
            '此消融重新在每個 training fold 的 flags 集合上做 RFE，與 T 的入選 flags 集合可能不同，因此不是逐一刪除固定紋理欄位的精確因果貢獻測量。T 的高分幾乎由非 radiomics 指標構成，不能據此宣稱核心 radiomics novelty 已驗證。']
    lines+=['','## 判讀限制','',
        '局部特徵保留 AAL3 區域身分，但 ROI 定義、強度處理與融合已不同於原稿，必須修正方法段。availability flags 與病灶占比本身攜帶病灶位置／大小訊息；高分不能單獨證明 texture 或對側差異有增量，需消融比較。',
        '每個 training fold 內才補值、縮放、RFE、SMOTE；ET/LR 不做跨 fold 特徵篩選。所有模型為固定比較，未使用 outer labels 調參。跨 V2–V6 的反覆研究仍使最高分具有選擇樂觀性。',
        '95% CI 是固定 OOF 預測的分層病人 bootstrap（2000 次），不涵蓋整個重訓與研究選擇，也未作多重比較校正。',
        '小病灶、MNI 重採樣與 native 厚切面仍限制紋理穩健性；本輪沒有 native 2D 重建，也沒有改動原始 outcome。','',
        '## 驗證','',
        '4 項抽取測試與 2 項 training-only 前處理測試通過。獨立核對 30 個 folds 的 train/test ID、每人三次 OOF、全部 AUC／平均／基準增量；V5 強基準在每個 fold 完全重現。詳見 independent_audit.json。']
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(checks,indent=2));print('Wrote',HERE/'RESULTS.md')
if __name__=='__main__':main()
