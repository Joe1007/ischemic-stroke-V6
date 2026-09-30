# V5 結果：AAL3 radiomics 重新抽取

主比較事前固定為 mRS 的 normal_tissue_unclipped vs clip_p99，不從 outer 結果挑 protocol。
同一 cohort、5-fold × 3 repeats、相同 V1-like LightGBM/RFE、training-only scaled SMOTE。
F 為 AAL3 whole-region radiomics；R 為 global whole-lesion radiomics；
B 為基本臨床 + 病灶體積 + texture availability + training-selected atlas burden。

所有結果均為已反覆研究的 cohort 內部估計，不是外部驗證。AUC 為病人 OOF 平均後計算；
另列三次 repeat 各自 AUC 的平均。CI 是固定 OOF 預測的 conditional patient bootstrap，未含重新訓練的不確定性。

## mRS：n=191，positive=80

| Method | Averaged-OOF AUC (95% conditional CI) | Mean repeat AUC | AUPRC | Brier |
|---|---:|---:|---:|---:|
| C | 0.650 (0.567–0.729) | 0.647 | 0.618 | 0.226 |
| B | 0.740 (0.666–0.810) | 0.729 | 0.704 | 0.200 |
| F__clip_p99 | 0.625 (0.543–0.703) | 0.616 | 0.596 | 0.243 |
| F__unclip_p99 | 0.658 (0.576–0.735) | 0.639 | 0.612 | 0.233 |
| F__normal_tissue_unclipped | 0.632 (0.549–0.713) | 0.610 | 0.609 | 0.235 |
| B_F__clip_p99 | 0.627 (0.545–0.703) | 0.616 | 0.601 | 0.245 |
| B_F__unclip_p99 | 0.665 (0.585–0.744) | 0.647 | 0.625 | 0.229 |
| B_F__normal_tissue_unclipped | 0.651 (0.574–0.731) | 0.629 | 0.606 | 0.229 |
| B_R__clip_p99 | 0.689 (0.610–0.767) | 0.677 | 0.671 | 0.217 |
| B_R__unclip_p99 | 0.758 (0.688–0.825) | 0.744 | 0.716 | 0.194 |
| B_R__normal_tissue_unclipped | 0.719 (0.642–0.791) | 0.705 | 0.677 | 0.208 |

### 預先指定的 paired comparisons

| Comparison | ΔAUC (95% conditional CI) |
|---|---:|
| F__normal_tissue_unclipped vs F__clip_p99 | +0.007 (-0.070–+0.090) |
| F__unclip_p99 vs F__clip_p99 | +0.033 (-0.005–+0.073) |
| B_F__normal_tissue_unclipped vs B_F__clip_p99 | +0.024 (-0.048–+0.103) |
| B_F__normal_tissue_unclipped vs B | -0.089 (-0.168–-0.014) |
| B_F__normal_tissue_unclipped vs B_R__normal_tissue_unclipped | -0.069 (-0.162–+0.022) |

Exact-feature mean pairwise Jaccard:
```json
{
  "clip_p99|F": 0.08235227005682556,
  "clip_p99|R": 0.28644620048732305,
  "normal_tissue_unclipped|F": 0.08091568356937662,
  "normal_tissue_unclipped|R": 0.30301500569741346,
  "unclip_p99|F": 0.10348604705084138,
  "unclip_p99|R": 0.32655183667606436
}
```

## NIHSS：n=192，positive=78

| Method | Averaged-OOF AUC (95% conditional CI) | Mean repeat AUC | AUPRC | Brier |
|---|---:|---:|---:|---:|
| C | 0.491 (0.402–0.578) | 0.495 | 0.427 | 0.264 |
| B | 0.672 (0.588–0.750) | 0.664 | 0.660 | 0.219 |
| F__clip_p99 | 0.683 (0.602–0.755) | 0.661 | 0.638 | 0.220 |
| F__unclip_p99 | 0.679 (0.595–0.753) | 0.669 | 0.649 | 0.217 |
| F__normal_tissue_unclipped | 0.694 (0.618–0.766) | 0.658 | 0.607 | 0.217 |
| B_F__clip_p99 | 0.683 (0.607–0.756) | 0.663 | 0.632 | 0.222 |
| B_F__unclip_p99 | 0.689 (0.607–0.760) | 0.673 | 0.660 | 0.215 |
| B_F__normal_tissue_unclipped | 0.685 (0.608–0.759) | 0.655 | 0.583 | 0.220 |
| B_R__clip_p99 | 0.620 (0.536–0.699) | 0.624 | 0.611 | 0.233 |
| B_R__unclip_p99 | 0.712 (0.634–0.785) | 0.698 | 0.663 | 0.208 |
| B_R__normal_tissue_unclipped | 0.720 (0.643–0.792) | 0.708 | 0.683 | 0.206 |

### 預先指定的 paired comparisons

| Comparison | ΔAUC (95% conditional CI) |
|---|---:|
| F__normal_tissue_unclipped vs F__clip_p99 | +0.012 (-0.059–+0.086) |
| F__unclip_p99 vs F__clip_p99 | -0.004 (-0.049–+0.043) |
| B_F__normal_tissue_unclipped vs B_F__clip_p99 | +0.002 (-0.070–+0.076) |
| B_F__normal_tissue_unclipped vs B | +0.013 (-0.073–+0.100) |
| B_F__normal_tissue_unclipped vs B_R__normal_tissue_unclipped | -0.036 (-0.124–+0.059) |

Exact-feature mean pairwise Jaccard:
```json
{
  "clip_p99|F": 0.07399544763123801,
  "clip_p99|R": 0.4124666780428685,
  "normal_tissue_unclipped|F": 0.08830830039852892,
  "normal_tissue_unclipped|R": 0.3707951154348413,
  "unclip_p99|F": 0.08612024659594575,
  "unclip_p99|R": 0.4225731930456322
}
```

## 範圍限制

- 本輪重抽 Original 的 93 個 intensity/texture descriptors；未抽 Wavelet/LoG，不能代表所有 radiomics protocol 的上限。
- 保留 MNI-space 3D 以隔離 normalization 因子；native-space 2D 尚未做。
- Pilot ICC 僅是 12 位病例、8 個 ROI、三方向位移的敏感度診斷，未用來篩選 modelling features。
- Clinical 未加入 admission/discharge severity，既有 NIHSS label/source discrepancies 仍未由病歷確認。
- 不與 V2/V4 的不同 cohort／模型直接作單一因素比較。

