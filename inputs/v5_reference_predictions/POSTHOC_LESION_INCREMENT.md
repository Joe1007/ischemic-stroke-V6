# Whole-lesion radiomics 的增量（事後描述性分析）

看到 mRS 的 B+R(unclip_p99) AUC 較高後，另外計算全部三個 R protocol 相對 B 的 paired CI。
没有重訓、挑特徵或改模型；這不是預先指定的主比較，也不是獨立確認。
B 已含臨床、病灶體積與 atlas lesion burden；R 是整個病灶的 radiomics，不是 AAL3 區域紋理。

| Task | Method | AUC | ΔAUC vs B (conditional 95% CI) |
|---|---|---:|---:|
| mRS | B_R__clip_p99 | 0.689 | -0.051 (-0.086–-0.015) |
| mRS | B_R__unclip_p99 | 0.758 | +0.018 (-0.027–+0.063) |
| mRS | B_R__normal_tissue_unclipped | 0.719 | -0.020 (-0.071–+0.034) |
| NIHSS | B_R__clip_p99 | 0.620 | -0.052 (-0.111–+0.003) |
| NIHSS | B_R__unclip_p99 | 0.712 | +0.040 (-0.035–+0.108) |
| NIHSS | B_R__normal_tissue_unclipped | 0.720 | +0.048 (-0.026–+0.121) |
