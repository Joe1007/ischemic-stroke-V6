# V6 結果：病灶相關 AAL3 radiomics

本輪完整執行 198 人特徵重抽、兩結局各 5 × 3 repeated CV。所有預先固定組合均列出。這是已反覆探索 cohort 的結果，並非外部或未碰觸測試集的效能。

## 特徵與核心 novelty

- L：逐區病灶一階統計與同源對側非病灶組織差異。
- T：逐區病灶與 AAL3 交集的 93 個 Original radiomics。
- W：完整 AAL3 腦區 radiomics 乘以該區病灶占比，最接近原稿 ROI。
- base：臨床 core + 體積 + AAL3 burden + 全病灶 radiomics；V1_base 完全重現 V5 B_R__unclip_p99。
- V1：原稿 LightGBM / staged RFE / SMOTE；ET：固定 ExtraTrees；LR：固定 L2 logistic。ET/LR 各自的 base 使用全部 training-preprocessed burden/global 特徵，必須以同模型比較增量。
- blend：固定 75% V1_base + 25% 區域單模態預測，未依 outer 表現選權重。

可用區域中位數：受影響 8、一階特徵 5、對側對照 5、紋理 1。
沒有任何可用局部紋理的個案：77/198；沒有局部一階特徵：24/198。

## 全部配對結果

AUC 是每人三次 OOF probability 平均後計分；mean-repeat AUC 為三個 repeat 分別計分後平均，兩者不可混稱。

| 方法 | mRS AUC | mRS mean-repeat | NIHSS AUC | NIHSS mean-repeat |
|---|---:|---:|---:|---:|
| ET_L | 0.559 | 0.551 | 0.630 | 0.620 |
| ET_T | 0.551 | 0.554 | 0.609 | 0.605 |
| ET_W | 0.561 | 0.555 | 0.607 | 0.603 |
| ET_base | 0.709 | 0.698 | 0.714 | 0.705 |
| ET_base_L | 0.693 | 0.682 | 0.706 | 0.696 |
| ET_base_T | 0.675 | 0.665 | 0.711 | 0.703 |
| ET_base_W | 0.691 | 0.680 | 0.695 | 0.684 |
| LR_L | 0.504 | 0.500 | 0.503 | 0.505 |
| LR_T | 0.564 | 0.564 | 0.492 | 0.496 |
| LR_W | 0.533 | 0.523 | 0.481 | 0.478 |
| LR_base | 0.667 | 0.655 | 0.667 | 0.661 |
| LR_base_L | 0.564 | 0.561 | 0.536 | 0.537 |
| LR_base_T | 0.566 | 0.568 | 0.545 | 0.549 |
| LR_base_W | 0.556 | 0.544 | 0.501 | 0.497 |
| V1_L | 0.619 | 0.607 | 0.603 | 0.600 |
| V1_T | 0.616 | 0.605 | 0.625 | 0.614 |
| V1_W | 0.639 | 0.621 | 0.611 | 0.605 |
| V1_base | 0.758 | 0.744 | 0.712 | 0.698 |
| V1_base_L | 0.755 | 0.742 | 0.707 | 0.692 |
| V1_base_T | 0.767 | 0.752 | 0.718 | 0.705 |
| V1_base_W | 0.753 | 0.736 | 0.702 | 0.687 |
| V1_blend_L | 0.758 | 0.743 | 0.701 | 0.693 |
| V1_blend_T | 0.757 | 0.741 | 0.703 | 0.697 |
| V1_blend_W | 0.761 | 0.744 | 0.707 | 0.697 |

## 預先指定主要比較

mRS V1_base_L 相對 V1_base：ΔAUC -0.003，conditional 95% CI [-0.019, +0.013]。

## 各模型融合特徵的配對增量

| 結局 | 方法 | 相同模型基準 | ΔAUC | conditional 95% CI |
|---|---|---|---:|---|
| mRS | V1_base_L | V1_base | -0.003 | [-0.019, +0.013] |
| mRS | V1_base_T | V1_base | +0.009 | [-0.005, +0.023] |
| mRS | V1_base_W | V1_base | -0.005 | [-0.021, +0.011] |
| mRS | V1_blend_L | V1_base | +0.000 | [-0.015, +0.015] |
| mRS | V1_blend_T | V1_base | -0.001 | [-0.016, +0.014] |
| mRS | V1_blend_W | V1_base | +0.002 | [-0.013, +0.017] |
| mRS | ET_base_L | ET_base | -0.016 | [-0.042, +0.008] |
| mRS | ET_base_T | ET_base | -0.035 | [-0.058, -0.013] |
| mRS | ET_base_W | ET_base | -0.018 | [-0.048, +0.011] |
| mRS | LR_base_L | LR_base | -0.103 | [-0.172, -0.039] |
| mRS | LR_base_T | LR_base | -0.101 | [-0.188, -0.018] |
| mRS | LR_base_W | LR_base | -0.111 | [-0.189, -0.034] |
| NIHSS | V1_base_L | V1_base | -0.005 | [-0.025, +0.014] |
| NIHSS | V1_base_T | V1_base | +0.006 | [-0.015, +0.025] |
| NIHSS | V1_base_W | V1_base | -0.011 | [-0.032, +0.010] |
| NIHSS | V1_blend_L | V1_base | -0.011 | [-0.029, +0.007] |
| NIHSS | V1_blend_T | V1_base | -0.009 | [-0.027, +0.009] |
| NIHSS | V1_blend_W | V1_base | -0.006 | [-0.023, +0.013] |
| NIHSS | ET_base_L | ET_base | -0.008 | [-0.034, +0.021] |
| NIHSS | ET_base_T | ET_base | -0.003 | [-0.027, +0.024] |
| NIHSS | ET_base_W | ET_base | -0.020 | [-0.051, +0.015] |
| NIHSS | LR_base_L | LR_base | -0.131 | [-0.194, -0.072] |
| NIHSS | LR_base_T | LR_base | -0.122 | [-0.193, -0.049] |
| NIHSS | LR_base_W | LR_base | -0.166 | [-0.230, -0.101] |

## 入選特徵的種類與事後消融

- mRS L：300 個跨 fold 入選名額中，159 個是病灶／availability 指標，141 個是 radiomics 衍生特徵。
- mRS T：300 個跨 fold 入選名額中，299 個是病灶／availability 指標，1 個是 radiomics 衍生特徵。
- mRS W：300 個跨 fold 入選名額中，29 個是病灶／availability 指標，271 個是 radiomics 衍生特徵。
- NIHSS L：525 個跨 fold 入選名額中，215 個是病灶／availability 指標，310 個是 radiomics 衍生特徵。
- NIHSS T：525 個跨 fold 入選名額中，449 個是病灶／availability 指標，76 個是 radiomics 衍生特徵。
- NIHSS W：525 個跨 fold 入選名額中，23 個是病灶／availability 指標，502 個是 radiomics 衍生特徵。

事後機制消融（並非新驗證）：mRS 相同 V1 baseline 加上僅有 flags 的 RFE 分支 AUC 0.761；T 分支 0.767，相對 flags-only 的 ΔAUC +0.005，conditional 95% CI [-0.004, +0.014]。
此消融重新在每個 training fold 的 flags 集合上做 RFE，與 T 的入選 flags 集合可能不同，因此不是逐一刪除固定紋理欄位的精確因果貢獻測量。T 的高分幾乎由非 radiomics 指標構成，不能據此宣稱核心 radiomics novelty 已驗證。

## 判讀限制

局部特徵保留 AAL3 區域身分，但 ROI 定義、強度處理與融合已不同於原稿，必須修正方法段。availability flags 與病灶占比本身攜帶病灶位置／大小訊息；高分不能單獨證明 texture 或對側差異有增量，需消融比較。
每個 training fold 內才補值、縮放、RFE、SMOTE；ET/LR 不做跨 fold 特徵篩選。所有模型為固定比較，未使用 outer labels 調參。跨 V2–V6 的反覆研究仍使最高分具有選擇樂觀性。
95% CI 是固定 OOF 預測的分層病人 bootstrap（2000 次），不涵蓋整個重訓與研究選擇，也未作多重比較校正。
小病灶、MNI 重採樣與 native 厚切面仍限制紋理穩健性；本輪沒有 native 2D 重建，也沒有改動原始 outcome。

## 驗證

4 項抽取測試與 2 項 training-only 前處理測試通過。獨立核對 30 個 folds 的 train/test ID、每人三次 OOF、全部 AUC／平均／基準增量；V5 強基準在每個 fold 完全重現。詳見 independent_audit.json。
