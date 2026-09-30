# V6 實驗紀錄（摘自專案總誌 EXPERIMENT_HISTORY.md 第 23 節）

> 原總誌涵蓋 V2–V8 迭代，已於 2026-09-30 刪除；此處只保留 V6 相關段落，內容未修改。

## 23. V6：病灶相關 AAL3 radiomics 重抽（2026-09-09）

使用者授權再次重抽特徵、保留原稿 atlas/radiomics 核心。建立 `ischemic stroke analyseV6/`，在讀取本輪模型效能前固定 protocol。

- L：lesion ∩ AAL3 的 8 個一階統計與同源對側整區非病灶參照差異；至少 8 voxel。
- T：lesion ∩ AAL3 的 93 個 Original radiomics；>64 voxel、至少 2 維有空間範圍。
- W：V5 unclip_p99 whole-AAL3 radiomics 乘對應腦區病灶占比，保留原稿完整腦區 ROI。
- 正規化固定 unclip_p99。無病灶零、有病灶但不可計算 NaN，另有 availability flags。對側依 `_L`/`_R` 解剖名稱配對，參照排除膨脹 2 voxel 的病灶，不是假定奇偶標籤，也不是鏡像病灶 ROI。
- 198/198 完成抽取；可用一階區域中位數 5、紋理區域中位數 1。828 個可算紋理 ROI 的自算一階統計與 PyRadiomics 對應結果一致（max absolute difference 1.96e-14）。4 項抽取、2 項 training-only preprocessing 測試通過，全部矩陣的 absence / missing 編碼核對通過。
- mRS/NIHSS 仍為 V5 同 cohort、相同 5×3 splits。V1 baseline 每 fold 必須完全重現 V5 B_R__unclip_p99。所有固定模型組合皆報告，主要比較為 mRS V1_base_L vs V1_base；另有固定 ET/LR 與 75% baseline +25% regional 的 V1 blend。
- 研究 cohort 已被反覆探索，結果不能視為外部驗證；高分也需要區分 availability／病灶位置與 radiomics 的貢獻。

本節建立時模型正在執行，完成狀態以 `ischemic stroke analyseV6/results/lesion_20260909/COMPLETE.json` 與 `RESULTS.md` 為準，結果完成後續補。

### 23.1 mRS 完成與事後機制核對（2026-09-10）

mRS 191 人的 15 folds 已完成，V1 強基準逐 fold 完全重現 V5。各固定方法與配對 CI 已保存 `results/lesion_20260909/mRS_summary.csv`、`mRS_matched_comparisons.csv`。

- V1_base AUC 0.758108；base+L 0.755293；base+T 0.766779；base+W 0.753378。
- 預先指定主比較 base+L vs base：Δ −0.002815，conditional 95% CI [−0.019147, +0.012838]。
- 本輪最高 base+T vs base：Δ +0.008671，CI [−0.005293, +0.023201]。三個 repeat AUC 分別 0.757376 / 0.727646 / 0.770608，mean-repeat 0.751877；0.766779 是平均 OOF probability 後的 AUC。
- 重要限制：T 跨 15 folds 的 300 個入選名額中，299 個是 affected / availability flags，只有 1 個是真正 radiomics。不能把這組高分當成 texture novelty 成功。
- 因看到上述選擇結果，新增明確標註 post-hoc 的 flags-only 機制消融（`ablate_flags.py`，`results/flags_ablation_20260910/`）。相同 V1 baseline + training-RFE flags 得 AUC 0.761374。T 相對 flags-only Δ +0.005405，CI [−0.003716, +0.014417]，未證明紋理增量。此消融重新選 flags，不是精確刪除固定紋理欄位的因果比較。
- ET/LR 沒有提高 mRS，完整列出所有 24 個固定方法，未只挑最高組報告。

續跑時發現先前程序未持續執行，停在 mRS folds 完成後、bootstrap 彙整前；已由原 checkpoints 恢復，不改 protocol / 參數 / splits。NIHSS 正在完成中。

### 23.2 V6 全部完成與獨立核對（2026-09-10）

mRS 191 人、NIHSS 192 人，各 15 folds × 24 個固定方法皆完成；`results/lesion_20260909/COMPLETE.json`、`independent_audit.json` 與 `ischemic stroke analyseV6/RESULTS.md` 已產生。此完成記錄取代本節先前「正在執行」狀態。

NIHSS 最高 V1_base_T 0.718061，相同 unclip_p99 基準 0.712438，Δ +0.005623，conditional 95% CI [−0.014735, +0.024972]；亦未超過 V5 另一 normalization 的既有 0.720423。T 跨 fold 525 個入選名額中，449 個 flags、76 個 radiomics。

結論：本輪重新抽取病灶局部與對側特徵、病灶加權 whole-AAL3、替換固定分類器，未取得明確穩定突破，更未達 0.85。mRS 探索最高 0.766779 的增益小且 CI 跨零，主要選入非 radiomics 指標，不可用來維持原稿「atlas radiomics 明顯優於傳統全病灶特徵」的數值主張。這也不是 0.75 的理論上限證明。

全部 30 folds 的 ID 隔離、每人三次 OOF、完整 AUC／平均／基準 delta 均由 report_v6.py 獨立重算通過，V5 baseline 逐 fold 完全重現。完整方法、程式、protocol、source snapshot、6 項測試、抽取 QC 與事後消融均保存在 V6。

