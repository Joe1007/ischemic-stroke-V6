# 與原稿核心貢獻的關係

原稿第 1-3 節的主張是腦區圖譜與 radiomics 的 Fusion Feature、相對全病灶傳統 radiomics 的改善，以及重要特徵映射回腦区。保留這些研究問題，不等於原稿的數值主張已獲得支持。

| 元素 | V6 狀態 | 論文需要修正 |
|---|---|---|
| DWI 與 AAL3 分區 | 保留 | 補充固定影像 QC、實際有效 cohort |
| 特徵保有明確腦區身分 | 保留 L/T/W | 可映射回 AAL3；不能把零重要度或不穩定選擇視作可靠生物發現 |
| 原始完整腦區 radiomics | W 沿用 V5 whole-region Original，但乘病灶占比 | 明述病灶加權、取消高強度 clipping、僅 Original 的差異 |
| 病灶局部 radiomics | L/T 為新增 | ROI 改成 lesion ∩ AAL3；補充最小 voxel 數、可用性與缺失處理 |
| 同源對側差異 | L 為新增 | 同源整區非病灶參照，非鏡像同形 ROI；雙側病灶排除範圍與未配對中線區要說明 |
| 原稿 LightGBM / RFE / SMOTE | 保留 V1 分支 | 補值、縮放、選擇與重採樣均限 training fold |
| 次要模型 | 新增固定 ET/LR | 若採用其結果必須更新 classifier 方法，不能宣稱原稿模型完全不變 |
| 相對傳統 radiomics 增量 | 必須重新驗證 | 對照包含同樣臨床／負荷資訊與 global radiomics；以 matched comparison 檢驗 |
| 原稿約 0.94 或至少 +0.1 | 不沿用 | 受 leakage 影響的舊數值不得保留為有效驗證 |
| 四項任務 | 本輪只重跑 mRS / NIHSS | 不把本輪推論延伸到 TE / LVO；其原稿主張另需驗證 |
| 外部效度 | 尚未驗證 | 反覆探索同一 cohort，最高 AUC 只能視為探索結果 |

availability 與病灶加權本身可能反映病灶位置及體積。若新模型有提升，仍需與僅用相同 availability／burden 的模型做消融，才能區分真正 radiomics 增量與更有效的位置表徵。以本輪最高分回選模型不构成未偏估的驗證。
