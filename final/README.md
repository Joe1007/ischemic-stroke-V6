# V6 核心與原稿結果對照重現

本目錄的新程式與結果不覆寫 2026-09-09 的歷史實驗。

- `output/core/`：原封不動呼叫 `model_v6.py`，重新訓練原 V6 全部 24 方法、兩任務、各 15 folds 的結果（2026-09-17；原本放在 `reproduction_20260917/`，2026-09-30 整併至此以免重複）。使用既有且 hash 驗證的影像特徵，不重新抽取 MRI。
- `audit_core.py`：核對新舊逐人預測、所有分割、OOF 平均與 AUC；輸出 `reproduction_audit.json` 與核心結果報告。
- `reproduce_tables.py`：四任務的論文比較結構，共 36 固定方法／任務。包含四 classifier、兩選擇法、T/Fwhole 組合、flags 對照，以及 held-out LightGBM TreeSHAP。這是採 V6／V5 特徵的**調整後重現**，不是原方法完全不變。
- `report_tables.py`：拒絕未完成／smoke 結果；完整核對 OOF、SHAP 還原預測，輸出表格、ROC、SHAP、腦區 CSV／PNG／NIfTI。
- `REPRODUCTION_SCOPE.md`：原稿頁碼、novelty 判定與交付範圍。

## 特徵名稱

C 僅保守臨床欄位，排除 mRS／NIHSS 及其衍生量；TE 另排除 AF。R 是 V5 全病灶 Original radiomics；S 為 AAL3 burden；T 是 V6 lesion∩AAL3 radiomics 加 availability；Fwhole 是 V5 完整腦區 Original unclip_p99。S 不含原稿其他 atlas，Fwhole 不含原稿全部 filters，均須修正方法描述。

原稿 notebook 稱 backward 的方法其實是 single-fit importance top-k。本程式照該操作實作，不把它改成真正逐步向後消除。

不同特徵組分別在 training fold 做 staged RFE，再串接入模型；所有 SMOTE 都在選擇後且只使用 training 病人。沒有依 outer 分數挑選模型。圖中的 SHAP 以 raw log odds 為單位，背景為訓練樹路徑（含 training-only SMOTE），不是因果解釋。

## 執行

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 MPLCONFIGDIR=/tmp/mpl-stroke-v6 PYTHONWARNINGS=ignore /home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/model_v6.py' --jobs 8 --output 'ischemic stroke analyseV6/final/output/core'
/home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/final/audit_core.py'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MPLCONFIGDIR=/tmp/mpl-stroke-v6 /home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/final/reproduce_tables.py' --jobs 2 --workers 16
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mpl-stroke-v6 /home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/final/report_tables.py'
```

`reproduce_tables.py` 只有在全部工作成功後才寫 COMPLETE；每個工作各自保留 log。重啟會檢查資料／程式 hash、環境及設定，跳過已完成方法。`--fold-limit` 僅供試跑，不能作為完整結果報告輸入。

歷史最佳分支為 V1_base_T。這是看過結果後指定的核心，不能改寫成事前 primary，亦不能聲稱重跑相同 cohort 是獨立驗證。真正接近原稿 F 的路線是 whole-AAL3；V6 需要病灶遮罩，因此不能保留「抽取無需病灶分割」的原稿主張。
