# V6：病灶相關的 AAL3 radiomics

> **資料說明**：本 repo 只包含程式碼、設定、文件與彙總結果（指標表、圖、腦區圖譜）。
> 病人層級資料（影像、臨床表、逐人特徵、逐人預測與 fold 名單）依研究倫理規範保留在實驗室主機，未公開；
> 排除規則見 `.gitignore`。套件版本見 `requirements.txt`。完整方法與結果見 `final/output/METHODS_V6_DETAILED.md`。

保留 DWI、AAL3 分區身分、可解釋的 radiomics，以及原稿的 LightGBM / RFE / SMOTE 分支。重新抽取病灶與腦區交集，處理全腦區特徵稀釋小病灶訊號的假設。

固定三組特徵：

- **L**：逐區病灶 8 個一階統計與同源對側差異，另附可計算性指標。對照為整個同源 AAL3 區域的非病灶組織，排除病灶及其 2 voxel 膨脹範圍；不是鏡像病灶形狀的 ROI。
- **T**：逐區病灶 93 個 Original radiomics，需超過 64 個有效 voxel 且至少 2 維有空間範圍。
- **W**：V5 完整腦區的 93 個 Original radiomics 乘該區病灶占比，保留原稿全腦區抽取方式，但特徵意義變成病灶加權版本。

小病灶至少 8 voxel 可取得 L。沒有病灶設零；有病灶但特徵不可計算设 NaN，加入 availability 指標，僅於 training fold 補值。由於影像已重採樣到 2 mm MNI，voxel 數不是獨立的 native 採樣點數，尤其不能把小 ROI 紋理視為已證實穩健。

正規化固定為 V5 `unclip_p99`。不改分組、標籤或排除個案，不使用結局挑 ROI。198 人影像 cohort、mRS / NIHSS 的 5 × 3 repeated stratified CV 與 V5 相同。

主要比較預先固定為 mRS 原版模型 `V1_base_L` 相對 `V1_base`。後者是 V5 的 `B_R__unclip_p99`；每個 fold 都驗證其預測完全重現。這個基準包含保守臨床變數、體積、AAL3 burden、全病灶 radiomics。不能只與弱 clinical-only 基準比較。

另外固定比較 L/T/W 單獨、直接融合，以及原版模型的 75% baseline + 25% regional probability blend。ExtraTrees 與 L2 logistic 為次要模型，使用各自對應的相同模型 baseline，training-only 處理全部特徵；不在 outer folds 選模型或調權重。所有組合完整報告，不能把探索中的最高分當成獨立驗證。

執行：

```bash
OPENBLAS_NUM_THREADS=1 ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 /home/joe1007/miniconda3/envs/radiomics311/bin/python 'ischemic stroke analyseV6/extract_v6.py' --jobs 6
/home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/assemble_v6.py'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 MPLCONFIGDIR=/tmp/mpl-stroke-v6 PYTHONWARNINGS=ignore /home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/model_v6.py' --jobs 8
/home/joe1007/miniconda3/envs/AD/bin/python 'ischemic stroke analyseV6/report_v6.py'
```

`config/protocol.json` 在本輪 outcome 模型執行前寫入。研究已反覆使用此 cohort，這不是全新、未碰觸的驗證集；conditional bootstrap 不涵蓋重訓與跨版本研究選擇不確定性。
