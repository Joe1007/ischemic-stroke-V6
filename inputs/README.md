# inputs/ — V6 的輸入資料

V6 以 V5 抽取的特徵作為強基準與對照，這些檔案原本放在已刪除的 `ischemic stroke analyseV5/`，
2026-09-30 原樣搬入（內容與雜湊值不變）。

| 檔案 | 原位置 | 用途 |
|---|---|---|
| `v5_features.npz` | V5/output/features.npz | 臨床以外的基準特徵：burden、全病灶 R、完整腦區 Fwhole、volume（model_v6、reproduce_tables、atlas_comparison） |
| `v5_checkpoints/` | V5/output/checkpoints/ | 每人 V5 抽取結果；extract_v6 核對影像來源、assemble_v6 與 validate_aal3 比對 |
| `v5_protocol.json` | V5/config/protocol.json | 影像正規化與 radiomics 參數（extract_v6、extract_atlas） |
| `v5_reference_predictions/` | V5/results/normalization_20260909/ | model_v6 逐 fold 斷言基準預測與 V5 完全一致；report_v6 對照 |

`v5_features.npz` SHA256：`043da2b1cf1c57d3e45d832559333418d2e452595fc7c4b19711a88585d14417`（= `output/features_manifest.json` 的 `v5_archive_sha256`）

## 程式碼

原本 import 自 V2–V5 的函式已併入 V6 根目錄（邏輯未改，只保留 V6 實際呼叫的定義）：

- `modeling_base.py`：V2 leakage_free_experiment / atlas_group_experiment、V4 run_v4、V5 model_v5
- `extraction_base.py`：V3 extract_radiomics_v3、V5 extract_v5

因程式檔位置與內容改變，舊結果資料夾內 `config.json`、`*.audit.json` 記錄的來源路徑與雜湊值
不再對得上；重跑請指定新的 `--output`，不要續寫舊結果資料夾。

## 驗證（V2–V5、upstream、source_snapshot 皆已刪除的狀態）

單元測試、全部腳本 import、`assemble_atlas.verify_aal3()`、`validate_aal3.py`（198 人／32868 區域完全一致）、
mRS split 0 模型重算（24 方法與原結果最大差 2.2e-16）、以及 1 位病人從原始 MRI 重抽（166 區域全部相同）均通過。
