# EGFR-PACC 初步数字细胞构建报告

生成时间：2026-06-18T02:58:39+00:00

## 构建策略

本轮已将 COSMIC、GENIE、DrugBank 从必需依赖改为 optional enhanced dependency。初步数字细胞使用开放替代源构建：cBioPortal、TCGA/GDC、CIViC、ClinVar、Cancer Hotspots、ChEMBL、PubChem、DailyMed、openFDA 和 manual literature seed。

## 输出文件

- 机器可读模型快照：`configs/open_fallback_digital_cell.yaml`
- mutation catalog：`data/processed/reference/egfr_open_mutation_catalog.csv`
- enhanced mutation catalog：`data/processed/reference/egfr_enhanced_mutation_catalog.csv`（如 COSMIC 已下载并解析）
- frequency prior：`data/processed/reference/open_egfr_pacc_frequency.csv`
- CNS modifier：`data/processed/reference/egfr_tki_cns_modifier_table.csv`
- source priority：`data/processed/reference/source_priority_map.csv`

## 完成度检查

| 项目 | 状态 | 行数 | 路径 |
| --- | --- | ---: | --- |
| egfr_open_mutation_catalog | complete | 35 | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/processed/reference/egfr_open_mutation_catalog.csv` |
| egfr_enhanced_mutation_catalog | complete | 53 | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/processed/reference/egfr_enhanced_mutation_catalog.csv` |
| open_egfr_pacc_frequency | complete | 24 | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/processed/reference/open_egfr_pacc_frequency.csv` |
| egfr_tki_cns_modifier_table | complete | 11 | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/processed/reference/egfr_tki_cns_modifier_table.csv` |
| source_priority_map | complete | 26 | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/processed/reference/source_priority_map.csv` |
| optional_enhanced_cosmic | available |  | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/raw/cosmic` |
| optional_enhanced_genie | warning_missing_optional |  | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/raw/genie` |
| optional_enhanced_drugbank | warning_missing_optional |  | `/Volumes/deaneu01/egfr_pacc_digital_cell/data/raw/drugbank` |
| model_config | complete |  | `/Volumes/deaneu01/egfr_pacc_digital_cell/configs/open_fallback_digital_cell.yaml` |

## 当前可用模型组件

- EGFR 重点突变目录：8 个 mutation bucket。
- EGFR 增强突变目录：8 个 mutation bucket。
- 开放队列频率：6 个 cohort。
- EGFR 药物 CNS modifier：11 个药物。

## 关键限制

- CNS modifier 目前包含 manual literature seed，尚未完成正式文献锁表。
- OncoKB、COSMIC、GENIE、DrugBank 均可增强模型，但缺失不会阻断 MVP。
- 当前模型可用于候选排序与流程验证，不能视为临床验证模型。

## 建议下一步

1. 锁定 Robichaux、LUX-Lung、ACHILLES/TORG1834、UNICORN 的 mutation x drug ORR/PFS 表。
2. 用锁表更新 `configs/priors.yaml`，保留 `source_priority_map.csv` 的字段级溯源。
3. 将 37 例临床数据先导入 feature-only 表，确认不含 PFS/ORR/OS。
4. 冻结 priors/baselines 后再打开 outcome，做 EPDRI vs Baseline-1 的探索性验证。
5. 后续有 OncoKB/COSMIC/GENIE/DrugBank 时作为 enhanced layer 增量比较，不回改已冻结 open baseline。
