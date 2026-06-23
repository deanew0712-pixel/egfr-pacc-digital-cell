# EGFR-PACC Digital Cell 项目完成记录

生成日期：2026-06-20 Asia/Shanghai

## 当前结论

本项目已经完成从公开知识库构建、文献先验冻结、贝叶斯部分池化、EPDRI 排名、文献 leave-one-out 回测，到本地 31 例临床治疗数据探索性验证的闭环。

当前最稳妥的论文定位是：

> prior-locked, interpretable EGFR uncommon/PACC digital-cell framework with literature back-testing and exploratory local clinical validation.

不建议写成：

> clinically validated treatment decision system.

## 已完成内容

1. 公开数据与可选增强数据层
   - CIViC、cBioPortal、GDC/TCGA、ClinVar、Cancer Hotspots、ChEMBL、PubChem、openFDA/FDA label。
   - COSMIC v104 和 OncoKB 已作为 optional enhanced sources 接入。
   - 亚洲/中国 uncommon EGFR/PACC 数据已整理。

2. 模型核心
   - Robichaux-style 四分类。
   - 复合突变规则。
   - frozen baseline。
   - PyMC 4-chain posterior prior table。
   - EPDRI 三层评分：locked prior、clinical modifier、uncertainty/evidence cap。

3. 文献回测
   - locked literature table 共 34 行。
   - ORR-evaluable folds 共 29 行。
   - model AUC 0.8787。
   - Baseline-0 AUC 0.8033。
   - Baseline-1 AUC 0.8602。
   - 相对 Robichaux-only baseline 的 AUC 改善有 bootstrap 支持。

4. 本地临床验证
   - 来源文件：`C:\Users\Administrator\Desktop\SupplementaryTables_SUBMISSION_FINAL.xlsx`。
   - 总病例数：31。
   - 主分析：24 例 model-native first-line EGFR-TKI monotherapy。
   - 敏感性分析 1：25 例，加入 BCP + afatinib 组合暴露。
   - 敏感性分析 2：29 例，加入 icotinib 作为一代 EGFR-TKI proxy。
   - 2 例 platinum chemotherapy 被标记为 unsupported treatment，不强行纳入 EGFR-TKI ranking 验证。

5. 基础结构 docking 补充分析
   - 使用 4LRM chain A 作为 EGFR kinase receptor template。
   - 使用 YUN 共晶配体定义 ATP pocket docking box。
   - 生成 WT、G719S、L861Q、S768I receptor models。
   - 对 afatinib、osimertinib、furmonertinib 做 AutoDock Vina 1.2.7 粗略 docking。
   - 生成 docking score matrix、pose PDBQT、heatmap 和 representative pose overview。
   - 已解决 Meeko receptor preparation 在当前 Windows/Python 运行时的 RDKit 残基价态解析失败问题：后续 receptor PDBQT 统一使用 OpenBabel 生成，Meeko 仅用于 ligand PDBQT。
   - 该分析明确标注为 coarse non-covalent docking，不作为 FEP/MD 或实验亲和力解释。
   - 已补充 receptor preparation/protonation workflow sensitivity：OpenBabel default、OpenBabel Gasteiger、PDBFixer pH 7.4 + OpenBabel Gasteiger 三条 workflow 下，关键 PACC receptor-drug 排名方向整体稳定。
   - 已补充 acrylamide warhead-Cys pre-covalent geometry screen，用于评估 afatinib/osimertinib/furmonertinib pose 中 warhead 对 ATP-pocket cysteine 的可及性；该部分作为共价 pose plausibility，不作为严格 covalent binding free energy。

## 本地验证主要结果

主分析 model-native monotherapy set（n = 24）：

- Spearman score vs PFS1 rho = 0.288，bootstrap 95% CI -0.184 至 0.687。
- Pairwise C-index = 0.595，bootstrap 95% CI 0.426 至 0.752。
- Top-1 coverage = 29.2%。
- Top-2 coverage = 50.0%。
- Top-3 coverage = 54.2%。
- 高于中位 EPDRI score 组的 median PFS1 比低分组高 8.75 个月。
- KM / log-rank：高 EPDRI 组 median PFS1 16.95 月，低 EPDRI 组 8.20 月，log-rank p = 0.207。

Proxy sensitivity set including icotinib（n = 29）：

- Spearman rho = 0.396，bootstrap 95% CI -0.001 至 0.711。
- Pairwise C-index = 0.633，bootstrap 95% CI 0.497 至 0.757。
- Top-1 coverage = 27.6%。
- Top-2 coverage = 44.8%。
- Top-3 coverage = 48.3%。
- KM / log-rank：高 EPDRI 组 median PFS1 16.80 月，低 EPDRI 组 7.35 月，log-rank p = 0.060。

n=24 与 n=29 正式一致性比较：

- Spearman 和 C-index 在两组中方向一致，均为正向。
- 两组 bootstrap CI 重叠。
- 因此可以写作“敏感性分析支持主分析方向稳健”，但仍应避免写成确定性验证。

亚组结果：

- PACC subgroup 在 n=29 proxy sensitivity set 中 rho = 0.394，C-index = 0.638。
- Compound mutation subgroup 在 n=29 proxy sensitivity set 中 rho = 0.524，C-index = 0.673。
- classical-like subgroup 仅 n=8，主分析 compound subgroup 仅 n=9，应标注为 descriptive only。
- PACC subgroup 和 n=29 compound subgroup 可写作 hypothesis-generating subgroup signal，不能写成独立确证性结论。
- 这与文献 LOO 中 PACC-only subset 表现更强的发现方向一致，是论文中最值得强调的方向性一致证据之一。

## 论文写法建议

可以写：

- The frozen EPDRI score showed a directionally positive association with PFS1 in an independent local 31-patient cohort.
- The signal was consistent across model-native and proxy sensitivity analyses, although confidence intervals were wide.
- The validation supports feasibility and calibration of the prior-locked framework.

不建议写：

- The model is prospectively validated.
- The model improves clinical outcomes.
- The model is superior as a top-1 recommender.

## 新增文件

- `scripts/run_clinical_validation.py`
- `data/clinical_holdout/clinical_features_31.csv`
- `data/clinical_holdout/clinical_outcomes_31.csv`
- `reports/clinical_validation/clinical_validation_report.md`
- `reports/clinical_validation/clinical_patient_predictions_31.csv`
- `reports/clinical_validation/clinical_drug_rankings_31.csv`
- `reports/clinical_validation/clinical_mutation_classifications_31.csv`
- `reports/clinical_validation/clinical_validation_metrics.csv`
- `reports/clinical_validation/clinical_validation_km_logrank_summary.csv`
- `reports/clinical_validation/clinical_validation_bootstrap_distributions.csv`
- `reports/clinical_validation/clinical_validation_consistency_n24_n29.csv`
- `reports/clinical_validation/clinical_validation_subgroup_metrics.csv`
- `reports/clinical_validation/clinical_validation_subgroup_summary.csv`
- `reports/clinical_validation/clinical_validation_score_vs_pfs.png`
- `reports/clinical_validation/clinical_validation_pfs_by_rank.png`
- `reports/clinical_validation/clinical_validation_km_high_low_epdri.png`
- `reports/clinical_validation/clinical_validation_bootstrap_distributions.png`
- `docs/clinical_validation_results_addendum_2026-06-20.md`
- `scripts/run_basic_pacc_docking.py`
- `requirements-docking.txt`
- `reports/docking_pacc_core/docking_scores.csv`
- `reports/docking_pacc_core/basic_pacc_docking_report.md`
- `reports/docking_pacc_core/basic_pacc_docking_heatmap.png`
- `reports/docking_pacc_core/basic_pacc_docking_pose_overview.png`
- `docs/docking_methods_results_addendum_2026-06-20.md`

## 复现命令

```bash
PYTHONPATH=src python scripts/run_clinical_validation.py --bootstrap 5000
```

在当前 Windows/Codex 运行时中使用的是 bundled Python。若在标准本地环境运行，请先安装：

```bash
pip install -r requirements.txt
```
