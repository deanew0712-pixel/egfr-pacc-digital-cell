# EGFR-PACC 数字细胞项目阶段性检查点

日期：2026-06-21  
用途：后续继续分析、修稿或投稿前，优先阅读本文件，避免从头重新梳理项目。  
当前定位：可进入 CSBJ 投稿准备；结构计算作为机制/可信度补充，不作为单独强结论。

## 一句话结论

本项目已经形成一个可投稿的完整闭环：冻结的 EGFR-PACC 数字细胞/EPDRI 药物排序模型、文献 leave-one-out back-testing、本地 31 例临床验证、KM/log-rank 与 bootstrap 稳健性展示、PACC/compound 亚组分析，以及基础 docking、receptor workflow sensitivity、pre-covalent geometry screen 和 Colab GPU OpenMM-relaxed docking 补充。

最合适的论文叙事不是“docking 证明药物结合”，而是：

> A structure-informed, literature-locked digital-cell framework for EGFR-PACC inhibitor prioritization, supported by literature back-testing, independent local clinical validation, and supplementary structural plausibility analyses.

## 当前完成度

### 1. 模型与文献 back-testing

已完成内容：

- EGFR uncommon mutation 分类：PACC、classical-like、T790M/resistance-like、exon20 insertion-like。
- 复合突变规则：保留原始 mutation set，同时按预设规则 collapse 到主分析 class。
- EPDRI 药物排序：整合 class-level prior、clinical modifiers、uncertainty discount、evidence cap。
- 文献 prior lock：本地临床结局没有进入模型训练或调参。
- 文献 leave-one-out back-testing：34 条 locked literature rows，其中 29 条 ORR-evaluable folds。

关键结果：

- Full model pooled drug AUC：0.8787。
- Baseline-0 AUC：0.8033；模型 ΔAUC +0.0753。
- Baseline-1 AUC：0.8602；模型 ΔAUC +0.0184。
- Top-1 hit rate：模型、Baseline-0、Baseline-1 均为 0.5862。
- PACC-only subset：模型 AUC 0.9024；相对 Baseline-0 ΔAUC +0.1151；相对 Baseline-1 ΔAUC +0.0540。

主要文件：

- `reports/backtest_loo/loo_metric_summary.csv`
- `reports/backtest_loo/loo_ablation_summary.csv`
- `reports/backtest_loo/loo_subgroup_rank_summary.csv`
- `reports/backtest_loo/backtest_loo_phase2_summary_2026-06-18.md`
- `reports/figures_publication/figure3_backtesting.*`

写作强度：

- 可以说模型提高了总体 drug-discrimination recovery，尤其相对 Robichaux-only structural baseline。
- 不要说模型已经证明 top-1 推荐优于所有 rule-based baseline，因为 top-1 hit rate 未提高。
- PACC-only subset 是本研究主线最强的方法学信号，但仍应写为 prespecified/intended-population signal 或 hypothesis-supporting evidence。

## 2. 本地 31 例临床验证

已完成内容：

- 本地治疗数据 31 例已整理并映射到模型输入。
- 模型-native monotherapy 主分析：n=24。
- model-native plus combo：n=25。
- proxy sensitivity including icotinib：n=29。
- unsupported platinum chemotherapy：2 例，不纳入核心模型药物评价。
- 已补齐 Spearman、pairwise C-index、bootstrap CI、KM curve、log-rank test、n=24/n=29 consistency、PACC 与 compound 亚组 CI。

主分析 n=24：

- Spearman rho：0.288，95% CI -0.184 至 0.687。
- Pairwise C-index：0.595，95% CI 0.426 至 0.752。
- Top-1/Top-2/Top-3 coverage：29.2%、50.0%、54.2%。
- High vs low EPDRI median PFS1：16.95 vs 8.20 months。
- Log-rank p：0.207。
- Median high-low EPDRI PFS difference：8.75 months。

Proxy sensitivity n=29：

- Spearman rho：0.396，95% CI -0.001 至 0.711。
- Pairwise C-index：0.633，95% CI 0.497 至 0.757。
- Top-1/Top-2/Top-3 coverage：27.6%、44.8%、48.3%。
- High vs low EPDRI median PFS1：16.80 vs 7.35 months。
- Log-rank p：0.060。
- Median high-low EPDRI PFS difference：9.45 months。

一致性解释：

- n=24 与 n=29 方向一致：Spearman、C-index、median PFS separation 均朝同一方向。
- CI 较宽是小样本回顾性验证的预期结果。
- n=29 proxy sensitivity 的信号更强，支持结果对 icotinib proxy 映射具有一定稳健性。
- 不应写成显著性确证，应写成 independent retrospective validation with directionally consistent discrimination and PFS separation。

亚组结果：

- PACC primary n=16：Spearman 0.319，95% CI -0.217 至 0.750；C-index 0.608，95% CI 0.421 至 0.785。
- PACC proxy n=21：Spearman 0.394，95% CI -0.039 至 0.740；C-index 0.638，95% CI 0.488 至 0.780。
- Compound primary n=9：Spearman 0.303，95% CI -0.493 至 0.858；C-index 0.597，95% CI 0.339 至 0.833。
- Compound proxy n=13：Spearman 0.524，95% CI -0.039 至 0.862；C-index 0.673，95% CI 0.486 至 0.836。
- Classical-like subgroup n=8 方向相反，提示模型主要适配 PACC target population，而非所有 uncommon EGFR。

亚组写法：

- PACC 与 compound 亚组样本量小，尤其 compound primary n=9。
- 文中应明确写为 hypothesis-generating subgroup signal。
- 论点可以是：文献 LOO 与本地病例均提示模型信号在 PACC target population 中更强，而不是亚组已独立验证。

主要文件：

- `reports/clinical_validation/clinical_patient_predictions_31.csv`
- `reports/clinical_validation/clinical_drug_rankings_31.csv`
- `reports/clinical_validation/clinical_validation_metrics.csv`
- `reports/clinical_validation/clinical_validation_subgroup_metrics.csv`
- `reports/clinical_validation/clinical_validation_subgroup_interpretation.csv`
- `reports/clinical_validation/clinical_validation_km_logrank_summary.csv`
- `reports/clinical_validation/clinical_validation_bootstrap_distributions.csv`
- `reports/clinical_validation/clinical_validation_figure_qc_2026-06-20.md`
- `reports/figures_publication/figure4_local_validation_score_vs_pfs.*`
- `reports/figures_publication/figure4b_local_validation_km_high_low_epdri.*`
- `reports/figures_publication/figureS_clinical_validation_bootstrap_distributions.*`

## 3. 结构 docking 与 sensitivity

### 基础 Vina docking

已完成内容：

- 对 WT、G719S、L861Q、S768I EGFR kinase domain 与 afatinib、osimertinib、furmonertinib 做基础非共价 Vina docking。
- Receptor template：4LRM chain A。
- Ligand/receptor 处理策略：ligand 使用 Meeko；receptor 因 Meeko Windows/Python 残基价态解析失败，改用 OpenBabel PDBQT。

基础 docking 排名：

- WT：osimertinib -9.684，furmonertinib -9.382，afatinib -9.348 kcal/mol。
- G719S：osimertinib -9.430，furmonertinib -9.141，afatinib -9.001 kcal/mol。
- L861Q：furmonertinib -9.591，afatinib -9.071，osimertinib -8.890 kcal/mol。
- S768I：furmonertinib -9.573，afatinib -9.358，osimertinib -9.134 kcal/mol。

主要文件：

- `scripts/run_basic_pacc_docking.py`
- `reports/docking_pacc_core/docking_scores.csv`
- `reports/docking_pacc_core/basic_pacc_docking_report.md`
- `reports/docking_pacc_core/receptor_preparation_resolution_2026-06-20.md`
- `reports/figures_publication/figureS_basic_pacc_docking_heatmap.*`
- `reports/figures_publication/figureS_basic_pacc_docking_pose_overview.*`

写作限制：

- 只能写 coarse non-covalent docking。
- 不要写成 binding free energy 或临床疗效证明。
- 可作为结构可解释性/补充机制合理性。

### Receptor protonation/workflow sensitivity 与 covalent screen

已完成内容：

- 三种 receptor workflow：
  - OpenBabel default。
  - OpenBabel Gasteiger。
  - PDBFixer pH 7.4 + OpenBabel Gasteiger。
- Pre-covalent warhead geometry screen：用 acrylamide beta-carbon 到 Cys797/模板 A:CYS800:SG 的距离与角度评估可接近性。

workflow 稳健性：

- G719S-osimertinib：3/3 workflows rank 1。
- S768I-furmonertinib：3/3 workflows rank 1。
- L861Q-furmonertinib：2/3 workflows rank 1。
- WT-osimertinib：2/3 workflows rank 1。

covalent screen：

- 大多数 ligand-mutant pair 有 near-attack pose。
- 例外或边缘：G719S-afatinib angle 68.5°，L861Q-furmonertinib angle 65.6°，低于 70° cutoff。

主要文件：

- `scripts/run_docking_sensitivity_and_covalent_screen.py`
- `reports/docking_workflow_sensitivity/receptor_workflow_sensitivity_scores.csv`
- `reports/docking_workflow_sensitivity/receptor_workflow_sensitivity_summary.csv`
- `reports/docking_workflow_sensitivity/covalent_warhead_geometry_screen.csv`
- `reports/docking_workflow_sensitivity/covalent_warhead_geometry_best_poses.csv`
- `reports/docking_workflow_sensitivity/workflow_sensitivity_and_covalent_screen_report.md`
- `reports/figures_publication/figureS_docking_workflow_and_covalent_sensitivity.*`

写作限制：

- 这是 pre-reaction geometry screen，不是真正 tethered covalent docking。
- 不报告为 covalent binding energy。

### Colab GPU OpenMM-relaxed docking

已完成内容：

- 已在 Google Colab GPU 上完成 OpenMM receptor restrained minimization + relaxed receptor docking。
- Colab 实际可用 OpenMM platform 为 OpenCL，不是 CUDA；GPU 设备为 Tesla T4。
- 输出 zip 已下载并导入本地项目。

本地导入位置：

- `reports/docking_md_relaxed_gpu/openmm_gpu_relaxed_outputs.zip`
- `reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/`

正式结果：

- WT：furmonertinib -8.636 rank 1；osimertinib -8.567 rank 2；afatinib -8.475 rank 3。
- G719S：furmonertinib -8.851 rank 1；afatinib -8.794 rank 2；osimertinib -8.626 rank 3。
- L861Q：furmonertinib -8.719 rank 1；osimertinib -8.456 rank 2；afatinib -8.392 rank 3。
- S768I：osimertinib -8.744 rank 1；furmonertinib -8.739 rank 2；afatinib -8.437 rank 3。

重要技术修正：

- 原始 `openmm_gpu_relaxed_receptor_metrics.csv` 中的 heavy_atom_rmsd_A 约 23 Å，是 atom-order mismatch 导致，不能引用。
- 应使用修正文件 `openmm_gpu_relaxed_receptor_metrics_corrected_rmsd.csv`。
- 修正后 heavy-atom RMSD 约 0.90-0.92 Å：
  - WT 0.904 Å。
  - G719S 0.905 Å。
  - L861Q 0.915 Å。
  - S768I 0.903 Å。

主要文件：

- `scripts/colab_md_relaxed_egfr_pacc.py`
- `scripts/package_colab_md_relaxed_inputs.py`
- `colab_egfr_pacc_md_relaxed_runner.ipynb`
- `docs/COLAB_GPU_MD_RELAXED_RUNBOOK_CN.md`
- `reports/colab_md_relaxed_package/colab_egfr_pacc_md_input.zip`
- `reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/openmm_gpu_relaxed_docking_scores.csv`
- `reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/openmm_gpu_relaxed_vs_unrelaxed_comparison.csv`
- `reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/openmm_gpu_relaxed_receptor_metrics_corrected_rmsd.csv`
- `reports/figures_publication/figureS_openmm_gpu_relaxed_docking_sensitivity.*`

写作解释：

- Relaxed docking 后分数整体变得不那么负，且部分排名改变，说明 docking 对 receptor conformation 敏感。
- L861Q-furmonertinib 在 unrelaxed 与 relaxed 中均保持强信号。
- WT/G719S/S768I 的 relaxed 排名改变提示不应把单一 docking 分数作为药物疗效结论。
- 可写为 receptor-conformation sensitivity analysis strengthens methodological caution and supports the use of docking as a supplementary plausibility layer。

## 4. 推荐论文结构

建议题目：

Structure-informed digital-cell ranking of EGFR-PACC inhibitor sensitivity with literature back-testing and retrospective clinical validation

主图建议：

- Figure 1：模型框架。
- Figure 2：文献 LOO back-testing 与 PACC-only subset。
- Figure 3：本地 31 例验证，包括 score vs PFS、C-index/Spearman、bootstrap。
- Figure 4：KM/log-rank high vs low EPDRI。

补充图建议：

- Figure S1：demo case ranking。
- Figure S2：PFS by rank 或 top-k coverage。
- Figure S3：basic docking heatmap/pose overview。
- Figure S4：workflow sensitivity + covalent geometry screen。
- Figure S5：OpenMM-relaxed docking sensitivity。

## 5. 投稿定位评估

CSBJ：可投，当前成熟度中等偏稳。

适合原因：

- 有方法学框架，不只是小样本临床观察。
- 有文献回测、本地独立验证和结构解释三层证据。
- PACC 是明确的结构生物学 target population。
- CSBJ 接受 computational + structural biotechnology + translational validation 的组合型工作。

主要风险：

- 本地验证 n=31，CI 宽。
- KM/log-rank 主分析 p=0.207，proxy p=0.060，未达到传统显著。
- Docking 是 Vina-level，不是 FEP/production MD。
- 亚组样本量小。

规避策略：

- 不声称 prospective clinical utility。
- 不写 statistically significant survival benefit。
- 使用 directionally consistent、retrospective validation、hypothesis-generating subgroup signal、supplementary structural plausibility。
- 明确 frozen model/no local retraining。

## 6. 后续不要重复跑的内容

除非输入数据或模型规则改变，否则不要重跑：

- 文献 LOO back-testing。
- 本地 31 例主分析。
- KM/log-rank。
- Bootstrap distributions。
- PACC/compound 亚组 CI。
- 基础 docking。
- receptor workflow sensitivity。
- pre-covalent geometry screen。
- Colab OpenMM-relaxed docking。

需要重跑的情况：

- 本地 clinical Excel 更新或病例数变化。
- EPDRI scoring logic、mutation classifier、drug list 或 priors 发生改变。
- docking box、template receptor、ligand protonation 或 receptor preparation protocol 改变。
- 需要 production-level MD 或真正 covalent docking/FEP，作为后续修订增强。

## 7. 当前最重要的写作原则

- 文章主结论：模型在 PACC target population 中提供可解释、冻结、方向一致的 inhibitor prioritization signal。
- 临床验证定位：small independent retrospective validation，不是确证性临床试验。
- 结构分析定位：supplementary plausibility and sensitivity analysis，不是疗效机制证明。
- 亚组定位：hypothesis-generating，尤其 compound subgroup。
- 最强论点：文献 back-testing 的 PACC-only 增强信号与本地 PACC/proxy 亚组方向一致。

