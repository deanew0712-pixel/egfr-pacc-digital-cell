# CSBJ 投稿就绪评估

日期：2026-06-21

## 结论

本项目可以投 Computational and Structural Biotechnology Journal (CSBJ)，但定位应是：

> structure-informed EGFR-PACC digital-cell drug-ranking framework with retrospective clinical validation and supplementary structural plausibility analyses.

不要定位成单纯 docking 文章，也不要写成已经具备前瞻性临床决策能力。

当前稳妥程度：中等偏稳。  
推荐：可以进入投稿准备；MD-relaxed receptor docking 可作为 Major Revision 阶段的补强材料，不需要现在卡住投稿。

## 已完成内容

1. 数字细胞模型

- EGFR PACC/classical-like/compound mutation 分类规则已实现。
- EPDRI 药物排序框架已冻结。
- 文献 leave-one-out back-testing 已完成。
- PACC-only subset 信号强于整体队列，是当前主线中最强的方法学证据。

2. 本地 31 例临床验证

- 主分析 n=24。
- proxy sensitivity n=29。
- n=24：Spearman rho 0.288；C-index 0.595。
- n=29：Spearman rho 0.396；C-index 0.633。
- KM/log-rank 已完成：
  - n=24：high vs low EPDRI median PFS1 16.95 vs 8.20 months，p=0.207。
  - n=29：16.80 vs 7.35 months，p=0.060。
- Bootstrap distribution、subgroup n、subgroup CI、n=24/n=29 consistency 说明均已补齐。

3. 结构补充

- WT/G719S/L861Q/S768I × afatinib/osimertinib/furmonertinib 基础 Vina docking 已完成。
- Meeko receptor preparation 失败问题已解决：receptor PDBQT 改用 OpenBabel，ligand PDBQT 保留 Meeko。
- receptor preparation/protonation workflow sensitivity 已完成：
  - OpenBabel default；
  - OpenBabel Gasteiger；
  - PDBFixer pH 7.4 + OpenBabel Gasteiger。
- covalent warhead geometry screen 已完成，用 acrylamide beta-carbon 到 ATP-pocket cysteine SG 的距离和角度评价 pre-covalent pose plausibility。

4. Colab/AWS GPU 补强模块

- 已准备 Google Colab/AWS 可执行脚本：
  - `scripts/colab_md_relaxed_egfr_pacc.py`
- 已准备本地打包脚本：
  - `scripts/package_colab_md_relaxed_inputs.py`
- 已生成 Colab 上传包：
  - `reports/colab_md_relaxed_package/colab_egfr_pacc_md_input.zip`
- 已生成操作说明：
  - `docs/COLAB_GPU_MD_RELAXED_RUNBOOK_CN.md`

## 为什么 MD-relaxed 不建议现在卡住投稿

本地 Windows CPU 上全 EGFR kinase receptor OpenMM restrained minimization 超时，即使把迭代数降到很低也没有及时输出 relaxed PDB。问题不是力场不可用，而是本地计算环境不适合这一步。

处理策略：

- 当前稿件：把已完成的 docking + workflow sensitivity + pre-covalent geometry 作为结构补充。
- 修稿阶段：用 Colab T4 或云 GPU 跑 restrained minimization + relaxed docking。
- 写法：OpenMM restrained receptor-relaxation sensitivity analysis。
- 不写成 production MD、FEP 或正式 covalent binding free-energy calculation。

## CSBJ 适配度

适合点：

- CSBJ 接受 computational/structural biology 与生物医学转化结合的工作。
- 本项目不是只有小样本临床回顾，而是有模型、文献 back-testing、本地验证和结构解释。
- EGFR-PACC 是明确目标人群，PACC-only subset 的增强信号有文章主线价值。

风险点：

- 本地验证 n=31 小，CI 宽。
- log-rank p 值未达到传统显著阈值。
- docking 是 coarse Vina，不是 FEP/MD。
- subgroup 只能写 hypothesis-generating。

规避方式：

- 全文使用 “retrospective validation / independent local validation / directionally consistent”。
- 不写 “prospective prediction”。
- docking 写作 supplementary structural plausibility。
- MD-relaxed 写作 revision-ready extension。

## 投稿建议

可以投 CSBJ。建议按 regular research article 准备。

主标题方向：

1. Structure-informed digital-cell ranking of EGFR-PACC inhibitor sensitivity with retrospective clinical validation
2. A structure-guided computational framework for prioritizing EGFR inhibitors in PACC and compound EGFR-mutant lung cancer
3. Digital-cell modeling of EGFR-PACC inhibitor response integrates literature back-testing, local clinical validation, and structural docking

## 投稿前最后清单

- Figure 1：数字细胞模型框架。
- Figure 2：literature LOO back-testing，突出 PACC-only subset。
- Figure 3：本地验证，score vs PFS、KM/log-rank、bootstrap distribution。
- Figure 4 或 Supplement：docking heatmap、workflow sensitivity、pre-covalent geometry。
- Supplementary Tables：31 例临床表、mutation class、drug ranking、bootstrap metrics、subgroup metrics、docking scores。
- Methods 明确 frozen model，无本地训练。
- Limitations 明确 retrospective small cohort、CI wide、Vina coarse docking、无 production MD/FEP。
