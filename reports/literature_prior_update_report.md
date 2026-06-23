# 文献锁表与 Priors 更新报告

生成时间：2026-06-15 Asia/Singapore

## 已完成

- 已审核 5 篇用户上传 PDF。
- 已复制纳入原始 PDF：`data/raw/literature_pdfs/`
- 已生成可检索文本：`reports/literature_pdf_text/`
- 已生成 PDF 审核报告：`reports/literature_pdf_review_2026-06-15.md`
- 已重写锁表：`data/processed/literature_orr_locked.csv`
- 已更新：`configs/priors.yaml`
- 已更新：`configs/drug_properties.yaml`
- 已重建：`configs/open_fallback_digital_cell.yaml`

## 锁表状态

当前锁表包含 34 行 mutation x drug、group x drug、structure-class 或
comparator evidence。

纳入来源：

- LUX-Lung 2/3/6 post-hoc pooled analysis：4 行，已由 primary PDF 原文锁定。
- Afatinib uncommon mutation updated database：5 行，已修正此前列位错位。
- Robichaux 2021 Nature：5 行，PACC 结构功能分类与 TKI class TTF/DOT evidence。
- UNICORN phase 2 first-line osimertinib uncommon EGFR：10 行。
- ACHILLES/TORG1834 randomized afatinib vs chemotherapy：10 行。

排除来源：

- Neoadjuvant pembrolizumab in MSI-H/dMMR solid tumors：不进入 EGFR-PACC TKI
  prior 锁表，仅作为背景免疫治疗文献保存。

## 关键更新

- Robichaux 支持 PACC second-generation TKI 优先：PACC 二代 TKI mTTF 21.7
  months，对比一代 10.0 months、三代 4.1 months。
- ACHILLES/TORG1834 提供随机证据：afatinib mPFS 10.6 months，对比化疗
  5.7 months，HR 0.421。
- LUX-Lung 2/3/6 primary PDF 提供 afatinib 经典少见突变证据：Group 1 ORR
  71.1%，mDOR 11.1 months，mPFS 10.7 months；G719X/L861Q/S768I 分别为
  ORR 77.8%/56.3%/100.0%，mPFS 13.8/8.2/14.7 months。
- ACHILLES compound uncommon subgroup：afatinib ORR 72.7%，mPFS 21.0 months。
- UNICORN overall uncommon EGFR：osimertinib ORR 55.0%，mPFS 9.4 months。
- UNICORN G719X：osimertinib ORR 45.0%，mPFS 5.1 months；solitary G719X 更低。
- UNICORN L861Q：osimertinib ORR 75.0%，mPFS 22.7 months。
- UNICORN S768I：osimertinib ORR 50.0%，mPFS 9.4 months。
- UNICORN E709X：osimertinib ORR 83.3%，但主文未锁定 PFS。

## Prior 更新逻辑

- `pacc.afatinib` 继续作为主 prior，来源从 LUX/数据库增强为 Robichaux +
  ACHILLES + afatinib database 的一致支持。
- `pacc.dacomitinib` 保留二代 TKI structure-class extrapolation，但仍缺直接
  PACC clinical ORR/PFS 锁表。
- `pacc.osimertinib` 保留为可选 prior，但对 G719X/PACC-like mutation 加 caution。
- CNS modifier 继续使用 FDA label/openFDA evidence，而不是 DrugBank。

## 当前验证结果

以下命令已成功运行：

```bash
python3 scripts/update_priors_from_locked_literature.py
python3 scripts/build_preliminary_digital_cell.py
python3 scripts/run_demo_case.py
python3 scripts/run_batch_cases.py
python3 -m py_compile scripts/*.py src/egfr_pacc/*.py
```

初步 digital-cell build 仍为 9 项检查通过。

## 限制

- LUX-Lung、Robichaux、ACHILLES/TORG1834、UNICORN 核心文献均已完成本地
  PDF 原文层面的 prior 锁表；进入本地 30-37 例结局前仍建议人工签字 freeze。
- UNICORN E709X 目前只锁 ORR，未锁 PFS。
- ACHILLES S768I single mutation n=3，不能过度解释。
- 本锁表仍是 clinical outcome review 前的 research-use prior candidate。
