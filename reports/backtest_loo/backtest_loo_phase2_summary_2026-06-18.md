# 第二阶段 Back-Testing 阶段总结

生成时间：2026-06-18 Asia/Singapore

## 目标

将 EGFR-PACC digital cell 从框架描述推进到可验证的方法学模型。

当前实现为 leave-one-out literature back-testing：

- 输入：`data/processed/literature_orr_locked.csv` 的 34 行锁定文献证据。
- 每次 held out 1 行 mutation-drug evidence。
- 使用其余行生成 mutation/drug prior score。
- held-out 行的 ORR 作为已知结局标签。
- 与 Baseline-0 和 Baseline-1 比较 pooled drug AUC、top-1 hit rate 和 mean target rank。

## 已实现脚本

脚本：

- `scripts/run_backtest_loo.py`

运行：

```bash
PYTHONPATH=src .venv/bin/python scripts/run_backtest_loo.py
```

输出目录：

- `reports/backtest_loo/`

三张核心输出表：

1. `loo_predictions.csv`
   - 每个 held-out fold × candidate drug 的模型分数、Baseline-0 分数、Baseline-1 分数。
   - 包含 `label_is_heldout_drug`，用于 pooled AUC。

2. `loo_fold_summary.csv`
   - 每个 fold × method 的 top-1 命中、held-out target rank、top1 drugs。
   - 34 行全部保留审计，但只有 ORR-evaluable fold 进入指标。

3. `loo_metric_summary.csv`
   - model、Baseline-0、Baseline-1 的总体指标。
   - 包含 `pooled_drug_auc`、`delta_auc_vs_baseline_0`、`delta_auc_vs_baseline_1`、
     `top1_hit_rate`、`mean_target_rank`。

补充报告：

- `loo_backtest_report.md`

## Fold 与指标定义

总锁定行：

- 34 行。

ORR-evaluable fold：

- 29 行。

非 ORR 行：

- Robichaux TTF/DOT 或 CNS subgroup PFS 等无 ORR 行保留在 fold audit。
- 不进入 AUC/top-1 指标，避免把非 ORR endpoint 当 ORR outcome。

候选药物数：

- 12 个。

AUC 定义：

- pooled one-positive-per-fold ROC AUC。
- 对每个 fold，held-out drug/member 标记为 positive，其余 candidate drugs 为 negative。
- drug-class 行如 `second_generation_TKI` 会展开为 afatinib/dacomitinib。

Top-1 hit 定义：

- 若 held-out target drug 位于最高分并列组中，记为命中。

## 模型与基线

### Model

使用剩余文献行生成 LOO prior：

- exact mutation-drug evidence 优先。
- class-drug evidence 次之。
- 若无 LOO 文献证据，则回退 locked config prior。

为避免直接文献证据被纯 fallback prior 轻微压过，采用 direct-evidence floor：

```text
exact score = max(config_prior, 0.25 * config_prior + 0.75 * exact_ORR)
class score = max(config_prior, 0.50 * config_prior + 0.50 * class_ORR)
```

### Baseline-0

Robichaux-style 四分类直接推荐：

- PACC：afatinib / dacomitinib。
- classical-like：osimertinib / afatinib / erlotinib / gefitinib。
- T790M-like：osimertinib / furmonertinib。
- ex20ins：amivantamab。

### Baseline-1

Baseline-0 + OncoKB level：

- OncoKB LEVEL_1/2/3/4 转为 actionability score。
- resistance level 不作为敏感证据。

## 当前结果

| method | pooled_drug_auc | ΔAUC vs Baseline-0 | ΔAUC vs Baseline-1 | top-1 hit rate | mean target rank |
| --- | ---: | ---: | ---: | ---: | ---: |
| model | 0.878662 | +0.075343 | +0.018430 | 0.586207 | 2.275862 |
| baseline_0 | 0.803319 | 0.000000 | -0.056913 | 0.586207 | 1.724138 |
| baseline_1 | 0.860231 | +0.056913 | 0.000000 | 0.586207 | 1.931034 |

解释：

- Model 的 pooled AUC 高于 Baseline-0 和 Baseline-1。
- Top-1 hit rate 与两个 baseline 当前相同，均为 0.586。
- Model 的 mean target rank 仍弱于 baseline，提示排序顶部仍受 class fallback 与候选药物集合影响。
- 因样本只有 29 个 ORR-evaluable fold，当前结果应作为方法学 proof-of-concept，不应过度解释为稳定临床优效。

## Bootstrap 95% CI

已对 29 个 ORR-evaluable folds 做 1000 次 fold-level bootstrap。

关键结果：

| method | metric | mean | 95% CI |
| --- | --- | ---: | --- |
| model | pooled AUC | 0.875861 | 0.795369-0.938331 |
| model | ΔAUC vs Baseline-0 | +0.074044 | +0.003021 to +0.160672 |
| model | ΔAUC vs Baseline-1 | +0.018564 | -0.026862 to +0.078878 |
| model | top-1 hit rate | 0.582069 | 0.379310-0.758621 |

解释：

- Model 相对 Baseline-0 的 ΔAUC 95% CI 下界 > 0，可支持“相对
  Robichaux-only baseline 的 AUC 改善”。
- Model 相对 Baseline-1 的 ΔAUC 95% CI 跨 0，因此相对 OncoKB-enhanced
  baseline 只能表述为 exploratory improvement / proof-of-concept signal。
- top-1 hit rate 的 CI 宽，且 model 与 baseline 点估计相同，不应声称 top-1
  命中率改善。

## Ablation 分析

已完成四个 ablation：

1. `no_oncokb`
2. `no_exact_evidence_floor`
3. `no_config_fallback`
4. `pacc_only_subset`

主要发现：

- `no_exact_evidence_floor` 使 model AUC 从 0.878662 降至 0.815588，
  top-1 hit rate 从 0.586207 降至 0.379310，说明 direct evidence floor 是当前
  设计里最关键的性能来源之一。
- `no_oncokb` 对 model 无影响，因为 OncoKB 当前只进入 Baseline-1，不进入 model
  scoring；但 Baseline-1 退化到 Baseline-0。
- `no_config_fallback` 使 model AUC 升至 0.906388、mean target rank 改善到
  1.655，但 top-1 hit rate 降至 0.413793。这说明 config fallback 有助于稳住
  top-1，但也会把部分无直接文献支撑的结构先验药物推得过高。
- `pacc_only_subset` 中 model AUC 为 0.902400，ΔAUC vs Baseline-0 为
  +0.115127，提示模型主要优势集中在 PACC 相关 folds。

## Mean Target Rank 分层解释

分层分析显示：

- PACC folds：model mean target rank 2.00，Baseline-0 为 1.72；top-1 hit rate
  均为 0.64。
- non-PACC folds：model mean target rank 4.00，Baseline-0 为 1.75；但仅 4 个
  folds，不能过度解释。
- compound/group folds：model mean target rank 2.21，Baseline-0 为 1.71；
  top-1 hit rate 均为 0.643。
- single/specific folds：model mean target rank 2.33，Baseline-0 为 1.73；
  top-1 hit rate 均为 0.533。

解释：

- Model 的 AUC 改善主要来自“把 held-out drug 整体排得更靠前”，但不是来自
  top-1/rank 的全面改善。
- Mean target rank 弱于 Baseline 的原因是 Baseline-0 是非常窄的规则推荐集，
  在 afatinib-heavy 文献锁表中天然有 top-rank 优势；model 则同时保留更宽候选
  药物和结构/文献 fallback，因此 rank 更保守。
- 论文中应强调 AUC 改善和可解释 ranking，而不是声称 top-1 ranking 已优于
  Baseline。

## 下一步

1. 增加 bootstrap CI：
   - 对 29 个 ORR-evaluable folds 做 bootstrap。
   - 输出 AUC、ΔAUC、top-1 hit rate 的 95% CI。

2. 增加 endpoint 分层：
   - mutation-drug rows。
   - group-drug rows。
   - compound rows。
   - randomized comparator rows。

3. 增加 ablation：
   - no OncoKB。
   - no exact evidence floor。
   - no config fallback。
   - PACC-only subset。

4. 为论文方法图增加 Back-testing panel：
   - leave-one-out split。
   - Baseline-0/Baseline-1/model 三路评分。
   - AUC/top-1 输出。
