# Leave-One-Out Back-Testing Report

- Total locked literature rows: 34
- ORR-evaluable folds used for metrics: 29
- Candidate drugs: 12
- Bootstrap replicates: 1000
- Outcome label: held-out mutation-drug ORR row; one positive candidate per fold.
- Non-ORR rows are retained in fold audit but excluded from AUC/top-1 metrics.

## Metric Summary

| method     |   n_total_locked_rows |   n_orr_evaluable_folds |   candidate_drug_count |   pooled_drug_auc |   delta_auc_vs_baseline_0 |   delta_auc_vs_baseline_1 |   top1_hit_rate |   mean_target_rank | metric_definition                                                                                            |
|:-----------|----------------------:|------------------------:|-----------------------:|------------------:|--------------------------:|--------------------------:|----------------:|-------------------:|:-------------------------------------------------------------------------------------------------------------|
| model      |                    34 |                      29 |                     12 |          0.878662 |                 0.0753432 |                 0.0184304 |        0.586207 |            2.27586 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |
| baseline_0 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                -0.0569128 |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |
| baseline_1 |                    34 |                      29 |                     12 |          0.860231 |                 0.0569128 |                 0         |        0.586207 |            1.93103 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |

## Bootstrap 95% CI

| method     | metric                  |   bootstrap_n |       mean |   ci95_lower |   ci95_upper |
|:-----------|:------------------------|--------------:|-----------:|-------------:|-------------:|
| baseline_0 | pooled_drug_auc         |          1000 |  0.801817  |   0.694326   |   0.899805   |
| baseline_0 | delta_auc_vs_baseline_0 |          1000 |  0         |   0          |   0          |
| baseline_0 | delta_auc_vs_baseline_1 |          1000 | -0.0554802 |  -0.109623   |  -0.00385634 |
| baseline_0 | top1_hit_rate           |          1000 |  0.582069  |   0.37931    |   0.758621   |
| baseline_1 | pooled_drug_auc         |          1000 |  0.857297  |   0.764731   |   0.932496   |
| baseline_1 | delta_auc_vs_baseline_0 |          1000 |  0.0554802 |   0.00385634 |   0.109623   |
| baseline_1 | delta_auc_vs_baseline_1 |          1000 |  0         |   0          |   0          |
| baseline_1 | top1_hit_rate           |          1000 |  0.582069  |   0.37931    |   0.758621   |
| model      | pooled_drug_auc         |          1000 |  0.875861  |   0.795369   |   0.938331   |
| model      | delta_auc_vs_baseline_0 |          1000 |  0.0740443 |   0.00302129 |   0.160672   |
| model      | delta_auc_vs_baseline_1 |          1000 |  0.0185641 |  -0.026862   |   0.078878   |
| model      | top1_hit_rate           |          1000 |  0.582069  |   0.37931    |   0.758621   |

## Ablation Summary

| scenario                | method     |   n_total_locked_rows |   n_orr_evaluable_folds |   candidate_drug_count |   pooled_drug_auc |   delta_auc_vs_baseline_0 |   delta_auc_vs_baseline_1 |   top1_hit_rate |   mean_target_rank | metric_definition                                                                                            |   delta_auc_vs_full_model |   delta_top1_vs_full_model |   delta_mean_rank_vs_full_model |
|:------------------------|:-----------|----------------------:|------------------------:|-----------------------:|------------------:|--------------------------:|--------------------------:|----------------:|-------------------:|:-------------------------------------------------------------------------------------------------------------|--------------------------:|---------------------------:|--------------------------------:|
| full_model              | model      |                    34 |                      29 |                     12 |          0.878662 |                 0.0753432 |                 0.0184304 |        0.586207 |            2.27586 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| full_model              | baseline_0 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                -0.0569128 |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| full_model              | baseline_1 |                    34 |                      29 |                     12 |          0.860231 |                 0.0569128 |                 0         |        0.586207 |            1.93103 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_oncokb               | model      |                    34 |                      29 |                     12 |          0.878662 |                 0.0753432 |                 0.0753432 |        0.586207 |            2.27586 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_oncokb               | baseline_0 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                 0         |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_oncokb               | baseline_1 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                 0         |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                -0.0569128 |                  0         |                     -0.206897   |
| no_exact_evidence_floor | model      |                    34 |                      29 |                     12 |          0.815588 |                 0.0122689 |                -0.0446438 |        0.37931  |            3       | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                -0.0630743 |                 -0.206897  |                      0.724138   |
| no_exact_evidence_floor | baseline_0 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                -0.0569128 |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_exact_evidence_floor | baseline_1 |                    34 |                      29 |                     12 |          0.860231 |                 0.0569128 |                 0         |        0.586207 |            1.93103 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_config_fallback      | model      |                    34 |                      29 |                     12 |          0.906388 |                 0.10307   |                 0.0461572 |        0.413793 |            1.65517 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0.0277267 |                 -0.172414  |                     -0.62069    |
| no_config_fallback      | baseline_0 |                    34 |                      29 |                     12 |          0.803319 |                 0         |                -0.0569128 |        0.586207 |            1.72414 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| no_config_fallback      | baseline_1 |                    34 |                      29 |                     12 |          0.860231 |                 0.0569128 |                 0         |        0.586207 |            1.93103 | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0         |                  0         |                      0          |
| pacc_only_subset        | model      |                    30 |                      25 |                     12 |          0.9024   |                 0.115127  |                 0.0540364 |        0.64     |            2       | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                 0.0237382 |                  0.0537931 |                     -0.275862   |
| pacc_only_subset        | baseline_0 |                    30 |                      25 |                     12 |          0.787273 |                 0         |                -0.0610909 |        0.64     |            1.72    | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                -0.0160458 |                  0.0537931 |                     -0.00413793 |
| pacc_only_subset        | baseline_1 |                    30 |                      25 |                     12 |          0.848364 |                 0.0610909 |                 0         |        0.64     |            1.96    | Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member. |                -0.0118677 |                  0.0537931 |                      0.0289655  |

## Subgroup Rank Summary

| subgroup_type   | subgroup           | method     |   n_folds |   mean_target_rank |   top1_hit_rate |   delta_mean_rank_vs_baseline_0 |   delta_top1_vs_baseline_0 |
|:----------------|:-------------------|:-----------|----------:|-------------------:|----------------:|--------------------------------:|---------------------------:|
| class_group     | PACC               | baseline_0 |        25 |            1.72    |        0.64     |                        0        |                          0 |
| class_group     | PACC               | baseline_1 |        25 |            1.96    |        0.64     |                        0.24     |                          0 |
| class_group     | PACC               | model      |        25 |            2       |        0.64     |                        0.28     |                          0 |
| class_group     | non_PACC           | baseline_0 |         4 |            1.75    |        0.25     |                        0        |                          0 |
| class_group     | non_PACC           | baseline_1 |         4 |            1.75    |        0.25     |                        0        |                          0 |
| class_group     | non_PACC           | model      |         4 |            4       |        0.25     |                        2.25     |                          0 |
| compound_group  | compound_or_group  | baseline_0 |        14 |            1.71429 |        0.642857 |                        0        |                          0 |
| compound_group  | compound_or_group  | baseline_1 |        14 |            2.14286 |        0.642857 |                        0.428571 |                          0 |
| compound_group  | compound_or_group  | model      |        14 |            2.21429 |        0.642857 |                        0.5      |                          0 |
| compound_group  | single_or_specific | baseline_0 |        15 |            1.73333 |        0.533333 |                        0        |                          0 |
| compound_group  | single_or_specific | baseline_1 |        15 |            1.73333 |        0.533333 |                        0        |                          0 |
| compound_group  | single_or_specific | model      |        15 |            2.33333 |        0.533333 |                        0.6      |                          0 |

## Outputs

- `loo_predictions.csv`
- `loo_fold_summary.csv`
- `loo_metric_summary.csv`
- `loo_bootstrap_replicates.csv`
- `loo_bootstrap_ci.csv`
- `loo_ablation_summary.csv`
- `loo_subgroup_rank_summary.csv`
