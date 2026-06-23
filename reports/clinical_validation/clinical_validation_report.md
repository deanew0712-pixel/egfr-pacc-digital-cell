# Local Clinical Validation Report

Generated from `C:\Users\Administrator\Desktop\SupplementaryTables_SUBMISSION_FINAL.xlsx` using frozen posterior prior table `F:\egfr_pacc_digital_cell\reports\pymc_phase1_prior_posterior_4chains_2000draws.csv`.

## Cohort and Analysis Sets

- Total local cases extracted: 31.
- Primary model-native monotherapy analysis set: 24.
- Model-native plus afatinib-combination sensitivity set: 25.
- Proxy sensitivity set including icotinib as first-generation TKI proxy: 29.
- Unsupported first-line treatments excluded from score-PFS validation: 2.

Unsupported treatments:

| patient_id | drug_received | treatment_mapping_note |
| --- | --- | --- |
| Pt-21 | Chemotherapy (platinum-based) | platinum chemotherapy is outside model drug candidates |
| Pt-26 | Chemotherapy (platinum-based) | platinum chemotherapy is outside model drug candidates |

## Primary Metrics

| analysis_set | metric | estimate | ci95_lower | ci95_upper | definition |
| --- | --- | --- | --- | --- | --- |
| model_native_monotherapy | n_evaluable | 24.0 |  |  | Patients with mappable first-line exposure and PFS1. |
| model_native_monotherapy | spearman_score_vs_pfs1 | 0.2877182325581308 | -0.18398650897405636 | 0.6866325523310298 | Spearman correlation between frozen EPDRI score for administered/proxy drug and PFS1 months. |
| model_native_monotherapy | pairwise_c_index | 0.5945454545454546 | 0.4263565891472868 | 0.7518796992481203 | Pairwise concordance of higher EPDRI score with longer PFS1; censoring not modeled. |
| model_native_monotherapy | median_pfs1_months | 15.7 |  |  | Median observed PFS1 in analysis set. |
| model_native_monotherapy | top1_coverage | 0.2916666666666667 |  |  | Fraction where administered/proxy drug is ranked first by frozen model. |
| model_native_monotherapy | top2_coverage | 0.5 |  |  | Fraction where administered/proxy drug is ranked in top 2. |
| model_native_monotherapy | top3_coverage | 0.5416666666666666 |  |  | Fraction where administered/proxy drug is ranked in top 3. |
| model_native_monotherapy | median_pfs_high_minus_low_score | 8.750000000000004 |  |  | Median PFS1 difference using within-set median EPDRI score split. |
| model_native_plus_combo | n_evaluable | 25.0 |  |  | Patients with mappable first-line exposure and PFS1. |
| model_native_plus_combo | spearman_score_vs_pfs1 | 0.2676073927751076 | -0.20836922940938093 | 0.6542256126689046 | Spearman correlation between frozen EPDRI score for administered/proxy drug and PFS1 months. |
| model_native_plus_combo | pairwise_c_index | 0.5836120401337793 | 0.42077464788732394 | 0.735191637630662 | Pairwise concordance of higher EPDRI score with longer PFS1; censoring not modeled. |
| model_native_plus_combo | median_pfs1_months | 15.7 |  |  | Median observed PFS1 in analysis set. |
| model_native_plus_combo | top1_coverage | 0.32 |  |  | Fraction where administered/proxy drug is ranked first by frozen model. |
| model_native_plus_combo | top2_coverage | 0.52 |  |  | Fraction where administered/proxy drug is ranked in top 2. |
| model_native_plus_combo | top3_coverage | 0.56 |  |  | Fraction where administered/proxy drug is ranked in top 3. |
| model_native_plus_combo | median_pfs_high_minus_low_score | 8.600000000000001 |  |  | Median PFS1 difference using within-set median EPDRI score split. |
| proxy_sensitivity_including_icotinib | n_evaluable | 29.0 |  |  | Patients with mappable first-line exposure and PFS1. |
| proxy_sensitivity_including_icotinib | spearman_score_vs_pfs1 | 0.39614361571403145 | -0.0005048553980595689 | 0.7112283272037574 | Spearman correlation between frozen EPDRI score for administered/proxy drug and PFS1 months. |
| proxy_sensitivity_including_icotinib | pairwise_c_index | 0.6333333333333333 | 0.49743589743589745 | 0.7570332480818415 | Pairwise concordance of higher EPDRI score with longer PFS1; censoring not modeled. |
| proxy_sensitivity_including_icotinib | median_pfs1_months | 11.6 |  |  | Median observed PFS1 in analysis set. |
| proxy_sensitivity_including_icotinib | top1_coverage | 0.27586206896551724 |  |  | Fraction where administered/proxy drug is ranked first by frozen model. |
| proxy_sensitivity_including_icotinib | top2_coverage | 0.4482758620689655 |  |  | Fraction where administered/proxy drug is ranked in top 2. |
| proxy_sensitivity_including_icotinib | top3_coverage | 0.4827586206896552 |  |  | Fraction where administered/proxy drug is ranked in top 3. |
| proxy_sensitivity_including_icotinib | median_pfs_high_minus_low_score | 9.450000000000001 |  |  | Median PFS1 difference using within-set median EPDRI score split. |

## Kaplan-Meier and Log-Rank

| analysis_set | split_variable | median_score_cutpoint | n_high | n_low | events_high | events_low | median_pfs_high | median_pfs_low | logrank_chi2 | logrank_p_value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| model_native_monotherapy | administered_drug_epdri_median | 0.679 | 14 | 10 | 10 | 8 | 16.950000000000003 | 8.2 | 1.5942290227947729 | 0.20672296352266062 |
| proxy_sensitivity_including_icotinib | administered_drug_epdri_median | 0.679 | 15 | 14 | 11 | 12 | 16.8 | 7.35 | 3.542580940726299 | 0.05981232339576183 |

## Bootstrap Distribution Consistency

| metric | primary_n24_estimate | primary_ci95 | sensitivity_n29_estimate | sensitivity_ci95 | direction_consistent | ci_overlap_interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| spearman_score_vs_pfs1 | 0.2877182325581308 | -0.18398650897405636 to 0.6866325523310298 | 0.39614361571403145 | -0.0005048553980595689 to 0.7112283272037574 | True | CIs overlap; sensitivity analysis preserves the positive direction. |
| pairwise_c_index | 0.5945454545454546 | 0.4263565891472868 to 0.7518796992481203 | 0.6333333333333333 | 0.49743589743589745 to 0.7570332480818415 | True | CIs overlap; sensitivity analysis preserves the positive direction. |

The n=24 primary model-native analysis and n=29 proxy sensitivity analysis are directionally consistent. Both Spearman and pairwise C-index remain positive, and their bootstrap confidence intervals overlap. This supports robustness of the local signal to inclusion of icotinib as a first-generation EGFR-TKI proxy, while preserving the primary analysis as the cleaner model-native estimate.

## Subgroup Metrics

| analysis_set | subgroup_family | subgroup | n | spearman_score_vs_pfs1 | spearman_ci95_lower | spearman_ci95_upper | pairwise_c_index | cindex_ci95_lower | cindex_ci95_upper | median_actual_model_score | median_pfs1_months | top3_coverage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| model_native_monotherapy | PACC vs non-PACC | classical_like | 8 | -0.6801842415181628 | -0.9058216273156765 | -0.17429786813223805 | 0.2777777777777778 | 0.15217391304347827 | 0.5 | 0.679 | 19.6 | 0.0 |
| model_native_monotherapy | PACC vs non-PACC | pacc | 16 | 0.31931298801289854 | -0.21710130434855482 | 0.749871241152025 | 0.6083333333333333 | 0.4212962962962963 | 0.7850877192982456 | 0.666 | 8.75 | 0.8125 |
| model_native_monotherapy | single vs compound mutation | Compound mutation | 9 | 0.30276503540974914 | -0.49319696191607193 | 0.8583950752789521 | 0.5972222222222222 | 0.3387096774193548 | 0.8333333333333334 | 0.749 | 15.9 | 0.8888888888888888 |
| model_native_monotherapy | single vs compound mutation | Single mutation | 15 | 0.22371654858937365 | -0.4760165546487071 | 0.774643448368334 | 0.5673076923076923 | 0.3177083333333333 | 0.7806122448979592 | 0.679 | 15.7 | 0.3333333333333333 |
| proxy_sensitivity_including_icotinib | PACC vs non-PACC | classical_like | 8 | -0.6801842415181628 | -0.9058216273156765 | -0.17541160386140583 | 0.2777777777777778 | 0.14705882352941177 | 0.5 | 0.679 | 19.6 | 0.0 |
| proxy_sensitivity_including_icotinib | PACC vs non-PACC | pacc | 21 | 0.39444258845230645 | -0.03889693200247054 | 0.7402492961657883 | 0.638095238095238 | 0.4876237623762376 | 0.78 | 0.666 | 8.7 | 0.6666666666666666 |
| proxy_sensitivity_including_icotinib | single vs compound mutation | Compound mutation | 13 | 0.5241424183609592 | -0.03863712414277084 | 0.8619050770310418 | 0.6730769230769231 | 0.4855072463768116 | 0.8356164383561644 | 0.666 | 11.3 | 0.6923076923076923 |
| proxy_sensitivity_including_icotinib | single vs compound mutation | Single mutation | 16 | 0.25118348964921217 | -0.3931147550474972 | 0.732927649323365 | 0.5756302521008403 | 0.35046728971962615 | 0.7616822429906542 | 0.679 | 13.649999999999999 | 0.3125 |

## Subgroup Interpretation Strength

| analysis_set | subgroup_family | subgroup | n | spearman_score_vs_pfs1 | spearman_ci95_lower | spearman_ci95_upper | pairwise_c_index | cindex_ci95_lower | cindex_ci95_upper | interpretation_strength | recommended_language |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| model_native_monotherapy | PACC vs non-PACC | classical_like | 8 | -0.6801842415181628 | -0.9058216273156765 | -0.17429786813223805 | 0.2777777777777778 | 0.15217391304347827 | 0.5 | very small; descriptive only | Report as descriptive only because subgroup n < 10. |
| model_native_monotherapy | PACC vs non-PACC | pacc | 16 | 0.31931298801289854 | -0.21710130434855482 | 0.749871241152025 | 0.6083333333333333 | 0.4212962962962963 | 0.7850877192982456 | exploratory; hypothesis-generating | Directionally positive subgroup signal; present as hypothesis-generating. |
| model_native_monotherapy | single vs compound mutation | Compound mutation | 9 | 0.30276503540974914 | -0.49319696191607193 | 0.8583950752789521 | 0.5972222222222222 | 0.3387096774193548 | 0.8333333333333334 | very small; descriptive only | Report as descriptive only because subgroup n < 10. |
| model_native_monotherapy | single vs compound mutation | Single mutation | 15 | 0.22371654858937365 | -0.4760165546487071 | 0.774643448368334 | 0.5673076923076923 | 0.3177083333333333 | 0.7806122448979592 | exploratory; hypothesis-generating | Directionally positive subgroup signal; present as hypothesis-generating. |
| proxy_sensitivity_including_icotinib | PACC vs non-PACC | classical_like | 8 | -0.6801842415181628 | -0.9058216273156765 | -0.17541160386140583 | 0.2777777777777778 | 0.14705882352941177 | 0.5 | very small; descriptive only | Report as descriptive only because subgroup n < 10. |
| proxy_sensitivity_including_icotinib | PACC vs non-PACC | pacc | 21 | 0.39444258845230645 | -0.03889693200247054 | 0.7402492961657883 | 0.638095238095238 | 0.4876237623762376 | 0.78 | exploratory; hypothesis-generating | Directionally positive subgroup signal; present as hypothesis-generating. |
| proxy_sensitivity_including_icotinib | single vs compound mutation | Compound mutation | 13 | 0.5241424183609592 | -0.03863712414277084 | 0.8619050770310418 | 0.6730769230769231 | 0.4855072463768116 | 0.8356164383561644 | small; hypothesis-generating | Directionally positive subgroup signal; present as hypothesis-generating. |
| proxy_sensitivity_including_icotinib | single vs compound mutation | Single mutation | 16 | 0.25118348964921217 | -0.3931147550474972 | 0.732927649323365 | 0.5756302521008403 | 0.35046728971962615 | 0.7616822429906542 | exploratory; hypothesis-generating | Directionally positive subgroup signal; present as hypothesis-generating. |

Subgroups with n < 10 should be described as descriptive only. Other local subgroup findings remain hypothesis-generating because confidence intervals are wide and the cohort was not powered for formal interaction testing.

## Interpretation

- This validation uses PFS1, not ORR; the original local table does not contain response categories.
- PFS1 event/censoring status is inferred only from PD type or second-line therapy fields, so formal survival modeling is not used as the primary analysis.
- The primary analysis is restricted to model-native first-line EGFR-TKI monotherapy to avoid post hoc model expansion.
- Icotinib is handled only in sensitivity analysis because it was not part of the frozen model candidate drug set.
- Platinum chemotherapy is outside the EGFR-TKI ranking model and is reported as unsupported rather than forced into the ranking.

## Outputs

- `clinical_patient_predictions_31.csv`
- `clinical_drug_rankings_31.csv`
- `clinical_mutation_classifications_31.csv`
- `clinical_validation_metrics.csv`
- `clinical_validation_bootstrap_distributions.csv`
- `clinical_validation_km_curves.csv`
- `clinical_validation_km_logrank_summary.csv`
- `clinical_validation_consistency_n24_n29.csv`
- `clinical_validation_subgroup_metrics.csv`
- `clinical_validation_subgroup_interpretation.csv`
- `clinical_validation_subgroup_summary.csv`
- `clinical_validation_score_vs_pfs.png/pdf/svg`
- `clinical_validation_pfs_by_rank.png/pdf/svg`
- `clinical_validation_km_high_low_epdri.png/pdf/svg`
- `clinical_validation_bootstrap_distributions.png/pdf/svg`
