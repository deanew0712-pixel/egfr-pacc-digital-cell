# Local Clinical Validation Addendum

Generated: 2026-06-20 Asia/Shanghai

## Validation Design

The frozen EGFR-PACC digital-cell model was applied to a de-identified local
31-patient treatment cohort from
`C:\Users\Administrator\Desktop\SupplementaryTables_SUBMISSION_FINAL.xlsx`.
The source workbook contains first-line treatment, EGFR mutation subtype, PFS1,
second-line treatment, PD type, and docking/MM-GBSA annotations. It does not
contain objective response categories, brain metastasis status, or co-mutation
fields. Therefore, local validation was prespecified as an exploratory PFS1
association analysis rather than an ORR or formal survival validation.

The workbook was parsed into physically separated files:

- `data/clinical_holdout/clinical_features_31.csv`
- `data/clinical_holdout/clinical_outcomes_31.csv`

Model predictions used the frozen posterior prior table
`reports/pymc_phase1_prior_posterior_4chains_2000draws.csv`. No local outcome
data were used to modify priors, weights, evidence caps, mutation rules, or drug
ranking logic.

## Analysis Sets

Three analysis sets were reported:

1. Primary model-native monotherapy set: 24 patients treated with a first-line
   drug represented directly in the frozen model candidate set.
2. Model-native plus combination sensitivity set: 25 patients, adding one
   BCP plus afatinib case mapped to afatinib as a sensitivity exposure.
3. Proxy sensitivity set: 29 patients, adding four icotinib-treated patients as
   first-generation EGFR-TKI proxy exposures, represented by the higher of the
   gefitinib and erlotinib model scores.

Two platinum-chemotherapy first-line cases were excluded from score-PFS
validation because chemotherapy is outside the EGFR-TKI ranking model.

## Results

In the primary model-native monotherapy set (n = 24), the frozen EPDRI score for
the administered first-line TKI showed a positive but imprecise association with
PFS1: Spearman rho 0.288, bootstrap 95% CI -0.184 to 0.687. Pairwise
concordance of higher EPDRI score with longer PFS1 was 0.595, bootstrap 95% CI
0.426 to 0.752. The administered drug was ranked first in 29.2% of patients,
within the top two in 50.0%, and within the top three in 54.2%. Median PFS1 was
15.7 months. Patients above the median administered-drug EPDRI score had a
median PFS1 advantage of 8.75 months over those below the median score.
Kaplan-Meier analysis by median administered-drug EPDRI score showed separation
between high- and low-score groups, with median PFS1 16.95 versus 8.20 months
and log-rank p = 0.207.

In the model-native plus afatinib-combination sensitivity set (n = 25),
Spearman rho was 0.268, bootstrap 95% CI -0.208 to 0.654. Pairwise concordance
was 0.584, bootstrap 95% CI 0.421 to 0.735. Top-1, top-2, and top-3 coverage
were 32.0%, 52.0%, and 56.0%, respectively.

In the proxy sensitivity set including icotinib (n = 29), Spearman rho increased
to 0.396, bootstrap 95% CI -0.001 to 0.711. Pairwise concordance was 0.633,
bootstrap 95% CI 0.497 to 0.757. Top-1, top-2, and top-3 coverage were 27.6%,
44.8%, and 48.3%, respectively. Kaplan-Meier separation was stronger in this
sensitivity set, with median PFS1 16.80 versus 7.35 months for high versus low
administered-drug EPDRI score and log-rank p = 0.060.

The n = 24 primary analysis and n = 29 proxy sensitivity analysis were formally
compared. Both Spearman and pairwise C-index estimates remained positive, and
their bootstrap confidence intervals overlapped. The sensitivity analysis
therefore supports robustness of the direction of effect while preserving the
model-native monotherapy set as the cleaner primary estimate.

Subgroup analyses were performed by PACC versus non-PACC primary class and by
single versus compound mutation status. In the proxy sensitivity set, the PACC
subgroup (n = 21) showed a positive score-PFS1 association, with Spearman rho
0.394 and pairwise C-index 0.638. The compound-mutation subgroup (n = 13) showed
an even stronger exploratory signal, with Spearman rho 0.524 and pairwise
C-index 0.673. These subgroup results are underpowered and should be described
as hypothesis-generating rather than independent confirmatory evidence. The
classical-like subgroup contained only eight patients, and the model-native
compound subgroup contained nine patients; these strata should be reported as
descriptive only. Taken together, the local subgroup pattern directionally
aligns with the literature leave-one-out finding that the framework is most
informative in the intended PACC-centered setting.

## Interpretation

The local validation provides an independent, outcome-unlocked feasibility test
of the frozen model. The association between frozen EPDRI score and PFS1 is
directionally positive across all analysis sets, with the strongest signal in
the proxy sensitivity set. However, confidence intervals remain wide because
the cohort is small, treatment exposures are heterogeneous, and event/censoring
status is not encoded as a formal survival endpoint. These results should be
reported as exploratory retrospective validation and calibration support, not
as definitive clinical utility validation.

The most defensible manuscript claim is that the model was frozen before local
outcome analysis and showed a consistent positive score-PFS1 signal in a small
local EGFR uncommon/PACC cohort, while retaining transparent handling of
unsupported therapies and non-model-native drugs.

## Generated Outputs

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
