# Clinical Validation Figure QC

Generated: 2026-06-20 Asia/Shanghai

## Scope

This QC pass reviewed the new local clinical validation figures requested for
publication use:

1. Kaplan-Meier curves with log-rank p values.
2. Spearman and pairwise C-index bootstrap distribution plots.
3. Existing score-PFS scatter plot.
4. Existing PFS-by-rank plot.

## Figure Checks

### KM Curve

Files:

- `reports/clinical_validation/clinical_validation_km_high_low_epdri.png`
- `reports/clinical_validation/clinical_validation_km_high_low_epdri.pdf`
- `reports/clinical_validation/clinical_validation_km_high_low_epdri.svg`
- `reports/figures_publication/figure4b_local_validation_km_high_low_epdri.png`
- `reports/figures_publication/figure4b_local_validation_km_high_low_epdri.pdf`
- `reports/figures_publication/figure4b_local_validation_km_high_low_epdri.svg`

QC result:

- Format is KM step-curve, not median/bar/box plot.
- Log-rank p values are annotated in each panel.
- P-value labels are placed in lower-left panel space with a light background
  box and do not overlap the survival curves.
- Axis labels are explicit: PFS1 months and progression-free probability.
- Legend is readable and does not obscure curves.
- Font sizes are suitable for supplement or multi-panel main figure use.

Recommended use:

- Main Figure 4 panel if local validation is central to the manuscript.
- Otherwise Supplementary Figure with Figure 4 scatter plot in main text.

### Bootstrap Distribution

Files:

- `reports/clinical_validation/clinical_validation_bootstrap_distributions.png`
- `reports/clinical_validation/clinical_validation_bootstrap_distributions.pdf`
- `reports/clinical_validation/clinical_validation_bootstrap_distributions.svg`
- `reports/figures_publication/figureS_clinical_validation_bootstrap_distributions.png`
- `reports/figures_publication/figureS_clinical_validation_bootstrap_distributions.pdf`
- `reports/figures_publication/figureS_clinical_validation_bootstrap_distributions.svg`

QC result:

- Includes both histogram and violin views.
- Shows both n=24 primary and n=29 proxy sensitivity distributions.
- Zero reference line is visible.
- Axes are labeled for Spearman rho and pairwise C-index.
- Legend and tick labels are readable.
- Distribution shape is visually interpretable; the n=29 distribution shifts
  modestly positive relative to n=24, consistent with the tabular sensitivity
  analysis.

Recommended use:

- Supplementary Figure. This is methodologically useful but visually dense for
  a main-text figure unless the target journal allows multi-panel statistical
  supplements in the main article.

### Score vs PFS Scatter

Files:

- `reports/clinical_validation/clinical_validation_score_vs_pfs.png`
- `reports/clinical_validation/clinical_validation_score_vs_pfs.pdf`
- `reports/clinical_validation/clinical_validation_score_vs_pfs.svg`
- `reports/figures_publication/figure4_local_validation_score_vs_pfs.png`
- `reports/figures_publication/figure4_local_validation_score_vs_pfs.pdf`
- `reports/figures_publication/figure4_local_validation_score_vs_pfs.svg`

QC result:

- Scatter plot is readable and useful for visualizing score-PFS association.
- Patient labels are visible; some label density remains but does not obscure
  the main trend.
- Recommended as a main-text local validation panel or as a companion to KM.

### PFS by Rank

Files:

- `reports/clinical_validation/clinical_validation_pfs_by_rank.png`
- `reports/clinical_validation/clinical_validation_pfs_by_rank.pdf`
- `reports/clinical_validation/clinical_validation_pfs_by_rank.svg`
- `reports/figures_publication/figureS_local_validation_pfs_by_rank.png`
- `reports/figures_publication/figureS_local_validation_pfs_by_rank.pdf`
- `reports/figures_publication/figureS_local_validation_pfs_by_rank.svg`

QC result:

- This is not a KM figure and should not be described as survival analysis.
- It is a descriptive box/scatter comparison of observed PFS1 by rank of the
  administered TKI.
- Recommended only as a supplementary descriptive figure.

## Subgroup Reporting Guardrail

The subgroup interpretation table is saved at:

- `reports/clinical_validation/clinical_validation_subgroup_interpretation.csv`

Key reporting constraints:

- classical-like subgroup: n = 8, descriptive only.
- model-native compound subgroup: n = 9, descriptive only.
- PACC subgroup: n = 16 in primary analysis and n = 21 in proxy sensitivity
  analysis; hypothesis-generating, not confirmatory.
- proxy-sensitivity compound subgroup: n = 13; hypothesis-generating, not
  confirmatory.

Recommended language:

> Subgroup findings were underpowered and are interpreted as
> hypothesis-generating. The directionally positive PACC and compound-mutation
> signals align with the independent literature leave-one-out result, but do
> not constitute definitive interaction evidence.

## Final QC Verdict

The KM and bootstrap distribution figures are publication-usable after the
current formatting pass. The preferred figure hierarchy is:

1. Main text: score-vs-PFS scatter plus KM curve.
2. Supplement: bootstrap distributions and PFS-by-rank descriptive plot.
