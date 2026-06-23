# Publication Package Phase 3 Summary

Generated: 2026-06-18 Asia/Singapore

## Scope

This package converts the current EGFR-PACC digital-cell MVP into a manuscript-ready analysis bundle:

1. Publication-quality Figure 1, Figure 2, and Figure 3 in PDF, SVG, and 300 DPI PNG.
2. A Methods and Results working draft based on the frozen prior, formal PyMC run, demo case ranking, and leave-one-out back-testing.
3. ACHILLES/TORG1834 and UNICORN source-PDF verification for the locked ORR rows used in back-testing.

## Generated Figures

Figure source script:

- `scripts/make_publication_figures.py`

Outputs:

- `reports/figures_publication/figure1_framework.pdf`
- `reports/figures_publication/figure1_framework.svg`
- `reports/figures_publication/figure1_framework.png`
- `reports/figures_publication/figure2_demo_case_ranking.pdf`
- `reports/figures_publication/figure2_demo_case_ranking.svg`
- `reports/figures_publication/figure2_demo_case_ranking.png`
- `reports/figures_publication/figure3_backtesting.pdf`
- `reports/figures_publication/figure3_backtesting.svg`
- `reports/figures_publication/figure3_backtesting.png`

Figure content:

- Figure 1: EGFR-PACC digital-cell framework, from locked public priors to mutation classification, PyMC partial pooling, EPDRI ranking, back-testing, and future clinical validation.
- Figure 2: demo case drug ranking for EGFR G719S plus L861Q with CNS disease and TP53 co-mutation, using the formal four-chain PyMC posterior prior table.
- Figure 3: leave-one-out back-testing, including AUC with bootstrap CI, delta AUC versus baselines, ablation comparison, and subgroup target-rank analysis.

Graphics settings:

- PDF/SVG vector outputs.
- PNG output at 300 DPI.
- Matplotlib font settings use TrueType-compatible PDF embedding and text-preserving SVG output.

## Manuscript Draft

Draft file:

- `docs/manuscript_methods_results_draft_2026-06-18.md`

Status:

- Working draft, approximately 2900 words.
- Includes Methods and Results sections with placeholders for Figure 1, Figure 2, and Figure 3.
- Uses current locked outputs only; no local clinical outcome data are used.

Core statements supported by current data:

- The framework is an auditable, prior-driven EGFR uncommon/PACC ranking model.
- The formal PyMC prior run completed with 4 chains, 2000 tuning iterations, 2000 posterior draws per chain, target acceptance 0.995, 0 divergences, maximum R-hat 1.0017, minimum bulk ESS 1567.5, and minimum tail ESS 1665.8.
- Leave-one-out back-testing across 29 ORR-evaluable literature rows showed model AUC 0.8787, Baseline-0 AUC 0.8033, and Baseline-1 AUC 0.8602.
- Bootstrap CI supports improvement versus Robichaux-only Baseline-0: mean delta AUC +0.0740, 95% CI +0.0030 to +0.1607.
- Delta AUC versus OncoKB-enhanced Baseline-1 crosses zero: mean +0.0186, 95% CI -0.0269 to +0.0789, so this comparison should be described as exploratory.
- Top-1 hit rate is not improved over baselines; current claim should focus on overall discrimination rather than superior top-1 recommendation.

## ACHILLES and UNICORN Verification

Verification report:

- `reports/literature_lock_achilles_unicorn_check_2026-06-18.md`

Local source PDFs/text:

- `data/raw/literature_pdfs/pragmatic_afatinib_vs_chemotherapy_nsclc.pdf`
- `data/raw/literature_pdfs/unicorn_osimertinib_uncommon_egfr.pdf`
- `reports/literature_pdf_text/pragmatic_afatinib_vs_chemotherapy_nsclc.txt`
- `reports/literature_pdf_text/unicorn_osimertinib_uncommon_egfr.txt`

Conclusion:

- Locked ORR values used for leave-one-out back-testing are supported by local ACHILLES/TORG1834 and UNICORN PDF text extraction.
- UNICORN uncommon/uncommon ORR should use 54.5% from the response table; one OCR passage reads 54.4%.
- ACHILLES subgroup PFS values are adequate for manuscript background but should be rechecked against supplement/source tables before being used as primary PFS endpoints.

## Reproducibility Commands

From project root:

```bash
PYTHONPATH=src .venv/bin/python scripts/make_publication_figures.py
PYTHONPATH=src .venv/bin/python -m py_compile scripts/make_publication_figures.py scripts/run_backtest_loo.py src/egfr_pacc/classifier.py src/egfr_pacc/bayes_model.py src/egfr_pacc/epdri.py src/egfr_pacc/ranking.py
```

## Immediate Next Steps

1. Human review of Figure 1-3 layout and labels before journal-specific formatting.
2. Decide whether to keep the manuscript draft in English or convert to Chinese outline first.
3. Obtain or verify ACHILLES supplementary tables if PFS subgroup values will be treated as primary locked outcomes.
4. Freeze the figure/table numbering before integrating the 30-patient local clinical validation cohort.
