# Excluded Materials

This public package excludes materials that should not be redistributed through
GitHub:

- Patient-level clinical cohort files:
  - `data/clinical_holdout/clinical_features_31.csv`
  - `data/clinical_holdout/clinical_outcomes_31.csv`
  - patient-level prediction/ranking tables under `reports/clinical_validation/`
- Raw restricted or credentialed data:
  - COSMIC raw downloads
  - OncoKB raw API outputs
  - DrugBank exports, if available
- Copyrighted source PDFs and conference slide/poster files.
- Large generated binary artifacts:
  - TIFF/PNG figure exports
  - docking pose/receptor coordinate intermediates
  - OpenMM relaxed output zip archives
  - local Vina executables

Where possible, aggregate outputs, open fallback tables, scripts, and
regeneration instructions are included instead.

