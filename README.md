# EGFR-PACC Digital Cell Publication Package

This repository contains the public reproducibility package for the EGFR-PACC
digital-cell framework described in the CSBJ manuscript:

**Structure-Informed Digital-Cell Ranking of EGFR-PACC Inhibitor Sensitivity**

The package includes source code, locked configuration files, open/derived
reference tables, literature back-testing outputs, aggregate clinical-validation
outputs, docking summaries, and publication figure files. It is intended for
research reproducibility, review, and methods inspection. It is not a clinical
decision system.

## What Is Included

- `src/egfr_pacc/`: mutation classification, Bayesian prior, EPDRI ranking, and
  report-generation modules.
- `scripts/`: entry points for demo cases, literature leave-one-out
  back-testing, local validation reporting, publication figures, basic docking,
  docking sensitivity, and OpenMM-relaxed redocking workflows.
- `configs/`: locked mutation-class rules, priors, baselines, drug properties,
  and fallback model configuration.
- `data/processed/`: open or derived non-sensitive tables used by the framework.
- `data/clinical_holdout/*_template.csv`: templates for local clinical cohort
  inputs. Patient-level clinical data are not redistributed here.
- `reports/backtest_loo/`: leave-one-out literature back-testing outputs.
- `reports/clinical_validation/`: aggregate clinical-validation metrics and
  summaries. Patient-level prediction tables are not redistributed here.
- `reports/docking_*`: docking score summaries, workflow sensitivity summaries,
  and method reports. Large pose/receptor coordinate files are not included.
- `reports/figures_publication/` and
  `reports/csbj_submission_figures_2026-06-21/`: editable SVG/PDF figure files.

## What Is Excluded

The following materials are intentionally excluded from this public GitHub
package:

- Patient-level clinical feature, outcome, prediction, and treatment-ranking
  tables. The de-identified clinical cohort is deposited separately under OMIX
  accession `OMIX017364`, subject to the access conditions of that repository.
- Credentialed or restricted third-party resources, including raw COSMIC,
  OncoKB, DrugBank, and any API-token-dependent downloads.
- Copyrighted source PDFs and conference slide/poster files.
- Large raw structure files, docking poses, receptor coordinate intermediates,
  generated zip archives, local virtual environments, and binary executables.

Scripts and source-status reports are retained where possible so that users with
appropriate third-party access can regenerate restricted components.

## Quick Start

Create a Python environment and install the core dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run a demonstration case:

```bash
PYTHONPATH=src python scripts/run_demo_case.py
```

Run literature leave-one-out back-testing:

```bash
PYTHONPATH=src python scripts/run_backtest_loo.py
```

Regenerate publication figures from available outputs:

```bash
PYTHONPATH=src python scripts/make_publication_figures.py
PYTHONPATH=src python scripts/make_csbj_submission_figures.py
```

Clinical validation requires the de-identified cohort files. Place approved
copies under `data/clinical_holdout/` using the provided templates, then run:

```bash
PYTHONPATH=src python scripts/run_clinical_validation.py --bootstrap 5000
```

Docking workflows require additional docking dependencies and a local Vina
binary. See `requirements-docking.txt` and the docking reports under
`reports/docking_pacc_core/` and `reports/docking_workflow_sensitivity/`.

## Reproducibility Notes

- The model is a frozen, literature-prior framework. Local clinical outcomes are
  not used to construct, select, or tune EPDRI.
- OncoKB-aware baseline regeneration requires user-provided OncoKB access.
  OncoKB-derived raw files are not redistributed.
- COSMIC-enhanced catalog regeneration requires licensed COSMIC files. The open
  fallback catalog is included.
- Docking outputs are rank-order structural plausibility analyses, not binding
  free-energy estimates.

## Suggested Citation

Please cite the associated manuscript when available. If this package is
archived on Zenodo, cite the versioned Zenodo DOI for exact reproducibility.

