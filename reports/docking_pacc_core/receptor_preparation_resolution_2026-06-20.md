# Receptor Preparation Resolution

## Problem

Meeko receptor preparation failed on the current Windows/Python runtime when parsing residue valence states for the 4LRM-derived EGFR receptor models. The observed failure was consistent with RDKit residue/bond perception producing invalid explicit valence assignments during receptor preparation.

## Resolution

The docking pipeline now prepares receptor PDBQT files with OpenBabel first. Meeko is retained for ligand PDBQT generation, where it worked reliably for afatinib, osimertinib, and furmonertinib.

Implementation:

- Script: `scripts/run_basic_pacc_docking.py`
- Primary receptor route: OpenBabel PDB to PDBQT conversion with rigid receptor output
- Fallback route: Meeko receptor preparation, retained only if OpenBabel is unavailable
- Ligand route: RDKit 3D embedding/minimization followed by Meeko PDBQT generation

## Verification

OpenBabel converted all receptor models successfully:

| receptor model | conversion log | status |
| --- | --- | --- |
| WT | `EGFR_4LRM_chainA_WT.obabel_receptor.log` | 1 molecule converted |
| G719S | `EGFR_4LRM_chainA_G719S.obabel_receptor.log` | 1 molecule converted |
| L861Q | `EGFR_4LRM_chainA_L861Q.obabel_receptor.log` | 1 molecule converted |
| S768I | `EGFR_4LRM_chainA_S768I.obabel_receptor.log` | 1 molecule converted |

AutoDock Vina accepted the OpenBabel-prepared receptor PDBQT files and completed all 12 receptor-ligand runs:

- 4 receptor models: WT, G719S, L861Q, S768I
- 3 ligands: afatinib, osimertinib, furmonertinib
- Output table: `reports/docking_pacc_core/docking_scores.csv`
- Pose files: `reports/docking_pacc_core/poses/*.pdbqt`

No Meeko receptor-preparation log files were generated in the final run, confirming that the OpenBabel receptor route was used.

## Boundary For Future Work

This fix resolves the receptor-preparation runtime failure and makes the basic docking workflow reproducible for supplementary structural screening. It does not replace higher-stringency receptor preparation, protonation benchmarking, covalent docking, MD relaxation, or FEP/MM-GBSA binding-energy estimation.
