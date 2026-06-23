# Docking Workflow Sensitivity and Covalent Geometry Addendum

## Rationale

To strengthen the structural supplement for the EGFR-PACC digital-cell manuscript, the basic Vina docking analysis was extended in two directions:

1. receptor preparation/protonation workflow sensitivity; and
2. acrylamide warhead geometry screening against the ATP-pocket cysteine.

This analysis is designed as a robustness and pose-plausibility supplement. It should not be described as FEP, MD, experimental affinity, or formal covalent binding-energy estimation.

## Receptor Workflow Sensitivity

AutoDock Vina was rerun for WT, G719S, L861Q, and S768I receptor models against afatinib, osimertinib, and furmonertinib using three receptor-preparation workflows:

- OpenBabel default receptor PDBQT conversion.
- OpenBabel receptor PDBQT conversion with Gasteiger partial charges.
- PDBFixer pH 7.4 protonation, nonpolar hydrogen removal, and OpenBabel Gasteiger receptor PDBQT conversion.

Across these workflows, the strongest median-rank patterns remained directionally consistent:

- G719S: osimertinib ranked first in all 3 workflows.
- L861Q: furmonertinib ranked first in 2 of 3 workflows.
- S768I: furmonertinib ranked first in all 3 workflows.
- WT: osimertinib ranked first in 2 of 3 workflows.

The largest workflow-driven score ranges were observed for L861Q-osimertinib and S768I-osimertinib, indicating that score magnitudes are sensitive to receptor preparation. The rank-level conclusion is therefore more appropriate than overinterpreting absolute Vina score differences.

Primary outputs:

- `reports/docking_workflow_sensitivity/receptor_workflow_sensitivity_scores.csv`
- `reports/docking_workflow_sensitivity/receptor_workflow_sensitivity_summary.csv`
- `reports/docking_workflow_sensitivity/docking_workflow_and_covalent_sensitivity.png/pdf/svg`

## Covalent Geometry Screen

Afatinib, osimertinib, and furmonertinib contain acrylamide warheads. For each ligand, the acrylamide motif was detected using the SMARTS pattern `C=CC(=O)N`. For each Vina pose, two pre-reaction geometric descriptors were calculated:

- distance from the ATP-pocket cysteine sulfur to the acrylamide beta-carbon; and
- SG-betaC-alphaC angle.

The ATP-pocket cysteine in the 4LRM-derived local template is CYS A800. This local template residue was used for the covalent geometry screen.

Most receptor-drug pairs contained at least one near-attack pose under the predefined loose geometric screen. Two exceptions were flagged:

- G719S-afatinib: best geometry-ranked pose had distance 3.76 A but angle 68.5 degrees, just outside the 70-degree lower bound.
- L861Q-furmonertinib: best geometry-ranked pose had distance 3.79 A but angle 65.6 degrees, just outside the 70-degree lower bound.

These exceptions should be interpreted conservatively. The geometry screen supports warhead accessibility and pose plausibility, but it does not quantify covalent reaction energetics.

Primary outputs:

- `reports/docking_workflow_sensitivity/covalent_warhead_geometry_screen.csv`
- `reports/docking_workflow_sensitivity/covalent_warhead_geometry_best_poses.csv`
- `reports/docking_workflow_sensitivity/workflow_sensitivity_and_covalent_screen_report.md`

## Manuscript Wording

Suggested concise wording:

> Receptor-preparation sensitivity analysis across three OpenBabel/PDBFixer workflows showed directionally stable rank patterns for key PACC models, with osimertinib favored in G719S and furmonertinib favored in L861Q and S768I. Because EGFR TKIs in this panel contain acrylamide warheads, we further screened Vina poses for pre-covalent geometry against the ATP-pocket cysteine in the 4LRM-derived template. This analysis supported warhead accessibility for most receptor-drug pairs, but was treated as pose-plausibility evidence rather than formal covalent binding free-energy estimation.

## Boundary

For a stricter future structural study, the next upgrade should use a validated tethered/flexible-side-chain covalent docking protocol and benchmark reaction-state ligand preparation, receptor protonation alternatives, and MD-relaxed mutant receptor conformations.
