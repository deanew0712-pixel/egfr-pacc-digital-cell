# Receptor Workflow Sensitivity and Covalent Geometry Screen

## Scope

This analysis extends the basic EGFR-PACC docking screen with receptor preparation sensitivity and a pre-covalent warhead-geometry screen.
It is intended as supplementary structural evidence, not as a formal covalent binding-energy calculation.

## Receptor Workflow Sensitivity

AutoDock Vina was rerun with exhaustiveness 8 across 3 receptor-preparation workflows:

- `obabel_default`: OpenBabel PDBQT conversion from heavy-atom PDB; default OpenBabel charges.
- `obabel_gasteiger`: OpenBabel PDBQT conversion with Gasteiger partial charges.
- `pdbfixer_ph7p4_obabel_gasteiger`: PDBFixer protonation at pH 7.4, nonpolar hydrogens removed, followed by OpenBabel Gasteiger PDBQT.

Top-rank counts across 3 workflows:

| variant_model | drug | top_rank_count_across_3_workflows |
| --- | --- | --- |
| G719S | osimertinib | 3 |
| L861Q | furmonertinib | 2 |
| L861Q | osimertinib | 1 |
| S768I | furmonertinib | 3 |
| WT | osimertinib | 2 |
| WT | afatinib | 1 |

Score-range summary:

| variant_model | drug | mean_score | sd_score | min_score | max_score | median_rank | best_rank_count | score_range_kcal_mol |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| G719S | afatinib | -9.030 | 0.050 | -9.087 | -9.001 | 3.000 | 0 | 0.086 |
| G719S | furmonertinib | -9.268 | 0.105 | -9.389 | -9.207 | 2.000 | 0 | 0.182 |
| G719S | osimertinib | -9.452 | 0.295 | -9.793 | -9.282 | 1.000 | 3 | 0.511 |
| L861Q | afatinib | -9.139 | 0.118 | -9.275 | -9.071 | 2.000 | 0 | 0.204 |
| L861Q | furmonertinib | -9.649 | 0.100 | -9.764 | -9.591 | 1.000 | 2 | 0.173 |
| L861Q | osimertinib | -9.195 | 0.528 | -9.804 | -8.890 | 3.000 | 1 | 0.914 |
| S768I | afatinib | -9.269 | 0.155 | -9.358 | -9.090 | 2.000 | 0 | 0.268 |
| S768I | furmonertinib | -9.493 | 0.218 | -9.745 | -9.367 | 1.000 | 3 | 0.378 |
| S768I | osimertinib | -9.331 | 0.341 | -9.724 | -9.134 | 3.000 | 0 | 0.590 |
| WT | afatinib | -9.352 | 0.007 | -9.360 | -9.348 | 2.000 | 1 | 0.012 |
| WT | furmonertinib | -9.179 | 0.079 | -9.224 | -9.088 | 3.000 | 0 | 0.136 |
| WT | osimertinib | -9.551 | 0.279 | -9.712 | -9.229 | 1.000 | 2 | 0.483 |

## Covalent Geometry Screen

The acrylamide warhead was detected from each ligand SMILES using the `C=CC(=O)N` SMARTS pattern.
For each Vina pose, the distance from the ATP-pocket cysteine sulfur to the acrylamide beta-carbon and the SG-betaC-alphaC angle were calculated.
The ATP-pocket cysteine in the 4LRM-derived template is reported by the structure as CYS A800; this is the local template residue used for the geometry screen.

Best geometry-ranked poses by receptor-drug pair:

| variant_model | drug | pose_model | vina_score_kcal_mol | template_covalent_cys | cys_sg_to_acrylamide_beta_c_distance_A | sg_beta_alpha_angle_deg | near_attack_pose |
| --- | --- | --- | --- | --- | --- | --- | --- |
| G719S | afatinib | 5 | -8.367 | A:CYS800:SG | 3.759 | 68.518 | False |
| G719S | furmonertinib | 8 | -8.237 | A:CYS800:SG | 4.707 | 78.190 | True |
| G719S | osimertinib | 8 | -8.734 | A:CYS800:SG | 4.769 | 80.702 | True |
| L861Q | afatinib | 2 | -8.913 | A:CYS800:SG | 3.722 | 80.719 | True |
| L861Q | furmonertinib | 6 | -8.940 | A:CYS800:SG | 3.791 | 65.611 | False |
| L861Q | osimertinib | 9 | -8.066 | A:CYS800:SG | 4.163 | 94.740 | True |
| S768I | afatinib | 8 | -8.108 | A:CYS800:SG | 3.792 | 82.835 | True |
| S768I | furmonertinib | 2 | -9.569 | A:CYS800:SG | 4.019 | 81.723 | True |
| S768I | osimertinib | 8 | -8.481 | A:CYS800:SG | 3.821 | 96.858 | True |
| WT | afatinib | 3 | -8.881 | A:CYS800:SG | 3.809 | 80.812 | True |
| WT | furmonertinib | 7 | -8.445 | A:CYS800:SG | 3.657 | 73.280 | True |
| WT | osimertinib | 5 | -8.999 | A:CYS800:SG | 4.328 | 74.672 | True |

## Interpretation Boundary

- Receptor protonation/workflow sensitivity is a robustness check for the non-covalent Vina score ranking.
- The covalent analysis is a pre-reaction geometry screen, not a tethered covalent AutoDock4 or FEP calculation.
- A strict covalent docking upgrade should use a validated flexible-side-chain/tethered protocol and should benchmark reaction-state ligand preparation.
