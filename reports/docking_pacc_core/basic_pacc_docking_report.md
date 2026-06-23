# Basic EGFR-PACC Docking Report

## Scope

This is a coarse AutoDock Vina 1.2.7 docking screen intended as a CSBJ-oriented supplementary structural analysis.
It is not a covalent docking, MD, FEP, or binding free-energy calculation.

## Inputs

- Receptor template: 4LRM chain A.
- Docking box reference ligand: YUN.
- Box center: [49.297, 360.413, 31.953].
- Box size: [24.0, 24.0, 24.0] Angstrom.
- Receptor variants: WT, G719S, L861Q, S768I.
- Ligands: afatinib, osimertinib, furmonertinib.
- Receptor point mutants were generated with PDBFixer and were not MD-relaxed.
- Receptor PDBQT files were generated with OpenBabel. This route avoids the RDKit residue-valence parsing failure observed with Meeko receptor preparation on the current Windows/Python runtime.
- Ligands were desalted, embedded in 3D with RDKit, minimized, and converted to PDBQT with Meeko.

## Score Summary

| variant_model | drug | vina_score_kcal_mol | rank_within_variant | receptor_template | mutation_model | pose_pdbqt | log_file | method_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| G719S | osimertinib | -9.43 | 1.0 | 4LRM_chainA | GLY-719-SER | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\G719S__osimertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\G719S__osimertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| G719S | furmonertinib | -9.141 | 2.0 | 4LRM_chainA | GLY-719-SER | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\G719S__furmonertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\G719S__furmonertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| G719S | afatinib | -9.001 | 3.0 | 4LRM_chainA | GLY-719-SER | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\G719S__afatinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\G719S__afatinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| L861Q | furmonertinib | -9.591 | 1.0 | 4LRM_chainA | LEU-861-GLN | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\L861Q__furmonertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\L861Q__furmonertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| L861Q | afatinib | -9.071 | 2.0 | 4LRM_chainA | LEU-861-GLN | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\L861Q__afatinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\L861Q__afatinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| L861Q | osimertinib | -8.89 | 3.0 | 4LRM_chainA | LEU-861-GLN | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\L861Q__osimertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\L861Q__osimertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| S768I | furmonertinib | -9.573 | 1.0 | 4LRM_chainA | SER-768-ILE | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\S768I__furmonertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\S768I__furmonertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| S768I | afatinib | -9.358 | 2.0 | 4LRM_chainA | SER-768-ILE | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\S768I__afatinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\S768I__afatinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| S768I | osimertinib | -9.134 | 3.0 | 4LRM_chainA | SER-768-ILE | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\S768I__osimertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\S768I__osimertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| WT | osimertinib | -9.684 | 1.0 | 4LRM_chainA | WT | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\WT__osimertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\WT__osimertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| WT | furmonertinib | -9.382 | 2.0 | 4LRM_chainA | WT | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\WT__furmonertinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\WT__furmonertinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |
| WT | afatinib | -9.348 | 3.0 | 4LRM_chainA | WT | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\poses\WT__afatinib__vina_pose.pdbqt | F:\egfr_pacc_digital_cell\reports\docking_pacc_core\logs\WT__afatinib__vina.log | Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel. |

## Interpretation Boundary

- EGFR TKIs such as afatinib, osimertinib, and furmonertinib can form covalent interactions; this Vina run models only non-covalent poses.
- PDBFixer-generated point-mutant side chains are not MD-relaxed.
- Receptor PDBQT files were prepared with OpenBabel; future higher-stringency studies should compare against alternative receptor preparation and protonation workflows.
- Scores should be interpreted as coarse pose/compatibility evidence and used only as supplementary support for the structure-informed ranking framework.
- The manuscript should report these values as approximate Vina docking scores, not as experimental affinity or DeltaG_bind.

## Outputs

- `docking_scores.csv`
- `basic_pacc_docking_heatmap.png/pdf/svg`
- `basic_pacc_docking_pose_overview.png/pdf/svg`
- `poses/*.pdbqt`
- `receptors/*.pdb` and `receptors/*.pdbqt`
- `ligands/*.sdf` and `ligands/*.pdbqt`
