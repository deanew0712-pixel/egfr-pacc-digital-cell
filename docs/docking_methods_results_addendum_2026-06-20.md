# Basic Docking Addendum for CSBJ-Oriented Submission

Generated: 2026-06-20 Asia/Shanghai

## Rationale

To strengthen the structural-computational component of the EGFR-PACC digital
cell manuscript, a basic AutoDock Vina docking screen was added for core EGFR
uncommon/PACC-relevant mutations and representative EGFR TKIs. This analysis is
intended as supplementary structural support rather than a high-precision
binding free-energy calculation.

## Methods

The EGFR kinase domain structure 4LRM chain A was used as the receptor template.
The docking grid was centered on the co-crystal ligand YUN in the ATP-binding
pocket, with a 24 A cubic search space centered at x = 49.297, y = 360.413,
z = 31.953. Receptor models were generated for WT, G719S, L861Q, and S768I.
Point-mutant side chains were generated with PDBFixer and were not MD-relaxed.

Afatinib, osimertinib, and furmonertinib/firmonertinib were prepared from
ChEMBL canonical SMILES. Salt fragments were removed by retaining the largest
molecular fragment. Three-dimensional conformers were generated with RDKit
ETKDG and force-field minimized. Ligand PDBQT files were generated with Meeko.

Docking was performed with AutoDock Vina 1.2.7 using exhaustiveness 16,
9 output modes, and a fixed random seed. Receptor PDBQT files were generated
with OpenBabel. This route was selected because Meeko receptor preparation
triggered RDKit residue-valence parsing failures for the 4LRM-derived receptor
on the current Windows/Python runtime. OpenBabel converted all WT and mutant
receptor PDB files successfully, and conversion logs are retained in
`reports/docking_pacc_core/receptors/*.obabel_receptor.log`.

## Results

The complete docking score matrix is saved in
`reports/docking_pacc_core/docking_scores.csv`. Vina scores ranged from
approximately -8.98 to -9.73 kcal/mol. Across the core mutation models,
osimertinib or furmonertinib showed the most favorable Vina scores, while
afatinib remained within the same broad affinity range but ranked lower in this
coarse non-covalent docking setup.

Best-scoring drugs by receptor model after OpenBabel receptor preparation were:

- WT: osimertinib, -9.68 kcal/mol.
- G719S: osimertinib, -9.43 kcal/mol.
- L861Q: furmonertinib, -9.59 kcal/mol.
- S768I: furmonertinib, -9.57 kcal/mol.

The heatmap and representative pose overview are saved as:

- `reports/docking_pacc_core/basic_pacc_docking_heatmap.png`
- `reports/docking_pacc_core/basic_pacc_docking_pose_overview.png`

Publication copies are saved under `reports/figures_publication/`.

## Interpretation Boundary

This docking analysis should be reported cautiously. Afatinib, osimertinib, and
furmonertinib can act as covalent EGFR inhibitors, but this analysis models only
non-covalent binding poses. The receptor point-mutant models were not relaxed by
MD, and the receptor preparation/protonation workflow was not benchmarked
against multiple alternatives. The scores should therefore be interpreted as
approximate Vina pose/compatibility signals, not as experimental affinity,
covalent binding energetics, or formal DeltaG_bind.

Recommended manuscript language:

> A supplementary coarse AutoDock Vina screen was performed to provide a
> structure-oriented sanity check of ATP-pocket compatibility across core EGFR
> uncommon/PACC mutation models. The analysis supported plausible pocket
> accommodation of the tested EGFR TKIs but was not used to refit EPDRI priors
> or override the locked clinical evidence layer.
