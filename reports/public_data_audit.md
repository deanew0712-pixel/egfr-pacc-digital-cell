# EGFR-PACC public data audit

## P0/P1/P2 download status

- CIViC EGFR evidence rows: 141; accepted predictive rows: 80.
- cBioPortal studies summarized: 4; EGFR-mutant sample rows: 119.
- GDC TCGA LUAD/LUSC masked somatic mutation file records: 1167.
- GDC TCGA LUAD/LUSC MAF files downloaded/status rows: 1167.
- RCSB PDB structures downloaded/status rows: 5.
- AlphaFold EGFR P00533 files downloaded/status rows: 4.
- ChEMBL ligand rows: 4.

## cBioPortal EGFR co-mutation prior candidates

- luad_tcga_pub: EGFR-mut n=33, PACC-like n=4, TP53=0.6061, PIK3CA=0.0303, MET_amp=0.0303, RB1=0.0909
- luad_tcga_pan_can_atlas_2018: EGFR-mut n=65, PACC-like n=6, TP53=0.6615, PIK3CA=0.0615, MET_amp=0.0154, RB1=0.1077
- lusc_tcga_pub: EGFR-mut n=7, PACC-like n=0, TP53=0.8571, PIK3CA=0.1429, MET_amp=0.0, RB1=0.1429
- lusc_tcga_pan_can_atlas_2018: EGFR-mut n=14, PACC-like n=0, TP53=0.9286, PIK3CA=0.0714, MET_amp=0.0, RB1=0.0714

## Restricted/manual sources not yet downloaded

- OncoKB: needs `ONCOKB_TOKEN`; use after token approval for Baseline-1 actionability.
- COSMIC: optional enhanced dependency; open fallback uses cBioPortal/TCGA/GDC/ClinVar/CIViC/Cancer Hotspots.
- DrugBank: optional enhanced dependency; open fallback uses ChEMBL/PubChem/DailyMed/openFDA/manual literature seed.
- GENIE: optional enhanced dependency; open fallback uses cBioPortal/TCGA/GDC.
- Literature ORR table: must be manually curated and locked before opening the 37-case outcomes.

## Immediate interpretation

The open fallback layer is sufficient to run the MVP and generate auditable prior candidates. Final freeze still requires manual ORR/PFS curation; OncoKB/COSMIC/GENIE/DrugBank can enhance but should not block the MVP.
