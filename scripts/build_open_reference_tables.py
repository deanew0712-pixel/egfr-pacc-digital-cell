from __future__ import annotations

import csv
import gzip
import json
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import yaml


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
REFERENCE = PROCESSED / "reference"
OPEN_RAW = RAW / "open_fallback"

USER_AGENT = "egfr-pacc-digital-cell/0.1"
GENERATED_AT = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

KEY_MUTATIONS = ["G719X", "S768I", "E709X", "L861Q", "C797S", "L718Q", "ex20ins", "compound_EGFR"]
PACC_PREFIXES = ("G719", "E709")
PACC_EXACT = {"S768I", "L747P", "G719A", "G719C", "G719D", "G719S", "E709A", "E709G", "E709K"}
BYPASS_GENES = {"TP53", "RB1"}

DRUGS = [
    "gefitinib",
    "erlotinib",
    "afatinib",
    "dacomitinib",
    "osimertinib",
    "lazertinib",
    "furmonertinib",
    "aumolertinib",
    "mobocertinib",
    "amivantamab",
    "patritumab deruxtecan",
]

LICENSE_NOTES = {
    "cBioPortal": "cBioPortal public API; source cohort licenses vary, review original study terms.",
    "TCGA/GDC": "GDC open-access masked somatic mutation MAF files; use per GDC data policies.",
    "CIViC": "CIViC open variant evidence; cite CIViC and source publications.",
    "ClinVar": "NCBI ClinVar public data via E-utilities.",
    "Cancer Hotspots": "Cancer Hotspots public download; cite Chang et al./cancerhotspots.org.",
    "ChEMBL": "ChEMBL public data under EMBL-EBI terms.",
    "PubChem": "PubChem public data from NCBI.",
    "DailyMed": "DailyMed public SPL labeling from NLM.",
    "openFDA": "openFDA public API; data are unvalidated and subject to openFDA terms.",
    "manual_literature_seed": "Manual literature curation seed; must be reviewed before prior freeze.",
    "optional_enhanced": "Optional enhanced source; absence must not block open pipeline.",
}


def warn(message: str) -> None:
    print(f"WARNING: {message}")


def request_json(url: str) -> dict | list | None:
    try:
        req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        with urlopen(req, timeout=45) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as exc:
        warn(f"Could not fetch {url}: {exc}")
        return None


def request_bytes(url: str) -> bytes | None:
    try:
        req = Request(url, headers={"User-Agent": USER_AGENT})
        with urlopen(req, timeout=90) as response:
            return response.read()
    except Exception as exc:
        warn(f"Could not fetch {url}: {exc}")
        return None


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def read_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_protein_change(value: str | None) -> str:
    if not value:
        return ""
    value = value.strip()
    value = re.sub(r"^p\.", "", value)
    value = re.sub(r"^[A-Z][a-z]{2}", lambda m: three_to_one(m.group(0)), value)
    return value


def three_to_one(token: str) -> str:
    mapping = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    }
    return mapping.get(token, token)


def mutation_bucket(protein_change: str, variant_type: str = "", variant_classification: str = "") -> str | None:
    pc = normalize_protein_change(protein_change)
    if not pc:
        return None
    if pc.startswith("G719"):
        return "G719X"
    if pc.startswith("E709"):
        return "E709X"
    if pc in {"S768I", "L861Q", "C797S", "L718Q"}:
        return pc
    text = f"{pc} {variant_type} {variant_classification}".lower()
    if "ins" in text and ("exon 20" in text or "ex20" in text or re.search(r"[a-z]\d+_[a-z]\d+ins", text)):
        return "ex20ins"
    return None


def is_pacc_variant(protein_change: str) -> bool:
    pc = normalize_protein_change(protein_change)
    return pc in PACC_EXACT or pc.startswith(PACC_PREFIXES)


def ensure_optional_dirs() -> None:
    for dirname in ["cosmic", "genie", "drugbank"]:
        path = RAW / dirname
        path.mkdir(parents=True, exist_ok=True)
        if not any(path.iterdir()):
            warn(f"Optional enhanced dependency not found or empty: {path}")


def download_cancer_hotspots() -> Path | None:
    out = OPEN_RAW / "cancer_hotspots" / "hotspots_v2.xls"
    if out.exists() and out.stat().st_size > 0:
        return out
    data = request_bytes("https://www.cancerhotspots.org/files/hotspots_v2.xls")
    if data:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        return out
    return None


def fetch_clinvar() -> dict[str, dict]:
    out_dir = OPEN_RAW / "clinvar"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    for mutation in KEY_MUTATIONS:
        if mutation == "compound_EGFR":
            continue
        query = mutation.replace("X", "")
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            + urlencode({"db": "clinvar", "term": f"EGFR[gene] AND {query}", "retmode": "json", "retmax": "20"})
        )
        data = request_json(url) or {}
        write_json(out_dir / f"{mutation}.json", data)
        results[mutation] = data
        time.sleep(0.1)
    return results


def fetch_pubchem() -> dict[str, dict]:
    out_dir = OPEN_RAW / "pubchem"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    props = "CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight,XLogP,TPSA"
    for drug in DRUGS:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(drug)}/property/{props}/JSON"
        data = request_json(url) or {}
        write_json(out_dir / f"{drug.replace(' ', '_')}.json", data)
        results[drug] = data
        time.sleep(0.1)
    return results


def fetch_dailymed_openfda() -> tuple[dict[str, dict], dict[str, dict]]:
    dailymed_dir = OPEN_RAW / "dailymed"
    openfda_dir = OPEN_RAW / "openfda"
    dailymed_dir.mkdir(parents=True, exist_ok=True)
    openfda_dir.mkdir(parents=True, exist_ok=True)
    dailymed_results = {}
    openfda_results = {}
    for drug in DRUGS:
        dm_url = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?" + urlencode(
            {"drug_name": drug, "page": 1, "pageSize": 3}
        )
        dm_data = request_json(dm_url) or {}
        write_json(dailymed_dir / f"{drug.replace(' ', '_')}.json", dm_data)
        dailymed_results[drug] = dm_data

        fda_url = "https://api.fda.gov/drug/label.json?" + urlencode({"search": f"openfda.generic_name:{drug}", "limit": 1})
        fda_data = request_json(fda_url) or {}
        write_json(openfda_dir / f"{drug.replace(' ', '_')}.json", fda_data)
        openfda_results[drug] = fda_data
        time.sleep(0.1)
    return dailymed_results, openfda_results


def civic_catalog_rows() -> list[dict]:
    path = PROCESSED / "civic_egfr_evidence.csv"
    if not path.exists():
        return []
    counts = Counter()
    sources = defaultdict(set)
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            bucket = mutation_bucket(row.get("requested_variant") or row.get("civic_variant") or "")
            if bucket:
                counts[bucket] += 1
                if row.get("source"):
                    sources[bucket].add(row["source"])
    rows = []
    for mutation, count in sorted(counts.items()):
        rows.append(
            reference_row(
                {
                    "mutation": mutation,
                    "gene": "EGFR",
                    "evidence_type": "variant_drug_evidence_count",
                    "record_count": count,
                    "clinical_significance": "",
                    "notes": "; ".join(sorted(sources[mutation]))[:500],
                },
                "CIViC",
                "CIViC EGFR evidence",
                "data/processed/civic_egfr_evidence.csv",
            )
        )
    return rows


def clinvar_catalog_rows(clinvar_data: dict[str, dict]) -> list[dict]:
    rows = []
    for mutation, data in clinvar_data.items():
        result = data.get("esearchresult", {}) if isinstance(data, dict) else {}
        rows.append(
            reference_row(
                {
                    "mutation": mutation,
                    "gene": "EGFR",
                    "evidence_type": "clinvar_query_count",
                    "record_count": int(result.get("count") or 0),
                    "clinical_significance": "",
                    "notes": ";".join(result.get("idlist", [])),
                },
                "ClinVar",
                "ClinVar EGFR query",
                "NCBI E-utilities ClinVar",
            )
        )
    return rows


def cbio_catalog_rows() -> list[dict]:
    rows = []
    path = PROCESSED / "cbioportal_egfr_sample_flags.csv"
    if not path.exists():
        return rows
    counts = Counter()
    compound = 0
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            variants = [item for item in row.get("egfr_variants", "").split(";") if item]
            if len(variants) > 1:
                compound += 1
            for variant in variants:
                bucket = mutation_bucket(variant)
                if bucket:
                    counts[bucket] += 1
    counts["compound_EGFR"] += compound
    for mutation in KEY_MUTATIONS:
        rows.append(
            reference_row(
                {
                    "mutation": mutation,
                    "gene": "EGFR",
                    "evidence_type": "cbioportal_sample_count",
                    "record_count": counts.get(mutation, 0),
                    "clinical_significance": "",
                    "notes": "Counts from downloaded cBioPortal EGFR sample flags.",
                },
                "cBioPortal",
                "cBioPortal downloaded cohorts",
                "data/processed/cbioportal_egfr_sample_flags.csv",
            )
        )
    return rows


def cancer_hotspots_catalog_rows(hotspots_path: Path | None) -> list[dict]:
    status = "downloaded" if hotspots_path and hotspots_path.exists() else "not_available"
    return [
        reference_row(
            {
                "mutation": mutation,
                "gene": "EGFR",
                "evidence_type": "hotspot_reference_file_status",
                "record_count": 1 if status == "downloaded" else 0,
                "clinical_significance": "",
                "notes": f"Cancer Hotspots file status: {status}; xls parsing intentionally deferred from MVP.",
            },
            "Cancer Hotspots",
            "Cancer Hotspots",
            "https://www.cancerhotspots.org/files/hotspots_v2.xls",
        )
        for mutation in KEY_MUTATIONS
    ]


def gdc_scan() -> tuple[list[dict], list[dict]]:
    maf_dir = RAW / "gdc" / "maf_files"
    if not maf_dir.exists():
        warn("GDC MAF directory missing; frequency table will use cBioPortal only.")
        return [], []
    file_to_project = {}
    manifest = RAW / "gdc" / "tcga_luad_lusc_masked_somatic_mutation_files.json"
    if manifest.exists():
        for hit in read_json(manifest).get("data", {}).get("hits", []):
            projects = {
                case.get("project", {}).get("project_id", "")
                for case in hit.get("cases", [])
                if case.get("project")
            }
            file_to_project[hit.get("file_name", "")] = ";".join(sorted(p for p in projects if p))
    sample_info: dict[str, dict] = {}
    mutation_counts = Counter()
    for path in sorted(maf_dir.glob("*.maf.gz")):
        if path.name.startswith("._"):
            continue
        project = file_to_project.get(path.name, "TCGA")
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            reader = None
            for line in handle:
                if line.startswith("#"):
                    continue
                if line.startswith("Hugo_Symbol"):
                    reader = csv.DictReader(chain([line], handle), delimiter="\t")
                    break
            if reader is None:
                continue
            for row in reader:
                sample = row.get("Tumor_Sample_Barcode", "")
                if not sample:
                    continue
                info = sample_info.setdefault(
                    sample,
                    {"cohort": project, "genes": set(), "egfr_variants": set(), "egfr_mutation_count": 0},
                )
                gene = row.get("Hugo_Symbol", "")
                if gene:
                    info["genes"].add(gene)
                if gene == "EGFR":
                    pc = normalize_protein_change(row.get("HGVSp_Short") or row.get("Protein_Change") or row.get("HGVSp") or "")
                    if pc:
                        info["egfr_variants"].add(pc)
                        info["egfr_mutation_count"] += 1
                        bucket = mutation_bucket(pc, row.get("Variant_Type", ""), row.get("Variant_Classification", ""))
                        if bucket:
                            mutation_counts[(project, bucket)] += 1
    catalog_rows = []
    for (cohort, mutation), count in sorted(mutation_counts.items()):
        catalog_rows.append(
            reference_row(
                {
                    "mutation": mutation,
                    "gene": "EGFR",
                    "evidence_type": "tcga_gdc_mutation_count",
                    "record_count": count,
                    "clinical_significance": "",
                    "notes": "Counted from downloaded open GDC MAF.gz files.",
                },
                "TCGA/GDC",
                cohort,
                "data/raw/gdc/maf_files",
            )
        )
    freq_rows = frequency_rows_from_samples(sample_info, "TCGA/GDC", "data/raw/gdc/maf_files")
    return catalog_rows, freq_rows


def frequency_rows_from_samples(sample_info: dict[str, dict], source: str, source_id: str) -> list[dict]:
    cohorts: dict[str, list[dict]] = defaultdict(list)
    for info in sample_info.values():
        cohorts[info["cohort"]].append(info)
    rows = []
    for cohort, samples in sorted(cohorts.items()):
        sample_count = len(samples)
        egfr_samples = [info for info in samples if info["egfr_variants"]]
        pacc_samples = [info for info in egfr_samples if any(is_pacc_variant(v) for v in info["egfr_variants"])]
        compound_samples = [info for info in egfr_samples if len(info["egfr_variants"]) > 1 or info["egfr_mutation_count"] > 1]
        tp53_samples = [info for info in egfr_samples if "TP53" in info["genes"]]
        rb1_samples = [info for info in egfr_samples if "RB1" in info["genes"]]
        metrics = [
            ("EGFR PACC frequency", len(pacc_samples), sample_count),
            ("compound EGFR frequency", len(compound_samples), max(len(egfr_samples), 1)),
            ("EGFR+TP53 frequency", len(tp53_samples), max(len(egfr_samples), 1)),
            ("EGFR+RB1 frequency", len(rb1_samples), max(len(egfr_samples), 1)),
        ]
        for metric, mutation_count, denominator in metrics:
            rows.append(
                reference_row(
                    {
                        "cancer_type": "LUAD/LUSC" if cohort.startswith("TCGA") else "",
                        "sample_count": denominator,
                        "mutation_count": mutation_count,
                        "frequency": round(mutation_count / denominator, 6) if denominator else "",
                        "metric": metric,
                        "ancestry_race_available": False,
                    },
                    source,
                    cohort,
                    source_id,
                )
            )
    return rows


def cbio_frequency_rows() -> list[dict]:
    path = PROCESSED / "cbioportal_egfr_sample_flags.csv"
    if not path.exists():
        return []
    sample_info = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            variants = {v for v in row.get("egfr_variants", "").split(";") if v}
            sample_info[f"{row['study_id']}::{row['sample_id']}"] = {
                "cohort": row["study_id"],
                "genes": {gene for gene, flag in [("TP53", row.get("tp53_mut")), ("RB1", row.get("rb1_mut"))] if str(flag).lower() == "true"},
                "egfr_variants": variants,
                "egfr_mutation_count": len(variants),
            }
    return frequency_rows_from_samples(sample_info, "cBioPortal", "data/processed/cbioportal_egfr_sample_flags.csv")


def reference_row(values: dict, source: str, cohort: str, source_url_or_id: str) -> dict:
    return {
        **values,
        "source": source,
        "cohort": cohort,
        "source_url_or_id": source_url_or_id,
        "generated_at": GENERATED_AT,
        "license_note": LICENSE_NOTES.get(source, ""),
    }


def load_fda_label_evidence() -> dict[str, dict]:
    path = REFERENCE / "fda_label_cns_evidence.csv"
    if not path.exists():
        warn(f"FDA label CNS evidence not found yet: {path}")
        return {}
    with path.open(encoding="utf-8") as handle:
        return {row["drug"].lower(): row for row in csv.DictReader(handle)}


def build_cns_table(pubchem: dict[str, dict], dailymed: dict[str, dict], openfda: dict[str, dict], fda_label_evidence: dict[str, dict]) -> list[dict]:
    chembl_path = PROCESSED / "chembl_ligands.csv"
    chembl_by_drug = {}
    if chembl_path.exists():
        with chembl_path.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                chembl_by_drug[row["drug"].lower()] = row

    manual_cns = {
        "gefitinib": ("low", -0.04),
        "erlotinib": ("low_to_moderate", -0.02),
        "afatinib": ("moderate", 0.0),
        "dacomitinib": ("low_to_moderate", -0.02),
        "osimertinib": ("high", 0.10),
        "lazertinib": ("high", 0.08),
        "furmonertinib": ("high", 0.10),
        "aumolertinib": ("high", 0.08),
        "mobocertinib": ("low", -0.06),
        "amivantamab": ("low", -0.08),
        "patritumab deruxtecan": ("low", -0.08),
    }
    rows = []
    for drug in DRUGS:
        pc = pubchem.get(drug, {}).get("PropertyTable", {}).get("Properties", [{}])[0]
        dm_records = dailymed.get(drug, {}).get("data", []) if isinstance(dailymed.get(drug), dict) else []
        fda_results = openfda.get(drug, {}).get("results", []) if isinstance(openfda.get(drug), dict) else []
        cns_activity, cns_modifier = manual_cns[drug]
        source_parts = ["manual_literature_seed", "PubChem"]
        source_ids = ["manual_literature_curation_pending", f"PubChem CID:{pc.get('CID', '')}"]
        if drug in chembl_by_drug:
            source_parts.append("ChEMBL")
            source_ids.append(chembl_by_drug[drug].get("chembl_id", ""))
        if dm_records:
            source_parts.append("DailyMed")
            source_ids.append(dm_records[0].get("setid", ""))
        if fda_results:
            source_parts.append("openFDA")
            source_ids.append("openFDA drug label")
        fda_evidence = fda_label_evidence.get(drug, {})
        fda_strength = fda_evidence.get("fda_label_evidence_strength", "")
        if fda_evidence:
            source_parts.append("FDA_label_CNS_evidence")
            source_ids.append(f"FDA label effective_time:{fda_evidence.get('effective_time', '')}")
        notes = "CNS class is a manual literature seed and must be reviewed before final prior freeze."
        if fda_strength:
            notes = f"FDA label CNS/PK evidence={fda_strength}. {fda_evidence.get('cns_label_interpretation', '')} Manual CNS class remains a prior pending clinical literature lock."
        rows.append(
            {
                "drug": drug,
                "cns_activity": cns_activity,
                "cns_modifier_score": cns_modifier,
                "pubchem_cid": pc.get("CID", ""),
                "molecular_weight": pc.get("MolecularWeight", ""),
                "xlogp": pc.get("XLogP", ""),
                "tpsa": pc.get("TPSA", ""),
                "chembl_id": chembl_by_drug.get(drug, {}).get("chembl_id", ""),
                "dailymed_setid": dm_records[0].get("setid", "") if dm_records else "",
                "openfda_label_found": bool(fda_results),
                "fda_label_evidence_strength": fda_strength,
                "fda_label_effective_time": fda_evidence.get("effective_time", ""),
                "notes": notes,
                "source": ";".join(source_parts),
                "cohort": "EGFR TKI / antibody",
                "source_url_or_id": ";".join(str(item) for item in source_ids if item),
                "generated_at": GENERATED_AT,
                "license_note": "; ".join(LICENSE_NOTES.get(item, "") for item in source_parts),
            }
        )
    return rows


def write_source_priority_map() -> None:
    rows = [
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 1, "cBioPortal", "public API/downloaded processed flags", "required_open", LICENSE_NOTES["cBioPortal"], "Somatic cohort support."),
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 2, "TCGA/GDC", "data/raw/gdc/maf_files", "required_open", LICENSE_NOTES["TCGA/GDC"], "Open TCGA MAF support."),
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 3, "CIViC", "data/processed/civic_egfr_evidence.csv", "required_open", LICENSE_NOTES["CIViC"], "Evidence/actionability support."),
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 4, "ClinVar", "NCBI E-utilities", "required_open", LICENSE_NOTES["ClinVar"], "Clinical variant archive query count."),
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 5, "Cancer Hotspots", "https://www.cancerhotspots.org/files/hotspots_v2.xls", "required_open", LICENSE_NOTES["Cancer Hotspots"], "Downloaded; xls parsing deferred."),
        ("egfr_open_mutation_catalog.csv", "mutation/catalog", 6, "OncoKB", "ONCOKB_TOKEN", "optional_enhanced", LICENSE_NOTES["optional_enhanced"], "Enhances actionability only."),
        ("egfr_enhanced_mutation_catalog.csv", "mutation/catalog", 1, "egfr_open_mutation_catalog.csv", "data/processed/reference/egfr_open_mutation_catalog.csv", "required_open_baseline", "Inherited from source rows.", "Baseline rows copied into enhanced catalog."),
        ("egfr_enhanced_mutation_catalog.csv", "mutation/catalog", 2, "COSMIC MutantCensus", "data/raw/cosmic/Cosmic_MutantCensus_v104_GRCh37.tsv.gz", "optional_enhanced", LICENSE_NOTES["optional_enhanced"], "Enhances mutation bucket counts and exact AA summaries."),
        ("egfr_enhanced_mutation_catalog.csv", "drug_resistance", 3, "COSMIC ResistanceMutations", "data/raw/cosmic/Cosmic_ResistanceMutations_v104_GRCh37.tsv.gz", "optional_enhanced", LICENSE_NOTES["optional_enhanced"], "Enhances C797S/L718Q resistance catalog and drug-response links."),
        ("open_egfr_pacc_frequency.csv", "frequency", 1, "TCGA/GDC", "data/raw/gdc/maf_files", "required_open", LICENSE_NOTES["TCGA/GDC"], "Primary open frequency fallback."),
        ("open_egfr_pacc_frequency.csv", "frequency", 2, "cBioPortal", "data/processed/cbioportal_egfr_sample_flags.csv", "required_open", LICENSE_NOTES["cBioPortal"], "Cross-check cohort frequency."),
        ("open_egfr_pacc_frequency.csv", "frequency", 3, "MSK-IMPACT published cohorts", "cBioPortal datahub if available", "optional_open", LICENSE_NOTES["cBioPortal"], "Not required for MVP."),
        ("open_egfr_pacc_frequency.csv", "frequency", 4, "TRACERx processed data", "study download if available", "optional_open", "Review original study terms.", "Not required for MVP."),
        ("open_egfr_pacc_frequency.csv", "frequency", 5, "NGDC/CNGBdb", "China/Asia cohorts if available", "optional_enhanced", "Review cohort terms.", "Adds Asia-specific priors."),
        ("egfr_tki_cns_modifier_table.csv", "drug_identity/chemistry", 1, "ChEMBL", "data/processed/chembl_ligands.csv", "required_open", LICENSE_NOTES["ChEMBL"], "Ligand IDs and structures."),
        ("egfr_tki_cns_modifier_table.csv", "drug_identity/chemistry", 2, "PubChem", "PUG REST", "required_open", LICENSE_NOTES["PubChem"], "CID and computed properties."),
        ("egfr_tki_cns_modifier_table.csv", "label", 3, "DailyMed", "DailyMed API", "required_open", LICENSE_NOTES["DailyMed"], "US label if available."),
        ("egfr_tki_cns_modifier_table.csv", "label", 4, "openFDA", "openFDA drug label", "required_open", LICENSE_NOTES["openFDA"], "FDA label if available."),
        ("egfr_tki_cns_modifier_table.csv", "cns_pk_label_evidence", 5, "FDA_label_CNS_evidence", "data/processed/reference/fda_label_cns_evidence.csv", "required_open", LICENSE_NOTES["openFDA"], "Preferred open replacement for DrugBank CNS/PK label evidence."),
        ("egfr_tki_cns_modifier_table.csv", "adverse_events", 6, "SIDER", "manual import", "optional_open", "Review SIDER terms.", "Enhances toxicity only."),
        ("egfr_tki_cns_modifier_table.csv", "pharmacology", 7, "IUPHAR", "manual/API import", "optional_open", "Review IUPHAR terms.", "Enhances target pharmacology."),
        ("egfr_tki_cns_modifier_table.csv", "pharmacogenomics", 8, "PharmGKB", "manual import", "optional_enhanced", "Review PharmGKB terms.", "Optional."),
        ("egfr_tki_cns_modifier_table.csv", "cns_modifier_score", 9, "manual_literature_seed", "manual_literature_curation_pending", "required_manual_before_freeze", LICENSE_NOTES["manual_literature_seed"], "Must be reviewed before clinical outcome analysis."),
        ("all", "enhanced_catalog/frequency/drug", 99, "COSMIC/GENIE/DrugBank", "data/raw/cosmic|genie|drugbank", "optional_enhanced", LICENSE_NOTES["optional_enhanced"], "Absence must warn only, not fail."),
    ]
    write_csv(
        REFERENCE / "source_priority_map.csv",
        [
            {
                "output_table": row[0],
                "field": row[1],
                "priority_rank": row[2],
                "source": row[3],
                "cohort": "source_priority_map",
                "source_url_or_id": row[4],
                "dependency_status": row[5],
                "license_note": row[6],
                "notes": row[7],
                "generated_at": GENERATED_AT,
            }
            for row in rows
        ],
        ["output_table", "field", "priority_rank", "source", "cohort", "source_url_or_id", "dependency_status", "license_note", "notes", "generated_at"],
    )


def main() -> None:
    ensure_optional_dirs()
    hotspots_path = download_cancer_hotspots()
    clinvar = fetch_clinvar()
    pubchem = fetch_pubchem()
    dailymed, openfda = fetch_dailymed_openfda()

    gdc_catalog, gdc_frequency = gdc_scan()
    catalog_rows = []
    catalog_rows.extend(cbio_catalog_rows())
    catalog_rows.extend(gdc_catalog)
    catalog_rows.extend(civic_catalog_rows())
    catalog_rows.extend(clinvar_catalog_rows(clinvar))
    catalog_rows.extend(cancer_hotspots_catalog_rows(hotspots_path))

    frequency_rows = []
    frequency_rows.extend(cbio_frequency_rows())
    frequency_rows.extend(gdc_frequency)

    fda_label_evidence = load_fda_label_evidence()
    cns_rows = build_cns_table(pubchem, dailymed, openfda, fda_label_evidence)

    catalog_fields = [
        "mutation",
        "gene",
        "evidence_type",
        "record_count",
        "clinical_significance",
        "notes",
        "source",
        "cohort",
        "source_url_or_id",
        "generated_at",
        "license_note",
    ]
    frequency_fields = [
        "cancer_type",
        "cohort",
        "metric",
        "sample_count",
        "mutation_count",
        "frequency",
        "ancestry_race_available",
        "source",
        "source_url_or_id",
        "generated_at",
        "license_note",
    ]
    cns_fields = [
        "drug",
        "cns_activity",
        "cns_modifier_score",
        "pubchem_cid",
        "molecular_weight",
        "xlogp",
        "tpsa",
        "chembl_id",
        "dailymed_setid",
        "openfda_label_found",
        "fda_label_evidence_strength",
        "fda_label_effective_time",
        "notes",
        "source",
        "cohort",
        "source_url_or_id",
        "generated_at",
        "license_note",
    ]
    write_csv(REFERENCE / "egfr_open_mutation_catalog.csv", catalog_rows, catalog_fields)
    write_csv(REFERENCE / "open_egfr_pacc_frequency.csv", frequency_rows, frequency_fields)
    write_csv(REFERENCE / "egfr_tki_cns_modifier_table.csv", cns_rows, cns_fields)
    write_source_priority_map()
    print(
        {
            "mutation_catalog_rows": len(catalog_rows),
            "frequency_rows": len(frequency_rows),
            "cns_modifier_rows": len(cns_rows),
            "reference_dir": str(REFERENCE),
        }
    )


if __name__ == "__main__":
    main()
