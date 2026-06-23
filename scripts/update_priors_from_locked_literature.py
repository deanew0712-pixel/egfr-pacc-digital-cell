from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
LOCKED_TABLE = ROOT / "data" / "processed" / "literature_orr_locked.csv"
CNS_TABLE = ROOT / "data" / "processed" / "reference" / "egfr_tki_cns_modifier_table.csv"
PRIORS_OUT = ROOT / "configs" / "priors.yaml"
DRUG_PROPS_OUT = ROOT / "configs" / "drug_properties.yaml"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def evidence_sources(rows: list[dict], mutation_filter: set[str] | None = None, drug_filter: set[str] | None = None) -> str:
    selected = []
    for row in rows:
        if mutation_filter and row["mutation_bucket"] not in mutation_filter:
            continue
        if drug_filter and row["drug"] not in drug_filter:
            continue
        selected.append(f"{row['source_study']}:{row['lock_id']}")
    return "; ".join(selected[:8])


def cns_map(rows: list[dict]) -> dict[str, dict]:
    return {row["drug"]: row for row in rows}


def drug_generation(drug: str) -> str:
    if drug in {"gefitinib", "erlotinib"}:
        return "first"
    if drug in {"afatinib", "dacomitinib"}:
        return "second"
    if drug in {"osimertinib", "lazertinib", "furmonertinib", "aumolertinib"}:
        return "third"
    if drug == "mobocertinib":
        return "ex20ins_tki"
    if drug == "amivantamab":
        return "antibody"
    if drug == "patritumab deruxtecan":
        return "antibody_drug_conjugate"
    return "other"


def build_drug_properties(cns_rows: list[dict]) -> dict:
    drugs = {}
    for drug, row in cns_map(cns_rows).items():
        drugs[drug] = {
            "generation": drug_generation(drug),
            "covalent": drug in {"afatinib", "dacomitinib", "osimertinib", "lazertinib", "furmonertinib", "aumolertinib", "mobocertinib"},
            "cns_activity": row["cns_activity"],
            "cns_modifier_score": float(row["cns_modifier_score"]),
            "pubchem_cid": row.get("pubchem_cid") or None,
            "chembl_id": row.get("chembl_id") or None,
            "openfda_label_found": str(row.get("openfda_label_found")).lower() == "true",
            "fda_label_evidence_strength": row.get("fda_label_evidence_strength") or None,
            "fda_label_effective_time": row.get("fda_label_effective_time") or None,
            "availability_note": "open fallback CNS modifier; review against local formulary and labels",
            "source": row["source"],
            "source_url_or_id": row["source_url_or_id"],
        }
    return {"drugs": drugs}


def prior(score: float, level: str, source: str) -> dict:
    return {"score": round(score, 3), "evidence_level": level, "source": source}


def build_priors(locked_rows: list[dict]) -> dict:
    pacc_sources = evidence_sources(
        locked_rows,
        {"G719X", "S768I", "E709X", "major_uncommon", "compound_uncommon", "PACC", "sensitizing_uncommon_non_ex20ins_non_T790M"},
    )
    osi_sources = evidence_sources(locked_rows, {"G719X", "S768I", "E709X", "L861Q", "compound_uncommon", "uncommon_non_ex20ins"}, {"osimertinib"})
    afatinib_sources = evidence_sources(locked_rows, drug_filter={"afatinib", "second_generation_TKI"})
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return {
        "metadata": {
            "status": "open_fallback_literature_locked_candidate",
            "generated_at": generated_at,
            "warning": "Research-use prior candidate. Do not review clinical holdout outcomes before final freeze.",
            "locked_literature_table": "data/processed/literature_orr_locked.csv",
            "open_reference_inputs": [
                "data/processed/reference/egfr_open_mutation_catalog.csv",
                "data/processed/reference/open_egfr_pacc_frequency.csv",
                "data/processed/reference/egfr_tki_cns_modifier_table.csv",
                "data/processed/reference/source_priority_map.csv",
            ],
            "optional_enhanced_dependencies": ["OncoKB", "COSMIC", "GENIE", "DrugBank"],
        },
        "class_drug_priors": {
            "pacc": {
                "afatinib": prior(0.80, "B", f"Locked literature supports afatinib for major/PACC uncommon EGFR: {afatinib_sources}"),
                "dacomitinib": prior(0.70, "C", f"Second-generation TKI structure-class extrapolation from Robichaux PACC TTF and open chemistry; direct ORR lock absent. {pacc_sources}"),
                "osimertinib": prior(0.58, "B", f"UNICORN supports mixed uncommon EGFR activity but G719X subgroup PFS is shorter; use with mutation-specific caution. {osi_sources}"),
                "furmonertinib": prior(0.56, "D", "CNS-active third-generation TKI open fallback; PACC-specific locked ORR absent."),
                "lazertinib": prior(0.55, "D", "CNS-active third-generation TKI open fallback; PACC-specific locked ORR absent."),
                "aumolertinib": prior(0.54, "D", "CNS-active third-generation TKI open fallback; PACC-specific locked ORR absent."),
                "erlotinib": prior(0.42, "C", "Earlier-generation uncommon EGFR evidence exists but lower PACC prior than afatinib."),
                "gefitinib": prior(0.40, "C", "Earlier-generation uncommon EGFR evidence exists but lower PACC prior than afatinib."),
                "mobocertinib": prior(0.30, "D", "Exon20ins-oriented TKI; not a PACC-preferred prior."),
            },
            "classical_like": {
                "osimertinib": prior(0.82, "A", "Standard EGFR sensitizing mutation prior; retained as Baseline-1 comparator."),
                "lazertinib": prior(0.80, "B", "Third-generation CNS-active classical-like prior from open label/chemistry fallback."),
                "furmonertinib": prior(0.78, "B", "Third-generation CNS-active classical-like prior from open label/chemistry fallback."),
                "aumolertinib": prior(0.76, "B", "Third-generation CNS-active classical-like prior from open label/chemistry fallback."),
                "afatinib": prior(0.70, "B", "Active EGFR TKI prior; also supported in uncommon EGFR literature."),
                "dacomitinib": prior(0.68, "B", "Active EGFR TKI prior."),
                "erlotinib": prior(0.64, "B", "First-generation sensitizing mutation prior."),
                "gefitinib": prior(0.64, "B", "First-generation sensitizing mutation prior."),
            },
            "t790m_like": {
                "osimertinib": prior(0.86, "A", "T790M third-generation TKI prior."),
                "furmonertinib": prior(0.80, "B", "Third-generation TKI prior."),
                "lazertinib": prior(0.78, "B", "Third-generation TKI prior."),
                "aumolertinib": prior(0.76, "B", "Third-generation TKI prior."),
                "afatinib": prior(0.35, "C", "T790M reduces earlier-generation activity."),
                "dacomitinib": prior(0.35, "C", "T790M reduces earlier-generation activity."),
            },
            "ex20ins": {
                "amivantamab": prior(0.78, "B", "Exon20ins specialist antibody prior with DailyMed/openFDA support."),
                "mobocertinib": prior(0.62, "C", "Exon20ins-oriented TKI prior; availability and label status require review."),
                "patritumab deruxtecan": prior(0.52, "D", "HER3-DXd exploratory/specialist context; not mutation-class-specific."),
                "osimertinib": prior(0.40, "D", "Exon20ins response heterogeneous."),
                "afatinib": prior(0.35, "D", "Exon20ins response heterogeneous."),
            },
            "unknown": {
                "afatinib": prior(0.50, "D", "Unknown variant fallback; literature review required."),
                "osimertinib": prior(0.50, "D", "Unknown variant fallback; literature review required."),
            },
        },
        "modifiers": {
            "cns_high_activity_bonus": 0.10,
            "cns_low_activity_penalty": -0.08,
            "bypass_risk_penalty": -0.08,
            "compound_conflict_penalty": -0.04,
            "later_line_penalty": -0.03,
        },
    }


def main() -> None:
    locked_rows = read_csv(LOCKED_TABLE)
    cns_rows = read_csv(CNS_TABLE)
    priors = build_priors(locked_rows)
    drug_properties = build_drug_properties(cns_rows)
    PRIORS_OUT.write_text(yaml.safe_dump(priors, sort_keys=False, allow_unicode=True), encoding="utf-8")
    DRUG_PROPS_OUT.write_text(yaml.safe_dump(drug_properties, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(
        {
            "priors": str(PRIORS_OUT),
            "drug_properties": str(DRUG_PROPS_OUT),
            "locked_rows": len(locked_rows),
            "drug_rows": len(cns_rows),
        }
    )


if __name__ == "__main__":
    main()
