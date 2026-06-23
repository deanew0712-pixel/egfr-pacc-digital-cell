from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "processed" / "reference"
CNS_TABLE = REFERENCE / "egfr_tki_cns_modifier_table.csv"
FDA_EVIDENCE = REFERENCE / "fda_label_cns_evidence.csv"
SOURCE_PRIORITY = REFERENCE / "source_priority_map.csv"
GENERATED_AT = datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_once(parts: list[str], value: str) -> list[str]:
    if value and value not in parts:
        parts.append(value)
    return parts


def integrate_cns_table() -> int:
    cns_rows = read_csv(CNS_TABLE)
    evidence_by_drug = {row["drug"].lower(): row for row in read_csv(FDA_EVIDENCE)}
    fieldnames = list(cns_rows[0].keys())
    for field in ["fda_label_evidence_strength", "fda_label_effective_time"]:
        if field not in fieldnames:
            insert_at = fieldnames.index("notes") if "notes" in fieldnames else len(fieldnames)
            fieldnames.insert(insert_at, field)

    updated = 0
    for row in cns_rows:
        evidence = evidence_by_drug.get(row["drug"].lower())
        if not evidence:
            row.setdefault("fda_label_evidence_strength", "")
            row.setdefault("fda_label_effective_time", "")
            continue
        row["fda_label_evidence_strength"] = evidence.get("fda_label_evidence_strength", "")
        row["fda_label_effective_time"] = evidence.get("effective_time", "")
        source_parts = [part for part in row.get("source", "").split(";") if part]
        source_ids = [part for part in row.get("source_url_or_id", "").split(";") if part]
        append_once(source_parts, "FDA_label_CNS_evidence")
        append_once(source_ids, f"FDA label effective_time:{evidence.get('effective_time', '')}")
        row["source"] = ";".join(source_parts)
        row["source_url_or_id"] = ";".join(source_ids)
        row["notes"] = (
            f"FDA label CNS/PK evidence={evidence.get('fda_label_evidence_strength', '')}. "
            f"{evidence.get('cns_label_interpretation', '')} "
            "Manual CNS class remains a prior pending clinical literature lock."
        ).strip()
        if "openFDA public API" not in row.get("license_note", ""):
            row["license_note"] = f"{row.get('license_note', '')}; openFDA public API; data are unvalidated and subject to openFDA terms.".strip("; ")
        updated += 1
    write_csv(CNS_TABLE, cns_rows, fieldnames)
    return updated


def integrate_source_priority() -> None:
    rows = read_csv(SOURCE_PRIORITY)
    exists = any(
        row.get("output_table") == "egfr_tki_cns_modifier_table.csv"
        and row.get("field") == "cns_pk_label_evidence"
        and row.get("source") == "FDA_label_CNS_evidence"
        for row in rows
    )
    if not exists:
        rows.append(
            {
                "output_table": "egfr_tki_cns_modifier_table.csv",
                "field": "cns_pk_label_evidence",
                "priority_rank": "5",
                "source": "FDA_label_CNS_evidence",
                "cohort": "source_priority_map",
                "source_url_or_id": "data/processed/reference/fda_label_cns_evidence.csv",
                "dependency_status": "required_open",
                "license_note": "openFDA public API; data are unvalidated and subject to openFDA terms.",
                "notes": "Preferred open replacement for DrugBank CNS/PK label evidence.",
                "generated_at": GENERATED_AT,
            }
        )
    fieldnames = list(rows[0].keys())
    write_csv(SOURCE_PRIORITY, rows, fieldnames)


def main() -> None:
    updated = integrate_cns_table()
    integrate_source_priority()
    print({"cns_rows_updated": updated, "cns_table": str(CNS_TABLE), "source_priority": str(SOURCE_PRIORITY)})


if __name__ == "__main__":
    main()
