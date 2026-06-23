from __future__ import annotations

import csv
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw" / "fda_labels"
REFERENCE = ROOT / "data" / "processed" / "reference"
REPORTS = ROOT / "reports"

USER_AGENT = "egfr-pacc-digital-cell/0.1"
GENERATED_AT = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

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

DIRECT_CNS_PATTERNS = [
    r"\bCNS\b",
    r"central nervous system",
    r"brain metast",
    r"cerebral metast",
    r"intracranial",
    r"leptomeningeal",
    r"blood-brain",
    r"cerebrospinal",
    r"\bCSF\b",
]

PK_PATTERNS = [
    r"volume of distribution",
    r"distribution",
    r"pharmacokinetic",
    r"exposure",
    r"steady-state",
    r"plasma",
]

TRANSPORTER_PATTERNS = [
    r"\bP-gp\b",
    r"P-glycoprotein",
    r"\bBCRP\b",
    r"\bABCB1\b",
    r"\bABCG2\b",
    r"transporter",
]

SECTIONS = [
    "indications_and_usage",
    "clinical_pharmacology",
    "pharmacokinetics",
    "clinical_studies",
    "warnings_and_cautions",
    "use_in_specific_populations",
    "description",
]


def request_json(url: str) -> dict:
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urlopen(req, timeout=45) as response:
        return json.loads(response.read().decode("utf-8"))


def safe_name(drug: str) -> str:
    return drug.replace(" ", "_").replace("/", "_")


def openfda_url(drug: str) -> str:
    return "https://api.fda.gov/drug/label.json?" + urlencode(
        {"search": f"openfda.generic_name:{drug}", "limit": 1}
    )


def download_labels() -> dict[str, dict]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    labels = {}
    for drug in DRUGS:
        out = RAW_DIR / f"{safe_name(drug)}_label.json"
        if out.exists() and out.stat().st_size > 0:
            labels[drug] = json.loads(out.read_text(encoding="utf-8"))
            continue
        try:
            data = request_json(openfda_url(drug))
        except Exception as exc:
            data = {"error": str(exc), "results": []}
        out.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
        labels[drug] = data
        time.sleep(0.1)
    return labels


def label_text(result: dict) -> dict[str, str]:
    texts = {}
    for section in SECTIONS:
        value = result.get(section, [])
        if isinstance(value, list):
            texts[section] = " ".join(str(item) for item in value)
        else:
            texts[section] = str(value)
    return texts


def count_patterns(text: str, patterns: list[str]) -> int:
    return sum(len(re.findall(pattern, text, flags=re.IGNORECASE)) for pattern in patterns)


def snippets(texts: dict[str, str], patterns: list[str], limit: int = 4) -> list[str]:
    hits = []
    joined = " ".join(f"[{section}] {text}" for section, text in texts.items())
    for pattern in patterns:
        for match in re.finditer(pattern, joined, flags=re.IGNORECASE):
            start = max(0, match.start() - 180)
            end = min(len(joined), match.end() + 220)
            snippet = re.sub(r"\s+", " ", joined[start:end]).strip()
            if snippet not in hits:
                hits.append(snippet)
            if len(hits) >= limit:
                return hits
    return hits


def classify_evidence(texts: dict[str, str]) -> tuple[str, str]:
    clinical_text = " ".join(texts.get(section, "") for section in ["indications_and_usage", "clinical_studies"])
    all_text = " ".join(texts.values())
    direct_hits = count_patterns(clinical_text, DIRECT_CNS_PATTERNS)
    all_direct_hits = count_patterns(all_text, DIRECT_CNS_PATTERNS)
    transporter_hits = count_patterns(all_text, TRANSPORTER_PATTERNS)
    pk_hits = count_patterns(all_text, PK_PATTERNS)

    if direct_hits:
        return (
            "direct_clinical_cns",
            "FDA label contains CNS/brain-metastasis language in indication or clinical studies sections.",
        )
    if all_direct_hits:
        return (
            "direct_label_cns_nonclinical_or_context",
            "FDA label contains CNS-related language outside primary indication/clinical studies sections.",
        )
    if transporter_hits and pk_hits:
        return (
            "indirect_pk_transporter",
            "FDA label contains transporter and pharmacokinetic/distribution evidence but no direct CNS efficacy claim.",
        )
    if transporter_hits:
        return (
            "indirect_transporter_only",
            "FDA label contains transporter evidence but no direct CNS efficacy claim.",
        )
    if pk_hits:
        return (
            "indirect_pk_only",
            "FDA label contains pharmacokinetic/distribution evidence but no direct CNS efficacy claim.",
        )
    return ("not_found", "No CNS, transporter, or distribution evidence detected in parsed label sections.")


def build_rows(labels: dict[str, dict]) -> list[dict]:
    rows = []
    for drug, data in labels.items():
        results = data.get("results", []) if isinstance(data, dict) else []
        result = results[0] if results else {}
        texts = label_text(result)
        all_text = " ".join(texts.values())
        evidence_strength, interpretation = classify_evidence(texts) if result else ("label_missing", "No openFDA label result.")
        direct_count = count_patterns(all_text, DIRECT_CNS_PATTERNS)
        pk_count = count_patterns(all_text, PK_PATTERNS)
        transporter_count = count_patterns(all_text, TRANSPORTER_PATTERNS)
        evidence_snippets = snippets(texts, DIRECT_CNS_PATTERNS + TRANSPORTER_PATTERNS + PK_PATTERNS)
        openfda = result.get("openfda", {}) if result else {}
        rows.append(
            {
                "drug": drug,
                "label_found": bool(result),
                "brand_names": ";".join(openfda.get("brand_name", []) or []),
                "generic_names": ";".join(openfda.get("generic_name", []) or []),
                "set_id": result.get("set_id", ""),
                "effective_time": result.get("effective_time", ""),
                "direct_cns_term_count": direct_count,
                "pk_term_count": pk_count,
                "transporter_term_count": transporter_count,
                "fda_label_evidence_strength": evidence_strength,
                "cns_label_interpretation": interpretation,
                "evidence_snippets": " || ".join(evidence_snippets),
                "source": "openFDA",
                "cohort": "FDA drug label",
                "source_url_or_id": openfda_url(drug),
                "generated_at": GENERATED_AT,
                "license_note": "openFDA public API; data are unvalidated and subject to openFDA terms.",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "drug",
        "label_found",
        "brand_names",
        "generic_names",
        "set_id",
        "effective_time",
        "direct_cns_term_count",
        "pk_term_count",
        "transporter_term_count",
        "fda_label_evidence_strength",
        "cns_label_interpretation",
        "evidence_snippets",
        "source",
        "cohort",
        "source_url_or_id",
        "generated_at",
        "license_note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    counts = {}
    for row in rows:
        key = row["fda_label_evidence_strength"]
        counts[key] = counts.get(key, 0) + 1
    lines = [
        "# FDA Label CNS Evidence Report",
        "",
        f"Generated at: {GENERATED_AT}",
        "",
        "## Summary",
        "",
    ]
    for key, value in sorted(counts.items()):
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Drug-Level Evidence", ""])
    for row in rows:
        lines.append(
            f"- {row['drug']}: {row['fda_label_evidence_strength']} | "
            f"direct_cns_terms={row['direct_cns_term_count']} | "
            f"pk_terms={row['pk_term_count']} | transporter_terms={row['transporter_term_count']} | "
            f"effective_time={row['effective_time'] or 'NA'}"
        )
    lines.extend(
        [
            "",
            "## Interpretation Rule",
            "",
            "- `direct_clinical_cns`: CNS/brain-metastasis language appears in indications or clinical studies.",
            "- `indirect_pk_transporter`: transporter plus PK/distribution evidence, without direct CNS efficacy language.",
            "- `indirect_transporter_only` or `indirect_pk_only`: weaker label-derived support.",
            "- These labels replace DrugBank as the preferred open source for label-level CNS/PK evidence.",
        ]
    )
    (REPORTS / "fda_label_cns_evidence_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    labels = download_labels()
    rows = build_rows(labels)
    write_csv(REFERENCE / "fda_label_cns_evidence.csv", rows)
    write_report(rows)
    print(
        {
            "raw_dir": str(RAW_DIR),
            "evidence_csv": str(REFERENCE / "fda_label_cns_evidence.csv"),
            "report": str(REPORTS / "fda_label_cns_evidence_report.md"),
            "rows": len(rows),
        }
    )


if __name__ == "__main__":
    main()
