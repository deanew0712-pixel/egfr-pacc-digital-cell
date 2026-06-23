from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from egfr_pacc.classifier import classify_variants
from egfr_pacc.ranking import CaseFeatures, rank_drugs
from egfr_pacc.reporting import write_markdown_report


OUTCOME_LIKE_COLUMNS = {
    "best_response",
    "orr",
    "response",
    "pfs",
    "pfs_months",
    "pfs_event",
    "os",
    "os_months",
    "os_event",
    "death",
    "progression",
}


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "positive", "mut", "mutated"}


def assert_no_outcome_leakage(columns: list[str]) -> None:
    lower = {col.lower().strip() for col in columns}
    leaked = sorted(lower & OUTCOME_LIKE_COLUMNS)
    if leaked:
        raise ValueError(
            "Clinical feature file contains outcome-like columns. "
            f"Move these to a physically separate outcome file before ranking: {', '.join(leaked)}"
        )


def case_from_row(row: pd.Series) -> CaseFeatures:
    variants = [item.strip() for item in str(row["egfr_variants"]).split(";") if item.strip()]
    return CaseFeatures(
        patient_id=str(row["patient_id"]),
        egfr_variants=variants,
        brain_metastasis=parse_bool(row.get("brain_metastasis", False)),
        leptomeningeal_metastasis=parse_bool(row.get("leptomeningeal_metastasis", False)),
        treatment_line=int(row.get("treatment_line", 1)),
        tp53_mut=parse_bool(row.get("tp53_mut", False)),
        rb1_loss=parse_bool(row.get("rb1_loss", False)),
        met_amplification=parse_bool(row.get("met_amplification", False)),
        pik3ca_mut=parse_bool(row.get("pik3ca_mut", False)),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run EGFR-PACC MVP ranking for a feature-only clinical table.")
    parser.add_argument(
        "--features",
        default=str(ROOT / "data" / "clinical_holdout" / "clinical_features_template.csv"),
        help="Feature-only CSV. Must not contain outcome columns.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(ROOT / "reports" / "batch_demo"),
        help="Output directory for ranking tables and reports.",
    )
    args = parser.parse_args()

    features = pd.read_csv(args.features)
    assert_no_outcome_leakage(list(features.columns))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ranking_rows = []
    class_rows = []
    for _, row in features.iterrows():
        case = case_from_row(row)
        classification = classify_variants(case.egfr_variants)
        ranking = rank_drugs(case, classification)
        ranking_rows.append(ranking)
        class_rows.append(
            {
                "patient_id": case.patient_id,
                "normalized_variants": ";".join(classification.normalized_variants),
                "classes": ";".join(classification.classes),
                "primary_class": classification.primary_class,
                "compound_conflict_flag": classification.compound_conflict_flag,
                "conflict_resolution_rule": classification.conflict_resolution_rule,
                "uncertainty_delta": classification.uncertainty_delta,
            }
        )
        write_markdown_report(
            case,
            classification,
            ranking,
            out_dir / f"{case.patient_id}_report.md",
        )

    all_rankings = pd.concat(ranking_rows, ignore_index=True) if ranking_rows else pd.DataFrame()
    all_rankings.to_csv(out_dir / "drug_rankings.csv", index=False)
    pd.DataFrame(class_rows).to_csv(out_dir / "mutation_classifications.csv", index=False)
    print(
        {
            "rankings": str(out_dir / "drug_rankings.csv"),
            "classifications": str(out_dir / "mutation_classifications.csv"),
            "n_cases": int(len(features)),
        }
    )


if __name__ == "__main__":
    main()

