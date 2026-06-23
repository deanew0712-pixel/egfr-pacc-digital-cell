from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from egfr_pacc.classifier import classify_variants
from egfr_pacc.ranking import CaseFeatures, rank_drugs
from egfr_pacc.reporting import write_markdown_report


def main() -> None:
    case = CaseFeatures(
        patient_id="DEMO_PACC_COMPOUND_CNS",
        egfr_variants=["G719S", "L861Q"],
        brain_metastasis=True,
        treatment_line=1,
        tp53_mut=True,
    )
    classification = classify_variants(case.egfr_variants)
    ranking = rank_drugs(case, classification)

    reports = ROOT / "reports"
    reports.mkdir(exist_ok=True)
    ranking_out = reports / "demo_case_ranking.csv"
    report_out = reports / "demo_case_report.md"
    ranking.to_csv(ranking_out, index=False)
    write_markdown_report(case, classification, ranking, report_out)
    print({"ranking": str(ranking_out), "report": str(report_out)})


if __name__ == "__main__":
    main()

