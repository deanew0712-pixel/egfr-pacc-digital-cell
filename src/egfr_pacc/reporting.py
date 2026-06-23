from __future__ import annotations

from pathlib import Path

import pandas as pd

from .classifier import ClassificationResult
from .ranking import CaseFeatures


def write_markdown_report(
    case: CaseFeatures,
    classification: ClassificationResult,
    ranking: pd.DataFrame,
    out: str | Path,
) -> None:
    lines = [
        f"# EGFR-PACC Digital Cell Report: {case.patient_id}",
        "",
        "Research use only. Not for clinical decision-making without prospective validation.",
        "",
        "## Mutation Classification",
        "",
        f"- EGFR variants: {', '.join(classification.normalized_variants)}",
        f"- Classes detected: {', '.join(classification.classes)}",
        f"- Primary class for MVP ranking: `{classification.primary_class}`",
        f"- Compound conflict: `{classification.compound_conflict_flag}`",
        f"- Conflict rule: `{classification.conflict_resolution_rule}`",
        "",
        "## Clinical Modifiers",
        "",
        f"- Brain metastasis: `{case.brain_metastasis}`",
        f"- Leptomeningeal metastasis: `{case.leptomeningeal_metastasis}`",
        f"- Treatment line: `{case.treatment_line}`",
        f"- TP53 mutated: `{case.tp53_mut}`",
        f"- RB1 loss: `{case.rb1_loss}`",
        f"- MET amplification: `{case.met_amplification}`",
        f"- PIK3CA mutated: `{case.pik3ca_mut}`",
        "",
        "## Drug Ranking",
        "",
        ranking.to_markdown(index=False),
        "",
        "## Interpretation Boundary",
        "",
        "- Ranking is prior-dominated and intended for MDT discussion.",
        "- Structural/docking/MD calculations are not part of this MVP ranking.",
        "- A held-out clinical cohort should be used only after baselines and priors are frozen.",
    ]
    Path(out).write_text("\n".join(lines) + "\n", encoding="utf-8")

