from __future__ import annotations

from dataclasses import dataclass, field


EVIDENCE_ORDER = {"A": 4, "B": 3, "C": 2, "D": 1}
GRADE_ORDER = {"A": 4, "B": 3, "C": 2, "D": 1}
ORDER_TO_GRADE = {value: key for key, value in GRADE_ORDER.items()}


@dataclass(frozen=True)
class EpDriInput:
    base_prior_score: float
    evidence_level: str = "D"
    cns_activity: str = "unknown"
    cns_present: bool = False
    bypass_risk: bool = False
    compound_conflict: bool = False
    treatment_line: int = 1
    classification_uncertainty_delta: int = 0
    modifiers: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class EpDriResult:
    epdri_score: float
    recommendation_grade: str
    uncertainty_tier: str
    evidence_level: str
    evidence_cap_grade: str
    uncertainty_discount: float
    modifier_notes: list[str]
    layer_scores: dict[str, float]


def clamp_score(value: float) -> float:
    return max(0.01, min(0.99, float(value)))


def uncertainty_tier(score: int) -> str:
    if score <= 0:
        return "low"
    if score == 1:
        return "moderate"
    return "high"


def grade_from_score(score: float) -> str:
    if score >= 0.78:
        return "A"
    if score >= 0.65:
        return "B"
    if score >= 0.45:
        return "C"
    return "D"


def cap_grade(grade: str, cap: str) -> str:
    return ORDER_TO_GRADE[min(GRADE_ORDER.get(grade, 1), GRADE_ORDER.get(cap, 1))]


def evidence_level_cap(evidence_level: str, uncertainty: str) -> str:
    level = evidence_level.upper()
    if level not in EVIDENCE_ORDER:
        level = "D"

    # Evidence levels cap the final recommendation. High uncertainty imposes a
    # further C cap even when score is high.
    cap = level
    if uncertainty == "moderate":
        cap = cap_grade(cap, "B")
    elif uncertainty == "high":
        cap = cap_grade(cap, "C")
    return cap


def uncertainty_discount(uncertainty: str) -> float:
    if uncertainty == "low":
        return 0.0
    if uncertainty == "moderate":
        return 0.03
    return 0.07


def compute_epdri(value: EpDriInput) -> EpDriResult:
    """Three-layer collapsed EPDRI.

    Layer 1: class-level hierarchical prior.
    Layer 2: clinical modifiers such as CNS, bypass risk, compound conflict.
    Layer 3: uncertainty discount plus evidence-level grade cap.
    """
    modifiers = value.modifiers or {}
    layer1 = clamp_score(value.base_prior_score)
    score = layer1
    modifier_notes: list[str] = []
    uncertainty_score = int(value.classification_uncertainty_delta)

    if value.cns_present:
        if value.cns_activity == "high":
            score += float(modifiers.get("cns_high_activity_bonus", 0.0))
            modifier_notes.append("CNS activity bonus")
        elif value.cns_activity in {"low", "low_to_moderate"}:
            score += float(modifiers.get("cns_low_activity_penalty", 0.0))
            modifier_notes.append("CNS activity penalty")

    if value.bypass_risk:
        score += float(modifiers.get("bypass_risk_penalty", 0.0))
        uncertainty_score += 1
        modifier_notes.append("co-mutation/bypass risk penalty")

    if value.compound_conflict:
        score += float(modifiers.get("compound_conflict_penalty", 0.0))
        modifier_notes.append("compound mutation conflict penalty")

    if value.treatment_line > 1:
        score += float(modifiers.get("later_line_penalty", 0.0))
        uncertainty_score += 1
        modifier_notes.append("later-line evidence penalty")

    layer2 = clamp_score(score)
    tier = uncertainty_tier(uncertainty_score)
    discount = uncertainty_discount(tier)
    layer3 = clamp_score(layer2 - discount)

    raw_grade = grade_from_score(layer3)
    cap = evidence_level_cap(value.evidence_level, tier)
    final_grade = cap_grade(raw_grade, cap)

    return EpDriResult(
        epdri_score=round(layer3, 3),
        recommendation_grade=final_grade,
        uncertainty_tier=tier,
        evidence_level=value.evidence_level.upper(),
        evidence_cap_grade=cap,
        uncertainty_discount=discount,
        modifier_notes=modifier_notes,
        layer_scores={
            "layer1_hierarchical_prior": round(layer1, 3),
            "layer2_clinical_modified": round(layer2, 3),
            "layer3_uncertainty_discounted": round(layer3, 3),
        },
    )
