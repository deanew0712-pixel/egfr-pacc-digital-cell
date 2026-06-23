from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .io import config_path, load_yaml


@dataclass(frozen=True)
class ClassificationResult:
    normalized_variants: list[str]
    classes: list[str]
    primary_class: str
    compound_conflict_flag: bool
    conflict_resolution_rule: str
    uncertainty_delta: int


def normalize_variant(raw: str) -> str:
    value = raw.strip()
    value = re.sub(r"^EGFR[:\s_-]*", "", value, flags=re.IGNORECASE)
    value = value.replace("p.", "").replace(" ", "")
    low = value.lower()
    if "exon19" in low or "ex19" in low or "del19" in low:
        return "exon19del"
    if "exon20" in low or "ex20" in low or "ins" in low and "20" in low:
        return "ex20ins"
    return value.upper()


def _variant_members(rules: dict) -> dict[str, set[str]]:
    members_by_class: dict[str, set[str]] = {}
    for class_name, class_rule in rules.get("classes", {}).items():
        members_by_class[class_name] = {
            normalize_variant(item) for item in class_rule.get("variants", [])
        }
    return members_by_class


def _class_for_variant(variant: str, members_by_class: dict[str, set[str]]) -> str:
    normalized = normalize_variant(variant)
    for class_name, members in members_by_class.items():
        if normalized in members:
            return class_name

    # Family-level aliases keep rules compact while allowing patient reports to
    # contain concrete substitutions such as G719R or E709T.
    if re.match(r"^G719[A-Z]$", normalized):
        return "pacc"
    if re.match(r"^E709[A-Z]$", normalized):
        return "pacc"
    if normalized.startswith("E709") and ">" in normalized:
        return "pacc"
    if normalized in {"DEL19", "EXON19DEL"}:
        return "classical_like"
    return "unknown"


def resolve_compound(classes: list[str], rules: dict | None = None) -> tuple[str, bool, str, int]:
    """Resolve compound EGFR classes into a single MVP primary class.

    The returned primary class is only for collapsed ranking. The full class list
    is preserved in `ClassificationResult.classes` so compound biology is not
    hidden from reports.
    """
    if rules is None:
        rules = load_yaml(config_path("mutation_class_rules.yaml"))

    unique_classes: list[str] = []
    for class_name in classes:
        if class_name not in unique_classes:
            unique_classes.append(class_name)

    if not unique_classes:
        return "unknown", False, "no_variant_detected", 2

    class_set = set(unique_classes)
    compound_rules = rules.get("compound_rules", {})
    conflict = len(class_set) > 1

    if not conflict and "unknown" not in class_set:
        return unique_classes[0], False, "single_class", 0

    if "ex20ins" in class_set and len(class_set) > 1:
        rule = compound_rules.get("any_plus_ex20ins", {})
        return (
            rule.get("primary_class", "ex20ins"),
            True,
            rule.get("recommendation", "ex20ins_specialist_review"),
            int(rule.get("uncertainty_delta", 2)),
        )

    if {"pacc", "t790m_like"} <= class_set:
        rule = compound_rules.get("pacc_plus_t790m_like", {})
        return (
            rule.get("primary_class", "t790m_like"),
            True,
            rule.get("recommendation", "acquired_resistance_context_required"),
            int(rule.get("uncertainty_delta", 2)),
        )

    if {"pacc", "classical_like"} <= class_set:
        rule = compound_rules.get("pacc_plus_classical_like", {})
        return (
            rule.get("primary_class", "pacc"),
            True,
            rule.get("recommendation", "preserve_both_classes_with_pacc_bias"),
            int(rule.get("uncertainty_delta", 1)),
        )

    if "unknown" in class_set:
        rule = rules.get("unknown_variant", {})
        if len(class_set) == 1:
            return (
                "unknown",
                False,
                rule.get("recommendation", "literature_review_required"),
                int(rule.get("uncertainty_delta", 2)),
            )
        priority = compound_rules.get(
            "priority_order", ["ex20ins", "t790m_like", "pacc", "classical_like", "unknown"]
        )
        primary = next((item for item in priority if item in class_set and item != "unknown"), "unknown")
        return (
            primary,
            True,
            rule.get("recommendation", "literature_review_required"),
            int(rule.get("uncertainty_delta", 2)),
        )

    priority = compound_rules.get(
        "priority_order", ["ex20ins", "t790m_like", "pacc", "classical_like", "unknown"]
    )
    primary = next((item for item in priority if item in class_set), unique_classes[0])
    return primary, conflict, "priority_order_fallback", 1 if conflict else 0


def classify_robichaux(
    variants: list[str],
    rules_path: str | Path | None = None,
) -> ClassificationResult:
    rules = load_yaml(rules_path or config_path("mutation_class_rules.yaml"))
    normalized = [normalize_variant(v) for v in variants if str(v).strip()]
    members_by_class = _variant_members(rules)
    class_hits: list[str] = []
    for variant in normalized:
        class_hits.append(_class_for_variant(variant, members_by_class))

    unique_classes = []
    for class_name in class_hits:
        if class_name not in unique_classes:
            unique_classes.append(class_name)

    if not unique_classes:
        unique_classes = ["unknown"]

    primary, conflict, resolution_rule, uncertainty_delta = resolve_compound(unique_classes, rules)

    return ClassificationResult(
        normalized_variants=normalized,
        classes=unique_classes,
        primary_class=primary,
        compound_conflict_flag=conflict,
        conflict_resolution_rule=resolution_rule,
        uncertainty_delta=uncertainty_delta,
    )


def classify_variants(
    variants: list[str],
    rules_path: str | Path | None = None,
) -> ClassificationResult:
    """Backward-compatible wrapper for existing scripts."""
    return classify_robichaux(variants, rules_path)
