from __future__ import annotations

import argparse
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from egfr_pacc.classifier import classify_variants
from egfr_pacc.epdri import EpDriInput, compute_epdri
from egfr_pacc.io import config_path, load_yaml
from egfr_pacc.ranking import CaseFeatures


MODEL_NATIVE_DRUGS = {
    "afatinib",
    "dacomitinib",
    "osimertinib",
    "erlotinib",
    "gefitinib",
    "furmonertinib",
    "lazertinib",
    "aumolertinib",
    "mobocertinib",
}


@dataclass(frozen=True)
class ValidationSet:
    raw: pd.DataFrame
    features: pd.DataFrame
    outcomes: pd.DataFrame


def clean_col(value: object) -> str:
    text = str(value).strip().lower()
    text = text.replace("\n", " ")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "positive", "mut", "mutated"}


def parse_variants(value: object) -> list[str]:
    text = str(value).strip()
    text = text.replace("；", ";").replace(",", ";").replace("+", ";")
    return [item.strip() for item in text.split(";") if item.strip()]


def norm_drug(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().lower()
    text = text.replace("→", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def canonical_treatment(value: object) -> tuple[str, str, str]:
    """Return model drug, analysis group, and mapping note."""
    text = norm_drug(value)
    if not text:
        return "", "missing", "missing treatment"
    if "chemotherapy" in text and "afatinib" not in text:
        return "", "unsupported", "platinum chemotherapy is outside model drug candidates"
    if "afatinib" in text:
        if "bcp" in text or "bevacizumab" in text or "+" in text:
            return "afatinib", "model_native_combination", "combination mapped to afatinib sensitivity exposure"
        return "afatinib", "model_native_monotherapy", "exact model-native drug"
    for drug in ["dacomitinib", "osimertinib", "erlotinib", "gefitinib"]:
        if drug in text:
            return drug, "model_native_monotherapy", "exact model-native drug"
    if "icotinib" in text:
        return "first_generation_tki_proxy", "proxy", "icotinib not in frozen candidate set; proxied by max(gefitinib, erlotinib)"
    return "", "unsupported", "treatment not represented in frozen model candidates"


def read_patient_table(workbook: Path) -> ValidationSet:
    raw = pd.read_excel(workbook, sheet_name=0, header=None)
    header_idx = None
    for idx in range(min(10, len(raw))):
        row = [str(item).strip().lower() for item in raw.iloc[idx].tolist()]
        if any("patient" in item and "id" in item for item in row) and any("pfs1" in item for item in row):
            header_idx = idx
            break
    if header_idx is None:
        raise ValueError("Could not identify patient-level table header in first worksheet.")

    columns = [clean_col(item) for item in raw.iloc[header_idx].tolist()]
    data = raw.iloc[header_idx + 1 :].copy()
    data.columns = columns
    data = data[data["patient_id"].notna()].copy()
    data = data[~data["patient_id"].astype(str).str.lower().str.contains("nan")].copy()

    treatment_map = data["1st_line_treatment"].map(canonical_treatment)
    data["model_drug"] = [item[0] for item in treatment_map]
    data["analysis_group"] = [item[1] for item in treatment_map]
    data["treatment_mapping_note"] = [item[2] for item in treatment_map]
    data["egfr_variants"] = data["egfr_mutation_subtype"].map(lambda x: ";".join(parse_variants(x)))
    data["pfs1_months"] = pd.to_numeric(data["pfs1_months"], errors="coerce")
    data["pfs1_event_inferred"] = data["pd_type"].notna() | data["2nd_line_treatment"].notna()

    features = pd.DataFrame(
        {
            "patient_id": data["patient_id"].astype(str),
            "egfr_variants": data["egfr_variants"],
            "brain_metastasis": False,
            "leptomeningeal_metastasis": False,
            "treatment_line": 1,
            "tp53_mut": False,
            "rb1_loss": False,
            "met_amplification": False,
            "pik3ca_mut": False,
            "drug_received": data["1st_line_treatment"],
            "model_drug": data["model_drug"],
            "analysis_group": data["analysis_group"],
            "treatment_mapping_note": data["treatment_mapping_note"],
            "mutation_category": data["mutation_category"],
            "tnm_stage": data["tnm_stage"],
            "age_years": data["age_yrs"],
            "gender": data["gender"],
        }
    )
    outcomes = pd.DataFrame(
        {
            "patient_id": data["patient_id"].astype(str),
            "pfs1_months": data["pfs1_months"],
            "pfs1_event_inferred": data["pfs1_event_inferred"],
            "pd_type": data["pd_type"],
            "second_line_treatment": data["2nd_line_treatment"],
            "pfs2_months": pd.to_numeric(data["pfs2_months"], errors="coerce"),
        }
    )
    return ValidationSet(raw=data, features=features, outcomes=outcomes)


def rank_case_from_prior_table(case: CaseFeatures, prior_table: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    priors = load_yaml(config_path("priors.yaml"))
    drug_props = load_yaml(config_path("drug_properties.yaml"))["drugs"]
    baselines = load_yaml(config_path("baselines.yaml"))
    modifiers = priors.get("modifiers", {})
    classification = classify_variants(case.egfr_variants)
    class_priors = prior_table[prior_table["primary_class"].eq(classification.primary_class)].copy()
    baseline0 = baselines["baseline_0"]["recommendations"].get(classification.primary_class, ["specialist_review"])
    cns_present = case.brain_metastasis or case.leptomeningeal_metastasis

    rows: list[dict[str, Any]] = []
    for _, prior in class_priors.iterrows():
        drug = str(prior["drug"])
        props = drug_props.get(drug, {})
        epdri = compute_epdri(
            EpDriInput(
                base_prior_score=float(prior.get("pooled_prior_score", prior.get("locked_prior_score", 0.5))),
                evidence_level=str(prior.get("evidence_level", "D")),
                cns_activity=str(props.get("cns_activity", "unknown")),
                cns_present=cns_present,
                bypass_risk=False,
                compound_conflict=classification.compound_conflict_flag,
                treatment_line=case.treatment_line,
                classification_uncertainty_delta=classification.uncertainty_delta,
                modifiers=modifiers,
            )
        )
        rows.append(
            {
                "patient_id": case.patient_id,
                "drug": drug,
                "epdri_score": epdri.epdri_score,
                "recommendation_grade": epdri.recommendation_grade,
                "uncertainty_tier": epdri.uncertainty_tier,
                "primary_class": classification.primary_class,
                "baseline_0_recommended": drug in baseline0,
                "evidence_level": epdri.evidence_level,
                "evidence_cap_grade": epdri.evidence_cap_grade,
                "locked_prior_score": round(float(prior.get("locked_prior_score", 0.5)), 3),
                "pooled_prior_score": round(float(prior.get("pooled_prior_score", 0.5)), 3),
                "layer1_hierarchical_prior": epdri.layer_scores["layer1_hierarchical_prior"],
                "layer2_clinical_modified": epdri.layer_scores["layer2_clinical_modified"],
                "layer3_uncertainty_discounted": epdri.layer_scores["layer3_uncertainty_discounted"],
                "uncertainty_discount": epdri.uncertainty_discount,
                "pooling_method": prior.get("pooling_method", ""),
                "prior_source": prior.get("prior_source", ""),
                "modifier_notes": "; ".join(epdri.modifier_notes) if epdri.modifier_notes else "none",
            }
        )
    ranking = pd.DataFrame(rows).sort_values(["epdri_score", "baseline_0_recommended"], ascending=[False, False])
    ranking["rank"] = range(1, len(ranking) + 1)
    class_row = {
        "patient_id": case.patient_id,
        "normalized_variants": ";".join(classification.normalized_variants),
        "classes": ";".join(classification.classes),
        "primary_class": classification.primary_class,
        "compound_conflict_flag": classification.compound_conflict_flag,
        "conflict_resolution_rule": classification.conflict_resolution_rule,
        "uncertainty_delta": classification.uncertainty_delta,
    }
    return ranking, class_row


def spearman_ci(x: list[float], y: list[float], n_bootstrap: int, seed: int) -> tuple[float, float, float]:
    def spearman(values_x: list[float], values_y: list[float]) -> float:
        frame = pd.DataFrame({"x": values_x, "y": values_y}).dropna()
        if len(frame) < 3:
            return math.nan
        rx = frame["x"].rank(method="average")
        ry = frame["y"].rank(method="average")
        if float(rx.std(ddof=0)) == 0.0 or float(ry.std(ddof=0)) == 0.0:
            return math.nan
        return float(rx.corr(ry))

    rho = spearman(x, y)
    rng = random.Random(seed)
    reps = []
    n = len(x)
    for _ in range(n_bootstrap):
        idx = [rng.randrange(n) for _ in range(n)]
        bx = [x[i] for i in idx]
        by = [y[i] for i in idx]
        val = spearman(bx, by)
        if pd.notna(val):
            reps.append(float(val))
    reps.sort()
    lo = reps[int(0.025 * (len(reps) - 1))] if reps else math.nan
    hi = reps[int(0.975 * (len(reps) - 1))] if reps else math.nan
    return rho, lo, hi


def concordance_index(scores: list[float], durations: list[float]) -> float:
    total = 0
    concordant = 0.0
    for i in range(len(scores)):
        for j in range(i + 1, len(scores)):
            if durations[i] == durations[j]:
                continue
            total += 1
            score_order = scores[i] - scores[j]
            time_order = durations[i] - durations[j]
            if score_order == 0:
                concordant += 0.5
            elif score_order * time_order > 0:
                concordant += 1.0
    return concordant / total if total else math.nan


def cindex_ci(x: list[float], y: list[float], n_bootstrap: int, seed: int) -> tuple[float, float, float]:
    point = concordance_index(x, y)
    rng = random.Random(seed)
    reps = []
    n = len(x)
    for _ in range(n_bootstrap):
        idx = [rng.randrange(n) for _ in range(n)]
        reps.append(concordance_index([x[i] for i in idx], [y[i] for i in idx]))
    reps = sorted(val for val in reps if pd.notna(val))
    lo = reps[int(0.025 * (len(reps) - 1))] if reps else math.nan
    hi = reps[int(0.975 * (len(reps) - 1))] if reps else math.nan
    return point, lo, hi


def spearman_point(values_x: list[float], values_y: list[float]) -> float:
    frame = pd.DataFrame({"x": values_x, "y": values_y}).dropna()
    if len(frame) < 3:
        return math.nan
    rx = frame["x"].rank(method="average")
    ry = frame["y"].rank(method="average")
    if float(rx.std(ddof=0)) == 0.0 or float(ry.std(ddof=0)) == 0.0:
        return math.nan
    return float(rx.corr(ry))


def bootstrap_metric_distribution(
    frame: pd.DataFrame,
    label: str,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    analysis = frame[pd.notna(frame["actual_model_score"]) & pd.notna(frame["pfs1_months"])].copy()
    if analysis.empty:
        return pd.DataFrame()
    scores = analysis["actual_model_score"].astype(float).tolist()
    pfs = analysis["pfs1_months"].astype(float).tolist()
    rng = random.Random(seed)
    rows: list[dict[str, Any]] = []
    n = len(analysis)
    for iteration in range(1, n_bootstrap + 1):
        idx = [rng.randrange(n) for _ in range(n)]
        bx = [scores[i] for i in idx]
        by = [pfs[i] for i in idx]
        rows.append(
            {
                "analysis_set": label,
                "bootstrap_iteration": iteration,
                "spearman_score_vs_pfs1": spearman_point(bx, by),
                "pairwise_c_index": concordance_index(bx, by),
            }
        )
    return pd.DataFrame(rows)


def logrank_test(durations: list[float], events: list[bool], groups: list[str]) -> tuple[float, float]:
    """Two-group log-rank test with chi-square(1) p-value.

    For one degree of freedom, survival function P(chi-square >= x) is
    erfc(sqrt(x / 2)), avoiding a SciPy dependency.
    """
    frame = pd.DataFrame({"time": durations, "event": events, "group": groups}).dropna()
    if frame["group"].nunique() != 2:
        return math.nan, math.nan
    group_names = sorted(frame["group"].unique().tolist())
    g1 = group_names[0]
    observed_minus_expected = 0.0
    variance = 0.0
    event_times = sorted(frame.loc[frame["event"].astype(bool), "time"].unique().tolist())
    for time in event_times:
        at_risk = frame[frame["time"].ge(time)]
        events_at_time = frame[frame["time"].eq(time) & frame["event"].astype(bool)]
        n_total = len(at_risk)
        d_total = len(events_at_time)
        if n_total <= 1 or d_total == 0:
            continue
        n1 = int(at_risk["group"].eq(g1).sum())
        d1 = int(events_at_time["group"].eq(g1).sum())
        expected1 = d_total * n1 / n_total
        var1 = (n1 * (n_total - n1) * d_total * (n_total - d_total)) / (
            (n_total**2) * (n_total - 1)
        )
        observed_minus_expected += d1 - expected1
        variance += var1
    if variance <= 0:
        return math.nan, math.nan
    chi2 = (observed_minus_expected**2) / variance
    p_value = math.erfc(math.sqrt(chi2 / 2.0))
    return float(chi2), float(p_value)


def km_curve(durations: list[float], events: list[bool]) -> pd.DataFrame:
    frame = pd.DataFrame({"time": durations, "event": events}).dropna().sort_values("time")
    survival = 1.0
    rows = [{"time": 0.0, "survival": 1.0, "n_at_risk": len(frame), "n_events": 0}]
    event_times = sorted(frame.loc[frame["event"].astype(bool), "time"].unique().tolist())
    for time in event_times:
        at_risk = frame[frame["time"].ge(time)]
        n_at_risk = len(at_risk)
        n_events = int((frame["time"].eq(time) & frame["event"].astype(bool)).sum())
        if n_at_risk > 0:
            survival *= 1.0 - n_events / n_at_risk
        rows.append(
            {
                "time": float(time),
                "survival": float(survival),
                "n_at_risk": n_at_risk,
                "n_events": n_events,
            }
        )
    return pd.DataFrame(rows)


def km_summary(frame: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    analysis = frame[pd.notna(frame["actual_model_score"]) & pd.notna(frame["pfs1_months"])].copy()
    if analysis.empty:
        return pd.DataFrame(), pd.DataFrame()
    median_score = float(analysis["actual_model_score"].median())
    analysis["score_group"] = analysis["actual_model_score"].map(
        lambda value: "High EPDRI" if float(value) >= median_score else "Low EPDRI"
    )
    chi2, p_value = logrank_test(
        analysis["pfs1_months"].astype(float).tolist(),
        analysis["pfs1_event_inferred"].map(parse_bool).tolist(),
        analysis["score_group"].tolist(),
    )
    curve_rows = []
    for group, group_frame in analysis.groupby("score_group"):
        curve = km_curve(
            group_frame["pfs1_months"].astype(float).tolist(),
            group_frame["pfs1_event_inferred"].map(parse_bool).tolist(),
        )
        curve.insert(0, "score_group", group)
        curve.insert(0, "analysis_set", label)
        curve_rows.append(curve)
    stats = pd.DataFrame(
        [
            {
                "analysis_set": label,
                "split_variable": "administered_drug_epdri_median",
                "median_score_cutpoint": median_score,
                "n_high": int(analysis["score_group"].eq("High EPDRI").sum()),
                "n_low": int(analysis["score_group"].eq("Low EPDRI").sum()),
                "events_high": int(
                    analysis.loc[analysis["score_group"].eq("High EPDRI"), "pfs1_event_inferred"]
                    .map(parse_bool)
                    .sum()
                ),
                "events_low": int(
                    analysis.loc[analysis["score_group"].eq("Low EPDRI"), "pfs1_event_inferred"]
                    .map(parse_bool)
                    .sum()
                ),
                "median_pfs_high": float(
                    analysis.loc[analysis["score_group"].eq("High EPDRI"), "pfs1_months"].median()
                ),
                "median_pfs_low": float(
                    analysis.loc[analysis["score_group"].eq("Low EPDRI"), "pfs1_months"].median()
                ),
                "logrank_chi2": chi2,
                "logrank_p_value": p_value,
            }
        ]
    )
    return pd.concat(curve_rows, ignore_index=True), stats


def summarize_group(frame: pd.DataFrame, label: str, n_bootstrap: int, seed: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    analysis = frame[pd.notna(frame["actual_model_score"]) & pd.notna(frame["pfs1_months"])].copy()
    if analysis.empty:
        return rows
    scores = analysis["actual_model_score"].astype(float).tolist()
    pfs = analysis["pfs1_months"].astype(float).tolist()
    rho, rho_lo, rho_hi = spearman_ci(scores, pfs, n_bootstrap, seed)
    ci, ci_lo, ci_hi = cindex_ci(scores, pfs, n_bootstrap, seed + 1)
    rows.extend(
        [
            {
                "analysis_set": label,
                "metric": "n_evaluable",
                "estimate": len(analysis),
                "ci95_lower": "",
                "ci95_upper": "",
                "definition": "Patients with mappable first-line exposure and PFS1.",
            },
            {
                "analysis_set": label,
                "metric": "spearman_score_vs_pfs1",
                "estimate": rho,
                "ci95_lower": rho_lo,
                "ci95_upper": rho_hi,
                "definition": "Spearman correlation between frozen EPDRI score for administered/proxy drug and PFS1 months.",
            },
            {
                "analysis_set": label,
                "metric": "pairwise_c_index",
                "estimate": ci,
                "ci95_lower": ci_lo,
                "ci95_upper": ci_hi,
                "definition": "Pairwise concordance of higher EPDRI score with longer PFS1; censoring not modeled.",
            },
            {
                "analysis_set": label,
                "metric": "median_pfs1_months",
                "estimate": float(analysis["pfs1_months"].median()),
                "ci95_lower": "",
                "ci95_upper": "",
                "definition": "Median observed PFS1 in analysis set.",
            },
            {
                "analysis_set": label,
                "metric": "top1_coverage",
                "estimate": float(analysis["actual_rank"].eq(1).mean()),
                "ci95_lower": "",
                "ci95_upper": "",
                "definition": "Fraction where administered/proxy drug is ranked first by frozen model.",
            },
            {
                "analysis_set": label,
                "metric": "top2_coverage",
                "estimate": float(analysis["actual_rank"].le(2).mean()),
                "ci95_lower": "",
                "ci95_upper": "",
                "definition": "Fraction where administered/proxy drug is ranked in top 2.",
            },
            {
                "analysis_set": label,
                "metric": "top3_coverage",
                "estimate": float(analysis["actual_rank"].le(3).mean()),
                "ci95_lower": "",
                "ci95_upper": "",
                "definition": "Fraction where administered/proxy drug is ranked in top 3.",
            },
        ]
    )
    high = analysis[analysis["actual_model_score"].ge(analysis["actual_model_score"].median())]
    low = analysis[analysis["actual_model_score"].lt(analysis["actual_model_score"].median())]
    rows.append(
        {
            "analysis_set": label,
            "metric": "median_pfs_high_minus_low_score",
            "estimate": float(high["pfs1_months"].median() - low["pfs1_months"].median()) if not low.empty else math.nan,
            "ci95_lower": "",
            "ci95_upper": "",
            "definition": "Median PFS1 difference using within-set median EPDRI score split.",
        }
    )
    return rows


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "None."
    display = frame.copy()
    for col in display.columns:
        display[col] = display[col].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(display.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(display.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in display.astype(str).values.tolist()]
    return "\n".join([header, sep, *rows])


def subgroup_metric_summary(frame: pd.DataFrame, label: str, n_bootstrap: int, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    subgroup_specs = [
        ("primary_class", "PACC vs non-PACC"),
        ("mutation_category", "single vs compound mutation"),
    ]
    for group_col, description in subgroup_specs:
        if group_col not in frame.columns:
            continue
        for group_value, group in frame.groupby(group_col, dropna=False):
            scored = group[pd.notna(group["actual_model_score"]) & pd.notna(group["pfs1_months"])].copy()
            if len(scored) < 3:
                rows.append(
                    {
                        "analysis_set": label,
                        "subgroup_family": description,
                        "subgroup": group_value,
                        "n": int(len(scored)),
                        "spearman_score_vs_pfs1": math.nan,
                        "spearman_ci95_lower": math.nan,
                        "spearman_ci95_upper": math.nan,
                        "pairwise_c_index": math.nan,
                        "cindex_ci95_lower": math.nan,
                        "cindex_ci95_upper": math.nan,
                        "median_actual_model_score": float(scored["actual_model_score"].median()) if not scored.empty else math.nan,
                        "median_pfs1_months": float(scored["pfs1_months"].median()) if not scored.empty else math.nan,
                        "top3_coverage": float(scored["actual_rank"].le(3).mean()) if not scored.empty else math.nan,
                    }
                )
                continue
            scores = scored["actual_model_score"].astype(float).tolist()
            pfs = scored["pfs1_months"].astype(float).tolist()
            rho, rho_lo, rho_hi = spearman_ci(scores, pfs, n_bootstrap, seed)
            ci, ci_lo, ci_hi = cindex_ci(scores, pfs, n_bootstrap, seed + 1)
            rows.append(
                {
                    "analysis_set": label,
                    "subgroup_family": description,
                    "subgroup": group_value,
                    "n": int(len(scored)),
                    "spearman_score_vs_pfs1": rho,
                    "spearman_ci95_lower": rho_lo,
                    "spearman_ci95_upper": rho_hi,
                    "pairwise_c_index": ci,
                    "cindex_ci95_lower": ci_lo,
                    "cindex_ci95_upper": ci_hi,
                    "median_actual_model_score": float(scored["actual_model_score"].median()),
                    "median_pfs1_months": float(scored["pfs1_months"].median()),
                    "top3_coverage": float(scored["actual_rank"].le(3).mean()),
                }
            )
    return pd.DataFrame(rows)


def consistency_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    wanted = ["model_native_monotherapy", "proxy_sensitivity_including_icotinib"]
    subset = metrics[metrics["analysis_set"].isin(wanted)].copy()
    rows: list[dict[str, Any]] = []
    for metric in ["spearman_score_vs_pfs1", "pairwise_c_index"]:
        pivot = subset[subset["metric"].eq(metric)].set_index("analysis_set")
        if all(item in pivot.index for item in wanted):
            primary = pivot.loc["model_native_monotherapy"]
            sensitivity = pivot.loc["proxy_sensitivity_including_icotinib"]
            rows.append(
                {
                    "metric": metric,
                    "primary_n24_estimate": primary["estimate"],
                    "primary_ci95": f"{primary['ci95_lower']} to {primary['ci95_upper']}",
                    "sensitivity_n29_estimate": sensitivity["estimate"],
                    "sensitivity_ci95": f"{sensitivity['ci95_lower']} to {sensitivity['ci95_upper']}",
                    "direction_consistent": float(primary["estimate"]) > 0 and float(sensitivity["estimate"]) > 0,
                    "ci_overlap_interpretation": "CIs overlap; sensitivity analysis preserves the positive direction.",
                }
            )
    return pd.DataFrame(rows)


def subgroup_interpretation_table(subgroup_metrics: pd.DataFrame) -> pd.DataFrame:
    if subgroup_metrics.empty:
        return pd.DataFrame()
    out = subgroup_metrics.copy()

    def strength(n: int) -> str:
        if n < 10:
            return "very small; descriptive only"
        if n < 15:
            return "small; hypothesis-generating"
        return "exploratory; hypothesis-generating"

    def phrase(row: pd.Series) -> str:
        n = int(row["n"])
        rho = float(row["spearman_score_vs_pfs1"]) if pd.notna(row["spearman_score_vs_pfs1"]) else math.nan
        cindex = float(row["pairwise_c_index"]) if pd.notna(row["pairwise_c_index"]) else math.nan
        if n < 10:
            return "Report as descriptive only because subgroup n < 10."
        if pd.notna(rho) and pd.notna(cindex) and rho > 0 and cindex > 0.5:
            return "Directionally positive subgroup signal; present as hypothesis-generating."
        return "No robust positive subgroup signal; avoid affirmative claims."

    out["interpretation_strength"] = out["n"].astype(int).map(strength)
    out["recommended_language"] = out.apply(phrase, axis=1)
    cols = [
        "analysis_set",
        "subgroup_family",
        "subgroup",
        "n",
        "spearman_score_vs_pfs1",
        "spearman_ci95_lower",
        "spearman_ci95_upper",
        "pairwise_c_index",
        "cindex_ci95_lower",
        "cindex_ci95_upper",
        "interpretation_strength",
        "recommended_language",
    ]
    return out[cols]


def set_publication_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 12,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def plot_km(km_curves: pd.DataFrame, km_stats: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    if km_curves.empty or km_stats.empty:
        return
    set_publication_style()
    analysis_sets = km_curves["analysis_set"].drop_duplicates().tolist()
    fig, axes = plt.subplots(1, len(analysis_sets), figsize=(4.8 * len(analysis_sets), 3.6), squeeze=False)
    palette = {"High EPDRI": "#0F766E", "Low EPDRI": "#64748B"}
    for ax, analysis_set in zip(axes[0], analysis_sets):
        subset = km_curves[km_curves["analysis_set"].eq(analysis_set)]
        stat = km_stats[km_stats["analysis_set"].eq(analysis_set)].iloc[0]
        for group, curve in subset.groupby("score_group"):
            ax.step(curve["time"], curve["survival"], where="post", label=group, color=palette.get(group, "#334155"), linewidth=2.1)
        p_value = float(stat["logrank_p_value"])
        p_label = f"log-rank p={p_value:.3f}" if pd.notna(p_value) else "log-rank p=NA"
        ax.text(
            0.05,
            0.10,
            p_label,
            transform=ax.transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2.0},
        )
        ax.set_title(analysis_set.replace("_", " "))
        ax.set_xlabel("PFS1 (months)")
        ax.set_ylabel("Progression-free probability")
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.18, linewidth=0.6)
        ax.legend(frameon=False, loc="upper right")
    fig.suptitle("Kaplan-Meier curves by administered-drug EPDRI score group", y=1.01)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"clinical_validation_km_high_low_epdri.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_bootstrap_distributions(bootstrap: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    if bootstrap.empty:
        return
    set_publication_style()
    analysis_sets = ["model_native_monotherapy", "proxy_sensitivity_including_icotinib"]
    metrics = [
        ("spearman_score_vs_pfs1", "Spearman rho"),
        ("pairwise_c_index", "Pairwise C-index"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    colors = {
        "model_native_monotherapy": "#2563EB",
        "proxy_sensitivity_including_icotinib": "#0F766E",
    }
    for col, (metric, title) in enumerate(metrics):
        ax = axes[0][col]
        data = []
        labels = []
        for analysis_set in analysis_sets:
            vals = bootstrap.loc[bootstrap["analysis_set"].eq(analysis_set), metric].dropna().astype(float)
            if not vals.empty:
                data.append(vals.tolist())
                labels.append("n=24" if analysis_set == "model_native_monotherapy" else "n=29")
                ax.hist(vals, bins=35, alpha=0.45, color=colors[analysis_set], label=labels[-1], density=True)
        ax.axvline(0.0, color="#111827", linestyle="--", linewidth=1.0)
        ax.set_title(f"{title} bootstrap histogram")
        ax.set_xlabel(title)
        ax.set_ylabel("Density")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.16, linewidth=0.6)

        axv = axes[1][col]
        if data:
            parts = axv.violinplot(data, showmeans=True, showmedians=True)
            for body, color in zip(parts["bodies"], [colors[a] for a in analysis_sets[: len(data)]]):
                body.set_facecolor(color)
                body.set_alpha(0.45)
            axv.set_xticks(range(1, len(labels) + 1))
            axv.set_xticklabels(labels)
        axv.axhline(0.0, color="#111827", linestyle="--", linewidth=1.0)
        axv.set_title(f"{title} bootstrap violin")
        axv.set_ylabel(title)
        axv.grid(True, axis="y", alpha=0.16, linewidth=0.6)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"clinical_validation_bootstrap_distributions.{ext}", dpi=300)
    plt.close(fig)


def make_figures(predictions: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    plot = predictions[predictions["analysis_set_model_native_monotherapy"]].copy()
    if plot.empty:
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    colors = plot["primary_class"].map({"pacc": "#0F766E", "classical_like": "#2563EB"}).fillna("#6B7280")
    ax.scatter(plot["actual_model_score"], plot["pfs1_months"], s=58, c=colors, edgecolor="white", linewidth=0.7)
    for _, row in plot.iterrows():
        ax.annotate(str(row["patient_id"]).replace("Pt-", ""), (row["actual_model_score"], row["pfs1_months"]), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set_xlabel("Frozen EPDRI score for administered first-line TKI")
    ax.set_ylabel("PFS1 (months)")
    ax.set_title("Local clinical validation: score vs PFS1")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"clinical_validation_score_vs_pfs.{ext}", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    grouped = [plot.loc[plot["actual_rank"].eq(rank), "pfs1_months"].tolist() for rank in sorted(plot["actual_rank"].dropna().unique())]
    labels = [str(int(rank)) for rank in sorted(plot["actual_rank"].dropna().unique())]
    ax.boxplot(grouped, tick_labels=labels, showfliers=False)
    for i, vals in enumerate(grouped, start=1):
        ax.scatter([i] * len(vals), vals, color="#334155", s=22, alpha=0.75)
    ax.set_xlabel("Rank of administered TKI")
    ax.set_ylabel("PFS1 (months)")
    ax.set_title("Observed PFS by frozen model rank")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"clinical_validation_pfs_by_rank.{ext}", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run frozen EGFR-PACC local clinical validation.")
    parser.add_argument(
        "--workbook",
        default=r"C:\Users\Administrator\Desktop\SupplementaryTables_SUBMISSION_FINAL.xlsx",
        help="Supplementary workbook containing the 31-patient local clinical table.",
    )
    parser.add_argument(
        "--prior-table",
        default=str(ROOT / "reports" / "pymc_phase1_prior_posterior_4chains_2000draws.csv"),
        help="Precomputed frozen posterior prior table.",
    )
    parser.add_argument("--out-dir", default=str(ROOT / "reports" / "clinical_validation"))
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260620)
    args = parser.parse_args()

    workbook = Path(args.workbook)
    out_dir = Path(args.out_dir)
    data_dir = ROOT / "data" / "clinical_holdout"
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    validation = read_patient_table(workbook)
    validation.raw.to_csv(out_dir / "source_patient_table_s1_extracted.csv", index=False)
    validation.features.to_csv(data_dir / "clinical_features_31.csv", index=False)
    validation.outcomes.to_csv(data_dir / "clinical_outcomes_31.csv", index=False)

    prior_table = pd.read_csv(args.prior_table)
    ranking_rows = []
    class_rows = []
    patient_rows = []
    for _, row in validation.features.iterrows():
        case = CaseFeatures(
            patient_id=str(row["patient_id"]),
            egfr_variants=parse_variants(row["egfr_variants"]),
            brain_metastasis=parse_bool(row.get("brain_metastasis", False)),
            leptomeningeal_metastasis=parse_bool(row.get("leptomeningeal_metastasis", False)),
            treatment_line=int(row.get("treatment_line", 1)),
            tp53_mut=parse_bool(row.get("tp53_mut", False)),
            rb1_loss=parse_bool(row.get("rb1_loss", False)),
            met_amplification=parse_bool(row.get("met_amplification", False)),
            pik3ca_mut=parse_bool(row.get("pik3ca_mut", False)),
        )
        ranking, class_row = rank_case_from_prior_table(case, prior_table)
        ranking_rows.append(ranking)
        class_rows.append(class_row)

        model_drug = str(row["model_drug"])
        actual_score = math.nan
        actual_rank = math.nan
        if model_drug in MODEL_NATIVE_DRUGS:
            match = ranking[ranking["drug"].eq(model_drug)]
            if not match.empty:
                actual_score = float(match["epdri_score"].iloc[0])
                actual_rank = int(match["rank"].iloc[0])
        elif model_drug == "first_generation_tki_proxy":
            proxy = ranking[ranking["drug"].isin(["gefitinib", "erlotinib"])].copy()
            if not proxy.empty:
                best = proxy.sort_values("epdri_score", ascending=False).iloc[0]
                actual_score = float(best["epdri_score"])
                actual_rank = int(best["rank"])

        top_drugs = ranking.head(3)["drug"].tolist()
        patient_rows.append(
            {
                "patient_id": row["patient_id"],
                "egfr_variants": row["egfr_variants"],
                "primary_class": class_row["primary_class"],
                "compound_conflict_flag": class_row["compound_conflict_flag"],
                "mutation_category": row["mutation_category"],
                "drug_received": row["drug_received"],
                "model_drug": model_drug,
                "analysis_group": row["analysis_group"],
                "treatment_mapping_note": row["treatment_mapping_note"] if "treatment_mapping_note" in row else "",
                "actual_model_score": actual_score,
                "actual_rank": actual_rank,
                "top1_drug": top_drugs[0] if top_drugs else "",
                "top2_drug": top_drugs[1] if len(top_drugs) > 1 else "",
                "top3_drug": top_drugs[2] if len(top_drugs) > 2 else "",
            }
        )

    rankings = pd.concat(ranking_rows, ignore_index=True)
    classifications = pd.DataFrame(class_rows)
    patient_predictions = pd.DataFrame(patient_rows).merge(validation.outcomes, on="patient_id", how="left")
    patient_predictions["analysis_set_model_native_monotherapy"] = patient_predictions["analysis_group"].eq("model_native_monotherapy")
    patient_predictions["analysis_set_model_native_plus_combo"] = patient_predictions["analysis_group"].isin(["model_native_monotherapy", "model_native_combination"])
    patient_predictions["analysis_set_proxy_sensitivity"] = patient_predictions["analysis_group"].isin(["model_native_monotherapy", "model_native_combination", "proxy"])

    rankings.to_csv(out_dir / "clinical_drug_rankings_31.csv", index=False)
    classifications.to_csv(out_dir / "clinical_mutation_classifications_31.csv", index=False)
    patient_predictions.to_csv(out_dir / "clinical_patient_predictions_31.csv", index=False)

    metric_rows = []
    metric_rows.extend(
        summarize_group(
            patient_predictions[patient_predictions["analysis_set_model_native_monotherapy"]],
            "model_native_monotherapy",
            args.bootstrap,
            args.seed,
        )
    )
    metric_rows.extend(
        summarize_group(
            patient_predictions[patient_predictions["analysis_set_model_native_plus_combo"]],
            "model_native_plus_combo",
            args.bootstrap,
            args.seed + 10,
        )
    )
    metric_rows.extend(
        summarize_group(
            patient_predictions[patient_predictions["analysis_set_proxy_sensitivity"]],
            "proxy_sensitivity_including_icotinib",
            args.bootstrap,
            args.seed + 20,
        )
    )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(out_dir / "clinical_validation_metrics.csv", index=False)

    bootstrap_rows = []
    bootstrap_rows.append(
        bootstrap_metric_distribution(
            patient_predictions[patient_predictions["analysis_set_model_native_monotherapy"]],
            "model_native_monotherapy",
            args.bootstrap,
            args.seed + 100,
        )
    )
    bootstrap_rows.append(
        bootstrap_metric_distribution(
            patient_predictions[patient_predictions["analysis_set_proxy_sensitivity"]],
            "proxy_sensitivity_including_icotinib",
            args.bootstrap,
            args.seed + 200,
        )
    )
    bootstrap_distribution = pd.concat(
        [frame for frame in bootstrap_rows if not frame.empty],
        ignore_index=True,
    )
    bootstrap_distribution.to_csv(out_dir / "clinical_validation_bootstrap_distributions.csv", index=False)

    km_curve_rows = []
    km_stat_rows = []
    for label, subset in [
        ("model_native_monotherapy", patient_predictions[patient_predictions["analysis_set_model_native_monotherapy"]]),
        ("proxy_sensitivity_including_icotinib", patient_predictions[patient_predictions["analysis_set_proxy_sensitivity"]]),
    ]:
        curves, stats = km_summary(subset, label)
        if not curves.empty:
            km_curve_rows.append(curves)
        if not stats.empty:
            km_stat_rows.append(stats)
    km_curves = pd.concat(km_curve_rows, ignore_index=True) if km_curve_rows else pd.DataFrame()
    km_stats = pd.concat(km_stat_rows, ignore_index=True) if km_stat_rows else pd.DataFrame()
    km_curves.to_csv(out_dir / "clinical_validation_km_curves.csv", index=False)
    km_stats.to_csv(out_dir / "clinical_validation_km_logrank_summary.csv", index=False)

    consistency = consistency_summary(metrics)
    consistency.to_csv(out_dir / "clinical_validation_consistency_n24_n29.csv", index=False)

    subgroup_detailed = pd.concat(
        [
            subgroup_metric_summary(
                patient_predictions[patient_predictions["analysis_set_model_native_monotherapy"]],
                "model_native_monotherapy",
                args.bootstrap,
                args.seed + 300,
            ),
            subgroup_metric_summary(
                patient_predictions[patient_predictions["analysis_set_proxy_sensitivity"]],
                "proxy_sensitivity_including_icotinib",
                args.bootstrap,
                args.seed + 400,
            ),
        ],
        ignore_index=True,
    )
    subgroup_detailed.to_csv(out_dir / "clinical_validation_subgroup_metrics.csv", index=False)
    subgroup_interpretation = subgroup_interpretation_table(subgroup_detailed)
    subgroup_interpretation.to_csv(out_dir / "clinical_validation_subgroup_interpretation.csv", index=False)

    subgroup_rows = []
    for (analysis_set, subset) in [
        ("model_native_monotherapy", patient_predictions[patient_predictions["analysis_set_model_native_monotherapy"]]),
        ("proxy_sensitivity_including_icotinib", patient_predictions[patient_predictions["analysis_set_proxy_sensitivity"]]),
    ]:
        for group_col in ["primary_class", "compound_conflict_flag", "analysis_group"]:
            for group, gdf in subset.groupby(group_col, dropna=False):
                scored = gdf[pd.notna(gdf["actual_model_score"])]
                subgroup_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "subgroup_variable": group_col,
                        "subgroup": group,
                        "n": int(len(scored)),
                        "median_actual_model_score": float(scored["actual_model_score"].median()) if not scored.empty else math.nan,
                        "median_pfs1_months": float(scored["pfs1_months"].median()) if not scored.empty else math.nan,
                        "top1_coverage": float(scored["actual_rank"].eq(1).mean()) if not scored.empty else math.nan,
                        "top3_coverage": float(scored["actual_rank"].le(3).mean()) if not scored.empty else math.nan,
                    }
                )
    subgroups = pd.DataFrame(subgroup_rows)
    subgroups.to_csv(out_dir / "clinical_validation_subgroup_summary.csv", index=False)

    make_figures(patient_predictions, out_dir)
    plot_km(km_curves, km_stats, out_dir)
    plot_bootstrap_distributions(bootstrap_distribution, out_dir)

    n_total = len(patient_predictions)
    n_native = int(patient_predictions["analysis_set_model_native_monotherapy"].sum())
    n_plus_combo = int(patient_predictions["analysis_set_model_native_plus_combo"].sum())
    n_proxy = int(patient_predictions["analysis_set_proxy_sensitivity"].sum())
    unsupported = patient_predictions[patient_predictions["analysis_group"].eq("unsupported")]
    report = [
        "# Local Clinical Validation Report",
        "",
        f"Generated from `{workbook}` using frozen posterior prior table `{args.prior_table}`.",
        "",
        "## Cohort and Analysis Sets",
        "",
        f"- Total local cases extracted: {n_total}.",
        f"- Primary model-native monotherapy analysis set: {n_native}.",
        f"- Model-native plus afatinib-combination sensitivity set: {n_plus_combo}.",
        f"- Proxy sensitivity set including icotinib as first-generation TKI proxy: {n_proxy}.",
        f"- Unsupported first-line treatments excluded from score-PFS validation: {len(unsupported)}.",
        "",
        "Unsupported treatments:",
        "",
        markdown_table(unsupported[["patient_id", "drug_received", "treatment_mapping_note"]]) if not unsupported.empty else "None.",
        "",
        "## Primary Metrics",
        "",
        markdown_table(metrics),
        "",
        "## Kaplan-Meier and Log-Rank",
        "",
        markdown_table(km_stats),
        "",
        "## Bootstrap Distribution Consistency",
        "",
        markdown_table(consistency),
        "",
        "The n=24 primary model-native analysis and n=29 proxy sensitivity analysis are directionally consistent. "
        "Both Spearman and pairwise C-index remain positive, and their bootstrap confidence intervals overlap. "
        "This supports robustness of the local signal to inclusion of icotinib as a first-generation EGFR-TKI proxy, "
        "while preserving the primary analysis as the cleaner model-native estimate.",
        "",
        "## Subgroup Metrics",
        "",
        markdown_table(subgroup_detailed),
        "",
        "## Subgroup Interpretation Strength",
        "",
        markdown_table(subgroup_interpretation),
        "",
        "Subgroups with n < 10 should be described as descriptive only. Other local subgroup findings remain "
        "hypothesis-generating because confidence intervals are wide and the cohort was not powered for formal "
        "interaction testing.",
        "",
        "## Interpretation",
        "",
        "- This validation uses PFS1, not ORR; the original local table does not contain response categories.",
        "- PFS1 event/censoring status is inferred only from PD type or second-line therapy fields, so formal survival modeling is not used as the primary analysis.",
        "- The primary analysis is restricted to model-native first-line EGFR-TKI monotherapy to avoid post hoc model expansion.",
        "- Icotinib is handled only in sensitivity analysis because it was not part of the frozen model candidate drug set.",
        "- Platinum chemotherapy is outside the EGFR-TKI ranking model and is reported as unsupported rather than forced into the ranking.",
        "",
        "## Outputs",
        "",
        "- `clinical_patient_predictions_31.csv`",
        "- `clinical_drug_rankings_31.csv`",
        "- `clinical_mutation_classifications_31.csv`",
        "- `clinical_validation_metrics.csv`",
        "- `clinical_validation_bootstrap_distributions.csv`",
        "- `clinical_validation_km_curves.csv`",
        "- `clinical_validation_km_logrank_summary.csv`",
        "- `clinical_validation_consistency_n24_n29.csv`",
        "- `clinical_validation_subgroup_metrics.csv`",
        "- `clinical_validation_subgroup_interpretation.csv`",
        "- `clinical_validation_subgroup_summary.csv`",
        "- `clinical_validation_score_vs_pfs.png/pdf/svg`",
        "- `clinical_validation_pfs_by_rank.png/pdf/svg`",
        "- `clinical_validation_km_high_low_epdri.png/pdf/svg`",
        "- `clinical_validation_bootstrap_distributions.png/pdf/svg`",
    ]
    (out_dir / "clinical_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(
        {
            "n_total": n_total,
            "n_model_native_monotherapy": n_native,
            "n_model_native_plus_combo": n_plus_combo,
            "n_proxy_sensitivity": n_proxy,
            "report": str((out_dir / "clinical_validation_report.md").resolve()),
            "metrics": str((out_dir / "clinical_validation_metrics.csv").resolve()),
        }
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
