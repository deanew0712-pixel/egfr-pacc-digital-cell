from __future__ import annotations

import argparse
import math
import random
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from egfr_pacc.classifier import classify_robichaux


LEVEL_SCORE = {
    "LEVEL_1": 1.00,
    "LEVEL_2": 0.85,
    "LEVEL_3A": 0.65,
    "LEVEL_3B": 0.60,
    "LEVEL_4": 0.45,
    "LEVEL_R1": 0.05,
    "LEVEL_R2": 0.00,
}

DRUG_ALIASES = {
    "second_generation_tki": ["afatinib", "dacomitinib"],
    "first_generation_tki": ["gefitinib", "erlotinib"],
    "third_generation_tki": ["osimertinib", "furmonertinib", "lazertinib", "aumolertinib"],
}


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def norm_drug(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_")


def display_drug(value: str) -> str:
    return value.replace("_", " ")


def bucket_to_variants(bucket: object) -> list[str]:
    value = str(bucket).strip()
    lower = value.lower()
    mapping = {
        "g719x": ["G719S"],
        "s768i": ["S768I"],
        "l861q": ["L861Q"],
        "e709x": ["E709A"],
        "pacc": ["G719S"],
        "major_uncommon": ["G719S", "S768I", "L861Q"],
        "compound_uncommon": ["G719S", "L861Q"],
        "classical_like_or_pacc": ["G719S", "L861Q"],
        "uncommon_non_ex20ins": ["G719S", "S768I", "L861Q"],
        "solitary_uncommon": ["G719S"],
        "sensitizing_uncommon_non_ex20ins_non_t790m": ["G719S", "S768I", "L861Q"],
        "common_plus_uncommon": ["L858R", "G719S"],
        "uncommon_plus_uncommon": ["G719S", "S768I"],
        "brain_metastasis_uncommon": ["G719S"],
        "other_uncommon": ["E709A"],
    }
    if lower in mapping:
        return mapping[lower]
    if lower.endswith("x") and len(value) >= 4:
        return [value[:-1].upper() + "S"]
    return [value]


@lru_cache(maxsize=None)
def primary_class_for_bucket(bucket: object) -> str:
    return classify_robichaux(bucket_to_variants(bucket)).primary_class


def expand_drug(drug: object) -> list[str]:
    key = norm_drug(drug)
    return DRUG_ALIASES.get(key, [key])


def weighted_mean_orr(frame: pd.DataFrame) -> float | None:
    subset = frame[pd.notna(frame["orr_percent"])].copy()
    if subset.empty:
        return None
    weights = pd.to_numeric(subset["n"], errors="coerce").fillna(1.0).clip(lower=1.0)
    values = pd.to_numeric(subset["orr_percent"], errors="coerce") / 100.0
    ok = pd.notna(values)
    if not ok.any():
        return None
    return float((values[ok] * weights[ok]).sum() / weights[ok].sum())


def auc_score(y_true: list[int], scores: list[float]) -> float:
    pairs = [(float(score), int(label)) for score, label in zip(scores, y_true)]
    positives = sum(label for _, label in pairs)
    negatives = len(pairs) - positives
    if positives == 0 or negatives == 0:
        return float("nan")
    pairs.sort(key=lambda item: item[0])
    ranks: list[float] = [0.0] * len(pairs)
    i = 0
    while i < len(pairs):
        j = i + 1
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    rank_sum_pos = sum(rank for rank, (_, label) in zip(ranks, pairs) if label == 1)
    return (rank_sum_pos - positives * (positives + 1) / 2.0) / (positives * negatives)


def base_prior_score(priors: dict[str, Any], class_name: str, drug: str) -> float:
    entry = priors.get("class_drug_priors", {}).get(class_name, {}).get(drug, {})
    return float(entry.get("score", 0.05))


def literature_model_score(
    train: pd.DataFrame,
    priors: dict[str, Any],
    heldout_bucket: str,
    class_name: str,
    candidate_drug: str,
    exact_floor: bool = True,
    config_fallback: bool = True,
) -> tuple[float, str]:
    base = base_prior_score(priors, class_name, candidate_drug) if config_fallback else 0.05
    train = train.copy()
    train["expanded_drugs"] = train["drug"].map(expand_drug)
    train = train.explode("expanded_drugs")
    train["expanded_drugs"] = train["expanded_drugs"].map(norm_drug)
    train["derived_class"] = train["mutation_bucket"].map(primary_class_for_bucket)

    exact = train[
        train["mutation_bucket"].astype(str).str.lower().eq(str(heldout_bucket).lower())
        & train["expanded_drugs"].eq(candidate_drug)
    ]
    exact_orr = weighted_mean_orr(exact)
    if exact_orr is not None:
        blended = 0.25 * base + 0.75 * exact_orr
        score = max(base, blended) if exact_floor else blended
        return max(0.01, min(0.99, score)), "exact_mutation_drug_loo_literature"

    class_drug = train[
        train["derived_class"].eq(class_name)
        & train["expanded_drugs"].eq(candidate_drug)
    ]
    class_orr = weighted_mean_orr(class_drug)
    if class_orr is not None:
        blended = 0.50 * base + 0.50 * class_orr
        score = max(base, blended) if exact_floor else blended
        return max(0.01, min(0.99, score)), "class_drug_loo_literature"

    return max(0.01, min(0.99, base)), "locked_config_prior_fallback"


def baseline0_score(baselines: dict[str, Any], class_name: str, drug: str) -> tuple[float, bool]:
    recs = [norm_drug(item) for item in baselines["baseline_0"]["recommendations"].get(class_name, [])]
    if drug not in recs:
        return 0.05, False
    rank = recs.index(drug)
    return max(0.10, 1.0 - 0.05 * rank), True


def oncokb_level_for_bucket(oncokb: pd.DataFrame, bucket: str, drug: str) -> tuple[str, float]:
    variants = [v.upper() for v in bucket_to_variants(bucket)]
    if oncokb.empty:
        return "not_available", 0.0
    matrix = oncokb.copy()
    matrix["alteration_norm"] = matrix["alteration"].astype(str).str.upper()
    matrix["drug_norm"] = matrix["drug"].map(norm_drug)
    subset = matrix[
        matrix["alteration_norm"].isin(variants)
        & matrix["drug_norm"].eq(drug)
        & ~matrix["is_resistance"].astype(str).str.lower().eq("true")
    ]
    if subset.empty:
        return "not_found", 0.0
    scored = subset.assign(level_score=subset["best_level"].map(LEVEL_SCORE).fillna(0.0))
    best = scored.sort_values("level_score", ascending=False).iloc[0]
    return str(best["best_level"]), float(best["level_score"])


def baseline1_score(
    baselines: dict[str, Any],
    oncokb: pd.DataFrame,
    bucket: str,
    class_name: str,
    drug: str,
) -> tuple[float, str]:
    b0, _ = baseline0_score(baselines, class_name, drug)
    level, level_score = oncokb_level_for_bucket(oncokb, bucket, drug)
    if level_score <= 0:
        return b0, level
    return max(b0, 0.30 * b0 + 0.70 * level_score), level


def rank_with_ties(scores: dict[str, float], target: str) -> tuple[int, bool, str]:
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    target_score = scores[target]
    rank = 1 + sum(score > target_score for _, score in ordered)
    top_score = ordered[0][1]
    top_ties = [drug for drug, score in ordered if score == top_score]
    return rank, target in top_ties, ";".join(top_ties)


def build_candidate_drugs(literature: pd.DataFrame, priors: dict[str, Any], baselines: dict[str, Any]) -> set[str]:
    candidate_drugs: set[str] = set()
    for class_map in priors.get("class_drug_priors", {}).values():
        candidate_drugs.update(norm_drug(drug) for drug in class_map)
    for drug in literature["drug"].dropna():
        candidate_drugs.update(expand_drug(drug))
    for recs in baselines["baseline_0"]["recommendations"].values():
        candidate_drugs.update(norm_drug(drug) for drug in recs if drug != "specialist_review")
    return {drug for drug in candidate_drugs if drug and drug != "specialist_review"}


def compound_group(bucket: object) -> str:
    lower = str(bucket).lower()
    variants = bucket_to_variants(bucket)
    if "compound" in lower or "_plus_" in lower or "classical_like_or_pacc" in lower or len(variants) > 1:
        return "compound_or_group"
    return "single_or_specific"


def class_group(class_name: str) -> str:
    return "PACC" if class_name == "pacc" else "non_PACC"


def compute_metric_summary(predictions: pd.DataFrame, folds: pd.DataFrame) -> pd.DataFrame:
    metric_rows: list[dict[str, Any]] = []
    metric_predictions = predictions[predictions["orr_evaluable"]].copy()
    metric_folds = folds[folds["is_metric_fold"]].copy()
    for method, score_col in [
        ("model", "model_score"),
        ("baseline_0", "baseline_0_score"),
        ("baseline_1", "baseline_1_score"),
    ]:
        auc = auc_score(
            metric_predictions["label_is_heldout_drug"].astype(int).tolist(),
            metric_predictions[score_col].astype(float).tolist(),
        )
        method_folds = metric_folds[metric_folds["method"].eq(method)]
        top1 = float(method_folds["top1_hit"].mean()) if not method_folds.empty else float("nan")
        mean_rank = float(method_folds["target_rank"].mean()) if not method_folds.empty else float("nan")
        metric_rows.append(
            {
                "method": method,
                "n_total_locked_rows": int(folds["fold_id"].nunique()),
                "n_orr_evaluable_folds": int(metric_predictions["fold_id"].nunique()),
                "candidate_drug_count": int(predictions["candidate_drug"].nunique()),
                "pooled_drug_auc": auc,
                "top1_hit_rate": top1,
                "mean_target_rank": mean_rank,
                "metric_definition": "Pooled one-positive-per-fold ROC AUC over candidate drugs; positive label is the held-out drug class/member.",
            }
        )

    metrics = pd.DataFrame(metric_rows)
    b0_auc = float(metrics.loc[metrics["method"].eq("baseline_0"), "pooled_drug_auc"].iloc[0])
    b1_auc = float(metrics.loc[metrics["method"].eq("baseline_1"), "pooled_drug_auc"].iloc[0])
    metrics["delta_auc_vs_baseline_0"] = metrics["pooled_drug_auc"] - b0_auc
    metrics["delta_auc_vs_baseline_1"] = metrics["pooled_drug_auc"] - b1_auc
    ordered = [
        "method",
        "n_total_locked_rows",
        "n_orr_evaluable_folds",
        "candidate_drug_count",
        "pooled_drug_auc",
        "delta_auc_vs_baseline_0",
        "delta_auc_vs_baseline_1",
        "top1_hit_rate",
        "mean_target_rank",
        "metric_definition",
    ]
    return metrics[ordered]


def run_backtest(
    literature: pd.DataFrame,
    priors: dict[str, Any],
    baselines: dict[str, Any],
    oncokb: pd.DataFrame,
    candidate_drugs: set[str],
    use_oncokb: bool = True,
    exact_floor: bool = True,
    config_fallback: bool = True,
    pacc_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prediction_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []

    for idx, heldout in literature.iterrows():
        train = literature.drop(index=idx)
        bucket = str(heldout["mutation_bucket"])
        class_name = primary_class_for_bucket(bucket)
        if pacc_only and class_name != "pacc":
            continue
        heldout_drugs = expand_drug(heldout["drug"])
        # A row with a drug class, such as second_generation_TKI, is counted
        # correct if any expanded member is top-ranked.
        target_drug = heldout_drugs[0]
        orr = pd.to_numeric(pd.Series([heldout["orr_percent"]]), errors="coerce").iloc[0]
        evaluable = bool(pd.notna(orr))

        score_maps = {"model": {}, "baseline_0": {}, "baseline_1": {}}
        details: dict[str, dict[str, str]] = {"model_source": {}, "baseline1_oncokb_level": {}}

        for drug in sorted(candidate_drugs):
            model_score, source = literature_model_score(
                train,
                priors,
                bucket,
                class_name,
                drug,
                exact_floor=exact_floor,
                config_fallback=config_fallback,
            )
            b0_score, b0_rec = baseline0_score(baselines, class_name, drug)
            b1_score, oncokb_level = (
                baseline1_score(baselines, oncokb, bucket, class_name, drug)
                if use_oncokb
                else (b0_score, "disabled")
            )

            score_maps["model"][drug] = model_score
            score_maps["baseline_0"][drug] = b0_score
            score_maps["baseline_1"][drug] = b1_score
            details["model_source"][drug] = source
            details["baseline1_oncokb_level"][drug] = oncokb_level

            prediction_rows.append(
                {
                    "fold_id": idx + 1,
                    "heldout_lock_id": heldout["lock_id"],
                    "heldout_source_study": heldout["source_study"],
                    "mutation_bucket": bucket,
                    "primary_class": class_name,
                    "class_group": class_group(class_name),
                    "compound_group": compound_group(bucket),
                    "heldout_drug_raw": heldout["drug"],
                    "heldout_target_drugs": ";".join(heldout_drugs),
                    "candidate_drug": drug,
                    "known_orr_percent": orr if evaluable else math.nan,
                    "orr_evaluable": evaluable,
                    "label_is_heldout_drug": int(drug in heldout_drugs),
                    "model_score": model_score,
                    "model_source": source,
                    "baseline_0_score": b0_score,
                    "baseline_0_recommended": b0_rec,
                    "baseline_1_score": b1_score,
                    "baseline_1_oncokb_level": oncokb_level,
                    "ablation_use_oncokb": use_oncokb,
                    "ablation_exact_floor": exact_floor,
                    "ablation_config_fallback": config_fallback,
                    "ablation_pacc_only": pacc_only,
                    "notes": heldout.get("notes", ""),
                }
            )

        for method, scores in score_maps.items():
            ranks = []
            top_hits = []
            top_ties = []
            for target in heldout_drugs:
                if target in scores:
                    rank, hit, ties = rank_with_ties(scores, target)
                    ranks.append(rank)
                    top_hits.append(hit)
                    top_ties.append(ties)
            fold_rows.append(
                {
                    "fold_id": idx + 1,
                    "heldout_lock_id": heldout["lock_id"],
                    "mutation_bucket": bucket,
                    "primary_class": class_name,
                    "class_group": class_group(class_name),
                    "compound_group": compound_group(bucket),
                    "heldout_drug_raw": heldout["drug"],
                    "heldout_target_drugs": ";".join(heldout_drugs),
                    "known_orr_percent": orr if evaluable else math.nan,
                    "orr_evaluable": evaluable,
                    "method": method,
                    "target_rank": min(ranks) if ranks else math.nan,
                    "top1_hit": bool(any(top_hits)),
                    "top1_drugs": top_ties[0] if top_ties else "",
                    "candidate_count": len(candidate_drugs),
                    "is_metric_fold": evaluable,
                    "ablation_use_oncokb": use_oncokb,
                    "ablation_exact_floor": exact_floor,
                    "ablation_config_fallback": config_fallback,
                    "ablation_pacc_only": pacc_only,
                }
            )

    predictions = pd.DataFrame(prediction_rows)
    folds = pd.DataFrame(fold_rows)
    metrics = compute_metric_summary(predictions, folds) if not predictions.empty else pd.DataFrame()
    return predictions, folds, metrics


def bootstrap_ci(
    predictions: pd.DataFrame,
    folds: pd.DataFrame,
    n_bootstrap: int = 1000,
    seed: int = 20260618,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = random.Random(seed)
    fold_ids = sorted(folds.loc[folds["is_metric_fold"], "fold_id"].unique().tolist())
    replicates: list[dict[str, Any]] = []
    for iteration in range(1, n_bootstrap + 1):
        sampled = [rng.choice(fold_ids) for _ in fold_ids]
        sampled_predictions = pd.concat(
            [predictions[predictions["fold_id"].eq(fold_id)] for fold_id in sampled],
            ignore_index=True,
        )
        sampled_folds = pd.concat(
            [folds[folds["fold_id"].eq(fold_id)] for fold_id in sampled],
            ignore_index=True,
        )
        metrics = compute_metric_summary(sampled_predictions, sampled_folds)
        for _, row in metrics.iterrows():
            replicates.append(
                {
                    "bootstrap_iteration": iteration,
                    "method": row["method"],
                    "pooled_drug_auc": row["pooled_drug_auc"],
                    "delta_auc_vs_baseline_0": row["delta_auc_vs_baseline_0"],
                    "delta_auc_vs_baseline_1": row["delta_auc_vs_baseline_1"],
                    "top1_hit_rate": row["top1_hit_rate"],
                    "mean_target_rank": row["mean_target_rank"],
                }
            )

    boot = pd.DataFrame(replicates)
    rows: list[dict[str, Any]] = []
    for method, group in boot.groupby("method"):
        for metric in ["pooled_drug_auc", "delta_auc_vs_baseline_0", "delta_auc_vs_baseline_1", "top1_hit_rate"]:
            rows.append(
                {
                    "method": method,
                    "metric": metric,
                    "bootstrap_n": n_bootstrap,
                    "mean": float(group[metric].mean()),
                    "ci95_lower": float(group[metric].quantile(0.025)),
                    "ci95_upper": float(group[metric].quantile(0.975)),
                }
            )
    return boot, pd.DataFrame(rows)


def ablation_summary(
    literature: pd.DataFrame,
    priors: dict[str, Any],
    baselines: dict[str, Any],
    oncokb: pd.DataFrame,
    candidate_drugs: set[str],
) -> pd.DataFrame:
    scenarios = [
        ("full_model", True, True, True, False),
        ("no_oncokb", False, True, True, False),
        ("no_exact_evidence_floor", True, False, True, False),
        ("no_config_fallback", True, True, False, False),
        ("pacc_only_subset", True, True, True, True),
    ]
    rows: list[pd.DataFrame] = []
    for name, use_oncokb, exact_floor, config_fallback, pacc_only in scenarios:
        _, _, metrics = run_backtest(
            literature,
            priors,
            baselines,
            oncokb,
            candidate_drugs,
            use_oncokb=use_oncokb,
            exact_floor=exact_floor,
            config_fallback=config_fallback,
            pacc_only=pacc_only,
        )
        metrics.insert(0, "scenario", name)
        rows.append(metrics)
    summary = pd.concat(rows, ignore_index=True)
    full = summary[summary["scenario"].eq("full_model")][
        ["method", "pooled_drug_auc", "top1_hit_rate", "mean_target_rank"]
    ].rename(
        columns={
            "pooled_drug_auc": "full_auc",
            "top1_hit_rate": "full_top1",
            "mean_target_rank": "full_mean_rank",
        }
    )
    summary = summary.merge(full, on="method", how="left")
    summary["delta_auc_vs_full_model"] = summary["pooled_drug_auc"] - summary["full_auc"]
    summary["delta_top1_vs_full_model"] = summary["top1_hit_rate"] - summary["full_top1"]
    summary["delta_mean_rank_vs_full_model"] = summary["mean_target_rank"] - summary["full_mean_rank"]
    return summary.drop(columns=["full_auc", "full_top1", "full_mean_rank"])


def subgroup_rank_summary(folds: pd.DataFrame) -> pd.DataFrame:
    metric = folds[folds["is_metric_fold"]].copy()
    rows: list[dict[str, Any]] = []
    for group_col in ["class_group", "compound_group"]:
        for group_value, group in metric.groupby(group_col):
            pivot_rank = group.pivot_table(index="fold_id", columns="method", values="target_rank", aggfunc="min")
            pivot_hit = group.pivot_table(index="fold_id", columns="method", values="top1_hit", aggfunc="max")
            for method, method_group in group.groupby("method"):
                rows.append(
                    {
                        "subgroup_type": group_col,
                        "subgroup": group_value,
                        "method": method,
                        "n_folds": int(method_group["fold_id"].nunique()),
                        "mean_target_rank": float(method_group["target_rank"].mean()),
                        "top1_hit_rate": float(method_group["top1_hit"].mean()),
                        "delta_mean_rank_vs_baseline_0": (
                            float(pivot_rank[method].mean() - pivot_rank["baseline_0"].mean())
                            if method in pivot_rank and "baseline_0" in pivot_rank
                            else float("nan")
                        ),
                        "delta_top1_vs_baseline_0": (
                            float(pivot_hit[method].mean() - pivot_hit["baseline_0"].mean())
                            if method in pivot_hit and "baseline_0" in pivot_hit
                            else float("nan")
                        ),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Leave-one-out back-testing for locked EGFR-PACC ORR priors.")
    parser.add_argument("--literature", default=str(ROOT / "data/processed/literature_orr_locked.csv"))
    parser.add_argument("--out-dir", default=str(ROOT / "reports/backtest_loo"))
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260618)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    literature = pd.read_csv(args.literature)
    priors = load_yaml(ROOT / "configs/priors.yaml")
    baselines = load_yaml(ROOT / "configs/baselines.yaml")
    oncokb_path = ROOT / "data/processed/reference/oncokb_egfr_drug_level_matrix.csv"
    oncokb = pd.read_csv(oncokb_path) if oncokb_path.exists() else pd.DataFrame()
    candidate_drugs = build_candidate_drugs(literature, priors, baselines)

    predictions, folds, metrics = run_backtest(
        literature,
        priors,
        baselines,
        oncokb,
        candidate_drugs,
        use_oncokb=True,
        exact_floor=True,
        config_fallback=True,
        pacc_only=False,
    )

    predictions.to_csv(out_dir / "loo_predictions.csv", index=False)
    folds.to_csv(out_dir / "loo_fold_summary.csv", index=False)
    metrics.to_csv(out_dir / "loo_metric_summary.csv", index=False)
    boot_reps, boot_ci = bootstrap_ci(predictions, folds, n_bootstrap=args.bootstrap, seed=args.seed)
    boot_reps.to_csv(out_dir / "loo_bootstrap_replicates.csv", index=False)
    boot_ci.to_csv(out_dir / "loo_bootstrap_ci.csv", index=False)
    ablations = ablation_summary(literature, priors, baselines, oncokb, candidate_drugs)
    ablations.to_csv(out_dir / "loo_ablation_summary.csv", index=False)
    subgroups = subgroup_rank_summary(folds)
    subgroups.to_csv(out_dir / "loo_subgroup_rank_summary.csv", index=False)

    report = [
        "# Leave-One-Out Back-Testing Report",
        "",
        f"- Total locked literature rows: {len(literature)}",
        f"- ORR-evaluable folds used for metrics: {int(predictions[predictions['orr_evaluable']]['fold_id'].nunique())}",
        f"- Candidate drugs: {len(candidate_drugs)}",
        f"- Bootstrap replicates: {args.bootstrap}",
        "- Outcome label: held-out mutation-drug ORR row; one positive candidate per fold.",
        "- Non-ORR rows are retained in fold audit but excluded from AUC/top-1 metrics.",
        "",
        "## Metric Summary",
        "",
        metrics.to_markdown(index=False),
        "",
        "## Bootstrap 95% CI",
        "",
        boot_ci.to_markdown(index=False),
        "",
        "## Ablation Summary",
        "",
        ablations.to_markdown(index=False),
        "",
        "## Subgroup Rank Summary",
        "",
        subgroups.to_markdown(index=False),
        "",
        "## Outputs",
        "",
        "- `loo_predictions.csv`",
        "- `loo_fold_summary.csv`",
        "- `loo_metric_summary.csv`",
        "- `loo_bootstrap_replicates.csv`",
        "- `loo_bootstrap_ci.csv`",
        "- `loo_ablation_summary.csv`",
        "- `loo_subgroup_rank_summary.csv`",
    ]
    (out_dir / "loo_backtest_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(
        {
            "predictions": str((out_dir / "loo_predictions.csv").resolve()),
            "fold_summary": str((out_dir / "loo_fold_summary.csv").resolve()),
            "metric_summary": str((out_dir / "loo_metric_summary.csv").resolve()),
            "report": str((out_dir / "loo_backtest_report.md").resolve()),
        }
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
