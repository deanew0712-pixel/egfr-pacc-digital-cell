from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "csbj_submission_figures_2026-06-21"
PUB = ROOT / "reports" / "figures_publication"

COLORS = {
    "model": "#1F77B4",
    "baseline0": "#777777",
    "baseline1": "#2CA25F",
    "pacc": "#1F77B4",
    "classical": "#D95F02",
    "compound": "#756BB1",
    "single": "#31A354",
    "high": "#1F77B4",
    "low": "#D95F02",
    "accent": "#D62728",
    "grid": "#D9D9D9",
    "text": "#222222",
}

METHOD_LABEL = {
    "baseline_0": "Baseline-0",
    "baseline_1": "Baseline-1",
    "model": "EPDRI model",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "figure.titlesize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "lines.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_all(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png", "tiff"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=600, bbox_inches="tight")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def clean_axis(ax: plt.Axes, axis: str = "y") -> None:
    ax.grid(axis=axis, color=COLORS["grid"], linewidth=0.5, alpha=0.8)
    ax.set_axisbelow(True)


def rounded_box(
    ax: plt.Axes,
    xy,
    wh,
    text: str,
    fc: str,
    ec: str = "#444444",
    fontsize: float = 7.0,
) -> None:
    x, y = xy
    w, h = wh
    patch = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.025,rounding_size=0.035",
        linewidth=0.9,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=COLORS["text"], fontsize=fontsize)


def arrow(ax: plt.Axes, p1, p2) -> None:
    ax.annotate(
        "",
        xy=p2,
        xytext=p1,
        arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#333333", shrinkA=2, shrinkB=2),
    )


def figure1_framework() -> None:
    fig = plt.figure(figsize=(7.1, 4.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 1.0], hspace=0.28)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax.set_xlim(0, 12)
    ax.set_ylim(2.05, 6)
    ax.axis("off")

    rounded_box(ax, (0.35, 4.25), (2.05, 0.9), "EGFR variants\ncompound status", "#EAF3FB")
    rounded_box(ax, (0.35, 2.55), (2.05, 0.9), "Clinical context\nCNS / co-alterations", "#F3F0FA")
    rounded_box(ax, (3.15, 4.25), (2.05, 0.9), "Structural class\nPACC and comparators", "#E5F5E0")
    rounded_box(ax, (3.15, 2.55), (2.05, 0.9), "Compound resolver\nuncertainty flag", "#FFF2CC")
    rounded_box(ax, (5.95, 4.25), (2.15, 0.9), "Locked priors\n34 rows; 29 ORR folds", "#FEE5D9")
    rounded_box(ax, (5.95, 2.55), (2.15, 0.9), "Bayesian pooling\nclass-level hierarchy", "#FDE0DD")
    rounded_box(ax, (8.9, 3.25), (2.55, 1.15), "EPDRI output\nranked inhibitors\nscore, grade, uncertainty", "#F0F0F0")

    arrow(ax, (2.3, 4.8), (3.1, 4.8))
    arrow(ax, (2.3, 3.1), (3.1, 3.1))
    arrow(ax, (5.2, 4.8), (6.0, 4.8))
    arrow(ax, (5.2, 3.1), (6.0, 3.1))
    arrow(ax, (8.2, 4.8), (9.0, 4.25))
    arrow(ax, (8.2, 3.1), (9.0, 3.75))
    arrow(ax, (1.3, 4.35), (1.3, 3.55))
    arrow(ax, (4.15, 4.35), (4.15, 3.55))
    arrow(ax, (7.1, 4.35), (7.1, 3.55))
    ax.text(0.0, 5.75, "A", fontsize=10, fontweight="bold", va="top")
    ax.text(0.35, 5.75, "Digital-cell scoring workflow", fontsize=9.5, fontweight="bold", va="top")

    ax2.axis("off")
    ax2.set_xlim(0, 12)
    ax2.set_ylim(0, 2)
    ax2.text(0, 1.8, "B", fontsize=10, fontweight="bold", va="top")
    ax2.text(0.35, 1.8, "Model evaluation design", fontsize=9.5, fontweight="bold", va="top")
    rounded_box(ax2, (0.45, 0.55), (2.4, 0.75), "Literature lock\nbefore validation", "#FEE5D9")
    rounded_box(ax2, (3.35, 0.55), (2.4, 0.75), "Leave-one-out\nback-testing", "#EAF3FB")
    rounded_box(ax2, (6.25, 0.55), (2.4, 0.75), "Independent local\n31-case validation", "#E5F5E0")
    rounded_box(ax2, (9.15, 0.55), (2.4, 0.75), "Structural sensitivity\nVina / OpenMM", "#F3F0FA")
    for x in (2.85, 5.75, 8.65):
        arrow(ax2, (x, 0.92), (x + 0.5, 0.92))

    save_all(fig, "Figure_1_framework")
    plt.close(fig)


def figure2_backtesting() -> None:
    metrics = pd.read_csv(ROOT / "reports/backtest_loo/loo_metric_summary.csv")
    ci = pd.read_csv(ROOT / "reports/backtest_loo/loo_bootstrap_ci.csv")
    ablation = pd.read_csv(ROOT / "reports/backtest_loo/loo_ablation_summary.csv")
    folds = pd.read_csv(ROOT / "reports/backtest_loo/loo_subgroup_rank_summary.csv")

    fig = plt.figure(figsize=(7.1, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    order = ["baseline_0", "baseline_1", "model"]
    labels = [METHOD_LABEL[x] for x in order]
    colors = [COLORS["baseline0"], COLORS["baseline1"], COLORS["model"]]
    vals = [float(metrics.loc[metrics.method.eq(x), "pooled_drug_auc"].iloc[0]) for x in order]
    lowers = [
        float(ci[(ci.method.eq(x)) & (ci.metric.eq("pooled_drug_auc"))]["ci95_lower"].iloc[0])
        for x in order
    ]
    uppers = [
        float(ci[(ci.method.eq(x)) & (ci.metric.eq("pooled_drug_auc"))]["ci95_upper"].iloc[0])
        for x in order
    ]
    ax1.bar(labels, vals, color=colors, edgecolor="#333333", linewidth=0.5)
    ax1.errorbar(labels, vals, yerr=[np.array(vals) - np.array(lowers), np.array(uppers) - np.array(vals)], fmt="none", ecolor="#222222", capsize=3, linewidth=0.8)
    ax1.set_ylim(0.65, 1.0)
    ax1.set_ylabel("Pooled drug-level AUC")
    ax1.tick_params(axis="x", rotation=25)
    clean_axis(ax1)
    panel_label(ax1, "A")

    delta_specs = [
        ("Model - B0", "model", "delta_auc_vs_baseline_0", COLORS["model"]),
        ("Model - B1", "model", "delta_auc_vs_baseline_1", COLORS["accent"]),
        ("B1 - B0", "baseline_1", "delta_auc_vs_baseline_0", COLORS["baseline1"]),
    ]
    y = np.arange(len(delta_specs))
    means, lo, hi = [], [], []
    for _, method, metric, _ in delta_specs:
        row = ci[(ci.method.eq(method)) & (ci.metric.eq(metric))].iloc[0]
        means.append(float(row["mean"]))
        lo.append(float(row["ci95_lower"]))
        hi.append(float(row["ci95_upper"]))
    ax2.axvline(0, color="#333333", linewidth=0.8)
    for yy, (label, _, _, color), mean, low, high in zip(y, delta_specs, means, lo, hi):
        ax2.errorbar(mean, yy, xerr=[[mean - low], [high - mean]], fmt="o", color=color, ecolor="#222222", capsize=3, linewidth=0.8)
    ax2.set_yticks(y, [x[0] for x in delta_specs])
    ax2.set_xlabel("Delta AUC (bootstrap 95% CI)")
    clean_axis(ax2, "x")
    panel_label(ax2, "B")

    abl = ablation[ablation.method.eq("model")].copy()
    scenario_label = {
        "full_model": "Full model",
        "no_oncokb": "No OncoKB",
        "no_exact_evidence_floor": "No exact-evidence floor",
        "no_config_fallback": "No config fallback",
        "pacc_only_subset": "PACC-only subset",
    }
    abl["label"] = abl.scenario.map(scenario_label)
    abl = abl.set_index("scenario").loc[
        ["full_model", "no_exact_evidence_floor", "no_config_fallback", "pacc_only_subset", "no_oncokb"]
    ].reset_index()
    ax3.barh(abl["label"], abl["pooled_drug_auc"], color="#9ECAE1", edgecolor="#333333", linewidth=0.5)
    ax3.axvline(float(metrics.loc[metrics.method.eq("baseline_0"), "pooled_drug_auc"].iloc[0]), color=COLORS["baseline0"], linestyle="--", linewidth=0.9)
    ax3.axvline(float(metrics.loc[metrics.method.eq("baseline_1"), "pooled_drug_auc"].iloc[0]), color=COLORS["baseline1"], linestyle="--", linewidth=0.9)
    ax3.set_xlim(0.76, 0.93)
    ax3.set_xlabel("Model AUC")
    clean_axis(ax3, "x")
    panel_label(ax3, "C")

    rank = folds[folds.subgroup_type.eq("class_group")]
    pivot = rank.pivot(index="subgroup", columns="method", values="mean_target_rank").loc[["PACC", "non_PACC"]]
    x = np.arange(len(pivot.index))
    width = 0.25
    for off, method, color in [(-width, "baseline_0", COLORS["baseline0"]), (0, "baseline_1", COLORS["baseline1"]), (width, "model", COLORS["model"])]:
        ax4.bar(x + off, pivot[method].values, width=width, color=color, edgecolor="#333333", linewidth=0.5, label=METHOD_LABEL[method])
    ax4.set_xticks(x, ["PACC\n(n=25 folds)", "non-PACC\n(n=4 folds)"])
    ax4.set_ylabel("Mean target rank")
    ax4.legend(frameon=False, ncol=1, loc="upper left")
    clean_axis(ax4)
    panel_label(ax4, "D")

    fig.suptitle("Literature leave-one-out back-testing", x=0.02, ha="left", fontweight="bold")
    save_all(fig, "Figure_2_backtesting")
    plt.close(fig)


def _metric(metrics: pd.DataFrame, analysis_set: str, metric: str) -> float:
    return float(metrics[(metrics.analysis_set.eq(analysis_set)) & (metrics.metric.eq(metric))]["estimate"].iloc[0])


def _metric_ci(metrics: pd.DataFrame, analysis_set: str, metric: str) -> tuple[float, float, float]:
    row = metrics[(metrics.analysis_set.eq(analysis_set)) & (metrics.metric.eq(metric))].iloc[0]
    return float(row["estimate"]), float(row["ci95_lower"]), float(row["ci95_upper"])


def figure3_clinical_validation() -> None:
    pred = pd.read_csv(ROOT / "reports/clinical_validation/clinical_patient_predictions_31.csv")
    metrics = pd.read_csv(ROOT / "reports/clinical_validation/clinical_validation_metrics.csv")
    subgroup = pd.read_csv(ROOT / "reports/clinical_validation/clinical_validation_subgroup_metrics.csv")

    native = pred[pred.analysis_set_model_native_monotherapy.astype(bool)].copy()
    proxy = pred[pred.analysis_set_proxy_sensitivity.astype(bool)].copy()

    fig = plt.figure(figsize=(7.1, 6.5))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.38)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    class_colors = native.primary_class.map({"pacc": COLORS["pacc"], "classical_like": COLORS["classical"]}).fillna("#999999")
    ax1.scatter(native.actual_model_score, native.pfs1_months, s=34, c=class_colors, edgecolor="#333333", linewidth=0.4)
    m, b = np.polyfit(native.actual_model_score, native.pfs1_months, 1)
    xx = np.linspace(native.actual_model_score.min() - 0.02, native.actual_model_score.max() + 0.02, 100)
    ax1.plot(xx, m * xx + b, color="#333333", linewidth=1.0)
    ax1.set_xlabel("Administered-drug EPDRI")
    ax1.set_ylabel("PFS1 (months)")
    ax1.text(0.03, 0.95, "n=24\nrho=0.288\nC-index=0.595", transform=ax1.transAxes, va="top")
    clean_axis(ax1)
    panel_label(ax1, "A")

    sets = [
        ("Model-native\n(n=24)", "model_native_monotherapy"),
        ("Proxy sensitivity\n(n=29)", "proxy_sensitivity_including_icotinib"),
    ]
    x = np.arange(len(sets))
    for i, metric in enumerate(["spearman_score_vs_pfs1", "pairwise_c_index"]):
        vals, low, high = [], [], []
        for _, aset in sets:
            v, l, u = _metric_ci(metrics, aset, metric)
            vals.append(v)
            low.append(l)
            high.append(u)
        offset = -0.12 if i == 0 else 0.12
        label = "Spearman rho" if i == 0 else "Pairwise C-index"
        color = COLORS["pacc"] if i == 0 else COLORS["compound"]
        ax2.errorbar(x + offset, vals, yerr=[np.array(vals) - np.array(low), np.array(high) - np.array(vals)], fmt="o", color=color, ecolor="#222222", capsize=3, label=label)
    ax2.axhline(0, color="#888888", linewidth=0.7)
    ax2.set_xticks(x, [s[0] for s in sets])
    ax2.set_ylim(-0.35, 0.9)
    ax2.set_ylabel("Estimate with bootstrap 95% CI")
    ax2.legend(frameon=False, loc="upper right")
    clean_axis(ax2)
    panel_label(ax2, "B")

    top_metrics = ["top1_coverage", "top2_coverage", "top3_coverage"]
    width = 0.32
    for j, (set_label, aset) in enumerate(sets):
        vals = [_metric(metrics, aset, m) * 100 for m in top_metrics]
        ax3.bar(np.arange(3) + (j - 0.5) * width, vals, width=width, label=set_label.replace("\n", " "), edgecolor="#333333", linewidth=0.5)
    ax3.set_xticks(np.arange(3), ["Top 1", "Top 2", "Top 3"])
    ax3.set_ylim(0, 70)
    ax3.set_ylabel("Coverage of administered drug (%)")
    ax3.legend(frameon=False, loc="upper right")
    clean_axis(ax3)
    panel_label(ax3, "C")

    sub = subgroup[subgroup.analysis_set.eq("proxy_sensitivity_including_icotinib")].copy()
    sub = sub[sub.subgroup_family.isin(["PACC vs non-PACC", "single vs compound mutation"])]
    sub["label"] = sub["subgroup"].replace({"pacc": "PACC", "classical_like": "classical-like"})
    sub = sub.iloc[::-1]
    ax4.axvline(0.5, color="#888888", linewidth=0.7, linestyle="--")
    y = np.arange(len(sub))
    vals = sub.pairwise_c_index.astype(float).values
    low = sub.cindex_ci95_lower.astype(float).values
    high = sub.cindex_ci95_upper.astype(float).values
    ax4.errorbar(vals, y, xerr=[vals - low, high - vals], fmt="o", color=COLORS["pacc"], ecolor="#222222", capsize=3)
    wrapped = {
        "classical-like": "classical-like\n(n=8)",
        "PACC": "PACC\n(n=21)",
        "Compound mutation": "Compound\nmutation\n(n=13)",
        "Single mutation": "Single\nmutation\n(n=16)",
    }
    ax4.set_yticks(y, [wrapped.get(r.label, f"{r.label}\n(n={int(r.n)})") for r in sub.itertuples()])
    ax4.set_xlabel("Subgroup C-index (proxy set)")
    ax4.set_xlim(0.1, 0.9)
    clean_axis(ax4, "x")
    panel_label(ax4, "D")

    fig.suptitle("Independent local clinical validation", x=0.02, ha="left", fontweight="bold")
    save_all(fig, "Figure_3_clinical_validation")
    plt.close(fig)


def figure4_km() -> None:
    km = pd.read_csv(ROOT / "reports/clinical_validation/clinical_validation_km_curves.csv")
    summary = pd.read_csv(ROOT / "reports/clinical_validation/clinical_validation_km_logrank_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.3), sharey=True)
    specs = [
        ("model_native_monotherapy", "Model-native monotherapy"),
        ("proxy_sensitivity_including_icotinib", "Proxy sensitivity"),
    ]
    for ax, (aset, title), label in zip(axes, specs, ["A", "B"]):
        for group, color in [("High EPDRI", COLORS["high"]), ("Low EPDRI", COLORS["low"])]:
            d = km[(km.analysis_set.eq(aset)) & (km.score_group.eq(group))].copy()
            ax.step(d.time.astype(float), d.survival.astype(float), where="post", label=group, color=color)
        row = summary[summary.analysis_set.eq(aset)].iloc[0]
        ax.set_title(title)
        ax.set_xlabel("PFS1 time (months)")
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, max(26, km[km.analysis_set.eq(aset)].time.astype(float).max() + 1))
        ax.text(
            0.04,
            0.08,
            f"High vs low median PFS1:\n{float(row.median_pfs_high):.2f} vs {float(row.median_pfs_low):.2f} mo\nlog-rank p={float(row.logrank_p_value):.3f}",
            transform=ax.transAxes,
            fontsize=7.5,
        )
        clean_axis(ax)
        panel_label(ax, label)
    axes[0].set_ylabel("Progression-free survival")
    axes[1].legend(frameon=False, loc="upper right")
    fig.suptitle("Kaplan-Meier analysis by EPDRI score stratum", x=0.02, ha="left", fontweight="bold")
    save_all(fig, "Figure_4_km_logrank")
    plt.close(fig)


def figure_s1_demo() -> None:
    ranking = pd.read_csv(ROOT / "reports/demo_case_ranking_pymc_formal_4chains_2000draws.csv")
    ranking = ranking.sort_values("epdri_score", ascending=True)
    fig, ax = plt.subplots(figsize=(7.1, 3.8))
    grade_colors = {"A": "#1A9850", "B": "#66BD63", "C": "#FDAE61", "D": "#D73027"}
    colors = ranking.recommendation_grade.map(grade_colors).fillna("#999999")
    ax.barh(ranking.drug.str.capitalize(), ranking.epdri_score, color=colors, edgecolor="#333333", linewidth=0.5)
    for y, row in enumerate(ranking.itertuples()):
        ax.text(float(row.epdri_score) + 0.012, y, f"{float(row.epdri_score):.3f} | grade {row.recommendation_grade}", va="center", fontsize=7)
    ax.set_xlim(0, 0.66)
    ax.set_xlabel("Collapsed EPDRI score")
    ax.set_title("Demo case: EGFR G719S + L861Q with CNS and co-alteration modifiers", loc="left", fontweight="bold")
    clean_axis(ax, "x")
    save_all(fig, "Figure_S1_demo_case_ranking")
    plt.close(fig)


def figure_s2_bootstrap() -> None:
    boot = pd.read_csv(ROOT / "reports/clinical_validation/clinical_validation_bootstrap_distributions.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.1), sharey=False)
    specs = [
        ("model_native_monotherapy", "spearman_score_vs_pfs1", "Model-native Spearman"),
        ("model_native_monotherapy", "pairwise_c_index", "Model-native C-index"),
        ("proxy_sensitivity_including_icotinib", "spearman_score_vs_pfs1", "Proxy Spearman"),
        ("proxy_sensitivity_including_icotinib", "pairwise_c_index", "Proxy C-index"),
    ]
    for ax, (aset, col, title), label in zip(axes.ravel(), specs, ["A", "B", "C", "D"]):
        vals = boot[boot.analysis_set.eq(aset)][col].astype(float).dropna()
        ax.hist(vals, bins=30, color="#9ECAE1", edgecolor="#333333", linewidth=0.3)
        ax.axvline(vals.quantile(0.025), color="#444444", linestyle="--", linewidth=0.8)
        ax.axvline(vals.quantile(0.975), color="#444444", linestyle="--", linewidth=0.8)
        ax.set_title(title)
        ax.set_xlabel("Bootstrap estimate")
        ax.set_ylabel("Frequency")
        clean_axis(ax)
        panel_label(ax, label)
    fig.suptitle("Bootstrap distributions for local validation metrics", x=0.02, ha="left", fontweight="bold")
    save_all(fig, "Figure_S2_bootstrap_distributions")
    plt.close(fig)


def figure_s3_docking() -> None:
    scores = pd.read_csv(ROOT / "reports/docking_pacc_core/docking_scores.csv")
    variants = ["WT", "G719S", "L861Q", "S768I"]
    drugs = ["afatinib", "osimertinib", "furmonertinib"]
    heat = scores.pivot(index="variant_model", columns="drug", values="vina_score_kcal_mol").loc[variants, drugs]
    ranks = scores.pivot(index="variant_model", columns="drug", values="rank_within_variant").loc[variants, drugs]

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.8), gridspec_kw={"width_ratios": [1.05, 1.1]})
    ax = axes[0]
    im = ax.imshow(heat.values, cmap="viridis_r", aspect="auto")
    ax.set_xticks(np.arange(len(drugs)), ["Afatinib", "Osimertinib", "Furmonertinib"], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(variants)), variants)
    for i in range(len(variants)):
        for j in range(len(drugs)):
            ax.text(j, i, f"{heat.iloc[i, j]:.2f}\nrank {int(ranks.iloc[i, j])}", ha="center", va="center", fontsize=7, color="white" if heat.iloc[i, j] < -9.4 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Vina score (kcal/mol)")
    ax.set_title("Basic non-covalent docking")
    panel_label(ax, "A")

    ax2 = axes[1]
    ax2.axis("off")
    ax2.set_title("Representative rank-1 pairs", loc="left")
    top = scores.sort_values(["variant_model", "rank_within_variant"]).groupby("variant_model").head(1)
    top = top.set_index("variant_model").loc[variants].reset_index()
    y0 = 0.82
    for idx, row in enumerate(top.itertuples()):
        y = y0 - idx * 0.18
        rounded_box(ax2, (0.05, y - 0.055), (0.9, 0.1), f"{row.variant_model}: {row.drug} ({row.vina_score_kcal_mol:.2f} kcal/mol)", "#F0F0F0")
    ax2.text(
        0.05,
        0.08,
        "Scores are rank-order Vina outputs from homology-derived mutant receptor models; they are not binding free energies.",
        transform=ax2.transAxes,
        fontsize=7.5,
        color="#444444",
        va="bottom",
    )
    panel_label(ax2, "B")
    save_all(fig, "Figure_S3_basic_docking")
    plt.close(fig)


def figure_s4_workflow_covalent() -> None:
    workflow = pd.read_csv(ROOT / "reports/docking_workflow_sensitivity/receptor_workflow_sensitivity_summary.csv")
    cov = pd.read_csv(ROOT / "reports/docking_workflow_sensitivity/covalent_warhead_geometry_best_poses.csv")
    variants = ["WT", "G719S", "L861Q", "S768I"]
    drugs = ["afatinib", "osimertinib", "furmonertinib"]

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.8))
    top = workflow.pivot(index="variant_model", columns="drug", values="best_rank_count").loc[variants, drugs]
    im = axes[0].imshow(top.values, cmap="Blues", vmin=0, vmax=3, aspect="auto")
    axes[0].set_xticks(np.arange(len(drugs)), ["Afatinib", "Osimertinib", "Furmonertinib"], rotation=35, ha="right")
    axes[0].set_yticks(np.arange(len(variants)), variants)
    for i in range(len(variants)):
        for j in range(len(drugs)):
            axes[0].text(j, i, f"{int(top.iloc[i, j])}/3", ha="center", va="center", fontsize=8)
    cb = fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)
    cb.set_label("Rank-1 workflows")
    axes[0].set_title("Receptor-workflow rank stability")
    panel_label(axes[0], "A")

    cov["near"] = cov.near_attack_pose.astype(str).str.lower().eq("true")
    for drug, color in zip(drugs, ["#D95F02", "#1F77B4", "#2CA25F"]):
        d = cov[cov.drug.eq(drug)]
        axes[1].scatter(
            d.cys_sg_to_acrylamide_beta_c_distance_A.astype(float),
            d.sg_beta_alpha_angle_deg.astype(float),
            s=45,
            label=drug.capitalize(),
            edgecolor="#333333",
            linewidth=0.4,
            color=color,
            alpha=0.9,
        )
        for r in d.itertuples():
            axes[1].text(float(r.cys_sg_to_acrylamide_beta_c_distance_A) + 0.02, float(r.sg_beta_alpha_angle_deg) + 0.8, r.variant_model, fontsize=6.5)
    axes[1].axvspan(0, 5.0, color="#E5F5E0", alpha=0.35)
    axes[1].axhspan(70, 120, color="#E5F5E0", alpha=0.35)
    axes[1].set_xlabel("Cys SG to warhead beta-C (Angstrom)")
    axes[1].set_ylabel("SG-beta-alpha angle (degrees)")
    axes[1].set_title("Pre-covalent geometry screen")
    axes[1].legend(frameon=False, loc="lower right")
    clean_axis(axes[1])
    panel_label(axes[1], "B")
    save_all(fig, "Figure_S4_workflow_covalent")
    plt.close(fig)


def figure_s5_openmm() -> None:
    comp = pd.read_csv(ROOT / "reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/openmm_gpu_relaxed_vs_unrelaxed_comparison.csv")
    rmsd = pd.read_csv(ROOT / "reports/docking_md_relaxed_gpu/colab_output/openmm_gpu_relaxed/openmm_gpu_relaxed_receptor_metrics_corrected_rmsd.csv")
    variants = ["WT", "G719S", "L861Q", "S768I"]
    drugs = ["afatinib", "osimertinib", "furmonertinib"]

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.6))
    for drug, color in zip(drugs, ["#D95F02", "#1F77B4", "#2CA25F"]):
        d = comp[comp.drug.eq(drug)]
        axes[0].scatter(d.unrelaxed_score_kcal_mol.astype(float), d.relaxed_score_kcal_mol.astype(float), s=46, color=color, edgecolor="#333333", linewidth=0.4, label=drug.capitalize())
        for r in d.itertuples():
            axes[0].text(float(r.unrelaxed_score_kcal_mol) + 0.015, float(r.relaxed_score_kcal_mol) + 0.015, r.variant_model, fontsize=6.5)
    mn = min(comp.unrelaxed_score_kcal_mol.min(), comp.relaxed_score_kcal_mol.min()) - 0.05
    mx = max(comp.unrelaxed_score_kcal_mol.max(), comp.relaxed_score_kcal_mol.max()) + 0.05
    axes[0].plot([mn, mx], [mn, mx], color="#777777", linewidth=0.9)
    axes[0].set_xlabel("Unrelaxed Vina score (kcal/mol)")
    axes[0].set_ylabel("OpenMM-relaxed Vina score (kcal/mol)")
    axes[0].legend(frameon=False, loc="lower right")
    clean_axis(axes[0])
    panel_label(axes[0], "A")

    delta = comp.pivot(index="variant_model", columns="drug", values="score_delta_relaxed_minus_unrelaxed").loc[variants, drugs]
    im = axes[1].imshow(delta.values, cmap="coolwarm", vmin=0.2, vmax=1.2, aspect="auto")
    axes[1].set_xticks(np.arange(len(drugs)), ["Afatinib", "Osimertinib", "Furmonertinib"], rotation=35, ha="right")
    axes[1].set_yticks(np.arange(len(variants)), variants)
    for i, v in enumerate(variants):
        for j, d in enumerate(drugs):
            mark = "*" if bool(comp[(comp.variant_model.eq(v)) & (comp.drug.eq(d))]["rank_changed"].iloc[0]) else ""
            axes[1].text(j, i, f"{delta.loc[v, d]:.2f}{mark}", ha="center", va="center", fontsize=7.5)
    cb = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    cb.set_label("Relaxed - unrelaxed score")
    note = "; ".join([f"{r.variant_model} RMSD {float(r.corrected_heavy_atom_rmsd_A):.2f} A" for r in rmsd.itertuples()])
    axes[1].set_title("Score shift after receptor relaxation")
    axes[1].text(0, -0.32, f"* rank changed. Corrected heavy-atom RMSD: {note}.", transform=axes[1].transAxes, fontsize=6.8, color="#444444")
    panel_label(axes[1], "B")
    save_all(fig, "Figure_S5_openmm_relaxed")
    plt.close(fig)


def write_manifest() -> None:
    rows = []
    for p in sorted(OUT.glob("*")):
        rows.append({"file": p.name, "bytes": p.stat().st_size})
    pd.DataFrame(rows).to_csv(OUT / "figure_manifest.csv", index=False)
    qc = OUT / "CSBJ_figure_QC.md"
    qc.write_text(
        "\n".join(
            [
                "# CSBJ figure QC",
                "",
                "Generated on 2026-06-21.",
                "",
                "Design target: two-column figures at approximately 180 mm width, minimum 600 dpi raster export, editable PDF/SVG vector export, panel letters A-D, >=7.5 pt visible font in plotted panels.",
                "",
                "CSBJ accepts PDF/SVG-derived vector artwork and 300 dpi or higher bitmap images; these outputs include PDF, SVG, PNG, and TIFF copies for each generated figure.",
                "",
                "Main figures:",
                "- Figure_1_framework",
                "- Figure_2_backtesting",
                "- Figure_3_clinical_validation",
                "- Figure_4_km_logrank",
                "",
                "Supplementary figures:",
                "- Figure_S1_demo_case_ranking",
                "- Figure_S2_bootstrap_distributions",
                "- Figure_S3_basic_docking",
                "- Figure_S4_workflow_covalent",
                "- Figure_S5_openmm_relaxed",
                "",
                "Interpretive cautions retained in figure design:",
                "- Subgroup plots label sample sizes and should be interpreted as hypothesis-generating.",
                "- Docking panels label Vina outputs as rank-order scores, not binding free energies.",
                "- OpenMM-relaxed panels use corrected matched-heavy-atom RMSD, not the invalid original atom-order RMSD.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    figure1_framework()
    figure2_backtesting()
    figure3_clinical_validation()
    figure4_km()
    figure_s1_demo()
    figure_s2_bootstrap()
    figure_s3_docking()
    figure_s4_workflow_covalent()
    figure_s5_openmm()
    write_manifest()
    print(f"CSBJ figures written to {OUT}")


if __name__ == "__main__":
    main()
