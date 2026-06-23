from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "reports" / "figures_publication"


COLORS = {
    "model": "#2C7FB8",
    "baseline_0": "#7F8C8D",
    "baseline_1": "#31A354",
    "accent": "#E6550D",
    "light_blue": "#EAF3FB",
    "light_green": "#EDF8E9",
    "light_orange": "#FFF4E6",
    "light_purple": "#F3E8FF",
    "light_gray": "#F4F4F4",
    "text": "#222222",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=300, bbox_inches="tight")


def draw_box(ax, xy, wh, text, facecolor, edgecolor="#333333", fontsize=8):
    x, y = xy
    w, h = wh
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.0,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color=COLORS["text"])


def arrow(ax, start, end):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="->", lw=1.1, color="#333333", shrinkA=3, shrinkB=3),
    )


def figure1_framework() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    draw_box(
        ax,
        (0.35, 4.15),
        (2.1, 0.9),
        "Locked evidence\n34 literature rows\nOpen reference priors",
        COLORS["light_blue"],
    )
    draw_box(
        ax,
        (0.35, 2.25),
        (2.1, 0.9),
        "Mutation layer\nRobichaux classes\ncompound resolver",
        COLORS["light_green"],
    )
    draw_box(
        ax,
        (3.1, 4.15),
        (2.1, 0.9),
        "Bayesian layer\nclass-level partial pooling\nlocked priors only",
        COLORS["light_orange"],
    )
    draw_box(
        ax,
        (3.1, 2.25),
        (2.1, 0.9),
        "Clinical modifiers\nCNS disease\nbypass/co-mutation",
        "#E8F4F2",
    )
    draw_box(
        ax,
        (5.85, 3.2),
        (2.0, 1.0),
        "EPDRI collapse\npooled prior\nmodifier discount\nevidence cap",
        "#FCE8E8",
    )
    draw_box(
        ax,
        (8.35, 3.2),
        (1.35, 1.0),
        "Output\nranking\nuncertainty\nreport",
        COLORS["light_gray"],
    )

    arrow(ax, (2.45, 4.6), (3.1, 4.6))
    arrow(ax, (2.45, 2.7), (3.1, 2.7))
    arrow(ax, (5.2, 4.6), (5.85, 3.9))
    arrow(ax, (5.2, 2.7), (5.85, 3.45))
    arrow(ax, (7.85, 3.7), (8.35, 3.7))
    arrow(ax, (1.4, 4.15), (1.4, 3.15))
    arrow(ax, (4.15, 4.15), (4.15, 3.15))

    ax.text(5.0, 5.55, "EGFR-PACC digital-cell framework", ha="center", va="center", fontsize=11, weight="bold")
    ax.text(
        5.0,
        0.65,
        "Local outcomes are held out from prior construction; back-testing evaluates locked literature rows by leave-one-out.",
        ha="center",
        va="center",
        fontsize=7.5,
        color="#444444",
    )
    save_figure(fig, "figure1_framework")
    plt.close(fig)


def figure2_demo_ranking() -> None:
    ranking_path = ROOT / "reports" / "demo_case_ranking_pymc_formal_4chains_2000draws.csv"
    ranking = pd.read_csv(ranking_path).head(9).sort_values("epdri_score", ascending=True)

    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    grade_colors = {"A": "#1A9850", "B": "#66BD63", "C": "#FDAE61", "D": "#D73027"}
    colors = ranking["recommendation_grade"].map(grade_colors).fillna("#999999")
    ax.barh(ranking["drug"].str.replace("_", " "), ranking["epdri_score"], color=colors, edgecolor="#333333", lw=0.5)
    for y, (_, row) in enumerate(ranking.iterrows()):
        ax.text(
            row["epdri_score"] + 0.01,
            y,
            f"{row['epdri_score']:.3f} | {row['recommendation_grade']} | {row['uncertainty_tier']}",
            va="center",
            fontsize=7,
        )
    ax.set_xlim(0, max(0.72, ranking["epdri_score"].max() + 0.16))
    ax.set_xlabel("Collapsed EPDRI score")
    ax.set_title("Demo case: G719S + L861Q, CNS disease, TP53 mutation", loc="left", weight="bold", pad=12)
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.set_axisbelow(True)
    ax.text(
        0,
        -0.18,
        "Scores use the formal 4-chain PyMC posterior prior table followed by clinical modifiers and evidence caps.",
        transform=ax.transAxes,
        fontsize=7,
        color="#444444",
    )
    save_figure(fig, "figure2_demo_case_ranking")
    plt.close(fig)


def ci_lookup(ci: pd.DataFrame, method: str, metric: str) -> pd.Series:
    row = ci[(ci["method"].eq(method)) & (ci["metric"].eq(metric))]
    if row.empty:
        raise ValueError(f"Missing CI for {method}/{metric}")
    return row.iloc[0]


def figure3_backtest() -> None:
    metrics = pd.read_csv(ROOT / "reports/backtest_loo/loo_metric_summary.csv")
    ci = pd.read_csv(ROOT / "reports/backtest_loo/loo_bootstrap_ci.csv")
    ablation = pd.read_csv(ROOT / "reports/backtest_loo/loo_ablation_summary.csv")
    subgroup = pd.read_csv(ROOT / "reports/backtest_loo/loo_subgroup_rank_summary.csv")

    fig = plt.figure(figsize=(7.2, 6.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.38, wspace=0.35)
    ax_auc = fig.add_subplot(gs[0, 0])
    ax_delta = fig.add_subplot(gs[0, 1])
    ax_ablation = fig.add_subplot(gs[1, 0])
    ax_rank = fig.add_subplot(gs[1, 1])

    method_order = ["baseline_0", "baseline_1", "model"]
    labels = ["Baseline-0", "Baseline-1", "Model"]
    auc_values = [float(metrics.loc[metrics["method"].eq(m), "pooled_drug_auc"].iloc[0]) for m in method_order]
    auc_err_low = [v - float(ci_lookup(ci, m, "pooled_drug_auc")["ci95_lower"]) for v, m in zip(auc_values, method_order)]
    auc_err_high = [float(ci_lookup(ci, m, "pooled_drug_auc")["ci95_upper"]) - v for v, m in zip(auc_values, method_order)]
    ax_auc.bar(labels, auc_values, color=[COLORS["baseline_0"], COLORS["baseline_1"], COLORS["model"]], edgecolor="#333333", lw=0.5)
    ax_auc.errorbar(labels, auc_values, yerr=[auc_err_low, auc_err_high], fmt="none", ecolor="#222222", capsize=3, lw=0.8)
    ax_auc.set_ylim(0.62, 0.98)
    ax_auc.set_ylabel("Pooled drug AUC")
    ax_auc.set_title("A. Leave-one-out AUC", loc="left", weight="bold")
    ax_auc.grid(axis="y", color="#DDDDDD", lw=0.6)

    delta_rows = [
        ("Model vs B0", ci_lookup(ci, "model", "delta_auc_vs_baseline_0"), COLORS["model"]),
        ("Model vs B1", ci_lookup(ci, "model", "delta_auc_vs_baseline_1"), COLORS["accent"]),
        ("B1 vs B0", ci_lookup(ci, "baseline_1", "delta_auc_vs_baseline_0"), COLORS["baseline_1"]),
    ]
    y = list(range(len(delta_rows)))
    means = [float(row["mean"]) for _, row, _ in delta_rows]
    low = [m - float(row["ci95_lower"]) for m, (_, row, _) in zip(means, delta_rows)]
    high = [float(row["ci95_upper"]) - m for m, (_, row, _) in zip(means, delta_rows)]
    ax_delta.axvline(0, color="#333333", lw=0.8)
    ax_delta.errorbar(means, y, xerr=[low, high], fmt="o", color="#222222", ecolor="#222222", capsize=3, lw=0.8)
    for yy, (label, _, color) in zip(y, delta_rows):
        ax_delta.scatter(means[yy], yy, color=color, s=28, zorder=3)
    ax_delta.set_yticks(y, [label for label, _, _ in delta_rows])
    ax_delta.set_xlabel("ΔAUC with bootstrap 95% CI")
    ax_delta.set_title("B. Incremental AUC", loc="left", weight="bold")
    ax_delta.grid(axis="x", color="#DDDDDD", lw=0.6)

    abl_model = ablation[ablation["method"].eq("model")].copy()
    scenario_labels = {
        "full_model": "Full",
        "no_oncokb": "No\nOncoKB",
        "no_exact_evidence_floor": "No exact\nevidence floor",
        "no_config_fallback": "No config\nfallback",
        "pacc_only_subset": "PACC-only",
    }
    abl_model["label"] = abl_model["scenario"].map(scenario_labels)
    abl_model = abl_model.iloc[::-1]
    ax_ablation.barh(
        abl_model["label"],
        abl_model["pooled_drug_auc"],
        color="#9ECAE1",
        edgecolor="#333333",
        lw=0.5,
    )
    ax_ablation.axvline(float(metrics.loc[metrics["method"].eq("baseline_0"), "pooled_drug_auc"].iloc[0]), color=COLORS["baseline_0"], ls="--", lw=1, label="B0 AUC")
    ax_ablation.axvline(float(metrics.loc[metrics["method"].eq("baseline_1"), "pooled_drug_auc"].iloc[0]), color=COLORS["baseline_1"], ls="--", lw=1, label="B1 AUC")
    ax_ablation.set_xlim(0.76, 0.93)
    ax_ablation.set_xlabel("Model AUC")
    ax_ablation.set_title("C. Design ablation", loc="left", weight="bold")
    ax_ablation.legend(frameon=False, loc="lower right")
    ax_ablation.grid(axis="x", color="#DDDDDD", lw=0.6)

    rank = subgroup[subgroup["method"].isin(["baseline_0", "baseline_1", "model"])].copy()
    keep = rank[rank["subgroup_type"].eq("class_group")]
    pivot = keep.pivot(index="subgroup", columns="method", values="mean_target_rank").loc[["PACC", "non_PACC"]]
    x = range(len(pivot.index))
    width = 0.25
    for offset, method, color, label in [
        (-width, "baseline_0", COLORS["baseline_0"], "B0"),
        (0, "baseline_1", COLORS["baseline_1"], "B1"),
        (width, "model", COLORS["model"], "Model"),
    ]:
        ax_rank.bar([i + offset for i in x], pivot[method], width=width, color=color, edgecolor="#333333", lw=0.5, label=label)
    ax_rank.set_xticks(list(x), pivot.index)
    ax_rank.set_ylabel("Mean target rank")
    ax_rank.set_title("D. Rank by mutation class", loc="left", weight="bold")
    ax_rank.legend(frameon=False, ncol=3, loc="upper left")
    ax_rank.grid(axis="y", color="#DDDDDD", lw=0.6)

    fig.suptitle("Back-testing performance of EGFR-PACC digital-cell priors", x=0.03, y=0.99, ha="left", fontsize=11, weight="bold")
    fig.text(
        0.03,
        0.01,
        "B0: Robichaux four-class baseline. B1: B0 plus OncoKB actionability. CI from 1000 fold-level bootstrap resamples.",
        fontsize=7,
        color="#444444",
    )
    save_figure(fig, "figure3_backtesting")
    plt.close(fig)


def main() -> None:
    setup_style()
    figure1_framework()
    figure2_demo_ranking()
    figure3_backtest()
    print(
        {
            "figure_dir": str(FIG_DIR.resolve()),
            "figures": [
                "figure1_framework.pdf/svg/png",
                "figure2_demo_case_ranking.pdf/svg/png",
                "figure3_backtesting.pdf/svg/png",
            ],
        }
    )


if __name__ == "__main__":
    main()
