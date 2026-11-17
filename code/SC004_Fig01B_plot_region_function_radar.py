from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"D:\桌面\PDE2A")
INPUT = ROOT / "results" / "Spatial_source_annotation_recovery_20260730"
OUTPUT = ROOT / "results" / "Spatial_source_annotation_function_redraw_20260805"
OUTPUT.mkdir(parents=True, exist_ok=True)

SEED = 20260805
BLUE = "#0072B2"
VERMILLION = "#D55E00"
PURPLE = "#7C3AED"
GRAY = "#7A7A7A"
LIGHT_GRAY = "#D9D9D9"

LABELS = {
    "1.Hypoxia": "Hypoxia",
    "2.Angiogenesis": "Angiogenesis",
    "3.Glycolysis": "Glycolysis",
    "4.Cell_cycle": "Cell cycle",
    "5.DNA_replication": "DNA replication",
    "6.P53_signaling_pathway": "p53 signaling",
    "7.Pyrimidine_metabolism": "Pyrimidine metabolism",
    "8.Stem_score": "Stemness",
    "9.EMT": "EMT",
}

EFFECT_ORDER = [
    "9.EMT",
    "4.Cell_cycle",
    "5.DNA_replication",
    "8.Stem_score",
    "2.Angiogenesis",
    "7.Pyrimidine_metabolism",
    "3.Glycolysis",
    "6.P53_signaling_pathway",
    "1.Hypoxia",
]

RADAR_ORDER = [
    "8.Stem_score",
    "9.EMT",
    "4.Cell_cycle",
    "5.DNA_replication",
    "7.Pyrimidine_metabolism",
    "3.Glycolysis",
    "2.Angiogenesis",
    "1.Hypoxia",
    "6.P53_signaling_pathway",
]


def significance(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def bootstrap_median_ci(values, rng, n_boot=20000):
    values = np.asarray(values, dtype=float)
    boot = np.median(
        rng.choice(values, size=(n_boot, values.size), replace=True), axis=1
    )
    return np.quantile(boot, [0.025, 0.975])


def load_and_summarize():
    sample = pd.read_csv(INPUT / "source_author_region_function_by_sample.csv")
    stats = pd.read_csv(INPUT / "source_author_region_function_paired_stats.csv")
    sample["label"] = sample["program"].map(LABELS)

    rng = np.random.default_rng(SEED)
    records = []
    for program in EFFECT_ORDER:
        values = sample.loc[sample["program"] == program, "edge_minus_core"].to_numpy()
        lo, hi = bootstrap_median_ci(values, rng)
        stat = stats.loc[stats["program"] == program].iloc[0]
        records.append(
            {
                "program": program,
                "label": LABELS[program],
                "n_pairs": int(stat["n_pairs"]),
                "median_edge_minus_core": float(np.median(values)),
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
                "wilcoxon_p": float(stat["p"]),
                "bh_q": float(stat["q_bh"]),
                "significance": significance(float(stat["q_bh"])),
            }
        )
    summary = pd.DataFrame(records)
    return sample, summary


def draw_effect_plot(sample, summary):
    fig = plt.figure(figsize=(11.6, 7.1), constrained_layout=False)
    grid = fig.add_gridspec(1, 2, width_ratios=[5.6, 1.05], wspace=0.03)
    ax = fig.add_subplot(grid[0, 0])
    ax_q = fig.add_subplot(grid[0, 1], sharey=ax)
    y = np.arange(len(EFFECT_ORDER))

    ax.axvline(0, color="#333333", linewidth=1.0, linestyle=(0, (3, 3)), zorder=0)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.7, zorder=0)

    for i, program in enumerate(EFFECT_ORDER):
        row = summary.loc[summary["program"] == program].iloc[0]
        vals = sample.loc[sample["program"] == program, "edge_minus_core"].to_numpy()
        color = VERMILLION if row["bh_q"] < 0.05 else GRAY
        jitter = np.linspace(-0.11, 0.11, len(vals))
        ax.scatter(
            vals,
            np.full(len(vals), i) + jitter,
            s=24,
            facecolors="white",
            edgecolors=color,
            linewidths=0.8,
            alpha=0.72,
            zorder=2,
        )
        ax.plot(
            [row["bootstrap_ci_low"], row["bootstrap_ci_high"]],
            [i, i],
            color=color,
            linewidth=2.0,
            solid_capstyle="round",
            zorder=3,
        )
        ax.scatter(
            row["median_edge_minus_core"],
            i,
            marker="D",
            s=63,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )

    ax.set_yticks(y, [LABELS[p] for p in EFFECT_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Edge − Core program score", labelpad=9)
    ax.set_ylabel("")
    ax.set_xlim(-0.34, 0.47)
    ax.set_title(
        "Original-author spatial regions: paired functional effects",
        loc="left",
        fontsize=15,
        fontweight="bold",
        pad=19,
    )
    ax.text(
        0,
        1.015,
        "Pathology SCC gate; section-level pairing (n = 12). Diamonds show medians; lines show bootstrap 95% CIs.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.6,
        color="#4D4D4D",
    )

    ax_q.set_xlim(0, 1)
    ax_q.text(0.12, -0.72, "BH q", fontsize=9.2, fontweight="bold", ha="left", va="center")
    for i, program in enumerate(EFFECT_ORDER):
        row = summary.loc[summary["program"] == program].iloc[0]
        color = VERMILLION if row["bh_q"] < 0.05 else GRAY
        ax_q.text(0.12, i, f"{row['bh_q']:.3f}", fontsize=9.0, ha="left", va="center", color=color)
        ax_q.text(0.92, i, row["significance"], fontsize=10, fontweight="bold", ha="right", va="center", color=color)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="y", length=0, pad=7)
    ax.tick_params(axis="x", width=0.8)
    ax_q.set_axis_off()
    fig.subplots_adjust(left=0.20, right=0.96, top=0.84, bottom=0.13)

    out = OUTPUT / "Fig_source_author_regions_paired_function_effects.png"
    fig.savefig(out, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def make_radar_data(sample):
    rows = []
    for program in RADAR_ORDER:
        subset = sample.loc[sample["program"] == program, ["sample", "core", "edge"]].copy()
        pooled = np.r_[subset["core"].to_numpy(), subset["edge"].to_numpy()]
        low, high = pooled.min(), pooled.max()
        span = high - low
        if span == 0:
            subset["core_scaled"] = 0.5
            subset["edge_scaled"] = 0.5
        else:
            subset["core_scaled"] = (subset["core"] - low) / span
            subset["edge_scaled"] = (subset["edge"] - low) / span
        for region in ["Core", "Edge"]:
            col = f"{region.lower()}_scaled"
            rows.append(
                {
                    "program": program,
                    "label": LABELS[program],
                    "region": region,
                    "median_scaled_score": subset[col].median(),
                    "minimum_used_for_scaling": low,
                    "maximum_used_for_scaling": high,
                }
            )
    return pd.DataFrame(rows)


def draw_radar(radar, summary):
    sig_map = summary.set_index("program")["significance"].to_dict()
    label_text = {
        "8.Stem_score": "Stemness",
        "9.EMT": "EMT",
        "4.Cell_cycle": "Cell cycle",
        "5.DNA_replication": "DNA replication",
        "7.Pyrimidine_metabolism": "Pyrimidine\nmetabolism",
        "3.Glycolysis": "Glycolysis",
        "2.Angiogenesis": "Angiogenesis",
        "1.Hypoxia": "Hypoxia",
        "6.P53_signaling_pathway": "p53 signaling",
    }
    n = len(RADAR_ORDER)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    closed_angles = np.r_[angles, angles[0]]

    fig = plt.figure(figsize=(10.6, 9.2))
    ax = fig.add_axes([0.20, 0.20, 0.60, 0.60], projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    for region, color, marker in [("Core", BLUE, "o"), ("Edge", VERMILLION, "s")]:
        values = (
            radar.loc[radar["region"] == region]
            .set_index("program")
            .loc[RADAR_ORDER, "median_scaled_score"]
            .to_numpy()
        )
        closed = np.r_[values, values[0]]
        ax.plot(
            closed_angles,
            closed,
            color=color,
            linewidth=2.25,
            marker=marker,
            markersize=5.3,
            markerfacecolor="white",
            markeredgewidth=1.4,
            label=region,
            zorder=3,
        )

    ax.set_xticks(angles)
    ax.set_xticklabels([])
    for angle, program in zip(angles, RADAR_ORDER):
        stars = sig_map[program] if sig_map[program] != "ns" else ""
        text = f"{label_text[program]}{f'  {stars}' if stars else ''}"
        horizontal = np.sin(angle)
        vertical = np.cos(angle)
        if horizontal > 0.18:
            ha = "left"
        elif horizontal < -0.18:
            ha = "right"
        else:
            ha = "center"
        if vertical > 0.28:
            va = "bottom"
        elif vertical < -0.28:
            va = "top"
        else:
            va = "center"
        ax.text(
            angle,
            1.17,
            text,
            ha=ha,
            va=va,
            fontsize=9.6,
            linespacing=1.08,
            clip_on=False,
        )
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=8.3, color="#666666")
    ax.set_rlabel_position(72)
    ax.yaxis.grid(True, color=LIGHT_GRAY, linewidth=0.8)
    ax.xaxis.grid(True, color="#E5E5E5", linewidth=0.7)
    ax.spines["polar"].set_color("#BDBDBD")
    ax.spines["polar"].set_linewidth(0.8)

    fig.text(
        0.5,
        0.965,
        "Core–Edge functional profile",
        ha="center",
        va="top",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.920,
        "Median across 12 sections; within-program min–max scaling",
        ha="center",
        va="center",
        fontsize=9.8,
        color="#4D4D4D",
    )
    fig.text(
        0.5,
        0.025,
        "Descriptive display only; inferential results are reported in the paired-effect plot.  * BH q < 0.05; ** q < 0.01; *** q < 0.001",
        ha="center",
        va="center",
        fontsize=8.6,
        color="#4D4D4D",
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.23),
        ncol=2,
        frameon=False,
        fontsize=10,
        handlelength=2.4,
    )
    out = OUTPUT / "Fig_source_author_regions_function_radar.png"
    fig.savefig(out, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.linewidth": 0.8,
            "axes.unicode_minus": False,
            "savefig.transparent": False,
        }
    )
    sample, summary = load_and_summarize()
    radar = make_radar_data(sample)

    sample.to_csv(OUTPUT / "source_author_function_paired_deltas.csv", index=False)
    summary.to_csv(OUTPUT / "source_author_function_effect_summary.csv", index=False)
    radar.to_csv(OUTPUT / "source_author_function_radar_scaled.csv", index=False)

    effect_path = draw_effect_plot(sample, summary)
    radar_path = draw_radar(radar, summary)
    print(effect_path)
    print(radar_path)


if __name__ == "__main__":
    main()
