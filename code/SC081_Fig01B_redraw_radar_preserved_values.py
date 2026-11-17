#!/usr/bin/env python3
"""Restyle the supplied Core/Edge radar without changing its 18 values.

The supplied raster used twice the archived within-program min-max scaled
medians.  This script reuses that authoritative table and changes styling only.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "results" / "Spatial_source_annotation_function_redraw_20260805"
INPUT_VALUES = SOURCE_DIR / "source_author_function_radar_scaled.csv"
INPUT_STATS = SOURCE_DIR / "source_author_function_effect_summary.csv"
OUTPUT_DIR = ROOT / "results" / "Spatial_radar_preserved_values_large_filled_20261107"

ORDER = [
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

LABELS = {
    "8.Stem_score": "Stemness",
    "9.EMT": "EMT",
    "4.Cell_cycle": "Cell cycle",
    "5.DNA_replication": "DNA replication",
    "7.Pyrimidine_metabolism": "Pyrimidine metabolism",
    "3.Glycolysis": "Glycolysis",
    "2.Angiogenesis": "Angiogenesis",
    "1.Hypoxia": "Hypoxia",
    "6.P53_signaling_pathway": "p53 signaling",
}

BLUE = "#0B63A9"
ORANGE = "#F36C0A"

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 24,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def load_values() -> tuple[pd.DataFrame, pd.DataFrame]:
    values = pd.read_csv(INPUT_VALUES)
    stats = pd.read_csv(INPUT_STATS)
    values["radius"] = 2.0 * values["median_scaled_score"]
    expected = {(program, region) for program in ORDER for region in ("Core", "Edge")}
    actual = set(zip(values["program"], values["region"]))
    if expected - actual:
        raise ValueError(f"Missing radar values: {sorted(expected - actual)}")
    return values, stats


def draw(values: pd.DataFrame, stats: pd.DataFrame) -> tuple[Path, Path, Path]:
    stars = stats.set_index("program")["significance"].to_dict()
    angles = np.linspace(0, 2 * np.pi, len(ORDER), endpoint=False)
    closed_angles = np.r_[angles, angles[0]]

    fig = plt.figure(figsize=(10.4, 8.2), facecolor="white")
    ax = fig.add_axes([0.19, 0.18, 0.62, 0.72], projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    for region, color in (("Core", BLUE), ("Edge", ORANGE)):
        radii = (
            values.loc[values["region"] == region]
            .set_index("program")
            .loc[ORDER, "radius"]
            .to_numpy()
        )
        ax.fill(
            closed_angles,
            np.r_[radii, radii[0]],
            color=color,
            alpha=0.18,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            closed_angles,
            np.r_[radii, radii[0]],
            color=color,
            linewidth=6.2,
            marker="o",
            markersize=13.5,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=2.4,
            solid_capstyle="round",
            solid_joinstyle="round",
            label=region,
            zorder=4,
        )

    ax.set_xticks(angles)
    ax.set_xticklabels([])
    for angle, program in zip(angles, ORDER):
        sig = stars.get(program, "ns")
        suffix = f"  $\\mathbf{{{sig}}}$" if sig and sig != "ns" else ""
        horizontal = np.sin(angle)
        vertical = np.cos(angle)
        ha = "left" if horizontal > 0.18 else "right" if horizontal < -0.18 else "center"
        va = "bottom" if vertical > 0.28 else "top" if vertical < -0.28 else "center"
        ax.text(
            angle,
            2.22,
            f"{LABELS[program]}{suffix}",
            ha=ha,
            va=va,
            fontsize=24,
            clip_on=False,
        )

    ax.set_ylim(0, 2.0)
    ax.set_yticks([0.5, 1.0, 1.5, 2.0])
    ax.set_yticklabels(["0.5", "1.0", "1.5", "2.0"], fontsize=19, color="#333333")
    ax.set_rlabel_position(67)
    ax.text(np.deg2rad(145), 0.18, "0", fontsize=19, color="#333333", ha="center", va="center")

    ax.yaxis.grid(True, color="#929292", linewidth=2.35, linestyle="--", dashes=(2.4, 2.4))
    ax.xaxis.grid(True, color="#929292", linewidth=2.25, linestyle="--", dashes=(2.4, 2.4))
    ax.spines["polar"].set_color("#777777")
    ax.spines["polar"].set_linewidth(2.8)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=2,
        frameon=False,
        fontsize=22,
        handlelength=2.0,
        handletextpad=0.35,
        columnspacing=0.3,
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT_DIR / "Fig_source_author_regions_function_radar_preserved_values_large_filled"
    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(svg, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(pdf, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return png, svg, pdf


def main() -> None:
    values, stats = load_values()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    audit = values.loc[:, ["program", "label", "region", "median_scaled_score", "radius"]].copy()
    audit["formula"] = "2 * median_scaled_score"
    audit.to_csv(OUTPUT_DIR / "radar_values_preserved.csv", index=False)
    pd.DataFrame({
        "element": ["program labels", "radial ticks", "legend", "series lines", "markers", "grid lines", "outer frame"],
        "size_or_width": [24.0, 19.0, 22.0, 6.2, 13.5, 2.3, 2.8],
        "unit": ["pt", "pt", "pt", "pt", "pt", "pt", "pt"],
    }).to_csv(OUTPUT_DIR / "STYLE_AUDIT.tsv", sep="\t", index=False)
    png, svg, pdf = draw(values, stats)
    pivot = audit.pivot(index="program", columns="region", values="radius").loc[ORDER]
    print(pivot.to_string())
    print(f"PNG: {png}")
    print(f"SVG: {svg}")
    print(f"PDF: {pdf}")


if __name__ == "__main__":
    main()
