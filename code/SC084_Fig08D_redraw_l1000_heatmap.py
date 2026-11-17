#!/usr/bin/env python3
"""Redraw the standalone Wide L1000CDS2 heatmap as a narrow, tall panel."""

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

from build_l1000_heatmap_panels_bc import (
    BLUE,
    NAVY,
    PURPLE,
    RED,
    build_matrix,
    display_name,
    load_query,
    make_cell_palette,
)


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "L1000CDS2_wide_narrow_tall_large_type_20261107"
STEM = OUTDIR / "Fig_L1000CDS2_wide_narrow_tall_large_type"


def draw() -> None:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 9.0,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    strict_conditions, _, _ = load_query("Strict")
    conditions, up_genes, down_genes = load_query("Wide")
    rows_up, rows_down, matrix = build_matrix(conditions, up_genes, down_genes)
    genes = rows_up + rows_down
    palette = make_cell_palette(
        [str(item["cell_id"]) for item in strict_conditions + conditions]
    )

    if matrix.shape != (96, 50):
        raise ValueError(f"Expected a 96 x 50 matrix, got {matrix.shape}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.4, 17.0), facecolor="white")
    grid = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.90, 0.24, 13.55],
        hspace=0.035,
        left=0.205,
        right=0.985,
        top=0.992,
        bottom=0.105,
    )
    ax_drug = fig.add_subplot(grid[0])
    ax_cell = fig.add_subplot(grid[1])
    ax = fig.add_subplot(grid[2])

    n_cols = len(conditions)
    ax_drug.set_xlim(-0.5, n_cols - 0.5)
    ax_drug.set_ylim(0, 1)
    for index, condition in enumerate(conditions):
        is_cgp = str(condition["pert_desc"]).lower() == "cgp-60474"
        ax_drug.text(
            index,
            0.02,
            display_name(condition["pert_desc"]),
            rotation=90,
            rotation_mode="anchor",
            ha="left",
            va="center",
            fontsize=8.6,
            linespacing=0.92,
            color=PURPLE if is_cgp else NAVY,
            weight="bold" if is_cgp else "normal",
            clip_on=True,
        )
    ax_drug.set_xticks([])
    ax_drug.set_yticks([])
    for spine in ax_drug.spines.values():
        spine.set_visible(False)

    cells = [str(item["cell_id"]) for item in conditions]
    rgb = np.array([[mpl.colors.to_rgb(palette[cell]) for cell in cells]])
    ax_cell.imshow(rgb, aspect="auto", interpolation="none")
    ax_cell.set_ylabel(
        "Cell line",
        rotation=0,
        ha="right",
        va="center",
        fontsize=10.5,
        color=NAVY,
        labelpad=8,
    )
    ax_cell.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_cell.grid(which="minor", axis="x", color="white", linewidth=0.85)
    ax_cell.tick_params(which="minor", bottom=False)
    ax_cell.set_xticks([])
    ax_cell.set_yticks([])
    for spine in ax_cell.spines.values():
        spine.set_linewidth(1.25)
        spine.set_color("#C9D0D6")

    cmap = ListedColormap([BLUE, "#FFFFFF", RED])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    ax.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="#E6EBEF", linewidth=0.78)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=9.3)
    ax.tick_params(axis="y", length=0, pad=5, colors=NAVY)
    for label in ax.get_yticklabels():
        if label.get_text() == "PDE2A":
            label.set_color(PURPLE)
            label.set_weight("bold")
    ax.axhline(len(rows_up) - 0.5, color=NAVY, linewidth=2.2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.25)
        spine.set_color("#B8C2CA")

    direction_handles = [
        Patch(facecolor=BLUE, edgecolor=NAVY, linewidth=0.8,
              label="Loss-up reversed down"),
        Patch(facecolor=RED, edgecolor="#8F3C3C", linewidth=0.8,
              label="Loss-down reversed up"),
    ]
    cell_handles = [
        Patch(facecolor=color, edgecolor="#666666", linewidth=0.55, label=cell)
        for cell, color in palette.items()
    ]
    legend1 = fig.legend(
        handles=direction_handles,
        loc="lower left",
        bbox_to_anchor=(0.20, 0.025),
        frameon=False,
        ncol=1,
        fontsize=9.8,
        labelspacing=0.65,
        handlelength=1.35,
        handleheight=0.95,
    )
    fig.add_artist(legend1)
    fig.legend(
        handles=cell_handles,
        title="Cell line",
        loc="lower right",
        bbox_to_anchor=(0.985, 0.020),
        frameon=False,
        ncol=4,
        fontsize=9.0,
        title_fontsize=9.4,
        columnspacing=0.65,
        handlelength=0.9,
        handletextpad=0.30,
        labelspacing=0.35,
    )

    for suffix in ("png", "svg", "pdf"):
        fig.savefig(
            STEM.with_suffix(f".{suffix}"),
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
            pad_inches=0.06,
            facecolor="white",
        )
    plt.close(fig)

    np.savetxt(OUTDIR / "reverse_overlap_matrix.tsv", matrix, fmt="%d", delimiter="\t")
    with (OUTDIR / "row_order.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["row_index", "gene", "direction"])
        for index, gene in enumerate(genes, start=1):
            direction = "loss_up_reversed_down" if index <= len(rows_up) else "loss_down_reversed_up"
            writer.writerow([index, gene, direction])
    with (OUTDIR / "condition_order.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["column_index", "pert_desc", "cell_id", "is_CGP_60474"])
        for index, condition in enumerate(conditions, start=1):
            name = str(condition["pert_desc"])
            writer.writerow([index, name, condition["cell_id"], name.lower() == "cgp-60474"])


if __name__ == "__main__":
    draw()
