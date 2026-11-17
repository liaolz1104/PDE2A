"""Render strict/wide L1000CDS2 heatmaps as clean publication panels B and C."""

import json
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "mechanism_drug_work" / "results"
OUTDIR = ROOT / "results" / "L1000CDS2_heatmap_panels_BC_20261002"
OUTPUT = OUTDIR / "Fig_L1000CDS2_heatmap_panels_BC.png"

CONFIG = {
    "Strict": {
        "conditions": RESULTS / "l1000cds2_reverse" / "l1000cds2_reverse_conditions.json",
        "up": RESULTS / "vendor_preliminary" / "pde2a_loss_up.grp",
        "down": RESULTS / "vendor_preliminary" / "pde2a_loss_down.grp",
    },
    "Wide": {
        "conditions": RESULTS / "l1000cds2_reverse_wide" / "l1000cds2_reverse_wide_conditions.json",
        "payload": RESULTS / "l1000cds2_reverse_wide" / "l1000cds2_reverse_wide_payload.json",
    },
}

PURPLE = "#7C3AED"
BLUE = "#2C7FB8"
RED = "#D95F5F"
WHITE = "#FFFFFF"
NAVY = "#17324D"


def load_query(query):
    config = CONFIG[query]
    conditions = json.loads(config["conditions"].read_text(encoding="utf-8"))
    if query == "Strict":
        up = config["up"].read_text(encoding="utf-8").split()
        down = config["down"].read_text(encoding="utf-8").split()
    else:
        payload = json.loads(config["payload"].read_text(encoding="utf-8"))
        up = payload["data"]["upGenes"]
        down = payload["data"]["dnGenes"]
    return conditions, up, down


def build_matrix(conditions, up_genes, down_genes):
    observed_up = {gene for item in conditions for gene in item.get("overlap", {}).get("up/dn", [])}
    observed_down = {gene for item in conditions for gene in item.get("overlap", {}).get("dn/up", [])}
    rows_up = [gene for gene in up_genes if gene in observed_up]
    rows_down = [gene for gene in down_genes if gene in observed_down]
    genes = rows_up + rows_down
    row_index = {gene: index for index, gene in enumerate(genes)}
    matrix = np.zeros((len(genes), len(conditions)), dtype=np.int8)
    for column, condition in enumerate(conditions):
        overlap = condition.get("overlap", {})
        for gene in overlap.get("up/dn", []):
            if gene in row_index:
                matrix[row_index[gene], column] = -1
        for gene in overlap.get("dn/up", []):
            if gene in row_index:
                matrix[row_index[gene], column] = 1
    return rows_up, rows_down, matrix


def display_name(name):
    clean = " ".join(str(name).replace('"', "").split())
    if len(clean) > 52:
        return "\n".join(textwrap.wrap(clean, width=52, break_long_words=False,
                                       break_on_hyphens=False))
    return clean


def make_cell_palette(all_cells):
    unique = list(dict.fromkeys(all_cells))
    colors = list(plt.get_cmap("tab20").colors)
    return {cell: colors[index % len(colors)] for index, cell in enumerate(unique)}


def draw_drug_header(ax, conditions):
    n = len(conditions)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    for index, condition in enumerate(conditions):
        is_cgp = str(condition["pert_desc"]).lower() == "cgp-60474"
        ax.text(index, 0.02, display_name(condition["pert_desc"]), rotation=90,
                rotation_mode="anchor", ha="left", va="center", fontsize=4.6,
                linespacing=0.88, color=PURPLE if is_cgp else NAVY,
                weight="bold" if is_cgp else "normal", clip_on=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_cell_strip(ax, conditions, palette):
    cells = [str(item["cell_id"]) for item in conditions]
    rgb = np.array([[mpl.colors.to_rgb(palette[cell]) for cell in cells]])
    ax.imshow(rgb, aspect="auto", interpolation="none")
    ax.set_ylabel("Cell line", rotation=0, ha="right", va="center",
                  fontsize=6.5, color=NAVY, labelpad=6)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_heatmap(ax, rows_up, rows_down, matrix):
    genes = rows_up + rows_down
    cmap = ListedColormap([BLUE, WHITE, RED])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    ax.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#F1F4F6", linewidth=0.25)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=5.3)
    ax.tick_params(axis="y", length=0, pad=3, colors=NAVY)
    for label in ax.get_yticklabels():
        if label.get_text() == "PDE2A":
            label.set_color(PURPLE)
            label.set_weight("bold")
    if rows_up and rows_down:
        ax.axhline(len(rows_up) - 0.5, color=NAVY, linewidth=0.85)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.unicode_minus": False,
    })
    strict_conditions, strict_up, strict_down = load_query("Strict")
    wide_conditions, wide_up, wide_down = load_query("Wide")
    strict_rows_up, strict_rows_down, strict_matrix = build_matrix(strict_conditions, strict_up, strict_down)
    wide_rows_up, wide_rows_down, wide_matrix = build_matrix(wide_conditions, wide_up, wide_down)

    all_cells = [str(item["cell_id"]) for item in strict_conditions + wide_conditions]
    palette = make_cell_palette(all_cells)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    grid = fig.add_gridspec(
        8, 1,
        height_ratios=[1.55, 0.13, 2.15, 0.22, 0.72, 0.13, 8.72, 0.48],
        hspace=0.06, left=0.12, right=0.992, top=0.99, bottom=0.025,
    )

    ax_b_drug = fig.add_subplot(grid[0])
    ax_b_cell = fig.add_subplot(grid[1])
    ax_b = fig.add_subplot(grid[2])
    ax_gap = fig.add_subplot(grid[3])
    ax_c_drug = fig.add_subplot(grid[4])
    ax_c_cell = fig.add_subplot(grid[5])
    ax_c = fig.add_subplot(grid[6])
    ax_legend = fig.add_subplot(grid[7])

    draw_drug_header(ax_b_drug, strict_conditions)
    draw_cell_strip(ax_b_cell, strict_conditions, palette)
    draw_heatmap(ax_b, strict_rows_up, strict_rows_down, strict_matrix)
    draw_drug_header(ax_c_drug, wide_conditions)
    draw_cell_strip(ax_c_cell, wide_conditions, palette)
    draw_heatmap(ax_c, wide_rows_up, wide_rows_down, wide_matrix)

    ax_gap.axis("off")
    ax_legend.axis("off")
    ax_b_drug.text(-0.105, 0.92, "B", transform=ax_b_drug.transAxes,
                   fontsize=10, weight="bold", color="black", ha="left", va="top")
    ax_b_drug.text(-0.065, 0.92, "Strict", transform=ax_b_drug.transAxes,
                   fontsize=7.0, weight="bold", color=NAVY, ha="left", va="top")
    ax_c_drug.text(-0.105, 0.92, "C", transform=ax_c_drug.transAxes,
                   fontsize=10, weight="bold", color="black", ha="left", va="top")
    ax_c_drug.text(-0.065, 0.92, "Wide", transform=ax_c_drug.transAxes,
                   fontsize=7.0, weight="bold", color=NAVY, ha="left", va="top")

    direction_handles = [
        Patch(facecolor=BLUE, edgecolor="none", label="Loss-up reversed down"),
        Patch(facecolor=RED, edgecolor="none", label="Loss-down reversed up"),
    ]
    cell_handles = [Patch(facecolor=color, edgecolor="none", label=cell)
                    for cell, color in palette.items()]
    legend1 = ax_legend.legend(
        handles=direction_handles, loc="upper left", bbox_to_anchor=(0.0, 0.95),
        frameon=False, ncol=3, fontsize=5.8, columnspacing=1.0, handlelength=1.2,
    )
    ax_legend.add_artist(legend1)
    ax_legend.legend(
        handles=cell_handles, title="Cell line", loc="upper right", bbox_to_anchor=(1.0, 0.98),
        frameon=False, ncol=8, fontsize=5.0, title_fontsize=5.4,
        columnspacing=0.55, handlelength=0.9, handletextpad=0.25,
    )

    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    draw()
