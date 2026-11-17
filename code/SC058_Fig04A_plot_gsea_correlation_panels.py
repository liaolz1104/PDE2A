#!/usr/bin/env python3
"""Render the audited GSEA plus PDE2A-correlation figure on the analysis server.

The script reads existing ranked genes, GSEA statistics and correlation source
tables. It does not rerun differential expression, GSEA permutations or the
correlation analysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


GSEA_GENES = {
    "Proliferation_cell_cycle": [
        "MKI67", "TOP2A", "PCNA", "MCM2", "MCM4", "MCM6", "CDK1",
        "CCNB1", "CCNB2", "BUB1", "AURKB", "TYMS", "UBE2C",
    ],
    "Stemness_CSC": [
        "SOX2", "KLF4", "ALDH1A1", "PROM1", "CD44", "ITGA6", "LGR5",
        "EPCAM", "MKI67", "TOP2A", "BIRC5", "HMGA2", "KRT5", "KRT14", "TP63",
    ],
}

GSEA_TITLES = {
    "Proliferation_cell_cycle": "Proliferation / cell cycle",
    "Stemness_CSC": "Stemness",
}

SUBSETS = ["all_bulk_samples", "tumor_IT_VT_only", "paired_IT_VT_patients"]
SUBSET_LABELS = {
    "all_bulk_samples": "All specimens\n(n=65)",
    "tumor_IT_VT_only": "Tumor IT/VT\n(n=35)",
    "paired_IT_VT_patients": "IT/VT-paired patients\n(n=51)",
}
DISPLAY_ORDER = [
    "stemness_basal",
    "stemness_CSC_surface",
    "proliferation",
    "EMT_invasion",
    "leading_edge",
    "tumor_core",
    "cAMP_axis_no_PDE2A",
    "MAPK_activity",
]
ROW_LABELS = {
    "stemness_basal": "Stemness basal",
    "stemness_CSC_surface": "Stemness CSC",
    "proliferation": "Proliferation",
    "EMT_invasion": "EMT/invasion",
    "leading_edge": "Leading-edge program",
    "tumor_core": "Tumor-core program",
    "cAMP_axis_no_PDE2A": "cAMP-related program",
    "MAPK_activity": "MAPK",
}


def configure_matplotlib() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.labelsize": 7.3,
        "axes.titlesize": 8.5,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#333333",
        "xtick.labelsize": 6.6,
        "ytick.labelsize": 6.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "legend.frameon": False,
        "savefig.facecolor": "white",
    })


def fdr_label(value: float) -> str:
    return "<0.001" if value == 0 or value < 0.001 else f"={value:.4f}"


def running_enrichment(rank: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    hits = rank["gene_upper"].isin({g.upper() for g in genes}).to_numpy()
    scores = rank["tstat_low_vs_high"].to_numpy(float)
    n_hits = int(hits.sum())
    if n_hits == 0 or n_hits == len(rank):
        raise ValueError("The gene set has no usable hits or covers the ranked list.")
    weights = np.abs(scores)
    hit_weights = np.where(hits, weights, 0.0)
    phit = np.cumsum(hit_weights / hit_weights.sum())
    pmiss = np.cumsum((~hits).astype(float) / (len(rank) - n_hits))
    return phit - pmiss, np.flatnonzero(hits)


def add_gsea_panel(fig, slot, rank: pd.DataFrame, stats: pd.DataFrame, term: str, label: str) -> None:
    inner = slot.subgridspec(4, 1, height_ratios=[2.45, 0.30, 0.16, 1.05], hspace=0.035)
    ax_es = fig.add_subplot(inner[0])
    ax_hits = fig.add_subplot(inner[1], sharex=ax_es)
    ax_strip = fig.add_subplot(inner[2], sharex=ax_es)
    ax_rank = fig.add_subplot(inner[3], sharex=ax_es)

    es, hit_idx = running_enrichment(rank, GSEA_GENES[term])
    x = np.arange(1, len(rank) + 1)
    score = rank["tstat_low_vs_high"].to_numpy(float)
    row = stats.loc[stats["Term"].eq(term)].iloc[0]

    purple = "#6F3FB5"
    pale_purple = "#B79BE3"
    low_color = "#C66A3D"
    high_color = "#2B7896"

    ax_es.plot(x, es, color=purple, lw=1.65, solid_capstyle="round", zorder=3)
    ax_es.fill_between(x, 0, es, color=pale_purple, alpha=0.13, linewidth=0, zorder=1)
    ax_es.axhline(0, color="#444444", lw=0.55)
    ax_es.set_ylabel("Running enrichment score (ES)", labelpad=4)
    ax_es.set_ylim(min(-0.05, float(es.min()) - 0.03), float(es.max()) + 0.12)
    ax_es.tick_params(axis="x", bottom=False, labelbottom=False, length=2)
    ax_es.text(
        0.985, 0.94,
        f"NES={float(row['NES']):.2f}, FDR q{fdr_label(float(row['FDR q-val']))}",
        transform=ax_es.transAxes, ha="right", va="top", fontsize=6.4,
    )
    ax_es.set_title(GSEA_TITLES[term], fontsize=8.6, weight="bold", pad=5)

    ax_hits.vlines(hit_idx + 1, 0, 1, color="#222222", lw=0.48)
    ax_hits.set_ylim(0, 1)
    ax_hits.axis("off")

    limit = max(float(np.nanpercentile(np.abs(score), 99)), 1e-9)
    cmap = LinearSegmentedColormap.from_list("rank_direction", [high_color, "#F7F7F7", low_color])
    ax_strip.imshow(
        np.clip(score, -limit, limit).reshape(1, -1), aspect="auto", cmap=cmap,
        vmin=-limit, vmax=limit, extent=[1, len(rank), 0, 1], interpolation="nearest",
    )
    ax_strip.axis("off")

    ax_rank.fill_between(x, score, 0, color="#B7B7B7", lw=0)
    ax_rank.axhline(0, color="#555555", lw=0.5)
    ax_rank.set_ylabel("Rank metric", labelpad=4)
    ax_rank.set_xlabel("Ranked genes", labelpad=3)
    ax_rank.tick_params(length=2, pad=2)
    for axis in (ax_es, ax_rank):
        axis.spines["left"].set_linewidth(0.65)
        axis.spines["bottom"].set_linewidth(0.65)

    ax_es.text(-0.18, 1.17, label, transform=ax_es.transAxes, fontsize=11, weight="bold", va="top")
    ax_rank.text(0.0, -0.62, "PDE2A-low", transform=ax_rank.transAxes, color=low_color,
                 fontsize=6.5, weight="bold", ha="left", va="top")
    ax_rank.text(1.0, -0.62, "PDE2A-high", transform=ax_rank.transAxes, color=high_color,
                 fontsize=6.5, weight="bold", ha="right", va="top")


def significance(q: float) -> str:
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def add_heatmap_panel(fig, slot, corr: pd.DataFrame) -> None:
    ax = fig.add_subplot(slot)
    selected = corr[corr["signature"].isin(DISPLAY_ORDER) & corr["subset"].isin(SUBSETS)].copy()
    rho = selected.pivot(index="signature", columns="subset", values="spearman_rho").loc[DISPLAY_ORDER, SUBSETS]
    qvals = selected.pivot(index="signature", columns="subset", values="fdr_within_subset_12_signatures").loc[DISPLAY_ORDER, SUBSETS]
    annotations = rho.copy().astype(object)
    for row in DISPLAY_ORDER:
        for col in SUBSETS:
            annotations.loc[row, col] = f"{rho.loc[row, col]:.2f}{significance(qvals.loc[row, col])}"
    rho.index = [ROW_LABELS[x] for x in rho.index]
    rho.columns = [SUBSET_LABELS[x] for x in rho.columns]
    annotations.index = rho.index
    annotations.columns = rho.columns

    hm = sns.heatmap(
        rho, cmap=sns.diverging_palette(245, 25, s=82, l=48, center="light", as_cmap=True),
        center=0, vmin=-1, vmax=1, annot=annotations, fmt="",
        annot_kws={"fontsize": 6.5}, linewidths=0.65, linecolor="white",
        cbar_kws={"label": "Spearman’s ρ (PDE2A expression)", "shrink": 0.78, "pad": 0.035}, ax=ax,
    )
    ax.set_title("PDE2A-associated expression programs", fontsize=8.6, weight="bold", pad=5)
    ax.set_xlabel("")
    ax.set_ylabel("Expression program", labelpad=5)
    ax.tick_params(axis="x", labelrotation=0, labelsize=5.7, length=0, pad=4)
    ax.tick_params(axis="y", labelrotation=0, labelsize=6.4, length=0, pad=3)
    for text in ax.get_yticklabels():
        text.set_ha("left")
        text.set_x(-0.02)
    cb = hm.collections[0].colorbar
    cb.set_ticks([-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=5.8, length=2.5, width=0.5)
    cb.outline.set_linewidth(0.5)
    ax.text(-0.18, 1.17, "c", transform=ax.transAxes, fontsize=11, weight="bold", va="top")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gsea-rank", type=Path, required=True)
    parser.add_argument("--gsea-stats", type=Path, required=True)
    parser.add_argument("--correlation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    rank = pd.read_csv(args.gsea_rank).dropna(subset=["gene", "tstat_low_vs_high"]).copy()
    rank["gene_upper"] = rank["gene"].astype(str).str.upper()
    rank = rank.sort_values("tstat_low_vs_high", ascending=False).reset_index(drop=True)
    gsea = pd.read_csv(args.gsea_stats)
    corr = pd.read_csv(args.correlation, sep="\t")
    if set(DISPLAY_ORDER) - set(corr["signature"]):
        raise ValueError("Correlation source is missing one or more requested programs.")

    fig = plt.figure(figsize=(11.2, 4.15), facecolor="white")
    outer = fig.add_gridspec(1, 3, width_ratios=[1.08, 1.08, 1.30], wspace=0.42)
    add_gsea_panel(fig, outer[0], rank, gsea, "Proliferation_cell_cycle", "a")
    add_gsea_panel(fig, outer[1], rank, gsea, "Stemness_CSC", "b")
    add_heatmap_panel(fig, outer[2], corr)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.93, bottom=0.22)

    stem = args.output / "PDE2A_GSEA_Spearman_publication"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white")
    plt.close(fig)

    (args.output / "source_manifest.json").write_text(json.dumps({
        "backend": "Python/matplotlib/seaborn",
        "inputs": {
            "GSEA_ranked_genes": str(args.gsea_rank),
            "GSEA_statistics": str(args.gsea_stats),
            "correlation_source": str(args.correlation),
        },
        "no_reanalysis": True,
        "panel_labels": {"a": "Proliferation / cell cycle GSEA", "b": "Stemness GSEA", "c": "PDE2A correlation heatmap"},
    }, indent=2), encoding="utf-8")
    (args.output / "figure_caption.md").write_text(
        "**PDE2A-low HNSCC tumors are associated with proliferation/cell-cycle and stemness-program enrichment.** "
        "Panels a and b show pre-ranked GSEA of GSE178537 bulk tumors stratified at the median PDE2A expression level; "
        "positive rank statistics indicate higher expression in PDE2A-low samples. The purple curves show the running "
        "enrichment score, black ticks indicate gene-set hits, the rank-metric strip is colored from PDE2A-low to PDE2A-high, "
        "and NES/FDR q values are shown in each panel. Panel c shows sample-level Spearman correlations between PDE2A "
        "expression and selected expression programs. P values were adjusted within each subset across all 12 prespecified "
        "signatures using the Benjamini-Hochberg method. *FDR q<0.05, **q<0.01, ***q<0.001. The IT/VT-paired-patient "
        "subset includes all specimens from patients with both sample types and is not a paired-difference analysis. "
        "Associations do not establish causality.\n",
        encoding="utf-8",
    )
    print(stem)


if __name__ == "__main__":
    main()
