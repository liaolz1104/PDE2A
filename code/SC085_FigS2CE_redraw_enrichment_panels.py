"""Arial 14 pt redraws from frozen ORA, GSEA and spatial annotation tables."""
from pathlib import Path
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
BULK = ROOT / "PDE2A_5patient_bulkRNAseq_pipeline/results/07_figures/manuscript_redraw_v4_Arial7pt"
plt.rcParams.update({"font.family":"Arial", "font.size":14, "axes.titlesize":14,
                     "axes.labelsize":14, "xtick.labelsize":14, "ytick.labelsize":14,
                     "legend.fontsize":14, "pdf.fonttype":42, "svg.fonttype":"none"})

def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)

def style(ax):
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", width=1.4, length=5)
    ax.grid(axis="x", color="#E5E5E5", lw=1.0)
    ax.set_axisbelow(True)

# Figure 1: two ORA bubble panels. Values and terms come from frozen source tables.
fig, axes = plt.subplots(1, 2, figsize=(18, 9.6))
fig.subplots_adjust(left=.18, right=.93, top=.93, bottom=.20, wspace=.73)
for ax, folder, title, cmap in zip(axes,
        ("Figure_10_ORA_EdgeDown", "Figure_11_ORA_EdgeUp"),
        ("Edge-downregulated DEGs", "Edge-upregulated DEGs"),
        (plt.cm.Blues, plt.cm.Reds)):
    d = pd.read_csv(BULK / folder / "source_data.tsv", sep="\t")
    y = np.arange(len(d))
    colors = cmap(.30 + .62 * np.clip(d.minus_log10_padj / 5, 0, 1))
    ax.hlines(y, 0, d.foldEnrichment, color=colors, lw=2.5)
    points = ax.scatter(d.foldEnrichment, y, c=d.minus_log10_padj,
                        cmap=cmap, vmin=0, vmax=5, s=30+d.overlapSize*26,
                        edgecolor="white", lw=.8, zorder=3)
    ax.set_yticks(y, [textwrap.fill(v, width=28) for v in d.pathway_display_name])
    ax.invert_yaxis(); ax.set_title(title, weight="bold", pad=13)
    ax.set_xlim(0, max(d.foldEnrichment)*1.10)
    ax.set_xlabel("Fold enrichment"); style(ax)
    cb = fig.colorbar(points, ax=ax, shrink=.56, pad=.025)
    cb.set_label("−log10(FDR q-value)", fontsize=14)
    cb.ax.tick_params(labelsize=14, width=1.2)
    sizes = [3, 5, 7, 9] if "Down" in folder else [1, 2, 3, 4, 5]
    handles = [Line2D([0], [0], marker="o", ls="", color="#777777",
                      markersize=np.sqrt(30+n*26), label=str(n)) for n in sizes]
    fig.legend(handles=handles, title="Overlapping genes", loc="lower center",
               bbox_to_anchor=(.30 if "Down" in folder else .78, .025),
               ncol=len(sizes), frameon=False, title_fontsize=14)
save(fig, "ORA_EdgeDown_EdgeUp_Arial14")

# Figure 2: two GSEA lollipop panels with the original NES and FDR values.
fig, axes = plt.subplots(1, 2, figsize=(17.5, 8.0))
fig.subplots_adjust(left=.19, right=.94, top=.90, bottom=.22, wspace=.70)
for ax, folder, title in zip(axes,
        ("Figure_08_GSEA_GO_BP", "Figure_09_GSEA_Reactome"),
        ("GO Biological Process GSEA", "Reactome pathway GSEA")):
    d = pd.read_csv(BULK / folder / "source_data.tsv", sep="\t")
    d = d.sort_values("NES", ascending=False).reset_index(drop=True)
    y = np.arange(len(d))
    col = np.where(d.NES > 0, "#E64B35", "#3C5488")
    for yy, row, c in zip(y, d.itertuples(), col):
        ax.plot([0, row.NES], [yy, yy], color=c, lw=3.0, zorder=2)
        ax.scatter(row.NES, yy, s=110+row.minus_log10_padj*3.4,
                   color=c, edgecolor="white", lw=.8, zorder=3)
    ax.axvline(0, color="#444444", lw=1.7)
    ax.set_yticks(y, [textwrap.fill(v, width=27) for v in d.pathway_display_name])
    ax.invert_yaxis(); ax.set_xlim(-3.2, 3.6)
    ax.set_title(title, weight="bold", pad=13)
    ax.set_xlabel("Normalized enrichment score (NES; Edge vs. Core)")
    style(ax)
fig.legend([Line2D([0],[0],color=c,lw=3,marker="o",markersize=9)
            for c in ("#3C5488", "#E64B35")],
           ["Core-enriched", "Edge-enriched"], loc="upper center",
           bbox_to_anchor=(.56, .995), ncol=2, frameon=False)
fig.legend(handles=[Line2D([0],[0],marker="o",ls="",color="#444444",
                   markersize=np.sqrt(110+n*3.4),label=str(n))
            for n in (10,20,30,40)],
           title="−log10(FDR q-value)", loc="lower center",
           bbox_to_anchor=(.57,.018), ncol=4, frameon=False, title_fontsize=14)
save(fig, "GSEA_GO_Reactome_Arial14")

# Figure 4: author region labels, retaining each original spot coordinate.
src = ROOT / "results/Spatial_source_annotation_recovery_20260730/spatial_source_annotation_merged.csv"
spots = pd.read_csv(src, usecols=["sample", "x_hires", "y_hires", "source_author_region"])
samples = [f"sample_{i}" for i in range(1,13)]
palette = {"core":"#4C78A8", "transitory":"#F2CF5B",
           "edge":"#E45756", "nc":"#BDBDBD"}
fig, axes = plt.subplots(2, 6, figsize=(24, 9), dpi=300)
fig.subplots_adjust(left=.025,right=.99,top=.92,bottom=.13,wspace=.10,hspace=.16)
for ax, sample in zip(axes.flat, samples):
    d = spots[spots["sample"].eq(sample)]
    for group, color in palette.items():
        z = d[d.source_author_region.eq(group)]
        ax.scatter(z.x_hires, z.y_hires, s=7, color=color, linewidths=0,
                   rasterized=True)
    ax.set_title(sample, fontsize=14, weight="bold", pad=6)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
handles=[Line2D([0],[0],marker="o",ls="",color=color,markersize=9,label=group.title())
         for group,color in palette.items()]
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,
           bbox_to_anchor=(.5,.035),fontsize=14)
save(fig,"Spatial_author_regions_Arial14")
print(OUT)
