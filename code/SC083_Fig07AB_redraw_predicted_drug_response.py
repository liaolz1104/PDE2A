"""Large-type redraw from the frozen HN6 prediction figure source tables."""
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "drug_figure_revision_20260908"
DRUGS = ["etoposide", "bleomycin", "bosutinib", "lenalidomide",
         "midostaurin", "pemetrexed", "pyrimethamine"]
BLUE, GREEN = "#2878A5", "#527D61"
CMAP = LinearSegmentedColormap.from_list("response", ["#316A91", "#FAFAFA", "#B45F48"])
plt.rcParams.update({"font.family": "Arial", "font.size": 14,
                     "axes.labelsize": 15, "xtick.labelsize": 14,
                     "ytick.labelsize": 14, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})

heat = {k: pd.read_csv(SOURCE / f"HN6_{k}_heatmap_source.tsv", sep="\t")
        for k in ("KD", "OE")}
forest = {k: pd.read_csv(SOURCE / f"HN6_{k}_model_effects_source.tsv", sep="\t")
          for k in ("KD", "OE")}
assert all(len(heat[k]) == 54 and len(forest[k]) == 9 for k in heat)
assert all(set(heat[k].drug) == set(DRUGS) for k in heat)
assert all((forest[k].low <= forest[k].effect).all() and
           (forest[k].effect <= forest[k].high).all() for k in forest)

def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)

def strip(ax, title, color):
    ax.add_patch(Rectangle((0, 1.02), 1, .14, transform=ax.transAxes,
                           color=color, clip_on=False))
    ax.text(.5, 1.09, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=16, color="white", weight="bold")

limit = max(float(np.abs(h.control_centered_response).max()) for h in heat.values())
fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.6))
fig.subplots_adjust(left=.15, right=.88, top=.79, bottom=.24, wspace=.48)
for ax, k, title, strip_color, groups in zip(
        axes, ("KD", "OE"), ("siCtrl vs siPDE2A", "oeCtrl vs oePDE2A"),
        (BLUE, GREEN), (("NC", "KD"), ("Vector", "OE"))):
    h = heat[k]
    assert h.groupby(["drug", "sample"]).size().max() <= 2
    samples = [f"HN6_{g}{suffix}" for g in groups for suffix in ("", "_r1", "_r2")]
    matrix = (h.groupby(["drug", "sample"]).control_centered_response.median()
              .unstack("sample").reindex(index=DRUGS, columns=samples))
    assert matrix.notna().all().all()
    image = ax.imshow(matrix, cmap=CMAP, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                      aspect="auto", interpolation="nearest")
    ax.set_yticks(range(7), [d.capitalize() for d in DRUGS])
    labels = ([f"{groups[0]}-{i}" for i in (1, 2, 3)] +
              [f"{groups[1]}-{i}" for i in (1, 2, 3)])
    ax.set_xticks(range(6), labels, rotation=45, ha="right")
    ax.axvline(2.5, color="white", lw=3)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    strip(ax, title, strip_color)
    matrix.rename_axis("drug").to_csv(HERE / f"heatmap_{k}_display_matrix.tsv", sep="\t")
cax = fig.add_axes([.905, .32, .018, .39])
cb = fig.colorbar(image, cax=cax)
cb.set_label("Control-centered predicted response\n(training-response SD units)", fontsize=13)
cb.ax.tick_params(labelsize=13)
save(fig, "HN6_predicted_response_heatmaps_large")

fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.8), sharey=True)
fig.subplots_adjust(left=.16, right=.98, top=.78, bottom=.20, wspace=.11)
for ax, k, title, strip_color in zip(
        axes, ("KD", "OE"), ("siCtrl vs siPDE2A", "oeCtrl vs oePDE2A"),
        (BLUE, GREEN)):
    f = forest[k]
    ax.set_ylim(6.5, -.5)
    ax.set_yticks(range(7), [d.capitalize() for d in DRUGS])
    ax.grid(axis="y", color="#E7E7E7", lw=1.1)
    ax.grid(axis="x", color="#F0F0F0", lw=.8)
    ax.set_axisbelow(True)
    ax.axvline(0, color="#777777", lw=1.6, ls=(0, (4, 3)), zorder=1)
    for db, color, offset in (("GDSC1", BLUE, -.13), ("CTRP2", GREEN, .13)):
        subset = f[f.database.eq(db)]
        yy = subset.drug.map({d: i for i, d in enumerate(DRUGS)}).to_numpy() + offset
        ax.hlines(yy, subset.low, subset.high, color=color, lw=2.4, zorder=3)
        ax.vlines(subset.low, yy-.065, yy+.065, color=color, lw=2.4, zorder=3)
        ax.vlines(subset.high, yy-.065, yy+.065, color=color, lw=2.4, zorder=3)
        ax.scatter(subset.effect, yy, marker="D", s=63, facecolor=color,
                   edgecolor="white", lw=.6, zorder=4)
    ax.tick_params(axis="both", length=4, width=1.2)
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"): ax.spines[side].set_linewidth(1.3)
    strip(ax, title, strip_color)
    ax.set_xlabel("Predicted response change / training SD")
    f[["database", "drug", "effect", "low", "high"]].to_csv(
        HERE / f"forest_{k}_display_data.tsv", sep="\t", index=False)
axes[0].set_xlim(-1.14, .30)
axes[1].set_xlim(-.36, .40)
axes[1].tick_params(labelleft=False)
fig.legend([Line2D([0], [0], color=c, marker="D", lw=2.4, markersize=7)
            for c in (BLUE, GREEN)], ["GDSC1", "CTRP2"],
           loc="upper center", bbox_to_anchor=(.57, .985), ncol=2,
           frameon=False, fontsize=14)
save(fig, "HN6_predicted_response_forest_large")

for k in ("KD", "OE"):
    for kind in ("heatmap", "model_effects"):
        shutil.copy2(SOURCE / f"HN6_{k}_{kind}_source.tsv", HERE)
print(HERE)
