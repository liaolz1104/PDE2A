"""Create three standalone manuscript PNGs from frozen PDE2A drug-screen outputs.

No model is refit here. TCGA-HNSC results are external transcriptional/
pharmacogenomic projections, not measured clinical drug responses.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "mechanism_drug_work/fda_repurposing_work"
TCGA_DIR = ROOT / "results/PDE2A_TCGA_drug_projection_20260807"
OUT_DIR = ROOT / "results/PDE2A_TCGA_HN6_drug_validation_figures_20260808"

TCGA_FILE = TCGA_DIR / "TCGA_HNSC_PDE2A_drug_projection_per_database.tsv"
MASTER_FILE = WORK / "cache/candidate_master_table.tsv"
PRED_FILE = WORK / "cache/all_sample_predictions.tsv"
MODEL_FILE = WORK / "cache/all_model_results.tsv"

PURPLE = "#7C3AED"
PURPLE_DARK = "#4C1D95"
PURPLE_PALE = "#F3E8FF"
BLUE = "#2563A6"
TEAL = "#16877B"
AMBER = "#D69000"
RED = "#C2415D"
INK = "#172033"
GREY = "#A8B1BF"
GREY_DARK = "#667085"
GRID = "#E5E7EB"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.labelsize": 9.5,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "savefig.dpi": 600,
    }
)


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / name, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure_tcga_cross_database(tcga: pd.DataFrame) -> None:
    # Some active ingredients have more than one GDSC assay record. Average
    # within database first so the cross-database statistic is not record-order dependent.
    wide = tcga.pivot_table(
        index="drug", columns="database", values="loss_score_spearman_rho", aggfunc="mean"
    ).dropna(subset=["GDSC1", "CTRP2"])
    strict = tcga.pivot_table(
        index="drug", columns="database", values="TCGA_support_strict", aggfunc="max"
    ).reindex(wide.index).fillna(False)
    concordant = strict.get("GDSC1", False).astype(bool) & strict.get("CTRP2", False).astype(bool)

    rho, pval = stats.spearmanr(wide["GDSC1"], wide["CTRP2"])
    fig, ax = plt.subplots(figsize=(6.5, 5.35))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axvspan(-0.65, 0, ymin=0, ymax=0.5, color=PURPLE_PALE, alpha=0.78, zorder=0)
    ax.axhline(0, color="#98A2B3", lw=0.8, ls=(0, (3, 3)))
    ax.axvline(0, color="#98A2B3", lw=0.8, ls=(0, (3, 3)))

    other = wide.loc[~concordant]
    hit = wide.loc[concordant]
    ax.scatter(other["GDSC1"], other["CTRP2"], s=45, c=GREY, alpha=0.75,
               edgecolors="white", linewidths=0.7, zorder=2)
    ax.scatter(hit["GDSC1"], hit["CTRP2"], s=72, c=AMBER, alpha=0.95,
               edgecolors="white", linewidths=0.9, zorder=4)

    offsets = {
        "afatinib": (8, 9), "belinostat": (8, 8), "bosutinib": (8, 6),
        "erlotinib": (8, -8), "gefitinib": (-10, -18), "gemcitabine": (-10, 6),
    }
    for drug, row in hit.iterrows():
        if drug == "etoposide":
            continue
        dx, dy = offsets.get(drug, (5, 5))
        ax.annotate(drug.capitalize(), (row["GDSC1"], row["CTRP2"]),
                    xytext=(dx, dy), textcoords="offset points",
                    ha="right" if dx < 0 else "left", va="center",
                    fontsize=7.4, color="#805A00", fontweight="semibold")

    etop = wide.loc["etoposide"]
    ax.scatter([etop["GDSC1"]], [etop["CTRP2"]], s=205, marker="*", c=PURPLE,
               edgecolors=PURPLE_DARK, linewidths=0.9, zorder=7)
    ax.annotate("Etoposide", (etop["GDSC1"], etop["CTRP2"]), xytext=(-10, 7),
                textcoords="offset points", ha="right", va="center", fontsize=9,
                color=PURPLE_DARK, fontweight="bold")

    ax.text(0.03, 0.05, "Concordant PDE2A-loss\nassociation", transform=ax.transAxes,
            color=PURPLE_DARK, fontsize=8.2, fontweight="bold", va="bottom")
    ax.text(0.97, 0.97, f"Cross-database Spearman r = {rho:.3f}\nP = {pval:.4f}  |  n = {len(wide)} drugs",
            transform=ax.transAxes, ha="right", va="top", color=INK, fontsize=8.3,
            bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": GRID, "lw": 0.8})

    ax.set_xlim(-0.62, 0.33)
    ax.set_ylim(-0.68, 0.33)
    ax.set_xlabel("TCGA association: GDSC1 projection\nSpearman r (PDE2A-loss score vs predicted response)", labelpad=8)
    ax.set_ylabel("TCGA association: CTRP2 projection\nSpearman r (PDE2A-loss score vs predicted response)", labelpad=8)
    ax.set_title("TCGA-HNSC supports a cross-database drug-response signal", loc="left", color=INK, pad=18)
    ax.text(0, 1.015, "Negative values indicate greater predicted sensitivity in PDE2A-loss-like tumors",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, color=GREY_DARK)
    ax.tick_params(length=3.2, color=GREY_DARK)
    ax.spines["left"].set_color(GREY_DARK)
    ax.spines["bottom"].set_color(GREY_DARK)
    ax.text(0, -0.235, "External pharmacogenomic projection; TCGA does not provide measured drug response.",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.2, color=GREY_DARK)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.22, top=0.86)
    save(fig, "Figure_TCGA_cross_database_drug_screen.png")


def fisher_ci(r: float, n: int = 520) -> tuple[float, float]:
    r = np.clip(r, -0.999999, 0.999999)
    z = np.arctanh(r)
    half = 1.96 / np.sqrt(n - 3)
    return tuple(np.tanh([z - half, z + half]))


def figure_pde2a_association(tcga: pd.DataFrame, master: pd.DataFrame) -> None:
    candidates = master.loc[master["final_tier"].isin(["Tier 1", "Tier 2"]), ["drug", "final_tier"]]
    data = tcga.merge(candidates, on="drug", how="inner", suffixes=("", "_master"))
    data = data.dropna(subset=["PDE2A_spearman_rho"]).copy()
    order = (data.groupby("drug")["PDE2A_spearman_rho"].max()
             .sort_values(ascending=True).index.tolist())
    ypos = {drug: i for i, drug in enumerate(order)}
    offsets = {"GDSC1": -0.18, "GDSC2": 0.0, "CTRP2": 0.18}
    colors = {"GDSC1": BLUE, "GDSC2": TEAL, "CTRP2": AMBER}

    fig, ax = plt.subplots(figsize=(7.2, 4.85))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    if "etoposide" in ypos:
        ax.axhspan(ypos["etoposide"] - 0.48, ypos["etoposide"] + 0.48,
                   color=PURPLE_PALE, alpha=0.75, zorder=0)
    ax.axvspan(0, 0.48, color="#F6F1FF", alpha=0.42, zorder=0)
    ax.axvline(0, color="#98A2B3", lw=0.9, ls=(0, (3, 3)))

    for row in data.itertuples():
        y = ypos[row.drug] + offsets.get(row.database, 0)
        lo, hi = fisher_ci(row.PDE2A_spearman_rho, int(row.n_TCGA))
        color = colors.get(row.database, GREY_DARK)
        ax.plot([lo, hi], [y, y], color=color, lw=1.25, alpha=0.88, zorder=2)
        significant = bool(row.PDE2A_spearman_FDR < 0.05)
        ax.scatter([row.PDE2A_spearman_rho], [y], s=58,
                   facecolors=color if significant else "white", edgecolors=color,
                   linewidths=1.15, zorder=4)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([d.capitalize() for d in order])
    for tick, drug in zip(ax.get_yticklabels(), order):
        if drug == "etoposide":
            tick.set_color(PURPLE_DARK)
            tick.set_fontweight("bold")
    ax.set_ylim(-0.65, len(order) - 0.35)
    ax.set_xlim(-0.38, 0.48)
    ax.set_xlabel("Spearman r (PDE2A expression vs predicted drug response)\nPositive r: lower PDE2A corresponds to greater predicted sensitivity", labelpad=8)
    ax.set_title("Selected drugs show distinct associations with PDE2A expression", loc="left", color=INK, pad=18)
    ax.text(0, 1.015, "TCGA-HNSC, n = 520 primary tumors; bars are approximate 95% Fisher-z intervals",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, color=GREY_DARK)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=c,
                      markeredgecolor=c, label=db) for db, c in colors.items()]
    handles += [Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=INK,
                       markeredgecolor=INK, label="FDR < 0.05"),
                Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor="white",
                       markeredgecolor=INK, label="FDR >= 0.05")]
    ax.legend(handles=handles, ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.31),
              handletextpad=0.35, columnspacing=1.05)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=3, color=GREY_DARK)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(GREY_DARK)
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.25, top=0.84)
    save(fig, "Figure_TCGA_PDE2A_drug_association.png")


def figure_hn6_predictions(master: pd.DataFrame, pred: pd.DataFrame, models: pd.DataFrame) -> None:
    candidates = master.loc[master["final_tier"].isin(["Tier 1", "Tier 2"]), ["drug", "final_tier"]].copy()
    candidates["drug"] = candidates["drug"].str.lower()
    passed = models.loc[models["model_pass"].eq(True), ["database", "drug_id", "drug"]].copy()
    passed["drug_id"] = passed["drug_id"].astype(str)
    passed["drug"] = passed["drug"].str.lower()
    pred = pred.copy()
    pred["drug_id"] = pred["drug_id"].astype(str)
    pred["drug"] = pred["drug"].str.lower()
    pred = pred.merge(passed.drop_duplicates(), on=["database", "drug_id", "drug"], how="inner")
    pred = pred.merge(candidates, on="drug", how="inner")

    pred["model_z"] = pred.groupby(["database", "drug", "drug_id"])["predicted_response"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
    )
    consensus = pred.groupby(["drug", "sample"], as_index=False)["model_z"].median()
    sample_order = [
        "HN6_NC", "HN6_NC_r1", "HN6_NC_r2",
        "HN6_KD", "HN6_KD_r1", "HN6_KD_r2",
        "HN6_Vector", "HN6_Vector_r1", "HN6_Vector_r2",
        "HN6_OE", "HN6_OE_r1", "HN6_OE_r2",
    ]
    drug_order = candidates.sort_values(["final_tier", "drug"])["drug"].tolist()
    matrix = consensus.pivot(index="drug", columns="sample", values="model_z").reindex(
        index=drug_order, columns=sample_order
    )
    if matrix.isna().any().any():
        missing = matrix.isna().stack()
        raise ValueError(f"Missing HN6 prediction cells: {missing[missing].index.tolist()}")

    cmap = LinearSegmentedColormap.from_list("sensitivity", [PURPLE_DARK, PURPLE, "#FAFAFC", "#F5C15D", RED])
    vmax = max(1.5, float(np.nanpercentile(np.abs(matrix.values), 97)))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(8.25, 4.75))
    fig.patch.set_facecolor("white")
    im = ax.imshow(matrix.values, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    labels = ["NC-1", "NC-2", "NC-3", "KD-1", "KD-2", "KD-3",
              "Vector-1", "Vector-2", "Vector-3", "OE-1", "OE-2", "OE-3"]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(drug_order)))
    ax.set_yticklabels([d.capitalize() for d in drug_order])
    for tick, drug in zip(ax.get_yticklabels(), drug_order):
        if drug == "etoposide":
            tick.set_color(PURPLE_DARK)
            tick.set_fontweight("bold")
    for x in [2.5, 5.5, 8.5]:
        ax.axvline(x, color="white", lw=2.2)
    for y in np.arange(len(drug_order) - 1) + 0.5:
        ax.axhline(y, color="white", lw=0.65, alpha=0.75)

    group_colors = [GREY_DARK, PURPLE, "#5B8DB8", AMBER]
    for i, (name, color) in enumerate(zip(["NC", "KD", "Vector", "OE"], group_colors)):
        start = i * 3 - 0.5
        ax.plot([start, start + 3], [-0.72, -0.72], color=color, lw=5,
                solid_capstyle="butt", clip_on=False)
        ax.text(start + 1.5, -1.05, name, ha="center", va="bottom", color=color,
                fontsize=8.2, fontweight="bold", clip_on=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.028)
    cbar.set_label("Within-drug standardized predicted response", rotation=90, labelpad=9)
    cbar.ax.text(0.5, -0.055, "more\nsensitive", transform=cbar.ax.transAxes,
                 ha="center", va="top", fontsize=7, color=PURPLE_DARK)
    cbar.ax.text(0.5, 1.03, "less\nsensitive", transform=cbar.ax.transAxes,
                 ha="center", va="bottom", fontsize=7, color=RED)
    cbar.outline.set_linewidth(0.6)

    ax.set_title("HN6 RNA-seq predicts PDE2A-state-dependent drug responses", loc="left", color=INK, pad=34)
    ax.text(0, 1.13, "Median standardized prediction across available quality-gated pharmacogenomic models",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, color=GREY_DARK)
    ax.set_xlabel("HN6 RNA-seq samples", labelpad=10)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(0, -0.24, "Lower standardized response indicates greater predicted sensitivity; values are comparable within each drug, not across drugs.",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.2, color=GREY_DARK)
    fig.subplots_adjust(left=0.18, right=0.89, bottom=0.25, top=0.80)
    save(fig, "Figure_HN6_predicted_drug_response.png")


def main() -> None:
    tcga = pd.read_csv(TCGA_FILE, sep="\t")
    master = pd.read_csv(MASTER_FILE, sep="\t")
    pred = pd.read_csv(PRED_FILE, sep="\t")
    models = pd.read_csv(MODEL_FILE, sep="\t")
    figure_tcga_cross_database(tcga)
    figure_pde2a_association(tcga, master)
    figure_hn6_predictions(master, pred, models)
    for path in sorted(OUT_DIR.glob("*.png")):
        print(path)


if __name__ == "__main__":
    main()
