#!/usr/bin/env python3
"""Redraw Fig. 3B in the scatter/regression template without changing statistics."""

from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT_MATRIX = ROOT / "results/Fig3_v3_audit_20260628/Fig3B_PDE2A_signature_correlation_matrix_used.csv"
RLD = ROOT / "data/GSE178537_bulk_data/GSE178537_HNSCC_Patient_rld.txt.gz"
FROZEN_STATS = ROOT / "results/Fig3B_selected_programs_heatmap_20260812/SourceData_PDE2A_selected_program_correlations.tsv"
OUT = ROOT / "results/Fig3B_template_scatter_Arial17_global_axes_20261106"

FONT_PT = 17
FONT_REGULAR = Path("C:/Windows/Fonts/arial.ttf")
FONT_BOLD = Path("C:/Windows/Fonts/arialbd.ttf")
REGULAR = FontProperties(fname=str(FONT_REGULAR), size=FONT_PT)
BOLD = FontProperties(fname=str(FONT_BOLD), size=FONT_PT)

EMT_GENES = [
    "VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "ITGA5", "ITGB1",
    "MMP2", "MMP9", "LAMC2", "COL1A1", "COL1A2", "SPARC", "TNC",
]
TESTED = [
    "stemness_basal", "stemness_CSC_surface", "proliferation", "leading_edge",
    "tumor_core", "MAPK_activity", "Rap1_axis", "cAMP_axis_no_PDE2A",
    "epithelial_purity_proxy", "stromal_purity_proxy", "immune_purity_proxy",
    "EMT_invasion",
]
PROGRAMS = [
    "stemness_basal", "stemness_CSC_surface", "proliferation", "cAMP_axis_no_PDE2A",
]
LABELS = {
    "stemness_basal": "Stemness\nbasal",
    "stemness_CSC_surface": "Stemness\nCSC",
    "proliferation": "Proliferation",
    "EMT_invasion": "EMT /\ninvasion",
    "leading_edge": "Leading-edge\nprogram",
    "tumor_core": "Tumor-core\nprogram",
    "cAMP_axis_no_PDE2A": "cAMP-related\nprogram",
    "MAPK_activity": "MAPK",
}
COLORS = {
    "stemness_basal": "#F06A5F",
    "stemness_CSC_surface": "#F29F80",
    "proliferation": "#4C9BD3",
    "EMT_invasion": "#14AFA4",
    "leading_edge": "#7B6DB0",
    "tumor_core": "#8A6F5A",
    "cAMP_axis_no_PDE2A": "#E3A51A",
    "MAPK_activity": "#6E9D64",
}
SUBSETS = ["all_bulk_samples", "tumor_IT_VT_only", "paired_IT_VT_patients"]
SUBSET_LABELS = {
    "all_bulk_samples": "All specimens\n(n = 65)",
    "tumor_IT_VT_only": "Tumor IT/VT\n(n = 35)",
    "paired_IT_VT_patients": "IT/VT-paired\npatients\n(n = 51 samples)",
}


def configure_style() -> None:
    if not FONT_REGULAR.exists() or not FONT_BOLD.exists():
        raise FileNotFoundError("Arial font files are required.")
    mpl.font_manager.fontManager.addfont(str(FONT_REGULAR))
    mpl.font_manager.fontManager.addfont(str(FONT_BOLD))
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": FONT_PT,
        "axes.titlesize": FONT_PT, "axes.labelsize": FONT_PT,
        "xtick.labelsize": FONT_PT, "ytick.labelsize": FONT_PT,
        "legend.fontsize": FONT_PT, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 300, "axes.unicode_minus": True,
    })


def read_rld() -> pd.DataFrame:
    frame = pd.read_csv(RLD, sep="\t")
    gene_col = frame.columns[0]
    frame.index = frame[gene_col].astype(str).str.replace(r"^ENSG\d+_", "", regex=True)
    frame = frame.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce")
    return frame.groupby(frame.index).mean()


def prepare_data() -> tuple[pd.DataFrame, list[str]]:
    analysis = pd.read_csv(INPUT_MATRIX)
    rld = read_rld()
    matched = [gene for gene in EMT_GENES if gene in rld.index]
    if not matched:
        raise ValueError("No EMT/invasion genes were found in the RLD matrix.")
    z = rld.sub(rld.mean(axis=1), axis=0).div(rld.std(axis=1).replace(0, np.nan), axis=0)
    emt = z.loc[matched].mean(axis=0).rename("EMT_invasion")
    analysis = analysis.merge(emt.rename_axis("sample").reset_index(), on="sample", how="left")
    if len(analysis) != 65 or analysis[PROGRAMS + ["PDE2A"]].isna().any().any():
        raise ValueError("Expected 65 complete samples after adding the EMT score.")
    return analysis, matched


def subset_frames(analysis: pd.DataFrame) -> dict[str, pd.DataFrame]:
    paired = (analysis.pivot_table(index="patient", columns="region_raw", values="sample", aggfunc="first")
              .dropna(subset=["IT", "VT"]).index)
    frames = {
        "all_bulk_samples": analysis,
        "tumor_IT_VT_only": analysis[analysis.region_raw.isin(["IT", "VT"])],
        "paired_IT_VT_patients": analysis[analysis.patient.isin(paired)],
    }
    expected = {"all_bulk_samples": 65, "tumor_IT_VT_only": 35, "paired_IT_VT_patients": 51}
    if {key: len(value) for key, value in frames.items()} != expected:
        raise ValueError("Analysis subset sizes do not match the frozen figure.")
    return frames


def calculate_stats(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for subset_name, subset in frames.items():
        for signature in TESTED:
            rho, p_value = stats.spearmanr(subset.PDE2A, subset[signature])
            rows.append({"subset": subset_name, "signature": signature,
                         "n_samples": len(subset), "spearman_rho": rho, "p_value": p_value})
    result = pd.DataFrame(rows)
    result["fdr_within_subset_12_signatures"] = result.groupby("subset")["p_value"].transform(
        lambda values: multipletests(values, method="fdr_bh")[1]
    )
    return result


def verify_frozen(stats_now: pd.DataFrame) -> pd.DataFrame:
    frozen = pd.read_csv(FROZEN_STATS, sep="\t")
    check = frozen.merge(stats_now, on=["subset", "signature"], suffixes=("_frozen", "_now"))
    check["rho_abs_diff"] = (check.spearman_rho_frozen - check.spearman_rho_now).abs()
    check["q_abs_diff"] = (check.fdr_within_subset_12_signatures_frozen -
                           check.fdr_within_subset_12_signatures_now).abs()
    check["match"] = (check.rho_abs_diff < 1e-12) & (check.q_abs_diff < 1e-12)
    if len(check) != 24 or not check.match.all():
        raise AssertionError("Recomputed selected-program statistics differ from the frozen heatmap.")
    return check


def significance(q_value: float) -> str:
    if q_value < .001: return "***"
    if q_value < .01: return "**"
    if q_value < .05: return "*"
    return "n.s."


def regression_band(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    order = np.argsort(x)
    x_sorted = x[order]
    slope, intercept, _, _, _ = stats.linregress(x, y)
    fitted = intercept + slope * x_sorted
    residual = y - (intercept + slope * x)
    dof = len(x) - 2
    s_err = np.sqrt(np.sum(residual ** 2) / dof)
    x_mean = x.mean()
    sxx = np.sum((x - x_mean) ** 2)
    ci = stats.t.ppf(.975, dof) * s_err * np.sqrt(1 / len(x) + (x_sorted - x_mean) ** 2 / sxx)
    ax.fill_between(x_sorted, fitted - ci, fitted + ci, color=color, alpha=.16, linewidth=0)
    ax.plot(x_sorted, fitted, color=color, lw=2.1)


def kde_fill(ax: plt.Axes, values: np.ndarray, color: str, vertical: bool = False,
             alpha: float = .28) -> None:
    grid = np.linspace(values.min(), values.max(), 250)
    density = stats.gaussian_kde(values)(grid)
    if vertical:
        ax.fill_betweenx(grid, 0, density, color=color, alpha=alpha)
        ax.plot(density, grid, color=color, lw=1.1)
    else:
        ax.fill_between(grid, 0, density, color=color, alpha=alpha)
        ax.plot(grid, density, color=color, lw=1.1)


def draw(analysis: pd.DataFrame, frames: dict[str, pd.DataFrame], results: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(15.5, 19.5), facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[5.8, 10.2], left=.19, right=.975,
                             bottom=.18, top=.95, hspace=.25)

    top_band = outer[0].subgridspec(1, 3, width_ratios=[.4, 3.2, .4], wspace=0)
    top = top_band[0, 1].subgridspec(2, 2, height_ratios=[1.0, 5.0],
                                     width_ratios=[6.5, .8], hspace=.03, wspace=.04)
    ax_x = fig.add_subplot(top[0, 0])
    ax_main = fig.add_subplot(top[1, 0], sharex=ax_x)
    ax_y = fig.add_subplot(top[1, 1], sharey=ax_main)
    representatives = ["stemness_basal", "proliferation", "cAMP_axis_no_PDE2A"]
    for program in representatives:
        color = COLORS[program]
        ax_main.scatter(analysis.PDE2A, analysis[program], s=28, color=color, alpha=.42,
                        edgecolors="none", label=LABELS[program].replace("\n", " "))
        kde_fill(ax_y, analysis[program].to_numpy(), color, vertical=True, alpha=.22)
    kde_fill(ax_x, analysis.PDE2A.to_numpy(), "#F2B705", alpha=.70)
    ax_main.set_xlabel("PDE2A expression", fontproperties=BOLD)
    ax_main.set_ylabel("Program score", fontproperties=BOLD)
    ax_main.grid(True, color="#E7E7E7", lw=.8)
    ax_main.set_axisbelow(True)
    ax_main.legend(loc="upper left", ncol=3, frameon=True, prop=REGULAR,
                   facecolor="white", edgecolor="none", framealpha=.88,
                   columnspacing=.6, handletextpad=.25, borderpad=.25)
    ax_x.axis("off"); ax_y.axis("off")
    for spine in ax_main.spines.values():
        spine.set_color("#555555"); spine.set_linewidth(1.0)

    grid = outer[1].subgridspec(3, len(PROGRAMS), hspace=.27, wspace=.14)
    stat_lookup = results.set_index(["subset", "signature"])
    x_margin = (analysis.PDE2A.max() - analysis.PDE2A.min()) * .04
    x_limits = (analysis.PDE2A.min() - x_margin, analysis.PDE2A.max() + x_margin)
    y_limits = {}
    for program in PROGRAMS:
        margin = (analysis[program].max() - analysis[program].min()) * .07
        y_limits[program] = (analysis[program].min() - margin, analysis[program].max() + margin)
    for row, subset_name in enumerate(SUBSETS):
        subset = frames[subset_name]
        for col, program in enumerate(PROGRAMS):
            ax = fig.add_subplot(grid[row, col])
            color = COLORS[program]
            x = subset.PDE2A.to_numpy(float)
            y = subset[program].to_numpy(float)
            ax.scatter(x, y, s=18, color=color, alpha=.43, edgecolors="none")
            regression_band(ax, x, y, color)
            ax.set_xlim(*x_limits)
            ax.set_ylim(*y_limits[program])
            stat_row = stat_lookup.loc[(subset_name, program)]
            marker = significance(float(stat_row.fdr_within_subset_12_signatures))
            ax.text(.04, .94, f"ρ = {stat_row.spearman_rho:.2f}", transform=ax.transAxes,
                    ha="left", va="top", fontproperties=BOLD)
            ax.text(.04, .80, marker, transform=ax.transAxes, ha="left", va="top",
                    fontproperties=BOLD)
            ax.grid(True, color="#E8E8E8", lw=.65)
            ax.set_axisbelow(True)
            ax.tick_params(length=2.5, width=.8, labelsize=FONT_PT)
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontproperties(REGULAR)
            for spine in ax.spines.values():
                spine.set_color("#5C5C5C"); spine.set_linewidth(.9)
            if row < 2:
                ax.set_xticklabels([])
            if col > 0:
                ax.set_yticklabels([])
            if row == 0:
                ax.add_patch(Rectangle((0, 1.02), 1, .27, transform=ax.transAxes,
                                       color=color, clip_on=False))
                ax.text(.5, 1.155, LABELS[program], transform=ax.transAxes,
                        ha="center", va="center", color="white", fontproperties=BOLD,
                        linespacing=1.18)
            if col == 0:
                ax.text(-.27, .5, SUBSET_LABELS[subset_name], transform=ax.transAxes,
                        ha="right", va="center", rotation=0, fontproperties=BOLD,
                        linespacing=1.05)

    fig.add_artist(FancyArrowPatch((.183, .18), (.183, .59), transform=fig.transFigure,
                                  arrowstyle="-|>", mutation_scale=22, lw=2.0,
                                  color="black", clip_on=False))
    fig.text(.164, .385, "Feature score", ha="center", va="center", rotation=90,
             fontproperties=BOLD)
    fig.add_artist(FancyArrowPatch((.19, .142), (.975, .142), transform=fig.transFigure,
                                  arrowstyle="-|>", mutation_scale=22, lw=2.0,
                                  color="black", clip_on=False))
    fig.text(.582, .108, "PDE2A expression", ha="center", va="center",
             fontproperties=BOLD)

    fig.suptitle("PDE2A-associated expression programs", x=.50, y=.984,
                 fontproperties=BOLD)
    fig.text(.018, .966, "B", ha="left", va="top", fontproperties=BOLD)
    fig.text(.50, .045,
             "Spearman ρ; BH-FDR across 12 prespecified signatures within each subset.\n"
             "Regression lines and 95% confidence intervals are descriptive; "
             "* q < 0.05, ** q < 0.01, *** q < 0.001.",
             ha="center", va="center", fontproperties=REGULAR)
    return fig


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    analysis, matched = prepare_data()
    frames = subset_frames(analysis)
    results = calculate_stats(frames)
    verification = verify_frozen(results)

    point_rows = []
    for subset_name, subset in frames.items():
        for program in PROGRAMS:
            block = subset[["sample", "patient", "region_raw", "PDE2A", program]].copy()
            block.insert(0, "subset", subset_name)
            block.insert(1, "program", program)
            block = block.rename(columns={program: "program_score"})
            point_rows.append(block)
    pd.concat(point_rows, ignore_index=True).to_csv(OUT / "SourceData_scatter_points.tsv", sep="\t", index=False)
    results.to_csv(OUT / "recomputed_all_12_signature_stats.tsv", sep="\t", index=False)
    verification.to_csv(OUT / "FROZEN_STATISTICS_VERIFICATION.tsv", sep="\t", index=False)
    pd.DataFrame({
        "font_family": ["Arial"], "font_pt": [FONT_PT],
        "statistics_changed": [False], "n_programs": [len(PROGRAMS)],
        "n_subsets": [3], "n_scatter_facets": [3 * len(PROGRAMS)],
    }).to_csv(OUT / "STYLE_AUDIT.tsv", sep="\t", index=False)
    pd.DataFrame({"requested_gene": EMT_GENES,
                  "used": [gene in matched for gene in EMT_GENES]}).to_csv(
        OUT / "EMT_invasion_signature_audit.tsv", sep="\t", index=False)

    fig = draw(analysis, frames, results)
    stem = OUT / "PDE2A_program_correlations_scatter_template_Arial17_global_axes"
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    shutil.copy2(INPUT_MATRIX, OUT / INPUT_MATRIX.name)
    shutil.copy2(FROZEN_STATS, OUT / FROZEN_STATS.name)
    shutil.copy2(__file__, OUT / "plotting_script.py")
    (OUT / "README.txt").write_text(
        "Template redraw of the frozen Fig. 3B statistics. No correlation or FDR value was changed.\n"
        "The top overlay is descriptive and shows three representative programs for readability.\n"
        "The 12 lower facets show the 4 requested programs across all 3 original analysis subsets.\n"
        "The paired-patient subset contains all specimens from patients with both IT and VT; it is not a paired-difference analysis.\n",
        encoding="utf-8",
    )
    print(f"Created: {stem}")


if __name__ == "__main__":
    main()
