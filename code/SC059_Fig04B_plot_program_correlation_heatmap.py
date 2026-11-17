#!/usr/bin/env python3
"""Draw a selected-program PDE2A correlation heatmap with an added EMT score."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


BASE = Path("/data/liaolz/PDE2A")
SRC = BASE / "Fig3_v3_audit_20260628"
RLD = BASE / "data/GSE178537_bulk_data/GSE178537_HNSCC_Patient_rld.txt.gz"
OUT = BASE / "Fig3B_selected_programs_heatmap_20260812"
OUT.mkdir(parents=True, exist_ok=True)

EMT_GENES = [
    "VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "ITGA5", "ITGB1",
    "MMP2", "MMP9", "LAMC2", "COL1A1", "COL1A2", "SPARC", "TNC",
]

TESTED_SIGNATURES = [
    "stemness_basal",
    "stemness_CSC_surface",
    "proliferation",
    "leading_edge",
    "tumor_core",
    "MAPK_activity",
    "Rap1_axis",
    "cAMP_axis_no_PDE2A",
    "epithelial_purity_proxy",
    "stromal_purity_proxy",
    "immune_purity_proxy",
    "EMT_invasion",
]

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

SUBSETS = ["all_bulk_samples", "tumor_IT_VT_only", "paired_IT_VT_patients"]
COL_LABELS = {
    "all_bulk_samples": "All specimens\nn = 65",
    "tumor_IT_VT_only": "Tumor IT/VT\nn = 35",
    "paired_IT_VT_patients": "IT/VT-paired patients\nn = 51 samples",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def read_rld(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    gene_col = frame.columns[0]
    frame.index = frame[gene_col].astype(str).str.replace(r"^ENSG\d+_", "", regex=True)
    frame = frame.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce")
    return frame.groupby(frame.index).mean()


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.sub(frame.mean(axis=1), axis=0).div(
        frame.std(axis=1).replace(0, np.nan), axis=0
    )


def add_emt_score(analysis: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    rld = read_rld(RLD)
    matched = [gene for gene in EMT_GENES if gene in rld.index]
    if not matched:
        raise ValueError("No EMT/invasion genes were found in the RLD matrix.")
    emt = zscore_rows(rld).loc[matched].mean(axis=0).rename("EMT_invasion")
    analysis = analysis.merge(emt.rename_axis("sample").reset_index(), on="sample", how="left")
    return analysis, matched


def calculate_stats(analysis: pd.DataFrame) -> pd.DataFrame:
    paired_ids = (
        analysis.pivot_table(
            index="patient", columns="region_raw", values="sample", aggfunc="first"
        )
        .dropna(subset=["IT", "VT"])
        .index
    )
    subset_frames = {
        "all_bulk_samples": analysis,
        "tumor_IT_VT_only": analysis[analysis["region_raw"].isin(["IT", "VT"])],
        "paired_IT_VT_patients": analysis[analysis["patient"].isin(paired_ids)],
    }

    rows = []
    for subset_name, subset in subset_frames.items():
        for signature in TESTED_SIGNATURES:
            ok = np.isfinite(subset["PDE2A"]) & np.isfinite(subset[signature])
            rho, p_value = stats.spearmanr(
                subset.loc[ok, "PDE2A"], subset.loc[ok, signature]
            )
            rows.append({
                "subset": subset_name,
                "signature": signature,
                "n_samples": int(ok.sum()),
                "spearman_rho": float(rho),
                "p_value": float(p_value),
            })

    result = pd.DataFrame(rows)
    result["fdr_within_subset_12_signatures"] = result.groupby("subset")["p_value"].transform(
        lambda values: multipletests(values, method="fdr_bh")[1]
    )
    return result


def stars(q_value: float) -> str:
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return ""


def main() -> None:
    analysis = pd.read_csv(SRC / "Fig3B_PDE2A_signature_correlation_matrix_used.csv")
    analysis, emt_genes_used = add_emt_score(analysis)
    all_stats = calculate_stats(analysis)
    all_stats.to_csv(OUT / "PDE2A_all_12_signature_correlation_stats.tsv", sep="\t", index=False)

    selected = all_stats[all_stats["signature"].isin(DISPLAY_ORDER)].copy()
    selected["signature_label"] = selected["signature"].map(ROW_LABELS)
    selected["subset_label"] = selected["subset"].map(COL_LABELS).str.replace("\n", " ", regex=False)
    selected.to_csv(OUT / "SourceData_PDE2A_selected_program_correlations.tsv", sep="\t", index=False)

    rho = selected.pivot(index="signature", columns="subset", values="spearman_rho").loc[
        DISPLAY_ORDER, SUBSETS
    ]
    qvals = selected.pivot(
        index="signature", columns="subset", values="fdr_within_subset_12_signatures"
    ).loc[DISPLAY_ORDER, SUBSETS]
    annotations = rho.copy().astype(object)
    for row in DISPLAY_ORDER:
        for col in SUBSETS:
            annotations.loc[row, col] = f"{rho.loc[row, col]:.2f}{stars(qvals.loc[row, col])}"

    rho.index = [ROW_LABELS[name] for name in rho.index]
    rho.columns = [COL_LABELS[name] for name in rho.columns]
    annotations.index = rho.index
    annotations.columns = rho.columns

    fig, ax = plt.subplots(figsize=(7.6, 4.4), facecolor="white")
    heatmap = sns.heatmap(
        rho,
        cmap=sns.diverging_palette(245, 25, s=82, l=48, center="light", as_cmap=True),
        center=0,
        vmin=-1,
        vmax=1,
        annot=annotations,
        fmt="",
        annot_kws={"fontsize": 8.2},
        linewidths=0.9,
        linecolor="white",
        cbar_kws={"label": "Spearman's rho with PDE2A expression", "shrink": 0.84, "pad": 0.035},
        ax=ax,
    )

    ax.set_title("PDE2A correlations with selected expression programs", fontsize=11, weight="bold", pad=12)
    ax.set_xlabel("Analysis subset", fontsize=9, labelpad=9)
    ax.set_ylabel("Expression program", fontsize=9, labelpad=10)
    ax.tick_params(axis="x", labelrotation=0, labelsize=7.5, length=0, pad=6)
    ax.tick_params(axis="y", labelrotation=0, labelsize=8.2, length=0, pad=5)
    for label in ax.get_xticklabels():
        label.set_ha("center")

    colorbar = heatmap.collections[0].colorbar
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    colorbar.ax.tick_params(labelsize=7, width=0.6, length=3)
    colorbar.outline.set_linewidth(0.6)

    fig.text(
        0.32,
        0.040,
        "Values are sample-level Spearman correlations; * FDR q<0.05, ** q<0.01, *** q<0.001.",
        ha="left",
        va="bottom",
        fontsize=6.8,
        color="#444444",
    )
    fig.text(
        0.32,
        0.017,
        "FDR was recalculated within each subset across all 12 prespecified signatures.",
        ha="left",
        va="bottom",
        fontsize=6.7,
        color="#555555",
    )
    fig.subplots_adjust(left=0.32, right=0.90, top=0.86, bottom=0.22)

    stem = OUT / "PDE2A_selected_program_correlation_heatmap"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    audit = pd.DataFrame({
        "signature": ["EMT_invasion"],
        "genes_requested": [";".join(EMT_GENES)],
        "genes_used": [";".join(emt_genes_used)],
        "n_genes_used": [len(emt_genes_used)],
        "score_method": ["mean of per-gene z scores across GSE178537 samples"],
    })
    audit.to_csv(OUT / "EMT_invasion_signature_audit.tsv", sep="\t", index=False)

    caption = (
        "**Associations between PDE2A expression and selected expression programs in GSE178537.** "
        "Heatmap values are sample-level Spearman correlations. EMT/invasion activity was calculated as "
        "the mean of per-gene standardized expression for a project-defined 15-gene signature. P values "
        "were adjusted within each analysis subset across all 12 prespecified signatures using the "
        "Benjamini-Hochberg method. The IT/VT-paired-patient subset includes all specimens from patients "
        "with both sample types and is not a paired-difference analysis. Correlations do not establish causality."
    )
    (OUT / "figure_caption.md").write_text(caption + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
