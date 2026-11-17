#!/usr/bin/env python3
"""
Figure 4: PDE2A restricts malignant progression in HNSCC.

This script reruns the Figure 4 evidence chain on the v3 GSE198315
annotation and GSE178537 paired bulk RNA-seq:
  A. PDE2A expression and PDE2A-high/low malignant cells on UMAP.
  B. PDE2A-low versus PDE2A-high malignant-cell differential volcano.
  C. Enrichment of PDE2A-loss activated genes.
  D. Bulk PDE2A-low versus PDE2A-high preranked GSEA.
  E. Bulk immune/TIDE-like signature scores by PDE2A group.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import anndata as ad
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

H5AD = Path("/data/liaolz/PDE2A/data/GSE198315_data/GSE198315_full_annotated_v3.h5ad")
RLD_FILE = Path("/data/liaolz/PDE2A/data/GSE178537_bulk_data/GSE178537_HNSCC_Patient_rld.txt.gz")
LABEL_FILE = Path("/data/liaolz/PDE2A/Figure3_output/bulk_sample_labels.csv")
OUT_DIR = Path("/data/liaolz/PDE2A/Figure4_v3_output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", font="DejaVu Sans")
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

GENE_SETS = {
    "cAMP_Rap1_signaling": [
        "ADCY1", "ADCY2", "ADCY3", "ADCY5", "ADCY6", "ADCY7", "ADCY9",
        "PRKACA", "PRKACB", "PRKAR1A", "PRKAR2A", "RAPGEF3", "RAPGEF4",
        "RAP1A", "RAP1B", "BRAF", "MAPK1", "MAPK3", "CREB1", "CREB3L1",
    ],
    "MAPK_ERK_activation": [
        "EGFR", "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1",
        "MAP2K2", "MAPK1", "MAPK3", "DUSP4", "DUSP5", "FOS", "JUN", "ETV4", "ETV5",
    ],
    "Stemness_CSC": [
        "SOX2", "KLF4", "ALDH1A1", "PROM1", "CD44", "ITGA6", "LGR5",
        "EPCAM", "MKI67", "TOP2A", "BIRC5", "HMGA2", "KRT5", "KRT14", "TP63",
    ],
    "EMT_invasion": [
        "VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "ITGA5", "ITGB1",
        "MMP2", "MMP9", "LAMC2", "COL1A1", "COL1A2", "SPARC", "TNC",
    ],
    "Proliferation_cell_cycle": [
        "MKI67", "TOP2A", "PCNA", "MCM2", "MCM4", "MCM6", "CDK1",
        "CCNB1", "CCNB2", "BUB1", "AURKB", "TYMS", "UBE2C",
    ],
    "Apoptosis": [
        "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "BCL2", "BCL2L1",
        "BAD", "BIRC5", "XIAP", "FAS", "TNFSF10",
    ],
}

IMMUNE_SIGNATURES = {
    "CD8_T_cell": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "CXCL9"],
    "CD4_T_cell": ["CD3D", "CD3E", "IL7R", "CCR7", "LTB"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Macrophage_M1": ["IL1B", "TNF", "CXCL9", "CXCL10", "CD80"],
    "Macrophage_M2": ["CD163", "MRC1", "MSR1", "CCL18", "MERTK"],
    "Myeloid": ["LST1", "S100A8", "S100A9", "FCGR3A", "LYZ"],
    "Dendritic": ["ITGAX", "CLEC9A", "FCER1A", "LILRA4", "IRF8"],
    "TIDE_CTL": ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9"],
    "TIDE_Exclusion": ["TGFB1", "COL1A1", "COL1A2", "ACTA2", "VIM"],
    "TIDE_Dysfunction": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX"],
}


def savefig(fig: plt.Figure, name: str, w: float | None = None, h: float | None = None) -> None:
    if w and h:
        fig.set_size_inches(w, h)
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def expr_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def read_bulk() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(RLD_FILE, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    symbols = expr.index.to_series().str.replace(r"^ENSG[0-9.]+_", "", regex=True)
    expr.index = symbols.astype(str)
    expr = expr.groupby(expr.index).mean(numeric_only=True)
    labels = pd.read_csv(LABEL_FILE).rename(columns={"Unnamed: 0": "sample"})
    labels["sample"] = labels["sample"].astype(str)
    labels = labels[labels["sample"].isin(expr.columns)].copy()
    return expr, labels


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    z = mat.sub(mat.mean(axis=1), axis=0)
    z = z.div(mat.std(axis=1).replace(0, np.nan), axis=0)
    return z


def score_signatures(expr: pd.DataFrame, signatures: dict[str, list[str]]) -> pd.DataFrame:
    z = zscore_rows(expr)
    rows = []
    upper = pd.Series(expr.index, index=expr.index).str.upper()
    for sig, genes in signatures.items():
        wanted = {g.upper() for g in genes}
        matched = [idx for idx in expr.index if str(idx).upper() in wanted]
        if not matched:
            continue
        score = z.loc[matched].mean(axis=0)
        rows.append(pd.DataFrame({"sample": score.index, "signature": sig, "score": score.values, "n_genes": len(matched)}))
    return pd.concat(rows, ignore_index=True)


def run_sc_malignant() -> pd.DataFrame:
    print("[1/6] Load v3 malignant cells")
    backed = sc.read_h5ad(H5AD, backed="r")
    obs = backed.obs.copy()
    mask = obs["cell_type"].astype(str).isin(["Malignant", "Cycling_Malignant"])
    mask &= obs["region"].astype(str).isin(["Edge", "Core"])
    mal = backed[mask.values, :].to_memory()
    backed.file.close()

    pde = expr_vector(mal, "PDE2A")
    positive = pde > 0
    group = np.repeat("PDE2A_low", len(pde)).astype(object)
    group[positive] = "PDE2A_high"
    if positive.sum() < 20:
        cutoff = float(np.quantile(pde, 0.90))
        group = np.repeat("PDE2A_low", len(pde)).astype(object)
        group[pde > cutoff] = "PDE2A_high"
    print(
        "    PDE2A grouping:",
        f"low={int((group == 'PDE2A_low').sum())}",
        f"high={int((group == 'PDE2A_high').sum())}",
    )
    mal.obs["PDE2A_expr"] = pde
    mal.obs["PDE2A_group"] = pd.Categorical(group, categories=["PDE2A_low", "PDE2A_high"])

    meta = mal.obs[["patient", "region", "sample", "cell_type", "PDE2A_expr", "PDE2A_group"]].copy()
    meta.to_csv(OUT_DIR / "Fig4_malignant_pde2a_cell_metadata.csv")

    print("[2/6] Plot UMAP PDE2A state")
    xy = pd.DataFrame(mal.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=mal.obs_names)
    xy = pd.concat([xy, meta], axis=1)
    rng = np.random.default_rng(7)
    xy_plot = xy.sample(min(22000, xy.shape[0]), random_state=7)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4))
    sc0 = axes[0].scatter(
        xy_plot["UMAP_1"], xy_plot["UMAP_2"], c=xy_plot["PDE2A_expr"],
        s=4, cmap="viridis", linewidths=0, rasterized=True,
    )
    axes[0].set_title("PDE2A expression")
    plt.colorbar(sc0, ax=axes[0], fraction=0.04, pad=0.02).set_label("log-normalized expression")
    palette = {"PDE2A_low": "#C44E52", "PDE2A_mid": "#BDBDBD", "PDE2A_high": "#4C72B0"}
    for name, sub in xy_plot.groupby("PDE2A_group", observed=True):
        display_name = {"PDE2A_low": "PDE2A negative", "PDE2A_high": "PDE2A detected"}.get(str(name), str(name))
        axes[1].scatter(sub["UMAP_1"], sub["UMAP_2"], s=4, label=display_name,
                        color=palette.get(str(name), "grey"), linewidths=0, rasterized=True)
    axes[1].set_title("PDE2A-negative/detected malignant states")
    axes[1].legend(frameon=False, markerscale=3, loc="best")
    for ax in axes:
        ax.set_xlabel("UMAP_1")
        ax.set_ylabel("UMAP_2")
        ax.grid(False)
    savefig(fig, "Fig4A_PDE2A_malignant_UMAP")

    print("[3/6] Differential analysis PDE2A-low vs PDE2A-high")
    keep = mal.obs["PDE2A_group"].isin(["PDE2A_low", "PDE2A_high"]).values
    de = mal[keep, :].copy()
    de.obs["PDE2A_group"] = de.obs["PDE2A_group"].astype(str).astype("category")
    sc.tl.rank_genes_groups(de, "PDE2A_group", groups=["PDE2A_low"], reference="PDE2A_high", method="wilcoxon")
    rg = de.uns["rank_genes_groups"]
    deg = pd.DataFrame({
        "gene": rg["names"]["PDE2A_low"],
        "score": rg["scores"]["PDE2A_low"],
        "logfoldchange_low_vs_high": rg["logfoldchanges"]["PDE2A_low"],
        "pval": rg["pvals"]["PDE2A_low"],
        "padj": rg["pvals_adj"]["PDE2A_low"],
    })
    deg["neglog10_padj"] = -np.log10(deg["padj"].replace(0, np.nan))
    deg["vkd_effect"] = "not_affected"
    deg.loc[(deg["padj"] < 0.05) & (deg["logfoldchange_low_vs_high"] > 0.25), "vkd_effect"] = "PDE2A_loss_activated"
    deg.loc[(deg["padj"] < 0.05) & (deg["logfoldchange_low_vs_high"] < -0.25), "vkd_effect"] = "PDE2A_loss_repressed"
    deg.to_csv(OUT_DIR / "Fig4_scRNA_PDE2A_low_vs_high_DEG.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.7, 4.8))
    pal = {"not_affected": "#B8B8B8", "PDE2A_loss_activated": "#C44E52", "PDE2A_loss_repressed": "#4C72B0"}
    for eff, sub in deg.groupby("vkd_effect"):
        ax.scatter(sub["logfoldchange_low_vs_high"], sub["neglog10_padj"], s=8,
                   color=pal[eff], alpha=0.75, linewidths=0, label=eff.replace("_", " "))
    label_genes = ["PDE2A", "MKI67", "TOP2A", "CD44", "ITGA5", "VIM", "LAMC2", "EGFR", "MAPK1", "RAP1A"]
    lab = deg[deg["gene"].isin(label_genes)].copy()
    for _, row in lab.iterrows():
        if np.isfinite(row["neglog10_padj"]):
            ax.text(row["logfoldchange_low_vs_high"], row["neglog10_padj"], row["gene"], fontsize=7)
    ax.axvline(0, color="black", lw=0.7)
    ax.axvline(0.25, color="#777777", lw=0.6, ls="--")
    ax.axvline(-0.25, color="#777777", lw=0.6, ls="--")
    ax.axhline(-np.log10(0.05), color="#777777", lw=0.6, ls="--")
    ax.set_xlabel("log fold-change: PDE2A-negative vs PDE2A-detected")
    ax.set_ylabel("-log10 adjusted p value")
    ax.set_title("Virtual PDE2A-loss malignant-cell effects")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    savefig(fig, "Fig4B_VKD_PDE2A_volcano")
    return deg


def enrich_de_genes(deg: pd.DataFrame) -> pd.DataFrame:
    print("[4/6] Enrichment of PDE2A-loss activated genes")
    genes = deg.loc[deg["vkd_effect"] == "PDE2A_loss_activated", "gene"].head(500).astype(str).tolist()
    enrich_tables = []
    if len(genes) >= 10:
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
                organism="human",
                outdir=None,
                cutoff=0.5,
            )
            enrich_tables.append(enr.results)
        except Exception as exc:
            print(f"    Enrichr unavailable, fallback to custom Fisher tests: {exc}")

    if enrich_tables:
        enr = pd.concat(enrich_tables, ignore_index=True)
        enr["source"] = enr["Gene_set"]
        enr = enr.rename(columns={"Adjusted P-value": "adj_p", "P-value": "pval", "Overlap": "overlap"})
        enr["score"] = -np.log10(enr["adj_p"].replace(0, np.nan))
    else:
        universe = set(deg["gene"].astype(str))
        activated = set(genes)
        rows = []
        for term, gset in GENE_SETS.items():
            gs = set(gset) & universe
            a = len(gs & activated)
            b = len(activated - gs)
            c = len(gs - activated)
            d = len(universe - activated - gs)
            if len(gs) == 0:
                continue
            _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({"source": "custom", "Term": term, "pval": p, "overlap": f"{a}/{len(gs)}", "Genes": ";".join(sorted(gs & activated))})
        enr = pd.DataFrame(rows)
        if not enr.empty:
            enr["adj_p"] = multipletests(enr["pval"], method="fdr_bh")[1]
            enr["score"] = -np.log10(enr["adj_p"].replace(0, np.nan))

    enr = enr.sort_values(["adj_p", "pval"]).reset_index(drop=True)
    enr.to_csv(OUT_DIR / "Fig4C_PDE2A_loss_activated_enrichment.csv", index=False)

    top = enr.head(12).copy()
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    sns.barplot(data=top, y="Term", x="score", hue="source", dodge=False, ax=ax, palette="mako")
    ax.set_xlabel("-log10 adjusted p value")
    ax.set_ylabel("")
    ax.set_title("Programs enriched in PDE2A-loss activated genes")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    savefig(fig, "Fig4C_PDE2A_loss_enrichment")
    return enr


def bulk_gsea_and_immune() -> tuple[pd.DataFrame, pd.DataFrame]:
    print("[5/6] Bulk PDE2A high/low GSEA")
    expr, labels = read_bulk()
    tumor = labels[labels["Location"].isin(["Edge", "Core", "mLN"])].copy()
    tumor = tumor[tumor["sample"].isin(expr.columns)].copy()
    pde = expr.loc["PDE2A", tumor["sample"]]
    med = float(pde.median())
    tumor["PDE2A_group"] = np.where(pde.loc[tumor["sample"]].values >= med, "PDE2A_high", "PDE2A_low")
    tumor.to_csv(OUT_DIR / "Fig4_bulk_PDE2A_group_labels.csv", index=False)

    high = tumor.loc[tumor["PDE2A_group"] == "PDE2A_high", "sample"].tolist()
    low = tumor.loc[tumor["PDE2A_group"] == "PDE2A_low", "sample"].tolist()
    ranks = []
    for gene, row in expr.iterrows():
        a = row[low].astype(float).values
        b = row[high].astype(float).values
        try:
            t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        except Exception:
            t, p = np.nan, np.nan
        ranks.append((gene, float(np.nanmean(a) - np.nanmean(b)), t, p))
    rank_df = pd.DataFrame(ranks, columns=["gene", "mean_low_minus_high", "tstat_low_vs_high", "pval"])
    rank_df = rank_df.dropna(subset=["tstat_low_vs_high"]).sort_values("tstat_low_vs_high", ascending=False)
    rank_df.to_csv(OUT_DIR / "Fig4_bulk_PDE2A_low_vs_high_ranked_genes.csv", index=False)
    rnk = rank_df[["gene", "tstat_low_vs_high"]]
    rnk.to_csv(OUT_DIR / "Fig4_bulk_PDE2A_low_vs_high_prerank.rnk", sep="\t", header=False, index=False)

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=GENE_SETS,
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=7,
        threads=4,
        outdir=None,
        verbose=False,
    )
    gsea = pre.res2d.copy()
    gsea.to_csv(OUT_DIR / "Fig4D_bulk_PDE2A_low_vs_high_GSEA.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    gsea_plot = gsea.sort_values("NES")
    colors = ["#4C72B0" if x < 0 else "#C44E52" for x in gsea_plot["NES"]]
    ax.barh(gsea_plot["Term"], gsea_plot["NES"], color=colors)
    for y, (_, row) in enumerate(gsea_plot.iterrows()):
        ax.text(row["NES"], y, f"  FDR={row['FDR q-val']:.3g}", va="center", fontsize=7)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("NES: PDE2A-low vs PDE2A-high")
    ax.set_ylabel("")
    ax.set_title("Bulk PDE2A-low high preranked GSEA")
    savefig(fig, "Fig4D_bulk_PDE2A_low_high_GSEA")

    print("[6/6] Bulk immune/TIDE-like scoring")
    scores = score_signatures(expr[tumor["sample"]], IMMUNE_SIGNATURES)
    scores = scores.merge(tumor[["sample", "Location", "PDE2A_group"]], on="sample", how="left")
    scores.to_csv(OUT_DIR / "Fig4E_bulk_immune_TIDE_scores.csv", index=False)

    stats_rows = []
    for sig, sub in scores.groupby("signature"):
        a = sub.loc[sub["PDE2A_group"] == "PDE2A_low", "score"]
        b = sub.loc[sub["PDE2A_group"] == "PDE2A_high", "score"]
        stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        stats_rows.append({
            "signature": sig,
            "mean_low": a.mean(),
            "mean_high": b.mean(),
            "delta_low_minus_high": a.mean() - b.mean(),
            "pval": p,
        })
    score_stats = pd.DataFrame(stats_rows)
    score_stats["padj"] = multipletests(score_stats["pval"], method="fdr_bh")[1]
    score_stats.to_csv(OUT_DIR / "Fig4E_bulk_immune_TIDE_stats.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), gridspec_kw={"width_ratios": [1.2, 1.0]})
    order = score_stats.sort_values("delta_low_minus_high", ascending=False)["signature"].tolist()
    sns.boxplot(data=scores, x="score", y="signature", hue="PDE2A_group", order=order,
                palette={"PDE2A_low": "#C44E52", "PDE2A_high": "#4C72B0"}, ax=axes[0], fliersize=0)
    sns.stripplot(data=scores, x="score", y="signature", hue="PDE2A_group", order=order,
                  dodge=True, palette={"PDE2A_low": "#C44E52", "PDE2A_high": "#4C72B0"},
                  ax=axes[0], size=2, alpha=0.55, legend=False)
    axes[0].set_xlabel("signature z-score")
    axes[0].set_ylabel("")
    axes[0].set_title("Immune/TIDE-like scores")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles[:2], labels[:2], frameon=False, fontsize=7)

    heat = score_stats.set_index("signature").loc[order, ["delta_low_minus_high"]]
    sns.heatmap(heat, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={"label": "Low - High"},
                ax=axes[1])
    axes[1].set_ylabel("")
    axes[1].set_title("PDE2A-low immune shift")
    savefig(fig, "Fig4E_bulk_immune_TIDE_scores")
    return gsea, score_stats


def write_report(deg: pd.DataFrame, enr: pd.DataFrame, gsea: pd.DataFrame, immune_stats: pd.DataFrame) -> None:
    activated = int((deg["vkd_effect"] == "PDE2A_loss_activated").sum())
    repressed = int((deg["vkd_effect"] == "PDE2A_loss_repressed").sum())
    top_terms = "; ".join(enr.head(5)["Term"].astype(str).tolist()) if not enr.empty else "none"
    top_gsea = gsea.sort_values("NES", ascending=False).head(3)
    top_immune = immune_stats.sort_values("delta_low_minus_high", ascending=False).head(3)

    report = f"""# Figure 4 v3 result interpretation

Input scRNA object: `{H5AD}`
Input bulk matrix: `{RLD_FILE}`
Output directory: `{OUT_DIR}`

## Main findings

1. v3 malignant cells were split into PDE2A-negative and PDE2A-positive/detected states because PDE2A expression is highly sparse at the single-cell level.
2. PDE2A-low versus PDE2A-high malignant-cell comparison identified {activated} activated genes and {repressed} repressed genes at FDR < 0.05 and |logFC| > 0.25.
3. The strongest enrichment terms among PDE2A-loss activated genes were: {top_terms}.
4. Bulk GSE178537 PDE2A-low versus PDE2A-high GSEA showed the strongest positive NES terms:
{top_gsea[['Term', 'NES', 'FDR q-val']].to_string(index=False)}

5. Bulk immune/TIDE-like signatures with the largest PDE2A-low minus high shifts were:
{top_immune[['signature', 'delta_low_minus_high', 'pval', 'padj']].to_string(index=False)}

## Interpretation

The single-cell PDE2A-negative versus detected comparison is limited by sparse PDE2A
expression and should be interpreted as a virtual-loss proxy rather than a balanced
knockout experiment. In this run, the single-cell PDE2A-loss activated genes were
dominated by translation/ribosome terms, while the bulk PDE2A-low tumors showed
clearer proliferation and stemness enrichment. Notably, cAMP/Rap1 signaling and EMT
were negatively enriched in PDE2A-low bulk tumors in this run, so Figure4 should not
be written as uniform activation of all proposed malignant pathways. The bulk validation
should be treated as cohort-level external support because GSE178537 is tissue-level
RNA-seq and mixes malignant, stromal and immune components.
"""
    (OUT_DIR / "Fig4_v3_interpretation_report.md").write_text(report, encoding="utf-8")

    caption = """**Figure 4. PDE2A restricts malignant progression in HNSCC.**

**A.** UMAP of v3 annotated GSE198315 Edge/Core malignant epithelial cells showing PDE2A expression and PDE2A-negative/positive malignant-cell states. Because PDE2A expression is sparse in scRNA-seq, any detected expression was used to define the PDE2A-positive state.

**B.** Virtual PDE2A-loss volcano plot comparing PDE2A-negative versus PDE2A-positive malignant cells. Red points indicate genes activated in the PDE2A-negative state, blue points indicate genes repressed in the PDE2A-negative state.

**C.** Functional enrichment of PDE2A-loss activated genes using Enrichr GO/KEGG/Reactome when available, with a curated pathway fallback.

**D.** GSE178537 bulk RNA-seq preranked GSEA comparing PDE2A-low versus PDE2A-high tumor samples. Positive NES indicates enrichment in PDE2A-low samples.

**E.** Immune and TIDE-like signature scores in GSE178537 tumor samples stratified by PDE2A expression.
"""
    (OUT_DIR / "Fig4_v3_figure_caption.md").write_text(caption, encoding="utf-8")


def main() -> None:
    deg = run_sc_malignant()
    enr = enrich_de_genes(deg)
    gsea, immune_stats = bulk_gsea_and_immune()
    write_report(deg, enr, gsea, immune_stats)
    print(f"Done. Output: {OUT_DIR}")


if __name__ == "__main__":
    main()
