#!/usr/bin/env python3
"""Draw standalone GSE208253 spatial maps on the original H&E sections.

Each output contains one section and one variable. No multi-panel figures are
created. Signature definitions are frozen to the current PDE2A Figure 3 audit.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree


BASE = Path("/data/liaolz/PDE2A")
SAMPLES_FILE = BASE / "Figure2/Fig2_3_spatial_stemness/Fig2_3_spatial_samples.csv"
ANNOTATION_FILE = BASE / "results/Spatial_source_annotation_recovery_20260730/spatial_source_annotation_merged.csv"
DECONV_FILE = BASE / "results/spatial_deconvolution_he_mask_20260727/spatial_deconvolution_he_scores_merged.tsv.gz"
OUT = BASE / "results/GSE208253_standalone_spatial_maps_20260812"

GENE_SETS = {
    "stemness_basal": (
        "Stemness basal",
        ["S100A2", "KRT17", "KRT14", "KRT5", "KRT6A", "TP63", "ITGA6", "CD44", "JAG1", "FABP5", "EGFR"],
    ),
    "stemness_CSC": (
        "Stemness CSC",
        ["CD44", "ITGA6", "ALDH1A1", "ALDH1A3", "PROM1", "EPCAM", "TP63"],
    ),
    "proliferation": (
        "Proliferation",
        ["MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNB2", "MCM2", "MCM5", "MCM6", "UBE2C"],
    ),
    "EMT_invasion": (
        "EMT / invasion",
        ["VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "CDH2", "ITGA5", "ITGAV", "MMP2", "MMP9", "LAMC2", "COL1A1", "COL1A2", "SPARC", "TGFB1"],
    ),
    "leading_edge": (
        "Leading-edge program",
        ["LAMC2", "VIM", "FN1", "ITGA5", "MMP2", "MMP9", "COL1A1", "COL1A2", "SPARC", "TNC", "CD44"],
    ),
    "tumor_core": (
        "Tumor-core program",
        ["CLDN4", "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "TACSTD2", "MUC1"],
    ),
    "cAMP_related": (
        "cAMP-related program",
        ["ADCY1", "ADCY3", "ADCY7", "PRKACA", "PRKACB", "CREB1", "CREM", "PDE4D"],
    ),
    "MAPK": (
        "MAPK",
        ["DUSP6", "DUSP4", "SPRY2", "SPRY4", "EGR1", "FOS", "JUN", "ETV4", "ETV5"],
    ),
}

REGION_ORDER = ["nc", "core", "transitory", "edge", "unannotated"]
REGION_COLORS = {
    "nc": "#BDBDBD",
    "core": "#3C78A8",
    "transitory": "#E6AB02",
    "edge": "#D55E00",
    "unannotated": "#F2F2F2",
}
REGION_DISPLAY = {"nc": "NC", "core": "Core", "transitory": "Transition", "edge": "Edge", "unannotated": "NA"}

DECONV_FEATURES = {
    "malignant_abundance": "Malignant epithelial abundance",
    "stem_like_malignant_abundance": "Basal / stem-like malignant abundance",
    "fibroblast_abundance": "Fibroblast abundance",
    "endothelial_abundance": "Endothelial abundance",
    "myeloid_abundance": "Myeloid abundance",
    "T_NK_abundance": "T / NK abundance",
    "B_plasma_abundance": "B / plasma abundance",
}

mpl.use("Agg")
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def load_image(path: str) -> Image.Image:
    with gzip.open(path, "rb") as handle:
        return Image.open(io.BytesIO(handle.read())).convert("RGB")


def score_signature(adata: sc.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.full(adata.n_obs, np.nan), []
    x = dense(adata[:, present].X).astype(np.float32)
    mean = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0)
    sd[sd == 0] = 1.0
    return np.nanmean((x - mean) / sd, axis=1), present


def crop_limits(frame: pd.DataFrame, image: Image.Image) -> tuple[float, float, float, float]:
    x = frame["x_hires"].to_numpy(float)
    y = frame["y_hires"].to_numpy(float)
    width, height = image.size
    padx = max(24.0, 0.045 * (x.max() - x.min()))
    pady = max(24.0, 0.045 * (y.max() - y.min()))
    return (
        max(0.0, x.min() - padx),
        min(float(width), x.max() + padx),
        max(0.0, y.min() - pady),
        min(float(height), y.max() + pady),
    )


def configure_axis(ax, image: Image.Image, limits):
    width, height = image.size
    ax.imshow(image, extent=[0, width, height, 0], origin="upper", interpolation="nearest")
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[3], limits[2])
    ax.set_aspect("equal")
    ax.set_axis_off()


def spot_area(ax, frame: pd.DataFrame) -> float:
    points = frame[["x_hires", "y_hires"]].dropna().to_numpy(float)
    nearest = cKDTree(points).query(points, k=2)[0][:, 1]
    diameter_data = 0.80 * float(np.nanmedian(nearest))
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((diameter_data, 0))
    diameter_pt = abs(p1[0] - p0[0]) * 72.0 / ax.figure.dpi
    return np.pi * (diameter_pt / 2.0) ** 2


def figure_size(limits) -> tuple[float, float]:
    width = limits[1] - limits[0]
    height = limits[3] - limits[2]
    fig_width = 3.35
    fig_height = float(np.clip(fig_width * height / width + 0.68, 2.85, 4.75))
    return fig_width, fig_height


def save_figure(fig, folder: Path, stem: str):
    folder.mkdir(parents=True, exist_ok=True)
    common = dict(bbox_inches="tight", pad_inches=0.03, facecolor="white")
    fig.savefig(folder / f"{stem}.png", dpi=600, **common)
    fig.savefig(folder / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"}, **common)
    fig.savefig(folder / f"{stem}.pdf", **common)
    plt.close(fig)


def draw_continuous(frame, image, sample_label, key, label, limits, cmap="magma", group="continuous"):
    crop = crop_limits(frame, image)
    fig = plt.figure(figsize=figure_size(crop))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.16, 0.84], hspace=0.02)
    header = fig.add_subplot(grid[0])
    ax = fig.add_subplot(grid[1])
    header.set_axis_off()
    configure_axis(ax, image, crop)
    values = frame[key].to_numpy(float)
    valid = np.isfinite(values)
    points = ax.scatter(
        frame.loc[valid, "x_hires"],
        frame.loc[valid, "y_hires"],
        c=values[valid],
        s=spot_area(ax, frame),
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
        alpha=0.78,
        linewidths=0,
        rasterized=True,
    )
    header.text(0.0, 0.94, f"{sample_label}  |  {label}", transform=header.transAxes,
                ha="left", va="top", fontsize=10, fontweight="bold")
    header.text(1.0, 0.38, f"n = {valid.sum():,} spots", transform=header.transAxes,
                ha="right", va="center", fontsize=6.5, color="#444444")
    cax = header.inset_axes([0.02, 0.12, 0.39, 0.12])
    colorbar = fig.colorbar(points, cax=cax, orientation="horizontal")
    colorbar.set_ticks([limits[0], limits[1]])
    colorbar.ax.tick_params(labelsize=6.5, length=0, pad=1)
    colorbar.outline.set_visible(False)
    save_figure(fig, OUT / group / key, f"GSE208253_{sample_label}_{key}")


def draw_regions(frame, image, sample_label):
    crop = crop_limits(frame, image)
    fig = plt.figure(figsize=figure_size(crop))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.16, 0.84], hspace=0.02)
    header = fig.add_subplot(grid[0])
    ax = fig.add_subplot(grid[1])
    header.set_axis_off()
    configure_axis(ax, image, crop)
    labels = frame["source_author_region"].fillna("unannotated").astype(str)
    codes = pd.Categorical(labels, categories=REGION_ORDER).codes
    cmap = ListedColormap([REGION_COLORS[name] for name in REGION_ORDER])
    ax.scatter(
        frame["x_hires"], frame["y_hires"], c=codes,
        s=spot_area(ax, frame), cmap=cmap, vmin=-0.5, vmax=len(REGION_ORDER) - 0.5,
        alpha=0.82, linewidths=0, rasterized=True,
    )
    header.text(0.0, 0.94, f"{sample_label}  |  Author molecular regions",
                transform=header.transAxes, ha="left", va="top", fontsize=10, fontweight="bold")
    handles = [mpl.lines.Line2D([], [], marker="o", linestyle="", markersize=4.5,
                               markerfacecolor=REGION_COLORS[name], markeredgecolor="none",
                               label=REGION_DISPLAY[name])
               for name in REGION_ORDER]
    header.legend(handles=handles, ncol=5, loc="lower left", bbox_to_anchor=(0.0, 0.02),
                  frameon=False, fontsize=5.8, handletextpad=0.20, columnspacing=0.55,
                  borderaxespad=0)
    header.text(1.0, 0.54, f"n = {len(frame):,} spots", transform=header.transAxes,
                ha="right", va="center", fontsize=6.5, color="#444444")
    save_figure(fig, OUT / "author_molecular_regions", f"GSE208253_{sample_label}_author_regions")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(SAMPLES_FILE)
    annotations = pd.read_csv(
        ANNOTATION_FILE,
        usecols=["sample", "barcode", "x_hires", "y_hires", "source_author_region"],
    )
    deconv_columns = [
        "sample_original", "barcode_original", "Malignant_basal_squamous_stem_like",
        "Malignant_other", "Fibroblast", "Endothelial", "Macrophage", "DC",
        "Neutrophil", "Mast", "T_cell", "CD8T", "Treg", "Cycling_T", "NK",
        "B_cell", "Plasma",
    ]
    deconv = pd.read_csv(DECONV_FILE, sep="\t", usecols=deconv_columns)
    deconv = deconv.rename(columns={"sample_original": "sample", "barcode_original": "barcode"})
    deconv["malignant_abundance"] = (
        deconv["Malignant_basal_squamous_stem_like"] + deconv["Malignant_other"]
    )
    deconv["stem_like_malignant_abundance"] = deconv["Malignant_basal_squamous_stem_like"]
    deconv["fibroblast_abundance"] = deconv["Fibroblast"]
    deconv["endothelial_abundance"] = deconv["Endothelial"]
    deconv["myeloid_abundance"] = deconv[["Macrophage", "DC", "Neutrophil", "Mast"]].sum(axis=1)
    deconv["T_NK_abundance"] = deconv[["T_cell", "CD8T", "Treg", "Cycling_T", "NK"]].sum(axis=1)
    deconv["B_plasma_abundance"] = deconv[["B_cell", "Plasma"]].sum(axis=1)
    deconv = deconv[["sample", "barcode", *DECONV_FEATURES.keys()]]
    score_tables = []
    coverage = []

    for row in samples.sort_values("sample", key=lambda s: s.str.extract(r"(\d+)")[0].astype(int)).itertuples(index=False):
        section = sc.read_10x_h5(row.h5)
        section.var_names_make_unique()
        sc.pp.normalize_total(section, target_sum=1e4)
        sc.pp.log1p(section)
        frame = pd.DataFrame({"barcode": section.obs_names.astype(str)})
        frame["sample"] = row.sample
        if "PDE2A" in section.var_names:
            frame["PDE2A"] = dense(section[:, ["PDE2A"]].X).reshape(-1)
        else:
            frame["PDE2A"] = np.nan
        for key, (_, genes) in GENE_SETS.items():
            frame[key], present = score_signature(section, genes)
            coverage.append({
                "sample": row.sample,
                "feature": key,
                "n_requested": len(genes),
                "n_present": len(present),
                "present_genes": ";".join(present),
            })
        frame = frame.merge(
            annotations[annotations["sample"] == row.sample],
            on=["sample", "barcode"], how="inner", validate="one_to_one",
        )
        score_tables.append(frame)

    all_scores = pd.concat(score_tables, ignore_index=True)
    all_scores = all_scores.merge(deconv, on=["sample", "barcode"], how="left", validate="one_to_one")
    all_scores.to_csv(OUT / "GSE208253_standalone_spatial_source_data.tsv.gz", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(OUT / "signature_gene_coverage.tsv", sep="\t", index=False)

    feature_limits = {}
    for key in ["PDE2A", *GENE_SETS.keys(), *DECONV_FEATURES.keys()]:
        values = all_scores[key].to_numpy(float)
        values = values[np.isfinite(values)]
        qlo = 0.0 if key == "PDE2A" or key in DECONV_FEATURES else float(np.quantile(values, 0.02))
        qhi = float(np.quantile(values, 0.995 if key == "PDE2A" else 0.98))
        feature_limits[key] = (qlo, qhi)

    for row in samples.sort_values("sample", key=lambda s: s.str.extract(r"(\d+)")[0].astype(int)).itertuples(index=False):
        frame = all_scores[all_scores["sample"] == row.sample].copy()
        image = load_image(row.image)
        sample_label = f"S{int(row.sample.split('_')[-1])}"
        draw_regions(frame, image, sample_label)
        draw_continuous(frame, image, sample_label, "PDE2A", "PDE2A expression", feature_limits["PDE2A"])
        for key, (label, _) in GENE_SETS.items():
            draw_continuous(frame, image, sample_label, key, label, feature_limits[key])
        for key, label in DECONV_FEATURES.items():
            draw_continuous(
                frame, image, sample_label, key, label, feature_limits[key],
                cmap="viridis", group="deconvolution",
            )

    manifest = []
    for file in sorted(OUT.rglob("*")):
        if file.is_file():
            manifest.append({"file": str(file.relative_to(OUT)), "bytes": file.stat().st_size})
    pd.DataFrame(manifest).to_csv(OUT / "output_manifest.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# GSE208253 standalone spatial maps\n\n"
        "Each file contains one Visium section and one variable; no multi-panel figures were generated. "
        "Continuous overlays use log-normalized expression (10,000 counts per spot followed by log1p). "
        "Program scores are the mean of gene-wise z scores within each section, using the frozen PDE2A "
        "Figure 3 gene sets. Color limits are fixed across all 12 sections at robust cohort-wide quantiles. "
        "Deconvolution overlays reuse the current GSE198315 source-guided non-negative abundance estimates "
        "and show major malignant, stromal, endothelial, myeloid, lymphoid, and B/plasma compartments. "
        "Author molecular regions (core/transitory/edge/nc) come from the released GSE208253 Seurat "
        "objects and are molecular cluster annotations, not pathologist-drawn tumor-core or leading-edge ROIs. "
        "The H&E images and spot coordinates are the original GEO GSE208253 files.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
