#!/usr/bin/env python3
"""Publication-grade spatial maps for the significant PDE2A analyses.

Inputs are the audited GSE208253 spot table and the original Visium images/H5
files on the analysis server. The maps are descriptive overlays; inferential
claims remain in the paired summary figures.
"""
from pathlib import Path
import gzip
import io
import json
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=FutureWarning)
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

SPOTS = Path("/data/liaolz/PDE2A/Figure2/Fig2_3_spatial_edge_core_panels/Fig2_3_spatial_edge_core_plot_data.csv")
SAMPLES = Path("/data/liaolz/PDE2A/Figure2/Fig2_3_spatial_stemness/Fig2_3_spatial_samples.csv")
OUT = Path("/data/liaolz/PDE2A/results/FigS_spatial_maps_significant_redraw")
REPRESENTATIVE = ["sample_9", "sample_6", "sample_7"]

STATE_COLORS = {"Other": "#D9D9D9", "Core-like": "#0072B2", "Edge-like": "#D55E00"}
STATE_ORDER = ["Other", "Core-like", "Edge-like"]


def dense_one(x):
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).reshape(-1)


def load_hires(path: str):
    with gzip.open(path, "rb") as fh:
        return Image.open(io.BytesIO(fh.read())).convert("RGB")


def read_pde2a(sample_rows):
    rows = []
    for _, r in sample_rows.iterrows():
        a = sc.read_h5ad(r["h5"]) if str(r["h5"]).endswith(".h5ad") else sc.read_10x_h5(r["h5"])
        a.var_names_make_unique()
        # Match the audited spatial analysis: library-size normalization + log1p.
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        vals = dense_one(a[:, ["PDE2A"]].X) if "PDE2A" in a.var_names else np.full(a.n_obs, np.nan)
        rows.append(pd.DataFrame({"sample": r.name, "barcode": a.obs_names.astype(str), "PDE2A_expr": vals}))
    return pd.concat(rows, ignore_index=True)


def scale_limits(df, col, q=(0.01, 0.99), symmetric=False):
    x = pd.to_numeric(df[col], errors="coerce").to_numpy()
    x = x[np.isfinite(x)]
    lo, hi = np.quantile(x, q)
    if symmetric:
        m = max(abs(lo), abs(hi))
        return -m, m
    return float(lo), float(hi)


def crop_limits(sample_df, image):
    x = sample_df["x_hires"].to_numpy(float)
    y = sample_df["y_hires"].to_numpy(float)
    w, h = image.size
    # Keep a compact field of view while preserving enough tissue context.
    pad_x = max(28.0, 0.07 * (x.max() - x.min()))
    pad_y = max(28.0, 0.07 * (y.max() - y.min()))
    return (max(0, x.min() - pad_x), min(w, x.max() + pad_x),
            max(0, y.min() - pad_y), min(h, y.max() + pad_y))


def spot_area(ax, d, fill_ratio=0.82):
    """Return a scatter area whose circle is slightly smaller than one spot pitch."""
    pts = d[["x_hires", "y_hires"]].dropna().to_numpy(float)
    nn = float(np.median(cKDTree(pts).query(pts, k=2)[0][:, 1]))
    x0, y0 = ax.transData.transform((ax.get_xlim()[0], ax.get_ylim()[0]))
    x1, _ = ax.transData.transform((ax.get_xlim()[0] + fill_ratio * nn, ax.get_ylim()[0]))
    diameter_pt = abs(x1 - x0) * 72.0 / ax.figure.dpi
    return np.pi * (diameter_pt / 2.0) ** 2


def scatter_spots(ax, d, col, vmin=None, vmax=None, cmap="viridis", alpha=0.9):
    x, y = d["x_hires"].to_numpy(float), d["y_hires"].to_numpy(float)
    vals = pd.to_numeric(d[col], errors="coerce").to_numpy(float)
    mask = np.isfinite(vals)
    return ax.scatter(x[mask], y[mask], c=vals[mask], s=spot_area(ax, d), cmap=cmap,
                      vmin=vmin, vmax=vmax, linewidths=0, alpha=alpha, rasterized=True)


def add_background(ax, image, lims):
    w, h = image.size
    ax.imshow(image, extent=[0, w, h, 0], origin="upper", interpolation="nearest", alpha=0.94)
    ax.set_xlim(lims[0], lims[1]); ax.set_ylim(lims[3], lims[2])
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def add_state_legend(ax):
    handles = [mpl.lines.Line2D([0], [0], marker="o", linestyle="", markersize=5,
                                markerfacecolor=STATE_COLORS[s], markeredgecolor="none", label=s)
               for s in STATE_ORDER]
    ax.legend(handles=handles, title="Spatial state", frameon=False, loc="upper left",
              bbox_to_anchor=(0, 1.02), ncol=3, fontsize=7, title_fontsize=7)


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_map_figure(df, samples, stem, value_cols, titles, cmaps, limits, state=False):
    nrow, ncol = len(REPRESENTATIVE), len(value_cols)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.05 * ncol, 2.75 * nrow), squeeze=False)
    mappables = []
    for i, sample in enumerate(REPRESENTATIVE):
        d = df[df["sample"] == sample].copy()
        image = load_hires(samples.loc[sample, "image"])
        lims = crop_limits(d, image)
        for j, (col, title, cmap) in enumerate(zip(value_cols, titles, cmaps)):
            ax = axes[i, j]
            add_background(ax, image, lims)
            if state and col == "edge_core_region":
                # Draw the neutral background first, then the two inferred states.
                for cat in STATE_ORDER:
                    q = d[d[col] == cat]
                    if len(q):
                        ax.scatter(q["x_hires"], q["y_hires"], s=spot_area(ax, d), c=STATE_COLORS[cat],
                                   linewidths=0, alpha=0.95, rasterized=True)
            else:
                mappables.append(scatter_spots(ax, d, col, limits[col][0], limits[col][1], cmap))
            if i == 0:
                ax.set_title(title, pad=8, fontweight="bold")
            if j == 0:
                ax.text(-0.02, 1.06, sample.replace("sample_", "S"), transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=9, fontweight="bold")
            if state and col == "edge_core_region" and i == 0:
                add_state_legend(ax)
    # One colorbar per continuous column, shared across samples.
    for j, col in enumerate(value_cols):
        if state and col == "edge_core_region":
            continue
        sm = mappables.pop(0)
        cbar = fig.colorbar(sm, ax=axes[:, j].tolist(), fraction=0.035, pad=0.02, aspect=22)
        cbar.ax.tick_params(labelsize=7, length=2)
        cbar.set_label({"stem_score": "Stemness score", "KRT17": "KRT17 (log1p)",
                        "S100A2": "S100A2 (log1p)", "PDE2A_expr": "PDE2A (log1p normalized)"}.get(col, col), fontsize=8)
    fig.subplots_adjust(wspace=0.06, hspace=0.18, left=0.02, right=0.95, top=0.91, bottom=0.04)
    fig.suptitle("Representative spatial maps: audited HNSCC Visium cohort", y=0.985,
                 fontsize=12, fontweight="bold")
    save(fig, stem)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    spots = pd.read_csv(SPOTS, low_memory=False)
    samples = pd.read_csv(SAMPLES).set_index("sample")
    samples = samples.loc[REPRESENTATIVE]
    df = spots[spots["sample"].isin(REPRESENTATIVE)].copy()
    pde = read_pde2a(samples)
    df = df.merge(pde, on=["sample", "barcode"], how="left")
    df.to_csv(OUT / "spatial_map_data_used.csv", index=False)

    # Fixed cohort-wide scales prevent visual exaggeration between samples.
    limits = {
        "stem_score": scale_limits(df, "stem_score", symmetric=True),
        "KRT17": scale_limits(df, "KRT17", q=(0.01, 0.995)),
        "S100A2": scale_limits(df, "S100A2", q=(0.01, 0.995)),
        "PDE2A_expr": (0.0, float(np.nanquantile(df["PDE2A_expr"], 0.995))),
    }
    make_map_figure(
        df, samples, "FigS_spatial_maps_stemness_marker",
        ["edge_core_region", "stem_score", "KRT17", "S100A2"],
        ["H&E + inferred state", "Stemness", "KRT17", "S100A2"],
        ["viridis", "coolwarm", "magma", "magma"], limits, state=True,
    )
    make_map_figure(
        df, samples, "FigS_spatial_maps_PDE2A_context",
        ["edge_core_region", "stem_score", "PDE2A_expr"],
        ["H&E + inferred state", "Stemness", "PDE2A expression"],
        ["viridis", "coolwarm", "viridis"], limits, state=True,
    )
    (OUT / "README.md").write_text(
        "# Spatial map redraws\n\n"
        "Representative samples S9, S6 and S7 from the audited GSE208253 Visium cohort. "
        "State labels are marker-defined LE-like/TC-like proxies (not manual pathology masks). "
        "PDE2A uses the same library-size normalization and log1p transform as the spatial audit. "
        "Coordinates and source files are recorded in spatial_map_data_used.csv.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
