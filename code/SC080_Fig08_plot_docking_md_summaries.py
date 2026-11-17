#!/usr/bin/env python3
"""Create standalone publication figures from the PDE2A-CGP60474 results."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import capped_distance
import numpy as np
from PIL import Image


PURPLE = "#7C3AED"
PURPLE_DARK = "#5B21B6"
GRAY = "#666666"
LIGHT_GRAY = "#B8B8B8"
BLACK = "#222222"


def parse_args() -> argparse.Namespace:
    desktop = Path.cwd().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-root",
        type=Path,
        default=desktop / "t........3" / "t........3",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd()
        / "results"
        / "PDE2A_CGP60474_publication_figures_20260807",
    )
    return parser.parse_args()


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "legend.frameon": False,
            "lines.linewidth": 1.15,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
        }
    )


def read_xvg(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments=("#", "@"))


def polish_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.7, color=BLACK)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.pdf")
    fig.savefig(output_dir / f"{stem}.png", dpi=600)
    plt.close(fig)


def window_stats(data: np.ndarray, start_ns: float = 20.0) -> dict[str, float]:
    values = data[data[:, 0] >= start_ns, 1]
    return {
        "n_frames": int(values.size),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "median": float(np.median(values)),
        "p05": float(np.percentile(values, 5)),
        "p95": float(np.percentile(values, 95)),
        "final": float(data[-1, 1]),
    }


def parse_docking_modes(path: Path) -> list[dict[str, float]]:
    modes = []
    pattern = re.compile(
        r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$"
    )
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            modes.append(
                {
                    "mode": int(match.group(1)),
                    "affinity_kcal_mol": float(match.group(2)),
                    "rmsd_lb_A": float(match.group(3)),
                    "rmsd_ub_A": float(match.group(4)),
                }
            )
    return modes


def analyze_ligand_trajectory(
    topology: Path,
    trajectory: Path,
    rmsf_complex: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    # The XTC C reader cannot open this dataset through a Unicode absolute path
    # on Windows, while relative paths from the same directory work correctly.
    os.chdir(topology.parent)
    universe = mda.Universe(topology.name, trajectory.name)
    ligand_selection = "resname MOL and not name H*"

    rmsd_analysis = rms.RMSD(
        universe,
        universe,
        select="backbone",
        groupselections=[ligand_selection],
        ref_frame=0,
    ).run(verbose=False)
    rmsd_result = rmsd_analysis.results.rmsd
    ligand_rmsd = np.column_stack((rmsd_result[:, 1] / 1000.0, rmsd_result[:, 3]))

    ligand = universe.select_atoms(ligand_selection)
    protein_heavy = universe.select_atoms("protein and not name H*")
    min_distance = np.empty(len(universe.trajectory), dtype=float)
    contact_frames: dict[int, int] = {}

    for frame_index, ts in enumerate(universe.trajectory):
        pairs, distances = capped_distance(
            ligand.positions,
            protein_heavy.positions,
            max_cutoff=4.0,
            box=ts.dimensions,
            return_distances=True,
        )
        min_distance[frame_index] = np.min(distances) if distances.size else np.nan
        if pairs.size:
            residue_indices = set(protein_heavy[pairs[:, 1]].resindices)
            for residue_index in residue_indices:
                contact_frames[residue_index] = contact_frames.get(residue_index, 0) + 1

    first_xvg_residue = int(rmsf_complex[0, 0])
    first_tpr_residue = int(universe.select_atoms("protein").residues[0].resid)
    residue_offset = first_xvg_residue - first_tpr_residue
    contacts = []
    frame_count = len(universe.trajectory)
    for residue_index, count in contact_frames.items():
        residue = universe.residues[residue_index]
        contacts.append(
            {
                "residue": f"{residue.resname}{int(residue.resid) + residue_offset}",
                "occupancy_fraction": count / frame_count,
                "frames": count,
            }
        )
    contacts.sort(key=lambda row: row["occupancy_fraction"], reverse=True)

    min_distance_data = np.column_stack((ligand_rmsd[:, 0], min_distance))
    return ligand_rmsd, min_distance_data, contacts


def plot_protein_rmsd(
    complex_data: np.ndarray, apo_data: np.ndarray, output_dir: Path
) -> None:
    complex_A = complex_data.copy()
    apo_A = apo_data.copy()
    complex_A[:, 1] *= 10.0
    apo_A[:, 1] *= 10.0
    complex_stats = window_stats(complex_A)
    apo_stats = window_stats(apo_A)

    fig, ax = plt.subplots(figsize=(3.54, 2.65))
    ax.plot(
        complex_A[:, 0],
        complex_A[:, 1],
        color=PURPLE,
        label=f"PDE2A-CGP60474 ({complex_stats['mean']:.2f} $\\AA$)",
    )
    ax.plot(
        apo_A[:, 0],
        apo_A[:, 1],
        color=GRAY,
        linestyle="--",
        label=f"Apo PDE2A ({apo_stats['mean']:.2f} $\\AA$)",
    )
    ax.set(xlabel="Time (ns)", ylabel=r"Backbone RMSD ($\AA$)", xlim=(0, 100))
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.legend(loc="lower right", title="Mean, 20-100 ns", title_fontsize=7.5)
    polish_axis(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_protein_backbone_RMSD_100ns")


def plot_ligand_rmsd(ligand_rmsd: np.ndarray, output_dir: Path) -> None:
    stats = window_stats(ligand_rmsd)
    fig, ax = plt.subplots(figsize=(3.54, 2.65))
    ax.plot(ligand_rmsd[:, 0], ligand_rmsd[:, 1], color=PURPLE)
    ax.axhline(
        stats["mean"],
        color=PURPLE_DARK,
        linestyle="--",
        linewidth=0.9,
        label=f"20-100 ns mean: {stats['mean']:.2f} $\\AA$",
    )
    ax.set(
        xlabel="Time (ns)",
        ylabel=r"Ligand heavy-atom RMSD ($\AA$)",
        xlim=(0, 100),
    )
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.legend(loc="lower right")
    polish_axis(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_ligand_pose_RMSD_100ns")


def plot_min_distance(min_distance: np.ndarray, output_dir: Path) -> None:
    stats = window_stats(min_distance)
    fig, ax = plt.subplots(figsize=(3.54, 2.65))
    ax.plot(min_distance[:, 0], min_distance[:, 1], color=PURPLE, linewidth=0.9)
    ax.axhline(
        stats["mean"],
        color=PURPLE_DARK,
        linestyle="--",
        linewidth=0.9,
        label=f"20-100 ns mean: {stats['mean']:.2f} $\\AA$",
    )
    ax.set(
        xlabel="Time (ns)",
        ylabel=r"Minimum protein-ligand distance ($\AA$)",
        xlim=(0, 100),
    )
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.legend(loc="upper right")
    polish_axis(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_protein_ligand_minimum_distance_100ns")


def plot_contacts(contacts: list[dict[str, float]], output_dir: Path) -> None:
    selected = [row for row in contacts if row["occupancy_fraction"] >= 0.30]
    selected = list(reversed(selected))
    labels = [row["residue"] for row in selected]
    values = [100.0 * row["occupancy_fraction"] for row in selected]

    height = max(2.8, 0.23 * len(selected) + 0.9)
    fig, ax = plt.subplots(figsize=(3.54, height))
    bars = ax.barh(labels, values, color=PURPLE, height=0.66)
    ax.axvline(30, color=GRAY, linestyle="--", linewidth=0.8)
    for bar, value in zip(bars, values):
        ax.text(
            min(value + 1.2, 98.5),
            bar.get_y() + bar.get_height() / 2,
            f"{value:.0f}%",
            va="center",
            ha="left" if value < 95 else "right",
            fontsize=7,
            color=BLACK,
        )
    ax.set(
        xlabel=r"Frames with heavy-atom contact $\leq$ 4.0 $\AA$ (%)",
        xlim=(0, 105),
    )
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.grid(axis="y", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.7)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_residue_contact_occupancy_100ns")


def plot_rmsf(
    complex_data: np.ndarray, apo_data: np.ndarray, output_dir: Path
) -> None:
    complex_A = complex_data.copy()
    apo_A = apo_data.copy()
    complex_A[:, 1] *= 10.0
    apo_A[:, 1] *= 10.0

    fig, ax = plt.subplots(figsize=(5.5, 2.65))
    ax.plot(
        complex_A[:, 0],
        complex_A[:, 1],
        color=PURPLE,
        label="PDE2A-CGP60474",
    )
    ax.plot(
        apo_A[:, 0],
        apo_A[:, 1],
        color=GRAY,
        linestyle="--",
        label="Apo PDE2A",
    )
    ax.set(xlabel="Residue number", ylabel=r"Residue RMSF ($\AA$)")
    ax.legend(loc="upper center", ncol=2)
    polish_axis(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_residue_RMSF_100ns")


def plot_hydrogen_bonds(hbond_data: np.ndarray, output_dir: Path) -> None:
    frame_interval_ns = float(np.median(np.diff(hbond_data[:, 0])))
    window = max(1, int(round(1.0 / frame_interval_ns)))
    kernel = np.ones(window) / window
    rolling = np.convolve(hbond_data[:, 1], kernel, mode="same")
    stats = window_stats(hbond_data)

    fig, ax = plt.subplots(figsize=(3.54, 2.65))
    ax.plot(
        hbond_data[:, 0],
        hbond_data[:, 1],
        color=LIGHT_GRAY,
        linewidth=0.55,
        alpha=0.55,
        label="Per frame",
    )
    ax.plot(
        hbond_data[:, 0],
        rolling,
        color=PURPLE,
        linewidth=1.25,
        label="1 ns moving mean",
    )
    ax.set(
        xlabel="Time (ns)",
        ylabel="Protein-ligand hydrogen bonds",
        xlim=(0, 100),
        ylim=(0, max(5.3, np.max(hbond_data[:, 1]) + 0.3)),
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        borderaxespad=0,
    )
    ax.text(
        0.98,
        0.95,
        f"20-100 ns mean: {stats['mean']:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
    )
    polish_axis(ax)
    fig.tight_layout()
    save_figure(fig, output_dir, "Fig_MD_protein_ligand_hydrogen_bonds_100ns")


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    source_root = args.source_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    set_style()

    data_root = source_root / "Dynamics" / "DATA"
    complex_data = data_root / "fuhewu"
    apo_data = data_root / "pure-protein"
    full_complex = source_root / "Dynamics" / "fuhewu"

    protein_rmsd_complex = read_xvg(complex_data / "rmsd.xvg")
    protein_rmsd_apo = read_xvg(apo_data / "rmsd.xvg")
    rmsf_complex = read_xvg(complex_data / "rmsf.xvg")
    rmsf_apo = read_xvg(apo_data / "rmsf.xvg")
    hbond_data = read_xvg(complex_data / "hbnum.xvg")

    metrics_path = output_dir / "ligand_pose_metrics.tsv"
    contacts_path = output_dir / "residue_contact_occupancy.tsv"
    if metrics_path.exists() and contacts_path.exists():
        cached_metrics = np.loadtxt(metrics_path, delimiter="\t", skiprows=1)
        ligand_rmsd = cached_metrics[:, :2]
        min_distance = cached_metrics[:, (0, 2)]
        with contacts_path.open(encoding="utf-8", newline="") as handle:
            contacts = []
            for row in csv.DictReader(handle, delimiter="\t"):
                contacts.append(
                    {
                        "residue": row["residue"],
                        "occupancy_fraction": float(row["occupancy_fraction"]),
                        "frames": int(row["frames"]),
                    }
                )
    else:
        ligand_rmsd, min_distance, contacts = analyze_ligand_trajectory(
            full_complex / "md_0_1.tpr",
            full_complex / "md_0_1_center.xtc",
            rmsf_complex,
        )

    plot_protein_rmsd(protein_rmsd_complex, protein_rmsd_apo, output_dir)
    plot_ligand_rmsd(ligand_rmsd, output_dir)
    plot_min_distance(min_distance, output_dir)
    plot_contacts(contacts, output_dir)
    plot_rmsf(rmsf_complex, rmsf_apo, output_dir)
    plot_hydrogen_bonds(hbond_data, output_dir)

    docking_image = (
        source_root
        / "Docking"
        / "Picture"
        / "Snipaste_2026-08-03_16-13-54.jpg"
    )
    with Image.open(docking_image) as image:
        image.save(
            output_dir / "Fig_Docking_predicted_binding_pose_CGP60474_PDE2A.png",
            dpi=(300, 300),
        )

    np.savetxt(
        output_dir / "ligand_pose_metrics.tsv",
        np.column_stack((ligand_rmsd, min_distance[:, 1])),
        delimiter="\t",
        header="time_ns\tligand_heavy_atom_RMSD_A\tminimum_distance_A",
        comments="",
        fmt="%.6f",
    )
    write_tsv(
        output_dir / "residue_contact_occupancy.tsv",
        contacts,
        ["residue", "occupancy_fraction", "frames"],
    )

    docking_modes = parse_docking_modes(
        source_root / "Docking" / "Results" / "PDE2A_2_CGP60474.txt"
    )
    write_tsv(
        output_dir / "docking_modes.tsv",
        docking_modes,
        ["mode", "affinity_kcal_mol", "rmsd_lb_A", "rmsd_ub_A"],
    )

    summary_rows = []
    datasets = {
        "protein_RMSD_complex_A": protein_rmsd_complex * np.array([1.0, 10.0]),
        "protein_RMSD_apo_A": protein_rmsd_apo * np.array([1.0, 10.0]),
        "ligand_RMSD_after_protein_fit_A": ligand_rmsd,
        "minimum_protein_ligand_distance_A": min_distance,
        "protein_ligand_hydrogen_bonds": hbond_data,
    }
    for metric, data in datasets.items():
        row = {"metric": metric, **window_stats(data)}
        summary_rows.append(row)
    write_tsv(
        output_dir / "summary_statistics_20_100ns.tsv",
        summary_rows,
        ["metric", "n_frames", "mean", "sd", "median", "p05", "p95", "final"],
    )

    print(f"Created publication figures in: {output_dir}")


if __name__ == "__main__":
    main()
