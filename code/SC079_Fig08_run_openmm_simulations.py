"""Run an auditable short OpenMM MD simulation for PDE2A.

This script deliberately uses the PDE2A 6C7E crystal structure as the
structural starting point. The current project contains a docking gate and
PDE2A structure/ligand evidence, but no saved docking pose file. Small-molecule
protein-ligand production MD therefore remains blocked until ligand parameters
are generated with a proper small-molecule force field.
"""

from __future__ import annotations

import argparse
import json
import math
import time
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis import align, rms
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import (
    CheckpointReporter,
    CutoffNonPeriodic,
    DCDReporter,
    ForceField,
    HBonds,
    Modeller,
    NoCutoff,
    PDBFile,
    PME,
    StateDataReporter,
)
from openmm.unit import (
    angstroms,
    bar,
    kelvin,
    kilojoules_per_mole,
    molar,
    nanometer,
    picosecond,
)
from pdbfixer import PDBFixer


PDB_ID = "6C7E"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"


def ps_to_steps(ps: float, timestep_fs: float) -> int:
    return max(1, int(round(ps * 1000.0 / timestep_fs)))


def choose_platform(preferred: str) -> tuple[Platform, dict[str, str]]:
    if preferred != "auto":
        return Platform.getPlatformByName(preferred), {}
    for name, props in (("CUDA", {"Precision": "mixed"}), ("OpenCL", {}), ("CPU", {})):
        try:
            return Platform.getPlatformByName(name), props
        except Exception:
            continue
    raise RuntimeError("No OpenMM platform is available.")


def download_pdb(outdir: Path) -> Path:
    out = outdir / f"{PDB_ID}.pdb"
    if not out.exists():
        urllib.request.urlretrieve(PDB_URL, out)
    return out


def extract_ligand_records(pdb_path: Path) -> list[dict[str, str]]:
    ligands = {}
    ignored = {"HOH", "WAT"}
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname in ignored:
                continue
            chain = line[21].strip()
            resid = line[22:26].strip()
            ligands[(resname, chain, resid)] = {
                "resname": resname,
                "chain": chain,
                "resid": resid,
            }
    return list(ligands.values())


def write_chain_pdb(pdb_path: Path, outdir: Path, keep_chain: str) -> Path:
    chain_path = outdir / f"{PDB_ID}_chain_{keep_chain}.pdb"
    keep_prefixes = ("HEADER", "TITLE ", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "JRNL", "REMARK", "CRYST1")
    with pdb_path.open() as src, chain_path.open("w") as dst:
        for line in src:
            if line.startswith(keep_prefixes):
                dst.write(line)
            elif line.startswith(("ATOM  ", "ANISOU")) and line[21].strip() == keep_chain:
                dst.write(line)
            elif line.startswith("TER") and line[21].strip() == keep_chain:
                dst.write(line)
        dst.write("END\n")
    return chain_path


def prepare_protein(pdb_path: Path, outdir: Path, ph: float, keep_chain: str) -> Path:
    fixed_path = outdir / f"{PDB_ID}_protein_fixed.pdb"
    chain_path = write_chain_pdb(pdb_path, outdir, keep_chain)
    fixer = PDBFixer(filename=str(chain_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with fixed_path.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    return fixed_path


def build_system(protein_pdb: Path, outdir: Path, padding_angstrom: float, environment: str):
    pdb = PDBFile(str(protein_pdb))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    if environment == "water":
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=padding_angstrom * angstroms,
            ionicStrength=0.15 * molar,
        )
        nonbonded_method = PME
        system_pdb = outdir / f"{PDB_ID}_protein_solvated.pdb"
    else:
        nonbonded_method = CutoffNonPeriodic
        system_pdb = outdir / f"{PDB_ID}_protein_vacuum_system.pdb"
    with system_pdb.open("w") as handle:
        PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    return modeller, system, system_pdb


def run_md(args: argparse.Namespace) -> dict[str, object]:
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = download_pdb(outdir)
    ligands = extract_ligand_records(pdb_path)
    protein_pdb = prepare_protein(pdb_path, outdir, args.ph, args.chain)
    modeller, system, system_pdb = build_system(protein_pdb, outdir, args.padding_angstrom, args.environment)

    integrator = LangevinMiddleIntegrator(
        args.temperature * kelvin,
        1 / picosecond,
        args.timestep_fs / 1000.0 * picosecond,
    )
    platform, properties = choose_platform(args.platform)
    from openmm.app import Simulation

    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    t0 = time.time()
    simulation.minimizeEnergy(
        tolerance=args.min_tolerance * kilojoules_per_mole / nanometer,
        maxIterations=args.min_steps,
    )
    min_state = simulation.context.getState(getEnergy=True, getPositions=True)
    minimized = outdir / f"{PDB_ID}_protein_minimized.pdb"
    with minimized.open("w") as handle:
        PDBFile.writeFile(simulation.topology, min_state.getPositions(), handle)

    report_interval = max(1, ps_to_steps(args.report_ps, args.timestep_fs))
    simulation.context.setVelocitiesToTemperature(args.temperature * kelvin)
    simulation.reporters = [
        StateDataReporter(
            str(outdir / "nvt_log.tsv"),
            report_interval,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            speed=True,
            separator="\t",
        ),
        DCDReporter(str(outdir / "nvt_traj.dcd"), report_interval),
    ]
    nvt_steps = ps_to_steps(args.nvt_ps, args.timestep_fs)
    simulation.step(nvt_steps)

    if args.environment == "water":
        system.addForce(MonteCarloBarostat(args.pressure_bar * bar, args.temperature * kelvin, 25))
        simulation.context.reinitialize(preserveState=True)
    simulation.reporters = [
        StateDataReporter(
            str(outdir / "production_log.tsv"),
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
            separator="\t",
        ),
        DCDReporter(str(outdir / "production_traj.dcd"), report_interval),
        CheckpointReporter(str(outdir / "production_checkpoint.chk"), report_interval),
    ]
    production_steps = ps_to_steps(args.production_ps, args.timestep_fs)
    simulation.step(production_steps)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_pdb = outdir / f"{PDB_ID}_protein_production_final.pdb"
    with final_pdb.open("w") as handle:
        PDBFile.writeFile(simulation.topology, final_state.getPositions(), handle)

    return {
        "pdb_id": PDB_ID,
        "pdb_url": PDB_URL,
        "ligands_in_crystal_after_water_exclusion": ligands,
        "md_scope": "protein-only apo simulation from 6C7E after removing heterogens",
        "chain": args.chain,
        "environment": args.environment,
        "ligand_md_status": "blocked: no saved docking pose and no ligand parameters generated in project",
        "force_field": "AMBER14 protein; TIP3P-FB water only for explicit-water mode",
        "platform": platform.getName(),
        "temperature_K": args.temperature,
        "pressure_bar": args.pressure_bar,
        "timestep_fs": args.timestep_fs,
        "nvt_ps": args.nvt_ps,
        "production_ps": args.production_ps,
        "nvt_steps": nvt_steps,
        "production_steps": production_steps,
        "initial_potential_energy": str(initial),
        "minimized_potential_energy": str(min_state.getPotentialEnergy()),
        "final_potential_energy": str(final_state.getPotentialEnergy()),
        "elapsed_seconds": round(time.time() - t0, 3),
        "outputs": {
            "source_pdb": str(pdb_path),
            "protein_fixed_pdb": str(protein_pdb),
            "system_pdb": str(system_pdb),
            "minimized_pdb": str(minimized),
            "nvt_log": str(outdir / "nvt_log.tsv"),
            "nvt_traj": str(outdir / "nvt_traj.dcd"),
            "production_log": str(outdir / "production_log.tsv"),
            "production_traj": str(outdir / "production_traj.dcd"),
            "final_pdb": str(final_pdb),
        },
    }


def analyze(outdir: Path, report_ps: float, environment: str) -> dict[str, object]:
    top = outdir / (
        f"{PDB_ID}_protein_solvated.pdb"
        if environment == "water"
        else f"{PDB_ID}_protein_vacuum_system.pdb"
    )
    traj = outdir / "production_traj.dcd"
    if not traj.exists() or traj.stat().st_size == 0:
        return {"analysis_status": "skipped: production trajectory missing"}

    universe = Universe(str(top), str(traj))
    protein = universe.select_atoms("protein and backbone")
    if protein.n_atoms == 0 or len(universe.trajectory) < 2:
        return {"analysis_status": "skipped: insufficient trajectory/protein atoms"}

    align.AlignTraj(universe, universe, select="protein and backbone", in_memory=True).run()
    rmsd_runner = rms.RMSD(universe, select="protein and backbone", ref_frame=0).run()
    rmsd = rmsd_runner.results.rmsd
    rmsd_path = outdir / "production_backbone_rmsd.tsv"
    np.savetxt(
        rmsd_path,
        rmsd,
        delimiter="\t",
        header="frame\ttime_ps\tbackbone_rmsd_A",
        comments="",
    )
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    ax.plot(rmsd[:, 1] / 1000.0, rmsd[:, 2], color="#7C3AED", linewidth=1.3)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Backbone RMSD (A)")
    ax.set_title("PDE2A 6C7E apo pilot MD")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    rmsd_png = outdir / "production_backbone_rmsd.png"
    fig.savefig(rmsd_png, dpi=200)
    plt.close(fig)

    return {
        "analysis_status": "complete",
        "trajectory_frames": len(universe.trajectory),
        "nominal_frame_interval_ps": report_ps,
        "backbone_rmsd_A_mean": float(np.mean(rmsd[:, 2])),
        "backbone_rmsd_A_last": float(rmsd[-1, 2]),
        "backbone_rmsd_table": str(rmsd_path),
        "backbone_rmsd_plot": str(rmsd_png),
    }


def write_report(outdir: Path, manifest: dict[str, object]) -> None:
    report = outdir / "PDE2A_md_audit_report.md"
    lines = [
        "# PDE2A MD audit report",
        "",
        "## Docking traceback",
        "",
        "- Local project search found a docking gate and PDE2A PDB/BindingDB evidence, but no saved docking pose file (`pdbqt`, `sdf`, `mol2`, Vina/Glide output, or protein-ligand complex from docking).",
        "- The current runnable structure is therefore the PDE2A 6C7E crystal complex, not a de novo docking result.",
        "- Protein-ligand production MD is not claimed here because ligand force-field parameters were not generated.",
        "",
        "## Completed simulation",
        "",
        f"- Scope: {manifest['md_scope']}",
        f"- PDB: {manifest['pdb_id']} ({manifest['pdb_url']})",
        f"- Force field: {manifest['force_field']}",
        f"- Platform: {manifest['platform']}",
        f"- Environment: {manifest['environment']}",
        f"- NVT: {manifest['nvt_ps']} ps; production: {manifest['production_ps']} ps",
        f"- Initial potential energy: {manifest['initial_potential_energy']}",
        f"- Minimized potential energy: {manifest['minimized_potential_energy']}",
        f"- Final potential energy: {manifest['final_potential_energy']}",
        "",
        "## Outputs",
        "",
    ]
    for key, value in manifest["outputs"].items():
        lines.append(f"- {key}: `{value}`")
    lines.extend(
        [
            "",
            "## Analysis",
            "",
            f"- Status: {manifest['analysis']['analysis_status']}",
        ]
    )
    if manifest["analysis"].get("analysis_status") == "complete":
        lines.extend(
            [
                f"- Frames: {manifest['analysis']['trajectory_frames']}",
                f"- Mean backbone RMSD: {manifest['analysis']['backbone_rmsd_A_mean']:.3f} A",
                f"- Last-frame backbone RMSD: {manifest['analysis']['backbone_rmsd_A_last']:.3f} A",
                f"- RMSD table: `{manifest['analysis']['backbone_rmsd_table']}`",
                f"- RMSD plot: `{manifest['analysis']['backbone_rmsd_plot']}`",
            ]
        )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results/PDE2A_MD_6C7E_openmm"))
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--pressure-bar", type=float, default=1.0)
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--padding-angstrom", type=float, default=8.0)
    parser.add_argument("--environment", choices=["water", "vacuum"], default="water")
    parser.add_argument("--platform", choices=["auto", "CUDA", "OpenCL", "CPU"], default="auto")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--min-steps", type=int, default=500)
    parser.add_argument("--min-tolerance", type=float, default=20.0)
    parser.add_argument("--nvt-ps", type=float, default=2.0)
    parser.add_argument("--production-ps", type=float, default=5.0)
    parser.add_argument("--report-ps", type=float, default=0.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = run_md(args)
    manifest["analysis"] = analyze(args.outdir, args.report_ps, args.environment)
    manifest_path = args.outdir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    write_report(args.outdir, manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
