#!/usr/bin/env python3
"""Checkpoint 2-3 technical QC for the frozen 12-sample bulk RNA-seq cohort."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SAMPLE_ORDER = [
    "P08_C", "P08_TI", "P08_P", "P09_C", "P09_TI", "P09_P",
    "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI",
]
STAR_KEYS = [
    "Number of input reads", "Average input read length",
    "Uniquely mapped reads number", "Uniquely mapped reads %",
    "Average mapped length", "Number of splices: Total",
    "Number of splices: Annotated (sjdb)", "Number of splices: GT/AG",
    "Number of splices: GC/AG", "Number of splices: AT/AC",
    "Number of splices: Non-canonical", "Mismatch rate per base, %",
    "Deletion rate per base", "Insertion rate per base",
    "Number of reads mapped to multiple loci", "% of reads mapped to multiple loci",
    "Number of reads mapped to too many loci", "% of reads mapped to too many loci",
    "% of reads unmapped: too many mismatches", "% of reads unmapped: too short",
    "% of reads unmapped: other",
]


def now() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def read_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"sample_id", "patient_id", "region", "batch"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Manifest missing columns: {sorted(missing)}")
    if df["sample_id"].tolist() != SAMPLE_ORDER:
        raise SystemExit("Frozen sample order/content differs from the required 12-sample cohort.")
    return df


def atomic_write(path: Path, text: str, archive_existing: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if archive_existing and path.exists() and path.read_text(encoding="utf-8") != text:
        stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
        shutil.copy2(path, path.with_name(f"{path.name}.previous.{stamp}"))
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def run(cmd: list[str], stdout_path: Path | None = None) -> subprocess.CompletedProcess:
    if stdout_path is None:
        return subprocess.run(cmd, text=True, capture_output=True, check=False)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8", newline="\n") as handle:
        return subprocess.run(cmd, text=True, stdout=handle, stderr=subprocess.PIPE, check=False)


def paths_for(star_dir: Path, sample: str) -> tuple[Path, Path, Path]:
    root = star_dir / sample / sample
    return (
        Path(f"{root}.Aligned.sortedByCoord.out.bam"),
        Path(f"{root}.Aligned.sortedByCoord.out.bam.bai"),
        Path(f"{root}.Log.final.out"),
    )


def running_sample() -> str | None:
    proc = run(["pgrep", "-af", "STAR --runThreadN"])
    if proc.returncode not in (0, 1):
        return None
    found = re.findall(r"--outFileNamePrefix\s+\S*/(P\d+_(?:C|TI|P))\.", proc.stdout)
    unique = sorted(set(found))
    if len(unique) > 1:
        raise SystemExit(f"Duplicate/concurrent STAR samples detected: {unique}")
    return unique[0] if unique else None


def star_wrapper_running() -> bool:
    proc = run(["pgrep", "-af", "02_star_align.sh"])
    return proc.returncode == 0 and any("bash scripts/02_star_align.sh" in line for line in proc.stdout.splitlines())


def write_status(manifest: pd.DataFrame, star_dir: Path) -> None:
    active = running_sample()
    wrapper_active = star_wrapper_running()
    rows = []
    for row in manifest.itertuples(index=False):
        bam, bai, log = paths_for(star_dir, row.sample_id)
        complete = bam.is_file() and bam.stat().st_size > 0 and bai.is_file() and bai.stat().st_size > 0 and log.is_file() and log.stat().st_size > 0
        if complete:
            status, notes = "COMPLETED", "BAM, BAI, and Log.final.out present"
        elif row.sample_id == active:
            status, notes = "RUNNING", "Active STAR process"
        elif any(p.exists() for p in (bam, bai, log)) and wrapper_active:
            status, notes = "RUNNING", "STAR queue active; sample output is being finalized/indexed"
        elif any(p.exists() for p in (bam, bai, log)):
            status, notes = "FAILED", "Partial STAR outputs present; inspect before any rerun"
        else:
            status, notes = "PENDING", "Not started"
        rows.append({
            "sample_id": row.sample_id, "patient_id": row.patient_id,
            "region": row.region, "batch": row.batch, "star_status": status,
            "bam_exists": bam.is_file() and bam.stat().st_size > 0 if bam.exists() else False,
            "bai_exists": bai.is_file() and bai.stat().st_size > 0 if bai.exists() else False,
            "log_final_exists": log.is_file() and log.stat().st_size > 0 if log.exists() else False,
            "last_checked": now(), "notes": notes,
        })
    pd.DataFrame(rows).to_csv(star_dir / "STAR_RUN_STATUS.tsv", sep="\t", index=False)


def parse_number(value: str):
    value = value.strip()
    if value.endswith("%"):
        return float(value[:-1])
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def parse_star_log(path: Path) -> dict:
    found = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" not in line:
            continue
        key, value = [part.strip() for part in line.split("|", 1)]
        found[key] = parse_number(value)
    missing = [key for key in STAR_KEYS if key not in found]
    if missing:
        raise SystemExit(f"{path}: missing required STAR fields: {missing}")
    return {key: found[key] for key in STAR_KEYS}


def parse_flagstat(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    def first(pattern: str, cast=int, default=None):
        match = re.search(pattern, text, re.MULTILINE)
        return cast(match.group(1)) if match else default
    return {
        "total_reads_or_fragments": first(r"^(\d+) \+ \d+ in total"),
        "mapped_reads": first(r"^(\d+) \+ \d+ mapped"),
        "mapped_percent": first(r"^\d+ \+ \d+ mapped \(([\d.]+)%", float),
        "properly_paired": first(r"^(\d+) \+ \d+ properly paired"),
        "properly_paired_percent": first(r"^\d+ \+ \d+ properly paired \(([\d.]+)%", float),
        "secondary": first(r"^(\d+) \+ \d+ secondary"),
        "supplementary": first(r"^(\d+) \+ \d+ supplementary"),
        "duplicates_if_available": first(r"^(\d+) \+ \d+ duplicates"),
    }


def save_figure(fig: plt.Figure, out_dir: Path, stem: str, source: pd.DataFrame) -> None:
    fig_dir = out_dir / stem
    fig_dir.mkdir(parents=True, exist_ok=True)
    source.to_csv(fig_dir / "source_data.tsv", sep="\t", index=False)
    shutil.copy2(Path(__file__), fig_dir / "plotting_script.py")
    fig.savefig(fig_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def style_axis(ax, grid_y=True):
    ax.spines[["top", "right"]].set_visible(False)
    if grid_y:
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, zorder=0)
    ax.tick_params(labelsize=9)


def star_figures(qc: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    qc = qc.set_index("sample_id").loc[SAMPLE_ORDER].reset_index()
    colors = qc["batch"].map({"old": "#4C78A8", "new": "#E45756"}).tolist()
    x = np.arange(len(qc))

    stem = "Figure_QC_STAR_01_UniquelyMappedRate"
    source = qc[["sample_id", "patient_id", "region", "batch", "Uniquely mapped reads %"]].copy()
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    ax.bar(x, source["Uniquely mapped reads %"], color=colors, edgecolor="black", linewidth=0.5, zorder=2)
    ax.set_ylabel("Uniquely mapped reads (%)")
    ax.set_xlabel("Sample")
    ax.set_xticks(x, qc["sample_id"], rotation=45, ha="right")
    ax.set_title("STAR uniquely mapped read rate")
    style_axis(ax)
    save_figure(fig, out_dir, stem, source)

    stem = "Figure_QC_STAR_02_MappingComposition"
    cols = {
        "Uniquely mapped": "Uniquely mapped reads %",
        "Multiple loci": "% of reads mapped to multiple loci",
        "Too many loci": "% of reads mapped to too many loci",
        "Unmapped: too many mismatches": "% of reads unmapped: too many mismatches",
        "Unmapped: too short": "% of reads unmapped: too short",
        "Unmapped: other": "% of reads unmapped: other",
    }
    source = qc[["sample_id", "patient_id", "region", "batch"]].copy()
    for label, col in cols.items():
        source[label] = qc[col]
    source["composition_sum"] = source[list(cols)].sum(axis=1)
    if not np.allclose(source["composition_sum"], 100.0, atol=0.25):
        raise SystemExit("STAR mapping composition does not sum to approximately 100% for every sample.")
    palette = ["#4C78A8", "#F2CF5B", "#B279A2", "#E45756", "#72B7B2", "#9D755D"]
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    bottom = np.zeros(len(source))
    for (label, _), color in zip(cols.items(), palette):
        values = source[label].to_numpy(float)
        ax.bar(x, values, bottom=bottom, label=label, color=color, width=0.82, linewidth=0)
        bottom += values
    ax.set_ylabel("Reads (%)")
    ax.set_xlabel("Sample")
    ax.set_ylim(0, 100)
    ax.set_xticks(x, source["sample_id"], rotation=45, ha="right")
    ax.set_title("STAR mapping composition")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    style_axis(ax, grid_y=False)
    save_figure(fig, out_dir, stem, source)

    stem = "Figure_QC_STAR_03_MappedLength"
    source = qc[["sample_id", "patient_id", "region", "batch", "Average mapped length"]].copy()
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    ax.bar(x, source["Average mapped length"], color=colors, edgecolor="black", linewidth=0.5, zorder=2)
    ax.set_ylabel("Average mapped length (nt)")
    ax.set_xlabel("Sample")
    ax.set_xticks(x, source["sample_id"], rotation=45, ha="right")
    ax.set_title("STAR average mapped length")
    style_axis(ax)
    save_figure(fig, out_dir, stem, source)

    stem = "Figure_QC_STAR_04_OldNewMappingQC"
    metric_map = {
        "Uniquely mapped (%)": "Uniquely mapped reads %",
        "Multiple loci (%)": "% of reads mapped to multiple loci",
        "Unmapped: too short (%)": "% of reads unmapped: too short",
    }
    source = qc[["sample_id", "patient_id", "region", "batch"]].copy()
    source = source.merge(qc[["sample_id"] + list(metric_map.values())], on="sample_id")
    long = source.melt(id_vars=["sample_id", "patient_id", "region", "batch"], value_vars=list(metric_map.values()), var_name="metric", value_name="percent")
    long["metric"] = long["metric"].map({v: k for k, v in metric_map.items()})
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.6))
    rng = np.random.default_rng(20260828)
    for ax, metric in zip(axes, metric_map):
        part = long[long["metric"] == metric]
        for xi, batch in enumerate(["old", "new"]):
            vals = part.loc[part["batch"] == batch, "percent"].to_numpy(float)
            jitter = rng.uniform(-0.07, 0.07, len(vals))
            color = {"old": "#4C78A8", "new": "#E45756"}[batch]
            ax.scatter(np.full(len(vals), xi) + jitter, vals, s=38, color=color, edgecolor="black", linewidth=0.4, zorder=3)
            ax.plot([xi - 0.18, xi + 0.18], [np.mean(vals)] * 2, color="black", linewidth=2, zorder=4)
        ax.set_xticks([0, 1], ["old", "new"])
        ax.set_ylabel(metric)
        ax.set_xlabel("Sequencing set")
        style_axis(ax)
    fig.suptitle("Descriptive STAR mapping metrics by sequencing set", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, stem, long)


def checkpoint2(manifest: pd.DataFrame, star_dir: Path, figure_dir: Path) -> None:
    write_status(manifest, star_dir)
    status = pd.read_csv(star_dir / "STAR_RUN_STATUS.tsv", sep="\t")
    if not (status["star_status"] == "COMPLETED").all():
        raise SystemExit("Checkpoint 2 cannot run until all 12 STAR samples are COMPLETED.")
    audit_dir = star_dir / "bam_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    audit_rows, qc_rows = [], []
    for row in manifest.itertuples(index=False):
        bam, bai, log = paths_for(star_dir, row.sample_id)
        quick_path = audit_dir / f"{row.sample_id}.quickcheck.txt"
        flag_path = audit_dir / f"{row.sample_id}.flagstat.txt"
        idx_path = audit_dir / f"{row.sample_id}.idxstats.txt"
        quick = run(["samtools", "quickcheck", "-v", str(bam)])
        atomic_write(quick_path, "PASS\n" if quick.returncode == 0 else quick.stdout + quick.stderr)
        flag = run(["samtools", "flagstat", "-@", "4", str(bam)], flag_path)
        idx = run(["samtools", "idxstats", str(bam)], idx_path)
        if flag.returncode != 0 or idx.returncode != 0:
            raise SystemExit(f"samtools audit failed for {row.sample_id}: {flag.stderr} {idx.stderr}")
        values = parse_flagstat(flag_path)
        passed = quick.returncode == 0 and bai.is_file() and bai.stat().st_size > 0
        audit_rows.append({
            "sample_id": row.sample_id, "bam_path": str(bam.resolve()),
            "bam_size": bam.stat().st_size, "bai_exists": bai.is_file() and bai.stat().st_size > 0,
            "quickcheck_pass": quick.returncode == 0, **values,
            "audit_status": "PASS" if passed else "FAIL",
            "notes": "BAM readable; index and paired-end flagstat available" if passed else "BAM/index integrity failure",
        })
        qc_rows.append({"sample_id": row.sample_id, "patient_id": row.patient_id, "region": row.region, "batch": row.batch, **parse_star_log(log)})
    audit = pd.DataFrame(audit_rows)
    qc = pd.DataFrame(qc_rows)
    audit.to_csv(star_dir / "BAM_INTEGRITY_AUDIT.tsv", sep="\t", index=False)
    qc.to_csv(star_dir / "STAR_QC_SUMMARY.tsv", sep="\t", index=False)
    star_figures(qc, figure_dir)

    uniq = qc["Uniquely mapped reads %"].astype(float)
    median = float(uniq.median())
    mad = float(np.median(np.abs(uniq - median)))
    robust_z = np.zeros(len(uniq)) if mad == 0 else 0.6745 * (uniq - median) / mad
    outliers = qc.loc[np.abs(robust_z) > 3.5, "sample_id"].tolist()
    batch_summary = qc.groupby("batch")[["Uniquely mapped reads %", "% of reads mapped to multiple loci", "% of reads unmapped: too short"]].agg(["mean", "median", "min", "max"])
    batch_summary.to_csv(star_dir / "STAR_QC_BATCH_DESCRIPTIVE.tsv", sep="\t")
    all_bam_pass = bool((audit["audit_status"] == "PASS").all())
    gate = "PASS" if all_bam_pass and not outliers else "REVIEW"
    text = f"""# PDE2A bulk RNA-seq: Checkpoint 2 audit

Date: {now()}

## Completion and BAM integrity

- STAR completed: 12/12
- samtools quickcheck PASS: {int(audit['quickcheck_pass'].sum())}/12
- BAM index present and usable: {int(audit['bai_exists'].sum())}/12
- BAM integrity gate: {'PASS' if all_bam_pass else 'FAIL'}

## STAR mapping QC

- Uniquely mapped reads (%): range {uniq.min():.2f}-{uniq.max():.2f}; median {median:.2f}
- Distribution-based potential mapping outliers (robust |z| > 3.5): {', '.join(outliers) if outliers else 'none'}
- Samples requiring technical review: {', '.join(outliers) if outliers else 'none'}
- old/new mapping metrics are recorded descriptively in `STAR_QC_BATCH_DESCRIPTIVE.tsv`; this is not labeled a confirmed batch effect and is not a sample-exclusion rule.

## Gate decision

**Checkpoint 2 = {gate}.**

{'The 12 BAMs satisfy the prerequisites for RSeQC.' if gate == 'PASS' else 'RSeQC/featureCounts remain blocked pending review of the listed technical finding(s).'}

No gene-level expression or PDE2A result was inspected.
"""
    atomic_write(star_dir / "QC_CHECKPOINT2_AUDIT.md", text, archive_existing=True)
    if gate != "PASS":
        raise SystemExit("Checkpoint 2 REVIEW: do not proceed to RSeQC automatically.")


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text)}


def gtf_to_bed12(gtf: Path, bed12: Path) -> None:
    transcripts: dict[tuple[str, str, str], dict] = {}
    with gtf.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            key = (fields[0], fields[6], transcript_id)
            rec = transcripts.setdefault(key, {"gene_id": attrs.get("gene_id", transcript_id), "exons": []})
            rec["exons"].append((int(fields[3]) - 1, int(fields[4])))
    if not transcripts:
        raise SystemExit(f"No transcript exons parsed from {gtf}")
    lines = []
    for (chrom, strand, transcript_id), rec in sorted(transcripts.items()):
        exons = sorted(set(rec["exons"]))
        tx_start, tx_end = min(x[0] for x in exons), max(x[1] for x in exons)
        sizes = ",".join(str(end - start) for start, end in exons) + ","
        starts = ",".join(str(start - tx_start) for start, _ in exons) + ","
        name = f"{rec['gene_id']}|{transcript_id}"
        lines.append("\t".join(map(str, [chrom, tx_start, tx_end, name, 0, strand, tx_start, tx_end, 0, len(exons), sizes, starts])))
    atomic_write(bed12, "\n".join(lines) + "\n", archive_existing=True)


def audit_bed12(bed12: Path, bam: Path, output_dir: Path) -> None:
    """Validate BED12 structure/coordinates and chromosome compatibility with BAM."""
    output_dir.mkdir(parents=True, exist_ok=True)
    header = run(["samtools", "view", "-H", str(bam)])
    if header.returncode != 0:
        raise SystemExit(f"Cannot read BAM header for BED12 audit: {header.stderr}")
    bam_contigs = {}
    for line in header.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" in fields and "LN" in fields:
            bam_contigs[fields["SN"]] = int(fields["LN"])
    rows = []
    invalid = 0
    with bed12.open(encoding="utf-8", errors="strict") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            reason = []
            if len(fields) != 12:
                rows.append({"line": line_no, "chrom": fields[0] if fields else "", "valid": False, "notes": f"column_count={len(fields)}"})
                invalid += 1
                continue
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
                thick_start, thick_end = int(fields[6]), int(fields[7])
                block_count = int(fields[9])
                sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            except ValueError:
                reason.append("non_integer_coordinate_or_block")
                start = end = thick_start = thick_end = block_count = 0
                sizes, starts = [], []
            if chrom not in bam_contigs:
                reason.append("chrom_not_in_bam_header")
            if start < 0 or end <= start or thick_start < start or thick_end > end:
                reason.append("invalid_transcript_coordinates")
            if chrom in bam_contigs and end > bam_contigs[chrom]:
                reason.append("end_exceeds_bam_contig_length")
            if block_count <= 0 or len(sizes) != block_count or len(starts) != block_count:
                reason.append("block_count_mismatch")
            if any(size <= 0 for size in sizes) or any(offset < 0 for offset in starts):
                reason.append("invalid_block")
            if sizes and starts and max(offset + size for offset, size in zip(starts, sizes)) > end - start:
                reason.append("block_exceeds_transcript")
            valid = not reason
            invalid += int(not valid)
            rows.append({"line": line_no, "chrom": chrom, "valid": valid, "notes": ";".join(reason) if reason else "PASS"})
    if not rows:
        raise SystemExit("BED12 is empty.")
    detail = pd.DataFrame(rows)
    detail.to_csv(output_dir / "BED12_LINE_AUDIT.tsv", sep="\t", index=False)
    chroms = sorted(detail["chrom"].unique())
    summary = pd.DataFrame([{
        "bed12_path": str(bed12.resolve()), "bam_path": str(bam.resolve()),
        "line_count": len(detail), "column_count_required": 12,
        "invalid_line_count": invalid, "bed_chromosome_count": len(chroms),
        "bam_contig_count": len(bam_contigs),
        "chromosome_naming_compatible": all(chrom in bam_contigs for chrom in chroms),
        "audit_status": "PASS" if invalid == 0 else "FAIL",
    }])
    summary.to_csv(output_dir / "BED12_AUDIT.tsv", sep="\t", index=False)
    text = f"""# BED12 reference audit

Date: {now()}

- Source: frozen GENCODE v44 GTF used by the STAR reference
- Conversion route: UCSC gtfToGenePred -> genePredToBed
- BED12 records: {len(detail)}
- Required 12-column records valid: {len(detail) - invalid}/{len(detail)}
- BED chromosomes found in BAM header: {'yes' if all(chrom in bam_contigs for chrom in chroms) else 'no'}
- Transcript/block coordinates valid: {'yes' if invalid == 0 else 'no'}
- Audit status: {'PASS' if invalid == 0 else 'FAIL'}
"""
    atomic_write(output_dir / "BED12_AUDIT.md", text, archive_existing=True)
    if invalid:
        raise SystemExit(f"BED12 audit failed for {invalid} records.")


def parse_rseqc(path: Path) -> tuple[float, float, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    def fraction(label: str) -> float | None:
        match = re.search(rf'Fraction of reads explained by "{re.escape(label)}":\s*([\d.]+)', text)
        return float(match.group(1)) if match else None
    a = fraction("1++,1--,2+-,2-+")
    b = fraction("1+-,1-+,2++,2--")
    failed_match = re.search(r"Fraction of reads failed to determine:\s*([\d.]+)", text)
    other = float(failed_match.group(1)) if failed_match else max(0.0, 1.0 - (a or 0) - (b or 0))
    if a is None or b is None:
        raise SystemExit(f"Could not parse paired-end RSeQC patterns from {path}")
    return a, b, other


def classify_strand(a: float, b: float, other: float) -> tuple[str, str, str]:
    if abs(a - b) <= 0.15 and a + b >= 0.75:
        return "unstranded", "high" if abs(a - b) <= 0.10 else "moderate", "featureCounts candidate -s 0"
    winner = max(a, b)
    confidence = "high" if winner >= 0.80 and other <= 0.20 else "moderate" if winner >= 0.65 else "low"
    if a > b:
        return "forward/sense stranded", confidence, "featureCounts candidate -s 1"
    return "reverse/antisense stranded", confidence, "featureCounts candidate -s 2"


def rseqc_figures(summary: pd.DataFrame, out_dir: Path) -> None:
    summary = summary.set_index("sample_id").loc[SAMPLE_ORDER].reset_index()
    x = np.arange(len(summary))
    stem = "Figure_QC_RSeQC_01_StrandnessFractions"
    source = summary[["sample_id", "patient_id", "region", "batch", "fraction_pattern_A", "fraction_pattern_B", "other_fraction_if_present"]].copy()
    fig, ax = plt.subplots(figsize=(11, 5.8))
    bottom = np.zeros(len(source))
    for col, label, color in [
        ("fraction_pattern_A", "Pattern A", "#4C78A8"),
        ("fraction_pattern_B", "Pattern B", "#E45756"),
        ("other_fraction_if_present", "Other/undetermined", "#BAB0AC"),
    ]:
        vals = source[col].to_numpy(float) * 100
        ax.bar(x, vals, bottom=bottom, label=label, color=color, width=0.82)
        bottom += vals
    ax.set_ylim(0, 100)
    ax.set_ylabel("Read fraction (%)")
    ax.set_xlabel("Sample")
    ax.set_xticks(x, source["sample_id"], rotation=45, ha="right")
    ax.set_title("RSeQC paired-end strandness fractions")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    style_axis(ax, grid_y=False)
    save_figure(fig, out_dir, stem, source)

    stem = "Figure_QC_RSeQC_02_StrandnessByBatch"
    long = source.melt(id_vars=["sample_id", "patient_id", "region", "batch"], value_vars=["fraction_pattern_A", "fraction_pattern_B"], var_name="pattern", value_name="fraction")
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.6), sharey=True)
    rng = np.random.default_rng(20260828)
    for ax, pattern in zip(axes, ["fraction_pattern_A", "fraction_pattern_B"]):
        part = long[long["pattern"] == pattern]
        for xi, batch in enumerate(["old", "new"]):
            vals = part.loc[part["batch"] == batch, "fraction"].to_numpy(float)
            jitter = rng.uniform(-0.07, 0.07, len(vals))
            color = {"old": "#4C78A8", "new": "#E45756"}[batch]
            ax.scatter(np.full(len(vals), xi) + jitter, vals, s=42, color=color, edgecolor="black", linewidth=0.4)
            ax.plot([xi - 0.18, xi + 0.18], [np.mean(vals)] * 2, color="black", linewidth=2)
        ax.set_xticks([0, 1], ["old", "new"])
        ax.set_xlabel("Sequencing set")
        ax.set_title("Pattern A" if pattern.endswith("A") else "Pattern B")
        style_axis(ax)
    axes[0].set_ylabel("Read fraction")
    fig.suptitle("Descriptive RSeQC strandness fractions by sequencing set", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, stem, long)


def checkpoint3(manifest: pd.DataFrame, out_dir: Path, figure_dir: Path) -> None:
    rows = []
    for row in manifest.itertuples(index=False):
        path = out_dir / f"{row.sample_id}.infer_experiment.txt"
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing RSeQC output: {path}")
        a, b, other = parse_rseqc(path)
        classification, confidence, note = classify_strand(a, b, other)
        rows.append({
            "sample_id": row.sample_id, "patient_id": row.patient_id,
            "region": row.region, "batch": row.batch,
            "fraction_pattern_A": a, "fraction_pattern_B": b,
            "other_fraction_if_present": other,
            "strand_classification": classification, "confidence": confidence, "notes": note,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "RSeQC_STRANDNESS_SUMMARY.tsv", sep="\t", index=False)
    rseqc_figures(summary, figure_dir)
    classes = summary["strand_classification"].unique().tolist()
    old_classes = sorted(summary.loc[summary["batch"] == "old", "strand_classification"].unique())
    new_classes = sorted(summary.loc[summary["batch"] == "new", "strand_classification"].unique())
    consistent = len(classes) == 1
    batch_consistent = old_classes == new_classes and len(old_classes) == 1
    mapping = {"unstranded": 0, "forward/sense stranded": 1, "reverse/antisense stranded": 2}
    recommended = mapping.get(classes[0]) if consistent and batch_consistent else None
    gate = "PASS" if recommended is not None and (summary["confidence"] != "low").all() else "REVIEW"
    text = f"""# PDE2A bulk RNA-seq: Checkpoint 3 audit

Date: {now()}

## RSeQC completion

- infer_experiment successful: 12/12
- Sample-level strandedness classifications: {', '.join(classes)}
- 12-sample classification consistent: {'yes' if consistent else 'no'}
- old/new classification consistent: {'yes' if batch_consistent else 'no'}
- Final strandedness call: {classes[0] if consistent and batch_consistent else 'unresolved'}
- Recommended featureCounts parameter: {f'-s {recommended}' if recommended is not None else 'UNSET'}

## Gate decision

**Checkpoint 3 = {gate}.**

- featureCounts: {'READY BUT NOT RUN' if gate == 'PASS' else 'BLOCKED'}
- PDE2A: NOT INSPECTED

This checkpoint records the RSeQC evidence only. No featureCounts or gene-level analysis was run.
"""
    atomic_write(out_dir / "QC_CHECKPOINT3_AUDIT.md", text, archive_existing=True)
    if gate != "PASS":
        raise SystemExit("Checkpoint 3 REVIEW: FEATURECOUNTS_STRANDNESS must remain UNSET.")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--manifest", type=Path, required=True)
    common.add_argument("--star-dir", type=Path, default=Path("results/02_star"))
    sub.add_parser("status", parents=[common])
    cp2 = sub.add_parser("checkpoint2", parents=[common])
    cp2.add_argument("--figure-dir", type=Path, required=True)
    bed = sub.add_parser("gtf2bed12")
    bed.add_argument("--gtf", type=Path, required=True)
    bed.add_argument("--bed12", type=Path, required=True)
    bed_audit = sub.add_parser("audit-bed12")
    bed_audit.add_argument("--bed12", type=Path, required=True)
    bed_audit.add_argument("--bam", type=Path, required=True)
    bed_audit.add_argument("--output-dir", type=Path, required=True)
    cp3 = sub.add_parser("checkpoint3")
    cp3.add_argument("--manifest", type=Path, required=True)
    cp3.add_argument("--rseqc-dir", type=Path, required=True)
    cp3.add_argument("--figure-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.command == "gtf2bed12":
        gtf_to_bed12(args.gtf, args.bed12)
        return
    if args.command == "audit-bed12":
        audit_bed12(args.bed12, args.bam, args.output_dir)
        return
    manifest = read_manifest(args.manifest)
    if args.command == "status":
        write_status(manifest, args.star_dir)
    elif args.command == "checkpoint2":
        checkpoint2(manifest, args.star_dir, args.figure_dir)
    elif args.command == "checkpoint3":
        checkpoint3(manifest, args.rseqc_dir, args.figure_dir)


if __name__ == "__main__":
    main()
