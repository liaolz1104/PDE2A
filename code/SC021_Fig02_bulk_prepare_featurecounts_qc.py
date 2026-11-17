#!/usr/bin/env python3
"""Validate featureCounts output and freeze raw matrices for Checkpoint 4 QC."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLE_ORDER = [
    "P08_C", "P08_TI", "P08_P", "P09_C", "P09_TI", "P09_P",
    "P34_C", "P34_TI", "P36_C", "P36_TI", "P38_C", "P38_TI",
]


def read_manifest(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str)
    required = {"sample_id", "patient_id", "region", "batch"}
    missing = required - set(meta.columns)
    if missing:
        raise SystemExit(f"Manifest missing required fields: {sorted(missing)}")
    if meta["sample_id"].tolist() != SAMPLE_ORDER:
        raise SystemExit("Manifest does not match the frozen 12-sample order.")
    if meta["sample_id"].duplicated().any():
        raise SystemExit("Duplicate sample_id in manifest.")
    return meta


def sample_from_bam_column(column: str) -> str:
    hits = [sample for sample in SAMPLE_ORDER if f"/{sample}.Aligned.sortedByCoord.out.bam" in column.replace("\\", "/")]
    if len(hits) != 1:
        raise SystemExit(f"Cannot map featureCounts BAM column uniquely: {column}")
    return hits[0]


def parse_gtf_annotation(gtf: Path) -> pd.DataFrame:
    records = []
    seen = set()
    with gtf.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            if not gid or gid.group(1) in seen:
                continue
            gene_id = gid.group(1)
            seen.add(gene_id)
            name = re.search(r'gene_name "([^"]+)"', attrs)
            gene_type = re.search(r'gene_type "([^"]+)"', attrs) or re.search(r'gene_biotype "([^"]+)"', attrs)
            records.append({
                "gene_id": gene_id,
                "gene_id_stripped": re.sub(r"\.\d+$", "", gene_id),
                "gene_name": name.group(1) if name else gene_id,
                "gene_type": gene_type.group(1) if gene_type else "",
            })
    if not records:
        raise SystemExit(f"No gene records parsed from frozen GTF: {gtf}")
    return pd.DataFrame(records)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--featurecounts", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--gtf", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    for path in (args.featurecounts, args.summary, args.manifest, args.gtf):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Required input missing or empty: {path}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    meta = read_manifest(args.manifest)

    fc = pd.read_csv(args.featurecounts, sep="\t", comment="#", dtype=str)
    if "Geneid" not in fc.columns or len(fc.columns) < 7:
        raise SystemExit("Malformed featureCounts output.")
    bam_cols = list(fc.columns[6:])
    if len(bam_cols) != 12 or len(set(bam_cols)) != 12:
        raise SystemExit(f"Expected 12 unique BAM count columns; found {len(bam_cols)}.")
    mapped_samples = [sample_from_bam_column(col) for col in bam_cols]
    if len(set(mapped_samples)) != 12 or set(mapped_samples) != set(SAMPLE_ORDER):
        raise SystemExit("featureCounts columns do not map one-to-one to the frozen manifest.")

    numeric = fc[bam_cols].apply(pd.to_numeric, errors="raise")
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise SystemExit("Counts must be finite non-negative integers.")
    numeric = numeric.astype(np.int64)
    rename = dict(zip(bam_cols, mapped_samples))
    matrix = pd.concat([fc[["Geneid"]].rename(columns={"Geneid": "gene_id"}), numeric.rename(columns=rename)], axis=1)
    if matrix["gene_id"].duplicated().any():
        raise SystemExit("Duplicate gene_id rows in featureCounts output.")
    matrix = matrix[["gene_id"] + SAMPLE_ORDER]
    matrix.to_csv(args.outdir / "count_matrix_gene_id.tsv", sep="\t", index=False)

    annotation = parse_gtf_annotation(args.gtf)
    annotation.to_csv(args.outdir / "gene_annotation.tsv", sep="\t", index=False)
    missing_annotation = sorted(set(matrix["gene_id"]) - set(annotation["gene_id"]))
    if missing_annotation:
        raise SystemExit(f"Count matrix has {len(missing_annotation)} gene IDs absent from frozen GTF gene records.")

    summary = pd.read_csv(args.summary, sep="\t", dtype=str)
    if "Status" not in summary.columns:
        raise SystemExit("featureCounts summary lacks Status column.")
    summary_cols = list(summary.columns[1:])
    if len(summary_cols) != 12 or len(set(summary_cols)) != 12:
        raise SystemExit("featureCounts summary must contain 12 unique BAM columns.")
    summary_samples = [sample_from_bam_column(col) for col in summary_cols]
    summary_values = summary[summary_cols].apply(pd.to_numeric, errors="raise").astype(np.int64)
    by_status = summary_values.copy()
    by_status.insert(0, "Status", summary["Status"])
    by_status = by_status.set_index("Status")

    rows = []
    explicit = {
        "assigned": "Assigned",
        "unassigned_unmapped": "Unassigned_Unmapped",
        "unassigned_multimapping": "Unassigned_MultiMapping",
        "unassigned_no_features": "Unassigned_NoFeatures",
        "unassigned_ambiguity": "Unassigned_Ambiguity",
        "unassigned_chimera": "Unassigned_Chimera",
        "unassigned_mapping_quality": "Unassigned_MappingQuality",
    }
    for col, sample in zip(summary_cols, summary_samples):
        counts = by_status[col]
        total = int(counts.sum())
        assigned = int(counts.get("Assigned", 0))
        row = meta.loc[meta["sample_id"] == sample].iloc[0]
        record = {
            "sample_id": sample, "patient_id": row["patient_id"],
            "region": row["region"], "batch": row["batch"],
            "bam_path": col, "total_fragments_considered": total,
            "assigned": assigned,
            "assigned_percent": 100.0 * assigned / total if total else np.nan,
        }
        used = set()
        for output_name, status_name in explicit.items():
            record[output_name] = int(counts.get(status_name, 0))
            used.add(status_name)
        record["other_unassigned_if_present"] = int(sum(value for status, value in counts.items() if status != "Assigned" and status not in used))
        record["audit_status"] = "PASS" if total > 0 else "FAIL"
        record["notes"] = "12-column mapping verified; non-negative integer counts" if total > 0 else "No fragments summarized"
        rows.append(record)
    audit = pd.DataFrame(rows).set_index("sample_id").loc[SAMPLE_ORDER].reset_index()
    rates = audit["assigned_percent"].to_numpy(float)
    median = float(np.median(rates))
    mad = float(np.median(np.abs(rates - median)))
    robust_z = np.zeros(len(rates)) if mad == 0 else 0.6745 * (rates - median) / mad
    audit["assignment_robust_z"] = robust_z
    audit.loc[np.abs(robust_z) > 3.5, "audit_status"] = "REVIEW"
    audit.loc[np.abs(robust_z) > 3.5, "notes"] = "Distribution-based assignment-rate outlier; technical review required"
    audit.to_csv(args.outdir / "FEATURECOUNTS_RUN_AUDIT.tsv", sep="\t", index=False)

    metadata = meta[["sample_id", "patient_id", "region", "batch"]].copy()
    metadata.to_csv(args.outdir / "metadata_all12.tsv", sep="\t", index=False)
    print(f"Validated featureCounts: {len(matrix)} genes x 12 samples; annotation rows={len(annotation)}")


if __name__ == "__main__":
    main()
