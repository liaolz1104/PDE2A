#!/usr/bin/env python3
"""Reconcile prior PDE2A values with the frozen pre-analysis audit outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


PATIENTS = ["P08", "P09", "P34", "P36", "P38"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("vendor_expression", type=Path)
    args = parser.parse_args()

    root = args.project_root.resolve()
    out = root / "results" / "05_pde2a_audit"
    out.mkdir(parents=True, exist_ok=True)
    required = [
        root / "manifest_combined_standardized.tsv",
        out / "PDE2A_patientwise_direction_audit.tsv",
        args.vendor_expression,
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty required input: {path}")

    manifest = pd.read_csv(required[0], sep="\t", dtype=str)
    current = pd.read_csv(required[1], sep="\t")
    vendor = pd.read_csv(args.vendor_expression, sep="\t", low_memory=False)
    prior_gene = vendor.loc[vendor["gene_symbol"].eq("PDE2A")]
    if len(prior_gene) != 1:
        raise RuntimeError(f"Prior vendor PDE2A row is not unique: {len(prior_gene)}")
    prior_gene = prior_gene.iloc[0]

    prior = {
        "P08": (float(prior_gene["fpkm_P8_T_C"]), float(prior_gene["fpkm_P8_T_I"])),
        "P09": (float(prior_gene["fpkm_P9_T_C"]), float(prior_gene["fpkm_P9_T_I"])),
    }
    rows = []
    source = str(args.vendor_expression.resolve())
    for patient in PATIENTS:
        cur = current.loc[current["patient_id"].eq(patient)]
        if len(cur) != 1:
            raise RuntimeError(f"Current direction row is not unique for {patient}")
        cur = cur.iloc[0]
        if patient in prior:
            previous_core, previous_edge = prior[patient]
            previous_delta = previous_edge - previous_core
            previous_direction = "Edge<Core" if previous_edge < previous_core else ("Edge>Core" if previous_edge > previous_core else "Equal")
            current_direction = "Edge<Core" if cur.edge_normalized < cur.core_normalized else ("Edge>Core" if cur.edge_normalized > cur.core_normalized else "Equal")
            concordant = "CONCORDANT" if previous_direction == current_direction else "REVIEW"
            previous_unit = "FPKM"
            previous_source_file = source
            previous_method = "Vendor gene_expression.xls FPKM; exact upstream reference/annotation/strandedness not documented in the recovered source table"
            notes = "Exact P8/P9 source row and sample columns recovered; labels map T_C to Core and T_I to Edge."
        else:
            previous_core = previous_edge = previous_delta = pd.NA
            previous_direction = previous_unit = "SOURCE_NOT_FOUND"
            previous_source_file = "SOURCE_NOT_FOUND"
            previous_method = "SOURCE_NOT_FOUND"
            current_direction = "Edge<Core" if cur.edge_normalized < cur.core_normalized else ("Edge>Core" if cur.edge_normalized > cur.core_normalized else "Equal")
            concordant = "SOURCE_NOT_FOUND"
            notes = "No numeric source data or plotting script for the prior P34/P36/P38 expression figure was recovered; image values were not inferred."
        rows.append({
            "patient_id": patient,
            "previous_core": previous_core,
            "previous_edge": previous_edge,
            "previous_delta": previous_delta,
            "previous_direction": previous_direction,
            "previous_unit": previous_unit,
            "previous_source_file": previous_source_file,
            "previous_method": previous_method,
            "current_raw_core": int(cur.core_raw),
            "current_raw_edge": int(cur.edge_raw),
            "current_normalized_core": float(cur.core_normalized),
            "current_normalized_edge": float(cur.edge_normalized),
            "current_delta": float(cur.delta_normalized_edge_minus_core),
            "current_direction": current_direction,
            "direction_concordant": concordant,
            "notes": notes,
        })
    pd.DataFrame(rows).to_csv(out / "PDE2A_previous_vs_current_reconciliation.tsv", sep="\t", index=False, na_rep="NA")

    selected = manifest.loc[manifest["sample_id"].isin([f"{p}_{r}" for p in PATIENTS for r in ("C", "TI")])].copy()
    selected["mapping_interpretation"] = selected["region"].map({"c": "Core", "ti": "Edge/TI"})
    selected[[
        "sample_id", "original_sample_id", "patient_id", "region", "mapping_interpretation", "batch",
        "read1_windows", "read2_windows", "include_main_edge_core",
    ]].to_csv(out / "PDE2A_sample_selection_provenance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
