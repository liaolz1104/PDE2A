#!/usr/bin/env python3
import argparse, os, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--featurecounts", required=True)
ap.add_argument("--manifest", required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

meta = pd.read_csv(args.manifest, sep="\t")
fc = pd.read_csv(args.featurecounts, sep="\t", comment="#")
if "Geneid" not in fc.columns:
    raise SystemExit("featureCounts output missing Geneid column.")

bam_cols = list(fc.columns[6:])
if len(bam_cols) != len(meta):
    raise SystemExit(f"BAM columns ({len(bam_cols)}) != manifest samples ({len(meta)}).")

rename = {}
for col in bam_cols:
    matches = [s for s in meta["sample_id"] if s in col]
    if len(matches) != 1:
        # fallback: rely on manifest order, featureCounts preserves BAM argument order
        matches = []
    if matches:
        rename[col] = matches[0]

if len(rename) != len(bam_cols):
    rename = dict(zip(bam_cols, meta["sample_id"].tolist()))

mat = fc[["Geneid"] + bam_cols].rename(columns=rename)
mat["Geneid"] = mat["Geneid"].str.replace(r"\.\d+$", "", regex=True)
mat.to_csv(args.output, sep="\t", index=False)
print(f"Wrote count matrix: {mat.shape[0]} genes x {mat.shape[1]-1} samples")
