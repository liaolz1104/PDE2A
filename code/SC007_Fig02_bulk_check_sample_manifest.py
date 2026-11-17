#!/usr/bin/env python3
import argparse, os, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--path-column", default="wsl", choices=["wsl", "windows"])
args = ap.parse_args()

df = pd.read_csv(args.manifest, sep="\t")
required = {
    "sample_id","patient_id","region","batch",
    "read1_windows","read2_windows","read1_wsl","read2_wsl"
}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Missing required columns: {sorted(missing)}")

if df["sample_id"].duplicated().any():
    raise SystemExit("Duplicate sample_id detected.")

allowed = {"c","ti","p"}
bad = sorted(set(df["region"]) - allowed)
if bad:
    raise SystemExit(f"Unexpected region labels: {bad}")

main = df[df["region"].isin(["c","ti"])].copy()
tab = main.groupby(["patient_id","region"]).size().unstack(fill_value=0)
expected_patients = ["P08","P09","P34","P36","P38"]
for p in expected_patients:
    if p not in tab.index:
        raise SystemExit(f"Missing main patient: {p}")
    if tab.loc[p].get("c",0) != 1 or tab.loc[p].get("ti",0) != 1:
        raise SystemExit(f"{p} does not have exactly one Core and one Edge sample.")

r1_col = "read1_wsl" if args.path_column == "wsl" else "read1_windows"
r2_col = "read2_wsl" if args.path_column == "wsl" else "read2_windows"
missing_paths = []
if args.path_column == "wsl":
    for _, r in df.iterrows():
        for col in (r1_col, r2_col):
            if not os.path.isfile(r[col]):
                missing_paths.append((r["sample_id"], col, r[col]))

print(f"Samples: {len(df)}")
print(f"Main paired Core/Edge samples: {len(main)}")
print(f"Paired patients: {main['patient_id'].nunique()}")
print("Regions:", df["region"].value_counts().to_dict())
print("Batches:", df["batch"].value_counts().to_dict())

if args.path_column == "wsl":
    if missing_paths:
        print("\nWARNING: paths not visible from this Linux/WSL environment:")
        for x in missing_paths:
            print("\t".join(map(str,x)))
        sys.exit(2)
    print("All WSL FASTQ paths exist.")

print("Manifest preflight passed.")
