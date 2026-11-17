#!/usr/bin/env python3
import argparse, gzip, re, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

opener = gzip.open if args.gtf.endswith(".gz") else open
seen = {}
with opener(args.gtf, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != 9 or f[2] != "gene":
            continue
        attrs = f[8]
        gid = re.search(r'gene_id "([^"]+)"', attrs)
        gname = re.search(r'gene_name "([^"]+)"', attrs)
        gtype = re.search(r'gene_type "([^"]+)"', attrs) or re.search(r'gene_biotype "([^"]+)"', attrs)
        if not gid:
            continue
        gene_id = gid.group(1).split(".")[0]
        seen[gene_id] = {
            "gene_id": gene_id,
            "gene_name": gname.group(1) if gname else gene_id,
            "gene_type": gtype.group(1) if gtype else "",
        }

pd.DataFrame(seen.values()).to_csv(args.output, sep="\t", index=False)
print(f"Wrote {len(seen)} genes to {args.output}")
