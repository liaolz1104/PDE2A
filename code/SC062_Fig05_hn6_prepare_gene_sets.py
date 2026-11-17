from __future__ import annotations

import csv
import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results" / "Edge_CSC_multialgorithm_source_author_20260805"
OUT = Path(__file__).resolve().parent / "frozen_gene_sets.gmt"

SEEDS = [
    "S100A2", "KRT17", "KRT14", "KRT5", "KRT6A", "TP63", "ITGA6",
    "CD44", "JAG1", "FABP5", "EGFR", "MYC", "KRT15", "NES", "BMI1",
]


def main() -> None:
    with gzip.open(SOURCE / "edge_CSC_high_low_gene_screen.tsv.gz", "rt") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        edge155 = sorted(row["gene"] for row in rows if row["passes_gate"] == "True")
    if len(edge155) != 155:
        raise RuntimeError(f"Expected 155 frozen Edge CSC genes, found {len(edge155)}")

    with (SOURCE / "consensus_coexpression_nodes.tsv").open() as handle:
        nodes = sorted(row["gene"] for row in csv.DictReader(handle, delimiter="\t"))
    if len(nodes) != 15:
        raise RuntimeError(f"Expected 15 strict network nodes, found {len(nodes)}")

    sets = {
        "EDGE_CSC_155": edge155,
        "STRICT_NETWORK_15": nodes,
        "ORIGINAL_CSC_SEEDS_15": SEEDS,
    }
    with OUT.open("w", newline="") as handle:
        for name, genes in sets.items():
            handle.write("\t".join([name, "frozen_2026-08-05", *genes]) + "\n")


if __name__ == "__main__":
    main()
