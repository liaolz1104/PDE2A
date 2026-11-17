"""Query L1000CDS2 for compounds predicted to reverse a PDE2A-loss signature."""

import argparse
import json
import urllib.request
from pathlib import Path

import pandas as pd


def summarize_l1000_results(top_meta: list[dict]) -> pd.DataFrame:
    """Aggregate L1000 conditions by compound without discarding condition-level evidence."""
    rows = []
    for item in top_meta:
        overlap = item.get("overlap", {})
        rows.append(
            {
                "compound": item.get("pert_desc"),
                "score": item.get("score"),
                "cell_line": item.get("cell_id"),
                "reverse_overlap_genes": sorted(
                    set(overlap.get("up/dn", []) + overlap.get("dn/up", []))
                ),
            }
        )
    conditions = pd.DataFrame(rows)
    grouped = conditions.groupby("compound", dropna=False, sort=False).agg(
        n_lincs_conditions=("score", "size"),
        best_reverse_score=("score", "max"),
        mean_reverse_score=("score", "mean"),
        cell_lines=("cell_line", lambda values: ";".join(sorted(set(values)))),
        reverse_overlap_genes=(
            "reverse_overlap_genes",
            lambda values: ";".join(sorted({gene for genes in values for gene in genes})),
        ),
    )
    return grouped.reset_index().sort_values(
        ["best_reverse_score", "n_lincs_conditions"], ascending=False
    ).reset_index(drop=True)


def query_l1000cds2(loss_up: list[str], loss_down: list[str]) -> dict:
    """Submit a reverse gene-set query to the public L1000CDS2 API."""
    payload = {
        "data": {"upGenes": loss_up, "dnGenes": loss_down},
        "config": {
            "aggravate": False,
            "searchMethod": "geneSet",
            "share": False,
            "combination": False,
            "db-version": "latest",
        },
        "metadata": [{"key": "tag", "value": "HN6_PDE2A_loss_bidirectional_core"}],
    }
    request = urllib.request.Request(
        "https://maayanlab.cloud/L1000CDS2/query",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an L1000CDS2 reverse query.")
    parser.add_argument("signature_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    loss_up = (args.signature_dir / "pde2a_loss_up.grp").read_text(encoding="utf-8").split()
    loss_down = (args.signature_dir / "pde2a_loss_down.grp").read_text(encoding="utf-8").split()
    result = query_l1000cds2(loss_up, loss_down)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "l1000cds2_reverse_raw.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame(result["topMeta"]).to_json(
        args.output_dir / "l1000cds2_reverse_conditions.json", orient="records", indent=2
    )
    summarize_l1000_results(result["topMeta"]).to_csv(
        args.output_dir / "l1000cds2_reverse_compounds.csv", index=False
    )
    print(f"Saved {len(result['topMeta'])} L1000 conditions to {args.output_dir}")


if __name__ == "__main__":
    main()
