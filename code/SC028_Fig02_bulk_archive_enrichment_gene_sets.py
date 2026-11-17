#!/usr/bin/env python3
"""Download and freeze the four prespecified standard Enrichr gene-set libraries."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen


LIBRARIES = {
    "Hallmark": "MSigDB_Hallmark_2020",
    "Reactome": "Reactome_2022",
    "GO_BP": "GO_Biological_Process_2023",
    "KEGG": "KEGG_2021_Human",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    out = args.project_root.resolve() / "results" / "06_enrichment" / "gene_sets"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for short_name, library_name in LIBRARIES.items():
        url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=" + quote(library_name)
        destination = out / f"{library_name}.gmt"
        if destination.exists() and destination.stat().st_size > 0:
            status = "REUSED_FROZEN"
        else:
            with urlopen(url, timeout=120) as response:
                payload = response.read()
            if not payload.strip():
                raise RuntimeError(f"Empty library download: {library_name}")
            destination.write_bytes(payload.rstrip(b"\r\n") + b"\n")
            status = "DOWNLOADED"
        # Normalize the terminal newline without changing any gene-set record.
        payload = destination.read_bytes()
        destination.write_bytes(payload.rstrip(b"\r\n") + b"\n")
        set_count = sum(1 for line in destination.read_text(encoding="utf-8").splitlines() if line.strip())
        rows.append({
            "collection": short_name,
            "library_name": library_name,
            "source_url": url,
            "frozen_file": destination.name,
            "gene_set_count": set_count,
            "download_status": status,
            "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        })
    with (out / "GENE_SET_SOURCE_MANIFEST.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
