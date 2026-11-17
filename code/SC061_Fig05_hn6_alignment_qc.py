from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = {
    "Number of input reads": "input_read_pairs",
    "Uniquely mapped reads %": "uniquely_mapped_percent",
    "% of reads mapped to multiple loci": "multimapped_percent",
    "% of reads unmapped: too short": "unmapped_too_short_percent",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("star_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = []
    for path in sorted(args.star_dir.glob("*/Log.final.out")):
        values = {"sample": path.parent.name}
        for line in path.read_text().splitlines():
            if "|" not in line:
                continue
            key, value = [part.strip() for part in line.split("|", 1)]
            if key in FIELDS:
                values[FIELDS[key]] = float(re.sub(r"[%\s]", "", value))
        rows.append(values)
    if len(rows) != 12:
        raise RuntimeError(f"Expected 12 STAR final logs, found {len(rows)}")
    columns = ["sample", *FIELDS.values()]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
