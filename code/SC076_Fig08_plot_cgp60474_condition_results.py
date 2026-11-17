"""Plot condition-level L1000CDS2 evidence for the CGP-60474 shortlist."""

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "mechanism_drug_work" / "results"
STRICT = RESULTS / "l1000cds2_reverse" / "l1000cds2_reverse_conditions.json"
WIDE = RESULTS / "l1000cds2_reverse_wide" / "l1000cds2_reverse_wide_conditions.json"
CONSENSUS = RESULTS / "l1000cds2_consensus" / "strict_wide_consensus_compounds.csv"
TRIAGE = RESULTS / "l1000cds2_consensus" / "consensus_candidate_triage.csv"
DSIGDB = RESULTS / "web_drug_evidence" / "DSigDB_table.txt"
OUTDIR = ROOT / "results" / "CGP60474_L1000_condition_level_20261001"

PURPLE = "#7C3AED"
TEAL = "#0F766E"
NAVY = "#17324D"
SLATE = "#5C6E7E"
GRID = "#D8E0E8"
HIGHLIGHT = "#F5EFFF"


def normalize_name(value):
    return value.strip().lower()


def load_records():
    with CONSENSUS.open(encoding="utf-8-sig", newline="") as handle:
        consensus = list(csv.DictReader(handle))
    shared = {normalize_name(row["compound_strict"]): row["compound_strict"] for row in consensus}

    with TRIAGE.open(encoding="utf-8-sig", newline="") as handle:
        triage = {normalize_name(row["compound"]): row for row in csv.DictReader(handle)}

    records = []
    for query, path in (("Strict", STRICT), ("Wide", WIDE)):
        conditions = json.loads(path.read_text(encoding="utf-8"))
        for condition in conditions:
            key = normalize_name(condition["pert_desc"])
            if key not in shared:
                continue
            overlap = condition.get("overlap", {})
            up_dn = overlap.get("up/dn", [])
            dn_up = overlap.get("dn/up", [])
            records.append({
                "compound": shared[key],
                "query": query,
                "score": float(condition["score"]),
                "cell_line": str(condition["cell_id"]),
                "dose": str(condition["pert_dose"]),
                "dose_unit": str(condition["pert_dose_unit"]),
                "time": str(condition["pert_time"]),
                "time_unit": str(condition["pert_time_unit"]),
                "sig_id": str(condition["sig_id"]),
                "overlap_n": len(up_dn) + len(dn_up),
                "up_dn_genes": ";".join(up_dn),
                "dn_up_genes": ";".join(dn_up),
                "status": triage.get(key, {}).get("consensus_status", "not_triaged"),
                "role": triage.get(key, {}).get("experimental_role", "not_triaged"),
            })

    counts = {}
    for record in records:
        counts[record["compound"]] = counts.get(record["compound"], 0) + 1
    status_order = {
        "retain_tier1": 0,
        "retain_tier2": 1,
        "hold_identity_resolution": 2,
        "exclude_first_round": 3,
    }
    compounds = sorted(
        counts,
        key=lambda compound: (
            -counts[compound],
            status_order.get(triage.get(normalize_name(compound), {}).get("consensus_status"), 9),
            compound.lower(),
        ),
    )
    return records, compounds


def dsigdb_audit():
    with DSIGDB.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for rank, row in enumerate(rows, start=1):
        if "cgp 60474" in row["Term"].lower():
            return {
                "rank": rank,
                "terms": len(rows),
                "term": row["Term"],
                "genes": row["Genes"],
                "p_value": row["P-value"],
                "adjusted_p_value": row["Adjusted P-value"],
                "odds_ratio": row["Odds Ratio"],
                "combined_score": row["Combined Score"],
            }
    raise RuntimeError("CGP 60474 was not found in the DSigDB table")


def offsets(n):
    if n == 1:
        return [0.0]
    span = 0.34
    return [(-span / 2) + span * i / (n - 1) for i in range(n)]


def configure_style():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.linewidth": 0.8,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def draw(records, compounds, audit):
    configure_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.6, 5.9), facecolor="white")
    fig.subplots_adjust(left=0.25, right=0.91, top=0.80, bottom=0.22)

    y_lookup = {compound: len(compounds) - 1 - index for index, compound in enumerate(compounds)}
    cgp_y = y_lookup["CGP-60474"]
    ax.axhspan(cgp_y - 0.43, cgp_y + 0.43, color=HIGHLIGHT, zorder=0)

    for compound in compounds:
        center = y_lookup[compound]
        for query, marker, color in (("Strict", "o", PURPLE), ("Wide", "s", TEAL)):
            subset = [row for row in records if row["compound"] == compound and row["query"] == query]
            subset.sort(key=lambda row: (row["score"], row["cell_line"], float(row["dose"])))
            for row, delta in zip(subset, offsets(len(subset))):
                size = 30 + row["overlap_n"] * 12
                ax.scatter(
                    row["score"], center + delta, s=size, marker=marker,
                    facecolor=color, edgecolor="white", linewidth=0.8,
                    alpha=1.0 if compound == "CGP-60474" else 0.68, zorder=3,
                )

        total = sum(row["compound"] == compound for row in records)
        ax.text(0.1245, center, f"n = {total}", ha="left", va="center",
                fontsize=9.2, color=PURPLE if compound == "CGP-60474" else SLATE,
                weight="bold" if compound == "CGP-60474" else "normal")

    ax.set_yticks([y_lookup[compound] for compound in compounds], compounds)
    for label in ax.get_yticklabels():
        if label.get_text() == "CGP-60474":
            label.set_color(PURPLE)
            label.set_weight("bold")
        else:
            label.set_color(NAVY)
    ax.tick_params(axis="y", length=0, pad=8)

    ax.set_xlim(0.066, 0.132)
    ax.set_xticks([0.07, 0.08, 0.09, 0.10, 0.11, 0.12])
    ax.set_xlabel("L1000CDS2 reverse score (query-specific)", color=NAVY, labelpad=10)
    ax.tick_params(axis="x", colors=SLATE, length=3)
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#9AA8B6")

    fig.text(
        0.25, 0.965, "Condition-level L1000CDS2 reversal evidence",
        ha="left", va="top", fontsize=15.0, color=NAVY, weight="bold",
    )
    fig.text(
        0.25, 0.917,
        "Each point represents one returned LINCS perturbation condition; point area indicates reversed-gene overlap",
        ha="left", va="top", fontsize=9.4, color=SLATE,
    )

    query_legend = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=PURPLE,
               markeredgecolor="white", markersize=8, label="Strict query"),
        Line2D([0], [0], marker="s", linestyle="", markerfacecolor=TEAL,
               markeredgecolor="white", markersize=8, label="Wide query"),
    ]
    size_legend = [
        plt.scatter([], [], s=30 + n * 12, facecolor="#98A6B3", edgecolor="white", label=f"{n} genes")
        for n in (5, 10, 15)
    ]
    legend1 = ax.legend(
        handles=query_legend, loc="upper left", bbox_to_anchor=(0.0, 1.13),
        frameon=False, ncol=2, fontsize=8.8, handletextpad=0.4, columnspacing=1.2,
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=size_legend, title="Reverse overlap", loc="upper right",
        bbox_to_anchor=(1.0, 1.145), frameon=False, ncol=3, fontsize=8.0,
        title_fontsize=8.2, handletextpad=0.3, columnspacing=0.7,
    )

    fig.text(
        0.25, 0.045,
        "Strict and wide scores arise from different query sizes and should not be compared as a common effect scale.\n"
        "CGP-60474 was prioritized by cross-query recurrence (2 strict + 4 wide conditions), not by DSigDB significance.",
        ha="left", va="bottom", fontsize=8.0, color=SLATE, linespacing=1.35,
    )

    for suffix in ("png", "pdf", "svg", "tiff"):
        fig.savefig(
            OUTDIR / f"Fig_CGP60474_L1000_condition_level_reversal.{suffix}",
            dpi=600 if suffix in {"png", "tiff"} else None,
            bbox_inches="tight", facecolor="white",
        )
    plt.close(fig)

    fields = [
        "compound", "query", "score", "cell_line", "dose", "dose_unit", "time",
        "time_unit", "sig_id", "overlap_n", "up_dn_genes", "dn_up_genes", "status", "role",
    ]
    with (OUTDIR / "SourceData_CGP60474_L1000_condition_level.tsv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

    audit_fields = list(audit)
    with (OUTDIR / "DSigDB_CGP60474_audit.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(audit)

    readme = f"""# Evidence decision

CGP-60474 was not selected by significant DSigDB enrichment. In the 44-gene DSigDB query, it was term rank {audit['rank']} of {audit['terms']}, overlapped only {audit['genes']}, and had raw P={float(audit['p_value']):.4f} and BH-adjusted P={float(audit['adjusted_p_value']):.4f}. A DSigDB Top-compound plot that presents CGP-60474 as significantly prioritized would therefore be unsupported.

The defensible source is L1000CDS2 reverse mode. CGP-60474 was returned in two strict-query conditions and four wide-query conditions. This figure displays every condition-level observation for the seven compounds shared by both queries. Point area represents the number of genes overlapping in the reverse direction. Strict and wide scores are query-specific because the input signatures differ in size.

## Figure caption

**Condition-level L1000CDS2 reversal evidence for consensus compounds.** Each point denotes one LINCS perturbation condition returned in reverse mode for a compound shared by the strict and wide PDE2A-loss transcriptional queries. Point area represents the number of query genes reversed in that condition. CGP-60474 recurred in six conditions spanning the strict (two conditions) and wide (four conditions) queries and was retained as a Tier 1 CDK research lead. The scores are query-specific and do not quantify drug potency, binding affinity, or direct PDE2A engagement.
"""
    (OUTDIR / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    records_, compounds_ = load_records()
    draw(records_, compounds_, dsigdb_audit())
