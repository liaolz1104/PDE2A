#!/usr/bin/env python3
"""Compact bracket-style redraw of the frozen patient-balanced zone summary."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont, fontManager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


parser = argparse.ArgumentParser()
parser.add_argument("--summary", type=Path, required=True)
parser.add_argument("--stats", type=Path, required=True)
parser.add_argument("--out-dir", type=Path, required=True)
args = parser.parse_args()
args.out_dir.mkdir(parents=True, exist_ok=True)

if Path("C:/Windows/Fonts/arial.ttf").is_file():
    fonts = [Path("C:/Windows/Fonts/arial.ttf"),
             Path("C:/Windows/Fonts/arialbd.ttf")]
else:
    root = Path("/data/liaolz/PDE2A/GSE268014_scRNA_spatial_integration/")
    fonts = [root / "results/11_invasive_edge_ecosystem_main_figure/fonts/Arial.ttf",
             root / "results/11_invasive_edge_ecosystem_main_figure/fonts/Arial-Bold.ttf"]
for font in fonts:
    if not font.is_file():
        raise FileNotFoundError(font)
    fontManager.addfont(str(font))
if Path(findfont("Arial", fallback_to_default=False)).name.lower() != "arial.ttf":
    raise RuntimeError("Arial TrueType not resolved")
plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": .8,
    "savefig.dpi": 300, "svg.fonttype": "none",
})

summary = pd.read_csv(args.summary, sep="\t")
stats = pd.read_csv(args.stats, sep="\t").set_index("program")
zones = ["Core", "Intermediate", "Edge"]
programs = [
    ("MP1 invasion/pEMT", "#205A72"),
    ("Cancer_EMT", "#70445F"),
    ("Fibroblast abundance", "#826341"),
]
if (len(summary) != 9 or
        set(summary.program) != {p for p, _ in programs} or
        summary.groupby(["program", "zone"]).size().ne(1).any() or
        set(summary.zone) != set(zones) or
        summary.patient_n.ne(8).any()):
    raise AssertionError("Expected the frozen 8-patient, 3-program, 3-zone summary")

fig, ax = plt.subplots(figsize=(5.4, 3.3), facecolor="white")
fig.subplots_adjust(left=.17, right=.68, top=.93, bottom=.24)
ax.axhline(0, color="#777777", linewidth=.75, zorder=0)
rows = []
handles = []
for program, color in programs:
    frame = summary.loc[summary.program.eq(program)].set_index("zone").loc[zones]
    frozen = stats.loc[program]
    p = float(frozen.paired_edge_core_p_two_sided)
    direction = int(frozen.edge_greater_core_n)
    if not np.allclose(frame.edge_core_wilcoxon_p_two_sided, p, atol=1e-12):
        raise AssertionError(f"Paired P changed: {program}")
    if (not frame.edge_greater_core_patients.eq(direction).all() or
            not np.allclose(frame.loc["Core", ["mean_paired_change_pct",
                         "bootstrap_95_ci_low", "bootstrap_95_ci_high"]], 0)):
        raise AssertionError(f"Frozen summary changed: {program}")
    stars = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
    y = frame.mean_paired_change_pct.to_numpy(float)
    lo = frame.bootstrap_95_ci_low.to_numpy(float)
    hi = frame.bootstrap_95_ci_high.to_numpy(float)
    x = np.arange(3)
    ax.errorbar(x, y, yerr=np.vstack((y-lo, hi-y)), color=color,
                ecolor="#3D3D3D", linewidth=1, elinewidth=.8,
                capsize=2.2, capthick=.8, marker="o", markersize=5,
                markeredgecolor="#343434", markeredgewidth=.6,
                zorder=4)
    handles.append(Line2D([], [], color=color, marker="o", markersize=5,
                          markeredgecolor="#343434", markeredgewidth=.6,
                          linewidth=1))
    for zone in zones:
        record = frame.loc[zone].to_dict()
        rows.append({"program": program, "zone": zone, **record,
                     "edge_core_stars": stars,
                     "paired_test": "two-sided paired Wilcoxon; unadjusted"})

for height, (program, color) in zip((61.5, 67, 72.5), programs):
    p = float(stats.loc[program, "paired_edge_core_p_two_sided"])
    stars = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
    ax.plot([0, 0, 2, 2], [height-1.8, height, height, height-1.8],
            color=color, linewidth=.8, solid_capstyle="butt", clip_on=False)
    ax.annotate(stars, xy=(1, height), xytext=(0, 0),
                textcoords="offset points", color=color, fontsize=7,
                ha="center", va="bottom")

ax.set_xlim(-.13, 2.22)
ax.set_ylim(-10, 76)
ax.set_xticks(range(3), zones)
ax.set_yticks([0, 20, 40, 60])
ax.set_xlabel("Author-defined malignant zone", labelpad=6)
ax.set_ylabel("Paired change from Core\n(% of cohort Core mean)", labelpad=5)
ax.tick_params(direction="out", length=2.5, width=.8)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(False)

fig.legend(handles, [p for p, _ in programs], title="Program",
           loc="upper left", bbox_to_anchor=(.69, .78), frameon=False,
           handlelength=1.7, handletextpad=.65, labelspacing=1.2,
           fontsize=7, title_fontsize=7)
fig.text(.17, .075,
         "Mean ± bootstrap 95% CI (n=8 patients); * P<0.05, ** P<0.01.\n"
         "Core–Edge: two-sided paired Wilcoxon, unadjusted.",
         ha="left", va="center", fontsize=7, color="#444444")

stem = args.out_dir / "GSE268014_patient_balanced_zone_gradient_centered_stars"
fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
fig.savefig(stem.with_suffix(".svg"), facecolor="white")
plt.close(fig)
pd.DataFrame(rows).to_csv(args.out_dir / f"{stem.name}_source_data.tsv",
                          sep="\t", index=False)
with Image.open(stem.with_suffix(".png")) as image:
    if (image.size != (1620, 990) or
            any(abs(v - 300) > 1 for v in image.info.get("dpi", (0, 0)))):
        raise AssertionError("PNG dimensions or resolution mismatch")
print(stem.with_suffix(".png"))
