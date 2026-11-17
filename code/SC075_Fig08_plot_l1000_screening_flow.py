"""Render the evidence-traceable screening path that prioritized CGP-60474."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


OUTDIR = Path(__file__).resolve().parents[1] / "results" / "CGP60474_screening_flow_20260926"

PURPLE = "#7C3AED"
PURPLE_LIGHT = "#F3E8FF"
NAVY = "#17324D"
BLUE = "#2563EB"
BLUE_LIGHT = "#EAF2FF"
TEAL = "#0F766E"
TEAL_LIGHT = "#E7F6F3"
GREEN = "#15803D"
GREEN_LIGHT = "#ECF8EF"
SLATE = "#526476"
LINE = "#CBD5E1"
PALE = "#F8FAFC"


def configure_fonts():
    candidates = ["Microsoft YaHei", "Noto Sans CJK SC", "SimHei", "Arial Unicode MS"]
    available = {font.name for font in mpl.font_manager.fontManager.ttflist}
    family = next((font for font in candidates if font in available), "DejaVu Sans")
    mpl.rcParams.update({
        "font.family": family,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def rounded_box(ax, x, y, w, h, face, edge, linewidth=1.5, radius=0.018):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=linewidth,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, x1, x2, y=0.53):
    ax.add_patch(FancyArrowPatch(
        (x1, y), (x2, y), arrowstyle="-|>", mutation_scale=15,
        linewidth=1.5, color="#94A3B8", transform=ax.transAxes,
    ))


def stage(ax, x, w, number, title, lines, face, edge):
    rounded_box(ax, x, 0.22, w, 0.57, face, edge)
    ax.text(
        x + 0.022, 0.744, number, ha="left", va="center", fontsize=9.5,
        color="white", weight="bold", transform=ax.transAxes,
        bbox=dict(boxstyle="circle,pad=0.34", facecolor=edge, edgecolor="none"),
    )
    ax.text(
        x + 0.025, 0.675, title, ha="left", va="top", fontsize=12.5,
        color=NAVY, weight="bold", transform=ax.transAxes,
    )
    y = 0.600
    for text, color, weight, size in lines:
        ax.text(
            x + 0.025, y, text, ha="left", va="top", fontsize=size,
            color=color, weight=weight, linespacing=1.38, transform=ax.transAxes,
        )
        y -= 0.075 if "\n" not in text else 0.115


def build_figure():
    configure_fonts()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16.0, 7.2), facecolor="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.035, 0.942, "CGP-60474 的转录逆转筛选路径",
        ha="left", va="top", fontsize=22, color=NAVY, weight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.035, 0.885,
        "HN6 PDE2A 双向遗传扰动  →  PDE2A-loss 签名  →  L1000CDS2 反向匹配  →  跨签名共识优先级",
        ha="left", va="top", fontsize=11.5, color=SLATE, transform=ax.transAxes,
    )

    xs = [0.035, 0.225, 0.415, 0.605, 0.795]
    w = 0.155

    stage(ax, xs[0], w, "1", "HN6双向扰动", [
        ("PDE2A KD vs NC", NAVY, "bold", 10.5),
        ("PDE2A OE vs Vector", NAVY, "bold", 10.5),
        ("每组 n = 3 生物学重复", SLATE, "normal", 9.2),
        ("检验同一基因在KD与OE中\n是否呈相反变化方向", SLATE, "normal", 9.0),
    ], PALE, NAVY)

    stage(ax, xs[1], w, "2", "方向相反取交集", [
        ("loss-up", PURPLE, "bold", 10.5),
        ("KD ↑  ∩  OE ↓", NAVY, "bold", 11.0),
        ("loss-down", TEAL, "bold", 10.5),
        ("KD ↓  ∩  OE ↑", NAVY, "bold", 11.0),
    ], PURPLE_LIGHT, PURPLE)

    stage(ax, xs[2], w, "3", "构建两套签名", [
        ("STRICT", BLUE, "bold", 10.5),
        ("13 loss-up + 31 loss-down\n共 44 个基因", NAVY, "normal", 9.5),
        ("WIDE", GREEN, "bold", 10.5),
        ("61 loss-up + 205 loss-down\n共 266 个基因", NAVY, "normal", 9.5),
    ], BLUE_LIGHT, BLUE)

    stage(ax, xs[3], w, "4", "L1000CDS2反向检索", [
        ("分别提交 strict / wide", NAVY, "bold", 10.0),
        ("检索与PDE2A-loss表达方向\n相反的药物扰动签名", SLATE, "normal", 9.1),
        ("两套查询共有 7 个化合物", TEAL, "bold", 10.0),
        ("再按重复条件、身份可解释性\n和机制用途进行项目分层", SLATE, "normal", 8.7),
    ], TEAL_LIGHT, TEAL)

    stage(ax, xs[4], w, "5", "CGP-60474", [
        ("STRICT：2个条件", PURPLE, "bold", 10.4),
        ("A549、HCC515", SLATE, "normal", 9.0),
        ("WIDE：4个条件", PURPLE, "bold", 10.4),
        ("HA1E、HT29、MCF7", SLATE, "normal", 9.0),
        ("共6个LINCS条件；Tier 1\nCDK机制研究候选", NAVY, "bold", 9.3),
    ], PURPLE_LIGHT, PURPLE)

    for left, right in zip(xs[:-1], xs[1:]):
        arrow(ax, left + w + 0.012, right - 0.012)

    # A single evidence strip keeps the output as one standalone figure rather than a panel composite.
    rounded_box(ax, 0.035, 0.075, 0.915, 0.105, "white", LINE, linewidth=1.0, radius=0.012)
    ax.add_patch(Rectangle((0.035, 0.075), 0.008, 0.105, facecolor=PURPLE,
                           edgecolor="none", transform=ax.transAxes))
    ax.text(
        0.055, 0.145, "结论边界", ha="left", va="center", fontsize=10.0,
        color=PURPLE, weight="bold", transform=ax.transAxes,
    )
    ax.text(
        0.055, 0.103,
        "CGP-60474因跨strict/wide签名重复命中而被优先保留；这不是药效、结合亲和力或PDE2A直接靶向排名。",
        ha="left", va="center", fontsize=10.0, color=NAVY, transform=ax.transAxes,
    )
    ax.text(
        0.035, 0.030,
        "注：当前签名来自vendor/rounded-count DESeq2敏感性分析；项目中无rescue组，因此表述为“双向支持的PDE2A-loss签名”。",
        ha="left", va="bottom", fontsize=8.4, color=SLATE, transform=ax.transAxes,
    )

    for suffix in ("png", "pdf", "svg"):
        fig.savefig(
            OUTDIR / f"Fig_CGP60474_screening_flow.{suffix}",
            dpi=600 if suffix == "png" else None,
            bbox_inches="tight", facecolor="white",
        )
    plt.close(fig)

    source = """item\tvalue\tsource\nKD_vs_NC_replicates\t3\tPDE2A_KD_OE_input_audit.md\nOE_vs_Vector_replicates\t3\tPDE2A_KD_OE_input_audit.md\nstrict_loss_up\t13\tpde2a_loss_up.grp\nstrict_loss_down\t31\tpde2a_loss_down.grp\nwide_loss_up\t61\tPDE2A_KD_OE_input_audit.md\nwide_loss_down\t205\tPDE2A_KD_OE_input_audit.md\nstrict_wide_shared_compounds\t7\tstrict_wide_consensus_compounds.csv\nCGP-60474_strict_conditions\t2\tstrict_wide_consensus_compounds.csv\nCGP-60474_strict_cell_lines\tA549;HCC515\tstrict_wide_consensus_compounds.csv\nCGP-60474_wide_conditions\t4\tstrict_wide_consensus_compounds.csv\nCGP-60474_wide_cell_lines\tHA1E;HT29;MCF7\tstrict_wide_consensus_compounds.csv\nCGP-60474_triage\tretain_tier1; CDK research lead; docking_eligible=No\tconsensus_candidate_triage.csv\n"""
    (OUTDIR / "SourceData_CGP60474_screening_flow.tsv").write_text(source, encoding="utf-8-sig")

    caption = """# Figure caption

**CGP-60474的转录逆转筛选路径。** 在HN6细胞中，分别比较PDE2A敲低与阴性对照、PDE2A过表达与空载体对照，并保留两种扰动中变化方向相反的基因，以构建双向支持的PDE2A-loss转录签名。strict签名包含13个loss-up和31个loss-down基因，wide敏感性签名包含61个loss-up和205个loss-down基因。两套签名分别提交L1000CDS2进行反向扰动检索，共得到7个跨查询共有化合物。CGP-60474在strict和wide查询中分别出现于2和4个LINCS条件，因此被列为Tier 1 CDK机制研究候选。该优先级表示转录逆转证据的重复性，不代表药效、PDE2A结合亲和力或直接靶向关系。

## Evidence note

当前签名基于vendor/rounded-count DESeq2敏感性分析，项目中没有rescue组。正式论文中应使用“bidirectionally supported PDE2A-loss transcriptional signature”，不应写成“rescue signature”。
"""
    (OUTDIR / "README.md").write_text(caption, encoding="utf-8")


if __name__ == "__main__":
    build_figure()
