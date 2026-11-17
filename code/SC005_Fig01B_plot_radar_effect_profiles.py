#!/usr/bin/env python3
"""Replace hard-to-read filled radar charts with interpretable effect profiles.

The radial panel preserves the requested radar-like overview, while the paired
effect panel carries direction, magnitude and BH-adjusted significance.
"""
from pathlib import Path
import argparse
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

warnings.filterwarnings('ignore')
mpl.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.titlesize': 10,
    'axes.labelsize': 8.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
})
EDGE = '#D55E00'; CORE = '#0072B2'; DARK = '#333333'; LIGHT = '#BDBDBD'

def bh(p):
    p = np.asarray(p, float); order = np.argsort(p); q = np.empty_like(p)
    q[order] = np.minimum.accumulate((p[order] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def program_label(x):
    x = x.split('.', 1)[-1].replace('_score', '').replace('_', ' ')
    return x.replace('P53 signaling pathway', 'P53 signaling')

def paired_effect(df, id_col, group_col, features, edge_name, core_name):
    rows = []
    for feat in features:
        wide = df.pivot(index=id_col, columns=group_col, values=feat).dropna()
        edge = wide[edge_name]; core = wide[core_name]
        test = wilcoxon(edge, core, alternative='two-sided', method='auto')
        rows.append({'program': feat, 'n_pairs': len(wide), 'delta_edge_minus_core': float((edge-core).mean()),
                     'p': float(test.pvalue), 'edge_values': edge, 'core_values': core})
    out = pd.DataFrame([{k:v for k,v in r.items() if k not in {'edge_values','core_values'}} for r in rows])
    out['q'] = bh(out['p'].to_numpy())
    return out, rows

def save(fig, out, stem):
    out.mkdir(parents=True, exist_ok=True)
    for ext, kw in [('png', {'dpi': 600}), ('tiff', {'dpi': 600}), ('pdf', {}), ('svg', {})]:
        fig.savefig(out / f'{stem}.{ext}', bbox_inches='tight', facecolor='white', **kw)
    plt.close(fig)

def draw_profile(eff, rows, out, stem, title, layer_note):
    labels = eff['program'].tolist(); n = len(labels)
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    max_abs = max(float(np.max(np.abs(eff['delta_edge_minus_core']))), 0.1)
    radial = 0.5 + 0.38 * eff['delta_edge_minus_core'].to_numpy() / max_abs
    radial = np.r_[radial, radial[0]]; ang = np.r_[angles, angles[0]]
    fig = plt.figure(figsize=(11.2, 5.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45], wspace=0.36)
    ax = fig.add_subplot(gs[0, 0], polar=True)
    ax.set_theta_offset(np.pi/2); ax.set_theta_direction(-1)
    ax.plot(ang, radial, color=EDGE, lw=2.1, marker='o', ms=4, zorder=3)
    ax.plot(np.r_[angles, angles[0]], np.repeat(0.5, n+1), color=LIGHT, lw=1.0, ls='--', zorder=1)
    ax.set_ylim(0.05, 0.95); ax.set_yticks([0.12, 0.31, 0.5, 0.69, 0.88]); ax.set_yticklabels(['−', '−', '0', '+', '+'], color='#666666', fontsize=7)
    ax.set_xticks(angles); ax.set_xticklabels([program_label(x) for x in labels], fontsize=7)
    ax.grid(color='#D9D9D9', lw=0.7); ax.spines['polar'].set_color('#777777')
    ax.set_title('Radial effect profile', pad=18, fontweight='bold')
    ax.text(0.5, -0.13, 'radius is scaled; dashed circle = Edge−Core = 0', transform=ax.transAxes, ha='center', fontsize=7, color='#666666')

    ax2 = fig.add_subplot(gs[0, 1])
    order = np.argsort(eff['delta_edge_minus_core'].to_numpy())
    y = np.arange(n); vals = eff['delta_edge_minus_core'].to_numpy()[order]; q = eff['q'].to_numpy()[order]
    labs = [program_label(labels[i]) for i in order]
    for yi, row in zip(y, [rows[i] for i in order]):
        pair = row['edge_values'] - row['core_values']
        ax2.scatter(pair, np.full(len(pair), yi), color=LIGHT, s=14, alpha=0.8, zorder=2)
    colors = [EDGE if v >= 0 else CORE for v in vals]
    ax2.barh(y, vals, color=colors, alpha=0.92, height=0.58, zorder=1)
    ax2.axvline(0, color=DARK, lw=0.8)
    ax2.set_yticks(y); ax2.set_yticklabels(labs); ax2.set_xlabel('Mean paired Edge − Core score'); ax2.set_title(f'Paired effect sizes ({layer_note})', fontweight='bold')
    xpad = 0.04 * max(1.0, float(np.max(np.abs(vals))))
    for yi, vi, qi in zip(y, vals, q):
        ax2.text(vi + (xpad if vi >= 0 else -xpad), yi, f'q={qi:.3g}', va='center', ha='left' if vi >= 0 else 'right', fontsize=7)
    ax2.spines[['top','right']].set_visible(False); ax2.grid(axis='x', color='#EEEEEE', lw=0.6)
    fig.suptitle(title, fontsize=13, fontweight='bold', y=0.99)
    save(fig, out, stem)
    eff.drop(columns=[], inplace=False).to_csv(out / f'{stem}_statistics.csv', index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spatial_dir', default='/data/liaolz/PDE2A/Figure2/Fig2_3_spatial_edge_core_panels')
    ap.add_argument('--cell_dir', default='/data/liaolz/PDE2A/Figure2/Fig2_7_8_edge_core_biology_cellchat')
    ap.add_argument('--out', default='/data/liaolz/PDE2A/results/Spatial_radar_redraw_20260713')
    a = ap.parse_args(); out = Path(a.out)

    sdir = Path(a.spatial_dir)
    s = pd.read_csv(sdir / 'Fig2_3_spatial_behavior_scores_by_sample_region.csv')
    sf = [c for c in s.columns if c.endswith('_z') and c.split('_z')[0] in [x.replace('.', '.') for x in ['1.Hypoxia','2.Angiogenesis','3.Glycolysis','4.Cell_cycle','5.DNA_replication','6.P53_signaling_pathway','7.Pyrimidine_metabolism','8.Stem_score','9.EMT']]]
    # Keep the audited z-score columns in their original order.
    sf = [f'{x}_z' for x in ['1.Hypoxia','2.Angiogenesis','3.Glycolysis','4.Cell_cycle','5.DNA_replication','6.P53_signaling_pathway','7.Pyrimidine_metabolism','8.Stem_score','9.EMT']]
    sp = s[s.edge_core_region.isin(['Edge-like','Core-like'])].copy()
    sp = sp.rename(columns={'edge_core_region':'region'}); sp['region'] = sp['region'].map({'Edge-like':'Edge','Core-like':'Core'})
    seff, srows = paired_effect(sp, 'sample', 'region', sf, 'Edge', 'Core')
    seff['program'] = [x.replace('_z','') for x in seff['program']]
    draw_profile(seff, srows, out, 'Fig2_3I_spatial_edge_core_behavior_v2', 'Spatial Edge/Core functional effect profile', 'spatial spots; n=12 pairs')

    cdir = Path(a.cell_dir)
    c = pd.read_csv(cdir / 'Fig2_7_patient_region_scores.csv')
    cf = [f'{i}_{name}_score_z' for i, name in [(1,'Hypoxia'),(2,'Angiogenesis'),(3,'Glycolysis'),(4,'Cell_cycle'),(5,'DNA_replication'),(6,'P53_signaling_pathway'),(7,'Pyrimidine_metabolism'),(8,'Stem_score'),(9,'EMT')]]
    ceff, crows = paired_effect(c, 'patient', 'region', cf, 'Edge', 'Core')
    ceff['program'] = [x.replace('_score_z','') for x in ceff['program']]
    draw_profile(ceff, crows, out, 'Fig2_7_malignant_biology_effect_profile_v4', 'Malignant-cell Edge/Core functional effect profile', 'malignant cells; n=9 pairs')
    (out / 'README.md').write_text('# Radar replacements\n\nFilled polygons were replaced by a radial effect profile plus paired effect-size plot. The effect-size panel is the inferential panel; q values are BH-adjusted within each layer. Spatial and malignant-cell layers are kept separate.\n', encoding='utf-8')

if __name__ == '__main__':
    main()
