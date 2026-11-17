"""Redraw frozen pharmacogenomic predictions; run only on the analysis server."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    read = lambda name: pd.read_csv(args.input / name, sep='\t')
    master = read('candidate_master_table.tsv')
    models = read('all_model_results.tsv')
    pred = read('all_sample_predictions.tsv')
    sample = read('hn6_sample_audit.tsv')
    selected = master[master.final_tier.isin(['Tier 1', 'Tier 2'])].sort_values(['final_tier', 'drug'])
    order = selected.drug.str.lower().tolist()
    for frame in (models, pred):
        frame['drug'] = frame.drug.str.lower()
        frame['drug_id'] = frame.drug_id.astype(str)
    keys = ['database', 'drug_id', 'drug']
    passed = models[models.model_pass.eq(True)].copy()
    assert not passed.duplicated(keys).any()
    data = pred.merge(passed[keys + ['response_sd']], on=keys, validate='many_to_one')
    assert not data.duplicated(keys + ['sample']).any()
    assert data.lower_is_more_sensitive.eq(True).all()
    assert data.groupby(keys).size().eq(12).all()
    assert (data.response_sd > 0).all()
    chosen = data[data.drug.isin(order)].copy()
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.titlesize': 10, 'axes.labelsize': 8, 'axes.linewidth': .65,
        'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    colors = {'GDSC1': '#2878A5', 'GDSC2': '#B36949', 'CTRP2': '#527D61'}
    cmap = LinearSegmentedColormap.from_list('response', ['#316A91', '#FAFAFA', '#B45F48'])
    manifest = []

    def save(fig, name, source, caption):
        assert len([a for a in fig.axes if a.get_label() != '<colorbar>']) == 1
        fig.tight_layout(pad=1.1)
        for ext in ('png', 'pdf', 'svg'):
            fig.savefig(out / f'{name}.{ext}', dpi=400, facecolor='white')
        plt.close(fig)
        source.to_csv(out / f'{name}_source.tsv', sep='\t', index=False)
        manifest.append({'figure': name, 'caption': caption})

    audit = chosen.groupby(keys + ['group']).predicted_response.mean().unstack('group')
    audit['NC_minus_Vector'] = audit.NC - audit.Vector
    audit.reset_index().to_csv(out / 'control_baseline_audit.tsv', sep='\t', index=False)
    sample.to_csv(out / 'sample_audit.tsv', sep='\t', index=False)
    matrices = {}
    centered_sources = {}
    for control, treatment in [('NC', 'KD'), ('Vector', 'OE')]:
        part = chosen[chosen.group.isin([control, treatment])].copy()
        baseline = part[part.group.eq(control)].groupby(keys).predicted_response.mean().rename('control_mean')
        part = part.merge(baseline, on=keys, validate='many_to_one')
        part['control_centered_response'] = (part.predicted_response - part.control_mean) / part.response_sd
        means = part[part.group.eq(control)].groupby(keys).control_centered_response.mean()
        assert np.allclose(means, 0, atol=1e-10)
        samples = [f'HN6_{g}{suffix}' for g in (control, treatment) for suffix in ('', '_r1', '_r2')]
        # Median is used only to display dimensionless model projections, not to infer clinical response.
        matrix = part.groupby(['drug', 'sample']).control_centered_response.median().unstack('sample').reindex(index=order, columns=samples)
        assert matrix.notna().all().all()
        matrices[treatment] = matrix
        centered_sources[treatment] = part
    limit = max(float(np.abs(m.values).max()) for m in matrices.values())
    for control, treatment in [('NC', 'KD'), ('Vector', 'OE')]:
        matrix = matrices[treatment]
        fig, ax = plt.subplots(figsize=(4.8, max(3, len(order)*.34)))
        im = ax.imshow(matrix, cmap=cmap, vmin=-limit, vmax=limit, aspect='auto')
        ax.set_yticks(range(len(order)), [d.capitalize() for d in order])
        ax.set_xticks(range(6), [f'{g}-{i}' for g in (control, treatment) for i in (1, 2, 3)], rotation=45, ha='right')
        ax.axvline(2.5, color='white', lw=2)
        ax.set_title(f'HN6: {treatment} versus {control}', pad=10)
        ax.tick_params(length=0)
        for spine in ax.spines.values(): spine.set_visible(False)
        fig.colorbar(im, ax=ax, shrink=.8, label='Control-centered predicted response\n(training-response SD units)')
        save(fig, f'HN6_{treatment}_heatmap', centered_sources[treatment],
             f'Three RNA-seq samples per group; replicate independence not re-established. Each model is centered on its own {control} mean and divided by training-response SD, then median-aggregated per drug/sample. Both heatmaps share a symmetric unclipped scale. Negative indicates lower predicted response relative to the corresponding control; no NC-versus-Vector sensitivity inference. Not measured drug response.')

        records = passed[passed.drug.isin(order)].copy()
        records['effect'] = records[f'std_delta_{treatment}']
        records['low'] = records[f'bootstrap_{treatment}_low'] / records.response_sd
        records['high'] = records[f'bootstrap_{treatment}_high'] / records.response_sd
        fig, ax = plt.subplots(figsize=(5.6, max(3, len(order)*.48)))
        for j, (db, color) in enumerate(colors.items()):
            subset = records[records.database.eq(db)]
            if subset.empty:
                continue
            yy = subset.drug.map({d:i for i,d in enumerate(order)}).to_numpy() + (j-1)*.2
            ax.hlines(yy, subset.low, subset.high, color=color, lw=1)
            ax.scatter(subset.effect, yy, color=color, s=19, label=db, zorder=3)
        ax.axvline(0, color='.6', lw=.7, linestyle='--')
        ax.set_yticks(range(len(order)), [d.capitalize() for d in order])
        ax.invert_yaxis()
        ax.set_xlabel(f'Predicted response change ({treatment} - {control}) / training SD')
        ax.set_title(f'HN6 {treatment}: model-specific effects')
        ax.legend(loc='best', fontsize=7)
        save(fig, f'HN6_{treatment}_model_effects', records,
             'Dots are frozen full-model effects. Lines are 2.5th-97.5th percentiles across 100 training-cell-line bootstrap refits, not patient confidence intervals. All passing model records for frozen Tier 1/2 candidates are retained; missing databases have no mark. No new p values or rankings are inferred.')

        landscape = master[master.standardized_delta_KD_GDSC.notna()].copy()
        effect_col = f'standardized_delta_{treatment}_GDSC'
        landscape = landscape.sort_values(effect_col)
        fig, ax = plt.subplots(figsize=(5.5, max(4, len(landscape)*.16)))
        ax.scatter(landscape[effect_col], range(len(landscape)), s=12, color='#2878A5')
        ax.axvline(0, color='.65', lw=.6)
        ax.set_yticks(range(len(landscape)), landscape.drug.str.capitalize(), fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel(f'Median GDSC standardized change ({treatment} - {control})')
        ax.set_title('Quality-gated FDA-drug screening landscape')
        save(fig, f'FDA_drug_landscape_{treatment}', landscape,
             'All frozen master-table compounds with passing GDSC estimates, ordered by displayed effect. Other screened compounds are retained in the source master table, not assigned zero effects. This is an exploratory prediction landscape, not an efficacy ranking.')

        for (db, drug_id, drug), part in chosen.groupby(keys):
            part = part[part.group.isin([control, treatment])]
            fig, ax = plt.subplots(figsize=(2.9, 3.1))
            for x, group in enumerate((control, treatment)):
                values = part.loc[part.group.eq(group)].sort_values('sample').predicted_response.to_numpy()
                assert len(values) == 3
                ax.scatter(x + np.array([-.065, 0, .065]), values, s=26,
                           color='.5' if x == 0 else '#2878A5', zorder=3)
                ax.plot([x-.17, x+.17], [np.mean(values)]*2, color='black', lw=1)
            ax.set_xlim(-.5, 1.5)
            ax.set_xticks([0, 1], [f'{control}\nn = 3', f'{treatment}\nn = 3'])
            ax.set_ylabel(f'Predicted response ({db} model scale)')
            ax.set_title(f'{drug.capitalize()} | {db}', pad=9)
            name = f'HN6_{drug.replace(" ", "_")}_{db}_{drug_id}_{treatment}'
            save(fig, name, part,
                 'Each point is one RNA-seq sample prediction; line is the arithmetic mean. No paired sample identity assumed. Original model response scale retained without relabeling as measured IC50/AUC. Lower is greater predicted sensitivity; no inferential p value is assigned to these three predictions.')

    pd.DataFrame(manifest).to_csv(out / 'figure_manifest.tsv', sep='\t', index=False)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.input.glob('*.tsv')}
    summary = {'n_figures': len(manifest), 'n_candidates': len(order), 'candidates': order,
               'n_screened_master': len(master), 'n_GDSC_landscape': len(landscape),
               'n_control_model_records': len(audit), 'n_NC_lower_than_Vector': int((audit.NC_minus_Vector < 0).sum()),
               'heatmap_symmetric_limit': limit, 'input_sha256': hashes,
               'new_model_training': False, 'render_host': __import__('socket').gethostname()}
    (out / 'run_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    (out / 'README.md').write_text(
        '# Standalone drug prediction figures\n\n'
        'Frozen source: FDA screen 2026-08-01. No retraining or candidate reclassification. Later validation outputs are not merged into this frozen screen.\n\n'
        'Original HN6 heatmap used within-model Z scores across all 12 samples before median aggregation. New heatmaps use separate matched-control centering and training SD scaling, with the same symmetric color limits. This changes the displayed estimand, not the raw predictions. NC and Vector are distinct controls; the cause of baseline differences remains unresolved.\n\n'
        'Original source scripts: integrate_and_plot.py (QC/landscape/evidence), plot_manuscript_conclusion.py (prioritization), plot_tcga_hn6_validation_figures.py (HN6 heatmap and TCGA projections). Original files are not overwritten. Screenshot order cannot be confirmed from a text-only attachment.\n\n'
        'Figure-level source data and captions are in figure_manifest.tsv. PNG 400 dpi, vector PDF and editable-text SVG are provided separately. No composite figures. Per-group sample count is three; independent biological replication must be confirmed before inferential use. Training bootstrap intervals quantify model-training uncertainty only.\n\n'
        'The old binary L1000FWD_support field is not used as negative evidence: false can mean no mapped hit rather than a tested negative. No mixed-scale evidence heatmap is recreated.\n\n'
        'Presentation reference: https://www.thno.org/v12p7267.htm#headingA5 (Fig. 6; ranking, response distributions, evidence separation). The published biological conclusions and numerical thresholds are not transferred to PDE2A.\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
