#!/usr/bin/env python3
"""Standardize existing KD/OE DESeq2 results and build reversible signatures."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
IN_DIR = ROOT / 'results' / 'rounded_count_deseq2'
OUT_DIR = ROOT / 'results' / 'reversible_signature'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def standardize(path: Path, contrast: str) -> pd.DataFrame:
    x = pd.read_csv(path, sep=None, engine='python')
    x = x.rename(columns={'gene_symbol': 'gene', 'log2FoldChange': 'logFC', 'pvalue': 'P.Value', 'padj': 'adj.P.Val', 'stat': 'statistic'})
    required = {'gene', 'logFC', 'P.Value', 'adj.P.Val', 'statistic'}
    missing = required - set(x.columns)
    if missing:
        raise ValueError(f'{path.name} missing {sorted(missing)}')
    x = x[['gene', 'logFC', 'P.Value', 'adj.P.Val', 'statistic']].copy()
    x['gene'] = x['gene'].astype(str).str.strip()
    for col in ['logFC', 'P.Value', 'adj.P.Val', 'statistic']:
        x[col] = pd.to_numeric(x[col], errors='coerce')
    x = x.dropna(subset=['gene', 'logFC', 'P.Value', 'adj.P.Val', 'statistic'])
    x['direction'] = x['logFC'].map(lambda v: 'up' if v > 0 else 'down' if v < 0 else 'neutral')
    x['contrast'] = contrast
    return x


kd = standardize(IN_DIR / 'KD_vs_NC.csv', 'PDE2A-KD vs NC')
oe = standardize(IN_DIR / 'OE_vs_Vector.csv', 'PDE2A-OE vs Vector')
kd.to_csv(OUT_DIR / 'PDE2A_KD_vs_Control_DEG.tsv', sep='\t', index=False)
oe.to_csv(OUT_DIR / 'PDE2A_OE_vs_Control_DEG.tsv', sep='\t', index=False)

merged = kd.merge(oe, on='gene', suffixes=('_KD', '_OE'), how='inner')
merged['loss_class'] = 'not_reversible'
wide = (merged['adj.P.Val_KD'] < 0.1) & (merged['adj.P.Val_OE'] < 0.1)
large = (merged['logFC_KD'].abs() > 0.25) & (merged['logFC_OE'].abs() > 0.25)
merged.loc[wide & large & (merged['logFC_KD'] > 0) & (merged['logFC_OE'] < 0), 'loss_class'] = 'loss_up'
merged.loc[wide & large & (merged['logFC_KD'] < 0) & (merged['logFC_OE'] > 0), 'loss_class'] = 'loss_down'
merged.to_csv(OUT_DIR / 'PDE2A_KD_OE_reversible_all_genes.tsv', sep='\t', index=False)

core = merged[merged['loss_class'].isin(['loss_up', 'loss_down'])].copy()
core.to_csv(OUT_DIR / 'PDE2A_KD_OE_reversible_signature.tsv', sep='\t', index=False)
for cls in ['loss_up', 'loss_down']:
    genes = core.loc[core['loss_class'] == cls, 'gene'].drop_duplicates()
    (OUT_DIR / f'PDE2A_{cls}.grp').write_text('\n'.join(genes) + ('\n' if len(genes) else ''), encoding='utf-8')

summary = pd.DataFrame({
    'contrast': ['KD vs NC', 'OE vs Vector'],
    'n_genes': [len(kd), len(oe)],
    'n_fdr_lt_0.1_abs_logFC_gt_0.25': [int(((kd['adj.P.Val'] < 0.1) & (kd['logFC'].abs() > 0.25)).sum()), int(((oe['adj.P.Val'] < 0.1) & (oe['logFC'].abs() > 0.25)).sum())],
    'source_type': ['DESeq2 result table; raw matrix not present in workspace', 'DESeq2 result table; raw matrix not present in workspace'],
})
summary.to_csv(OUT_DIR / 'PDE2A_KD_OE_input_summary.tsv', sep='\t', index=False)
print(summary.to_string(index=False))
print('reversible_loss_up', int((core.loss_class == 'loss_up').sum()))
print('reversible_loss_down', int((core.loss_class == 'loss_down').sum()))
