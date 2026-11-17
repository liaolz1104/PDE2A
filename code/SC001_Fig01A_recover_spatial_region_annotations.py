from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.spatial import cKDTree

BASE = Path('/data/liaolz/PDE2A')
OUT = BASE/'results'/'Spatial_source_annotation_recovery_20260730'
OUT.mkdir(parents=True, exist_ok=True)

plot_path = BASE/'Figure2'/'Fig2_3_spatial_edge_core_panels'/'Fig2_3_spatial_edge_core_plot_data.csv'
anno_path = BASE/'source_annotation_recovery'/'source_pathologist_annotations_by_spot.csv'
plot = pd.read_csv(plot_path)
anno = pd.read_csv(anno_path)
anno = anno.rename(columns={'pathologist_anno.x':'source_pathologist_annotation', 'cluster_annotations':'source_author_region'})
cols = ['sample','barcode','source_pathologist_annotation','source_author_region']
anno = anno[cols].drop_duplicates(['sample','barcode'])
plot['barcode_core'] = plot['barcode'].str.replace(r'-\d+$','',regex=True)
anno['barcode_core'] = anno['barcode'].str.replace(r'-\d+$','',regex=True)
merged = plot.merge(anno.drop(columns=['barcode']), on=['sample','barcode_core'], how='left', validate='m:1')
merged['source_pathology_group'] = merged['source_pathologist_annotation'].map({
    'SCC':'SCC', 'Lymphocyte Positive Stroma':'stroma_lymphocyte_positive',
    'Lymphocyte Negative Stroma':'stroma_lymphocyte_negative', 'Glandular Stroma':'stroma_glandular',
    'Muscle':'muscle', 'Non-cancerous Mucosa':'normal_mucosa', 'Keratin':'keratin',
    'Artery/Vein':'vessel', 'Artifact':'artifact', 'Cautery':'artifact', 'Fold':'artifact', 'Edge Effects':'artifact',
}).fillna('unannotated')
merged['source_tumor_scc'] = merged['source_pathology_group'].eq('SCC')
merged.to_csv(OUT/'spatial_source_annotation_merged.csv', index=False)

audit = merged.groupby('sample', sort=False).agg(
    n_spots=('barcode','size'), n_source_annotation=('source_pathologist_annotation',lambda x:x.notna().sum()),
    n_scc=('source_tumor_scc','sum'), n_author_core=('source_author_region',lambda x:(x=='core').sum()),
    n_author_transitory=('source_author_region',lambda x:(x=='transitory').sum()),
    n_author_edge=('source_author_region',lambda x:(x=='edge').sum()),
    n_author_nc=('source_author_region',lambda x:(x=='nc').sum()),
).reset_index()
audit['annotation_rate'] = audit['n_source_annotation']/audit['n_spots']
audit.to_csv(OUT/'source_annotation_mapping_audit_by_sample.csv', index=False)
pd.crosstab(merged['source_pathology_group'], merged['source_author_region'], dropna=False).to_csv(OUT/'source_pathology_vs_author_region_crosstab.csv')

mapping = pd.DataFrame({'field':['source_series','source_article','doi','pmid','organism','cohort','technology','n_sections','n_unique_patients','pathologist','pathology_labels','author_region_labels','source_annotation_status'], 'value':
['GSE208253','Spatial transcriptomics reveals distinct and conserved tumor core and edge architectures that predict survival and targeted therapy response','10.1038/s41467-023-40271-4','37596273','Homo sapiens','HPV-negative oral squamous cell carcinoma','10x Genomics Visium, fresh-frozen 10 um H&E sections',12,10,'M.H.','SCC; Lymphocyte Positive Stroma; Lymphocyte Negative Stroma; Non-cancerous Mucosa; Glandular Stroma; Muscle; Keratin; Artery/Vein; Artifact; plus Cautery/Fold/Edge Effects in released objects','core; transitory; edge; nc','spot-level annotations recovered from author Seurat Robj files; TC/LE are author molecular cluster labels, not pathologist-drawn TC/LE ROIs']})
mapping.to_csv(OUT/'source_article_and_annotation_metadata.csv', index=False)

colors = {'SCC':'#D73027','stroma_lymphocyte_positive':'#4575B4','stroma_lymphocyte_negative':'#74ADD1','stroma_glandular':'#FDAE61','muscle':'#984EA3','normal_mucosa':'#66A61E','keratin':'#E6AB02','vessel':'#1B9E77','artifact':'#7F7F7F','unannotated':'#D9D9D9'}
labels = {'SCC':'SCC','stroma_lymphocyte_positive':'Lymph+ stroma','stroma_lymphocyte_negative':'Lymph− stroma','stroma_glandular':'Glandular stroma','muscle':'Muscle','normal_mucosa':'Normal mucosa','keratin':'Keratin','vessel':'Artery/vein','artifact':'Artifact / technical','unannotated':'Unannotated'}

def pitch_size(x, y):
    xy = np.c_[x, y]
    if len(xy)<3: return 25
    d = cKDTree(xy).query(xy, k=2)[0][:,1]
    return float(np.nanmedian(d))*0.82

samples = list(merged['sample'].drop_duplicates())
fig, axes = plt.subplots(2, 6, figsize=(24, 8), dpi=220)
for ax, s in zip(axes.flat, samples):
    d = merged[merged['sample']==s]
    size = pitch_size(d['x_hires'], d['y_hires'])
    for g, col in colors.items():
        q = d[d['source_pathology_group']==g]
        if len(q): ax.scatter(q['x_hires'], q['y_hires'], s=size, c=col, marker='o', linewidths=0, rasterized=True)
    ax.set_title(s, fontsize=10, weight='bold'); ax.set_aspect('equal'); ax.invert_yaxis(); ax.axis('off')
fig.suptitle('GSE208253: original pathologist spot annotations recovered from author Seurat objects', fontsize=16, weight='bold', y=0.99)
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=7, label=labels[g]) for g,c in colors.items() if g!='unannotated' or merged['source_pathology_group'].eq(g).any()]
fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=9)
fig.tight_layout(rect=[0,0.08,1,0.96]); fig.savefig(OUT/'FigS_spatial_source_pathologist_annotations.png', bbox_inches='tight'); fig.savefig(OUT/'FigS_spatial_source_pathologist_annotations.pdf', bbox_inches='tight'); plt.close(fig)

areg = {'core':'#4C78A8','transitory':'#F2CF5B','edge':'#E45756','nc':'#BDBDBD'}
fig, axes = plt.subplots(2, 6, figsize=(24, 8), dpi=220)
for ax, s in zip(axes.flat, samples):
    d = merged[merged['sample']==s]
    size = pitch_size(d['x_hires'], d['y_hires'])
    for g, col in areg.items():
        q=d[d['source_author_region']==g]
        if len(q): ax.scatter(q['x_hires'],q['y_hires'],s=size,c=col,marker='o',linewidths=0,rasterized=True)
    ax.set_title(s, fontsize=10, weight='bold'); ax.set_aspect('equal'); ax.invert_yaxis(); ax.axis('off')
fig.suptitle('GSE208253: author molecular region labels (Seurat metadata)', fontsize=16, weight='bold', y=0.99)
handles=[plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markersize=7,label=g.title()) for g,c in areg.items()]
fig.legend(handles=handles,loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(0.5,-0.01),fontsize=9)
fig.tight_layout(rect=[0,0.08,1,0.96]); fig.savefig(OUT/'FigS_spatial_source_author_regions.png', bbox_inches='tight'); fig.savefig(OUT/'FigS_spatial_source_author_regions.pdf', bbox_inches='tight'); plt.close(fig)

readme = f'''# Source annotation recovery — GSE208253 (2026-07-30)

## Source
Arora et al., *Nature Communications* (2023), “Spatial transcriptomics reveals distinct and conserved tumor core and edge architectures that predict survival and targeted therapy response”. DOI: 10.1038/s41467-023-40271-4; PMID: 37596273; GEO: GSE208253.

The source paper states that H&E sections were annotated by pathologist M.H. and exported from Loupe Browser v6. The released Seurat objects on Figshare (DOI 10.6084/m9.figshare.20304456.v1) contain the spot-level field `pathologist_anno.x` and the author molecular region field `cluster_annotations`.

## What was recovered
- 12 sections / 24,876 author Seurat spots, mapped by `sample` + `barcode` to the existing spatial plot table.
- Pathology labels: SCC, lymphocyte-positive/negative stroma, non-cancerous mucosa, glandular stroma, muscle, keratin, artery/vein, artifact; released objects additionally contain cautery/fold/edge-effects labels.
- Molecular region labels: core, transitory, edge, nc.

## Important boundary
The author paper does **not** provide pathologist-drawn TC/transition/LE ROIs. TC/transition/LE are molecular cluster annotations after Louvain clustering (CLDN4/SPR1B for TC; LAMC2/ITGA5 for LE). Therefore the old geometry/quartile labels remain derived proxies and are not overwritten. Use `source_pathologist_annotation` for pathology and `source_author_region` for the original molecular regions.
'''
(OUT/'README.md').write_text(readme, encoding='utf-8')
print('wrote', OUT)
print(audit.to_string(index=False))
