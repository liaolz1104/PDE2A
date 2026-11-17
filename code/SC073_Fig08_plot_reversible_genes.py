#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == 'mechanism_drug_work' else HERE
OUT = ROOT / 'results' / 'PDE2A_reversible_figures'
OUT.mkdir(parents=True, exist_ok=True)
SIG_CANDIDATES = [ROOT / 'mechanism_drug_work' / 'results' / 'reversible_signature_outputs',
                  ROOT / 'results' / 'reversible_signature']
SIG = next(p for p in SIG_CANDIDATES if p.exists())
PROJ_CANDIDATES = [ROOT / 'results' / 'PDE2A_reversible_projection', ROOT / 'PDE2A_reversible_projection']
PROJ = next(p for p in PROJ_CANDIDATES if p.exists())

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
                     'axes.labelsize':9,'legend.fontsize':8,'figure.dpi':150,
                     'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False})
BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#009E73'; GREY='#8C8C8C'

kd = pd.read_csv(SIG/'PDE2A_KD_vs_Control_DEG.tsv', sep='\t')
oe = pd.read_csv(SIG/'PDE2A_OE_vs_Control_DEG.tsv', sep='\t')
kd['gene']=kd['gene'].astype(str); oe['gene']=oe['gene'].astype(str)
merged=kd[['gene','logFC','adj.P.Val']].rename(columns={'logFC':'kd_logFC','adj.P.Val':'kd_fdr'}).merge(
    oe[['gene','logFC','adj.P.Val']].rename(columns={'logFC':'oe_logFC','adj.P.Val':'oe_fdr'}),on='gene')
rows=[]
for fdr,fc in [(0.05,0.5),(0.1,0.5),(0.1,0.25),(0.2,0.25)]:
    up=((merged.kd_fdr<fdr)&(merged.oe_fdr<fdr)&(merged.kd_logFC>fc)&(merged.oe_logFC<-fc)).sum()
    down=((merged.kd_fdr<fdr)&(merged.oe_fdr<fdr)&(merged.kd_logFC<-fc)&(merged.oe_logFC>fc)).sum()
    rows.append({'FDR':fdr,'abs_logFC':fc,'loss_up':int(up),'loss_down':int(down),'total':int(up+down)})
sens=pd.DataFrame(rows); sens.to_csv(OUT/'PDE2A_reversible_threshold_sensitivity.tsv',sep='\t',index=False)

g178=pd.read_csv(PROJ/'GSE178537/GSE178537_IT_Edge_vs_Core_stats.tsv',sep='\t')
g198=pd.read_csv(PROJ/'GSE198315/GSE198315_malignant_edge_core_score_stats.tsv',sep='\t')
g208=pd.read_csv(PROJ/'GSE208253/GSE208253_LE_TC_sample_level_stats.tsv',sep='\t')

fig,ax=plt.subplots(2,2,figsize=(7.2,5.8),constrained_layout=True)
# A sensitivity
a=ax[0,0]; labels=[f"FDR<{r.FDR:g}, |LFC|>{r.abs_logFC:g}" for r in sens.itertuples()]
x=np.arange(len(labels)); w=.36
a.bar(x-w/2,sens.loss_up,w,color=BLUE,label='loss-up'); a.bar(x+w/2,sens.loss_down,w,color=ORANGE,label='loss-down')
a.set_xticks(x,labels,rotation=35,ha='right'); a.set_ylabel('Number of genes'); a.set_title('A  Reversible signature sensitivity'); a.legend(frameon=False)
# B GSE178537
a=ax[0,1]; sub=g178[g178.feature.isin(['PDE2A','PDE2A_loss_up_score','PDE2A_loss_down_score'])].copy(); x=np.arange(len(sub));
a.bar(x-.18,sub.edge_mean,.36,color=BLUE,label='Edge/IT'); a.bar(x+.18,sub.core_mean,.36,color=GREY,label='Core/VT'); a.set_xticks(x,['PDE2A','loss-up','loss-down'],rotation=25); a.set_ylabel('Mean score'); a.set_title('B  GSE178537 bulk projection'); a.legend(frameon=False)
# C GSE198315
a=ax[1,0]; sub=g198[g198.feature.isin(['PDE2A','PDE2A_loss_up_score','PDE2A_loss_down_score','stemness_score','MAPK_like_score'])].copy(); piv=sub.pivot(index='feature',columns='region',values='mean').reindex(['PDE2A','PDE2A_loss_up_score','PDE2A_loss_down_score','stemness_score','MAPK_like_score']); x=np.arange(len(piv));
a.bar(x-.18,piv['Edge'],.36,color=BLUE,label='Malignant Edge'); a.bar(x+.18,piv['Core'],.36,color=GREY,label='Malignant Core'); a.set_xticks(x,['PDE2A','loss-up','loss-down','stemness','MAPK-like'],rotation=30); a.set_ylabel('Mean score'); a.set_title('C  GSE198315 malignant cells'); a.legend(frameon=False)
# D spatial sample-level
a=ax[1,1]; sub=g208[g208.feature.isin(['PDE2A','PDE2A_loss_up_score','PDE2A_loss_down_score','proliferation_score','leading_edge_score'])].copy(); x=np.arange(len(sub));
a.bar(x-.18,sub.edge_mean,.36,color=BLUE,label='Edge/LE'); a.bar(x+.18,sub.core_mean,.36,color=GREY,label='Core/TC'); a.set_xticks(x,['PDE2A','loss-up','loss-down','prolif.','leading-edge'],rotation=30); a.set_ylabel('Sample-level mean'); a.set_title('D  GSE208253 spatial projection'); a.legend(frameon=False)
fig.suptitle('PDE2A-KD/OE reversible program: cross-cohort projection',fontsize=12,fontweight='bold')
fig.savefig(OUT/'Figure_PDE2A_reversible_projection.pdf'); fig.savefig(OUT/'Figure_PDE2A_reversible_projection.png'); plt.close(fig)
print('Wrote', OUT/'Figure_PDE2A_reversible_projection.pdf')
print(sens.to_string(index=False))
