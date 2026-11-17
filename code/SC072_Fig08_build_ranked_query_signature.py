from pathlib import Path
import pandas as pd
base=Path('/data/liaolz/PDE2A')
x=pd.read_csv(base/'results/PDE2A_ranked_signatures/PDE2A_ranked_reversible_all_genes.tsv',sep='\t')
ann=pd.read_csv(base/'results/PDE2A_reversible_enrichment/PDE2A_loss_up_module_annotation.tsv',sep='\t')
drop=set(ann.loc[ann.module=='inflammatory','gene'])
drop |= {'IGH','IGKC','IGLC','MS4A1','MZB1','JCHAIN','CD3D','CD3E','LYZ','S100A8','S100A9','FCGR3B','COL1A1','COL1A2','PECAM1','VWF'}
keep=~x.gene.astype(str).str.upper().str.startswith(('IGH','IGK','IGL','TRAV','TRBV','TRGV','TRDV'))
keep &= ~x.gene.isin(drop)
y=x[keep].copy()
out=base/'results/PDE2A_ranked_signatures'
y.to_csv(out/'PDE2A_clean_ranked_reversible_all_genes.tsv',sep='\t',index=False)
for n in [50,100,250]:
    y.head(n).gene.to_csv(out/f'PDE2A_clean_loss_up_top{n}.txt',index=False,header=False)
    y.tail(n).sort_values('reversible_score').gene.to_csv(out/f'PDE2A_clean_loss_down_top{n}.txt',index=False,header=False)
print('removed',len(x)-len(y),'remaining',len(y)); print(y.head(10)[['gene','reversible_score']].to_string(index=False))
