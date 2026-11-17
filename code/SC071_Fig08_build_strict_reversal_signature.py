from pathlib import Path
import pandas as pd
base=Path('/data/liaolz/PDE2A')
sigdir=base/'results/reversible_signature'
kd=pd.read_csv(sigdir/'PDE2A_KD_vs_Control_DEG.tsv',sep='\t')
oe=pd.read_csv(sigdir/'PDE2A_OE_vs_Control_DEG.tsv',sep='\t')
x=kd[['gene','logFC','adj.P.Val']].rename(columns={'logFC':'kd_logFC','adj.P.Val':'kd_fdr'}).merge(oe[['gene','logFC','adj.P.Val']].rename(columns={'logFC':'oe_logFC','adj.P.Val':'oe_fdr'}),on='gene')
x=x[(x.kd_fdr<0.1)&(x.oe_fdr<0.1)&(((x.kd_logFC>0.5)&(x.oe_logFC< -0.5))|((x.kd_logFC< -0.5)&(x.oe_logFC>0.5)))].copy()
x['loss_class']=x.apply(lambda r:'loss_up' if r.kd_logFC>0 else 'loss_down',axis=1)
x.to_csv(sigdir/'PDE2A_KD_OE_reversible_signature_strict.tsv',sep='\t',index=False)
print(x.loss_class.value_counts().to_string())
