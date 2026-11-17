"""Figure 2F: parallel coordinates from the frozen 35-candidate Pareto table."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/publication_Figure2EF_ComplexUpset_Pareto_20261016'
OUT.mkdir(exist_ok=True)
SRC = ROOT / 'results/11_candidate_prioritization'
front = pd.read_csv(SRC/'CANDIDATE_PARETO_FRONTS.tsv',sep='\t')
lopo = pd.read_csv(SRC/'LOPO_CANDIDATE_STABILITY.tsv',sep='\t')
sens = pd.read_csv(SRC/'PARETO_SENSITIVITY.tsv',sep='\t')
assert len(front)==35 and front.gene.is_unique
first = ['CCL21','ADCY4','EREG','CAMK2B','PDE2A']
assert set(front.loc[front.pareto_front==1,'gene'])==set(first)
assert int(sens.set_index('gene').loc['PDE2A','pareto_front1_frequency'])==4
assert (front.criterion_edge_down_effect + front.full_shrunken_LFC).abs().max()<1e-10
for raw,criterion in [('LOPO_negative_count','criterion_LOPO_robustness'),
                      ('edge_down_leading_edge_count','criterion_leading_edge_recurrence'),
                      ('functional_module_count','criterion_module_breadth'),
                      ('database_support_count','criterion_database_breadth')]:
    assert (front[raw]-front[criterion]).abs().max()<1e-10
df=front.merge(lopo[['gene','LOPO_available_count']],on='gene',validate='one_to_one')
assert len(df)==35
assert df.set_index('gene').loc['WNK4','LOPO_available_count']==4
assert df.set_index('gene').loc['WNK4','LOPO_negative_count']==4
raw_cols=['full_shrunken_LFC','LOPO_negative_count','edge_down_leading_edge_count',
          'functional_module_count','database_support_count']
benefit=pd.DataFrame({
    'edge_down_effect':-df.full_shrunken_LFC,
    'LOPO_negative_count':df.LOPO_negative_count,
    'Core_enriched_memberships':df.edge_down_leading_edge_count,
    'functional_module_breadth':df.functional_module_count,
    'database_coverage':df.database_support_count,
})
benefit_cols=list(benefit)
minima=benefit.min(); maxima=benefit.max()
assert (maxima>minima).all()
normalized=(benefit-minima)/(maxima-minima)
assert normalized.notna().all().all()
for col in benefit_cols:
    df['plot_0_to_1_'+col]=normalized[col]
df.to_csv(OUT/'Figure2F_complete_35_candidate_parallel_source.tsv',sep='\t',index=False,na_rep='NA')
front.to_csv(OUT/'Figure2F_frozen_Pareto_front_assignments.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame({'axis':benefit_cols,'favorable_raw_min':minima.values,
              'favorable_raw_max':maxima.values,
              'original_metric':raw_cols}).to_csv(OUT/'Figure2F_display_normalization.tsv',sep='\t',index=False)

font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
font_manager.fontManager.addfont('C:/Windows/Fonts/arialbd.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none',
                     'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
fig,ax=plt.subplots(figsize=(8.1,4.35))
fig.subplots_adjust(left=.095,right=.85,top=.79,bottom=.24)
ax.set_xlim(-.17,4.75);ax.set_ylim(-.03,1.03)
ax.axis('off')
x=np.arange(5)
for j in x:
    ax.plot([j,j],[0,1],color='#86929C',lw=.65,zorder=0)
    for y in [0,.5,1]:ax.plot([j-.028,j+.028],[y,y],color='#86929C',lw=.55,zorder=1)
    raw=front[raw_cols[j]]
    # All axis labels are on the original measurement scale.
    low,high=float(minima.iloc[j]),float(maxima.iloc[j])
    values=[low,(low+high)/2,high]
    if j==0:values=[-v for v in values]
    for y,v in zip([0,.5,1],values):
        label=f'{v:.3f}' if j==0 else (f'{v:g}')
        ax.text(j-.045,y,label,ha='right',va='center',fontsize=6.3,color='#4D5A65')
    title=['Shrunken log2FC','Negative LOPO\ncount / 5*','Core-enriched\npathway memberships',
           'Functional\nmodule breadth','Database\ncoverage / 4'][j]
    ax.text(j,1.12,title,ha='center',va='bottom',fontsize=7,weight='bold')
    ax.text(j,1.06,'more negative favorable' if j==0 else 'higher favorable',
            ha='center',va='bottom',fontsize=6.2,color='#64727D')
cols=['plot_0_to_1_'+c for c in benefit_cols]
for _,r in df.loc[df.pareto_front!=1].iterrows():
    ax.plot(x,r[cols].to_numpy(dtype=float),color='#A7B1B8',lw=.55,alpha=.38,zorder=2)
colors={'CCL21':'#596C7A','ADCY4':'#668789','EREG':'#A17870',
        'CAMK2B':'#829071','PDE2A':'#7C3AED'}
for gene in first:
    r=df.set_index('gene').loc[gene]
    ax.plot(x,r[cols].to_numpy(dtype=float),color=colors[gene],
            lw=1.9 if gene=='PDE2A' else 1.35,zorder=4,
            marker='o',markersize=2.3 if gene=='PDE2A' else 1.7)
# Direct labels at final axis; deterministic spacing separates equal endpoints.
endpoints={g:float(df.set_index('gene').loc[g,cols[-1]]) for g in first}
groups={}
for g,y in endpoints.items():groups.setdefault(y,[]).append(g)
for y,group in groups.items():
    label_positions=(np.linspace(.78,.96,len(group)) if y>.85 else
                     np.linspace(max(.04,y-.08),min(.96,y+.08),len(group)))
    for yy,g in zip(label_positions,group):
        ax.plot([4.01,4.18],[y,yy],color=colors[g],lw=.6)
        ax.text(4.22,yy,g,ha='left',va='center',fontsize=7,
                weight='bold' if g=='PDE2A' else 'normal',color=colors[g])
fig.text(.035,.945,'F',fontsize=9,weight='bold')
fig.text(.075,.945,'Candidate prioritization · parallel coordinates',fontsize=8,weight='bold')
fig.text(.075,.888,'35 candidates  |  5 nondominated first-front candidates',fontsize=7)
fig.text(.095,.135,'Axes are min–max scaled to 0–1 across all 35 candidates; tick labels show original values.',fontsize=6.7)
fig.text(.095,.105,'* WNK4: 4/4 available LOPO effects negative; the fifth estimate was unavailable.',fontsize=6.7)
fig.text(.095,.075,'First-front membership is the frozen Pareto assignment; line crossings show trade-offs, not a new ranking.',fontsize=6.7)
for ext in ['svg','pdf','png']:
    fig.savefig(OUT/f'Figure2F_Pareto_ParallelCoordinates.{ext}',dpi=300,facecolor='white')
plt.close(fig)
print('F: 35 candidates, 5 front-1; frozen criteria and direction checks passed.')
print('Observed benefit ranges:',list(zip(benefit_cols,minima,maxima)))
