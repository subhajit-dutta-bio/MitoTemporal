from pathlib import Path
import sys,json
r=Path(__file__).parent;sys.path.insert(0,str(r.parent/'python_packages'))
import pandas as pd,numpy as np
from scipy.stats import wilcoxon
gene_file=r.parent/'source/hummel2/gene_sets/mito_modules.json'
if not gene_file.exists():gene_file=r/'discovery/mito_modules.json'
mods=json.loads(gene_file.read_text())
manifest=pd.read_csv(r/'analysis/additional_clinical_sample_manifest.csv',dtype={'patient':str,'timepoint':str})
gene_tables={acc:pd.read_csv(r/'analysis/gene_expression'/f'{acc}_module_genes.csv.gz',index_col=0) for acc in ['GSE20181','GSE111563','GSE147271']}
rows=[];gene_rows=[]
for acc,pre,post,mod in [('GSE20181','A','C','Import_TIM_TOM'),('GSE111563','1','4','Import_TIM_TOM'),('GSE147271','Pretreatment','Tamoxifen','Mitophagy')]:
 a=gene_tables[acc];gs=[g for g in mods[mod] if g in a.index];md=manifest[manifest.accession==acc].set_index('sample')
 common=[g for g in mods[mod] if all(g in gene_tables[k].index for k in ['GSE20181','GSE111563'])] if mod=='Import_TIM_TOM' else gs
 for spec,genes in [('all_measured',gs),('common_letrozole_genes',common)]+[('omit_'+g,[v for v in gs if v!=g]) for g in gs]:
  z=a.loc[genes].mean().rename('score').to_frame().join(md[['patient','timepoint']]);pv=z.groupby(['patient','timepoint']).score.mean().unstack()[[pre,post]].dropna();d=pv[post]-pv[pre]
  rows.append({'accession':acc,'module':mod,'sensitivity':spec,'genes':';'.join(genes),'n_genes':len(genes),'n_pairs':len(d),'mean_delta':d.mean(),'p_unadjusted':wilcoxon(d).pvalue})
 for g in gs:
  z=a.loc[g].rename('score').to_frame().join(md[['patient','timepoint']]);pv=z.groupby(['patient','timepoint']).score.mean().unstack()[[pre,post]].dropna();d=pv[post]-pv[pre]
  gene_rows.append({'accession':acc,'module':mod,'gene':g,'n_pairs':len(d),'mean_delta':d.mean(),'n_positive':(d>0).sum(),'n_negative':(d<0).sum()})
pd.DataFrame(rows).to_csv(r/'analysis/exploratory_gene_robustness.csv',index=False)
pd.DataFrame(gene_rows).to_csv(r/'analysis/exploratory_gene_directions.csv',index=False)
print(pd.DataFrame(rows).groupby(['accession','module']).agg(n_checks=('mean_delta','size'),min_delta=('mean_delta','min'),max_delta=('mean_delta','max')).to_string())
print(pd.DataFrame(rows).query('sensitivity=="common_letrozole_genes"').to_string(index=False))
print(pd.DataFrame(gene_rows).groupby(['accession','module']).mean_delta.apply(lambda x:f'{(x<0).sum()} down / {(x>0).sum()} up').to_string())
