from pathlib import Path
import sys,json
root=Path(__file__).parent
sys.path.insert(0,str(root.parent/'python_packages'))
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
out=root/'analysis';out.mkdir(exist_ok=True)
def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
rng=np.random.default_rng(20261024)
rows=[];pair_records=[];coverage=[];library_samples=[]
def analyse(acc,d,name='primary',score='CIV',mt=20,mincells=10,demux='hashsolo',gate='strict',exclude_s8=False):
 a=d.copy();a=a[a.pct_mt<mt]
 if demux=='hashsolo':a=a[a.hashsolo_pass.astype(str).str.lower().eq('true')]
 else:
  ratio=5 if demux=='dom5' else 10
  a=a[a.hash_ratio.ge(ratio)&a.hash_top_count.ge(10)];a['patient']=a.patient_dom;a['timepoint']=a.timepoint_dom
 a=a[a.patient.ne('Unassigned')]
 if acc=='GSE205589':a=a[a.patient.ne('P4')]
 if exclude_s8:a=a[a.library.ne('S8')]
 if gate=='strict':a=a[a.lineage.isin(['CD4','CD8'])]
 elif gate=='CD8B':a=a[a.CD8B_positive.astype(str).str.lower().eq('true')];a['lineage']='CD8B'
 else:a=a[a.T_marker_positive.astype(str).str.lower().eq('true')];a['lineage']='All_T'
 agg=a.groupby(['patient','timepoint','lineage'],observed=True)[score].agg(['mean','count']).reset_index()
 agg['accession']=acc;agg['analysis']=name
 if name=='primary':
  agg.to_csv(out/(acc+'_patient_timepoint_scores.csv'),index=False)
  library_samples.extend(a.groupby(['library','patient','timepoint','lineage'],observed=True).size().rename('n_cells').reset_index().assign(accession=acc).to_dict('records'))
 pre,post=('Pre','Post') if acc=='GSE205589' else ('Week1','Week15')
 for lineage,g in agg.groupby('lineage'):
  g=g[g['count']>=mincells];p=g.pivot(index='patient',columns='timepoint',values='mean')
  if pre not in p or post not in p:continue
  pairs=p[[pre,post]].dropna();delta=(pairs[post]-pairs[pre]).to_numpy();n=len(delta)
  pval=float(wilcoxon(delta,alternative='two-sided',method='auto').pvalue) if n>=3 and np.any(delta) else (1. if n>=3 else np.nan)
  ci=np.quantile(rng.choice(delta,size=(10000,n),replace=True).mean(axis=1),[.025,.975]) if n>=3 else [np.nan,np.nan]
  row={'accession':acc,'analysis':name,'lineage':lineage,'contrast':post+'_minus_'+pre,'endpoint':'scheduled_on_treatment_not_progression','score':score,'n_pairs':n,'mean_pre':pairs[pre].mean(),'mean_post':pairs[post].mean(),'mean_delta':delta.mean() if n else np.nan,'median_delta':np.median(delta) if n else np.nan,'n_positive':int((delta>0).sum()),'n_negative':int((delta<0).sum()),'p':pval,'ci_low':ci[0],'ci_high':ci[1]}
  rows.append(row)
  if name=='primary':
   for patient,r in pairs.iterrows():pair_records.append({'accession':acc,'lineage':lineage,'patient':patient,'pre':r[pre],'post':r[post],'delta':r[post]-r[pre]})
for acc in ['GSE205589','GSE300475']:
 path=root/'results'/acc/'cell_scores.csv.gz'
 if not path.exists():print('MISSING',acc);continue
 d=pd.read_csv(path,index_col=0)
 analyse(acc,d)
 for name,kwargs in [('min50',{'mincells':50}),('mt10',{'mt':10}),('mt40',{'mt':40}),('nuclear',{'score':'CIV_nuclear'}),('mean_log',{'score':'CIV_mean_log'}),('nuclear_mean_log',{'score':'CIV_nuclear_mean_log'}),('demux_dom5',{'demux':'dom5'}),('demux_dom10',{'demux':'dom10'}),('CD8B_gate',{'gate':'CD8B'}),('all_T',{'gate':'all_T'})]:analyse(acc,d,name=name,**kwargs)
 if acc=='GSE300475':analyse(acc,d,name='exclude_PT15_rep_S8',exclude_s8=True)
 a=d[d.hashsolo_pass.astype(str).str.lower().eq('true')&d.pct_mt.lt(20)]
 coverage.extend(a.groupby(['library','patient','timepoint','lineage'],observed=True).size().rename('n_cells').reset_index().assign(accession=acc).to_dict('records'))
 if acc=='GSE300475':
  x=set(i.split('|')[-1] for i in a[a.library.eq('S5')].index);y=set(i.split('|')[-1] for i in a[a.library.eq('S8')].index)
  print('PT15_REPLICATE_BARCODE_OVERLAP',len(x),len(y),len(x&y))
res=pd.DataFrame(rows);res['q']=np.nan
for name,idx in res.groupby('analysis').groups.items():
 valid=res.loc[idx,'p'].dropna().index
 if len(valid):res.loc[valid,'q']=bh(res.loc[valid,'p'])
res.to_csv(out/'replication_all_results.csv',index=False)
res[res.analysis.eq('primary')].to_csv(out/'replication_primary_results.csv',index=False)
pd.DataFrame(pair_records).to_csv(out/'replication_patient_deltas.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'replication_cell_coverage.csv',index=False)
pd.DataFrame(library_samples).to_csv(out/'replication_marker_subset_coverage.csv',index=False)
print('PRIMARY\n'+res[res.analysis.eq('primary')].to_string(index=False))
print('SENSITIVITY\n'+res[['accession','analysis','lineage','n_pairs','mean_delta','n_positive','p','q']].to_string(index=False))
