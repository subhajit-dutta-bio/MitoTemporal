from pathlib import Path
import sys,json,itertools
root=Path(__file__).parent
sys.path.insert(0,str(root/'python_packages'))
import pandas as pd,numpy as np
from scipy.stats import wilcoxon
j=json.loads((root/'hummel_patient_metadata.json').read_text())
mods=list(j['gene_modules'])
base=pd.DataFrame(j['patient_scores'])
def bh(p):
 p=np.array(p);q=np.empty(len(p));o=np.argsort(p);q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
rng=np.random.default_rng(20261024)
def test(d):
 if len(d)<3:return np.nan
 return float(wilcoxon(d,zero_method='wilcox',alternative='two-sided',method='auto').pvalue) if np.any(d!=0) else 1.
def ci(d):
 if len(d)<3:return np.nan,np.nan
 x=rng.choice(d,size=(10000,len(d)),replace=True).mean(axis=1)
 return tuple(np.quantile(x,[.025,.975]))
def analyse(df,mincells=10):
 d=df[df.n_cells>=mincells];treatment=[];progression=[];coverage=[];triples=[]
 for (comp,ct),g in d.groupby(['compartment','cell_type']):
  for mod in mods:
   piv=g.pivot(index='patient',columns='timepoint',values='mito_'+mod)
   if {'Pre','Post'}.issubset(piv.columns):
    p=piv[['Pre','Post']].dropna();delta=(p.Post-p.Pre).to_numpy();lo,hi=ci(delta)
    treatment.append(dict(compartment=comp,cell_type=ct,module=mod,n_pairs=len(p),mean_delta=float(delta.mean()) if len(delta) else np.nan,p=test(delta),ci_low=lo,ci_high=hi,pct_positive=float((delta>0).mean()) if len(delta) else np.nan))
   if {'Pre','Post','Prog'}.issubset(piv.columns):
    p=piv[['Pre','Post','Prog']].dropna()
    if mod==mods[0]:coverage.append(dict(compartment=comp,cell_type=ct,n_triples=len(p),patients=';'.join(p.index)))
    for pat,row in p.iterrows():triples.append(dict(compartment=comp,cell_type=ct,module=mod,patient=pat,Pre=row.Pre,Post=row.Post,Prog=row.Prog))
    for a,b in [('Pre','Post'),('Post','Prog'),('Pre','Prog')]:
     delta=(p[b]-p[a]).to_numpy();lo,hi=ci(delta)
     progression.append(dict(compartment=comp,cell_type=ct,module=mod,contrast=b+'_minus_'+a,n_pairs=len(p),mean_delta=float(delta.mean()) if len(delta) else np.nan,p=test(delta),ci_low=lo,ci_high=hi,pct_positive=float((delta>0).mean()) if len(delta) else np.nan,mean_Pre=float(p.Pre.mean()),mean_Post=float(p.Post.mean()),mean_Prog=float(p.Prog.mean())))
 ts=pd.DataFrame(treatment);pr=pd.DataFrame(progression)
 for frame in [ts,pr]:
  frame['q']=np.nan
  for comp,idx in frame.groupby('compartment').groups.items():
   valid=frame.loc[idx,'p'].dropna().index
   if len(valid):frame.loc[valid,'q']=bh(frame.loc[valid,'p'])
 return ts,pr,pd.DataFrame(coverage),pd.DataFrame(triples)
ts,pr,cov,tr=analyse(base)
for name,d in [('treatment_patient_compartment',ts),('progression_patient_compartment',pr),('progression_coverage',cov),('progression_individual_trajectories',tr)]:d.to_csv(root/(name+'.csv'),index=False)
summary={'source':j['source'],'min_cells':10,'bootstrap_seed':20261024,'bootstrap_iterations':10000,'treatment':{},'progression':{}}
for name,frame in [('treatment',ts),('progression',pr)]:
 for comp,g in frame.groupby('compartment'):
  valid=g.dropna(subset=['p']);summary[name][comp]={'tests':len(valid),'n_range':[int(valid.n_pairs.min()),int(valid.n_pairs.max())] if len(valid) else None,'significant_q':int((valid.q<.05).sum()),'min_p':float(valid.p.min()) if len(valid) else None,'min_q':float(valid.q.min()) if len(valid) else None}
print('SUMMARY',json.dumps(summary))
print('COVERAGE\n',cov.to_string(index=False))
print('IV PROGRESSION\n',pr[(pr.module=='OXPHOS_Complex_IV')&(pr.contrast!='Post_minus_Pre')].to_string(index=False))
print('TOP TREATMENT\n',ts.sort_values('q').head(15).to_string(index=False))
sens=[]
for threshold,records in [('50_cells',None)]+list(j['qc_sensitivity'].items()):
 df=base if records is None else pd.DataFrame(records)
 s,p,c,t=analyse(df,mincells=50 if records is None else 10)
 for comp in ['blood','tissue']:
  v=s[s.compartment.eq(comp)];w=p[p.compartment.eq(comp)]
  sens.append(dict(analysis=threshold,compartment=comp,treatment_tests=int(v.p.notna().sum()),treatment_q_lt_005=int((v.q<.05).sum()),progression_tests=int(w.p.notna().sum()),progression_q_lt_005=int((w.q<.05).sum())))
 for name,d in [('treatment',s),('progression',p)]:d.to_csv(root/('sensitivity_'+threshold+'_'+name+'.csv'),index=False)
pd.DataFrame(sens).to_csv(root/'patient_sensitivity_summary.csv',index=False)
summary['sensitivity']=sens
(root/'patient_analysis_summary.json').write_text(json.dumps(summary,indent=2))
print('SENSITIVITY',json.dumps(sens))
