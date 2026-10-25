from pathlib import Path
import sys,json,re,io,gzip,csv,hashlib
r=Path(__file__).resolve().parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'))
import pandas as pd,numpy as np
import statsmodels.api as sm
from scipy import stats
h=pd.read_csv(root/'revision6/metadata/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False)
h=h[h.status=='Approved'];approved=set(h.symbol);aliases={}
for row in h.itertuples():
 for col in ['prev_symbol','alias_symbol']:
  for g in getattr(row,col).split('|'):
   if g:aliases.setdefault(g.upper(),set()).add(row.symbol)
def canon(g):
 g=str(g).upper()
 return g if g in approved else next(iter(aliases[g])) if g in aliases and len(aliases[g])==1 else g
mods=json.loads((root/'revision6/analysis/canonical_modules.json').read_text());mods['CIV_nuclear']=[g for g in mods['OXPHOS_Complex_IV'] if not g.startswith('MT-')]
controls=json.loads((root/'revision4/analysis/control_gene_sets.json').read_text());controls={k:[canon(g) for g in v] for k,v in controls.items()}
sets=dict(mods,**controls,NRF1=['NRF1'],TFAM=['TFAM'],ESR1=['ESR1'],PGR=['PGR'],GREB1=['GREB1'])
assert not set(mods['Import_TIM_TOM']) & set(controls['Proliferation'])
assert not set(mods['Import_TIM_TOM']) & set(controls['ER_programme'])
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def canonicalize(a):
 a=a.copy();a.index=[canon(g) for g in a.index];return a.groupby(level=0).median()
def scores(a,groups=sets):
 s=pd.DataFrame(index=a.columns)
 for name,gs in groups.items():
  hit=[g for g in gs if g in a.index];minimum=1 if len(gs)==1 else max(3,int(np.ceil(len(gs)/2)))
  s[name]=a.loc[hit].mean() if len(hit)>=minimum else np.nan
  if len(hit):s.loc[a.loc[hit].notna().sum()<minimum,name]=np.nan
 return s
def metadata(acc):
 values={}
 with gzip.open(r/'metadata'/f'{acc}_series_matrix.txt.gz','rt') as f:
  for l in f:
   if l.startswith('!Sample_'):
    row=next(csv.reader([l],delimiter='\t'));k=row[0][8:]
    if k in ['geo_accession','title']:values[k]=row[1:]
 return pd.DataFrame(values).rename(columns={'geo_accession':'sample'}).set_index('sample')
def tdiff(d):
 d=np.asarray(d,float);n=len(d);se=d.std(ddof=1)/np.sqrt(n);delta=d.mean();ci=stats.t.ppf(.975,n-1)*se
 return dict(n=n,estimate=delta,low=delta-ci,high=delta+ci,p=stats.ttest_1samp(d,0).pvalue)
def freeze(a,name):a.to_csv(out/name,index=False)
