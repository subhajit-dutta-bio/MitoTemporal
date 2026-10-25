from pathlib import Path
import sys,json,re,hashlib
r=Path(__file__).resolve().parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'));sys.path.insert(0,str(root/'revision10'))
import common as inherited
import pandas as pd,numpy as np,statsmodels.api as sm
from scipy import stats
canon=inherited.canon;approved=inherited.approved;mods=inherited.mods;controls=inherited.controls
comparators=json.loads((root/'revision10/fixed_comparator_gene_sets.json').read_text())
sets={'Import_TIM_TOM':mods['Import_TIM_TOM'],**comparators,'ER_programme':controls['ER_programme']}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);valid=np.flatnonzero(np.isfinite(p));v=p[valid];o=np.argsort(v);a=np.minimum(1,np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);q[valid[o]]=a;return q
def tidy_matrix(a):
 a=a.copy();a.index=[canon(str(g).split(' (')[0]) for g in a.index];a=a[a.index.isin(approved)];return a.groupby(level=0).median()
def module_score(a,genes,minimum=None):
 hit=[g for g in genes if g in a.index];minimum=minimum or max(3,int(np.ceil(len(genes)/2)))
 if len(hit)<minimum:return pd.Series(np.nan,index=a.columns)
 z=a.loc[hit];return z.mean().where(z.notna().sum()>=max(minimum,int(np.ceil(.8*len(hit)))))
def namekey(s):return re.sub('[^A-Z0-9]','',str(s).upper())
def fit_ols(df,y,predictor,covariates=(),lineage=False,standardize=True,min_n=10):
 columns=[y,predictor,*covariates]+(['lineage'] if lineage else [])
 z=df.loc[:,list(dict.fromkeys(columns))].replace([np.inf,-np.inf],np.nan).dropna().copy()
 if len(z)<min_n:return dict(n=len(z),status='unavailable',reason='below prespecified n gate')
 X=z[[predictor,*covariates]].astype(float)
 if standardize:
  sd=X.std(ddof=1)
  if (sd<1e-12).any():return dict(n=len(z),status='unavailable',reason='constant covariate/predictor')
  X=(X-X.mean())/sd
 if lineage:
  counts=z.lineage.value_counts();labels=z.lineage.where(z.lineage.map(counts)>=10,'Other');X=pd.concat([X,pd.get_dummies(labels,prefix='lineage',drop_first=True,dtype=float)],axis=1)
 X=sm.add_constant(X,has_constant='add').astype(float)
 rank=np.linalg.matrix_rank(X);cond=np.linalg.cond(X)
 if rank<X.shape[1]:return dict(n=len(z),status='unavailable',reason='rank-deficient fixed model',rank=int(rank),parameters=X.shape[1])
 model=sm.OLS(z[y].to_numpy(float),X).fit(cov_type='HC3',use_t=True)
 low,high=model.conf_int().loc[predictor];return dict(n=len(z),estimate=model.params[predictor],low=low,high=high,p=model.pvalues[predictor],df=model.df_resid,parameters=X.shape[1],rank=int(rank),condition_number=cond,status='tested')
def frozen_family(table,by):
 for _,ix in table.groupby(by,dropna=False).groups.items():table.loc[ix,'q']=bh(table.loc[ix,'p']);table.loc[ix,'family_slots']=len(ix);table.loc[ix,'family_tested']=table.loc[ix,'p'].notna().sum()
 return table
