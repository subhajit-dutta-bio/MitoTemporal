from common import *
comparators=json.loads((r/'fixed_comparator_gene_sets.json').read_text())
groups={'Import_TIM_TOM':mods['Import_TIM_TOM'],**comparators}
def welch(x,y):
 x=pd.Series(x).dropna();y=pd.Series(y).dropna();rec=dict(n_treated=len(x),n_control=len(y),estimate=np.nan,low=np.nan,high=np.nan,p=np.nan,df=np.nan)
 if min(len(x),len(y))<2:return rec
 a=x.var(ddof=1)/len(x);b=y.var(ddof=1)/len(y);se=np.sqrt(a+b);d=x.mean()-y.mean()
 if se==0:
  rec['estimate']=d
  return rec
 df=(a+b)**2/(a*a/(len(x)-1)+b*b/(len(y)-1));ci=stats.t.ppf(.975,df)*se
 rec.update(estimate=d,low=d-ci,high=d+ci,p=2*stats.t.sf(abs(d/se),df),df=df);return rec
def onesample(d):
 d=pd.Series(d).dropna();rec=dict(n_treated=len(d),n_control=0,estimate=np.nan,low=np.nan,high=np.nan,p=np.nan,df=np.nan)
 if len(d)<2:return rec
 est=d.mean();se=d.std(ddof=1)/np.sqrt(len(d));ci=stats.t.ppf(.975,len(d)-1)*se;rec.update(estimate=est,low=est-ci,high=est+ci,p=2*stats.t.sf(abs(est/se),len(d)-1),df=len(d)-1);return rec
def family(df,by,expected=None):
 for _,ix in df.groupby(by,dropna=False).groups.items():
  if expected is not None:assert len(ix)==expected,(by,len(ix),expected)
  df.loc[ix,'q']=bh(df.loc[ix,'p'].fillna(1));df.loc[ix,'family_size']=len(ix);df.loc[df.p.isna(),'q']=np.nan
 return df
def paired_scores(a,pairs,gs=groups):
 rows=[]
 for v in pairs.itertuples():
  z=dict(unit=str(v.unit),baseline=v.baseline,on=v.on,arm=getattr(v,'arm','paired'))
  for name,genes in gs.items():
   hit=[g for g in genes if g in a.index and pd.notna(a.at[g,v.baseline]) and pd.notna(a.at[g,v.on])];z[name]=(a.loc[hit,v.on]-a.loc[hit,v.baseline]).mean() if len(hit)>=max(3,int(np.ceil(len(genes)/2))) else np.nan;z[name+'_n_genes']=len(hit)
  rows.append(z)
 return pd.DataFrame(rows)
def clinical(acc):
 if acc in ['GSE20181','GSE111563']:
  a=pd.read_csv(root/'revision7/analysis'/f'{acc}_full_gene_expression.csv.gz',index_col=0)
  p=pd.read_csv(root/'revision7/analysis'/f'{acc}_mechanism_pair_crosswalk.csv',dtype=str);pre,post=('A','C') if acc=='GSE20181' else ('1','4');p=p.rename(columns={'patient':'unit',pre:'baseline',post:'on'});p['arm']='AI'
 else:
  a=pd.read_csv(root/'revision8/analysis/GSE105777_full_observed_expression.csv.gz',index_col=0);p=pd.read_csv(root/'revision8/analysis/Ki67_patient_crosswalk.csv',dtype=str).query('accession=="GSE105777"').rename(columns={'clinical_key':'unit'})
 return a,p[['unit','baseline','on','arm']]
