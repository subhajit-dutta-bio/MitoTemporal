from pathlib import Path
import sys,json,ast,io,re,gzip,csv
r=Path(__file__).parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'))
import pandas as pd,numpy as np
h=pd.read_csv(root/'revision6/metadata/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);h=h[h.status=='Approved'];approved=set(h.symbol.str.upper());prev={};alias={}
for row in h.itertuples():
 for col,d in [('prev_symbol',prev),('alias_symbol',alias)]:
  for t in getattr(row,col).split('|'):
   if t:d.setdefault(t.upper(),set()).add(row.symbol.upper())
def canon(g):
 g=str(g).upper()
 if g in approved:return g
 for d in [prev,alias]:
  if g in d:return next(iter(d[g])) if len(d[g])==1 else g
 return g
def canonicalize(a):
 a=a.copy();a.index=[canon(g) for g in a.index];return a.groupby(level=0).median()
mods=json.loads((root/'revision6/analysis/canonical_modules.json').read_text());mods['CIV_nuclear']=[g for g in mods['OXPHOS_Complex_IV'] if not g.startswith('MT-')]
controls=json.loads((root/'revision4/analysis/control_gene_sets.json').read_text());controls={k:[canon(g) for g in gs] for k,gs in controls.items()}
def bh(p):
 p=np.array(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def matrix(acc):
 path=r/'metadata'/f'{acc}_series_matrix.txt.gz';meta={};ch={};lines=[];inside=False
 with gzip.open(path,'rt') as f:
  for l in f:
   if l.startswith('!series_matrix_table_begin'):inside=True
   elif l.startswith('!series_matrix_table_end'):break
   elif inside:lines.append(l)
   elif l.startswith('!Sample_'):
    v=next(csv.reader([l],delimiter='\t'));key=v[0][8:]
    if key.startswith('characteristics'):
     k=v[1].split(': ',1)[0];ch[k]=[z.split(': ',1)[-1] for z in v[1:]]
    else:meta[key]=v[1:]
 a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0);a.index=a.index.astype(str)
 md=pd.DataFrame(ch,index=meta['geo_accession']);md.index.name='sample';md['title']=meta['title'];md.to_csv(out/f'{acc}_sample_metadata.csv');return a,md,meta
def oldmatrix(acc):
 cache=out/f'{acc}_full_gene_expression.pkl'
 if cache.exists():return pd.read_pickle(cache)
 frozen=out/f'{acc}_full_gene_expression.csv.gz'
 if frozen.exists():return pd.read_csv(frozen,index_col=0,float_precision='round_trip')
 if acc in ['GSE153470','GSE105777']:a=pd.read_pickle(root/'revision4/analysis'/f'{acc}_gene_expression.pkl')
 else:
  src=(root/'revision2/prepare_scores.py').read_text();nodes=[v for v in ast.parse(src).body if isinstance(v,ast.FunctionDef) and v.name in ['symbol','soft']];ns=dict(pd=pd,np=np,io=io,re=re,r=root/'revision2',out=out)
  exec(compile(ast.Module(body=nodes,type_ignores=[]),'preserved_parser','exec'),ns);a,_,_=ns['soft'](acc)
  if acc=='GSE20181':a=np.log2(a.clip(lower=.01))
 a=canonicalize(a);a.to_pickle(cache);return a
def scores(a,sets=None):
 sets=sets or dict(mods,**controls,**{g:[g] for g in ['NRF1','TFAM','PPARGC1A']})
 result=pd.DataFrame(index=a.columns)
 for name,gs in sets.items():
  hit=[g for g in gs if g in a.index]
  if len(hit)>=max(1 if len(gs)==1 else 3,np.ceil(len(gs)/2)):result[name]=a.loc[hit].mean()
 return result
