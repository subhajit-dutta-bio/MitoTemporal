from common import *
def estimate(a,canonical=True):
 cg=pd.read_csv(r/'metadata/estimate_common_genes.txt',sep='\t',dtype=str);gs={}
 for l in (r/'metadata/estimate_SI_geneset.gmt').read_text().splitlines():
  z=l.split('\t');gs[z[0]]=[canon(g) if canonical else g for g in z[2:]]
 # R merge sorts its character key; preserve original input in reference validation.
 names=[canon(g) if canonical else g for g in cg.GeneSymbol];names=sorted(set(names)&set(a.index));a=a.loc[names].dropna(axis=0,how='any');rr=a.rank(axis=0,method='average').to_numpy()*10000/len(a);result={}
 for key,genes in gs.items():
  isin=np.isin(a.index,genes);v=[]
  for j in range(a.shape[1]):
   order=np.argsort(-rr[:,j],kind='stable');hit=isin[order];w=rr[order,j]**.25;fn=np.cumsum(hit*w)/np.sum(hit*w);f0=np.cumsum(~hit)/np.sum(~hit);v.append(np.sum(fn-f0))
  result[key]=v
 z=pd.DataFrame(result,index=a.columns);z['ESTIMATEScore']=z.iloc[:,:2].sum(axis=1);z.attrs['common_genes']=len(a);return z
if __name__=='__main__':
 x=pd.read_csv(r/'metadata/estimate_sample_input.gct',sep='\t',skiprows=2,index_col=0).drop(columns='Description')
 y=pd.read_csv(r/'metadata/estimate_sample_estimate.gct',sep='\t',skiprows=2,index_col=0).drop(columns='Description').T
 z=estimate(x,False);error=float(np.max(np.abs(z.to_numpy()-y.to_numpy())));print(z);print('MAX_ERROR',error)
 (out/'ESTIMATE_validation.json').write_text(json.dumps(dict(max_absolute_error=error,n_samples=len(z),reference='Original ESTIMATE package bundled example',canonical_annotation_used_for_reference=False),indent=2));assert error<.01
