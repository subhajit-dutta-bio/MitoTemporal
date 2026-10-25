from common import *
from scipy.stats import hypergeom
raw=pd.read_csv(r/'metadata/GSE27463_RefSeq.all.tsv.gz',sep='\t');raw['gene']=raw.GeneSymbol.map(canon);cols=['VEH','E2_10m','E2_40m','E2_160m'];a=raw.groupby('gene')[cols].median();a.to_csv(out/'GROseq_gene_profiles.csv.gz',compression='gzip')
# Deposited values have no replicate columns; temporal module profiles are descriptive.
rows=[]
for mod,gs in {k:mods[k] for k in ['Import_TIM_TOM','CIV_nuclear']}.items():
 hit=sorted(set(gs)&set(a.index));z=np.log2((a.loc[hit]+1).div(a.loc[hit,'VEH']+1,axis=0))
 for time in cols:rows.append(dict(programme=mod,time=time,n_genes=len(hit),mean_log2_response=z[time].mean(),median_log2_response=z[time].median()))
for g in ['NRF1','TFAM','PPARGC1A']:
 if g in a.index:
  for time in cols:rows.append(dict(programme=g,time=time,n_genes=1,mean_log2_response=np.log2((a.loc[g,time]+1)/(a.loc[g,'VEH']+1)),median_log2_response=np.nan))
pd.DataFrame(rows).to_csv(out/'GROseq_descriptive_timecourse.csv',index=False)
bed=pd.read_csv(r/'metadata/NRF1_MCF7_GRCh38_optimal.bed.gz',sep='\t',header=None);peaks={c:(np.sort(z[1].to_numpy()),np.sort(z[2].to_numpy())) for c,z in bed.groupby(0)};tss=[]
with gzip.open(r/'metadata/gencode.v47.annotation.gtf.gz','rt') as f:
 for l in f:
  if l.startswith('#'):continue
  z=l.rstrip().split('\t')
  if z[2]!='transcript':continue
  g=canon(re.search(r'gene_name "([^"]+)"',z[8])[1]);pos=int(z[3])-1 if z[6]=='+' else int(z[4])-1;start=pos-1000;end=pos+1001
  bound=False
  if z[0] in peaks:
   starts,ends=peaks[z[0]];bound=np.searchsorted(starts,end,side='left')-np.searchsorted(ends,start,side='right')>0
  tss.append(dict(gene=g,chrom=z[0],tss=pos,strand=z[6],NRF1_bound=bool(bound)))
tss=pd.DataFrame(tss).drop_duplicates();tss.to_csv(out/'NRF1_GRCh38_TSS_binding.csv.gz',index=False,compression='gzip');bound=tss.groupby('gene').NRF1_bound.any();expressed=set(a.index[a.VEH>=1])&set(bound.index);universe=sorted(expressed);flags=bound.reindex(universe).astype(int);rng=np.random.default_rng(20261024);bins=pd.qcut(a.loc[universe,'VEH'].rank(method='first'),10,labels=False);pools={b:np.array(bins.index[bins==b]) for b in sorted(bins.unique())};results=[];geneout=[]
for mod in ['Import_TIM_TOM','CIV_nuclear']:
 hit=sorted(set(mods[mod])&expressed);k=int(bound.loc[hit].sum());N=len(universe);K=int(flags.sum());n=len(hit);pv=hypergeom.sf(k-1,N,K,n);sim=np.zeros(10000)
 for g in hit:sim+=bound.loc[rng.choice(pools[bins[g]],10000,replace=True)].to_numpy(dtype=int)
 results.append(dict(module=mod,universe='MCF7_GROseq_VEH_at_least_1_and_GENCODE47_TSS',N=N,K=K,n_module_expressed=n,n_bound=k,expected=n*K/N,fold_enrichment=(k/n)/(K/N),p=pv,expression_matched_permutation_p=(1+np.sum(sim>=k))/10001,n_permutations=10000))
 geneout.extend(dict(module=mod,gene=g,baseline_GRO_value=a.loc[g,'VEH'] if g in a.index else np.nan,expressed=g in expressed,NRF1_bound=bool(bound[g]) if g in bound.index else False) for g in mods[mod])
z=pd.DataFrame(results);z['q']=bh(z.p);z['matched_q']=bh(z.expression_matched_permutation_p);z.to_csv(out/'NRF1_binding_enrichment.csv',index=False);pd.DataFrame(geneout).to_csv(out/'NRF1_module_gene_binding.csv',index=False);print(z.to_string(index=False));print(pd.DataFrame(rows).to_string(index=False))
