from v14_common import *
import gzip,bisect,collections
raw=r/'raw';diagnostics=[]
# A small, transparent UCSC-chain endpoint mapper. Chain target is hg19 input;
# query is hg38 output; nonunique/haplotype mappings are rejected, never chosen.
blocks=collections.defaultdict(list)
with gzip.open(raw/'hg19ToHg38.over.chain.gz','rt') as stream:
 for line in stream:
  if line.startswith('chain '):
   z=line.split();tc,ts,te,qc,qs,qstrand,qpos,cid=z[2],int(z[5]),int(z[6]),z[7],int(z[8]),z[9],int(z[10]),z[12];tpos=ts
  elif line.strip():
   q=list(map(int,line.split()));size=q[0]
   if tc in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY'] and qc in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:blocks[tc].append((tpos,tpos+size,qc,qpos,qs,qstrand,cid))
   if len(q)==3:tpos+=size+q[1];qpos+=size+q[2]
for chrom in blocks:blocks[chrom].sort()
starts={c:np.array([v[0] for v in bs]) for c,bs in blocks.items()};endsmax={c:np.maximum.accumulate([v[1] for v in bs]) for c,bs in blocks.items()}
def mapbase(chrom,pos):
 if chrom not in starts:return None
 bs=blocks[chrom];stop=np.searchsorted(starts[chrom],pos,side='right');i=stop-1;hits=set()
 while i>=0 and endsmax[chrom][i]>pos:
  a,b,qc,qp,qs,strand,cid=bs[i]
  if a<=pos<b:
   p=qp+pos-a;p=p if strand=='+' else qs-p-1;hits.add((qc,p,strand,cid))
  i-=1
 return next(iter(hits)) if len(hits)==1 else None
def lift_peaks(path):
 a=pd.read_csv(path,sep='\t',header=None,usecols=[0,1,2]);mapped=[];reason=collections.Counter()
 for chrom,start,end in a.itertuples(index=False,name=None):
  x=mapbase(chrom,int(start));y=mapbase(chrom,int(end)-1)
  if x is None or y is None:reason['ambiguous_or_unmapped_endpoint']+=1;continue
  if x[0]!=y[0] or x[2]!=y[2] or x[3]!=y[3]:reason['inconsistent_chain_or_orientation']+=1;continue
  lo,hi=sorted([x[1],y[1]]);hi+=1
  if not .8<=((hi-lo)/(end-start))<=1.2:reason['length_discrepancy']+=1;continue
  mapped.append((x[0],lo,hi))
 diagnostics.append(dict(file=path.name,input_peaks=len(a),mapped_peaks=len(mapped),**reason));return pd.DataFrame(mapped,columns=['chrom','start','end']).drop_duplicates()
peaks={}
for sample,condition,rep in [('GSM3595410','E2',1),('GSM3595411','Vehicle',1),('GSM3595412','E2',2),('GSM3595413','Vehicle',2)]:
 path=next(raw.glob(sample+'*.narrowPeak.gz'));peaks[(condition,rep)]=lift_peaks(path);peaks[(condition,rep)].to_csv(out/f'ER_{condition}_rep{rep}_GRCh38_peaks.csv.gz',index=False,compression='gzip')

# Merge intervals for efficient, exact Boolean overlap queries.
def interval_index(a):
 result={}
 for c,z in a.groupby('chrom',sort=False):
  arr=z[['start','end']].sort_values('start').to_numpy(int);merged=[]
  for s,e in arr:
   if merged and s<=merged[-1][1]:merged[-1][1]=max(e,merged[-1][1])
   else:merged.append([s,e])
  arr=np.asarray(merged,int);result[c]=(arr[:,0],arr[:,1])
 return result
indices={k:interval_index(v) for k,v in peaks.items()}
def overlaps(index,chrom,s,e):
 if chrom not in index:return False
 a,b=index[chrom];j=np.searchsorted(a,e,side='left')-1;return bool(j>=0 and b[j]>s)
def vector_overlap(index,chrom,start,end):
 if chrom not in index:return np.zeros(len(start),bool)
 a,b=index[chrom];j=np.searchsorted(a,np.asarray(end),side='left')-1;return (j>=0)&(b[np.maximum(j,0)]>np.asarray(start))

h=inherited.h;protein=set(h.loc[h.locus_group=='protein-coding gene','symbol']);tss=pd.read_csv(root/'revision7/analysis/NRF1_GRCh38_TSS_binding.csv.gz');tss=tss[tss.gene.isin(protein)].copy();tss.tss=tss.tss.astype(int)
rna=tidy_matrix(pd.read_csv(root/'revision8/analysis/GSE78167_TPM.csv.gz',index_col=0));cross=pd.read_csv(root/'revision8/analysis/GSE78167_sample_crosswalk.csv');baseline=rna.loc[:,cross.loc[cross.hour==0,'sample']].mean(axis=1);universe=sorted(set(baseline.index[baseline>=1])&set(tss.gene));imp=[g for g in mods['Import_TIM_TOM'] if g in universe];assert len(imp)>=12
tss=tss[tss.gene.isin(universe)].reset_index(drop=True)
for width in [1000,2000]:
 for key,index in indices.items():
  flags=np.zeros(len(tss),bool)
  for chrom,ix in tss.groupby('chrom').groups.items():flags[ix]=vector_overlap(index,chrom,tss.loc[ix,'tss']-width,tss.loc[ix,'tss']+width+1)
  tss[f'{key[0]}_rep{key[1]}_promoter_{width}']=flags
for condition in ['E2','Vehicle']:
 tss[condition+'_promoter']=tss[f'{condition}_rep1_promoter_1000']&tss[f'{condition}_rep2_promoter_1000'];tss[condition+'_promoter_union']=tss[f'{condition}_rep1_promoter_1000']|tss[f'{condition}_rep2_promoter_1000']
gene=tss.groupby('gene').agg(n_TSS=('tss','size'),NRF1_promoter=('NRF1_bound','any'),E2_promoter=('E2_promoter','any'),Vehicle_promoter=('Vehicle_promoter','any'),E2_promoter_union=('E2_promoter_union','any'),Vehicle_promoter_union=('Vehicle_promoter_union','any')).reindex(universe)
gene['baseline_TPM']=baseline.reindex(universe);gene['is_import']=gene.index.isin(imp)
tss.to_csv(out/'ER_NRF1_TSS_annotation.csv.gz',index=False,compression='gzip')

tsindex={c:(z.tss.to_numpy(int),z.gene.to_numpy()) for c,zz in tss.groupby('chrom') for z in [zz.sort_values('tss')]}
gene_ts={g:z.tss.to_numpy(int) for g,z in tss.groupby('gene')}
def linked_regions(path):
 loops=pd.read_csv(path,sep='\t',header=None,usecols=range(7),names=['chrom','s1','e1','chrom2','s2','e2','PET']);N=len(loops)
 loops=loops[(loops.chrom==loops.chrom2)&(((loops.s1+loops.e1-loops.s2-loops.e2).abs()/2)<=1_000_000)].copy();regions=collections.defaultdict(set)
 for chrom,z in loops.groupby('chrom',sort=False):
  if chrom not in tsindex:continue
  pos,gs=tsindex[chrom]
  for a,b,c,d in z[['s1','e1','s2','e2']].itertuples(index=False,name=None):
   for ss,ee,ds,de in [(a,b,c,d),(c,d,a,b)]:
    lo=np.searchsorted(pos,ss-1000);hi=np.searchsorted(pos,ee+1000,side='right')
    for g in set(gs[lo:hi]):
     tt=gene_ts[g]
     if ((tt>=ds-1000)&(tt<=de+1000)).any():continue
     regions[g].add((chrom,int(ds),int(de)))
 diagnostics.append(dict(file=path.name,input_loops=N,intrachromosomal_le1Mb=len(loops),genes_with_links=len(regions)))
 return regions
link1=linked_regions(raw/'ENCFF377SXL.bedpe.gz');link2=linked_regions(raw/'ENCFF597SQA.bedpe.gz')
common=[]
for g in sorted(set(link1)&set(link2)):
 b=sorted(link2[g]);ib=interval_index(pd.DataFrame(b,columns=['chrom','start','end']))
 for chrom,s,e in sorted(link1[g]):
  if chrom not in ib:continue
  starts2,ends2=ib[chrom];lo=np.searchsorted(ends2,s,side='right');hi=np.searchsorted(starts2,e,side='left')
  for j in range(lo,hi):common.append((g,chrom,max(s,starts2[j]),min(e,ends2[j])))
links=pd.DataFrame(common,columns=['gene','chrom','start','end']).drop_duplicates()
for condition in ['E2','Vehicle']:
 both=[];union=[]
 for row in links.itertuples():
  a=overlaps(indices[(condition,1)],row.chrom,row.start,row.end);b=overlaps(indices[(condition,2)],row.chrom,row.start,row.end);both.append(a and b);union.append(a or b)
 links[condition+'_distal']=both;links[condition+'_distal_union']=union
 flags=links.groupby('gene')[condition+'_distal'].any();gene[condition+'_distal']=flags.reindex(gene.index).fillna(False).astype(bool)
 gene[condition+'_promoter_or_distal']=gene[condition+'_promoter']|gene[condition+'_distal']
 gene[condition+'_NRF1_cooccur']=gene[condition+'_promoter']&gene.NRF1_promoter
 gene[condition+'_same_TSS_cooccur']=tss.assign(co=tss[condition+'_promoter']&tss.NRF1_bound).groupby('gene').co.any().reindex(gene.index).fillna(False)
links.to_csv(out/'ER_replicate_supported_linked_distal_regions.csv.gz',index=False,compression='gzip');gene.to_csv(out/'ER_cistrome_gene_annotation.csv.gz',compression='gzip')

# Re-use the identical 10,000 matched panels across all occupancy endpoints.
gene['expression_decile']=pd.qcut(gene.baseline_TPM.rank(method='first'),10,labels=False);gene['TSS_bin']=pd.cut(gene.n_TSS,[0,1,2,5,np.inf],labels=False)
null=np.empty((10000,len(imp)),int);rng=np.random.default_rng(20261025);rowindex={g:i for i,g in enumerate(gene.index)}
for stratum,targets in gene.loc[imp].groupby(['expression_decile','TSS_bin']):
 pool=gene.index[(gene.expression_decile==stratum[0])&(gene.TSS_bin==stratum[1])&~gene.is_import];assert len(pool)>=len(targets)
 positions=[imp.index(g) for g in targets.index];choices=np.array([rng.choice(len(pool),len(targets),replace=False) for _ in range(10000)]);null[:,positions]=np.array([rowindex[g] for g in pool])[choices]
np.savez_compressed(out/'ER_matched_null_gene_panels.npz',panels=null,universe=np.array(gene.index),import_genes=np.array(imp))
tests=[];nullrows=[]
for test,condition,endpoint in [('C1','Vehicle','promoter'),('C2','E2','promoter'),('C3','Vehicle','promoter_or_distal'),('C4','E2','promoter_or_distal'),('C5','Vehicle','NRF1_cooccur'),('C6','E2','NRF1_cooccur')]:
 flags=gene[condition+'_'+endpoint].to_numpy(bool);observed=int(gene.loc[imp,condition+'_'+endpoint].sum());counts=flags[null].sum(axis=1);p=(1+(counts>=observed).sum())/10001
 bg=gene.loc[~gene.is_import,condition+'_'+endpoint];table=[[observed,len(imp)-observed],[int(bg.sum()),len(bg)-int(bg.sum())]];odds,fisher=stats.fisher_exact(table,alternative='greater')
 tests.append(dict(test=test,condition=condition,endpoint=endpoint,n_import=len(imp),n_background=len(gene),observed_bound=observed,matched_expected=counts.mean(),enrichment=observed/counts.mean() if counts.mean()>0 else np.nan,p=p,odds_ratio=odds,Fisher_context_p=fisher,status='tested'))
 nullrows.extend(dict(test=test,iteration=i,count=int(v)) for i,v in enumerate(counts))
bound=[g for g in imp if gene.at[g,'E2_promoter_or_distal']];unbound=[g for g in imp if g not in bound];early=[];onsets=[]
lrna=np.log2(rna+1)
for rep,md in cross.groupby('replicate'):
 base=md.loc[md.hour==0,'sample'].iloc[0]
 for hour in [1,2]:
  current=md.loc[md.hour==hour,'sample'].iloc[0];delta=lrna[current]-lrna[base]
  early.append(dict(replicate=rep,hour=hour,n_bound=len(bound),n_unbound=len(unbound),bound_mean=delta.reindex(bound).mean(),unbound_mean=delta.reindex(unbound).mean(),difference=delta.reindex(bound).mean()-delta.reindex(unbound).mean()))
early=pd.DataFrame(early);early.to_csv(out/'ER_bound_import_early_response.csv',index=False)
if min(len(bound),len(unbound))>=4:
 d=early.loc[early.hour==1,'difference'];se=d.std(ddof=1)/np.sqrt(3);e=d.mean();lo,hi=e+np.array([-1,1])*stats.t.ppf(.975,2)*se
 tests.append(dict(test='C7',condition='E2',endpoint='one-hour bound-minus-unbound induction',n_import=len(imp),n_trajectories=3,n_bound=len(bound),n_unbound=len(unbound),estimate=e,low=lo,high=hi,p=stats.ttest_1samp(d,0).pvalue,status='tested'))
else:tests.append(dict(test='C7',condition='E2',endpoint='one-hour bound-minus-unbound induction',n_bound=len(bound),n_unbound=len(unbound),p=np.nan,status='unavailable',reason='fewer than4 import genes in a binding group'))
tests=pd.DataFrame(tests);q=bh(tests.p.fillna(1));tests['q']=np.where(tests.p.notna(),q,np.nan);tests['family_slots']=7;tests.to_csv(out/'ER_cistrome_primary_tests.csv',index=False);pd.DataFrame(nullrows).to_csv(out/'ER_expression_TSS_matched_null_counts.csv.gz',index=False,compression='gzip');pd.DataFrame(diagnostics).to_csv(out/'ER_cistrome_input_QC.csv',index=False)
gene.loc[imp].to_csv(out/'ER_import_gene_annotation.csv');gene[['expression_decile','TSS_bin','baseline_TPM','n_TSS']].to_csv(out/'ER_matching_covariates.csv.gz',compression='gzip')
print('Background genes',len(gene),'expressed import',len(imp),'E2-linked bound/unbound',len(bound),len(unbound));print(tests.to_string(index=False))
