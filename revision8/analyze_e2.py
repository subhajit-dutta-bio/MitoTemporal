from common import *
import tarfile
md=metadata('GSE78167');md['replicate']=md.title.str.extract(r'Rep(\d+)',expand=False).astype(int);md['hour']=md.title.str.extract(r'_(\d+)hrs',expand=False).astype(int)
parts=[]
with tarfile.open(r/'metadata/GSE78167_RAW.tar') as tar:
 for item in tar.getmembers():
  if not item.isfile():continue
  sample=item.name.split('_')[0];b=tar.extractfile(item).read();b=gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b
  x=pd.read_csv(io.BytesIO(b),sep='\t');x['symbol']=x.gene_id.str.split('|').str[0];x=x[x.symbol!='?']
  s=x.groupby('symbol').TPM.median();s.name=sample;parts.append(s)
tpm=canonicalize(pd.concat(parts,axis=1));tpm.to_csv(out/'GSE78167_TPM.csv.gz',compression='gzip');a=np.log2(1+tpm)
coverage=[dict(module=k,n_measured=len(set(v)&set(a.index)),n_fixed=len(v),genes=';'.join(g for g in v if g in a.index)) for k,v in sets.items()]
freeze(pd.DataFrame(coverage),'GSE78167_gene_coverage.csv')
samples=[];rows=[];repd=[]
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v).join(md)
 for module in ['Import_TIM_TOM','CIV_nuclear','NRF1','TFAM','ER_programme','Proliferation']:
  for rep,z in sc.groupby('replicate'):
   baseline=z.loc[z.hour==0,module].iloc[0]
   for row in z.itertuples():samples.append(dict(sample=row.Index,replicate=rep,hour=row.hour,scale=scale,module=module,score=getattr(row,module),delta=getattr(row,module)-baseline))
  for window,times in [('late',[6,8,12,24]),('early',[1,2,3])]:
   dd=[]
   for rep,z in sc.groupby('replicate'):
    delta=z.loc[z.hour.isin(times),module].mean()-z.loc[z.hour==0,module].iloc[0];dd.append(delta);repd.append(dict(scale=scale,module=module,window=window,replicate=rep,delta=delta))
   z=tdiff(dd);z.update(scale=scale,module=module,window=window,family=('import_primary' if module=='Import_TIM_TOM' else 'regulators' if module in ['NRF1','TFAM'] else 'comparators') if window=='late' else 'early_secondary' if module in ['Import_TIM_TOM','CIV_nuclear','NRF1','TFAM'] else 'early_context')
   rows.append(z)
results=pd.DataFrame(rows)
for _,ix in results.groupby(['scale','family']).groups.items():results.loc[ix,'q']=bh(results.loc[ix,'p']);results.loc[ix,'family_size']=len(ix)
freeze(pd.DataFrame(samples),'GSE78167_module_trajectories.csv');freeze(pd.DataFrame(repd),'GSE78167_replicate_window_changes.csv');freeze(results,'GSE78167_tests.csv')
freeze(md.reset_index(),'GSE78167_sample_crosswalk.csv')
print(results[(results.scale=='source')].to_string(index=False))
