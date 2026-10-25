from pathlib import Path
import sys,json,re,tarfile,io,gzip,gc,time
import numpy as np,pandas as pd
import scanpy as sc,scanpy.external as sce,anndata as ad
from scipy import io as sio,sparse
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)
ROOT=Path(__import__('os').environ.get('MITOTEMPORAL_REPLICATION_ROOT',str(Path(__file__).resolve().parent)))
ACC=sys.argv[1];maps=json.load(open(ROOT/'scripts/sample_maps.json'))[ACC]
genes=json.load(open(ROOT/'scripts/mito_modules.json'))['OXPHOS_Complex_IV']
out=ROOT/'results'/ACC;out.mkdir(parents=True,exist_ok=True)
sc.settings.n_jobs=4
seqs={'GTCAACTCTTTAGCG':'C0251','TGATGGCCTATTGGG':'C0252','TTCCGCCTCTCTTTG':'C0253','AGTAAGTTCAGCGTA':'C0254'}
def log(x):print(time.strftime('%H:%M:%S'),x,flush=True)
def normbc(s):return re.sub(r'-\d+$','',str(s))
def gztxt(t,n):return gzip.decompress(t.extractfile(n).read()).decode().splitlines()
def readmtx(t,n):
 with gzip.GzipFile(fileobj=t.extractfile(n)) as f:return sio.mmread(f).tocsr().astype(np.float32)
def feature_tag(x):
 for seq,tag in seqs.items():
  if seq in str(x):return tag
 m=re.search(r'C0?(25[1-4])',str(x))
 return 'C0'+m.group(1) if m else None
arrays=[];logs=[]
with tarfile.open(ROOT/'downloads'/(ACC+'_RAW.tar')) as tar:
 names=tar.getnames()
 for lib,mp in maps.items():
  log('Reading '+lib)
  if ACC=='GSE205589':
   name=next(n for n in names if re.search('_'+lib+'_sample_feature_bc_matrix.h5$',n))
   p=out/Path(name).name
   if not p.exists():p.write_bytes(tar.extractfile(name).read())
   a=sc.read_10x_h5(p);a.var_names_make_unique();a.X=a.X.astype(np.float32)
   htname=next(n for n in names if '_'+lib+'_hto_' in n)
   with tarfile.open(fileobj=io.BytesIO(tar.extractfile(htname).read()),mode='r:gz') as ht:
    hn=ht.getnames();features=gztxt(ht,next(n for n in hn if n.endswith('features.tsv.gz')));bcs=gztxt(ht,next(n for n in hn if n.endswith('barcodes.tsv.gz')))
    hmat=readmtx(ht,next(n for n in hn if n.endswith('matrix.mtx.gz'))).T.tocsr()
   valid=[i for i,x in enumerate(features) if feature_tag(x)]
   htcols=[feature_tag(features[i]) for i in valid]
   htdf=pd.DataFrame(hmat[:,valid].toarray(),index=[normbc(x) for x in bcs],columns=htcols)
   assert htdf.index.is_unique
   overlap=np.mean([normbc(x) in htdf.index for x in a.obs_names]);assert overlap>.95,(lib,overlap)
   hto=htdf.reindex([normbc(x) for x in a.obs_names]).fillna(0).to_numpy()
   raw_tags=htcols
  else:
   prefix=next(n[:-len('matrix.mtx.gz')] for n in names if re.search('_'+lib+'_matrix.mtx.gz$',n))
   f=pd.DataFrame([x.split('\t') for x in gztxt(tar,prefix+'features.tsv.gz')]);bcs=gztxt(tar,prefix+'barcodes.tsv.gz')
   mat=readmtx(tar,prefix+'matrix.mtx.gz').T.tocsr()
   gm=f[2].eq('Gene Expression').to_numpy();hm=f[2].eq('Antibody Capture').to_numpy()
   a=ad.AnnData(mat[:,gm].copy(),obs=pd.DataFrame(index=bcs),var=pd.DataFrame(index=f.loc[gm,1].to_numpy()));a.var_names_make_unique()
   raw_tags=[feature_tag(' '.join(row)) for row in f.loc[hm,[0,1]].itertuples(index=False,name=None)];assert all(raw_tags)
   hto=mat[:,hm].toarray();del mat;gc.collect()
  counts=np.asarray(a.X.sum(axis=1)).ravel();ng=np.diff(a.X.indptr)
  mtmask=a.var_names.str.startswith('MT-');mt=np.asarray(a.X[:,mtmask].sum(axis=1)).ravel()/np.maximum(counts,1)*100
  q=(ng>=200)&(mt<40)
  n_initial=a.n_obs;a=a[q].copy();hto=hto[q];mt=mt[q];counts=counts[q];ng=ng[q]
  for i,tag in enumerate(raw_tags):a.obs[tag]=hto[:,i]
  sce.pp.hashsolo(a,raw_tags,priors=(.01,.8,.19),number_of_noise_barcodes=max(1,len(raw_tags)-2))
  a.obs['hashsolo_class']=a.obs['Classification'].astype(str)
  a.obs['hashsolo_p']=a.obs['singlet_hypothesis_probability'].to_numpy()
  sort=np.sort(hto,axis=1);top=np.argmax(hto,axis=1);topcount=sort[:,-1];second=sort[:,-2]
  a.obs['hash_ratio']=topcount/np.maximum(second,1);a.obs['hash_top_count']=topcount
  a.obs['dom_tag']=[raw_tags[i] for i in top]
  a.obs['patient']=[mp.get(x,('Unassigned','Unassigned'))[0] for x in a.obs['hashsolo_class']]
  a.obs['timepoint']=[mp.get(x,('Unassigned','Unassigned'))[1] for x in a.obs['hashsolo_class']]
  a.obs['patient_dom']=[mp.get(x,('Unassigned','Unassigned'))[0] for x in a.obs['dom_tag']]
  a.obs['timepoint_dom']=[mp.get(x,('Unassigned','Unassigned'))[1] for x in a.obs['dom_tag']]
  a.obs['pct_mt']=mt;a.obs['n_genes']=ng;a.obs['total_counts']=counts
  primary=a.obs.patient.ne('Unassigned')&a.obs.hashsolo_p.ge(.9)
  dom=a.obs.patient_dom.ne('Unassigned')&a.obs.hash_ratio.ge(5)&a.obs.hash_top_count.ge(10)
  # Retain source cells for primary and independently specified demultiplexing sensitivity.
  keep=primary|dom
  if ACC=='GSE205589':keep &= a.obs.patient.ne('P4') & a.obs.patient_dom.ne('P4')
  a.obs['hashsolo_pass']=primary
  for marker in ['CD3D','CD3E','TRAC','CD4','CD8A','CD8B']:
   a.obs['UMI_'+marker]=a.X[:,a.var_names.get_loc(marker)].toarray().ravel() if marker in a.var_names else 0
  a.obs['library']=lib;a.obs['accession']=ACC
  a.obs_names=[lib+'|'+x for x in a.obs_names]
  a=a[keep].copy()
  dropcols=list(raw_tags)+['most_likely_hypothesis','cluster_feature','negative_hypothesis_probability','singlet_hypothesis_probability','doublet_hypothesis_probability','Classification']
  a.obs=a.obs.drop(columns=dropcols,errors='ignore')
  logs.append({'library':lib,'initial_cells':n_initial,'qc40_cells':len(q[q]),'retained_hs_or_dominance':a.n_obs,'hashsolo_pass':int(a.obs.hashsolo_pass.sum()),'primary_sample_counts':a.obs.loc[a.obs.hashsolo_pass].groupby(['patient','timepoint'],observed=True).size().to_dict().__str__()})
  log(logs[-1]);arrays.append(a)
log('Concatenating')
a=ad.concat(arrays,join='inner');del arrays;gc.collect()
# Fixed gene pool across all retained cells in this independent cohort.
pool=a.var_names.tolist()
present=[g for g in genes if g in a.var_names];nuclear=[g for g in present if not g.startswith('MT-')]
assert len(present)>=24,(ACC,present)
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
log('Scoring frozen Complex IV module')
sc.tl.score_genes(a,present,score_name='CIV',ctrl_size=50,n_bins=25,random_state=0,gene_pool=pool,use_raw=False)
sc.tl.score_genes(a,nuclear,score_name='CIV_nuclear',ctrl_size=50,n_bins=25,random_state=0,gene_pool=pool,use_raw=False)
a.obs['CIV_mean_log']=np.asarray(a[:,present].X.mean(axis=1)).ravel()
a.obs['CIV_nuclear_mean_log']=np.asarray(a[:,nuclear].X.mean(axis=1)).ravel()
t=(a.obs.UMI_CD3D+a.obs.UMI_CD3E+a.obs.UMI_TRAC)>0
cd4=a.obs.UMI_CD4>0;cd8=(a.obs.UMI_CD8A+a.obs.UMI_CD8B)>0
a.obs['lineage']=np.select([t&cd4&~cd8,t&cd8&~cd4],['CD4','CD8'],default='Unclassified')
a.obs['T_marker_positive']=t;a.obs['CD8B_positive']=t&a.obs.UMI_CD8B.gt(0)&~cd4
a.obs.to_csv(out/'cell_scores.csv.gz',compression='gzip')
pd.DataFrame(logs).to_csv(out/'library_qc.csv',index=False)
report={'accession':ACC,'score_cells':a.n_obs,'gene_pool_size':len(pool),'module_genes':present,'missing_module_genes':[g for g in genes if g not in present],'nuclear_genes':nuclear,'scanpy':sc.__version__,'lineage_counts':a.obs.lineage.value_counts().to_dict(),'methods':'HashSolo posterior>=0.9; conservative RNA lineage gates; patient pairing performed separately after metadata review'}
(out/'processing_summary.json').write_text(json.dumps(report,indent=2));log(report)
log('DONE')
