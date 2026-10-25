from common import *
raw=pd.read_csv(r/'metadata/SRP093386_gene_sums.G026.gz',sep='\t',comment='#',index_col=0)
md=pd.read_csv(r/'metadata/SRP093386_sra.MD.gz',sep='\t').set_index('external_id')
geo=metadata('GSE89888');assert set(md.sample_name)==set(geo.index) and len(md)==48
assert (md.sample_title.to_numpy()==geo.loc[md.sample_name,'title'].to_numpy()).all()
md[['cell_line','treatment','genotype','replicate']]=md.sample_title.str.split('-',expand=True).values
assert (md.groupby(['cell_line','genotype','treatment']).size()==4).all()
# Sum gene coverage prior to dropping noncoding or unmapped genes; do not relabel as TPM.
cp=raw.divide(raw.sum(axis=0),axis=1)*1e6
names={}
with gzip.open(r/'metadata/recount_G026.gtf.gz','rt') as f:
 for l in f:
  if l.startswith('#'):continue
  fields=l.rstrip().split('\t');attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[-1]));names[attrs['gene_id']]=canon(attrs['gene_name'])
cp.index=[names.get(g,g) for g in cp.index];cp=cp.groupby(level=0).median();a=np.log2(1+cp);a.to_csv(out/'GSE89888_log2_coverage_CPM.csv.gz',compression='gzip')
freeze(md.reset_index(),'GSE89888_sample_crosswalk.csv')
freeze(pd.DataFrame([dict(module=k,n_measured=len(set(v)&set(a.index)),n_fixed=len(v),genes=';'.join(g for g in v if g in a.index)) for k,v in sets.items()]),'GSE89888_gene_coverage.csv')
allsc=[];results=[];context=[]
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v).join(md[['cell_line','treatment','genotype','replicate','sample_name']]);allsc.append(sc.reset_index().assign(scale=scale))
 for line,z in sc.groupby('cell_line'):
  X=pd.DataFrame({'const':1.,'Y537S':(z.genotype=='Y537S').astype(float),'D538G':(z.genotype=='D538G').astype(float),'E2':(z.treatment=='E2').astype(float)},index=z.index)
  X['Y537S_E2']=X.Y537S*X.E2;X['D538G_E2']=X.D538G*X.E2
  for module in ['Import_TIM_TOM','CIV_nuclear','ER_programme','Proliferation']:
   fit=sm.OLS(z[module],X).fit(cov_type='HC3')
   for name in ['E2','Y537S','D538G','Y537S_E2','D538G_E2']:
    ci=fit.conf_int().loc[name];rec=dict(scale=scale,cell_line=line,module=module,contrast=name,estimate=fit.params[name],low=ci.iloc[0],high=ci.iloc[1],p=fit.pvalues[name],n_cultures=len(z),n_per_cell=4)
    (results if module in ['Import_TIM_TOM','CIV_nuclear'] else context).append(rec)
result=pd.DataFrame(results)
for _,ix in result.groupby(['scale','module']).groups.items():result.loc[ix,'q']=bh(result.loc[ix,'p']);result.loc[ix,'family_size']=len(ix)
freeze(result,'GSE89888_primary_tests.csv');freeze(pd.DataFrame(context),'GSE89888_control_tests_exploratory.csv');freeze(pd.concat(allsc,ignore_index=True),'GSE89888_sample_scores.csv')
qc=pd.read_csv(r/'metadata/SRP093386_recount_qc.MD.gz',sep='\t');qc.to_csv(out/'GSE89888_recount_QC.csv',index=False)
print(result[(result.scale=='source')&(result.module=='Import_TIM_TOM')].to_string(index=False))
