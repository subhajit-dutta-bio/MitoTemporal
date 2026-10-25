from analysis_utils import *
tests=[];module=[];samples=[];coverage=[]
def record(context,scale,sc,treat=None,ctl=None,paired=False):
 for name in groups:
  result=onesample(sc[name]) if paired else welch(sc.loc[treat,name],sc.loc[ctl,name]);module.append(dict(context=context,scale=scale,module=name,**result))
 for name in comparators:
  d=sc.Import_TIM_TOM-sc[name];result=onesample(d) if paired else welch(d.loc[treat],d.loc[ctl]);tests.append(dict(context=context,scale=scale,comparator=name,**result))
 samples.append(sc.reset_index().assign(context=context,scale=scale))
def cover(acc,a):
 for name,genes in groups.items():coverage.append(dict(accession=acc,module=name,n_fixed=len(genes),n_measured=len(set(genes)&set(a.index))))
for acc in ['GSE20181','GSE111563','GSE105777']:
 a,p=clinical(acc);cover(acc,a)
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  sc=paired_scores(v,p).set_index('unit');record(acc,scale,sc,sc.arm=='AI',sc.arm=='Untreated',paired=acc!='GSE105777')
a=np.log2(1+pd.read_csv(root/'revision8/analysis/GSE78167_TPM.csv.gz',index_col=0));md=pd.read_csv(root/'revision8/analysis/GSE78167_sample_crosswalk.csv').set_index('sample');cover('GSE78167',a)
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v,groups).join(md);rows=[]
 for rep,z in sc.groupby('replicate'):
  late=z[z.hour.isin([6,8,12,24])][list(groups)].mean();base=z.loc[z.hour==0,list(groups)].iloc[0];rows.append((late-base).rename(str(rep)))
 record('GSE78167_late',scale,pd.DataFrame(rows),paired=True)
a=pd.read_csv(root/'revision8/analysis/GSE89888_log2_coverage_CPM.csv.gz',index_col=0);md=pd.read_csv(root/'revision8/analysis/GSE89888_sample_crosswalk.csv').set_index('external_id');cover('GSE89888',a)
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v,groups).join(md[['cell_line','genotype','treatment']]);sc=sc[sc.genotype=='WT']
 for line,z in sc.groupby('cell_line'):record('GSE89888_'+line,scale,z,z.treatment=='E2',z.treatment=='veh')
for acc in ['GSE152323','GSE231397']:
 a=pd.read_csv(root/'revision9/analysis'/f'{acc}_log2_CPM.csv.gz',index_col=0);md=pd.read_csv(root/'revision9/analysis/ER_loss_sample_scores.csv').query('accession==@acc and scale=="source"').set_index('index');cover(acc,a)
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  sc=scores(v,groups).join(md[['group']]);t,c=('ESR1','CTL') if acc=='GSE152323' else ('Fulv+E2','E2');sc=sc[sc.group.isin([t,c])];record(acc,scale,sc,sc.group==t,sc.group==c)
a=pd.read_csv(root/'revision9/analysis/PDX_log2_CPM.csv.gz',index_col=0);md=pd.read_csv(root/'revision9/analysis/PDX_sample_scores.csv').query('scale=="source"').set_index('index');cover('GSE246365',a)
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v,groups).join(md[['model','dose','animal_id']])
 for model,z in sc.groupby('model'):
  for dose in [1,80]:record(f'GSE246365_{model}_{dose}',scale,z,z.dose==dose,z.dose==0)
res=family(pd.DataFrame(tests),['scale'],48);freeze(res,'specificity_primary_tests.csv');freeze(pd.DataFrame(module),'specificity_module_effects_context.csv');freeze(pd.concat(samples,ignore_index=True),'specificity_unit_scores.csv');freeze(pd.DataFrame(coverage),'specificity_gene_coverage.csv')
print(res.groupby(['scale','comparator']).agg(tests=('p','size'),corrected=('q',lambda v:(v<.05).sum())).to_string())
print(res[(res.scale=='source')&(res.comparator=='Nonimport_mito')][['context','estimate','low','high','q']].to_string(index=False))
