from common import *
acc='GSE246365';md=metadata(acc);md['model']=md.title.str.extract(r'(HCI-\d+)',expand=False);md['animal_id']=md.title.str.extract(r'#(\d+)',expand=False);md['dose']=md.title.apply(lambda t:0 if 'Vehicle' in t else 80 if '80 mg/kg' in t else 1 if '1 mg/kg' in t else np.nan)
assert md.animal_id.is_unique and md.dose.notna().all() and md.groupby(['model','dose']).size().eq(4).all()
raw=pd.read_csv(r/'metadata/GSE246365_RNAseq_pdx_processed_counts.tsv.gz',sep='\t',index_col=0)
lookup={re.match(r'(SAM\d+)_',z.title)[1]:sample for sample,z in md.iterrows()};raw.columns=[lookup[re.search(r'SAM\d+',c)[0]] for c in raw.columns];assert set(raw.columns)==set(md.index)
raw=raw.apply(pd.to_numeric);assert raw.ge(0).all().all();raw.sum().rename('library_gene_counts').to_csv(out/'PDX_library_sizes.csv')
expr=raw.divide(raw.sum(),axis=1)*1e6;ensg={row.ensembl_gene_id:row.symbol for row in h.itertuples() if row.ensembl_gene_id};expr.index=[ensg.get(str(g).split('.')[0],canon(g)) for g in expr.index];expr=canonicalize(expr);expr=np.log2(1+expr);expr.to_csv(out/'PDX_log2_CPM.csv.gz',compression='gzip')
samples=[];tests=[];coverage=[]
for module,genes in sets.items():coverage.append(dict(module=module,n_fixed=len(genes),n_measured=len(set(genes)&set(expr.index))))
for scale,v in [('source',expr),('rank',expr.rank(axis=0,pct=True))]:
 sc=scores(v);sc['ER_targets_no_ESR1']=scores(v,{'ER_no_ESR1':[g for g in controls['ER_programme'] if g!='ESR1']}).ER_no_ESR1;sc=sc.join(md);samples.append(sc.reset_index().assign(scale=scale))
 for model,z in sc.groupby('model'):
  for module in ['Import_TIM_TOM','CIV_nuclear','ER_targets_no_ESR1','Proliferation','ESR1']:
   for dose in [1,80]:
    x=z.loc[z.dose==dose,module];y=z.loc[z.dose==0,module];v1=x.var(ddof=1)/4;v0=y.var(ddof=1)/4;se=np.sqrt(v1+v0);df=(v1+v0)**2/(v1*v1/3+v0*v0/3);delta=x.mean()-y.mean();ci=stats.t.ppf(.975,df)*se
    tests.append(dict(model=model,module=module,scale=scale,dose_mg_kg=dose,estimate=delta,low=delta-ci,high=delta+ci,p=stats.ttest_ind(x,y,equal_var=False).pvalue,df=df,n_treated=4,n_control=4,scope='post_primary_exploratory'))
res=pd.DataFrame(tests)
for _,ix in res.groupby(['scale','module']).groups.items():assert len(ix)==4;res.loc[ix,'q']=bh(res.loc[ix,'p']);res.loc[ix,'family_size']=4
freeze(res,'PDX_exploratory_tests.csv');freeze(pd.concat(samples,ignore_index=True),'PDX_sample_scores.csv');freeze(md.reset_index(),'PDX_sample_crosswalk.csv');freeze(pd.DataFrame(coverage),'PDX_gene_coverage.csv')
print(res[res.module=='Import_TIM_TOM'].to_string(index=False))
