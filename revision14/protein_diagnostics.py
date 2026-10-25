from v14_common import *
rng=np.random.default_rng(20261025);rows=[];groups=[]
for cohort in ['CCLE','Asleh']:
 d=pd.read_csv(out/(cohort+'_protein_analysis_data.csv'))
 for scale,z in d.groupby('scale'):
  for field in (['ER_status'] if cohort=='CCLE' else ['PAM50 Subtype','ER_status','era']):
   for label,a in z.groupby(field):
    v=a.Import_TIM_TOM.dropna();groups.append(dict(cohort=cohort,scale=scale,group_field=field,group=str(label),n=len(v),mean=v.mean(),median=v.median(),SD=v.std(ddof=1)))
  if cohort=='Asleh':
   for era,a in z.groupby('era'):
    for field in ['luminal','ER_status']:
     counts=a.dropna(subset=['Import_TIM_TOM',field,'Nonimport_mito']).groupby(field).size()
     result=fit_ols(a,'Import_TIM_TOM',field,['Nonimport_mito'],standardize=False) if len(counts)==2 and counts.min()>=10 else dict(status='unavailable',reason='era-specific group below10')
     rows.append(dict(cohort=cohort,scale=scale,era=era,endpoint=field,**result))
pd.DataFrame(groups).to_csv(out/'Independent_proteome_group_descriptives.csv',index=False);pd.DataFrame(rows).to_csv(out/'Asleh_fixed_era_sensitivities.csv',index=False)

z=pd.read_csv(out/'CCLE_protein_dependency_matched_models.csv');results=[]
for y in ['Import_TIM_TOM_GE','import_minus_other_GE']:
 a=z[['Import_TIM_TOM',y]].dropna();v=[]
 for _ in range(5000):
  b=a.iloc[rng.integers(0,len(a),len(a))];v.append(stats.spearmanr(b.iloc[:,0],b.iloc[:,1]).statistic)
 lo,hi=np.nanquantile(v,[.025,.975]);results.append(dict(endpoint=y,n=len(a),rho=stats.spearmanr(a.iloc[:,0],a.iloc[:,1]).statistic,low=lo,high=hi,status='prespecified contextual description, not independent replication'))
pd.DataFrame(results).to_csv(out/'CCLE_protein_dependency_context.csv',index=False)

# Check the constant Asleh primary panel also meets the deposited ER-label groups.
cover=pd.read_csv(out/'Independent_proteome_gene_coverage.csv');a=pd.read_csv(out/'Asleh_fixed_protein_quantities.csv.gz',index_col=0);d=pd.read_csv(out/'Asleh_protein_analysis_data.csv');d=d[d.scale=='source'].set_index('Study_ID');qc=[]
for row in cover[cover.assay=='Asleh_luminal_vs_basal'].itertuples():
 for g in row.genes.split(';'):
  for label,z in d.dropna(subset=['ER_status']).groupby('ER_status'):
   ids=z.index.intersection(a.columns);frac=a.loc[g,ids].notna().mean();qc.append(dict(module=row.module,gene=g,ER_status=label,n=len(ids),fraction_observed=frac,passes=frac>=.8))
qc=pd.DataFrame(qc);qc.to_csv(out/'Asleh_ER_group_panel_coverage_QC.csv',index=False);assert qc.passes.all(),'Asleh ER-specific panel coverage rule not satisfied'
print('Coverage verified in both Asleh ER groups; CCLE binary counts',pd.read_csv(out/'CCLE_protein_analysis_data.csv').query("scale=='source'").ER_status.value_counts().to_dict());print(pd.DataFrame(results).to_string(index=False))
