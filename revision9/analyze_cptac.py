from common import *
rng=np.random.default_rng(20261024)
protein=pd.read_csv(r/'metadata/CPTAC_data_protein_quantification.txt',sep='\t',low_memory=False)
rna=pd.read_csv(r/'metadata/CPTAC_data_mrna_seq_fpkm.txt',sep='\t',low_memory=False)
clin=pd.read_csv(r/'metadata/CPTAC_data_clinical_patient.txt',sep='\t',comment='#');samp=pd.read_csv(r/'metadata/CPTAC_data_clinical_sample.txt',sep='\t',comment='#')
assert clin.PATIENT_ID.is_unique and samp.SAMPLE_ID.is_unique and samp.PATIENT_ID.is_unique
md=samp.merge(clin,on='PATIENT_ID',validate='one_to_one').set_index('SAMPLE_ID')
assert 'Composite.Element.REF' in protein and 'Hugo_Symbol' in rna
protein=protein.rename(columns={'Composite.Element.REF':'Hugo_Symbol'})
protein['Hugo_Symbol']=protein.Hugo_Symbol.str.split('|').str[0]
def expression(df):
 samples=[s for s in md.index if s in df.columns];a=df.set_index('Hugo_Symbol')[samples].apply(pd.to_numeric,errors='coerce');return canonicalize(a)
p=expression(protein);a=expression(rna);patients=md.index.intersection(p.columns).intersection(a.columns)
# Require actual measured modalities, rather than assuming a table column supplies data.
patients=[s for s in patients if p[s].notna().any() and a[s].notna().any()]
md=md.loc[patients];p=p[patients];a=a[patients]
assert len(patients)==len(set(md.PATIENT_ID))
md['ER_positive']=md.ER_UPDATED_CLINICAL_STATUS.eq('Positive');md['ER_known']=md.ER_UPDATED_CLINICAL_STATUS.isin(['Positive','Negative'])
md['Stromal']=pd.to_numeric(md.ESTIMATE_STROMAL_SCORE,errors='coerce');md['Immune']=pd.to_numeric(md.ESTIMATE_IMMUNE_SCORE,errors='coerce')
freeze(md.reset_index(),'CPTAC_sample_crosswalk.csv');a.to_csv(out/'CPTAC_RNA_log2_FPKM.csv.gz',compression='gzip');p.to_csv(out/'CPTAC_protein_log2_ratio.csv.gz',compression='gzip')
mito=sorted(set(g for gs in mods.values() for g in gs)-set(mods['Import_TIM_TOM'])-set(controls['ER_programme'])-set(controls['Proliferation']));mito=[g for g in mito if not g.startswith('MT-')]
covgroups=dict(Import_TIM_TOM=mods['Import_TIM_TOM'],CIV_nuclear=mods['CIV_nuclear'],ER_programme=controls['ER_programme'],Proliferation=controls['Proliferation'],Nonimport_mito=mito)
coverage=[];eligible={}
for kind,expr in [('protein',p),('RNA',a)]:
 for module,genes in covgroups.items():
  good=[g for g in genes if g in expr.index and expr.loc[g].notna().mean()>=.8];eligible[kind,module]=good
  for g in genes:coverage.append(dict(modality=kind,module=module,gene=g,n_measured=int(expr.loc[g].notna().sum()) if g in expr.index else 0,n_tumors=len(md),reliable=g in good))
freeze(pd.DataFrame(coverage),'CPTAC_gene_coverage.csv')
def getscore(expr,genes,fixed):
 threshold=max(3,int(np.ceil(len(fixed)/2)),int(np.ceil(.8*len(genes))))
 if len(genes)<max(3,int(np.ceil(len(fixed)/2))):return pd.Series(np.nan,index=expr.columns)
 z=expr.loc[genes].mean();z[expr.loc[genes].notna().sum()<threshold]=np.nan;return z
def fit(z,xname):
 names=[xname,'Proliferation_RNA','Nonimport_mito_protein','Stromal','Immune'];z=z.dropna(subset=['Import_protein']+names)
 if len(z)<30:return None,None,z
 x=z[names].copy();x=(x-x.mean())/x.std(ddof=1);x=sm.add_constant(x);y=(z.Import_protein-z.Import_protein.mean())/z.Import_protein.std(ddof=1)
 if not np.isfinite(x).all().all() or np.linalg.matrix_rank(x)!=x.shape[1]:return None,None,z
 ff=sm.OLS(y,x).fit(cov_type='HC3',use_t=True);return ff,x,z
scores_all=[];tests=[];missing=[];models=[];extras=[];loo=[];concordance=[]
for scale in ['source','rank','gene_zscore_sensitivity']:
 pp=p if scale=='source' else p.rank(axis=0,pct=True) if scale=='rank' else p.sub(p.mean(axis=1),axis=0).divide(p.std(axis=1,ddof=1),axis=0)
 aa=a if scale=='source' else a.rank(axis=0,pct=True) if scale=='rank' else a.sub(a.mean(axis=1),axis=0).divide(a.std(axis=1,ddof=1),axis=0)
 sc=md[['PATIENT_ID','ER_positive','ER_known','Stromal','Immune']].copy()
 for module,genes in covgroups.items():
  sc[module+'_protein']=getscore(pp,eligible['protein',module],genes);sc[module+'_RNA']=getscore(aa,eligible['RNA',module],genes)
 sc['ESR1_protein']=pp.loc['ESR1'].reindex(sc.index) if 'ESR1' in pp.index else np.nan
 sc['Import_protein']=sc.Import_TIM_TOM_protein
 for module in ['Import_TIM_TOM','CIV_nuclear']:
  shared=[g for g in covgroups[module] if g in eligible['RNA',module] and g in eligible['protein',module]]
  pr=getscore(pp,shared,covgroups[module]);rr=getscore(aa,shared,covgroups[module]);z=pd.DataFrame({'RNA':rr,'protein':pr}).dropna()
  if len(z)>=30:
   rho,pv=stats.spearmanr(z.RNA,z.protein);boot=[]
   for i in range(5000):
    ix=rng.integers(0,len(z),len(z));boot.append(stats.spearmanr(z.RNA.to_numpy()[ix],z.protein.to_numpy()[ix]).statistic)
   ci=np.quantile(boot,[.025,.975]);tests.append(dict(scale=scale,module=module,test='matched_RNA_protein',estimate=rho,low=ci[0],high=ci[1],p=pv,n=len(z),genes=len(shared),method='Spearman_patient_bootstrap'))
   for sample,row in z.iterrows():concordance.append(dict(sample=sample,scale=scale,module=module,RNA=row.RNA,protein=row.protein,ER_positive=md.loc[sample,'ER_positive']))
  else:tests.append(dict(scale=scale,module=module,test='matched_RNA_protein',estimate=np.nan,low=np.nan,high=np.nan,p=np.nan,n=len(z),genes=len(shared),method='unavailable'))
  for xname,test in [('ER_programme_RNA','ER_RNA_adjusted_ERpos'),('ESR1_protein','ESR1_protein_adjusted_ERpos')]:
   zz=sc[sc.ER_positive].copy();zz['Import_protein']=sc[module+'_protein'];ff,x,z=fit(zz,xname)
   if ff is None:tests.append(dict(scale=scale,module=module,test=test,estimate=np.nan,low=np.nan,high=np.nan,p=np.nan,n=len(z),method='unavailable'));continue
   ci=ff.conf_int().loc[xname];tests.append(dict(scale=scale,module=module,test=test,estimate=ff.params[xname],low=ci.iloc[0],high=ci.iloc[1],p=ff.pvalues[xname],n=len(z),genes=len(eligible['protein',module]),df=ff.df_resid,condition_number=np.linalg.cond(x),method='HC3_residual_t'))
   for term in x:
    ci=ff.conf_int().loc[term];models.append(dict(scale=scale,module=module,test=test,term=term,estimate=ff.params[term],low=ci.iloc[0],high=ci.iloc[1],p=ff.pvalues[term],n=len(z)))
   if scale=='source':x.corr().to_csv(out/f'CPTAC_{module}_{test}_predictor_correlations.csv')
   if scale=='source' and module=='Import_TIM_TOM' and xname=='ER_programme_RNA':
    for gene in eligible['protein',module]:
     subset=[g for g in eligible['protein',module] if g!=gene];zz=sc.loc[z.index].copy();zz['Import_protein']=getscore(pp,subset,covgroups[module]).reindex(zz.index);ff2,_,zz=fit(zz,xname)
     if ff2 is not None:loo.append(dict(omitted_gene=gene,n=len(zz),ER_RNA_coefficient=ff2.params[xname]))
    core=[g for g in covgroups[module] if g in p.index and p.loc[g].notna().all()];zz=sc[sc.ER_positive].copy();zz['Import_protein']=getscore(pp,core,covgroups[module]);ff2,_,zz=fit(zz,xname)
    if ff2 is not None:ci=ff2.conf_int().loc[xname];extras.append(dict(scale=scale,test='complete_import_protein_core',estimate=ff2.params[xname],low=ci.iloc[0],high=ci.iloc[1],p=ff2.pvalues[xname],n=len(zz),genes=len(core)))
  # Retain complete comparator result family without a significance-versus-nonsignificance contrast.
 if scale in ['source','rank']:
  pos=sc.loc[sc.ER_positive,'Import_protein'].dropna();neg=sc.loc[sc.ER_known&~sc.ER_positive,'Import_protein'].dropna();v1=pos.var(ddof=1)/len(pos);v0=neg.var(ddof=1)/len(neg);se=np.sqrt(v1+v0);df=(v1+v0)**2/(v1*v1/(len(pos)-1)+v0*v0/(len(neg)-1));delta=pos.mean()-neg.mean();ci=stats.t.ppf(.975,df)*se
  extras.append(dict(scale=scale,test='ERpos_minus_ERneg_context',estimate=delta,low=delta-ci,high=delta+ci,p=stats.ttest_ind(pos,neg,equal_var=False).pvalue,n=len(pos)+len(neg),n_ERpos=len(pos),n_ERneg=len(neg)))
  ff,x,z=fit(sc,'ER_programme_RNA')
  if ff is not None:ci=ff.conf_int().loc['ER_programme_RNA'];extras.append(dict(scale=scale,test='all_tumor_ER_RNA_adjusted_context',estimate=ff.params.ER_programme_RNA,low=ci.iloc[0],high=ci.iloc[1],p=ff.pvalues.ER_programme_RNA,n=len(z)))
 scores_all.append(sc.reset_index().assign(scale=scale))
res=pd.DataFrame(tests)
for _,ix in res.groupby(['scale','module']).groups.items():
 assert len(ix)==3;res.loc[ix,'q']=bh(res.loc[ix,'p'].fillna(1));res.loc[ix,'family_size']=3;res.loc[res.p.isna(),'q']=np.nan
freeze(res,'CPTAC_primary_tests.csv');freeze(pd.concat(scores_all,ignore_index=True),'CPTAC_sample_scores.csv');freeze(pd.DataFrame(models),'CPTAC_model_coefficients.csv');freeze(pd.DataFrame(extras),'CPTAC_context_sensitivity.csv');freeze(pd.DataFrame(loo),'CPTAC_leave_one_import_gene_out.csv');freeze(pd.DataFrame(concordance),'CPTAC_matched_module_concordance.csv')
gene_results=[]
for gene in mods['Import_TIM_TOM']:
 z=pd.DataFrame({'RNA':a.loc[gene] if gene in a.index else np.nan,'protein':p.loc[gene] if gene in p.index else np.nan},index=md.index).dropna()
 rho,pv=stats.spearmanr(z.RNA,z.protein) if len(z)>=20 else (np.nan,np.nan);gene_results.append(dict(gene=gene,n=len(z),rho=rho,p=pv))
g=pd.DataFrame(gene_results);g['q']=bh(g.p.fillna(1));g.loc[g.p.isna(),'q']=np.nan;freeze(g,'CPTAC_import_gene_RNA_protein_correlations.csv')
summary=dict(source='Krug_Cell_2020 via commit-pinned cBioPortal',n_clinical_patients=len(clin),n_protein_columns=p.shape[1],n_RNA_columns=a.shape[1],n_unique_matched_tumors=len(md),n_ER_positive=int(md.ER_positive.sum()),n_ER_negative=int((md.ER_known&~md.ER_positive).sum()),reliable_import_proteins=eligible['protein','Import_TIM_TOM'],reliable_import_RNA=eligible['RNA','Import_TIM_TOM'],n_nonimport_mito_fixed=len(mito),n_nonimport_mito_protein_reliable=len(eligible['protein','Nonimport_mito']),source_commit=json.loads((r/'metadata/CPTAC_commit_record.json').read_text())[0]['sha'])
(out/'CPTAC_population_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(res[res.module=='Import_TIM_TOM'].to_string(index=False))
