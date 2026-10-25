from common import *
import warnings
cross=[];tables=[];coverage=[]
for acc in ['GSE105777','GSE87411']:
 frozen=out/f'{acc}_full_observed_expression.csv.gz'
 if frozen.exists():a=pd.read_csv(frozen,index_col=0,float_precision='round_trip')
 elif acc=='GSE87411':
  # The old plotting cache retains only genes complete across all 218 arrays.
  # Reconstruct the source matrix without that global complete-gene filter.
  lines=[];inside=False
  with gzip.open(root/'revision4/metadata/GSE87411_series_matrix.txt.gz','rt') as f:
   for line in f:
    if line.startswith('!series_matrix_table_begin'):inside=True
    elif line.startswith('!series_matrix_table_end'):break
    elif inside:lines.append(line)
  a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0)
  mapping=pd.read_csv(root/'revision4/analysis/GSE87411_probe_to_gene.csv',index_col=0).iloc[:,0]
  a['symbol']=mapping.reindex(a.index);a=a.dropna(subset=['symbol']).groupby('symbol').median();a=canonicalize(a)
  a.to_csv(out/'GSE87411_full_observed_expression.csv.gz',compression='gzip')
 else:
  a=canonicalize(pd.read_pickle(root/'revision4/analysis'/f'{acc}_gene_expression.pkl'))
  a.to_csv(frozen,compression='gzip')
 md=pd.read_csv(root/'revision4/analysis'/f'{acc}_sample_metadata.csv',dtype=str).set_index('sample')
 if acc=='GSE105777':
  clinical=pd.read_excel(r/'metadata/POETIC_2019_MOESM2.xlsx',sheet_name='TableS4',header=2)
  clinical=clinical.rename(columns={clinical.columns[1]:'clinical_key','Baseline.Ki67':'ki67_pre','Surgery.Ki67':'ki67_on','Patient.ID':'patient_id'})
  assert clinical.clinical_key.nunique()==254
  pairs={}
  for sample,row in md.iterrows():
   mt=re.fullmatch(r'AI\.Treated\.(\d+)(B|S)',row.title)
   if mt:key='Treated.'+mt[1];time='baseline' if mt[2]=='B' else 'on';arm='AI'
   else:
    mt=re.fullmatch(r'(\d+)(D|S) \[reanalysis\]',row.title)
    assert mt,row.title
    key='Control.'+mt[1];time='baseline' if mt[2]=='D' else 'on';arm='Untreated'
   assert time not in pairs.setdefault(key,{}),key
   pairs[key][time]=sample;pairs[key]['arm']=arm
  pair=pd.DataFrame.from_dict(pairs,orient='index').rename_axis('clinical_key').reset_index().merge(clinical[['clinical_key','patient_id','ki67_pre','ki67_on','Group','Paired or Baseline.Single']],on='clinical_key',validate='one_to_one')
  assert len(pair)==213 and (pair['Paired or Baseline.Single']=='Paired').all()
  assert ((pair.arm=='AI')==(pair.Group=='Peri AI')).all()
  for cc in ['ki67_pre','ki67_on']:
   pair[cc+'_source']=pair[cc].astype(str);pair[cc]=pd.to_numeric(pair[cc],errors='coerce')
   assert (pair[cc].dropna()>=0).all()
  pair['endpoint_time']='two_week_surgery_IHC'
 else:
  pair=pd.read_csv(root/'revision4/analysis/GSE87411_pair_crosswalk.csv',dtype=str).rename(columns={'patient':'patient_id'})
  for label,col in [('ki67_pre','characteristics_ch2_4'),('ki67_on','characteristics_ch2_5')]:
   pair[label+'_source']=[md.loc[s,col].split(':',1)[-1].strip().replace('%','') for s in pair.baseline]
   pair[label]=pd.to_numeric(pair[label+'_source'],errors='coerce')
   assert all(md.loc[b,col]==md.loc[o,col] for b,o in zip(pair.baseline,pair.on))
  pair['arm']='AI';pair['endpoint_time']='source_EOT_early_2_to_4_week_Z1031B';pair['source_response']=[md.loc[s,'characteristics_ch2_6'].split(': ')[1] for s in pair.baseline]
  valid=pair.ki67_on.notna() & pair.source_response.isin(['sensitive','resistant'])
  assert ((pair.loc[valid,'ki67_on']>10)==(pair.loc[valid,'source_response']=='resistant')).all()
 for col in ['ki67_pre','ki67_on']:pair[col]=pd.to_numeric(pair[col],errors='coerce')
 pair['accession']=acc;cross.append(pair)
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  sc=scores(v)
  for record in pair.to_dict('records'):
   z=dict(record,scale=scale)
   b,o=record['baseline'],record['on']
   for module in ['Import_TIM_TOM','CIV_nuclear','ER_programme','Proliferation','Epithelial']:
    z[module+'_pre']=sc.loc[b,module];z[module+'_on']=sc.loc[o,module]
    hit=[g for g in sets[module] if g in v.index and pd.notna(v.loc[g,b]) and pd.notna(v.loc[g,o])]
    threshold=max(3,int(np.ceil(len(sets[module])/2)))
    z[module+'_delta']=(v.loc[hit,o]-v.loc[hit,b]).mean() if len(hit)>=threshold else np.nan
    if scale=='source':coverage.append(dict(accession=acc,patient_id=record['patient_id'],module=module,pre_genes=sum(g in v.index and pd.notna(v.loc[g,b]) for g in sets[module]),paired_genes=len(hit),fixed_genes=len(sets[module])))
   tables.append(z)
allsc=pd.DataFrame(tables);freeze(pd.concat(cross,ignore_index=True),'Ki67_patient_crosswalk.csv');freeze(pd.DataFrame(coverage),'Ki67_gene_coverage.csv');freeze(allsc,'Ki67_patient_scores.csv')
results=[];pred=[];allcoefs=[];delt=[];missing=[]
def fit_model(z,pseudo=.1,full=True,interaction=False):
 z=z.copy();z['baseline_Ki67']=np.log2(z.ki67_pre+pseudo);y=np.log2(z.ki67_on+pseudo)
 names=['Import_TIM_TOM_pre','baseline_Ki67']+(['ER_programme_pre','Proliferation_pre'] if full else [])
 x=z[names].copy();expr=[g for g in names if g!='baseline_Ki67'];x[expr]=(x[expr]-x[expr].mean())/x[expr].std(ddof=1)
 if interaction:x['AI']=(z.arm=='AI').astype(float);x['Import_x_AI']=x.Import_TIM_TOM_pre*x.AI
 x=sm.add_constant(x);fit=sm.OLS(y,x).fit(cov_type='HC3');return fit,x,y
for scale in ['source','rank']:
 for acc in ['GSE87411','GSE105777','GSE153470']:
  zz=allsc[(allsc.scale==scale)&(allsc.accession==acc)&(allsc.arm=='AI')]
  z=zz.dropna(subset=['ki67_pre','ki67_on','Import_TIM_TOM_pre','ER_programme_pre','Proliferation_pre'])
  if len(z)<30:
   missing.append(dict(accession=acc,scale=scale,n_available=len(z),status='unavailable_patient_linked_IHC' if acc=='GSE153470' else 'insufficient_complete_cases'));continue
  for model,pseudo,full in [('primary_full',.1,True),('baseline_Ki67_only',.1,False),('pseudocount_1_sensitivity',1.,True)]:
   fit,x,y=fit_model(z,pseudo,full);ci=fit.conf_int().loc['Import_TIM_TOM_pre'];rec=dict(accession=acc,scale=scale,model=model,n=len(z),estimate=fit.params.Import_TIM_TOM_pre,low=ci.iloc[0],high=ci.iloc[1],p=fit.pvalues.Import_TIM_TOM_pre,condition_number=np.linalg.cond(x),r_squared=fit.rsquared,df_resid=fit.df_resid)
   rec['p_t_sensitivity']=2*stats.t.sf(abs(fit.tvalues.Import_TIM_TOM_pre),fit.df_resid);results.append(rec)
   for term in x:
    cc=fit.conf_int().loc[term];allcoefs.append(dict(accession=acc,scale=scale,model=model,term=term,estimate=fit.params[term],low=cc.iloc[0],high=cc.iloc[1],p=fit.pvalues[term]))
   if model=='primary_full':
    base=sm.OLS(y,x.drop(columns='Import_TIM_TOM_pre')).fit();aug=sm.OLS(y,x).fit();partial=(base.ssr-aug.ssr)/base.ssr
    for label,ff in [('covariates_only',base),('plus_import',aug)]:
     errors=ff.resid/(1-ff.get_influence().hat_matrix_diag);pred.append(dict(accession=acc,scale=scale,model=label,n=len(z),LOOCV_RMSE=np.sqrt(np.mean(errors**2)),partial_R2_import=partial))
  rho,pv=stats.spearmanr(z.Import_TIM_TOM_delta,np.log2((z.ki67_on+.1)/(z.ki67_pre+.1)),nan_policy='omit');delt.append(dict(accession=acc,scale=scale,n=z.Import_TIM_TOM_delta.notna().sum(),rho=rho,p=pv))
  # Fixed response threshold, secondary endpoint only.
  fit,x,y=fit_model(z);response=(z.ki67_on<=10).astype(int)
  try:
   ff=sm.GLM(response,x,family=sm.families.Binomial()).fit(cov_type='HC3');ci=ff.conf_int().loc['Import_TIM_TOM_pre'];results.append(dict(accession=acc,scale=scale,model='Ki67_le10_logistic_secondary',n=len(z),estimate=ff.params.Import_TIM_TOM_pre,low=ci.iloc[0],high=ci.iloc[1],p=ff.pvalues.Import_TIM_TOM_pre,n_response=int(response.sum())))
  except Exception as exc:missing.append(dict(accession=acc,scale=scale,status='binary_model_failed',error=str(exc)))
 # POETIC treated-versus-control interaction in the mapped subset.
 z=allsc[(allsc.scale==scale)&(allsc.accession=='GSE105777')].dropna(subset=['ki67_pre','ki67_on','Import_TIM_TOM_pre','ER_programme_pre','Proliferation_pre'])
 fit,x,y=fit_model(z,interaction=True);ci=fit.conf_int().loc['Import_x_AI'];results.append(dict(accession='GSE105777',scale=scale,model='POETIC_import_treatment_interaction_secondary',n=len(z),estimate=fit.params.Import_x_AI,low=ci.iloc[0],high=ci.iloc[1],p=fit.pvalues.Import_x_AI))
res=pd.DataFrame(results)
for (scale,model),ix in res.groupby(['scale','model']).groups.items():
 denominator=1 if model.startswith('POETIC_') else 3;p=res.loc[ix,'p'].tolist()+[1.]*(denominator-len(ix));res.loc[ix,'q']=bh(p)[:len(ix)];res.loc[ix,'family_size']=denominator
 if 'p_t_sensitivity' in res and res.loc[ix,'p_t_sensitivity'].notna().all():res.loc[ix,'q_t_sensitivity']=bh(res.loc[ix,'p_t_sensitivity'].tolist()+[1.]*(denominator-len(ix)))[:len(ix)]
dd=pd.DataFrame(delt)
for _,ix in dd.groupby('scale').groups.items():dd.loc[ix,'q']=bh(dd.loc[ix,'p'].tolist()+[1.])[:len(ix)]
freeze(res,'Ki67_model_tests.csv');freeze(pd.DataFrame(pred),'Ki67_LOOCV_prediction.csv');freeze(pd.DataFrame(allcoefs),'Ki67_all_model_coefficients.csv');freeze(dd,'Ki67_delta_correlations.csv');freeze(pd.DataFrame(missing),'Ki67_unavailable_tests.csv')
print(res[(res.scale=='source')].to_string(index=False));print(pd.DataFrame(pred).query('scale == "source"').to_string(index=False));print(dd.to_string(index=False))
