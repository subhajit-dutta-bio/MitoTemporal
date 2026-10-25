from common import *
from scipy.stats import wilcoxon,permutation_test,spearmanr
from statsmodels.duration.hazard_regression import PHReg
def genemap(a,p,field):
 p=p.copy();p['gene']=p[field].map(lambda x:canon(x) if isinstance(x,str) and x and not any(z in x for z in ['///','//',';']) else None);mp=p.set_index('ID').gene;a=a.copy();a['gene']=mp.reindex(a.index);return a.dropna(subset=['gene']).groupby('gene').median()
# NeoPalAna: exact source windows, paired observed genes because deposited two-colour arrays contain missing observations.
a,md,_=matrix('GSE93204');mp=pd.read_csv(out/'GPL6480_probe_gene_map.csv',dtype=str);a=genemap(a,mp,'GENE_SYMBOL');a.to_pickle(out/'GSE93204_full_gene_expression.pkl');pair=md.reset_index().pivot(index='patient id',columns='treatment',values='sample');pair.to_csv(out/'NeoPalAna_pair_crosswalk.csv');res=[];ds=[];cover=[]
for pre,post,window in [('Baseline','C1D1','anastrozole_run_in'),('C1D1','C1D15','added_palbociclib'),('Baseline','Surgery','combined_treatment_and_washout_context')]:
 pp=pair[[pre,post]].dropna()
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  for mod in ['Import_TIM_TOM','CIV_nuclear','Mitophagy']:
   matched=sorted(set(mods[mod])&set(a.index));dd=[]
   for pt,z in pp.iterrows():
    x=v.loc[matched,[z[pre],z[post]]].dropna();ok=len(x)>=max(3,np.ceil(len(mods[mod])/2))
    if scale=='source':cover.append(dict(window=window,patient=pt,module=mod,n_measured=len(x),eligible=bool(ok),genes=';'.join(x.index)))
    if ok:
     d=(x[z[post]]-x[z[pre]]).mean();dd.append(d);ds.append(dict(accession='GSE93204',patient=pt,window=window,module=mod,scale=scale,pre=x[z[pre]].mean(),post=x[z[post]].mean(),delta=d))
   if len(dd)>=3:res.append(dict(accession='GSE93204',window=window,module=mod,scale=scale,n=len(dd),mean_delta=np.mean(dd),n_positive=sum(d>0 for d in dd),p=wilcoxon(dd).pvalue))
z=pd.DataFrame(res)
for scale,ix in z.groupby('scale').groups.items():z.loc[ix,'q']=bh(z.loc[ix,'p']);z.loc[ix,'family_size']=len(ix)
z.to_csv(out/'NeoPalAna_all_tests.csv',index=False);pd.DataFrame(ds).to_csv(out/'NeoPalAna_patient_deltas.csv',index=False);pd.DataFrame(cover).to_csv(out/'NeoPalAna_pair_gene_coverage.csv',index=False);print(z.to_string(index=False))
# GSE75971 is Martin 2016, NOT Nguyen/Magnani 2015 (attachment's PMCID was incorrect).
a,md,_=matrix('GSE75971');mp=pd.read_csv(root/'revision2/analysis/GSE111563_probe_gene_map.csv',dtype=str);a=genemap(a,mp,'gene');md['condition']=md.title.str.extract(r'\.((?:9E2|1wk|\d+wk))\.')[0];md['model']=md.title.str.replace(r'\.(?:9E2|1wk|\d+wk)\.\d+$','',regex=True);md.to_csv(out/'GSE75971_model_metadata.csv');res=[];coverage=[];scoreparts=[]
for model,m in md.groupby('model'):
 base=m.index[m.condition=='9E2'];late=m.index[~m.condition.isin(['9E2','1wk'])]
 s=scores(a.loc[:,list(base)+list(late)].dropna(axis=0,how='any'));s['model']=model;scoreparts.append(s)
 for mod in ['Import_TIM_TOM','CIV_nuclear','Mitophagy']:
  hit=set(mods[mod])&set(a.loc[:,list(base)+list(late)].dropna(axis=0,how='any').index);coverage.append(dict(model=model,module=mod,n_measured=len(hit),eligible=mod in s,genes=';'.join(sorted(hit))))
  if mod not in s:
   res.append(dict(model=model,module=mod,n_parental=len(base),n_late=len(late),p=np.nan,status='insufficient_complete_gene_coverage'));continue
  x=s.loc[base,mod].to_numpy();y=s.loc[late,mod].to_numpy();pv=permutation_test((y,x),lambda y,x:np.mean(y)-np.mean(x),n_resamples=np.inf,alternative='two-sided').pvalue;pooled=np.sqrt(((len(x)-1)*np.var(x,ddof=1)+(len(y)-1)*np.var(y,ddof=1))/(len(x)+len(y)-2))
  res.append(dict(model=model,module=mod,n_parental=len(x),n_late=len(y),mean_parental=x.mean(),mean_late=y.mean(),mean_difference=y.mean()-x.mean(),standardized_difference=(y.mean()-x.mean())/pooled,p=pv))
z=pd.DataFrame(res);z['q']=bh(z.p.fillna(1));z.loc[z.p.isna(),'q']=np.nan;z['family_size']=len(z);z.to_csv(out/'GSE75971_LTED_tests.csv',index=False);pd.DataFrame(coverage).to_csv(out/'GSE75971_gene_coverage.csv',index=False);pd.concat(scoreparts).to_csv(out/'GSE75971_module_scores.csv');print(z.to_string(index=False))
# Metastatic cohort: H is the 97 endocrine-treated patients; C/other treatments excluded.
a,md,meta=matrix('GSE124647');mp=pd.read_csv(root/'revision2/analysis/GSE20181_probe_gene_map.csv',dtype=str);a=genemap(a,mp,'gene').dropna(axis=0,how='any');s=scores(a);md=md.join(s);assert md['protocol number'].is_unique
for c in ['pfs months','pfs.event','oas months','os event','set er_pr rnaseq']:md[c]=pd.to_numeric(md[c],errors='coerce')
md.to_csv(out/'GSE124647_survival_analysis_data.csv');results=[];diagnostics=[]
for endpoint,tc,ec in [('PFS','pfs months','pfs.event'),('OS','oas months','os event')]:
 for mod in ['Import_TIM_TOM','CIV_nuclear']:
  for model,cvs in [('unadjusted',[]),('proliferation_adjusted',['Proliferation']),('SET_RNAseq_complete_case',['Proliferation','set er_pr rnaseq'])]:
   z=md[md['protocol treatment']=='H'][[tc,ec,mod]+cvs].dropna();z=z[z[tc]>0];assert z[ec].isin([0,1]).all();X=z[[mod]+cvs];X=(X-X.mean())/X.std(ddof=1);fit=PHReg(z[tc],X,status=z[ec],ties='efron').fit();ci=fit.conf_int()[0];results.append(dict(endpoint=endpoint,module=mod,model=model,n=len(z),events=int(z[ec].sum()),hazard_ratio=np.exp(fit.params[0]),low=np.exp(ci[0]),high=np.exp(ci[1]),p=fit.pvalues[0],converged_score_maxabs=float(np.max(np.abs(fit.model.score(fit.params)))),SET_coverage='deposited_RNAseq_subset' if cvs and 'set er_pr rnaseq' in cvs else 'not_adjusted_for_SET'))
   sc=fit.schoenfeld_residuals
   for j,c in enumerate(X.columns):
    keep=(z[ec]==1).to_numpy()&np.isfinite(sc[:,j]);rho,pv=spearmanr(np.log(z[tc].to_numpy()[keep]),sc[keep,j]);diagnostics.append(dict(endpoint=endpoint,module=mod,model=model,covariate=c,n_events=int(keep.sum()),rho=rho,p=pv,diagnostic='Schoenfeld_residual_log_time_rank_correlation_not_formal_G_T_test'))
z=pd.DataFrame(results)
for _,ix in z.groupby(['endpoint','model']).groups.items():z.loc[ix,'q']=bh(z.loc[ix,'p']);z.loc[ix,'family_size']=len(ix)
z.to_csv(out/'GSE124647_Cox_tests.csv',index=False);dg=pd.DataFrame(diagnostics);dg['q']=bh(dg.p);dg.to_csv(out/'GSE124647_PH_diagnostics.csv',index=False);print(z.to_string(index=False))
