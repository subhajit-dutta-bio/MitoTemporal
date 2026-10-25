from common import *
from estimate_port import estimate
from scipy.stats import wilcoxon,spearmanr,rankdata
import statsmodels.api as sm
manifest=pd.read_csv(root/'revision2/analysis/additional_clinical_sample_manifest.csv',dtype=str);rows=[];corr=[];adjusted=[];allpairs=[];coverage=[]
for acc,pre,post in [('GSE20181','A','C'),('GSE111563','1','4'),('GSE147271','Pretreatment','Tamoxifen'),('GSE153470','baseline','on')]:
 a=oldmatrix(acc).dropna(axis=0,how='any');rr=a.rank(axis=0,pct=True)
 if acc=='GSE153470':pair=pd.read_csv(root/'revision4/analysis/GSE153470_pair_crosswalk.csv',dtype=str).set_index('patient')
 else:pair=manifest[manifest.accession==acc].pivot(index='patient',columns='timepoint',values='sample')
 pair=pair[[pre,post]].dropna();pair.to_csv(out/f'{acc}_mechanism_pair_crosswalk.csv')
 es=estimate(a);es.to_csv(out/f'{acc}_ESTIMATE_scores.csv');coverage.append(dict(accession=acc,common_genes=es.attrs['common_genes'],n_pairs=len(pair)))
 if acc=='GSE20181':
  es['AffymetrixPurity']=np.cos(.6049872018+.0001467884*es.ESTIMATEScore);es.loc[es.AffymetrixPurity<0,'AffymetrixPurity']=np.nan
 for scale,v in [('source',a),('rank',rr)]:
  s=scores(v).join(es);dd=s.loc[pair[post]].set_axis(pair.index)-s.loc[pair[pre]].set_axis(pair.index)
  for mod in ['Import_TIM_TOM','CIV_nuclear','Mitophagy']:
   for pt in pair.index:allpairs.append(dict(accession=acc,patient=pt,scale=scale,module=mod,pre=s.loc[pair.loc[pt,pre],mod],post=s.loc[pair.loc[pt,post],mod],delta=dd.loc[pt,mod]))
   d=dd[mod];rows.append(dict(accession=acc,module=mod,scale=scale,n=len(d),mean_delta=d.mean(),median_delta=d.median(),n_positive=int((d>0).sum()),p=wilcoxon(d).pvalue))
  if acc not in ['GSE20181','GSE111563']:continue
  # Both platforms: validated signature components, not a universal calibrated purity.
  cv=['Proliferation','StromalSignature','ImmuneSignature'];z=dd[cv];z=z/z.std(ddof=1);X=sm.add_constant(z)
  for mod in ['Import_TIM_TOM','CIV_nuclear']:
   fit=sm.OLS(dd[mod],X).fit(cov_type='HC3');ci=fit.conf_int().loc['const'];vif=[]
   for c in z:
    rf=sm.OLS(z[c],sm.add_constant(z.drop(columns=c))).fit();vif.append(1/(1-rf.rsquared))
   adjusted.append(dict(accession=acc,module=mod,scale=scale,n=len(z),estimate=fit.params['const'],low=ci.iloc[0],high=ci.iloc[1],p=fit.pvalues['const'],model='delta_proliferation_stromal_immune_HC3',condition_number=np.linalg.cond(X),max_VIF=max(vif),zero_covariate_intercept=True,covariate_ranges=';'.join(f'{c}:{dd[c].min():.4g},{dd[c].max():.4g}' for c in cv)))
   for g in ['NRF1','TFAM','PPARGC1A']:
    if g not in dd:
     for kind in ['unadjusted','partial_rank_proliferation_stromal_immune']:corr.append(dict(accession=acc,module=mod,gene=g,scale=scale,n=0,rho=np.nan,p=np.nan,kind=kind,status='not_eligible_complete_gene_coverage'))
     continue
    rho,pv=spearmanr(dd[g],dd[mod]);corr.append(dict(accession=acc,module=mod,gene=g,scale=scale,n=len(dd),rho=rho,p=pv,kind='unadjusted'))
    xs=sm.add_constant(z.apply(rankdata));rx=sm.OLS(rankdata(dd[g]),xs).fit().resid;ry=sm.OLS(rankdata(dd[mod]),xs).fit().resid
    # Partial Spearman t approximation, explicitly separate from unadjusted family.
    rho=np.corrcoef(rx,ry)[0,1];from scipy.stats import t
    pv=2*t.sf(abs(rho)*np.sqrt((len(dd)-5)/(1-rho*rho)),len(dd)-5);corr.append(dict(accession=acc,module=mod,gene=g,scale=scale,n=len(dd),rho=rho,p=pv,kind='partial_rank_proliferation_stromal_immune'))
for name,rs,groups in [('mechanism_paired_tests',rows,['scale']),('mechanism_adjusted_tests',adjusted,['scale']),('mechanism_regulator_correlations',corr,['scale','kind'])]:
 df=pd.DataFrame(rs)
 for _,ix in df.groupby(groups).groups.items():
  df.loc[ix,'q']=bh(df.loc[ix,'p'].fillna(1));df.loc[ix,'family_size']=len(ix)
  df.loc[df.p.isna(),'q']=np.nan
 df.to_csv(out/f'{name}.csv',index=False);print(name);print(df.to_string(index=False))
pd.DataFrame(allpairs).to_csv(out/'mechanism_patient_scores.csv',index=False);pd.DataFrame(coverage).to_csv(out/'ESTIMATE_gene_coverage.csv',index=False)
