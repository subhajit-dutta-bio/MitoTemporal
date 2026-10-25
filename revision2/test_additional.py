from pathlib import Path
import sys,json,re,itertools
r=Path(__file__).parent;out=r/'analysis'
sys.path.insert(0,str(r.parent/'python_packages'))
import pandas as pd,numpy as np
from scipy.stats import wilcoxon,mannwhitneyu,ttest_ind
gene_file=r.parent/'source/hummel2/gene_sets/mito_modules.json'
if not gene_file.exists():gene_file=r/'discovery/mito_modules.json'
mods=json.loads(gene_file.read_text())
mods['CIV_nuclear']=[g for g in mods['OXPHOS_Complex_IV'] if not g.startswith('MT-')]

def bh(p):
 p=np.asarray(p);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def permutation(a,b):
 v=np.r_[a,b];n=len(a);obs=abs(np.mean(b)-np.mean(a));diff=[]
 for ii in itertools.combinations(range(len(v)),n):
  mask=np.zeros(len(v),bool);mask[list(ii)]=True;diff.append(abs(v[~mask].mean()-v[mask].mean()))
 return np.mean(np.array(diff)>=obs-1e-12)

rows=[];pairs=[];patientmeta=[];descriptive=[]
for acc in ['GSE20181','GSE111563','GSE147271','GSE138459','GSE296966']:
 for scale,file in [('mean_log','module_scores'),('rank','rank_scores')]:
  scores=pd.read_csv(out/f'{acc}_{file}.csv',index_col=0)
  if acc in ['GSE20181','GSE111563','GSE147271']:
   md=pd.read_csv(out/f'{acc}_sample_metadata.csv').set_index('sample');parsed=[]
   for sample,x in md.iterrows():
    title=x.title
    if acc=='GSE20181':
     m=re.match(r'^(\d+)([ABC])\b',title);assert m,title
     response='nonresponder' if re.search(r'\bnonresponder\b',title,re.I) else ('responder' if re.search(r'\bresponder\b',title,re.I) else 'unknown')
     patient,time=m.groups()
    elif acc=='GSE111563':
     m=re.search(r'PatientX(\d+)[A-Za-z]*_Timepoint(\d+)',title);assert m,title
     patient,time=m.groups();response='unavailable'
    else:
     m=re.search(r'Patient(\d+)_(Tamoxifen|Pretreatment)',title);assert m,title
     patient,time=m.groups();response='unavailable'
    parsed.append({'sample':sample,'patient':patient,'timepoint':time,'response':response,'title':title})
   md=pd.DataFrame(parsed).set_index('sample')
   if scale=='mean_log':patientmeta.extend(md.reset_index().assign(accession=acc).to_dict('records'))
   a=scores.join(md[['patient','timepoint']]);agg=a.groupby(['patient','timepoint']).mean(numeric_only=True)
   contrasts=[('A','B'),('A','C')] if acc=='GSE20181' else ([('1','4')] if acc=='GSE111563' else [('Pretreatment','Tamoxifen')])
   responsemap={}
   if acc=='GSE20181':
    for pt,g in md.groupby('patient'):
     states=set(g.response)-{'unknown'};responsemap[pt]=next(iter(states)) if len(states)==1 else 'discordant_or_unknown'
   for pre,post in contrasts:
    for module in scores.columns:
     pv=agg[module].unstack('timepoint')[[pre,post]].dropna();delta=pv[post]-pv[pre];n=len(delta)
     if n<3:continue
     p=float(wilcoxon(delta,alternative='two-sided').pvalue) if np.any(delta) else 1.
     rows.append({'accession':acc,'contrast':post+'_minus_'+pre,'module':module,'scale':scale,'test':'paired_Wilcoxon','n1':n,'n2':n,'mean_delta':delta.mean(),'n_positive':int((delta>0).sum()),'p':p,'welch_sensitivity_p':np.nan})
     if scale=='mean_log':
      for patient,v in pv.iterrows():pairs.append({'accession':acc,'contrast':post+'_minus_'+pre,'module':module,'patient':patient,'pre':v[pre],'post':v[post],'delta':v[post]-v[pre],'response':responsemap.get(patient,'unavailable')})
     if acc=='GSE20181':
      groups=pd.Series(delta.index.map(responsemap),index=delta.index)
      responders=delta[groups=='responder'];non=delta[groups=='nonresponder']
      if min(len(responders),len(non))>=3:
       rows.append({'accession':acc,'contrast':'nonresponder_vs_responder_delta_'+post+'_minus_'+pre,'module':module,'scale':scale,'test':'Mann_Whitney_on_patient_changes','n1':len(responders),'n2':len(non),'mean_delta':non.mean()-responders.mean(),'n_positive':np.nan,'p':mannwhitneyu(non,responders,alternative='two-sided').pvalue,'welch_sensitivity_p':np.nan})
  else:
   comparisons=[('polyT_ac','polyT_chr'),('OTI_ac','OTI_chr'),('veh_ac','veh_chr'),('veh_chr','nac_chr')] if acc=='GSE138459' else [('DMSO','Enza')]
   for pre,post in comparisons:
    aa=scores.loc[scores.index.str.startswith(pre+'_')];bb=scores.loc[scores.index.str.startswith(post+'_')]
    assert len(aa)==len(bb)==3
    for module in scores.columns:
     v=aa[module].to_numpy();w=bb[module].to_numpy()
     rows.append({'accession':acc,'contrast':post+'_minus_'+pre,'module':module,'scale':scale,'test':'exact_unpaired_label_permutation','n1':len(v),'n2':len(w),'mean_delta':w.mean()-v.mean(),'n_positive':np.nan,'p':permutation(v,w),'welch_sensitivity_p':ttest_ind(w,v,equal_var=False).pvalue})

res=pd.DataFrame(rows);res['family']='';res['q']=np.nan;res['screen_q']=np.nan;res['welch_sensitivity_q']=np.nan
for ix,x in res.iterrows():
 fam=('CIV_primary' if x.module=='OXPHOS_Complex_IV' else ('CIV_nuclear' if x.module=='CIV_nuclear' else 'other_modules_screen'))+'_'+x.scale
 res.loc[ix,'family']=fam
for fam,idx in res.groupby('family').groups.items():
 res.loc[idx,'q']=bh(res.loc[idx,'p'])
 vv=res.loc[idx,'welch_sensitivity_p'].dropna().index
 if len(vv):res.loc[vv,'welch_sensitivity_q']=bh(res.loc[vv,'welch_sensitivity_p'])
for scale in ['mean_log','rank']:
 idx=res.index[(res.scale==scale)&res.module.ne('CIV_nuclear')]
 res.loc[idx,'screen_q']=bh(res.loc[idx,'p'])
 other=idx[res.loc[idx,'module'].ne('OXPHOS_Complex_IV')]
 res.loc[other,'q']=res.loc[other,'screen_q']
res.to_csv(out/'additional_all_tests.csv',index=False)
primary=res[res.family=='CIV_primary_mean_log'];primary.to_csv(out/'additional_CIV_primary_tests.csv',index=False)
pd.DataFrame(pairs).to_csv(out/'additional_patient_pairs.csv',index=False)
pd.DataFrame(patientmeta).to_csv(out/'additional_clinical_sample_manifest.csv',index=False)

# Descriptive model contrasts, without treating clones or engraftment as independent resistance replications.
raw=pd.read_pickle(out/'clone_matrix_staging.pkl');raw.columns=raw.iloc[5].tolist();raw=raw.iloc[6:]
a=raw[['GENE','Parental','CLONE_1','CLONE_8','SALINE_CH','CLONE_A','CLONE_D','CLONE_F']].copy()
a=a.dropna(subset=['GENE']);a['GENE']=a.GENE.astype(str).str.upper();a=a.set_index('GENE').apply(pd.to_numeric,errors='coerce').groupby(level=0).median().dropna()
a.to_csv(out/'gene_expression/GSE138233_group_gene_expression.csv.gz',compression='gzip')
cov=[];dscores=pd.DataFrame(index=a.columns)
for module,genes in mods.items():
 gs=[g.upper() for g in genes if g.upper() in a.index];ok=len(gs)>=max(3,np.ceil(len(genes)*.5));cov.append({'accession':'GSE138233','module':module,'n_original':len(genes),'n_mapped':len(gs),'eligible':bool(ok),'mapped_genes':';'.join(gs),'missing_genes':';'.join(g for g in genes if g.upper() not in a.index)})
 if not ok:continue
 v=a.loc[gs].mean();dscores[module]=v
 for pre,post in [('Parental','CLONE_1'),('Parental','CLONE_8'),('SALINE_CH','CLONE_A'),('SALINE_CH','CLONE_D'),('SALINE_CH','CLONE_F')]:
  descriptive.append({'accession':'GSE138233','contrast':post+'_minus_'+pre,'module':module,'delta':v[post]-v[pre],'p':'not_tested_clone_group_means'})
dscores.to_csv(out/'GSE138233_module_scores.csv',index_label='group')
cov0=pd.read_csv(out/'gene_coverage.csv');pd.concat([cov0[cov0.accession!='GSE138233'],pd.DataFrame(cov)]).to_csv(out/'gene_coverage.csv',index=False)
md=pd.read_csv(out/'GSE146661_sample_metadata.csv').set_index('sample');a=pd.read_csv(out/'GSE146661_module_scores.csv',index_col=0)
for pt in ['118','124','134','142']:
 sel=md[md.title.str.contains('patient '+pt)];pre=sel[sel.title.str.startswith('Breast tumor')].index[0];post=sel[sel.title.str.startswith('PDX')].index[0]
 for module in a.columns:descriptive.append({'accession':'GSE146661','contrast':'PDX_minus_primary_patient_'+pt,'module':module,'delta':a.loc[post,module]-a.loc[pre,module],'p':'not_tested_site_engraftment_confounded'})
pd.DataFrame(descriptive).to_csv(out/'additional_descriptive_models.csv',index=False)
print('PRIMARY COMPLEX IV\n'+primary.to_string(index=False))
print('SECONDARY SIGNIFICANT MEAN LOG\n'+res[(res.scale=='mean_log')&(res.q<.05)].to_string(index=False))
print('DESCRIPTIVE CIV\n'+pd.DataFrame(descriptive).query('module=="OXPHOS_Complex_IV"').to_string(index=False))
