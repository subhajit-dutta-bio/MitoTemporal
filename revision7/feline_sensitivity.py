from common import *
from scipy.stats import wilcoxon
p=pd.read_csv(out/'FELINE_pseudobulk_scores.csv');counts=pd.read_csv(out/'FELINE_patient_time_scores.csv')[['patient','time','platform','n_cells']];p=p.merge(counts,on=['patient','time','platform']);results=[]
for arm in ['A','B','C']:
 for mod in ['Import_TIM_TOM','CIV_nuclear']:
  for post in ['M','E']:
   z=p[(p.platform=='10x')&(p.ARM==arm)&(p.n_cells>=10)].pivot(index='patient',columns='time',values=mod)[['S',post]].dropna();d=z[post]-z.S
   if len(d)<3:continue
   results.append(dict(arm=arm,module=mod,contrast=post+'_minus_S',n=len(d),mean_delta=d.mean(),p=wilcoxon(d).pvalue if (d!=0).any() else 1.,family='A_primary_sensitivity' if arm=='A' else 'BC_exploratory_sensitivity'))
z=pd.DataFrame(results)
for _,ix in z.groupby('family').groups.items():z.loc[ix,'q']=bh(z.loc[ix,'p']);z.loc[ix,'family_size']=len(ix)
z.to_csv(out/'FELINE_pseudobulk_tests.csv',index=False);print(z.to_string(index=False))
manifest=pd.read_csv(root/'revision2/analysis/additional_clinical_sample_manifest.csv',dtype=str);rows=[]
for acc,pre,post in [('GSE20181','A','C'),('GSE147271','Pretreatment','Tamoxifen')]:
 x=manifest[manifest.accession==acc];p=x.pivot(index='patient',columns='timepoint',values='sample');p['analysis_pair']=p[pre].notna()&p[post].notna();p.to_csv(out/f'{acc}_all_identifier_audit.csv');rows.append(dict(accession=acc,records=len(x),unique_patients=x.patient.nunique(),pairs=int(p.analysis_pair.sum()),timepoint_counts=x.timepoint.value_counts().to_dict()))
(out/'patient_count_audit.json').write_text(json.dumps(rows,indent=2));print(rows)
