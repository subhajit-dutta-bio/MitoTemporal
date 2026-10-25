from pathlib import Path
import sys,json,re,hashlib,itertools
r=Path(__file__).parent;old=r.parent/'revision2'
sys.path.insert(0,str(r.parent/'python_packages'))
import pandas as pd,numpy as np
from scipy.stats import wilcoxon

# Full processed-array checks for duplicate public records; no expression-based donor identity inference.
f=open(old/'metadata/GSE147271_soft.txt',encoding='utf-8');current=None;signatures={};meta=[];md={}
for line in f:
 if line.startswith('^SAMPLE = '):
  if current:meta.append(md)
  current=line.strip().split(' = ')[1];md={'sample':current}
 elif current and line.startswith('!Sample_title = '):md['title']=line.strip().split(' = ')[1]
 elif line.startswith('!sample_table_begin'):
  header=next(f);h=hashlib.sha256();n=0
  for ll in f:
   if ll.startswith('!sample_table_end'):break
   h.update(ll.strip().encode());n+=1
  signatures[current]={'sha256_processed_table':h.hexdigest(),'n_rows':n}
if current:meta.append(md)
md=pd.DataFrame(meta);md['patient']=md.title.str.extract(r'Patient(\d+)');md['treatment']=md.title.str.extract(r'_(Tamoxifen|Pretreatment)')
assert len(md)==102 and md[['patient','treatment']].duplicated().sum()==0
both=md.pivot(index='patient',columns='treatment',values='sample')
complete=both.dropna();assert len(complete)==28
pre=set(md.loc[md.treatment=='Pretreatment','patient']);post=set(md.loc[md.treatment=='Tamoxifen','patient'])
hashes=pd.DataFrame(signatures).T;duplicates=hashes[hashes.sha256_processed_table.duplicated(keep=False)]
assert len(duplicates)==0
md.join(hashes,on='sample').to_csv(r/'analysis/tamoxifen_all_array_audit.csv',index=False)
complete.to_csv(r/'analysis/tamoxifen_28_pair_crosswalk.csv')

# Public-label audit: identifiers are local metadata, not a verified global donor crosswalk.
cm=pd.read_csv(old/'analysis/additional_clinical_sample_manifest.csv',dtype=str)
labelsets={acc:set(cm[cm.accession==acc].patient) for acc in ['GSE20181','GSE111563']}
for acc in ['GSE55374','GSE59515']:
 t=(r/'metadata'/f'{acc}_soft.txt').read_text(encoding='utf-8')
 labelsets[acc]=set(re.findall(r'patient(\d+)[A-Z]?',' '.join(re.findall(r'!Sample_title = (.*)',t))))
labelrows=[]
for x,y in itertools.combinations(labelsets,2):
 for label in sorted(labelsets[x]&labelsets[y],key=int):labelrows.append({'series_1':x,'series_2':y,'shared_numeric_label':label,'interpretation':'label collision only; donor match not proven'})
pd.DataFrame(labelrows).to_csv(r/'analysis/letrozole_public_label_collisions.csv',index=False)
tab=pd.read_excel(old/'downloads/Selli_supplement.xlsx',sheet_name='TableS1',header=1)
tab.to_csv(r/'analysis/Selli_TableS1_sample_accession_crosswalk.csv',index=False)
sets={str(acc):set(g['Patient number'].dropna().astype(str)) for acc,g in tab.groupby('Public dataset')}
supp_collision=[]
for x,y in itertools.combinations(sets,2):
 for label in sorted(sets[x]&sets[y]):supp_collision.append({'series_1':x,'series_2':y,'same_supplement_label':label})
pd.DataFrame(supp_collision,columns=['series_1','series_2','same_supplement_label']).to_csv(r/'analysis/Selli_TableS1_shared_patient_labels.csv',index=False)

# Explicit conservative scenarios, fixed from metadata, retaining both complete test families.
scenarios={'exclude_label_50':['50'],'exclude_legacy_labels_32_50':['32','50'],'exclude_all_labels_shared_with_related_series':sorted(labelsets['GSE111563']&(labelsets['GSE20181']|labelsets['GSE55374']|labelsets['GSE59515']),key=int)}
(r/'analysis/overlap_sensitivity_scenarios.json').write_text(json.dumps(scenarios,indent=2))
def bh(p):
 p=np.asarray(p);order=np.argsort(p);out=np.empty(len(p));out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out
base=pd.read_csv(old/'analysis/additional_all_tests.csv');records=[]
manifest=cm[cm.accession=='GSE111563'].set_index('sample')
for name,excluded in scenarios.items():
 aa=base.copy()
 for scale,filename in [('mean_log','module_scores'),('rank','rank_scores')]:
  z=pd.read_csv(old/'analysis'/f'GSE111563_{filename}.csv',index_col=0).join(manifest[['patient','timepoint']]);z=z[~z.patient.isin(excluded)]
  agg=z.groupby(['patient','timepoint']).mean(numeric_only=True)
  for module in agg:
   pv=agg[module].unstack()[['1','4']].dropna();delta=pv['4']-pv['1']
   idx=aa.index[(aa.accession=='GSE111563')&(aa.contrast=='4_minus_1')&(aa.module==module)&(aa.scale==scale)]
   assert len(idx)==1
   aa.loc[idx,'n1']=len(delta);aa.loc[idx,'n2']=len(delta);aa.loc[idx,'mean_delta']=delta.mean();aa.loc[idx,'n_positive']=(delta>0).sum();aa.loc[idx,'p']=wilcoxon(delta).pvalue
 for family,ix in aa.groupby('family').groups.items():aa.loc[ix,'q']=bh(aa.loc[ix,'p'])
 for scale in ['mean_log','rank']:
  ix=aa.index[(aa.scale==scale)&aa.module.ne('CIV_nuclear')];aa.loc[ix,'screen_q']=bh(aa.loc[ix,'p']);other=ix[aa.loc[ix,'module'].ne('OXPHOS_Complex_IV')];aa.loc[other,'q']=aa.loc[other,'screen_q']
 aa['scenario']=name;aa['excluded_GSE111563_labels']=';'.join(excluded);records.extend(aa.to_dict('records'))
allres=pd.DataFrame(records);allres.to_csv(r/'analysis/overlap_sensitivity_all_tests.csv',index=False)
focus=allres[(allres.accession=='GSE111563')&allres.module.isin(['OXPHOS_Complex_IV','Import_TIM_TOM'])]
focus.to_csv(r/'analysis/overlap_sensitivity_key_results.csv',index=False)
summary={'tamoxifen':{'pre':len(pre),'post':len(post),'unique_public_patient_identifiers':len(pre|post),'pairs':len(complete),'duplicate_patient_timepoint_records':0,'exact_duplicate_processed_arrays':len(duplicates),'resolution':'Use 28 deposited pairs, supported by the later primary trial publication PMC9909892. The source 27 is an unreconciled reporting inconsistency, not a basis to discard a patient; no specific extra patient or author correction can be established.'},'letrozole':{'numeric_labels_20181_111563':sorted(labelsets['GSE20181']&labelsets['GSE111563']),'direct_GSM_overlap':[],'donor_overlap_status':'Not identifiable conclusively from public records. Shared programme and related-series reuse are documented; numeric label collisions do not prove donor matches. Do not call the two series independent cohorts or pool them as unique patients.','sensitivity_scenarios':scenarios},'primary_analysis':'Original 28-pair tamoxifen and 37-pair GSE111563 analyses retained. New exclusions are provenance sensitivity scenarios, not post-hoc replacement primary analyses.'}
(r/'analysis/provenance_resolution.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print(focus[['scenario','module','scale','n1','mean_delta','p','q']].to_string(index=False));print('SUPPLEMENT_SHARED_LABELS',supp_collision)
