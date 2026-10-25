from pathlib import Path
r=Path(__file__).parent
source=(r/'analyze_expansion.py').read_text();exec(compile(source.split('# Metadata eligibility is fixed')[0],str(r/'analyze_expansion.py'),'exec'))
from scipy.stats import mannwhitneyu
a,md=matrix('GSE105777');a=array_genes(a,'GPL10558','Symbol','GSE105777')
md['arm']=np.where(md.title.str.startswith('AI.Treated.'),'AI','Control')
md['patient_local']=md.title.str.extract(r'(\d+)[BSD]')[0];md['patient']=md.arm+'_'+md.patient_local
md['timepoint']=np.where(md.title.str.contains(r'\d+[BD](?: |$)'),'baseline','on')
assert md[['patient','timepoint']].duplicated().sum()==0
md.to_csv(out/'POETIC_sample_manifest.csv',index=False)
score('GSE105777',a,'deposited_VST_RSN','clinical_controlled_followup')
xx=md.pivot(index='patient',columns='timepoint',values='sample')[['baseline','on']].dropna();xx['arm']=np.where(xx.index.str.startswith('AI_'),'AI','Control');xx.to_csv(out/'POETIC_pair_crosswalk.csv');assert (xx.arm=='AI').sum()==157 and (xx.arm=='Control').sum()==56
tests=[];within=[];rows=[]
for scale in ['source','rank']:
 ss=score_frames['GSE105777',scale];dd=ss.loc[xx.on].set_axis(xx.index)-ss.loc[xx.baseline].set_axis(xx.index);dd['arm']=xx.arm
 dd.to_csv(out/f'POETIC_{scale}_patient_deltas.csv')
 for module in ss:
  aa=dd.loc[dd.arm=='AI',module];bb=dd.loc[dd.arm=='Control',module]
  means=np.mean(rng.choice(aa.values,(5000,len(aa))),axis=1)-np.mean(rng.choice(bb.values,(5000,len(bb))),axis=1);ci=np.quantile(means,[.025,.975]);p=mannwhitneyu(aa,bb,alternative='two-sided').pvalue
  tests.append({'module':module,'scale':scale,'family':'mito' if module in mods else 'controls','n_AI':len(aa),'n_control':len(bb),'mean_AI':aa.mean(),'mean_control':bb.mean(),'mean_difference':aa.mean()-bb.mean(),'ci_low':ci[0],'ci_high':ci[1],'p':p,'test':'two_sided_Mann_Whitney_on_patient_changes'})
  for arm,v in [('AI',aa),('Control',bb)]:within.append({'module':module,'scale':scale,'arm':arm,'family':'mito' if module in mods else 'controls','n_pairs':len(v),'mean_delta':v.mean(),'n_positive':int((v>0).sum()),'p':wilcoxon(v).pvalue,'test':'paired_Wilcoxon_within_arm'})
tt=pd.DataFrame(tests);ww=pd.DataFrame(within)
for _,ix in tt.groupby(['scale','family']).groups.items():tt.loc[ix,'q']=bh(tt.loc[ix,'p'])
for _,ix in ww.groupby(['scale','family']).groups.items():ww.loc[ix,'q']=bh(ww.loc[ix,'p'])
tt.to_csv(out/'POETIC_between_arm_all_tests.csv',index=False);ww.to_csv(out/'POETIC_within_arm_all_tests.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'POETIC_gene_coverage.csv',index=False)
_,other=matrix('GSE126870');overlap=[]
for patient,z in other.iterrows():
 original=[str(z[c]).split(': ',1)[-1] for c in other if c.startswith('characteristics') and str(z[c]).startswith('accession: ')]
 for gsm in original:
  overlap.append({'GSE126870_sample':patient,'title':z.title,'original_GSM':gsm,'in_GSE105777':gsm in md.index})
pd.DataFrame(overlap).to_csv(out/'POETIC_reanalysis_crosswalk.csv',index=False)
print(tt[tt.module.isin(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV','Proliferation'])].to_string(index=False));print(ww[ww.module=='Import_TIM_TOM'].to_string(index=False));print('INFERENTIAL_BETWEEN',len(tt),'WITHIN',len(ww),'MODULES',ss.columns.tolist());print('GSE126870 baseline-only reanalyses/new samples excluded from paired contrasts; direct accession mappings',len(overlap))
