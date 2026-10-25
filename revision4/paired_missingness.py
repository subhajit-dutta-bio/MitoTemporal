import runpy,ast
from pathlib import Path
# Load only imports, definitions and initial constants, without rerunning primary analyses.
r=Path(__file__).parent
source=(r/'analyze_expansion.py').read_text();prefix=source.split('# Metadata eligibility is fixed')[0]
exec(compile(prefix,str(r/'analyze_expansion.py'),'exec'))
a,md=matrix('GSE87411');a=array_genes(a,'GPL6480','GENE_SYMBOL','GSE87411')
md['patient']=md.title.str.extract(r'ID_(\d+)')[0];md['timepoint']=np.where(md.title.str.contains(r'_(?:W2|2wk|M)(?:_| |$)',case=False),'on','baseline')
x=md.pivot(index='patient',columns='timepoint',values='sample')[['baseline','on']].dropna()
assert len(x)==109
rows=[];coverage=[];tests=[]
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 for name,gs in {**mods,**controls}.items():
  measured=[g for g in gs if g in a.index];dd=[]
  for patient,z in x.iterrows():
   pp=v.loc[measured,[z.baseline,z.on]].dropna();eligible=len(pp)>=max(3,np.ceil(len(gs)/2))
   if scale=='source':coverage.append({'patient':patient,'module':name,'n_original':len(gs),'n_paired_observed':len(pp),'eligible':eligible,'paired_observed_genes':';'.join(pp.index)})
   if eligible:
    d=(pp[z.on]-pp[z.baseline]).mean();dd.append(d);rows.append({'accession':'GSE87411','patient':patient,'module':name,'scale':scale,'delta':d,'n_genes':len(pp),'analysis':'paired_observed_gene_sensitivity'})
  if len(dd)>=3:
   ci=boot(dd);tests.append({'accession':'GSE87411','module':name,'scale':scale,'n_pairs':len(dd),'mean_delta':np.mean(dd),'median_delta':np.median(dd),'n_positive':sum(v>0 for v in dd),'p':wilcoxon(dd).pvalue,'ci_low':ci[0],'ci_high':ci[1],'family':'mito' if name in mods else 'controls','analysis':'paired_observed_gene_sensitivity'})
t=pd.DataFrame(tests)
for _,ix in t.groupby(['family','scale']).groups.items():t.loc[ix,'q']=bh(t.loc[ix,'p'])
t.to_csv(out/'GSE87411_paired_observed_tests.csv',index=False)
pd.DataFrame(rows).to_csv(out/'GSE87411_paired_observed_deltas.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'GSE87411_pair_gene_coverage.csv',index=False)
print(t[t.module.isin(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV','Proliferation'])].to_string(index=False))
tt=pd.read_csv(out/'expansion_all_module_tests.csv');print(tt[(tt.accession=='GSE153470')&tt.module.isin(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV','Proliferation'])][['module','scale','n0','mean_delta','p','q']].to_string(index=False))
print(pd.read_csv(out/'clinical_control_associations.csv').to_string(index=False))
