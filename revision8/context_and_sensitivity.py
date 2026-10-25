from common import *
res=pd.read_csv(out/'GSE89888_primary_tests.csv')
# Conservative small-sample reference distribution, added after primary inspection.
res['p_t18_sensitivity']=2*stats.t.sf(abs((res.estimate/(res.high-res.low)*2*stats.norm.ppf(.975))),18)
for _,ix in res.groupby(['scale','module']).groups.items():res.loc[ix,'q_t18_sensitivity']=bh(res.loc[ix,'p_t18_sensitivity'])
res.to_csv(out/'GSE89888_primary_tests.csv',index=False)
source=pd.read_csv(root/'revision7/deliverables/source_data/Figure_3_f.csv.gz',float_precision='round_trip')
print('Figure3f columns',source.columns.tolist())
rs=[]
for acc,z in source.groupby('accession'):
 for scale in (['source','rank'] if 'scale' in z else ['source']):
  zz=z[z.scale==scale] if 'scale' in z else z
  for c in controls:
   if c not in zz:continue
   rho,pv=stats.spearmanr(zz.Import_TIM_TOM,zz[c],nan_policy='omit');rs.append(dict(accession=acc,scale=scale,control=c,n=len(zz),rho=rho,p=pv))
 # The old panel table has source-scale deltas only. Derive rank deltas independently.
 a=pd.read_csv(root/'revision7/analysis'/f'{acc}_full_gene_expression.csv.gz',index_col=0,float_precision='round_trip').dropna(axis=0,how='any')
 pair=pd.read_csv(root/'revision7/analysis'/f'{acc}_mechanism_pair_crosswalk.csv',dtype=str).set_index('patient');pre,post=('A','C') if acc=='GSE20181' else ('1','4')
 sc=scores(a.rank(axis=0,pct=True));dd=sc.loc[pair[post]].set_axis(pair.index)-sc.loc[pair[pre]].set_axis(pair.index)
 for c in controls:
  rho,pv=stats.spearmanr(dd.Import_TIM_TOM,dd[c]);rs.append(dict(accession=acc,scale='rank',control=c,n=len(dd),rho=rho,p=pv))
if rs:
 a=pd.DataFrame(rs)
 for _,ix in a.groupby('scale').groups.items():a.loc[ix,'q']=bh(a.loc[ix,'p']);a.loc[ix,'family_size']=len(ix)
 freeze(a,'Letrozole_control_correlations.csv')
qc=pd.read_csv(out/'GSE89888_recount_QC.csv')
fields=[c for c in qc if 'uniquely_mapped' in c and '%' in c]
print('QC fields',fields)
print('Import coverage')
for acc in ['GSE89888','GSE78167']:
 z=pd.read_csv(out/f'{acc}_gene_coverage.csv');print(acc,z[z.module.isin(['Import_TIM_TOM','CIV_nuclear','ER_programme','Proliferation'])].to_string(index=False))
print('Rank ESR1',res.query('scale=="rank" and module=="Import_TIM_TOM"')[['cell_line','contrast','estimate','q','q_t18_sensitivity']].to_string(index=False))
