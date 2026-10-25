from pathlib import Path
import sys,json,gzip,io,csv,re,hashlib
r=Path(__file__).resolve().parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'revision10'));sys.path.insert(0,str(root/'python_packages'))
import common as inherited
import pandas as pd,numpy as np,statsmodels.api as sm
from scipy import stats
text=gzip.open(r/'metadata/GSE8597_series_matrix.txt.gz','rt').read()
meta={}
for line in text.splitlines():
 if line.startswith(('!Sample_title\t','!Sample_geo_accession\t')):
  row=next(csv.reader([line],delimiter='\t'));meta[row[0]]=row[1:]
md=pd.DataFrame({'sample':meta['!Sample_geo_accession'],'title':meta['!Sample_title']})
md['E2']=md.title.str.contains('_E2_').astype(int);md['CHX']=md.title.str.contains('_CHX_').astype(int);md['E2_CHX']=md.E2*md.CHX;md['replicate_label']=md.title.str.extract(r'rep(\d+)',expand=False).astype(int)
md['group']=np.select([(md.E2==0)&(md.CHX==0),(md.E2==1)&(md.CHX==0),(md.E2==0)&(md.CHX==1)],['Vehicle','E2','CHX'],default='E2+CHX')
assert len(md)==16 and md.groupby(['E2','CHX']).size().eq(4).all() and md['sample'].is_unique
md.to_csv(out/'CHX_sample_crosswalk.csv',index=False)
matrix=pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0]),sep='\t',index_col=0)
assert list(matrix.columns)==md['sample'].tolist() and matrix.shape==(54675,16)
annot=(root/'revision4/metadata/GPL570_annotation.soft.txt').read_text(encoding='utf-8')
annot=pd.read_csv(io.StringIO(annot.split('!platform_table_begin\n')[1].split('!platform_table_end')[0]),sep='\t',dtype=str,keep_default_na=False)
def mapping(v):
 gs={inherited.canon(g.strip()) for g in re.split(r'\s*///\s*',v) if g.strip()}
 return next(iter(gs)) if len(gs)==1 and next(iter(gs)) in inherited.approved else None
annot['symbol']=annot['Gene Symbol'].map(mapping);keys=annot.set_index('ID').symbol
matrix['symbol']=matrix.index.map(keys);valid=matrix.dropna(subset=['symbol']);a=valid.groupby('symbol').median();a=a.dropna(axis=0,how='any')
assert a.shape[1]==16 and a.notna().all().all();a.to_csv(out/'GSE8597_gene_RMA_log2.csv.gz',compression='gzip')
groups=json.loads((root/'revision10/fixed_comparator_gene_sets.json').read_text());groups={'Import_TIM_TOM':inherited.mods['Import_TIM_TOM'],**groups,'ER_programme':inherited.controls['ER_programme']}
coverage=[]
for module,gs in groups.items():
 measured=[g for g in gs if g in a.index];coverage.append(dict(module=module,n_fixed=len(gs),n_measured=len(measured),minimum=max(3,int(np.ceil(len(gs)/2))),genes=';'.join(measured)))
pd.DataFrame(coverage).to_csv(out/'CHX_gene_coverage.csv',index=False)
assert len(set(groups['Import_TIM_TOM'])&set(a.index))>=12
X=sm.add_constant(md[['E2','CHX','E2_CHX']]).astype(float);X.index=md['sample']
contrasts={'E2_without_CHX':np.array([0,1,0,0]),'E2_during_CHX':np.array([0,1,0,1]),'E2_by_CHX':np.array([0,0,0,1])}
rows=[];sample_scores=[];specific=[]
def fitcontrasts(y,module,scale):
 model=sm.OLS(y.loc[X.index].to_numpy(float),X).fit(cov_type='HC3',use_t=True);classic=sm.OLS(y.loc[X.index].to_numpy(float),X).fit();result=[]
 for name,c in contrasts.items():
  effect=float(c@model.params);se=float(np.sqrt(c@model.cov_params()@c));df=float(model.df_resid);p=float(2*stats.t.sf(abs(effect/se),df));lo,hi=effect+np.array([-1,1])*stats.t.ppf(.975,df)*se
  loose=[]
  for i in range(len(y)):
   keep=np.arange(len(y))!=i;sub=sm.OLS(y.loc[X.index].to_numpy(float)[keep],X.iloc[keep]).fit();loose.append(float(c@sub.params))
  se0=float(np.sqrt(c@classic.cov_params()@c));p0=float(2*stats.t.sf(abs(effect/se0),classic.df_resid))
  result.append(dict(scale=scale,module=module,contrast=name,estimate=effect,low=lo,high=hi,p=p,df=df,n=16,n_per_group=4,classic_OLS_p=p0,leave_one_sample_out_min=min(loose),leave_one_sample_out_max=max(loose)))
 return result
for scale,mat in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 s=inherited.scores(mat,groups);joined=md.set_index('sample').join(s);joined['scale']=scale;sample_scores.append(joined.reset_index())
 for module in groups:rows.extend(fitcontrasts(s[module],module,scale))
 for comparator in ['Nonimport_mito','Cytosolic_ribosome','Translation_factors','Proliferation']:
  ss=fitcontrasts(s.Import_TIM_TOM-s[comparator],comparator,scale)
  for rec in ss:rec['comparator']=rec.pop('module')
  specific.extend(ss)
tests=pd.DataFrame(rows)
for _,ix in tests.groupby(['scale','module']).groups.items():tests.loc[ix,'q']=inherited.bh(tests.loc[ix,'p']);tests.loc[ix,'family_size']=3
spec=pd.DataFrame(specific)
for _,ix in spec.groupby('scale').groups.items():spec.loc[ix,'q']=inherited.bh(spec.loc[ix,'p']);spec.loc[ix,'family_size']=12
tests.to_csv(out/'CHX_all_module_tests.csv',index=False);tests[tests.module=='Import_TIM_TOM'].to_csv(out/'CHX_primary_tests.csv',index=False);spec.to_csv(out/'CHX_specificity_tests.csv',index=False);pd.concat(sample_scores).to_csv(out/'CHX_sample_scores.csv',index=False)
gene_rows=[]
for gene in groups['Import_TIM_TOM']+['NRF1','TFAM']:
 for chx in [0,1]:
  value=a.loc[gene,md.loc[(md.CHX==chx)&(md.E2==1),'sample']].mean()-a.loc[gene,md.loc[(md.CHX==chx)&(md.E2==0),'sample']].mean() if gene in a.index else np.nan
  gene_rows.append(dict(gene=gene,CHX=chx,E2_mean_response=value))
pd.DataFrame(gene_rows).to_csv(out/'CHX_gene_descriptive_responses.csv',index=False)
audit={'input_probes':len(matrix),'unambiguously_mapped_probes':len(valid),'complete_unique_approved_genes':len(a),'normalisation':'deposited RMA log2; no re-log','inference':'independent deposited RNA samples; cross-arm pairing and independent culture batches not established','plan_sha256':hashlib.sha256((r/'analysis_plan_CHX.md').read_bytes()).hexdigest(),'platform_sha256':hashlib.sha256((root/'revision4/metadata/GPL570_annotation.soft.txt').read_bytes()).hexdigest()}
(out/'CHX_processing_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(tests[(tests.module=='Import_TIM_TOM')].to_string(index=False));print('coverage',coverage[0]);print('specificity interaction');print(spec[(spec.scale=='source')&(spec.contrast=='E2_by_CHX')].to_string(index=False))
