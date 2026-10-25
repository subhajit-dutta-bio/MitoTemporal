from pathlib import Path
import sys,json,gzip,csv,re,hashlib
r=Path(__file__).resolve().parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'));sys.path.insert(0,str(root/'revision10'))
import common as inherited
import pandas as pd,numpy as np,statsmodels.api as sm
from scipy import stats
groups={'Import_TIM_TOM':inherited.mods['Import_TIM_TOM'],**json.loads((root/'revision10/fixed_comparator_gene_sets.json').read_text()),'ER_programme':inherited.controls['ER_programme']}
comparators=['Nonimport_mito','Cytosolic_ribosome','Translation_factors','Proliferation']
rescue_groups=['siNT','siNRF1','siNT+NRF1-FL','siNRF1+NRF1-FL','siNT+NRF1-dE7','siNRF1+NRF1-dE7']
platforms=['SaCas9_ABE','SpCas9_ABE','SpCas9_CBE']
rescue_contrasts={
 'KD_vs_control':[-1,1,0,0,0,0],
 'FL_rescue_interaction':[1,-1,-1,1,0,0],
 'DE7_vs_FL_under_KD':[0,0,0,-1,0,1],
 'FL_rescue_vs_KD':[0,-1,0,1,0,0],
 'DE7_rescue_vs_KD':[0,-1,0,0,0,1],
 'FL_rescue_vs_control':[-1,0,0,1,0,0],
 'DE7_rescue_vs_control':[-1,0,0,0,0,1],
 'FL_background_vs_control':[-1,0,1,0,0,0],
 'DE7_background_vs_control':[-1,0,0,0,1,0],
 'DE7_rescue_interaction':[1,-1,0,0,-1,1],
 'isoform_interaction_difference':[0,0,1,-1,-1,1]}
rows=[];spec=[];scores=[];coverage=[];crosswalk=[];audit=[];gene_rows=[];normalisation=[]
raw1=gzip.open(r/'metadata/GSE292988_raw_counts.txt.gz','rb').read()
raw2=gzip.open(r/'metadata/GSE292986_raw_counts.txt.gz','rb').read()
assert raw1==raw2

def fit(y,md,levels,contrasts):
 X=pd.get_dummies(pd.Categorical(md.group,categories=levels),dtype=float)
 y=np.asarray(y,float);fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);classic=sm.OLS(y,X).fit()
 block=pd.get_dummies(md.replicate_label.astype(str),drop_first=True,dtype=float).reset_index(drop=True)
 xb=pd.concat([X,block],axis=1);blocked=sm.OLS(y,xb).fit(cov_type='HC3',use_t=True)
 result=[]
 for name,co in contrasts.items():
  c=np.asarray(co,float);est=float(c@fit.params);se=float(np.sqrt(c@fit.cov_params()@c));df=float(fit.df_resid);p=float(2*stats.t.sf(abs(est/se),df));low,high=est+np.array([-1,1])*stats.t.ppf(.975,df)*se
  loo=[]
  for i in range(len(y)):
   keep=np.arange(len(y))!=i;loo.append(float(c@sm.OLS(y[keep],X.loc[keep]).fit().params))
  se0=float(np.sqrt(c@classic.cov_params()@c));cb=np.r_[c,0,0];eb=float(cb@blocked.params);seb=float(np.sqrt(cb@blocked.cov_params()@cb))
  result.append(dict(contrast=name,estimate=est,low=low,high=high,p=p,df=df,n=len(y),n_per_group=3,classic_OLS_p=float(2*stats.t.sf(abs(est/se0),classic.df_resid)),unverified_block_estimate=eb,unverified_block_p=float(2*stats.t.sf(abs(eb/seb),blocked.df_resid)),leave_one_sample_out_min=min(loo),leave_one_sample_out_max=max(loo)))
 return result

for acc in ['GSE292988','GSE292986']:
 meta={}
 for line in gzip.open(r/'metadata'/f'{acc}_series_matrix.txt.gz','rt'):
  if line.startswith(('!Sample_title\t','!Sample_geo_accession\t','!Sample_description\t')):
   row=next(csv.reader([line],delimiter='\t'));meta[row[0]]=row[1:]
 md=pd.DataFrame({'sample':meta['!Sample_geo_accession'],'title':meta['!Sample_title'],'description':meta['!Sample_description']})
 md['replicate_label']=md.title.str.extract(r'Replicate\s*(\d+)',expand=False).astype(int)
 if acc=='GSE292988':
  md['group']=md.title.str.extract(r'with (.*), Replicate',expand=False).str.replace('DE7','dE7',regex=False)
  md['column']=md.group+'_Replicate'+md.replicate_label.astype(str);levels=rescue_groups;contrasts=rescue_contrasts
 else:
  md['platform']=md.title.map(lambda s:('SaCas9' if 'SaCas9' in s or 'SaCAs9' in s else 'SpCas9')+('_ABE' if 'ABE' in s else '_CBE'))
  md['target']=np.where(md.title.str.contains('intergenic'),'intergenic','NRF1');md['group']=md.platform+'_'+md.target
  md['column']=md.group+'_Replicate'+md.replicate_label.astype(str);levels=[p+'_'+t for p in platforms for t in ['intergenic','NRF1']]
  contrasts={p+'_targeted_vs_intergenic':[(-1 if j==2*i else 1 if j==2*i+1 else 0) for j in range(6)] for i,p in enumerate(platforms)}
 assert len(md)==18 and md['column'].is_unique and set(md.group)==set(levels) and md.groupby('group').size().eq(3).all()
 raw=pd.read_csv(r/'metadata'/f'{acc}_raw_counts.txt.gz',sep='\t',index_col=0)
 assert raw.shape[1]==36 and set(md.column)<=set(raw.columns)
 counts=raw.loc[:,md.column].copy();assert np.isfinite(counts).all().all() and (counts>=0).all().all()
 symbols=[inherited.canon(x.split('_',1)[1]) for x in counts.index]
 mapping=pd.DataFrame({'input_id':counts.index,'symbol':symbols});mapping['approved']=mapping.symbol.isin(inherited.approved);mapping.to_csv(out/f'{acc}_gene_mapping.csv',index=False)
 counts.index=symbols;counts=counts[counts.index.isin(inherited.approved)].groupby(level=0).sum()
 library=counts.sum(axis=0);cpm=counts.div(library,axis=1)*1e6;keep=(cpm>=1).sum(axis=1)>=3;counts=counts.loc[keep]
 positive=counts[(counts>0).all(axis=1)];geo=np.exp(np.log(positive).mean(axis=1));sf=positive.div(geo,axis=0).median(axis=0);sf=sf/np.exp(np.log(sf).mean())
 a=np.log2(counts.div(sf,axis=1)+1);assert a.notna().all().all();a.to_csv(out/f'{acc}_gene_log2_normalized_counts.csv.gz',compression='gzip')
 for col,samp,lib,size in zip(md.column,md['sample'],library,sf):normalisation.append(dict(accession=acc,column=col,sample=samp,approved_mapped_library_abundance=float(lib),median_ratio_size_factor=size))
 for module,gs in groups.items():
  hit=[g for g in gs if g in a.index];minimum=max(3,int(np.ceil(len(gs)/2)));coverage.append(dict(accession=acc,module=module,n_fixed=len(gs),n_measured=len(hit),minimum=minimum,genes=';'.join(hit)))
  assert len(hit)>=minimum,(acc,module,len(hit),minimum)
 assert not set(md['column']) & {x['column'] for x in crosswalk};crosswalk.extend(md.assign(accession=acc).to_dict('records'))
 for scale,mat in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  s=inherited.scores(mat,groups);s=s.loc[md.column].reset_index(drop=True);scores.append(pd.concat([md.reset_index(drop=True),s],axis=1).assign(accession=acc,scale=scale))
  for module in groups:
   for rec in fit(s[module],md,levels,contrasts):
    rec.update(accession=acc,module=module,scale=scale,family='primary' if acc=='GSE292986' or rec['contrast'] in list(contrasts)[:3] else 'secondary');rows.append(rec)
  for comp in comparators:
   for rec in fit(s.Import_TIM_TOM-s[comp],md,levels,{k:v for k,v in contrasts.items() if acc=='GSE292986' or k in list(contrasts)[:3]}):rec.update(accession=acc,comparator=comp,scale=scale);spec.append(rec)
 for g in groups['Import_TIM_TOM']+['NRF1','TFAM']:
  for lev in levels:gene_rows.append(dict(accession=acc,gene=g,group=lev,mean_log2_normalized_count=a.loc[g,md.loc[md.group==lev,'column']].mean() if g in a.index else np.nan,measured=g in a.index))
 audit.append(dict(accession=acc,raw_rows=len(raw),combined_input_libraries=36,selected_libraries=18,approved_expressed_genes=len(a),all_positive_size_factor_genes=len(positive),normalisation='median of ratios, log2(normalized count+1)',decompressed_input_sha256=hashlib.sha256(raw1).hexdigest(),plan_sha256=hashlib.sha256((r/'analysis_plan_NRF1.md').read_bytes()).hexdigest()))
tests=pd.DataFrame(rows)
for _,ix in tests.groupby(['accession','scale','module','family']).groups.items():tests.loc[ix,'q']=inherited.bh(tests.loc[ix,'p']);tests.loc[ix,'family_size']=len(ix)
specific=pd.DataFrame(spec)
for _,ix in specific.groupby('scale').groups.items():specific.loc[ix,'q']=inherited.bh(specific.loc[ix,'p']);specific.loc[ix,'family_size']=24
tests.to_csv(out/'NRF1_all_module_tests.csv',index=False);tests[tests.module=='Import_TIM_TOM'].to_csv(out/'NRF1_import_tests.csv',index=False);specific.to_csv(out/'NRF1_specificity_tests.csv',index=False)
pd.concat(scores).to_csv(out/'NRF1_sample_scores.csv',index=False);pd.DataFrame(crosswalk).to_csv(out/'NRF1_sample_crosswalk.csv',index=False);pd.DataFrame(coverage).to_csv(out/'NRF1_gene_coverage.csv',index=False);pd.DataFrame(normalisation).to_csv(out/'NRF1_normalisation.csv',index=False);pd.DataFrame(gene_rows).to_csv(out/'NRF1_gene_descriptive_responses.csv',index=False)
(out/'NRF1_processing_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(tests[(tests.module=='Import_TIM_TOM')&(tests.family=='primary')][['accession','scale','contrast','estimate','low','high','p','q']].to_string(index=False))
print('secondary source');print(tests[(tests.module=='Import_TIM_TOM')&(tests.family=='secondary')&(tests.scale=='source')][['contrast','estimate','p','q']].to_string(index=False))
