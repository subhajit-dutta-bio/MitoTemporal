from v14_common import *
import io,hashlib
cache=Path(__import__('os').environ.get('MITOTEMPORAL_DEPMAP_DIR',str(r/'raw/depmap')))
inventory={v['name']:v for v in json.loads((r/'metadata/DepMap24Q2_full_inventory.json').read_text())}
inputs={'CRISPRGeneEffect.csv':cache/'CRISPRGeneEffect.csv','Model.csv':cache/'Model.csv','OmicsExpressionProteinCodingGenesTPMLogp1.csv':cache/'OmicsExpressionTPMLogp1.csv'}
audit=[]
for name,p in inputs.items():
 b=p.read_bytes();md5=hashlib.md5(b).hexdigest();ref=inventory[name];assert len(b)==ref['size'] and md5==ref['computed_md5'],(name,md5)
 audit.append(dict(release='DepMap Public24Q2',source_name=name,local_path=str(p),bytes=len(b),md5=md5,sha256=hashlib.sha256(b).hexdigest(),url=ref['download_url'],plan_sha256=sha(r/'analysis_plan_v14.md')))
(r/'metadata/DepMap_cached_input_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
meta=pd.read_csv(inputs['Model.csv']).set_index('ModelID');meta['lineage']=meta.OncotreeLineage
meta['exclusion']='';meta.loc[meta.OncotreePrimaryDisease.str.contains('Non-Cancer|Normal',case=False,na=False),'exclusion']='non-cancerous model'
eng=meta.EngineeredModel.fillna('').astype(str).str.lower();meta.loc[~eng.isin(['','false','no','0','nan']),'exclusion']='explicit engineered model'
eligible=meta[meta.exclusion==''].sort_index().copy();eligible['patient_key']=eligible.PatientID.fillna(pd.Series(eligible.index,index=eligible.index));dup=eligible.patient_key.duplicated(keep='first');meta.loc[eligible.index[dup],'exclusion']='later ModelID sharing PatientID';eligible=eligible[~dup]
labels=pd.read_csv(r/'metadata/Neve2006_receptor_table.csv');labels['key']=labels['Cell line'].map(namekey);labels['ER_status']=labels.ER.map({'+':1,'−':0,'-':0})
labelmap=labels.set_index('key').ER_status
meta['name_key']=meta.StrippedCellLineName.map(namekey);meta['ER_status']=meta.name_key.map(labelmap);meta['ER_source']='Neve2006 explicit non-bracketed receptor annotation';meta.to_csv(out/'DepMap_model_eligibility_and_ER_crosswalk.csv')
def load(p):
 x=pd.read_csv(p,index_col=0);x=x.select_dtypes(include=[np.number]).astype('float32');return tidy_matrix(x.T)
ge=load(inputs['CRISPRGeneEffect.csv']);rna=load(inputs['OmicsExpressionProteinCodingGenesTPMLogp1.csv'])
ge.to_csv(out/'DepMap24Q2_selected_gene_effects.csv.gz',columns=ge.columns.intersection(eligible.index),compression='gzip') if False else None
ids=ge.columns.intersection(eligible.index);md=meta.loc[ids].copy();rows=[];coverage=[]
for name,genes in sets.items():
 hit=[g for g in genes if g in ge.index];coverage.append(dict(assay='CRISPR gene effect',module=name,n_fixed=len(genes),n_measured=len(hit),genes=';'.join(hit)))
 md[name+'_GE']=module_score(ge.loc[:,ids],genes)
for name in ['ER_programme','Proliferation']:
 md[name+'_RNA']=module_score(rna,sets[name]).reindex(ids)
rankrna=rna.rank(axis=0,pct=True);md['ER_programme_RANK']=module_score(rankrna,sets['ER_programme']).reindex(ids)
md['ESR1_RNA']=rna.loc['ESR1'].reindex(ids);md['ESR1_GE']=ge.loc['ESR1'].reindex(ids)
for gene in ['TOMM20','TOMM40','TIMM23']:md[gene+'_GE']=ge.loc[gene].reindex(ids) if gene in ge.index else np.nan
md['global_GE_median']=ge.loc[~ge.index.isin(mods['Import_TIM_TOM']+['ESR1']),ids].median()
md['import_minus_other_GE']=md.Import_TIM_TOM_GE-md.Nonimport_mito_GE
breast=md[md.lineage=='Breast'].copy();breast.to_csv(out/'DepMap_breast_analysis_data.csv');md.to_csv(out/'DepMap_pan_cancer_analysis_data.csv')
pd.DataFrame(coverage).to_csv(out/'DepMap_gene_coverage.csv',index=False)
g=breast.dropna(subset=['ER_status','import_minus_other_GE']);ns=g.groupby('ER_status').size()
rec=fit_ols(g,'import_minus_other_GE','ER_status',standardize=False,min_n=10) if len(ns)==2 and ns.min()>=5 else dict(status='unavailable',reason='ER groups below n=5',n=len(g))
rows.append(dict(test='D1',endpoint='breast ER-positive minus ER-negative differential import dependency',hypothesis='negative',**rec))
rows.append(dict(test='D2',endpoint='breast ER-programme / differential dependency',hypothesis='negative',**fit_ols(breast,'import_minus_other_GE','ER_programme_RNA',['global_GE_median','Proliferation_RNA'],min_n=20)))
for test,y in zip(['D3','D4','D5','D6'],['Import_TIM_TOM_GE','TOMM20_GE','TOMM40_GE','TIMM23_GE']):rows.append(dict(test=test,endpoint=y+' / ESR1 gene effect',hypothesis='positive',**fit_ols(md,y,'ESR1_GE',['global_GE_median'],lineage=True,min_n=100)))
for test,predictor in [('D7','ESR1_RNA'),('D8','ER_programme_RNA')]:rows.append(dict(test=test,endpoint='pan-cancer differential dependency / '+predictor,hypothesis='negative',**fit_ols(md,'import_minus_other_GE',predictor,['global_GE_median','Proliferation_RNA'],lineage=True,min_n=100)))
tests=pd.DataFrame(rows);tests['family']='DepMap_8_primary';tests['p']=tests.get('p',np.nan);tests=frozen_family(tests,['family']);tests.to_csv(out/'DepMap_primary_tests.csv',index=False)
specific=[]
for gene in mods['Import_TIM_TOM']:
 sub=breast[['ER_status','Nonimport_mito_GE']].copy();sub['difference']=ge.loc[gene].reindex(sub.index)-sub.Nonimport_mito_GE if gene in ge.index else np.nan
 z=sub.dropna();nn=z.groupby('ER_status').size();res=fit_ols(z,'difference','ER_status',standardize=False,min_n=10) if len(nn)==2 and nn.min()>=5 else dict(n=len(z),status='unavailable',reason='coverage or n gate')
 specific.append(dict(gene=gene,**res))
sp=pd.DataFrame(specific);sp['family']='24_fixed_import_genes';sp['p']=sp.get('p',np.nan);frozen_family(sp,['family']).to_csv(out/'DepMap_gene_specific_ER_tests.csv',index=False)
sens=[]
sens.append(dict(test='D2_rank_ER',**fit_ols(breast,'import_minus_other_GE','ER_programme_RANK',['global_GE_median','Proliferation_RNA'],min_n=20)))
for test,y in zip(['D3','D4','D5','D6'],['Import_TIM_TOM_GE','TOMM20_GE','TOMM40_GE','TIMM23_GE']):
 sens.append(dict(test=test+'_adjust_other_mito',**fit_ols(md,y,'ESR1_GE',['global_GE_median','Nonimport_mito_GE'],lineage=True,min_n=100)))
 sens.append(dict(test=test+'_exclude_breast',**fit_ols(md[md.lineage!='Breast'],y,'ESR1_GE',['global_GE_median'],lineage=True,min_n=100)))
for test,predictor in [('D7','ESR1_RNA'),('D8','ER_programme_RNA')]:sens.append(dict(test=test+'_exclude_breast',**fit_ols(md[md.lineage!='Breast'],'import_minus_other_GE',predictor,['global_GE_median','Proliferation_RNA'],lineage=True,min_n=100)))
pd.DataFrame(sens).to_csv(out/'DepMap_prespecified_sensitivities.csv',index=False)
loo=[]
for id in breast.index:
 sub=breast.drop(index=id);res=fit_ols(sub,'import_minus_other_GE','ER_programme_RNA',['global_GE_median','Proliferation_RNA'],min_n=20);loo.append(dict(omitted=id,**res))
pd.DataFrame(loo).to_csv(out/'DepMap_breast_leave_one_model_out.csv',index=False)
context=[]
for label,gg in breast.dropna(subset=['ER_status']).groupby('ER_status'):
 for endpoint in ['Import_TIM_TOM_GE','Nonimport_mito_GE','import_minus_other_GE','ESR1_GE','ER_programme_RNA']:
  vals=gg[endpoint].dropna();context.append(dict(ER_status=label,endpoint=endpoint,n=len(vals),mean=vals.mean(),SD=vals.std(ddof=1),median=vals.median()))
pd.DataFrame(context).to_csv(out/'DepMap_breast_group_descriptives.csv',index=False)
des=[]
for y in ['Import_TIM_TOM_GE','TOMM20_GE','TOMM40_GE','TIMM23_GE']:
 zz=md[[y,'ESR1_GE']].dropna();des.append(dict(endpoint=y,predictor='ESR1_GE',n=len(zz),rho=stats.spearmanr(zz[y],zz.ESR1_GE).statistic))
pd.DataFrame(des).to_csv(out/'DepMap_unadjusted_descriptives.csv',index=False)
geneexport=ge.loc[[g for g in mods['Import_TIM_TOM']+comparators['Nonimport_mito']+['ESR1'] if g in ge.index],ids];geneexport.to_csv(out/'DepMap_fixed_gene_effects.csv.gz',compression='gzip')
rna.loc[[g for g in set(mods['Import_TIM_TOM']+comparators['Nonimport_mito']+sets['ER_programme']+sets['Proliferation']+['ESR1']) if g in rna.index]].to_csv(out/'CCLE_fixed_RNA_log2TPM.csv.gz',compression='gzip')
print('Eligible screened models',len(md),'breast',len(breast),'ER groups',ns.to_dict());print(tests[['test','n','estimate','low','high','p','q','status']].to_string(index=False))
