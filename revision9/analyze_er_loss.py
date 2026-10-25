from common import *
from scipy.stats import ttest_ind
ensg={row.ensembl_gene_id:row.symbol for row in h.itertuples() if row.ensembl_gene_id}
def richmeta(acc):
 rows={};chars={}
 with gzip.open(r/'metadata'/f'{acc}_series_matrix.txt.gz','rt') as f:
  for l in f:
   if not l.startswith('!Sample_'):continue
   row=next(csv.reader([l],delimiter='\t'));key=row[0][8:]
   if key.startswith('characteristics'):
    field=row[1].split(':',1)[0];chars[field]=[v.split(':',1)[-1].strip() for v in row[1:]]
   elif key in ['geo_accession','title','description','treatment_protocol_ch1']:rows.setdefault(key,row[1:])
 z=pd.DataFrame(rows).rename(columns={'geo_accession':'sample'}).set_index('sample')
 for key,values in chars.items():z[key]=values
 return z
def normalize(raw,acc):
 raw=raw.apply(pd.to_numeric,errors='raise');assert raw.notna().all().all() and raw.ge(0).all().all()
 fractional=int(np.sum(np.abs(raw.to_numpy()-np.round(raw.to_numpy()))>1e-8))
 audit.append(dict(accession=acc,n_rows=len(raw),n_samples=raw.shape[1],n_fractional=fractional,minimum=raw.min().min(),maximum=raw.max().max(),normalization='log2(CPM+1)',library_min=raw.sum().min(),library_max=raw.sum().max()))
 raw.sum().rename('library_gene_count').to_csv(out/f'{acc}_library_sizes.csv')
 cp=raw.divide(raw.sum(),axis=1)*1e6
 cp.index=[ensg.get(str(g).split('.')[0],canon(g)) for g in cp.index];cp=canonicalize(cp)
 a=np.log2(1+cp);a.to_csv(out/f'{acc}_log2_CPM.csv.gz',compression='gzip')
 return a
coverage=[];samples=[];primary=[];contexts=[];descriptive=[];audit=[];mapping=[];targets=[]
for acc in ['GSE152323','GSE231397','GSE246363','GSE246364']:
 md=richmeta(acc)
 if acc=='GSE152323':
  raw=pd.read_csv(r/'metadata/GSE152323_rawcount.txt.gz',sep='\t',index_col=0)
  look={}
  for sample,z in md.iterrows():
   group={'CTL':'siCTL','ESR1':'siESR1','RECQ1':'siRecQ1'}[z.title.split(' ')[0]];rep=re.search(r'rep (\d)',z.title)[1];look[group+'_'+rep]=sample
  raw.columns=[look[re.search(r'Sample_([^/]+)',c)[1]] for c in raw.columns]
  md['group']=md.title.str.split().str[0]
 elif acc=='GSE231397':
  x=pd.read_csv(r/'metadata/GSE231397_2022-07-30_Annotated_read_counts.csv.gz');raw=x.filter(regex=r'^counts_').copy();raw.index=x.hgnc_symbol
  # Use the original source sample description, not its position in the file.
  fullcols={}
  with gzip.open(r/'metadata/GSE231397_series_matrix.txt.gz','rt') as f:
   for l in f:
    if l.startswith('!Sample_description'):
     rr=next(csv.reader([l],delimiter='\t'))
     if all(v.startswith('counts_') for v in rr[1:]):fullcols=dict(zip(md.index,rr[1:]));break
  assert len(fullcols)==12
  raw=raw.rename(columns={v:k for k,v in fullcols.items()});assert set(raw.columns)==set(md.index)
  md['group']=md.title.str.extract(r'MCF7-(.*?) biol',expand=False)
 else:
  file='GSE246363_RNAseq_mcf7_crispr_knockouts_processed_count.tsv.gz' if acc=='GSE246363' else 'GSE246364_RNAseq_mcf7_siRNA_processed_counts.tsv.gz'
  raw=pd.read_csv(r/'metadata'/file,sep='\t',index_col=0)
  lookup={re.match(r'(SAM\d+)_',z.title)[1]:sample for sample,z in md.iterrows()}
  raw.columns=[lookup[re.search(r'SAM\d+',c)[0]] for c in raw.columns]
  md['group']=md.treatment
 assert set(raw.columns)==set(md.index) and md.index.is_unique
 md=md.loc[raw.columns];md['accession']=acc;mapping.append(md.reset_index())
 a=normalize(raw,acc)
 for name,gs in sets.items():coverage.append(dict(accession=acc,module=name,n_fixed=len(gs),n_measured=len(set(gs)&set(a.index)),measured_genes=';'.join(g for g in gs if g in a.index)))
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  sc=scores(v);sc['ER_targets_no_ESR1']=scores(v,{'target':[g for g in controls['ER_programme'] if g!='ESR1']}).target
  sc=sc.join(md[['group','title']])
  modules=['Import_TIM_TOM','CIV_nuclear','ER_programme','ER_targets_no_ESR1','Proliferation','ESR1']
  if acc=='GSE152323':
   for module in modules:
    ctl=sc.loc[sc.group=='CTL',module];kd=sc.loc[sc.group=='ESR1',module];assert len(ctl)==len(kd)==3
    v1=kd.var(ddof=1)/3;v0=ctl.var(ddof=1)/3;se=np.sqrt(v1+v0);df=(v1+v0)**2/(v1*v1/2+v0*v0/2);estimate=kd.mean()-ctl.mean();ci=stats.t.ppf(.975,df)*se
    rec=dict(accession=acc,scale=scale,module=module,contrast='ESR1_siRNA_minus_CTL_48h',estimate=estimate,low=estimate-ci,high=estimate+ci,p=ttest_ind(kd,ctl,equal_var=False).pvalue,df=df,n_treated=3,n_control=3,inference='Welch')
    (primary if module in ['Import_TIM_TOM','CIV_nuclear'] else targets).append(rec)
  elif acc=='GSE231397':
   e=sc.group.isin(['E2','Fulv+E2']).astype(float);f=sc.group.isin(['Fulv','Fulv+E2']).astype(float)
   X=pd.DataFrame({'const':1.,'E2':e,'Fulv':f,'interaction':e*f},index=sc.index)
   for module in modules:
    fit=sm.OLS(sc[module],X).fit(cov_type='HC3',use_t=True)
    for name,L in [('Fulv_plus_E2_minus_E2',[0,0,1,1]),('E2_by_Fulv_interaction',[0,0,0,1]),('E2_minus_vehicle_context',[0,1,0,0]),('Fulv_minus_vehicle_context',[0,0,1,0])]:
     L=np.asarray(L);est=float(L@fit.params);se=float(np.sqrt(L@fit.cov_params()@L));ci=stats.t.ppf(.975,fit.df_resid)*se
     rec=dict(accession=acc,scale=scale,module=module,contrast=name,estimate=est,low=est-ci,high=est+ci,p=2*stats.t.sf(abs(est/se),fit.df_resid),df=fit.df_resid,n_treated=3,n_control=3,inference='factorial_HC3_t8')
     if module in ['Import_TIM_TOM','CIV_nuclear'] and not name.endswith('_context'):primary.append(rec)
     else:contexts.append(rec)
  else:
   # The same project contributes these small-n contrasts; do not infer independent derivations.
   if acc=='GSE246363':pairs=[('ER KO','NTC','CRISPR_ESR1_4d'),('NTC + GDC-9545 (300 nM)','NTC','giredestrant_4d')]
   else:
    # Source titles encode treatment and time; their exact mappings are retained.
    for sample,z in md.iterrows():
     text=z.title
     group=('siESR1_72h' if 'siRNA ER, 72' in text else 'siESR1_48h' if 'siRNA ER, 48' in text else 'NTC_72h' if 'NTC, 72' in text else 'NTC_48h' if 'NTC, 48' in text else 'GDC9545' if 'mock, GDC-9545' in text else 'GDC0927' if 'mock, GDC-0927' in text else 'vehicle' if 'mock, DMSO' in text else None)
     assert group,text;sc.loc[sample,'group']=group
    pairs=[('siESR1_48h','NTC_48h','siESR1_48h'),('siESR1_72h','NTC_72h','siESR1_72h'),('GDC9545','vehicle','giredestrant_mock'),('GDC0927','vehicle','GDC0927_mock')]
   for treat,ctl,name in pairs:
    for module in modules:
     x=sc.loc[sc.group==treat,module];y=sc.loc[sc.group==ctl,module];assert len(x)==len(y)==2,(acc,treat,sc.group.tolist())
     descriptive.append(dict(accession=acc,scale=scale,module=module,contrast=name,mean_treated=x.mean(),mean_control=y.mean(),estimate=x.mean()-y.mean(),n_treated=2,n_control=2,inference='descriptive_no_P'))
  samples.append(sc.reset_index().assign(accession=acc,scale=scale))
res=pd.DataFrame(primary)
for _,ix in res.groupby(['scale','module']).groups.items():
 assert len(ix)==3;res.loc[ix,'q']=bh(res.loc[ix,'p']);res.loc[ix,'family_size']=3
freeze(res,'ER_loss_primary_tests.csv');freeze(pd.DataFrame(contexts),'ER_loss_factorial_context.csv');freeze(pd.DataFrame(targets),'ER_loss_KD_target_checks.csv');freeze(pd.DataFrame(descriptive),'ER_loss_small_n_context.csv');freeze(pd.concat(samples,ignore_index=True),'ER_loss_sample_scores.csv');freeze(pd.concat(mapping,ignore_index=True),'ER_loss_sample_crosswalk.csv');freeze(pd.DataFrame(coverage),'ER_loss_gene_coverage.csv');freeze(pd.DataFrame(audit),'ER_loss_input_units_audit.csv')
print(res[res.module=='Import_TIM_TOM'].to_string(index=False))
