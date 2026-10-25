from analysis_utils import *
md=metadata('GSE241943');md['experiment']=md.title.str.extract(r' ([PR])\d+$',expand=False);md['group']=md.title.str.replace(r' [PR]\d+$','',regex=True);parts=[];audit=[]
for sample in md.index:
 file=next((r/'metadata').glob(sample+'_*counts.csv.gz'));x=pd.read_csv(file);counts=pd.to_numeric(x['Total gene reads'],errors='raise');assert counts.ge(0).all() and counts.notna().all();counts.index=x.Name.map(canon);lib=counts.sum();parts.append((1e6*counts/lib).groupby(level=0).median().rename(sample));audit.append(dict(sample=sample,library_counts=lib,rows=len(x),source_column='Total gene reads'));md.loc[sample,'animal_source_id']=re.search(r'-(M\d+)-',file.name)[1]
a=np.log2(1+pd.concat(parts,axis=1));a.to_csv(out/'GSE241943_log2_CPM.csv.gz',compression='gzip');freeze(pd.DataFrame(audit),'GSE241943_library_audit.csv');freeze(md.reset_index(),'GSE241943_sample_crosswalk.csv');freeze(pd.DataFrame([dict(module=k,n_fixed=len(gs),n_measured=len(set(gs)&set(a.index))) for k,gs in groups.items()]),'GSE241943_gene_coverage.csv')
allgroups=dict(groups,ER_programme=controls['ER_programme']);primary=[];context=[];specific=[];samples=[]
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v,allgroups).join(md);samples.append(sc.reset_index().assign(scale=scale))
 for experiment,z in sc.groupby('experiment'):
  cdk='Palbociclib 25' if experiment=='P' else 'Ribociclib 25';combo='OP-1250 1 '+cdk;op='OP-1250 1'
  for mod in allgroups:
   rec=welch(z.loc[z.group==op,mod],z.loc[z.group==cdk,mod]);primary.append(dict(scale=scale,experiment=experiment,module=mod,contrast='OP1_minus_CDK25',**rec))
   for treat in [op]+(['OP-1250 10'] if experiment=='P' else [])+[cdk,combo]:context.append(dict(scale=scale,experiment=experiment,module=mod,contrast=treat+'_minus_Vehicle',**welch(z.loc[z.group==treat,mod],z.loc[z.group=='Vehicle',mod])))
   zz=z[z.group.isin(['Vehicle',op,cdk,combo])];o=zz.group.isin([op,combo]).astype(float);c=zz.group.isin([cdk,combo]).astype(float);X=pd.DataFrame({'const':1.,'OP':o,'CDK':c,'interaction':o*c});ff=sm.OLS(zz[mod],X).fit(cov_type='HC3',use_t=True);ci=ff.conf_int().loc['interaction'];context.append(dict(scale=scale,experiment=experiment,module=mod,contrast='OP_by_CDK_interaction',estimate=ff.params.interaction,low=ci.iloc[0],high=ci.iloc[1],p=ff.pvalues.interaction,df=ff.df_resid,n_treated=len(zz),n_control=0))
  for comp in comparators:
   d=z.Import_TIM_TOM-z[comp];specific.append(dict(scale=scale,experiment=experiment,comparator=comp,**welch(d[z.group==op],d[z.group==cdk])))
freeze(family(pd.DataFrame(primary),['scale','module'],2),'CDK_control_primary_tests.csv');freeze(family(pd.DataFrame(context),['scale','module'],9),'CDK_control_context_tests.csv');freeze(family(pd.DataFrame(specific),['scale'],8),'CDK_control_specificity_tests.csv');freeze(pd.concat(samples,ignore_index=True),'CDK_control_sample_scores.csv')
print(pd.read_csv(out/'CDK_control_primary_tests.csv').query('scale=="source" and module in ["Import_TIM_TOM","Proliferation","ER_programme"]')[['experiment','module','estimate','low','high','q']].to_string(index=False))
