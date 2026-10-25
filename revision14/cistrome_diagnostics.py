from v14_common import *
gene=pd.read_csv(out/'ER_cistrome_gene_annotation.csv.gz',index_col=0);tss=pd.read_csv(out/'ER_NRF1_TSS_annotation.csv.gz');links=pd.read_csv(out/'ER_replicate_supported_linked_distal_regions.csv.gz');imp=[g for g in mods['Import_TIM_TOM'] if g in gene.index];rows=[]
for condition in ['E2','Vehicle']:
 for width in [1000,2000]:
  a=tss[f'{condition}_rep1_promoter_{width}'];b=tss[f'{condition}_rep2_promoter_{width}']
  for support,flags in [('both_replicates',a&b),('either_replicate',a|b)]:
   by=tss.assign(flag=flags).groupby('gene').flag.any();rows.append(dict(condition=condition,width=width,support=support,n_import=len(imp),import_promoter_bound=int(by.reindex(imp).sum()),background_bound=int(by.sum()),n_background=len(by)))
pd.DataFrame(rows).to_csv(out/'ER_cistrome_fixed_sensitivity_counts.csv',index=False)
positive=['GREB1','TFF1','PGR','PDZK1','CA12','XBP1','NRIP1'];gene.reindex(positive).to_csv(out/'ER_cistrome_positive_control_genes.csv')
early=pd.read_csv(out/'ER_bound_import_early_response.csv');z=early.groupby('replicate').difference.mean();se=z.std(ddof=1)/np.sqrt(len(z));lo,hi=z.mean()+np.array([-1,1])*stats.t.ppf(.975,len(z)-1)*se
pd.DataFrame([dict(endpoint='mean1h2h_bound_minus_unbound',n_trajectories=len(z),estimate=z.mean(),low=lo,high=hi,p=stats.ttest_1samp(z,0).pvalue,status='fixed secondary sensitivity')]).to_csv(out/'ER_bound_import_early_sensitivity.csv',index=False)
rna=tidy_matrix(pd.read_csv(root/'revision8/analysis/GSE78167_TPM.csv.gz',index_col=0));cross=pd.read_csv(root/'revision8/analysis/GSE78167_sample_crosswalk.csv');log=np.log2(rna+1);rows=[]
for g in imp:
 effects=[]
 for hour in sorted(cross.hour.unique()):
  delta=[]
  for rep,z in cross.groupby('replicate'):
   now=z.loc[z.hour==hour,'sample'].iloc[0];base=z.loc[z.hour==0,'sample'].iloc[0];delta.append(log.at[g,now]-log.at[g,base])
  effects.append((hour,float(np.mean(delta))))
 onset=next((effects[i][0] for i in range(len(effects)-1) if effects[i][1]>=.2 and effects[i+1][1]>=.2),np.nan)
 rows.append(dict(gene=g,E2_bound=gene.at[g,'E2_promoter_or_distal'],onset_hour=onset,definition='first>=0.2 log2 followed by another>=0.2 timepoint; descriptive'))
pd.DataFrame(rows).to_csv(out/'ER_import_gene_descriptive_onsets.csv',index=False)
print('Import NRF1 promoter support',int(gene.loc[imp].NRF1_promoter.sum()),'/',len(imp));print(gene.reindex(positive)[['E2_promoter','E2_distal','E2_promoter_or_distal']].to_string());print(pd.DataFrame(rows).to_string(index=False))
