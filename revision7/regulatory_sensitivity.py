from common import *
tss=pd.read_csv(out/'NRF1_GRCh38_TSS_binding.csv.gz');a=pd.read_csv(out/'GROseq_gene_profiles.csv.gz',index_col=0);bound=tss.groupby('gene').NRF1_bound.any();nt=tss.groupby('gene').size();universe=sorted(set(a.index[a.VEH>=1])&set(bound.index));eb=pd.qcut(a.loc[universe,'VEH'].rank(method='first'),10,labels=False);tb=pd.cut(nt.loc[universe],[0,1,3,7,100000],labels=False);strata=pd.DataFrame(dict(expression_decile=eb,TSS_bin=tb));rng=np.random.default_rng(20261025);rows=[]
for mod in ['Import_TIM_TOM','CIV_nuclear']:
 hit=sorted(set(mods[mod])&set(universe));k=int(bound.loc[hit].sum());sim=np.zeros(10000)
 for g in hit:
  z=strata.loc[g];pool=strata.index[(strata.expression_decile==z.expression_decile)&(strata.TSS_bin==z.TSS_bin)];sim+=bound.loc[rng.choice(pool,10000,replace=True)].to_numpy(dtype=int)
 rows.append(dict(module=mod,n_module=len(hit),n_bound=k,mean_matched_bound=sim.mean(),p=(1+np.sum(sim>=k))/10001,matching='baseline_expression_decile_and_number_of_distinct_TSS_bins_1_2to3_4to7_8plus',status='post_primary_sensitivity'))
z=pd.DataFrame(rows);z['q']=bh(z.p);z.to_csv(out/'NRF1_expression_TSS_matched_sensitivity.csv',index=False);print(z.to_string(index=False))
