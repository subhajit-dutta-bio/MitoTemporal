from analysis_utils import *
subgroups={'TOM': [g for g in mods['Import_TIM_TOM'] if g.startswith('TOMM')],'TIM_PAM':[g for g in mods['Import_TIM_TOM'] if not g.startswith('TOMM')]}
gene_groups=dict(groups,**subgroups);cov=[];audit=[];long=[];tests=[];spec=[];summaries=[]
def prepare(file,genecol,ratio_cols,orientation):
 x=pd.read_csv(file,sep='\t',low_memory=False);original=len(x);bad=np.zeros(len(x),bool)
 for c in ['Reverse','Only identified by site','Contaminant','Potential contaminant']:
  if c in x:bad |= x[c].fillna('').astype(str).str.strip().eq('+').to_numpy()
 def gene(s):
  names=re.findall(r'\bGN=([^\s;]+)',str(s)) if genecol=='Fasta headers' else re.split(r';|,',str(s));names={canon(g.strip()) for g in names if canon(g.strip()) in approved};return next(iter(names)) if len(names)==1 else None
 x['approved_symbol']=x[genecol].map(gene);amb=x.approved_symbol.isna();audit.append(dict(file=file.name,source_rows=original,flagged_rows=int(bad.sum()),ambiguous_or_unmapped_rows=int(amb.sum()),retained_rows=int((~bad&~amb).sum()),ratio_columns=';'.join(ratio_cols),orientation=orientation))
 y=x.loc[~bad&~amb].copy();a=y[ratio_cols].apply(pd.to_numeric,errors='coerce');a=a.where(a>0);a=orientation*np.log2(a);a.index=y.approved_symbol;a=a.groupby(level=0).median();return a
for hour in [6,36]:
 file=next((r/'metadata').glob(f'PXD002797_MaxQuant_results_{hour}h*proteinGroups.txt'));cols=['Ratio H/L Normalized '+rep for rep in ['A','B','C']];a=prepare(file,'Gene Names',cols,1);a.columns=['A','B','C'];a.to_csv(out/f'PXD002797_{hour}h_log2_E2_vehicle.csv.gz',compression='gzip')
 for mod,gs in gene_groups.items():
  hits=[g for g in gs if g in a.index and a.loc[g].notna().all()];required=max(3,int(np.ceil(len(gs)/2)));available=len(hits)>=required
  cov.append(dict(study='PXD002797',condition=str(hour)+'h',module=mod,n_fixed=len(gs),n_measured_all_replicates=len(hits),minimum=required,eligible=available,genes=';'.join(hits)))
  vals=a.loc[hits].mean() if available else pd.Series(np.nan,index=a.columns)
  for rep in a:long.append(dict(study='PXD002797',condition=str(hour)+'h',replicate=rep,module=mod,score=vals[rep],n_genes=len(hits)))
  result=onesample(vals) if available else dict(n_treated=3,n_control=0,estimate=np.nan,low=np.nan,high=np.nan,p=np.nan,df=2)
  tests.append(dict(study='PXD002797',hour=hour,module=mod,family='primary_import' if mod=='Import_TIM_TOM' else 'secondary_submodule' if mod in subgroups else 'context_'+mod,availability='measured' if available else 'insufficient_fixed_gene_coverage',**result))
 sc=pd.DataFrame([v for v in long if v['study']=='PXD002797' and v['condition']==str(hour)+'h']).pivot(index='replicate',columns='module',values='score')
 for gene in mods['Import_TIM_TOM']:
  for rep in a:long.append(dict(study='PXD002797_gene',condition=str(hour)+'h',replicate=rep,module=gene,score=a.at[gene,rep] if gene in a.index else np.nan,n_genes=1))
 for comp in comparators:
  d=sc.Import_TIM_TOM-sc[comp];spec.append(dict(hour=hour,comparator=comp,**onesample(d)))
for frac in ['Cyto','Nuclear','membrane']:
 file=r/'metadata'/f'PXD049343_proteinGroups_{frac}.txt';a=prepare(file,'Fasta headers',['Ratio H/L normalized'],-1);a.columns=['summary'];a.to_csv(out/f'PXD049343_{frac}_log2_E2_vehicle.csv.gz',compression='gzip')
 for mod,gs in gene_groups.items():
  hits=[g for g in gs if g in a.index and pd.notna(a.at[g,'summary'])];required=max(3,int(np.ceil(len(gs)/2)));available=len(hits)>=required
  cov.append(dict(study='PXD049343',condition=frac,module=mod,n_fixed=len(gs),n_measured_all_replicates=len(hits),minimum=required,eligible=available,genes=';'.join(hits)));summaries.append(dict(study='PXD049343',fraction=frac,module=mod,n_measured=len(hits),estimate=a.loc[hits,'summary'].mean() if available else np.nan,measured_subset_descriptive_mean=a.loc[hits,'summary'].mean() if hits else np.nan,p=np.nan,reason='single_fraction_summary_no_biological_replicate_ratios'))
  for g in hits:long.append(dict(study='PXD049343_gene',condition=frac,replicate='summary',module=g,score=a.at[g,'summary'],n_genes=1))
res=family(pd.DataFrame(tests),['family']);sp=pd.DataFrame(spec);sp['q']=bh(sp.p.fillna(1));sp.loc[sp.p.isna(),'q']=np.nan;sp['family_size']=8
freeze(res,'perturbational_proteomics_tests.csv');freeze(pd.DataFrame(cov),'perturbational_proteomics_coverage.csv');freeze(pd.DataFrame(audit),'perturbational_proteomics_input_audit.csv');freeze(pd.DataFrame(long),'perturbational_proteomics_replicate_scores.csv');freeze(sp,'perturbational_proteomics_specificity.csv');freeze(pd.DataFrame(summaries),'PXD049343_descriptive_scores.csv')
print(res[res.module.isin(['Import_TIM_TOM','TOM','TIM_PAM'])].to_string(index=False));print(pd.DataFrame(cov).query('module=="Import_TIM_TOM"').to_string(index=False))
