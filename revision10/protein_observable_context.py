from analysis_utils import *
rows=[];genes=[]
for hour in [6,36]:
 a=pd.read_csv(out/f'PXD002797_{hour}h_log2_E2_vehicle.csv.gz',index_col=0);hit=[g for g in mods['Import_TIM_TOM'] if g in a.index and a.loc[g].notna().all()]
 for rep in a:
  rows.append(dict(study='PXD002797',hour=hour,replicate=rep,n_measured=len(hit),detected_subset_mean=a.loc[hit,rep].mean(),scope='descriptive_below_prespecified_full_module_coverage_no_P'))
 for gene in mods['Import_TIM_TOM']:
  for rep in a:genes.append(dict(study='PXD002797',condition=f'{hour}h',gene=gene,replicate=rep,log2_E2_vehicle=a.at[gene,rep] if gene in a.index else np.nan,eligible_common_core=gene in hit))
for frac in ['Cyto','Nuclear','membrane']:
 a=pd.read_csv(out/f'PXD049343_{frac}_log2_E2_vehicle.csv.gz',index_col=0)
 for gene in mods['Import_TIM_TOM']:genes.append(dict(study='PXD049343',condition=frac,gene=gene,replicate='source_summary',log2_E2_vehicle=a.at[gene,'summary'] if gene in a.index else np.nan,eligible_common_core=False))
freeze(pd.DataFrame(rows),'proteomics_detected_subset_context.csv');freeze(pd.DataFrame(genes),'proteomics_import_gene_ratios.csv')
print(pd.DataFrame(rows).to_string(index=False))
