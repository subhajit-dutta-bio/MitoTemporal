from analysis_utils import *
rng=np.random.default_rng(20261025);results=[];sets_out=[];diagnostics=[];null=[]
pc=set(h.loc[h.locus_group=='protein-coding gene','symbol']);exclude=set(mods['Import_TIM_TOM'])|set(controls['ER_programme'])|set(controls['Proliferation'])
for acc in ['GSE20181','GSE111563','GSE105777']:
 a,p=clinical(acc);base=a.loc[:,p.baseline];valid=(base.notna().mean(axis=1)>=.8)&(base.std(axis=1)>0);features=pd.DataFrame({'mean':base.mean(axis=1),'sd':base.std(axis=1)});imp=[g for g in mods['Import_TIM_TOM'] if g in a.index and valid.loc[g]];assert len(imp)>=12
 actual={}
 for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
  ds=paired_scores(v,p,{'Import_TIM_TOM':imp});d=ds.Import_TIM_TOM;actual[scale]=d[ds.arm=='AI'].mean()-(d[ds.arm=='Untreated'].mean() if acc=='GSE105777' else 0)
 for universe in ['all_coding','nonimport_mito']:
  allowed=pc if universe=='all_coding' else set(comparators['Nonimport_mito']);candidates=[g for g in a.index if g in allowed-exclude and valid.loc[g]]
  feat=features.loc[candidates+imp];mu=feat.mean();sd=feat.std();F=(feat-mu)/sd
  nearest={g:((F.loc[candidates]-F.loc[g])**2).sum(axis=1).sort_values(kind='stable').index[:20].tolist() for g in imp};union=sorted(set(imp)|set(q for gs in nearest.values() for q in gs));cor=base.loc[union].T.corr(min_periods=max(3,int(np.ceil(.8*len(p)))));target=cor.loc[imp,imp].to_numpy()[np.triu_indices(len(imp),1)].mean();proposals=[];tries=0
  while len(proposals)<5000 and tries<20000:
   tries+=1;chosen=[]
   for g in imp:
    options=[q for q in nearest[g] if q not in chosen]
    if not options:break
    chosen.append(str(rng.choice(options)))
   if len(chosen)!=len(imp):continue
   cc=cor.loc[chosen,chosen].to_numpy()[np.triu_indices(len(imp),1)];avg=np.nanmean(cc)
   if not np.isfinite(avg):continue
   distance=float(np.mean([np.sqrt(((F.loc[g]-F.loc[q])**2).sum()) for g,q in zip(imp,chosen)]));proposals.append((abs(avg-target),avg,distance,chosen))
  selected=sorted(proposals,key=lambda x:x[0])[:1000];eligible=len(selected)>=500;diagnostics.append(dict(accession=acc,universe=universe,n_baseline_patients=len(p),n_candidates=len(candidates),n_import=len(imp),proposals=len(proposals),selected=len(selected),target_baseline_correlation=target,median_selected_correlation=np.median([v[1] for v in selected]),median_absolute_correlation_mismatch=np.median([v[0] for v in selected]),median_expression_SD_distance=np.median([v[2] for v in selected])))
  for i,(dist,corr,fd,genes) in enumerate(selected):sets_out.append(dict(accession=acc,universe=universe,set_id=i,baseline_correlation=corr,correlation_mismatch=dist,expression_SD_distance=fd,genes=';'.join(genes)))
  union=sorted(set(g for v in selected for g in v[3]));ix={g:i for i,g in enumerate(union)};arrs={}
  for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
   deltas=v.loc[union,p.on].to_numpy()-v.loc[union,p.baseline].to_numpy();vals=[]
   for i,entry in enumerate(selected):
    data=deltas[[ix[g] for g in entry[3]]];count=np.isfinite(data).sum(axis=0);mean=np.nanmean(data,axis=0);mean[count<max(3,int(np.ceil(len(imp)/2)))]=np.nan;effect=np.nanmean(mean[p.arm.eq('AI')])-(np.nanmean(mean[p.arm.eq('Untreated')]) if acc=='GSE105777' else 0);vals.append(effect);null.append(dict(accession=acc,universe=universe,scale=scale,set_id=i,contrast=effect))
   vals=np.array(vals);obs=actual[scale];pv=min(1,2*min((np.sum(vals<=obs)+1)/(len(vals)+1),(np.sum(vals>=obs)+1)/(len(vals)+1))) if eligible else np.nan;results.append(dict(accession=acc,universe=universe,scale=scale,observed_import_contrast=obs,n_import=len(imp),import_genes=';'.join(imp),null_median=np.nanmedian(vals),null_low=np.nanquantile(vals,.025),null_high=np.nanquantile(vals,.975),percentile=np.mean(vals<=obs),p=pv,n_selected=len(vals),inference='competitive_gene_set_context_not_patient_replication'))
  print(acc,universe,'matched',len(selected),flush=True)
res=family(pd.DataFrame(results),['scale'],6);freeze(res,'matched_gene_benchmark_tests.csv');pd.DataFrame(sets_out).to_csv(out/'matched_gene_sets.csv.gz',index=False,compression='gzip');freeze(pd.DataFrame(diagnostics),'matched_gene_matching_diagnostics.csv');pd.DataFrame(null).to_csv(out/'matched_gene_null_distributions.csv.gz',index=False,compression='gzip');print(res.to_string(index=False))
