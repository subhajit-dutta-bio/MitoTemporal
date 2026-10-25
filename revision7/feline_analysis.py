from common import *
from scipy.stats import wilcoxon,mannwhitneyu,spearmanr
src=pd.read_csv(r/'metadata/43018_2021_215_MOESM6_ESM.csv');src['patient']='P'+src['Sample'].str.extract(r'FEL0*(\d+)P')[0].str.zfill(2)
cols=['patient','ARM','Treatment','dynamic_class3'];clinical=src[cols].drop_duplicates();assert not clinical.patient.duplicated().any()
ext=pd.read_csv(r/'metadata/43018_2021_215_MOESM10_ESM.csv');ext['patient']='P'+ext['Sample'].str.extract(r'FEL0*(\d+)P')[0].str.zfill(2);ext=ext.rename(columns={'Dynamic.Biological.Response':'dynamic_class3','Clinical.Response':'clinical_response'})
both=clinical.merge(ext[['patient','ARM','dynamic_class3']],on='patient',suffixes=('_figure3','_extended1'));assert (both.ARM_figure3==both.ARM_extended1).all() and (both.dynamic_class3_figure3==both.dynamic_class3_extended1).all()
clinical=pd.concat([clinical,ext[cols]],ignore_index=True).drop_duplicates();assert not clinical.patient.duplicated().any();clinical=clinical.merge(ext[['patient','clinical_response','prop_change']],on='patient',how='left');clinical.to_csv(out/'FELINE_clinical_crosswalk.csv',index=False)
tm=src[['Timepoint','Day']].drop_duplicates();assert dict(zip(tm.Timepoint,tm.Day))==dict(S=0,M=14,E=180);tm.to_csv(out/'FELINE_timepoint_crosswalk.csv',index=False)
selected=set(sum([mods[k] for k in ['Import_TIM_TOM','CIV_nuclear','Mitophagy']],[])+sum(controls.values(),[])+['NRF1','TFAM','PPARGC1A']);cover=[];parts=[];pseudoparts=[]
for platform in ['10x','icell8']:
 cache=out/f'FELINE_{platform}_selected_counts.npz';label=pd.read_csv(r/'metadata'/f'GSE158724_{platform}.cell_types.txt.gz',sep='\t').set_index('Cell')
 if not cache.exists():
  with gzip.open(r/'metadata'/f'GSE158724_{platform}.counts.txt.gz','rt') as f:
   cells=next(f).rstrip().split('\t')[1:];ix=np.array([i for i,c in enumerate(cells) if c in label.index]);cells=np.array(cells)[ix];lib=np.zeros(len(cells));gene=[];count=[];pseudo=[];pg=[]
   groups=pd.Series(cells).str.extract(r'^(P\d+_[SME])')[0];ug=sorted(groups.unique());codes=pd.Categorical(groups,categories=ug).codes
   for j,line in enumerate(f):
    g,vs=line.rstrip().split('\t',1);g=canon(g);v=np.fromstring(vs,sep='\t',dtype=float)[ix];lib+=v
    pg.append(g);pseudo.append(np.bincount(codes,weights=v,minlength=len(ug)))
    if g in selected:gene.append(g);count.append(v.astype('float32'))
    if j%5000==0:print(platform,j,'genes',len(cells),'source-labelled cancer cells',flush=True)
  np.savez_compressed(cache,cells=cells,genes=np.array(gene),counts=np.array(count),library=lib)
  z=pd.DataFrame(pseudo,index=pg,columns=ug).groupby(level=0).sum();z.to_csv(out/f'FELINE_{platform}_pseudobulk_counts.csv.gz',compression='gzip')
 z=np.load(cache);a=pd.DataFrame(z['counts'],index=z['genes'],columns=z['cells']).groupby(level=0).sum();lib=pd.Series(z['library'],index=z['cells']);log=np.log1p(a.div(lib,axis=1)*1e6);s=scores(log);s['cell']=s.index;s['patient']=s.cell.str.extract(r'^(P\d+)')[0];s['time']=s.cell.str.extract(r'^P\d+_([SME])')[0];s['platform']=platform;s['library']=lib;s['n_cells']=1
 zscore=s.drop(columns='cell').groupby(['patient','time','platform']).agg({**{c:'mean' for c in s.columns if c not in ['cell','patient','time','platform','n_cells']},'n_cells':'sum'}).reset_index();parts.append(zscore)
 for mod in ['Import_TIM_TOM','CIV_nuclear','Mitophagy']:cover.append(dict(platform=platform,module=mod,n_module=len(mods[mod]),n_measured=len(set(mods[mod])&set(a.index)),genes=';'.join(sorted(set(mods[mod])&set(a.index)))))
 p=pd.read_csv(out/f'FELINE_{platform}_pseudobulk_counts.csv.gz',index_col=0);v=np.log1p(p.div(p.sum(),axis=1)*1e6);ps=scores(v);ps['sample']=ps.index;ps['patient']=ps['sample'].str.extract(r'^(P\d+)')[0];ps['time']=ps['sample'].str.extract(r'^P\d+_([SME])')[0];ps['platform']=platform;ps['scale']='pseudobulk_logCPM';pseudoparts.append(ps)
tab=pd.concat(parts,ignore_index=True).merge(clinical,on='patient',how='left',validate='many_to_one');assert tab.ARM.notna().all();tab.to_csv(out/'FELINE_patient_time_scores.csv',index=False);pd.DataFrame(cover).to_csv(out/'FELINE_gene_coverage.csv',index=False)
pd.concat(pseudoparts,ignore_index=True).merge(clinical,on='patient',validate='many_to_one').to_csv(out/'FELINE_pseudobulk_scores.csv',index=False)
rows=[];deltas=[];response=[]
for platform in ['10x','icell8']:
 for arm in ['A','B','C']:
  eligible=tab[(tab.platform==platform)&(tab.ARM==arm)&(tab.n_cells>=10)]
  for mod in ['Import_TIM_TOM','CIV_nuclear']:
   p=eligible.pivot(index='patient',columns='time',values=mod)
   for post in ['M','E']:
    if not {'S',post}<=set(p.columns):continue
    v=p[['S',post]].dropna();d=v[post]-v.S
    for pt,val in d.items():deltas.append(dict(patient=pt,platform=platform,arm=arm,module=mod,contrast=post+'_minus_S',pre=v.loc[pt,'S'],post=v.loc[pt,post],delta=val,response=clinical.set_index('patient').loc[pt,'dynamic_class3']))
    if len(d)<3:continue
    rows.append(dict(platform=platform,arm=arm,module=mod,contrast=post+'_minus_S',n=len(d),mean_delta=d.mean(),median_delta=d.median(),n_positive=int((d>0).sum()),p=wilcoxon(d).pvalue if np.any(d!=0) else 1.,family='primary_10x_letrozole' if platform=='10x' and arm=='A' else 'other_arm_platform_exploratory'))
    if platform=='10x' and arm=='A' and post=='E':
     labels=clinical.set_index('patient').loc[d.index,'dynamic_class3'];a=d[labels=='Response'];b=d[labels=='Non-response'];response.append(dict(module=mod,n_sensitive=len(a),n_resistant=len(b),mean_sensitive=a.mean(),mean_resistant=b.mean(),difference=a.mean()-b.mean(),p=mannwhitneyu(a,b,alternative='two-sided').pvalue if min(len(a),len(b))>=2 else np.nan))
tests=pd.DataFrame(rows)
for _,ix in tests.groupby('family').groups.items():tests.loc[ix,'q']=bh(tests.loc[ix,'p']);tests.loc[ix,'family_size']=len(ix)
tests.to_csv(out/'FELINE_all_tests.csv',index=False);pd.DataFrame(deltas).to_csv(out/'FELINE_patient_deltas.csv',index=False)
response=pd.DataFrame(response);response['q']=bh(response.p);response.to_csv(out/'FELINE_response_tests.csv',index=False)
print(tests.to_string(index=False));print(response.to_string(index=False))
