from pathlib import Path
import sys,json,re,io,gzip,csv,itertools,warnings
r=Path(__file__).parent;root=r.parent;old=root/'revision2';out=r/'analysis'
sys.path.insert(0,str(root/'python_packages'))
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,ttest_1samp,spearmanr
warnings.filterwarnings('ignore',category=RuntimeWarning)
mods=json.loads((root/'source/hummel2/gene_sets/mito_modules.json').read_text())
controls={
 'Proliferation':'MKI67 PCNA TOP2A UBE2C CCNB1 CCNB2 CDC20 CDK1 MCM2 MCM3 MCM4 MCM5 MCM6 MCM7',
 'ER_programme':'ESR1 PGR GREB1 TFF1 PDZK1 CA12 XBP1',
 'Leukocyte':'PTPRC CD3D CD3E CD8A CD79A MS4A1 LST1 TYROBP',
 'Epithelial':'EPCAM KRT8 KRT18 KRT19',
 'Hypoxia':'CA9 VEGFA SLC2A1 LDHA PDK1 BNIP3 NDRG1 EGLN3',
 'Interferon':'STAT1 IRF1 ISG15 IFIT1 IFIT2 IFIT3 MX1 OAS1'}
overlap=set(sum(mods.values(),[]));controls={k:[g for g in v.split() if g not in overlap] for k,v in controls.items()}
(r/'analysis/control_gene_sets.json').write_text(json.dumps(controls,indent=2))
coverage=[];diagnostics=[];matrices={};metas={};score_frames={};tests=[];deltas=[];genes=[];robust=[];associations=[];eligibility=[]
rng=np.random.default_rng(20261024)

def bh(p):
 p=np.asarray(p,float);p=np.nan_to_num(p,nan=1);ix=np.argsort(p);q=np.empty(len(p));q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def matrix(acc):
 lines=[];info={};inside=False
 with gzip.open(r/'metadata'/f'{acc}_series_matrix.txt.gz','rt',encoding='utf-8') as f:
  for l in f:
   if l.startswith('!Sample_'):
    v=next(csv.reader([l],delimiter='\t'));key=v[0][8:]
    if key in info:key=key+'_'+str(sum(k.startswith(key) for k in info))
    info[key]=v[1:]
   elif l.startswith('!series_matrix_table_begin'):inside=True
   elif l.startswith('!series_matrix_table_end'):break
   elif inside:lines.append(l)
 md=pd.DataFrame(info);md=md.rename(columns={'geo_accession':'sample'}).set_index('sample',drop=False)
 a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0) if len(lines)>1 else pd.DataFrame()
 md.to_csv(out/f'{acc}_sample_metadata.csv',index=False);metas[acc]=md
 return a,md
def platform(acc):
 compact=r/'metadata'/f'{acc}_annotation.soft.txt';large=r/'metadata'/f'{acc}_family.soft.gz'
 f=compact.open(encoding='utf-8') if compact.exists() and compact.stat().st_size>10000 else gzip.open(large,'rt',encoding='utf-8')
 lines=[]
 with f:
  for l in f:
   if l.startswith('!platform_table_begin'):
    for ll in f:
     if ll.startswith('!platform_table_end'):break
     lines.append(ll)
    break
 assert len(lines)>1,acc
 return pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str,keep_default_na=False,low_memory=False)
def array_genes(a,gpl,col,acc):
 b=platform(gpl);s=b[col].astype(str).str.strip();valid=~s.str.contains(r'///|//|;|\s',regex=True)&s.ne('')&s.ne('---')
 mp=b.loc[valid,['ID',col]].drop_duplicates('ID').set_index('ID')[col].str.upper()
 mp.to_csv(out/f'{acc}_probe_to_gene.csv');a.index=a.index.astype(str);a['gene']=mp.reindex(a.index);return a.dropna(subset=['gene']).groupby('gene').median()
def score(acc,a,transform,category):
 a=a.apply(pd.to_numeric,errors='coerce');a.index=a.index.astype(str).str.upper();a=a.groupby(level=0).median()
 diagnostics.append({'accession':acc,'transform':transform,'category':category,'genes_before_complete_filter':len(a),'samples':a.shape[1],'min':float(a.min().min()),'median':float(np.nanmedian(a.values)),'max':float(a.max().max())})
 if transform=='log2p1_TPM':a=np.log2(a+1)
 elif transform=='log2p1_CPM':a=np.log2(a.div(a.sum(axis=0),axis=1)*1e6+1)
 a=a.dropna(axis=0,how='any');matrices[acc]=a
 a.to_pickle(out/f'{acc}_gene_expression.pkl')
 # Gene centering affects displays but cancels from pair differences.
 centered=a.sub(a.mean(axis=1),axis=0);ranks=a.rank(axis=0,pct=True)
 for scale,v in [('source',centered),('rank',ranks)]:
  scores=pd.DataFrame(index=a.columns)
  for name,gs in {**mods,**controls}.items():
   mapped=[g.upper() for g in gs if g.upper() in a.index];eligible=len(mapped)>=max(3,np.ceil(.5*len(gs)))
   if scale=='source':coverage.append({'accession':acc,'module':name,'n_original':len(gs),'n_mapped':len(mapped),'eligible':bool(eligible),'mapped_genes':';'.join(mapped),'missing_genes':';'.join(g for g in gs if g.upper() not in a.index)})
   if eligible:scores[name]=v.loc[mapped].mean()
  scores.index.name='sample';score_frames[acc,scale]=scores;scores.to_csv(out/f'{acc}_{scale}_scores.csv')
 print('SCORED',acc,a.shape,'eligible',len(score_frames[acc,'source'].columns),flush=True)
def boot(d):return np.quantile(rng.choice(np.asarray(d),size=(5000,len(d)),replace=True).mean(axis=1),[.025,.975])
def clinical(acc,md):
 assert md[['patient','timepoint']].duplicated().sum()==0
 x=md.pivot(index='patient',columns='timepoint',values='sample')[['baseline','on']].dropna();x.to_csv(out/f'{acc}_pair_crosswalk.csv')
 for scale in ['source','rank']:
  scores=score_frames[acc,scale];dd=scores.loc[x.on].set_axis(x.index)-scores.loc[x.baseline].set_axis(x.index)
  dd.to_csv(out/f'{acc}_{scale}_paired_deltas.csv')
  for module in dd:
   d=dd[module];ci=boot(d);p=wilcoxon(d).pvalue if (d!=0).any() else 1.
   tests.append({'accession':acc,'contrast':'on_minus_baseline','module':module,'scale':scale,'family':'clinical_mito' if module in mods else 'clinical_controls','n0':len(d),'n1':len(d),'mean_delta':d.mean(),'median_delta':d.median(),'n_positive':int((d>0).sum()),'p':p,'ci_low':ci[0],'ci_high':ci[1],'test':'paired_Wilcoxon'})
   deltas.extend({'accession':acc,'patient':patient,'module':module,'scale':scale,'delta':v} for patient,v in d.items())
   if module in mods:
    for patient in dd.index:robust.append({'accession':acc,'module':module,'scale':scale,'omitted_patient':patient,'remaining_n':len(d)-1,'mean_delta':d.drop(patient).mean()})
  if scale=='source':
   for target in ['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV']:
    if target not in dd:continue
    for control in controls:
     if control not in dd:continue
     rho,p=spearmanr(dd[target],dd[control]);associations.append({'accession':acc,'target':target,'control':control,'n':len(dd),'rho':rho,'p':p,'interpretation':'paired-change expression association; not a causal or cell-proportion estimate'})
   # Partial rank association with proliferation, controlling epithelial/leukocyte marker deltas when eligible.
   if 'Import_TIM_TOM' in dd and 'Proliferation' in dd:
    covars=[v for v in ['Epithelial','Leukocyte'] if v in dd]
    rr=dd[['Import_TIM_TOM','Proliferation']+covars].rank();design=np.column_stack([np.ones(len(rr)),rr[covars].values])
    a=rr['Import_TIM_TOM'].to_numpy(copy=True);b=rr['Proliferation'].to_numpy(copy=True);a-=design@np.linalg.lstsq(design,a,rcond=None)[0];b-=design@np.linalg.lstsq(design,b,rcond=None)[0]
    rho=np.corrcoef(a,b)[0,1];ps=(1+sum(abs(np.corrcoef(a,rng.permutation(b))[0,1])>=abs(rho) for _ in range(4999)))/5000
    associations.append({'accession':acc,'target':'Import_TIM_TOM','control':'Proliferation_partial_'+','.join(covars),'n':len(rr),'rho':rho,'p':ps,'interpretation':'exploratory residual-rank permutation; marker adjustment does not establish cellular origin'})
 a=matrices[acc];d=a[x.on].set_axis(x.index,axis=1)-a[x.baseline].set_axis(x.index,axis=1);stat=ttest_1samp(d.values,0,axis=1,nan_policy='omit');p=np.nan_to_num(stat.pvalue,nan=1)
 g=pd.DataFrame({'gene':a.index,'mean_delta':d.mean(axis=1).values,'n_positive':(d>0).sum(axis=1).values,'n_pairs':len(x),'p_paired_t':p,'q_genomewide':bh(p)});g['accession']=acc;g.to_csv(out/f'{acc}_all_gene_paired_tests.csv',index=False)
 for module,gs in mods.items():
  z=g[g.gene.isin(gs)].copy();z['module']=module;genes.extend(z.to_dict('records'))
 return len(x)
def exact(a,b):
 v=np.r_[a,b];n=len(a);obs=abs(np.mean(b)-np.mean(a));count=0;total=0
 for ix in itertools.combinations(range(len(v)),n):
  mask=np.zeros(len(v),bool);mask[list(ix)]=True;d=abs(v[~mask].mean()-v[mask].mean());count+=d>=obs-1e-12;total+=1
 return count/total,total
def experiment(acc,groups,contrasts,infer=True):
 pd.DataFrame({'sample':list(groups),'group':list(groups.values())}).to_csv(out/f'{acc}_group_manifest.csv',index=False)
 for control,treatment in contrasts:
  c=[k for k,v in groups.items() if v==control];t=[k for k,v in groups.items() if v==treatment]
  for scale in ['source','rank']:
   scores=score_frames[acc,scale]
   for module in scores.columns.intersection(mods):
    a=scores.loc[c,module].values;b=scores.loc[t,module].values
    p,total=exact(a,b) if infer and min(len(a),len(b))>=3 else (np.nan,0)
    tests.append({'accession':acc,'contrast':treatment+'_minus_'+control,'module':module,'scale':scale,'family':'experimental_mito','n0':len(a),'n1':len(b),'mean_delta':np.mean(b)-np.mean(a),'median_delta':np.median(b)-np.median(a),'p':p,'permutations':total,'test':'unpaired_exact_permutation' if total else 'descriptive_only'})

# Metadata eligibility is fixed before the following scores are inspected.
eligibility=[
 {'accession':'GSE153470','decision':'clinical paired','reason':'source-defined anastrozole baseline/surgery; pairs mapped by numeric title; response/Ki67 labels not deposited in matrix'},
 {'accession':'GSE87411','decision':'clinical paired','reason':'Z1031B pre/2-4-week AI; no new counting of reused GSE29442/GSE35186/GSE35191 arrays as separate cohorts; common-reference log2 ratios retained'},
 {'accession':'GSE155192','decision':'experimental context','reason':'source-defined stimulation/hypoxia, terminally exhausted TIL from Ndufs4-deficient B16 tumours, and PGC1a-overexpression CD8 TIL; PGC sample names mapped by deposited Sample_description, not expression'},
 {'accession':'GSE67916','decision':'resistance context','reason':'four TamR derivatives, shared parental MCF7; n3/3/2/2 versus eight parental replicates; TamR7/8 descriptive, not pooled pseudoreplication'},
 {'accession':'GSE125738','decision':'resistance context','reason':'three resistant/three parental T47D libraries; deposited per-replicate integer columns documented counts; CPM log2p1, no use of fold-change/FDR columns'},
 {'accession':'GSE118713','decision':'resistance context','reason':'three each MCF7/TAMR/FASR libraries; sum transcript TPM per gene; no selection of isoforms by effect'},
 {'accession':'GSE129544','decision':'excluded from inference','reason':'one parental and one resistant T47D sample; insufficient replication; no independent validation claim'},
 {'accession':'GSE21618','decision':'excluded from expansion inference','reason':'complex ligand/time experiment with six resistant derivative controls but only one matched WT Control; no pooling timepoints or perturbations to create independent replicates'},
 {'accession':'GSE41998','decision':'excluded longitudinal','reason':'all deposited tumour biopsies taken before treatment'}]
pd.DataFrame(eligibility).to_csv(out/'candidate_eligibility_log.csv',index=False)

for acc,gpl,col in [('GSE153470','GPL13376','Symbol'),('GSE87411','GPL6480','GENE_SYMBOL'),('GSE67916','GPL570','Gene Symbol')]:
 a,md=matrix(acc);a=array_genes(a,gpl,col,acc)
 if acc=='GSE153470':
  s=md.title.str.extract(r'(\d+)([AB])$');md['patient']=s[0];md['timepoint']=s[1].map({'A':'baseline','B':'on'})
 elif acc=='GSE87411':
  md['patient']=md.title.str.extract(r'ID_(\d+)')[0];md['timepoint']=np.where(md.title.str.contains(r'_(?:W2|2wk|M)(?:_| |$)',case=False),'on','baseline')
 md.to_csv(out/f'{acc}_sample_metadata.csv',index=False)
 score(acc,a,'deposited_VST_RSN' if acc=='GSE153470' else 'deposited_log2', 'clinical' if acc!='GSE67916' else 'experiment')
 if acc!='GSE67916':print('PAIRS',acc,clinical(acc,md),flush=True)
 else:
  groups={k:re.match(r'(TamR\d|MCF7)',v).group(1) if re.match(r'(TamR\d|MCF7)',v) else 'parental' for k,v in md.title.items()}
  groups={k:'parental' if v=='MCF7' else v for k,v in groups.items()}
  experiment(acc,groups,[('parental',v) for v in ['TamR1','TamR4','TamR7','TamR8']])

a,md=matrix('GSE125738');raw=pd.read_csv(r/'metadata/GSE125738_fpkm_T47D-TR_vs_T47D.txt.gz',sep='\t',index_col=0)
columns=['T47D_1','T47D_2','T47D_3','T47D-TR_1','T47D-TR_2','T47D-TR_3'];a=raw[columns];assert np.allclose(a.values,np.rint(a.values))
score('GSE125738',a,'log2p1_CPM','experiment');experiment('GSE125738',{k:'resistant' if '-TR_' in k else 'parental' for k in a},[('parental','resistant')])
a,md=matrix('GSE118713');raw=pd.read_csv(r/'metadata/GSE118713_TRANSCRIPT_TAMR_FASR_MCF7_expression_data.tsv.gz',sep='\t');columns=[k for k in raw if re.match(r'Rep[123]\.(TAMR|FASR|MCF7)\.TPM',k)]
a=raw.dropna(subset=['gene.symbol']).groupby('gene.symbol')[columns].sum();score('GSE118713',a,'log2p1_TPM','experiment');experiment('GSE118713',{k:k.split('.')[1] for k in a},[('MCF7','TAMR'),('MCF7','FASR')])

_,md=matrix('GSE155192')
for name in ['ChronicStimHypoxia','Ndufs4KO','PGC1aOE']:
 acc='GSE155192_'+name;a=pd.read_csv(r/'metadata'/f'GSE155192_TPM_{name}.csv.gz',index_col=0)
 if name=='ChronicStimHypoxia':
  groups={k:re.sub(r'_M\d_.*$','',k).replace('Chronic1.10','Chronic') for k in a};contrasts=[('Acute_Norm','Acute_Hypoxia'),('Acute_Norm','Chronic_Norm'),('Chronic_Norm','Chronic_Hypoxia'),('Acute_Hypoxia','Chronic_Hypoxia')]
 elif name=='Ndufs4KO':
  # The source perturbation is tumour Ndufs4 deficiency, NOT a CD8-cell knockout.
  a=a[[k for k in a if k.endswith('_Texh')]];groups={k:'Ndufs4_def_tumour' if k.startswith('Ndufs4_ko') else 'WT_tumour' for k in a};contrasts=[('WT_tumour','Ndufs4_def_tumour')]
 else:
  rows=md[md.title.str.startswith(('Empty Vector','PGC1a'))];desc=[k for k in md if k.startswith('description')][0];groups={}
  cross=[]
  for k in a:
   z=rows[rows[desc].str.startswith(k+'_')];assert len(z)==1,(k,z)
   group='PGC1a_OE' if z.title.iloc[0].startswith('PGC1a') else 'EV';groups[k]=group;cross.append({'processed_column':k,'GSM':z.index[0],'title':z.title.iloc[0],'source_description':z[desc].iloc[0],'group':group})
  pd.DataFrame(cross).to_csv(out/'PGC1a_verified_sample_crosswalk.csv',index=False);contrasts=[('EV','PGC1a_OE')]
 score(acc,a,'log2p1_TPM','experiment');experiment(acc,groups,contrasts)

tt=pd.DataFrame(tests)
for _,ix in tt[tt.p.notna()].groupby(['family','scale']).groups.items():tt.loc[ix,'q']=bh(tt.loc[ix,'p'])
tt.to_csv(out/'expansion_all_module_tests.csv',index=False)
pd.DataFrame(deltas).to_csv(out/'clinical_patient_deltas.csv',index=False)
pd.DataFrame(genes).to_csv(out/'clinical_module_gene_tests.csv',index=False)
pd.DataFrame(robust).to_csv(out/'leave_one_patient_out.csv',index=False)
aa=pd.DataFrame(associations);aa['q']=bh(aa.p);aa.to_csv(out/'clinical_control_associations.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'gene_coverage.csv',index=False)
(out/'processing_diagnostics.json').write_text(json.dumps(diagnostics,indent=2))
print(tt[(tt.module.isin(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV']))&(tt.scale=='source')][['accession','contrast','module','n0','n1','mean_delta','p','q']].to_string(index=False))
print('ALL_TESTS',len(tt),'inferential',tt.p.notna().sum(),'CONTROL_ASSOC',len(aa),flush=True)
