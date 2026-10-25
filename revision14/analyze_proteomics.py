from v14_common import *
import gzip,openpyxl,io
raw=r/'raw';rng=np.random.default_rng(20261025);rows=[];cover=[]
def symbol(v):
 if not isinstance(v,str):return None
 gs={canon(x) for x in re.split(r'[;,|\s]+',v.strip()) if x}
 return next(iter(gs)) if len(gs)==1 and next(iter(gs)) in approved else None
def protein_tidy(frame):
 frame=frame.copy();frame.index=[symbol(x) for x in frame.index];frame=frame.loc[frame.index.notna()];return frame.groupby(level=0).median()
def panel(a,group,assay):
 result={}
 for name,genes in {'Import_TIM_TOM':mods['Import_TIM_TOM'],'Nonimport_mito':comparators['Nonimport_mito']}.items():
  hit=[g for g in genes if g in a.index];minimum=max(3,int(np.ceil(len(genes)/2)))
  fixed=[g for g in hit if all(a.loc[g,group.index[group==label]].notna().mean()>=.8 for label in group.dropna().unique())]
  cover.append(dict(assay=assay,module=name,n_fixed=len(genes),n_in_matrix=len(hit),n_constant_panel=len(fixed),minimum=minimum,genes=';'.join(fixed)))
  result[name]=fixed
 return result
def scores_for(a,fixed):
 s=pd.DataFrame(index=a.columns)
 for name,genes in fixed.items():s[name]=module_score(a,genes,minimum=max(12 if name=='Import_TIM_TOM' else int(np.ceil(len(comparators['Nonimport_mito'])/2)),int(np.ceil(.8*len(genes))))) if genes else np.nan
 return s
def family(z):
 z=pd.DataFrame(z);z['p']=z.get('p',np.nan);z['q']=np.where(z.p.notna(),bh(z.p.fillna(1)),np.nan);z['family_slots']=3;return z
def conditional_test(df,y,x,cov=(),binary=False,n=10):
 if binary:
  counts=df.dropna(subset=[y,x,*cov]).groupby(x).size()
  if len(counts)!=2 or counts.min()<5:return dict(status='unavailable',reason='binary groups below prespecified n',n=int(counts.sum()))
 return fit_ols(df,y,x,cov,standardize=not binary,min_n=n)

# CCLE: source quantities are already normalized log2 abundance ratios.
path=Path(__import__('os').environ.get('MITOTEMPORAL_CCLE_PROTEOME',str(r/'raw/ccle_prot/protein_quant_current_normalized.csv.gz')))
audit={'source':'Nusinow2020 Gygi CCLE normalized proteome','url':'https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz','local_path':str(path),'bytes':path.stat().st_size,'sha256':sha(path),'plan_sha256':sha(r/'analysis_plan_v14.md')}
(r/'metadata/CCLE_proteome_cached_source_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
p=pd.read_csv(path,low_memory=False);columns=[c for c in p if '_TenPx' in c and not c.endswith('_Peptides')];p=p.set_index('Gene_Symbol')[columns].apply(pd.to_numeric,errors='coerce');a=protein_tidy(p);a.columns=[c.rsplit('_TenPx',1)[0] for c in a.columns];a=a.T.groupby(level=0).mean().T
meta=pd.read_csv(out/'DepMap_model_eligibility_and_ER_crosswalk.csv').set_index('ModelID');meta=meta[meta.exclusion.isna()|meta.exclusion.eq('')];lookup=meta.dropna(subset=['CCLEName']).reset_index().drop_duplicates('CCLEName',keep='first').set_index('CCLEName')
cross=pd.DataFrame({'CCLEName':a.columns});cross['ModelID']=cross.CCLEName.map(lookup.ModelID);cross['mapped']=cross.ModelID.notna();cross.to_csv(out/'CCLE_protein_model_crosswalk.csv',index=False)
good=cross.dropna(subset=['ModelID']);a=a.loc[:,good.CCLEName];a.columns=good.ModelID;md=meta.loc[a.columns].copy();breastids=md.index[md.lineage=='Breast'];b=a.loc[:,breastids];md=md.loc[breastids];fixed=panel(b,md.ER_status,'CCLE_breast');data=[];alltests=[];rawcontrast=None
rna=pd.read_csv(out/'CCLE_fixed_RNA_log2TPM.csv.gz',index_col=0)
for scale,mat in [('source',b),('rank',b.rank(axis=0,pct=True))]:
 sc=scores_for(mat,fixed);z=md[['CellLineName','PatientID','ER_status']].join(sc);z['ESR1_protein']=mat.loc['ESR1'].reindex(z.index) if 'ESR1' in mat.index else np.nan
 rs=[]
 rs.append(dict(test='P1',endpoint='ER status / import protein adjusted other mitochondrial',**conditional_test(z,'Import_TIM_TOM','ER_status',['Nonimport_mito'],binary=True)))
 rs.append(dict(test='P2',endpoint='ESR1 protein / import protein adjusted other mitochondrial',**conditional_test(z,'Import_TIM_TOM','ESR1_protein',['Nonimport_mito'])))
 common=[g for g in fixed['Import_TIM_TOM'] if g in rna.index]
 z['import_RNA']=module_score(rna,common,minimum=max(12,int(np.ceil(.8*len(common))))).reindex(z.index) if len(common)>=12 else np.nan
 z['matched_import_protein']=module_score(mat,common,minimum=max(12,int(np.ceil(.8*len(common))))).reindex(z.index) if len(common)>=12 else np.nan
 zz=z.dropna(subset=['import_RNA','matched_import_protein'])
 if len(zz)>=10:
  rho,pv=stats.spearmanr(zz.import_RNA,zz.matched_import_protein);bs=[]
  for _ in range(5000):
   k=rng.integers(0,len(zz),len(zz));bs.append(stats.spearmanr(zz.import_RNA.to_numpy()[k],zz.matched_import_protein.to_numpy()[k]).statistic)
  lo,hi=np.nanquantile(bs,[.025,.975]);rs.append(dict(test='P3',endpoint='matched import RNA/protein Spearman',n=len(zz),n_genes=len(common),estimate=rho,low=lo,high=hi,p=pv,status='tested'))
 else:rs.append(dict(test='P3',endpoint='matched import RNA/protein Spearman',n=len(zz),status='unavailable',reason='coverage or complete-case n gate',p=np.nan))
 for row in rs:row['scale']=scale
 alltests.append(family(rs));z['scale']=scale;data.append(z.reset_index())
 if scale=='source':
  groups=z.dropna(subset=['ER_status','Import_TIM_TOM']).groupby('ER_status').Import_TIM_TOM.mean();rawcontrast=groups.get(1,np.nan)-groups.get(0,np.nan)
  deps=pd.read_csv(out/'DepMap_pan_cancer_analysis_data.csv',index_col=0);joined=z.join(deps[['Import_TIM_TOM_GE','Nonimport_mito_GE','import_minus_other_GE']],how='inner');joined.to_csv(out/'CCLE_protein_dependency_matched_models.csv')
pd.concat(alltests).to_csv(out/'CCLE_protein_primary_tests.csv',index=False);pd.concat(data).to_csv(out/'CCLE_protein_analysis_data.csv',index=False);b.loc[[g for g in set(fixed['Import_TIM_TOM']+fixed['Nonimport_mito']+['ESR1']) if g in b.index]].to_csv(out/'CCLE_fixed_protein_quantities.csv.gz',compression='gzip')

# Asleh: read the complete protein source sheet, not any selected variable list.
source=raw/'Asleh2022_41467_2022_28524_MOESM9_ESM.xlsx';cache=out/'Asleh_all_gene_protein_quantities.csv.gz'
if cache.exists():a=pd.read_csv(cache,index_col=0);source_rows=9088
else:
 w=openpyxl.load_workbook(source,read_only=True,data_only=True);sheet=w['Protein Data'];source_rows=sheet.max_row-2;it=sheet.iter_rows(min_row=2,values_only=True);header=next(it);headers=[str(x).strip() if x is not None else '' for x in header];idxsym=headers.index('Symbol');quant=[i for i,h in enumerate(headers) if h.startswith('set_')];cols=[headers[i] for i in quant];vals=[];symbols=[]
 for row in it:
  sym=symbol(row[idxsym])
  if sym:
   symbols.append(sym);vals.append([row[i] if isinstance(row[i],(int,float)) else np.nan for i in quant])
 w.close();a=pd.DataFrame(vals,index=symbols,columns=cols).groupby(level=0).median();a.to_csv(cache,compression='gzip')
clinical=pd.read_excel(raw/'Asleh2022_41467_2022_28524_MOESM4_ESM.xlsx',sheet_name='Supplementary Data 1c',header=1);clinical.columns=[str(c).strip() for c in clinical];clinical=clinical.dropna(subset=['Study_ID']).drop_duplicates('Study_ID');clinical['study_id']=clinical.Study_ID.astype(str)
design=pd.read_excel(raw/'Asleh2022_41467_2022_28524_MOESM4_ESM.xlsx',sheet_name='Supplementary Data 1a',header=2);design.columns=[str(c).strip() for c in design];design=design.dropna(subset=['Study_ID','Sample_ set ID']);design['Study_ID']=design.Study_ID.astype(str);design['Sample_ set ID']=design['Sample_ set ID'].astype(str).str.strip()
cross=design[design['Sample_ set ID'].isin(a.columns)&design.Study_ID.isin(clinical.study_id)][['Study_ID','Sample_ set ID']].drop_duplicates();assert not cross['Sample_ set ID'].duplicated().any();cross.to_csv(out/'Asleh_TMT_patient_crosswalk.csv',index=False)
a=a.loc[:,cross['Sample_ set ID']];a.columns=cross.Study_ID;a=a.T.groupby(level=0).mean().T
md=clinical.set_index('study_id').reindex(a.columns);md['luminal']=md['PAM50 Subtype'].map({'Luminal A':1,'Luminal B':1,'Basal-like':0});md['ER_status']=md['ER_status (clinical)'].astype(str).str.lower().map({'positive':1,'negative':0,'pos':1,'neg':0,'+':1,'-':0});md['era']=md.Cohort.astype(str).str.contains('86').astype(int);fixedA=panel(a,md.luminal,'Asleh_luminal_vs_basal');ad=[];at=[];aslehcontrast=None
for scale,mat in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 z=md.join(scores_for(mat,fixedA));rs=[]
 rs.append(dict(test='A1',endpoint='luminal versus basal adjusted import protein',**conditional_test(z,'Import_TIM_TOM','luminal',['Nonimport_mito','era'],binary=True)))
 counts=z.dropna(subset=['Import_TIM_TOM','ER_status','Nonimport_mito','era']).groupby('ER_status').size()
 res=conditional_test(z,'Import_TIM_TOM','ER_status',['Nonimport_mito','era'],binary=True) if len(counts)==2 and counts.min()>=10 else dict(status='unavailable',reason='clinical ER group below10 or coverage',p=np.nan,n=int(counts.sum()))
 rs.append(dict(test='A2',endpoint='clinical ER status adjusted import protein',**res));rs.append(dict(test='A3',endpoint='luminal breast-cancer-specific survival',status='unavailable',reason='individual time/event data controlled access, not publicly deposited',p=np.nan))
 for rec in rs:rec['scale']=scale
 at.append(family(rs));z['scale']=scale;ad.append(z.reset_index(drop=True))
 if scale=='source':
  means=z.dropna(subset=['luminal','Import_TIM_TOM']).groupby('luminal').Import_TIM_TOM.mean();aslehcontrast=means.get(1,np.nan)-means.get(0,np.nan)
pd.concat(at).to_csv(out/'Asleh_protein_primary_tests.csv',index=False);pd.concat(ad).to_csv(out/'Asleh_protein_analysis_data.csv',index=False)
a.loc[[g for g in set(fixedA['Import_TIM_TOM']+fixedA['Nonimport_mito']+['ESR1']) if g in a.index]].to_csv(out/'Asleh_fixed_protein_quantities.csv.gz',compression='gzip');pd.DataFrame(cover).to_csv(out/'Independent_proteome_gene_coverage.csv',index=False)
trigger=bool(np.isfinite(rawcontrast) and np.isfinite(aslehcontrast) and rawcontrast*aslehcontrast>0 and len(fixed['Import_TIM_TOM'])>=12 and len(fixedA['Import_TIM_TOM'])>=12)
decision=dict(CCLE_unadjusted_ER_difference=rawcontrast,Asleh_unadjusted_luminal_difference=aslehcontrast,Johansson_triggered=trigger,rule='both meet >=12 fixed protein coverage and same sign of unadjusted contrast, regardless of significance',Asleh_source_protein_rows=source_rows,Asleh_unique_clinical_patients=len(md))
(out/'Johansson_conditional_gate.json').write_text(json.dumps(decision,indent=2),encoding='utf-8')
print('CCLE');print(pd.concat(alltests)[['scale','test','n','estimate','p','q','status']].to_string(index=False));print('ASLEH');print(pd.concat(at).to_string(index=False));print('CONDITIONAL',json.dumps(decision));print(pd.DataFrame(cover).to_string(index=False))
