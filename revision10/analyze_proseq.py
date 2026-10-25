from analysis_utils import *
import itertools
md=metadata('GSE241764');md['line']=md.title.str.extract(r'\b(MCF7|CAMA1)\b',expand=False);md['hour']=md.title.str.extract(r'(15\s*min|1hr|6hr|24hr)',expand=False).map(lambda v:.25 if '15' in v else float(v[:-2]));md['group']=md.title.str.replace(r'^\d+ (MCF7|CAMA1) ','',regex=True).str.replace(r'\s*(15\s*mins?|1hr|6hr|24hr)$','',regex=True).str.strip();md['group']=md.group.replace({'24 HOUR ESTRADIOL + OP-1250':'E2_pre24h_OP','Estradiol':'E2','DMSO':'Vehicle','Fulvestrant':'Fulv','OP-1250':'OP1250','Tamoxifen':'Tam'})
assert md.groupby(['line','hour','group']).size().eq(2).all()
parts=[];audit=[];reference=None;refgenes=None
for sample in md.index:
 file=next((r/'metadata').glob(sample+'_*counts.txt.gz'));x=pd.read_csv(file,sep='\t',header=None,names=['chr','start','end','gene','counts']);key=x[['chr','start','end','gene']].astype(str).agg(':'.join,axis=1)
 if reference is None:reference=key.tolist();refgenes=x.gene.map(canon).tolist()
 assert key.is_unique and set(key)==set(reference);counts=pd.to_numeric(x.counts,errors='raise');assert counts.ge(0).all() and counts.notna().all();counts.index=key;parts.append(counts.reindex(reference).rename(sample));audit.append(dict(sample=sample,all_rows=len(x),library_counts=int(counts.sum()),input_file=file.name,row_alignment='exact_coordinate_symbol_key'))
raw=pd.concat(parts,axis=1);cp=raw.divide(raw.sum(),axis=1)*1e6;cp.index=refgenes;a=np.log2(1+cp.groupby(level=0).median());a.to_csv(out/'GSE241764_log2_CPM.csv.gz',compression='gzip');freeze(pd.DataFrame(audit),'GSE241764_library_audit.csv');freeze(md.reset_index(),'GSE241764_sample_crosswalk.csv')
allgroups=dict(groups,ER_programme=controls['ER_programme'],ESR1=['ESR1'],NRF1=['NRF1'],TFAM=['TFAM'])
freeze(pd.DataFrame([dict(module=k,n_fixed=len(gs),n_measured=len(set(gs)&set(a.index)),genes=';'.join(g for g in gs if g in a.index)) for k,gs in allgroups.items()]),'GSE241764_gene_coverage.csv')
tests=[];specific=[];samples=[]
def exact(x,y):
 vals=np.r_[x,y];obs=np.mean(x)-np.mean(y);null=[]
 for ix in itertools.combinations(range(4),2):
  mask=np.zeros(4,bool);mask[list(ix)]=True;null.append(vals[mask].mean()-vals[~mask].mean())
 return np.mean(np.abs(null)>=abs(obs)-1e-12)
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 sc=scores(v,allgroups).join(md);samples.append(sc.reset_index().assign(scale=scale))
 for (line,hour),z in sc.groupby(['line','hour']):
  for treat in ['E2','Fulv','OP1250']:
   for mod in allgroups:
    xx=z.loc[z.group==treat,mod].dropna();yy=z.loc[z.group=='Vehicle',mod].dropna();test=welch(xx,yy);test['conditional_Welch_p']=test.pop('p');test['conditional_permutation_p']=exact(xx,yy) if len(xx)==len(yy)==2 else np.nan;test['p']=np.nan
    tests.append(dict(line=line,hour=hour,contrast=treat+'_minus_Vehicle',module=mod,scale=scale,scope='primary_1h' if hour==1 else 'secondary_other_times',replication_status='two_source_libraries_biological_independence_unverified',**test))
   for comp in comparators:
    d=z.Import_TIM_TOM-z[comp];xx=d[z.group==treat].dropna();yy=d[z.group=='Vehicle'].dropna();test=welch(xx,yy);test['conditional_Welch_p']=test.pop('p');test['conditional_permutation_p']=exact(xx,yy) if len(xx)==len(yy)==2 else np.nan;test['p']=np.nan;specific.append(dict(line=line,hour=hour,contrast=treat+'_minus_Vehicle',comparator=comp,scale=scale,scope='primary_1h' if hour==1 else 'secondary_other_times',**test))
res=family(pd.DataFrame(tests),['scale','module','scope']);sp=family(pd.DataFrame(specific),['scale','scope'])
for df,by in [(res,['scale','module','scope']),(sp,['scale','scope'])]:
 for _,ix in df.groupby(by).groups.items():df.loc[ix,'conditional_Welch_q']=bh(df.loc[ix,'conditional_Welch_p'].fillna(1));df.loc[ix,'conditional_permutation_q']=bh(df.loc[ix,'conditional_permutation_p'].fillna(1))
freeze(res,'PROseq_all_tests.csv');freeze(sp,'PROseq_specificity_tests.csv');freeze(pd.concat(samples,ignore_index=True),'PROseq_sample_scores.csv')
print(res.query('module=="Import_TIM_TOM" and scale=="source" and hour==1')[['line','contrast','estimate','low','high','conditional_Welch_p','conditional_Welch_q','conditional_permutation_p']].to_string(index=False))
