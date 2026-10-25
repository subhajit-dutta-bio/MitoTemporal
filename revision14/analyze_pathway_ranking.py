from v14_common import *
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)

raw=r/'raw';inputs=[];members={};aliases={}
for file in ['h.all.v2025.1.Hs.symbols.gmt','c5.go.bp.v2025.1.Hs.symbols.gmt','c2.cp.reactome.v2025.1.Hs.symbols.gmt']:
 p=raw/file;inputs.append(dict(path=str(p.relative_to(root)),sha256=sha(p),bytes=p.stat().st_size))
 for line in p.read_text(encoding='utf-8').splitlines():
  fields=line.split('\t');genes=tuple(sorted({canon(g) for g in fields[2:] if canon(g) in approved}))
  if not 15<=len(genes)<=500:continue
  aliases.setdefault(genes,[]).append(fields[0])
for i,(genes,names) in enumerate(sorted(aliases.items(),key=lambda x:sorted(x[1])[0])):
 members['MSIG_%05d'%i]=genes
alias=pd.DataFrame([dict(module=k,display_name=sorted(aliases[g])[0],aliases=';'.join(sorted(aliases[g])),n_fixed=len(g),import_overlap=len(set(g)&set(mods['Import_TIM_TOM'])),genes=';'.join(g)) for k,g in members.items()])
members['Import_TIM_TOM']=tuple(mods['Import_TIM_TOM'])
alias=pd.concat([alias,pd.DataFrame([dict(module='Import_TIM_TOM',display_name='Canonical TIM/TOM import',aliases='fixed canonical panel',n_fixed=24,import_overlap=24,genes=';'.join(mods['Import_TIM_TOM']))])],ignore_index=True).set_index('module')
alias.to_csv(out/'Pathway_fixed_membership_and_aliases.csv.gz',compression='gzip')
allresults=[];summary=[];scores=[];crosswalk=[]

def load(path,kind='csv'):
 p=root/path;inputs.append(dict(path=path,sha256=sha(p),bytes=p.stat().st_size))
 a=pd.read_pickle(p) if kind=='pkl' else pd.read_csv(p,index_col=0)
 return tidy_matrix(a)

def incidence(index):
 lookup={g:i for i,g in enumerate(index)};rr=[];cc=[];names=list(members)
 for i,name in enumerate(names):
  for g in members[name]:
   if g in lookup:rr.append(i);cc.append(lookup[g])
 return names,csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(names),len(index)))

def score(a,names,M):
 values=a.to_numpy(float);valid=np.isfinite(values);counts=M@valid.astype(float);sums=M@np.where(valid,values,0)
 minimum=np.array([12 if n=='Import_TIM_TOM' else max(10,int(np.ceil(len(members[n])/2))) for n in names])[:,None]
 return np.where(counts>=minimum,sums/np.maximum(counts,1),np.nan),counts

def paired(a,cross,context,pre='baseline',post='on',patient='patient',direction=-1):
 names,M=incidence(a.index);ds=[];dr=[];counts=[]
 for row in cross.to_dict('records'):
  b=row[pre];o=row[post];x=a[[b,o]].copy();joint=x.notna().all(axis=1);x.loc[~joint,:]=np.nan
  s,n=score(x,names,M);rk,_=score(x.rank(axis=0,pct=True),names,M)
  ds.append(s[:,1]-s[:,0]);dr.append(rk[:,1]-rk[:,0]);counts.append(n[:,0]);crosswalk.append(dict(context=context,unit=str(row[patient]),baseline=b,on=o,arm=row.get('arm','AI')))
 return names,np.array(ds).T,np.array(dr).T,np.array(counts).T

def finish(context,names,source,rank,counts,direction,group=None,primary=True):
 minimum_n=int(np.ceil(.8*source.shape[1]));eligible=np.isfinite(rank).sum(axis=1)>=minimum_n
 if group is not None:
  eligible=np.ones(len(names),bool)
  for g in [0,1]:eligible&=np.isfinite(rank[:,group==g]).sum(axis=1)>=int(np.ceil(.8*(group==g).sum()))
 result=alias.loc[names].reset_index();result['context']=context;result['primary_context']=primary;result['direction']=direction
 result['original_units']=source.shape[1];result['n_observed']=np.isfinite(rank).sum(axis=1);result['minimum_genes_observed']=np.nanmin(np.where(np.isfinite(rank),counts,np.nan),axis=1);result['eligible']=eligible
 for scale,values in [('source',source),('rank',rank)]:
  if group is None:
   effect=np.nanmean(values,axis=1);sd=np.nanstd(values,axis=1,ddof=1);n=np.isfinite(values).sum(axis=1);se=sd/np.sqrt(n);df=n-1;stat=effect/se;p=2*stats.t.sf(abs(stat),df)
  else:
   b=values[:,group==0];t=values[:,group==1];n0=np.isfinite(b).sum(axis=1);n1=np.isfinite(t).sum(axis=1);v0=np.nanvar(b,axis=1,ddof=1)/n0;v1=np.nanvar(t,axis=1,ddof=1)/n1
   effect=np.nanmean(t,axis=1)-np.nanmean(b,axis=1);se=np.sqrt(v0+v1);df=(v0+v1)**2/(v0*v0/(n0-1)+v1*v1/(n1-1));stat=effect/se;p=2*stats.t.sf(abs(stat),df);result['n_control']=n0;result['n_exposed']=n1
  crit=stats.t.ppf(.975,df);result[scale+'_effect']=effect;result[scale+'_low']=effect-crit*se;result[scale+'_high']=effect+crit*se;result[scale+'_t']=stat;result[scale+'_p']=np.where(eligible,p,np.nan);result[scale+'_q']=bh(result[scale+'_p'])
 result.to_csv(out/(context+'_all_module_ranking.csv.gz'),index=False,compression='gzip');allresults.append(result)
 target=result[result.module.eq('Import_TIM_TOM')].iloc[0];background=result[result.eligible&~result.module.eq('Import_TIM_TOM')]
 for benchmark in ['all','size15to40','exclude_half_import_overlap']:
  z=background if benchmark=='all' else background[background.n_fixed.between(15,40)] if benchmark=='size15to40' else background[background.import_overlap<12]
  row=target.to_dict();row['benchmark']=benchmark;row['N_comparison_modules']=len(z)
  for method,column in [('rank_effect','rank_effect'),('rank_t','rank_t'),('source_effect','source_effect')]:
   v=direction*target[column];bs=direction*z[column];position=1+(bs>v).sum()+.5*(bs==v).sum();row[method+'_position']=position;row[method+'_top_percent']=100*position/(len(z)+1)
  summary.append(row)
 for j in range(rank.shape[1]):scores.append(dict(context=context,unit=j,source_import=source[names.index('Import_TIM_TOM'),j],rank_import=rank[names.index('Import_TIM_TOM'),j],n_import_genes=counts[names.index('Import_TIM_TOM'),j],group=None if group is None else int(group[j])))
 print(context,'Nsets',len(background),'units',source.shape[1],'import_rank_tail',round(summary[-3]['rank_effect_top_percent'],2),'effect',target['source_effect'],flush=True)

# E2: each trajectory has its own baseline and fixed late-time mean.
a=load('revision8/analysis/GSE78167_TPM.csv.gz');a=np.log2(a+1);meta=pd.read_csv(root/'revision8/analysis/GSE78167_sample_crosswalk.csv');names,M=incidence(a.index);ss=[];rs=[];cs=[]
for rep,z in meta.groupby('replicate'):
 z=z.set_index('hour');cols=[z.loc[h,'sample'] for h in [0,6,8,12,24]];x=a[cols].copy();x.loc[~x.notna().all(axis=1),:]=np.nan;s,n=score(x,names,M);rk,_=score(x.rank(axis=0,pct=True),names,M);ss.append(s[:,1:].mean(axis=1)-s[:,0]);rs.append(rk[:,1:].mean(axis=1)-rk[:,0]);cs.append(n[:,0]);crosswalk.append(dict(context='E2_GSE78167',unit=str(rep),baseline=cols[0],on=';'.join(cols[1:]),arm='E2 fixed late'))
finish('E2_GSE78167',names,np.array(ss).T,np.array(rs).T,np.array(cs).T,1)

# Four ligand-free mutant retention contexts.
a=load('revision8/analysis/GSE89888_log2_coverage_CPM.csv.gz');meta=pd.read_csv(root/'revision8/analysis/GSE89888_sample_crosswalk.csv');names,M=incidence(a.index)
for line in ['MCF7','T47D']:
 for genotype in ['Y537S','D538G']:
  z=meta[meta.cell_line.eq(line)&meta.treatment.eq('veh')&meta.genotype.isin(['WT',genotype])].copy();cols=z.external_id.tolist();group=z.genotype.ne('WT').to_numpy(int);s,n=score(a[cols],names,M);rk,_=score(a[cols].rank(axis=0,pct=True),names,M);context='Mutant_'+line+'_'+genotype
  for row in z.to_dict('records'):crosswalk.append(dict(context=context,unit=row['external_id'],baseline='',on=row['external_id'],arm=row['genotype']))
  finish(context,names,s,rk,n,1,group)

a=load('revision9/analysis/GSE231397_log2_CPM.csv.gz');meta=pd.read_csv(root/'revision9/analysis/ER_loss_sample_crosswalk.csv');z=meta[meta.accession.eq('GSE231397')&meta.group.isin(['E2','Fulv+E2'])];cols=z['index'].tolist();group=z.group.eq('Fulv+E2').to_numpy(int);names,M=incidence(a.index);s,n=score(a[cols],names,M);rk,_=score(a[cols].rank(axis=0,pct=True),names,M)
finish('Fulvestrant_GSE231397',names,s,rk,n,-1,group)
for row in z.to_dict('records'):crosswalk.append(dict(context='Fulvestrant_GSE231397',unit=row['index'],baseline='',on=row['index'],arm=row['group']))

for acc,pre,post,patient in [('GSE20181','A','C','patient'),('GSE111563','1','4','patient'),('GSE153470','baseline','on','patient')]:
 a=load('revision7/analysis/'+acc+'_full_gene_expression.csv.gz');cross=pd.read_csv(root/('revision7/analysis/'+acc+'_mechanism_pair_crosswalk.csv'));cross=cross.dropna(subset=[pre,post]);names,s,rk,n=paired(a,cross,acc,pre,post,patient);finish(acc,names,s,rk,n,-1)
a=load('revision8/analysis/GSE87411_full_observed_expression.csv.gz');cross=pd.read_csv(root/'revision4/analysis/GSE87411_pair_crosswalk.csv');names,s,rk,n=paired(a,cross,'GSE87411');finish('GSE87411',names,s,rk,n,-1)
a=load('revision7/analysis/GSE93204_full_gene_expression.pkl','pkl');cross=pd.read_csv(root/'revision7/analysis/NeoPalAna_pair_crosswalk.csv').dropna(subset=['Baseline','C1D1']);names,s,rk,n=paired(a,cross,'GSE93204','Baseline','C1D1','patient id');finish('GSE93204',names,s,rk,n,-1)
a=load('revision8/analysis/GSE105777_full_observed_expression.csv.gz');cross=pd.read_csv(root/'revision8/analysis/Ki67_patient_crosswalk.csv');cross=cross[cross.accession.eq('GSE105777')].copy();names,s,rk,n=paired(a,cross,'POETIC_all',patient='patient_id');treated=cross.arm.eq('AI').to_numpy();finish('GSE105777',names,s[:,treated],rk[:,treated],n[:,treated],-1);finish('POETIC_controlled',names,s,rk,n,-1,treated.astype(int),False)

pd.DataFrame(summary).to_csv(out/'Pathway_import_ranking_summary.csv',index=False);pd.DataFrame(scores).to_csv(out/'Pathway_import_unit_scores.csv',index=False);pd.DataFrame(crosswalk).to_csv(out/'Pathway_sample_unit_crosswalk.csv',index=False);pd.concat(allresults,ignore_index=True).to_csv(out/'Pathway_all_contexts_complete_results.csv.gz',index=False,compression='gzip')
(r/'metadata/Pathway_input_hashes.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
print('DONE',len(members)-1,'distinct backgrounds; 12 primary contexts + controlled POETIC')
