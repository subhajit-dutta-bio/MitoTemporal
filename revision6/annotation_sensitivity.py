from pathlib import Path
import ast,sys,json,csv,io,gzip,hashlib
r=Path(__file__).parent;root=r.parent;out=r/'analysis';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'))
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,mannwhitneyu

hpath=r/'metadata/hgnc_complete_set.txt'
h=pd.read_csv(hpath,sep='\t',dtype=str,keep_default_na=False);h=h[h.status=='Approved']
approved=set(h.symbol.str.upper());prev={};alias={}
for row in h.itertuples():
 for name,target in [('prev_symbol',prev),('alias_symbol',alias)]:
  for token in getattr(row,name).split('|'):
   if token:target.setdefault(token.upper(),set()).add(row.symbol.upper())
def canon(g):
 g=str(g).upper()
 if g in approved:return g,'approved'
 for label,table in [('previous',prev),('alias',alias)]:
  if g in table:
   return (next(iter(table[g])),label) if len(table[g])==1 else (g,'ambiguous_'+label)
 return g,'unresolved'

oldmods=json.loads((root/'source/hummel2/gene_sets/mito_modules.json').read_text());mods={};mr=[]
for name,gs in oldmods.items():
 current=[]
 for g in gs:
  new,status=canon(g);include=not(name=='OXPHOS_Complex_I' and g=='NDUFA4')
  duplicate=new in current
  mr.append({'module':name,'old_symbol':g,'canonical_symbol':new,'mapping_status':status,'retained':include and not duplicate,'reason':'Complex IV assignment' if not include else ('duplicate gene' if duplicate else 'retained')})
  if include and not duplicate:current.append(new)
 mods[name]=current
assert canon('NDUFA4')[0] not in mods['OXPHOS_Complex_I']
assert canon('NDUFA4')[0] in mods['OXPHOS_Complex_IV']
assert mods['mtUPR'].count('DDIT3')==1
pd.DataFrame(mr).to_csv(out/'canonical_module_gene_audit.csv',index=False)
(out/'canonical_modules.json').write_text(json.dumps(mods,indent=2))

def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
maps=[];coverage=[];rows=[];deltas=[]
def canonicalize(a,acc):
 mapping={g:canon(g) for g in a.index}
 maps.extend({'accession':acc,'source_symbol':g,'canonical_symbol':v[0],'status':v[1]} for g,v in mapping.items() if v[1]!='approved')
 a=a.copy();a.index=[mapping[g][0] for g in a.index];return a.groupby(level=0).median()
def score(a,acc):
 a=canonicalize(a,acc).dropna(axis=0,how='any');rr=a.rank(axis=0,pct=True);ss={}
 for scale,v in [('source',a),('rank',rr)]:
  s=pd.DataFrame(index=a.columns)
  for name,gs in mods.items():
   matched=[g for g in gs if g in a.index];eligible=len(matched)>=max(3,np.ceil(len(gs)/2))
   if scale=='source':coverage.append({'accession':acc,'module':name,'n_original_unique':len(gs),'n_measured':len(matched),'eligible':bool(eligible),'genes':';'.join(matched)})
   if eligible:s[name]=v.loc[matched].mean()
  ss[scale]=s
 return ss

# Import only pure parser functions from the preserved source, never execute its analyses.
src=(root/'revision2/prepare_scores.py').read_text();nodes=[v for v in ast.parse(src).body if isinstance(v,ast.FunctionDef) and v.name in ['symbol','soft']]
ns={'pd':pd,'np':np,'io':io,'re':__import__('re'),'r':root/'revision2','out':out}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(root/'revision2/prepare_scores.py'),'exec'),ns)
manifest=pd.read_csv(root/'revision2/analysis/additional_clinical_sample_manifest.csv',dtype=str)
for acc,contrasts in [('GSE20181',[('A','B'),('A','C')]),('GSE111563',[('1','4')]),('GSE147271',[('Pretreatment','Tamoxifen')]),('GSE153470',[('baseline','on')])]:
 if acc=='GSE153470':
  a=pd.read_pickle(root/'revision4/analysis/GSE153470_gene_expression.pkl');pair=pd.read_csv(root/'revision4/analysis/GSE153470_pair_crosswalk.csv',dtype=str).set_index('patient')
 else:
  a,_,_=ns['soft'](acc)
  if acc=='GSE20181':a=np.log2(a.clip(lower=.01))
  md=manifest[manifest.accession==acc];pair=md.pivot(index='patient',columns='timepoint',values='sample')
 ss=score(a,acc)
 for scale,s in ss.items():
  for pre,post in contrasts:
   pp=pair[[pre,post]].dropna();dd=s.loc[pp[post]].set_axis(pp.index)-s.loc[pp[pre]].set_axis(pp.index)
   for name in dd:
    d=dd[name].dropna();p=wilcoxon(d,alternative='two-sided').pvalue if (d!=0).any() else 1.
    rows.append({'accession':acc,'contrast':post+'_minus_'+pre,'module':name,'scale':scale,'family':'uncontrolled_clinical_all_contrasts','n':len(d),'mean_delta':d.mean(),'median_delta':d.median(),'n_positive':int((d>0).sum()),'p':p})
    deltas.extend({'accession':acc,'contrast':post+'_minus_'+pre,'scale':scale,'module':name,'patient':pt,'delta':v} for pt,v in d.items())
 print('canonical scores',acc,flush=True)

acc='GSE105777';a=pd.read_pickle(root/'revision4/analysis/GSE105777_gene_expression.pkl');ss=score(a,acc);pp=pd.read_csv(root/'revision4/analysis/POETIC_pair_crosswalk.csv',dtype=str).set_index('patient')
for scale,s in ss.items():
 dd=s.loc[pp.on].set_axis(pp.index)-s.loc[pp.baseline].set_axis(pp.index)
 for name in dd:
  ai=dd.loc[pp.arm=='AI',name];co=dd.loc[pp.arm=='Control',name]
  rows.append({'accession':acc,'contrast':'AI_minus_control_patient_changes','module':name,'scale':scale,'family':'controlled_POETIC','n':len(ai),'n_control':len(co),'mean_delta':ai.mean()-co.mean(),'mean_AI':ai.mean(),'mean_control':co.mean(),'p':mannwhitneyu(ai,co,alternative='two-sided').pvalue})
  deltas.extend({'accession':acc,'contrast':pp.loc[pt,'arm']+'_on_minus_baseline','scale':scale,'module':name,'patient':pt,'delta':v} for pt,v in dd[name].items())

acc='GSE87411';lines=[];inside=False
with gzip.open(root/'revision4/metadata/GSE87411_series_matrix.txt.gz','rt',encoding='utf-8') as f:
 for l in f:
  if l.startswith('!series_matrix_table_begin'):inside=True
  elif l.startswith('!series_matrix_table_end'):break
  elif inside:lines.append(l)
a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0);a.index=a.index.astype(str);mp=pd.read_csv(root/'revision4/analysis/GSE87411_probe_to_gene.csv',index_col=0,dtype=str).iloc[:,0];a['gene']=mp.reindex(a.index);a=a.dropna(subset=['gene']).groupby('gene').median();a=canonicalize(a,acc)
pp=pd.read_csv(root/'revision4/analysis/GSE87411_pair_crosswalk.csv',dtype=str).set_index('patient');pcover=[]
for scale,v in [('source',a),('rank',a.rank(axis=0,pct=True))]:
 for name,gs in mods.items():
  matched=[g for g in gs if g in a.index];ds=[]
  for pt,z in pp.iterrows():
   x=v.loc[matched,[z.baseline,z.on]].dropna();ok=len(x)>=max(3,np.ceil(len(gs)/2))
   if scale=='source':pcover.append({'patient':pt,'module':name,'n_original_unique':len(gs),'n_measured':len(x),'eligible':bool(ok),'genes':';'.join(x.index)})
   if ok:
    d=(x[z.on]-x[z.baseline]).mean();ds.append(d);deltas.append({'accession':acc,'contrast':'on_minus_baseline','scale':scale,'module':name,'patient':pt,'delta':d})
  if len(ds)>=3:
   d=np.asarray(ds);rows.append({'accession':acc,'contrast':'paired_observed_on_minus_baseline','module':name,'scale':scale,'family':'Z1031_technical','n':len(d),'mean_delta':d.mean(),'median_delta':np.median(d),'n_positive':int((d>0).sum()),'p':wilcoxon(d).pvalue if (d!=0).any() else 1.})
pd.DataFrame(pcover).to_csv(out/'canonical_Z1031_pair_gene_coverage.csv',index=False)
t=pd.DataFrame(rows)
for _,ix in t.groupby(['family','scale']).groups.items():t.loc[ix,'q']=bh(t.loc[ix,'p']);t.loc[ix,'family_size']=len(ix)
t.to_csv(out/'canonical_annotation_all_tests.csv',index=False);pd.DataFrame(deltas).to_csv(out/'canonical_annotation_patient_deltas.csv',index=False);pd.DataFrame(coverage).to_csv(out/'canonical_annotation_gene_coverage.csv',index=False);pd.DataFrame(maps).drop_duplicates().to_csv(out/'canonical_expression_symbol_mapping.csv',index=False)
summary=t[t.module.isin(['Import_TIM_TOM','Mitophagy'])];summary.to_csv(out/'canonical_key_results.csv',index=False)
source={'url':'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt','retrieved':'2026-10-24','sha256':hashlib.sha256(hpath.read_bytes()).hexdigest(),'bytes':hpath.stat().st_size,'approved_records':len(h),'protocol':'annotation_sensitivity_protocol.txt'}
(r/'metadata/HGNC_source_inventory.json').write_text(json.dumps(source,indent=2));print(summary.to_string(index=False));print('ALL_TESTS',len(t))
