from pathlib import Path
import sys,json,re,io
r=Path(__file__).parent
sys.path.insert(0,str(r.parent/'python_packages'))
import pandas as pd,numpy as np
gene_file=r.parent/'source/hummel2/gene_sets/mito_modules.json'
if not gene_file.exists():gene_file=r/'discovery/mito_modules.json'
mods=json.loads(gene_file.read_text())
mods['CIV_nuclear']=[g for g in mods['OXPHOS_Complex_IV'] if not g.startswith('MT-')]
out=r/'analysis';(out/'gene_expression').mkdir(exist_ok=True)
coverage=[];diagnostics=[]

def symbol(row,acc):
 if acc=='GSE146661':
  ss={x.split(' // ')[1].strip() for x in str(row.get('gene_assignment','')).split(' /// ') if len(x.split(' // '))>1}
 elif acc=='GSE147271':ss={str(row.get('ORF','')).strip()}
 else:
  name={'GSE111563':'Symbol','GSE20181':'Gene Symbol','GSE138233':'GENE_SYMBOL'}[acc]
  ss=set(re.split(r'\s*///\s*|\s*//\s*|\s*;\s*',str(row.get(name,''))))
 ss.discard('');ss.discard('nan');ss.discard('---')
 return next(iter(ss)).upper() if len(ss)==1 else None

def soft(acc):
 f=open(r/'metadata'/f'{acc}_soft.txt',encoding='utf-8');platform=[];samples={};meta=[];current=None;md={}
 for line in f:
  if line.startswith('!platform_table_begin'):
   lines=[]
   for l in f:
    if l.startswith('!platform_table_end'):break
    lines.append(l)
   platform=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str,keep_default_na=False,low_memory=False)
  elif line.startswith('^SAMPLE'):
   if current:meta.append(md)
   current=line.strip().split(' = ')[1];md={'sample':current}
  elif current and line.startswith('!Sample_'):
   k,_,v=line.strip().partition(' = ')
   if k in ['!Sample_title','!Sample_source_name_ch1','!Sample_characteristics_ch1','!Sample_data_processing']:
    key=k[8:];md[key]=md.get(key,'')+(' | ' if key in md else '')+v
  elif line.startswith('!sample_table_begin'):
   lines=[]
   for l in f:
    if l.startswith('!sample_table_end'):break
    lines.append(l)
   df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=['ID_REF','VALUE'],dtype={'ID_REF':str})
   samples[current]=df.set_index('ID_REF').VALUE
 if current:meta.append(md)
 platform['gene']=[symbol(row,acc) for row in platform.to_dict('records')]
 platform[['ID','gene']].to_csv(out/f'{acc}_probe_gene_map.csv',index=False)
 pd.DataFrame(meta).to_csv(out/f'{acc}_sample_metadata.csv',index=False)
 a=pd.DataFrame(samples)
 if len(a):
  a['gene']=platform.set_index('ID').gene.reindex(a.index)
  a=a.dropna(subset=['gene']).groupby('gene').median()
 return a,pd.DataFrame(meta),platform

def score(acc,a,transform):
 a=a.apply(pd.to_numeric,errors='coerce');a=a.groupby(a.index.str.upper()).median()
 diagnostics.append({'accession':acc,'transform':transform,'n_genes':len(a),'n_samples':a.shape[1],'before_min':float(a.min().min()),'before_median':float(np.nanmedian(a.to_numpy())),'before_max':float(a.max().max())})
 if transform=='log2_linear':
  assert (a>=0).all().all();a=np.log2(a.clip(lower=.01))
 elif transform=='log2_RPKM_plus1':a=np.log2(a+1)
 elif transform=='log2_CPM_plus1':a=np.log2(a.div(a.sum(axis=0),axis=1)*1e6+1)
 # Missing any sample means gene is not eligible; no imputation.
 a=a.dropna(axis=0,how='any')
 ranks=a.rank(axis=0,pct=True)
 centered=a.sub(a.mean(axis=1),axis=0)
 res=pd.DataFrame(index=a.columns);rankres=res.copy()
 used=set()
 for name,genes in mods.items():
  gs=[g.upper() for g in genes if g.upper() in a.index]
  eligible=len(gs)>=max(3,np.ceil(.5*len(genes)))
  coverage.append({'accession':acc,'module':name,'n_original':len(genes),'n_mapped':len(gs),'eligible':bool(eligible),'mapped_genes':';'.join(gs),'missing_genes':';'.join(g for g in genes if g.upper() not in a.index)})
  if eligible:res[name]=centered.loc[gs].mean();rankres[name]=ranks.loc[gs].mean()
  used.update(gs)
 for g in ['CRAT','AR','TCF7','TOX','PDCD1','PPARGC1A']:
  if g in a.index:used.add(g)
 a.loc[sorted(used)].to_csv(out/'gene_expression'/f'{acc}_module_genes.csv.gz',compression='gzip')
 res.index.name='sample';rankres.index.name='sample'
 res.to_csv(out/f'{acc}_module_scores.csv');rankres.to_csv(out/f'{acc}_rank_scores.csv')
 print(acc,'scored',a.shape,'CIV coverage',[v['n_mapped'] for v in coverage if v['accession']==acc and v['module']=='OXPHOS_Complex_IV'],flush=True)

for acc,transform in [('GSE111563','already_log2'),('GSE20181','log2_linear'),('GSE147271','already_log2'),('GSE146661','already_log2')]:
 a,md,p=soft(acc);score(acc,a,transform)

# Parse clone platform metadata now; expression is supplied separately as a workbook.
a,md,p=soft('GSE138233')
book=pd.ExcelFile(r/'downloads/GSE138233_Normalized_Matrix.xlsx')
print('CLONE_WORKBOOK',book.sheet_names,flush=True)
raw=pd.read_excel(book,sheet_name=0)
print('CLONE_COLUMNS',raw.columns.tolist(),flush=True)
raw.to_pickle(out/'clone_matrix_staging.pkl')

for acc,file,transform in [('GSE138459','GSE138459_rpkm.txt.gz','log2_RPKM_plus1'),('GSE296966','GSE296966_counts.txt.gz','log2_CPM_plus1')]:
 a=pd.read_csv(r/'downloads'/file,sep='\t',index_col=0);score(acc,a,transform)
pd.DataFrame(coverage).to_csv(out/'gene_coverage.csv',index=False)
(out/'processing_diagnostics.json').write_text(json.dumps(diagnostics,indent=2))
