from pathlib import Path
import urllib.request,gzip,json,re,concurrent.futures,hashlib
r=Path(__file__).parent
accs=['GSE153470','GSE87411','GSE155192','GSE67916','GSE21618','GSE125738','GSE118713','GSE129544']
def fetch(acc):
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
 row={'accession':acc,'source_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}'}
 for kind,url in [('brief',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),('suppl',base+'/suppl/'),('matrix_listing',base+'/matrix/')]:
  p=r/'metadata'/f'{acc}_{kind}.txt'
  try:
   if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=60).read())
   t=p.read_text(encoding='utf-8');row[kind+'_bytes']=len(t)
   if kind=='brief':row['series']=[l for l in t.splitlines() if re.match(r'!Series_(title|summary|overall_design|pubmed_id|platform_id)',l)];row['n_samples']=len(re.findall(r'!Series_sample_id',t))
   else:row[kind+'_files']=re.findall(r'href="([^"]+)"',t)
  except Exception as e:row[kind+'_error']=str(e)
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:result=list(pool.map(fetch,accs))
(r/'metadata/candidate_inventory.json').write_text(json.dumps(result,indent=2))
