from pathlib import Path
import urllib.request,json,gzip,re,concurrent.futures
r=Path(__file__).parent
jobs=[('PMC9909892.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9909892/fullTextXML')]
for pmid in ['26033813','25100562']:
 u=f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json'
 j=json.load(urllib.request.urlopen(u,timeout=60));(r/'metadata'/f'{pmid}_reference.json').write_text(json.dumps(j,indent=2))
 x=j['resultList']['result'][0];print(pmid,x.get('pmcid'),x.get('doi'),flush=True)
 if x.get('pmcid'):
  for suf,url in [('.xml','fullTextXML'),('_supplement.zip','supplementaryFiles')]:jobs.append((x['pmcid']+suf,'https://www.ebi.ac.uk/europepmc/webservices/rest/'+x['pmcid']+'/'+url))
for acc in ['GSE55374','GSE59515']:
 jobs.append((acc+'_soft.txt',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz'))
def get(job):
 n,u=job;p=r/'metadata'/n
 try:
  if not p.exists():
   b=urllib.request.urlopen(u,timeout=90).read();p.write_bytes(gzip.decompress(b) if n.endswith('_soft.txt') else b)
  return {'file':n,'url':u,'bytes':p.stat().st_size}
 except Exception as e:return {'file':n,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(get,jobs))
(r/'metadata/fetch_log.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
