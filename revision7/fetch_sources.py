from pathlib import Path
import sys,concurrent.futures,json,hashlib
r=Path(__file__).parent;sys.path.insert(0,str(r.parent/'python_packages'))
import urllib.request
out=r/'metadata';out.mkdir(exist_ok=True)
jobs=[]
for acc in ['GSE158724','GSE284956','GSE124647','GSE27463']:
 jobs.extend([(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_suppl.html',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/'),(acc+'_matrix.html',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/matrix/')])
for pmc in ['PMC8547038','PMC4888666','PMC5555232','PMC2262171','PMC4895922','PMC3058701']:
 jobs.append((pmc+'.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'))
jobs.append(('ENCSR135ANT.json','https://www.encodeproject.org/experiments/ENCSR135ANT/?format=json'))
def fetch(job):
 name,url=job;p=out/name
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'MitoTemporal secondary research source audit','Accept':'application/json' if name.endswith('.json') else '*/*'})
  with urllib.request.urlopen(req,timeout=90) as x:content=x.read()
  p.write_bytes(content)
  return {'file':name,'url':url,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'status':'downloaded'}
 except Exception as e:return {'file':name,'url':url,'status':'failed','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(out/'source_download_manifest.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
