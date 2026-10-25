from pathlib import Path
import urllib.request,hashlib,json,concurrent.futures,datetime
r=Path(__file__).resolve().parent;dest=r/'metadata';dest.mkdir(exist_ok=True)
def fetch(job):
 name,url=job;p=dest/name
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 MitoTemporal public data analysis'})
   with urllib.request.urlopen(req,timeout=90) as f:b=f.read()
   p.write_bytes(b)
  b=p.read_bytes();return dict(file=name,url=url,status='downloaded',bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 except Exception as e:return dict(file=name,url=url,status='failed',error=str(e))
jobs=[]
for acc in ['GSE231397','GSE246363','GSE246364']:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
 jobs.extend([(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_suppl.html',base+'/suppl/'),(acc+'_series_matrix.txt.gz',base+'/matrix/'+acc+'_series_matrix.txt.gz')])
jobs += [('RECQ1_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8088126/fullTextXML'),('CPTAC_datahub_listing.json','https://api.github.com/repos/cBioPortal/datahub/contents/public/brca_cptac_2020'),('CPTAC_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8077737/fullTextXML')]
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
 (dest/'metadata_download_manifest.json').write_text(json.dumps(rows,indent=2))
 for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
