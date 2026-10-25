from pathlib import Path
import urllib.request,json,hashlib,datetime,concurrent.futures
r=Path(__file__).parent;dest=r/'metadata';dest.mkdir(parents=True,exist_ok=True)
def fetch(job):
 name,url=job;p=dest/name
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'MitoTemporal academic public-data analysis/1.0'})
  with urllib.request.urlopen(req,timeout=90) as f:b=f.read()
  p.write_bytes(b)
  return dict(file=name,url=url,status='downloaded',bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 except Exception as e:return dict(file=name,url=url,status='failed',error=str(e))
if __name__=='__main__':
 jobs=[]
 for acc in ['GSE241764','GSE241943']:
  base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
  jobs.extend([(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_series_matrix.txt.gz',base+'/matrix/'+acc+'_series_matrix.txt.gz'),(acc+'_suppl.html',base+'/suppl/')])
 jobs.extend([('Drabovich2016.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5083091/fullTextXML'),('Zhou2024.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11536429/fullTextXML'),('PXD049343_files.json','https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=PXD049343&pageSize=100'),('PXD049343_project.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD049343')])
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
 (dest/'metadata_download_manifest.json').write_text(json.dumps(rows,indent=2))
 for v in rows:print(v['file'],v['status'],v.get('bytes',v.get('error')))
