from pathlib import Path
import json,gzip,re,urllib.request,hashlib,concurrent.futures
r=Path(__file__).resolve().parent;raw=r/'raw';m=r/'metadata'
jobs=[]
for part in re.split(r'\^SAMPLE = ',gzip.open(raw/'GSE99680_family.soft.gz','rt').read())[1:]:
 acc=part.splitlines()[0]
 if acc in ['GSM3595410','GSM3595411','GSM3595412','GSM3595413']:
  for line in part.splitlines():
   if line.startswith('!Sample_supplementary_file'):
    name=line.split('/')[-1];jobs.append((acc,name))
def get(item):
 acc,name=item;u=f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file&file={name}'
 try:
  with urllib.request.urlopen(u,timeout=55) as resp:b=resp.read()
  assert b[:2]==b'\x1f\x8b',(name,b[:80]);(raw/name).write_bytes(b);return dict(accession=acc,file=name,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='downloaded')
 except Exception as e:return dict(accession=acc,file=name,url=u,status='error',error=str(e))
rec=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(get,jobs));(m/'ER_peak_download_retry_manifest.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec,indent=2))
