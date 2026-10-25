from pathlib import Path
import urllib.request,json,hashlib,gzip,re,datetime,concurrent.futures
r=Path(__file__).resolve().parent;m=r/'metadata';raw=r/'raw'
jobs={'DepMap24Q2_full_inventory.json':'https://api.figshare.com/v2/articles/25880521/files?page_size=1000',
 'DepMap24Q2_article.json':'https://api.figshare.com/v2/articles/25880521',
 'Asleh2022_PRIDE_file_inventory.json':'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD024322/files?pageSize=100'}
text=gzip.open(raw/'GSE99680_family.soft.gz','rt').read()
for part in re.split(r'\^SAMPLE = ',text)[1:]:
 acc=part.splitlines()[0]
 if acc in ['GSM3595410','GSM3595411','GSM3595412','GSM3595413']:
  for line in part.splitlines():
   if line.startswith('!Sample_supplementary_file'):
    u=line.split(' = ',1)[1].replace('ftp://','https://');jobs[u.split('/')[-1]]=u
def get(item):
 name,u=item;folder=raw if name.endswith('.gz') else m
 try:
  with urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'MitoTemporal public data research','Accept':'application/json' if name.endswith('.json') else '*/*'}),timeout=90) as resp:b=resp.read()
  (folder/name).write_bytes(b);return dict(file=name,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='downloaded')
 except Exception as e:return dict(file=name,url=u,status='error',error=str(e))
records=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(get,jobs.items()));(m/'remaining_input_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8');print(json.dumps(records,indent=2))
