from pathlib import Path
import urllib.request,hashlib,json,datetime
r=Path(__file__).resolve().parent;m=r/'metadata';m.mkdir(exist_ok=True)
planhash=hashlib.sha256((r/'analysis_plan_NRF1.md').read_bytes()).hexdigest()
records=[]
for acc in ['GSE292988','GSE292986']:
 for folder,name in [('matrix',acc+'_series_matrix.txt.gz'),('suppl',acc+'_raw_counts.txt.gz')]:
  u=f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/{acc}/{folder}/{name}'
  with urllib.request.urlopen(u,timeout=90) as response:b=response.read()
  assert b[:2]==b'\x1f\x8b';(m/name).write_bytes(b)
  records.append(dict(accession=acc,url=u,downloaded_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),plan_sha256=planhash))
  print(acc,folder,len(b))
(m/'NRF1_source_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
