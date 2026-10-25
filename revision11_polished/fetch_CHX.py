from pathlib import Path
import urllib.request,hashlib,json,datetime
r=Path(__file__).resolve().parent;m=r/'metadata';m.mkdir(exist_ok=True)
p=r/'analysis_plan_CHX.md'
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8597/matrix/GSE8597_series_matrix.txt.gz'
with urllib.request.urlopen(u,timeout=90) as response:b=response.read()
assert b[:2]==b'\x1f\x8b';(m/'GSE8597_series_matrix.txt.gz').write_bytes(b)
rec={'accession':'GSE8597','url':u,'downloaded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'plan_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'browser_verified':['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8597','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM213318','https://pmc.ncbi.nlm.nih.gov/articles/PMC2248750/']}
(m/'GSE8597_source_manifest.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print('download bytes',len(b))
