from pathlib import Path
import urllib.request,concurrent.futures,json,hashlib,time
r=Path(__file__).parent; out=r/'metadata'; jobs=[]
for acc,kind,files in [
 ('GSE158724','suppl',['GSE158724_10x.cell_types.txt.gz','GSE158724_10x.counts.txt.gz','GSE158724_icell8.cell_types.txt.gz','GSE158724_icell8.counts.txt.gz']),
 ('GSE124647','matrix',['GSE124647_series_matrix.txt.gz']),
 ('GSE27463','suppl',['GSE27463_RefSeq.all.tsv.gz']),
 ('GSE93204','matrix',['GSE93204_series_matrix.txt.gz'])]:
 for name in files:jobs.append((name,f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/{kind}/{name}'))
jobs.append(('GSE93204_record.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93204&targ=self&form=text&view=full'))
def fetch(job):
 name,url=job;p=out/name;tmp=p.with_suffix(p.suffix+'.part');t=time.time()
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'MitoTemporal academic secondary analysis'})
   with urllib.request.urlopen(req,timeout=180) as x,tmp.open('wb') as f:
    while True:
     b=x.read(1048576)
     if not b:break
     f.write(b)
   tmp.replace(p)
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1048576),b''):h.update(b)
  row=dict(file=name,url=url,bytes=p.stat().st_size,sha256=h.hexdigest(),status='downloaded',seconds=round(time.time()-t,1))
 except Exception as e:row=dict(file=name,url=url,status='failed',error=str(e))
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(out/'extension_download_manifest.json').write_text(json.dumps(rows,indent=2))
