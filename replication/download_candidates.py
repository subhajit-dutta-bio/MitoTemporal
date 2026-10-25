from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,hashlib,json,tarfile,time
root=Path(__file__).parent
(root/'downloads').mkdir(exist_ok=True)
def download(acc):
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar'
 dest=root/'downloads'/(acc+'_RAW.tar')
 if not dest.exists():
  h=hashlib.sha256();total=0;start=time.time()
  with urllib.request.urlopen(url,timeout=90) as r,open(dest.with_suffix('.partial'),'wb') as w:
   size=int(r.headers.get('Content-Length',0));milestone=0
   while True:
    b=r.read(4*1024*1024)
    if not b:break
    w.write(b);h.update(b);total+=len(b)
    if total//(100*1024*1024)>milestone:
     milestone=total//(100*1024*1024);print(acc,'downloaded_MB',round(total/1e6),flush=True)
  dest.with_suffix('.partial').rename(dest)
  meta={'url':url,'bytes':total,'sha256':h.hexdigest(),'seconds':round(time.time()-start,1)}
  (root/'metadata'/(acc+'_download.json')).write_text(json.dumps(meta,indent=2))
 with tarfile.open(dest) as tar:listing=[{'name':m.name,'bytes':m.size} for m in tar.getmembers() if m.isfile()]
 (root/'metadata'/(acc+'_archive_manifest.json')).write_text(json.dumps(listing,indent=2))
 print(acc,'COMPLETE',dest.stat().st_size,'FILES',len(listing),json.dumps(listing)[:6000],flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(download,['GSE205589','GSE300475']))
