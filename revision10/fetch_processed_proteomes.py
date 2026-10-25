from fetch_metadata import fetch,dest
import urllib.request,io,zipfile,json,hashlib,concurrent.futures
rows=[]
cat=json.loads((dest/'PXD049343_files_path.json').read_text())
jobs=[]
for v in cat:
 if v['fileName'].endswith(('.txt','.tsv')):
  url=next(q['value'] for q in v['publicFileLocations'] if q['name']=='FTP Protocol').replace('ftp://','https://',1);jobs.append(('PXD049343_'+v['fileName'],url))
jobs.append(('Palazestrant_OA_record.xml','https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC10911704'))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows.extend(ex.map(fetch,jobs))

class RemoteZip(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.cache={};self.transferred=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  if n==0:return b''
  if n>100_000_000:raise ValueError('Unexpectedly large remote ZIP read')
  start=self.pos;end=min(self.size,start+n);key=(start,end)
  if key not in self.cache:
   req=urllib.request.Request(self.url,headers={'User-Agent':'MitoTemporal public data analysis','Range':f'bytes={start}-{end-1}'})
   with urllib.request.urlopen(req,timeout=120) as f:
    if f.status!=206:raise ValueError('Server ignored byte range; full archive not downloaded')
    b=f.read(end-start+1)
   assert len(b)==end-start,(len(b),end-start);self.transferred+=len(b);self.cache[key]=b
  self.pos=end;return self.cache[key]

for v in json.loads((dest/'PXD002797_files_path.json').read_text()):
 if not v['fileName'].endswith('.zip'):continue
 url=next(q['value'] for q in v['publicFileLocations'] if q['name']=='FTP Protocol').replace('ftp://','https://',1);remote=RemoteZip(url,v['fileSizeBytes']);rec=dict(file=v['fileName'],url=url,archive_bytes=v['fileSizeBytes'],method='HTTP_range_ZIP_CRC_checked_extraction')
 try:
  with zipfile.ZipFile(remote) as z:
   listing=[dict(name=q.filename,bytes=q.file_size,compressed_bytes=q.compress_size,CRC=q.CRC) for q in z.infolist()];(dest/('PXD002797_'+v['fileName']+'_listing.json')).write_text(json.dumps(listing,indent=2))
   members=[q for q in z.infolist() if q.filename.endswith('proteinGroups.txt')]
   assert members,[(q.filename,q.file_size) for q in z.infolist()][:10]
   files=[]
   for member in members:
    b=z.read(member);name='PXD002797_'+v['fileName'].replace('.zip','')+'_'+member.filename.replace('/','_').replace('\\','_');(dest/name).write_bytes(b);files.append(dict(file=name,member=member.filename,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),CRC=member.CRC))
   rec.update(status='extracted',files=files,transferred_bytes=remote.transferred)
 except Exception as e:rec.update(status='failed',error=str(e),transferred_bytes=remote.transferred)
 rows.append(rec)
 print(rec,flush=True)
(dest/'processed_proteomics_manifest.json').write_text(json.dumps(rows,indent=2))
for v in rows:
 if v.get('method') is None:print(v['file'],v['status'],v.get('bytes',v.get('error')))
