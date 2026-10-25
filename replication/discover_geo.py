from pathlib import Path
import urllib.request,gzip,json,re
root=Path(__file__).parent
(root/'metadata').mkdir(parents=True,exist_ok=True)
for acc in ['GSE300475','GSE205589','GSE291055']:
 prefix=acc[:-3]+'nnn';base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}'
 urls={'soft':f'{base}/soft/{acc}_family.soft.gz','supplementary_listing':f'{base}/suppl/'}
 for kind,url in urls.items():
  target=root/'metadata'/f'{acc}_{kind}.txt'
  try:
   b=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MitoTemporal-academic-reanalysis/1.0'}),timeout=60).read()
   text=gzip.decompress(b).decode() if kind=='soft' else b.decode()
   target.write_text(text,encoding='utf-8')
   print(acc,kind,'bytes',len(b),flush=True)
   if kind=='supplementary_listing':print(text[:7000],flush=True)
  except Exception as e:print(acc,kind,type(e).__name__,str(e),flush=True)
