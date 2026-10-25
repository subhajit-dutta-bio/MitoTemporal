from pathlib import Path
import urllib.request,gzip,json,zipfile,io,xml.etree.ElementTree as E
r=Path(__file__).parent
acc='GSE289825';base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/{acc}'
for kind,url in [('soft',f'{base}/soft/{acc}_family.soft.gz'),('suppl',f'{base}/suppl/')]:
 try:
  b=urllib.request.urlopen(url,timeout=40).read();txt=gzip.decompress(b).decode() if kind=='soft' else b.decode();(r/'metadata'/f'{acc}_{kind}.txt').write_text(txt)
  if kind=='suppl':print(txt)
  else:print('\n'.join(x for x in txt.splitlines() if x.startswith(('!Series_title','!Series_overall_design','!Sample_title','!Sample_characteristics_ch1','!Series_supplementary_file'))))
 except Exception as e:print(type(e).__name__,str(e))
z=json.loads((r/'metadata/zenodo10359740.json').read_text())
for f in z['files']:
 b=urllib.request.urlopen(f['links']['self'],timeout=45).read();(r/'downloads'/'TNBC_2024_public_code.zip').write_bytes(b)
 with zipfile.ZipFile(io.BytesIO(b)) as zz:
  items=[{'name':x.filename,'bytes':x.file_size} for x in zz.infolist() if not x.is_dir()]
  (r/'metadata/TNBC_2024_code_manifest.json').write_text(json.dumps(items,indent=2));print('CODE_FILES',json.dumps(items)[:9000])
