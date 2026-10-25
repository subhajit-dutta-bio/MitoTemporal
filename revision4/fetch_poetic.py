from pathlib import Path
import urllib.request,json,re,concurrent.futures,gzip,hashlib
r=Path(__file__).parent
def get(acc):
 b=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}';res={'accession':acc}
 for name,url in [(acc+'_brief.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_suppl.txt',b+'/suppl/'),(acc+'_series_matrix.txt.gz',b+'/matrix/'+acc+'_series_matrix.txt.gz')]:
  p=r/'metadata'/name
  if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=100).read())
  if name.endswith('brief.txt'):res['series']=[l for l in p.read_text(encoding='utf-8').splitlines() if re.match(r'!Series_(title|summary|overall_design|platform_id)',l)]
  if name.endswith('suppl.txt'):res['files']=re.findall(r'href="([^"]+)"',p.read_text())
  res[name]={'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 print(json.dumps(res),flush=True);return res
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:res=list(ex.map(get,['GSE105777','GSE126870']))
(r/'metadata/POETIC_download_log.json').write_text(json.dumps(res,indent=2))
