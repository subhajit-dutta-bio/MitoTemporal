from pathlib import Path
import urllib.request,json,hashlib,datetime,concurrent.futures
from lxml import html
r=Path(__file__).resolve().parent;m=r/'metadata';raw=r/'raw';raw.mkdir(exist_ok=True)
plan=hashlib.sha256((r/'analysis_plan_v14.md').read_bytes()).hexdigest()
jobs={
 'GSE99680_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99680/soft/GSE99680_family.soft.gz',
 'hg19ToHg38.over.chain.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
 'Nusinow2020_sample_information.xlsx':'https://pmc.ncbi.nlm.nih.gov/articles/instance/7339254/bin/NIHMS1569324-supplement-TS1.xlsx'}
s=html.fromstring((m/'Asleh2022_article.html').read_text(encoding='utf-8'))
for a in s.xpath('//a[@href]'):
 u=a.get('href')
 if 'MOESM' in u and u.endswith('.xlsx'):jobs['Asleh2022_'+u.split('/')[-1]]=u
e=json.loads((m/'ENCODE_ENCSR059HDE.json').read_text())
for f in e['files']:
 if f.get('output_type')=='loops' and f.get('assembly')=='GRCh38':jobs[f['accession']+'.bedpe.gz']='https://www.encodeproject.org'+f['href']
for col in ['h.all','c5.go.bp','c2.cp.reactome']:jobs[col+'.v2025.1.Hs.symbols.gmt']=f'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/{col}.v2025.1.Hs.symbols.gmt'
def get(item):
 name,u=item;p=raw/name
 try:
  with urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'MitoTemporal public data research'}),timeout=120) as resp:b=resp.read()
  p.write_bytes(b);return dict(file=name,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),plan_sha256=plan,status='downloaded')
 except Exception as e:return dict(file=name,url=u,status='error',error=str(e))
records=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(get,jobs.items()))
(m/'assay_input_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8');print(json.dumps(records,indent=2))
