from pathlib import Path
import urllib.request,hashlib,json,concurrent.futures
r=Path(__file__).parent
urls={
'GSE138459_rpkm.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138459/suppl/GSE138459_ac_chr_stim_CD8T_NAC_norm_rpkm.txt.gz',
'GSE296966_counts.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296966/suppl/GSE296966_gene_counts.txt.gz',
'PMC6323855.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6323855/fullTextXML',
'PMC7442749.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7442749/fullTextXML',
'PMC12679916.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12679916/fullTextXML'}
def fetch(item):
 name,url=item;p=r/('metadata' if name.endswith('xml') else 'downloads')/name
 if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=90).read())
 b=p.read_bytes();return {'name':name,'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(fetch,urls.items()))
(r/'metadata/selected_downloads.json').write_text(json.dumps(res,indent=2))
print(json.dumps(res,indent=2))
