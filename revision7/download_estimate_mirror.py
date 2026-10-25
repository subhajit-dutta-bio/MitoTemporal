from pathlib import Path
import urllib.request,hashlib,json,concurrent.futures
out=Path(__file__).parent/'metadata';base='https://raw.githubusercontent.com/r-forge/estimate/abd60e4979c50a88407f7bfb9211687cb443abc6/pkg/estimate/'
names=['inst/extdata/'+n for n in ['common_genes.txt','SI_geneset.gmt','sample_input.gct','sample_estimate.gct','sample_purity.gct']]+['R/estimateScore.R','R/filterCommonGenes.R']
def fetch(n):
 u=base+n;p=out/('estimate_'+n.split('/')[-1])
 with urllib.request.urlopen(u,timeout=120) as f:b=f.read()
 p.write_bytes(b);return dict(file=p.name,url=u,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,names))
(out/'estimate_source_manifest.json').write_text(json.dumps(rows,indent=2));print(rows)
