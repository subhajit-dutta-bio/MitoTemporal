from fetch_inputs import fetch,d
import json,concurrent.futures
base='https://duffel.rail.bio/recount3/human/'
jobs=[('SRP093386_gene_sums.G026.gz',base+'data_sources/sra/gene_sums/86/SRP093386/sra.gene_sums.SRP093386.G026.gz'),('SRP093386_sra.MD.gz',base+'data_sources/sra/metadata/86/SRP093386/sra.sra.SRP093386.MD.gz'),('SRP093386_recount_qc.MD.gz',base+'data_sources/sra/metadata/86/SRP093386/sra.recount_qc.SRP093386.MD.gz'),('recount_G026.gtf.gz',base+'annotations/gene_sums/human.gene_sums.G026.gtf.gz')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: rows=list(ex.map(fetch,jobs))
(d/'source_manifest_recount.json').write_text(json.dumps(rows,indent=2))
for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
