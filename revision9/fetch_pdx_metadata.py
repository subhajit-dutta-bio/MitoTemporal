from fetch_metadata import fetch,dest
import json,concurrent.futures
acc='GSE246365';base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
jobs=[(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_series_matrix.txt.gz',base+'/matrix/'+acc+'_series_matrix.txt.gz'),(acc+'_suppl.html',base+'/suppl/')]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(fetch,jobs))
(dest/'PDX_metadata_manifest.json').write_text(json.dumps(rows,indent=2))
for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
