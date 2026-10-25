from pathlib import Path
import urllib.request, hashlib, json, concurrent.futures, datetime
r=Path(__file__).parent
d=r/'metadata'; d.mkdir(exist_ok=True)
jobs=[]
for acc in ['GSE89888','GSE78167']:
    base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
    jobs += [(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_suppl.html',base+'/suppl/'),(acc+'_series_matrix.txt.gz',base+'/matrix/'+acc+'_series_matrix.txt.gz')]
jobs += [('POETIC_2019_article.html','https://link.springer.com/article/10.1186/s13058-019-1223-z'),('POETIC_2019_MOESM2.xlsx','https://static-content.springer-cdn.com/esm/art%3A10.1186%2Fs13058-019-1223-z/MediaObjects/13058_2019_1223_MOESM2_ESM.xlsx')]
def fetch(job):
    name,url=job; p=d/name
    try:
        if not p.exists():
            req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 MitoTemporal public-data secondary research'})
            with urllib.request.urlopen(req,timeout=90) as f: b=f.read()
            p.write_bytes(b)
        b=p.read_bytes()
        return dict(file=name,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='downloaded',retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    except Exception as e:return dict(file=name,url=url,status='failed',error=str(e))
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: rows=list(ex.map(fetch,jobs))
    (d/'source_manifest_initial.json').write_text(json.dumps(rows,indent=2))
    for row in rows:print(row['file'],row['status'],row.get('bytes',row.get('error')))
