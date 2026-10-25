from pathlib import Path
import urllib.request,json,concurrent.futures,gzip
r=Path(__file__).parent
def get(acc):
 p=r/'metadata'/f'{acc}_annotation.soft.txt'
 u=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'
 if not p.exists():p.write_bytes(urllib.request.urlopen(u,timeout=120).read())
 t=p.read_text(encoding='utf-8');assert '!platform_table_begin' in t
 print(acc,len(t),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(get,['GPL13376','GPL6480','GPL570','GPL10371','GPL10558']))
