from fetch_inputs import fetch,d
import json,urllib.request,concurrent.futures
queries=[]
for p in d.glob('figshare_exact_*.json'):queries+=json.loads(p.read_text())
jobs=[]
for item in queries:
 if not any(k in item['title'] for k in ['Table S','Table 2','Table 1','File S6']):continue
 url=item['url'];name=f'figshare_{item["id"]}.json'
 row=fetch((name,url))
 if row['status']!='downloaded':continue
 data=json.loads((d/name).read_text())
 for f in data.get('files',[]):jobs.append((f['name'],f['download_url']))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(d/'FAIMOS_supplement_manifest.json').write_text(json.dumps(rows,indent=2))
for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
