from fetch_inputs import fetch,d
import json,urllib.request,re,concurrent.futures
jobs=[]
text=(d/'PLK1_2020_article.html').read_text()
for url in sorted(set(re.findall(r'https[^" ]+MOESM\d+_ESM.xlsx',text))):jobs.append((url.split('/')[-1],url))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: rows=list(ex.map(fetch,jobs))
for query in ['"Effect of Aromatase Inhibition on Functional"','"Molecular Profiling of Aromatase Inhibitor"']:
 try:
  req=urllib.request.Request('https://api.figshare.com/v2/articles/search',data=json.dumps({'search_for':query,'limit':100}).encode(),headers={'Content-Type':'application/json'})
  with urllib.request.urlopen(req,timeout=60) as f:b=f.read()
  p=d/('figshare_exact_'+str(len(rows))+'.json');p.write_bytes(b);found=json.loads(b)
  print('EXACT',query,[(x['id'],x['title']) for x in found]);rows.append(dict(query=query,file=p.name,status='searched',n=len(found)))
 except Exception as e:print(e)
(d/'FAIMOS_followup_manifest.json').write_text(json.dumps(rows,indent=2))
for z in rows[:len(jobs)]:print(z['file'],z['status'],z.get('bytes',z.get('error')))
