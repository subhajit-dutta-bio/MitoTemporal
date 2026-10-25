from fetch_inputs import fetch,d
import json,urllib.request,concurrent.futures
jobs=[('FAIMOS2014_article.html','https://aacrjournals.org/clincancerres/article/20/9/2485/78894/Effect-of-Aromatase-Inhibition-on-Functional-Gene'),('FAIMOS2013_article.html','https://aacrjournals.org/clincancerres/article/19/10/2775/77640/Molecular-Profiling-of-Aromatase-Inhibitor-Treated'),('PLK1_2020_article.html','https://www.nature.com/articles/s41467-020-17697-1')]
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(fetch,jobs))
for query in ['CCR-13-2602','Molecular Profiling of Aromatase Inhibitor','Effect of Aromatase Inhibition on Functional']:
 try:
  req=urllib.request.Request('https://api.figshare.com/v2/articles/search',data=json.dumps({'search_for':query,'limit':100}).encode(),headers={'Content-Type':'application/json'})
  with urllib.request.urlopen(req,timeout=60) as f:b=f.read()
  p=d/('figshare_search_'+str(len(rows))+'.json');p.write_bytes(b)
  z=json.loads(b);print(query,[(a['id'],a['title']) for a in z][:15])
  rows.append(dict(query=query,file=p.name,status='searched',n=len(z)))
 except Exception as e:rows.append(dict(query=query,status='failed',error=str(e)))
(d/'FAIMOS_search_manifest.json').write_text(json.dumps(rows,indent=2))
for z in rows[:3]:print(z['file'],z['status'],z.get('bytes',z.get('error')))
