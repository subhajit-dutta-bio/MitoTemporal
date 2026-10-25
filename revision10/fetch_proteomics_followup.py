from fetch_metadata import fetch,dest
import urllib.request,concurrent.futures,json,re
jobs=[('PXD002797_project.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD002797'),('PXD002797_files_v3.json','https://www.ebi.ac.uk/pride/ws/archive/v3/files/byProject?accession=PXD002797&pageSize=100'),('PXD049343_files_v3.json','https://www.ebi.ac.uk/pride/ws/archive/v3/files/byProject?accession=PXD049343&pageSize=100'),('PXD002797_files_path.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD002797/files'),('Zhou2024_figshare.html','https://acs.figshare.com/doi/suppl/10.1021/acs.jproteome.4c00102'),('Drabovich2016_xlsx_alt.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/PMC5083091/bin/supp_M115.057257_SuppTables.xlsx'),('Palazestrant2024_methods.docx','https://pmc.ncbi.nlm.nih.gov/articles/PMC10911704/bin/mct-23-0351_supplementary_methods_1_suppsupp.docx')]
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
redirects=[]
for suffix in [2,3,4]:
 url=f'https://doi.org/10.1021/acs.jproteome.4c00102.s{suffix:03}'
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MitoTemporal academic analysis'}),timeout=60) as f:final=f.url;b=f.read()
  name=f'Zhou2024_s{suffix:03}_landing.html';(dest/name).write_bytes(b);redirects.append(dict(doi=url,url=final,file=name));mt=re.search(r'/(\d+)(?:/\d+)?(?:\?|$)',final)
  if mt:rows.append(fetch((f'Zhou2024_figshare_{mt[1]}.json',f'https://api.figshare.com/v2/articles/{mt[1]}')))
 except Exception as e:redirects.append(dict(doi=url,error=str(e)))
(dest/'proteomics_followup_manifest.json').write_text(json.dumps(rows,indent=2));(dest/'Zhou2024_doi_redirects.json').write_text(json.dumps(redirects,indent=2))
for v in rows:print(v['file'],v['status'],v.get('bytes',v.get('error')))
print('DOI redirects',redirects)
