from fetch_metadata import fetch,dest
import gzip,csv,concurrent.futures,json
jobs=[]
for acc in ['GSE241764','GSE241943']:
 with gzip.open(dest/(acc+'_series_matrix.txt.gz'),'rt') as f:
  for l in f:
   if l.startswith('!Sample_supplementary_file'):
    urls=next(csv.reader([l],delimiter='\t'))[1:]
    for url in urls:
     if url.startswith('ftp://'):jobs.append((url.split('/')[-1],url.replace('ftp://','https://',1)))
jobs.extend([('Drabovich2016_SuppTables.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/instance/5083091/bin/supp_M115.057257_SuppTables.xlsx'),('Drabovich2016_sciencedirect.html','https://www.sciencedirect.com/science/article/pii/S1535947620335593'),('Drabovich2016_PMC.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC5083091/'),('Palazestrant2024.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10911704/fullTextXML')])
for i in [3,4,5]:
 name=f'pr4c00102_si_{i:03}.xlsx';jobs.append((name,'https://pmc.ncbi.nlm.nih.gov/articles/instance/11536429/bin/'+name))
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:rows=list(ex.map(fetch,jobs))
(dest/'input_download_manifest.json').write_text(json.dumps(rows,indent=2))
print('Downloaded',sum(v['status']=='downloaded' for v in rows),'/',len(rows))
for v in rows:
 if v['status']!='downloaded' or not v['file'].startswith('GSM'):print(v['file'],v['status'],v.get('bytes',v.get('error')))
