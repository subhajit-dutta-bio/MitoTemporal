from pathlib import Path
import urllib.request,json,concurrent.futures,hashlib,gzip
r=Path(__file__).parent; jobs=[]
for a in ['GSE153470','GSE87411','GSE155192','GSE67916','GSE21618','GSE125738','GSE118713','GSE129544']:
 jobs.append((a+'_series_matrix.txt.gz',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}/matrix/{a}_series_matrix.txt.gz'))
for a,n in [('GSE155192','TPM_ChronicStimHypoxia.csv.gz'),('GSE155192','TPM_Ndufs4KO.csv.gz'),('GSE155192','TPM_PGC1aOE.csv.gz'),('GSE125738','fpkm_T47D-TR_vs_T47D.txt.gz'),('GSE118713','TRANSCRIPT_TAMR_FASR_MCF7_expression_data.tsv.gz'),('GSE129544','T47DR_vs_T47D_mRNA.txt.gz')]:
 jobs.append((a+'_'+n,f'https://ftp.ncbi.nlm.nih.gov/geo/series/{a[:-3]}nnn/{a}/suppl/{a}_{n}'))
for a in ['GPL13376','GPL6480','GPL570','GPL10371']:
 jobs.append((a+'_annotation.soft.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}&targ=self&form=text&view=full'))
def fetch(job):
 n,u=job;p=r/'metadata'/n
 try:
  if not p.exists():
   with urllib.request.urlopen(u,timeout=120) as f,p.open('wb') as out:
    while True:
     b=f.read(1024*1024)
     if not b:break
     out.write(b)
  data=p.read_bytes()
  if n.endswith('.gz'):gzip.decompress(data)
  else:assert b'!platform_table_begin' in data
  result={'file':n,'url':u,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'format_valid':True}
 except Exception as e:
  result={'file':n,'url':u,'error':str(e)}
 print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:result=list(ex.map(fetch,jobs))
(r/'metadata/download_log.json').write_text(json.dumps(result,indent=2))
