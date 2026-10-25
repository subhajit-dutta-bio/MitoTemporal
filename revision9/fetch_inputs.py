from fetch_metadata import fetch,dest
import urllib.request,json,re,concurrent.futures
jobs=[]
for acc in ['GSE152323']:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
 jobs.extend([(acc+'_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'),(acc+'_suppl.html',base+'/suppl/'),(acc+'_series_matrix.txt.gz',base+'/matrix/'+acc+'_series_matrix.txt.gz')])
for acc,name in [('GSE152323','GSE152323_rawcount.txt.gz'),('GSE231397','GSE231397_2022-07-30_Annotated_read_counts.csv.gz'),('GSE246363','GSE246363_RNAseq_mcf7_crispr_knockouts_processed_count.tsv.gz'),('GSE246364','GSE246364_RNAseq_mcf7_siRNA_processed_counts.tsv.gz')]:
 jobs.append((name,f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}'))
# Pin the documented cBioPortal transformation to its last modifying commit.
req=urllib.request.Request('https://api.github.com/repos/cBioPortal/datahub/commits?path=public/brca_cptac_2020&per_page=1',headers={'User-Agent':'MitoTemporal public data analysis'})
with urllib.request.urlopen(req,timeout=60) as f:b=f.read()
(dest/'CPTAC_commit_record.json').write_bytes(b);commit=json.loads(b)[0]['sha']
listing=json.loads((dest/'CPTAC_datahub_listing.json').read_text());names={v['name'] for v in listing}
for name in ['data_clinical_patient.txt','data_clinical_sample.txt','data_mrna_seq_fpkm.txt','data_protein_quantification.txt','meta_mrna_seq_fpkm.txt','meta_protein_quantification.txt','README.md','LICENSE']:
 assert name in names
 domain='https://media.githubusercontent.com/media/' if name.startswith('data_') else 'https://raw.githubusercontent.com/'
 jobs.append(('CPTAC_'+name,domain+f'cBioPortal/datahub/{commit}/public/brca_cptac_2020/{name}'))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(dest/'input_download_manifest.json').write_text(json.dumps(rows,indent=2))
for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
