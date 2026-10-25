from pathlib import Path
import urllib.request,concurrent.futures,json,hashlib
r=Path(__file__).parent;out=r/'metadata'; jobs=[
('estimate_1.0.13.tar.gz','https://download.r-forge.r-project.org/src/contrib/estimate_1.0.13.tar.gz'),
('FELINE_Supplementary_Tables.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/instance/8547038/bin/NIHMS1697171-supplement-Supplementary_Tables.xlsx'),
('FELINE_Extended1.csv','https://pmc.ncbi.nlm.nih.gov/articles/instance/8547038/bin/NIHMS1697171-supplement-Numerical_Source_Data_for_Extended_Data_Figure_1_.csv'),
('FELINE_Figure1.csv','https://pmc.ncbi.nlm.nih.gov/articles/instance/8547038/bin/NIHMS1697171-supplement-Numerical_Source_Data_for_Figure_1_.csv'),
('GSE75971_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75971/matrix/GSE75971_series_matrix.txt.gz'),
('GPL6480_soft.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480&targ=self&form=text&view=full'),
('NRF1_MCF7_GRCh38_optimal.bed.gz','https://www.encodeproject.org/files/ENCFF269RME/@@download/ENCFF269RME.bed.gz'),
('gencode.v47.annotation.gtf.gz','https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz'),
('FELINE_article.xml','https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC8547038'),
('GSE75971_record.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75971&targ=self&form=text&view=full'),
('GSE75971_suppl.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75971/suppl/'),
('GSE93204_suppl.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93204/suppl/'),
]
for name in ['common_genes.txt','SI_geneset.gmt','sample_input.txt','sample_input.gct','sample_estimate.gct','sample_purity.gct']:
 jobs.append(('estimate_'+name,'https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/estimate/inst/extdata/'+name+'?root=estimate'))
def fetch(j):
 name,url=j;p=out/name
 try:
  if not p.exists():
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MitoTemporal academic research'}),timeout=120) as f:p.write_bytes(f.read())
  row=dict(file=name,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='downloaded')
 except Exception as e:row=dict(file=name,url=url,status='failed',error=str(e))
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(out/'regulatory_download_manifest.json').write_text(json.dumps(rows,indent=2))
