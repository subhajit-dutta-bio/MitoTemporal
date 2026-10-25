from pathlib import Path
import sys,gzip,tarfile
r=Path(__file__).parent;sys.path.insert(0,str(r.parent/'python_packages'))
import pandas as pd
for acc in ['GSE89888','GSE78167']:
 print('\nMETADATA',acc)
 with gzip.open(r/'metadata'/f'{acc}_series_matrix.txt.gz','rt') as f:
  for line in f:
   if line.startswith(('!Sample_title','!Sample_characteristics','!Sample_data_processing','!Sample_supplementary')):print(line[:1600].strip())
print('\nESR1 DE HEAD')
print(pd.read_csv(r/'metadata/GSE89888_resWT_E2vVeh_MCF7.csv.gz').head(2).to_string(index=False))
with tarfile.open(r/'metadata/GSE78167_RAW.tar') as t:
 fs=[v for v in t.getmembers() if v.isfile()];print('E2 FILES',len(fs),[v.name for v in fs[:4]])
 f=t.extractfile(fs[0]);b=f.read();b=gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b;print('E2 firstfile',b[:2000].decode())
x=pd.ExcelFile(r/'metadata/POETIC_2019_MOESM2.xlsx');print('POETIC SHEETS',x.sheet_names)
print(pd.read_excel(x,sheet_name='TableS4',header=None).head(8).to_string(index=False,header=False))
for acc in ['GSE105777','GSE87411']:
 f=r.parent/'revision4/analysis'/f'{acc}_sample_metadata.csv';a=pd.read_csv(f);print(acc,a.shape);print(a[['sample','title']].head(3).to_string(index=False))
