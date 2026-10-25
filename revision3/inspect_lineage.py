from pathlib import Path
import re,json,sys,xml.etree.ElementTree as E
import pandas as pd
r=Path(__file__).parent;old=r.parent/'revision2'
new=set(pd.read_csv(old/'analysis/additional_clinical_sample_manifest.csv',dtype=str).query('accession == "GSE111563"').patient)
rec=[]
for acc in ['GSE55374','GSE59515']:
 t=(r/'metadata'/f'{acc}_soft.txt').read_text(encoding='utf-8')
 titles=re.findall(r'!Sample_title = (.*)',t);pts=set(re.findall(r'patient(\d+[A-Z]?)',' '.join(titles)))
 print(acc,'samples',len(titles),'patients',len(pts),'ids',sorted(pts),'numeric shared with111563',sorted({re.sub('[A-Z]','',p) for p in pts}&new))
 print([l for l in t.splitlines() if l.startswith('!Series_overall_design')])
 for title in titles:
  m=re.search(r'patient(\d+[A-Z]?)',title);rec.append({'accession':acc,'title':title,'patient_label':m.group(1) if m else None})
pd.DataFrame(rec).to_csv(r/'analysis/related_letrozole_sample_labels.csv',index=False)
t=E.parse(r/'metadata/PMC9909892.xml')
for x in t.iter('p'):
 s=''.join(x.itertext())
 if 'Twenty-eight' in s or 'twenty-eight' in s or 'GSE147271' in s:print('RADIO_PRIMARY',s)
meta=t.find('.//article-meta')
print('RADIO_BIB',E.tostring(meta,encoding='unicode')[:5500])
