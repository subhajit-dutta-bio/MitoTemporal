from fetch_inputs import fetch,d
import re,json,concurrent.futures
jobs=[('POETIC_2019_MOESM2.xlsx','https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13058-019-1223-z/MediaObjects/13058_2019_1223_MOESM2_ESM.xlsx'),('GSE78167_RAW.tar','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78167/suppl/GSE78167_RAW.tar')]
for acc in ['GSE89888','GSE78167']:
    jobs.append((acc+'_full_record.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=full'))
for name in re.findall('href="([^"]+csv.gz)"',(d/'GSE89888_suppl.html').read_text()):
    jobs.append((name,'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89888/suppl/'+name))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(d/'source_manifest_processed.json').write_text(json.dumps(rows,indent=2))
for z in rows:print(z['file'],z['status'],z.get('bytes',z.get('error')))
