from fetch_metadata import fetch,dest
import json
name='GSE246365_RNAseq_pdx_processed_counts.tsv.gz'
z=fetch((name,'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246365/suppl/'+name))
(dest/'PDX_input_manifest.json').write_text(json.dumps([z],indent=2));print(z['file'],z['status'],z.get('bytes',z.get('error')))
