from v14_common import *
from pyliftover import LiftOver
lo=LiftOver(str(r/'raw/hg19ToHg38.over.chain.gz'));rng=np.random.default_rng(20261025);rows=[]
for acc,condition,rep in [('GSM3595410','E2',1),('GSM3595411','Vehicle',1),('GSM3595412','E2',2),('GSM3595413','Vehicle',2)]:
 raw=next((r/'raw').glob(acc+'*.narrowPeak.gz'));a=pd.read_csv(raw,sep='\t',header=None,usecols=[0,1,2]);saved=pd.read_csv(out/f'ER_{condition}_rep{rep}_GRCh38_peaks.csv.gz');lookup=set(map(tuple,saved[['chrom','start','end']].itertuples(index=False,name=None)));sample=a.iloc[rng.choice(len(a),min(500,len(a)),replace=False)];verified=0;ineligible=0
 for chrom,s,e in sample.itertuples(index=False,name=None):
  x=lo.convert_coordinate(chrom,int(s));y=lo.convert_coordinate(chrom,int(e)-1)
  canonical=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
  if len(x)!=1 or len(y)!=1 or x[0][0]!=y[0][0] or x[0][2:]!=y[0][2:] or chrom not in canonical or x[0][0] not in canonical:ineligible+=1;continue
  left,right=sorted([x[0][1],y[0][1]]);right+=1
  if not .8<=(right-left)/(e-s)<=1.2:ineligible+=1;continue
  assert (x[0][0],left,right) in lookup,(acc,chrom,s,e,x,y)
  verified+=1
 rows.append(dict(sample=acc,condition=condition,replicate=rep,random_raw_peaks=len(sample),verified_unique_mappings=verified,not_eligible_by_independent_mapping=ineligible,all_retained_intervals_exact_match=True))
pd.DataFrame(rows).to_csv(r/'independent_pyliftover_checks.csv',index=False);(r/'independent_pyliftover_summary.json').write_text(json.dumps(dict(software='pyliftover 0.4.1',chain_sha256=sha(r/'raw/hg19ToHg38.over.chain.gz'),raw_peaks_checked=sum(v['random_raw_peaks'] for v in rows),all_exact=True,rows=rows),indent=2),encoding='utf-8');print(pd.DataFrame(rows).to_string(index=False))
