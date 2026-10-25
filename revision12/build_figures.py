from pathlib import Path
import json,re,shutil,subprocess,io
import sys
root=Path(__file__).resolve().parents[1];r=root/'revision12'
sys.path.insert(0,str(root/'python_packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pypdf import PdfReader,PdfWriter
from pypdf.generic import ContentStream,NameObject
import xml.etree.ElementTree as ET
from scipy.stats import spearmanr
from reportlab.pdfgen.canvas import Canvas

dest=r/'deliverables';figdir=dest/'figures';src=dest/'source_data'
figdir.mkdir(parents=True,exist_ok=True);src.mkdir(exist_ok=True)
out=r/'analysis';out.mkdir(exist_ok=True)
old=root/'revision10/deliverables';manifest=[];provenance=[]
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';RED='#AF5464';INK='#20343F';GREY='#77868E';palette=[BLUE,TEAL,ORANGE,RED]
rng=np.random.default_rng(20261025)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def data(name):return pd.read_csv(old/'source_data'/(name+'.csv.gz'))
def canvas(n,title,sub):
 rows=3 if n in ['Figure_4','Supplementary_Figure_S17'] else 4
 f,axs=plt.subplots(rows,2,figsize=(7.2,8.6 if rows==3 else 9.3));f.subplots_adjust(left=.15,right=.96,bottom=.07,top=.95,hspace=1.03,wspace=.68)
 return f,axs.ravel()
def panel(ax,n,letter,title,z,origin='recomputed_from_existing_tables'):
 ax.set_title(title,loc='left',pad=13,weight='bold',color=INK);pos=ax.get_position();ax.figure.text(pos.x0-.048,pos.y1+.017,letter,fontsize=11,weight='bold')
 file=f'{n}_{letter}.csv.gz';z.to_csv(src/file,index=False,compression='gzip');manifest.append(dict(figure=n,panel=letter,title=title,source_data=file));provenance.append(dict(figure=n,panel=letter,revision10_source=origin))
def note(ax,t):ax.text(0,1.015,t,transform=ax.transAxes,fontsize=5.9,color=GREY,va='bottom')
def zero(ax):ax.axhline(0,c='#CAD3D8',lw=.65,zorder=0)
def save(f,n):
 for ext in ['pdf','svg','png']:f.savefig(figdir/f'{n}.{ext}',dpi=350,bbox_inches='tight',pad_inches=.13)
 plt.close(f)
def forest(ax,z,labels,x='estimate',ref=0):
 y=np.arange(len(z));ax.errorbar(z[x],y,xerr=np.array([z[x]-z.low,z.high-z[x]]),fmt='o',c=BLUE,ms=3,capsize=2,lw=.85);ax.axvline(ref,c=GREY,lw=.65);ax.set_yticks(y,labels);ax.invert_yaxis()
def rain(ax,groups,ylabel):
 for i,(label,values) in enumerate(groups.items()):
  vals=np.asarray(values,float);vals=vals[np.isfinite(vals)];col=palette[i%4]
  if len(vals)>2 and vals.std()>0:
   yy=np.linspace(vals.min(),vals.max(),120);dd=gaussian_kde(vals)(yy);dd=.28*dd/dd.max();ax.fill_betweenx(yy,i-.08-dd,i-.08,color=col,alpha=.22,lw=0)
  ax.scatter(i+.04+rng.uniform(0,.13,len(vals)),vals,s=7,c=col,alpha=.55,edgecolors='none');ax.plot([i-.08,i+.24],[vals.mean()]*2,c=INK,lw=1.5)
 ax.set_xticks(range(len(groups)),list(groups));ax.set_xlim(-.45,len(groups)-.5+.35);ax.set_ylabel(ylabel);zero(ax)
def slope(ax,z):
 for row in z.itertuples():ax.plot([0,1],[row.pre,row.post],c=BLUE,alpha=.2,lw=.65)
 ax.plot([0,1],[z.pre.mean(),z.post.mean()],'-o',c=INK,ms=3,lw=1.5);ax.set_xticks([0,1],['Baseline','On treatment']);ax.set_xlim(-.2,1.2);ax.set_ylabel('Mean source-scale expression')
def dots(ax,z,groups):
 for i,g in enumerate(groups):
  v=z.loc[z.group==g,'Import_TIM_TOM'];ax.scatter(i+rng.uniform(-.065,.065,len(v)),v,s=16,c=palette[i%4]);ax.plot([i-.15,i+.15],[v.mean()]*2,c=INK,lw=1)
 ax.set_xticks(range(len(groups)),groups,rotation=15);ax.set_ylabel('Mean log2(CPM + 1)')
def cultures(ax,z):
 for j,(treatment,col) in enumerate([('veh',BLUE),('E2',ORANGE)]):
  means=[]
  for i,gen in enumerate(['WT','Y537S','D538G']):
   vals=z.loc[(z.genotype==gen)&(z.treatment==treatment),'Import_TIM_TOM'];ax.scatter(i+(j-.5)*.25+rng.uniform(-.035,.035,len(vals)),vals,c=col,s=12);means.append(vals.mean())
  ax.plot(np.arange(3)+(j-.5)*.25,means,'-_',c=col,lw=1,ms=8,label='Vehicle' if treatment=='veh' else 'E2')
 ax.set_xticks(range(3),['WT','Y537S','D538G']);ax.set_ylabel('Mean log2(coverage CPM + 1)');ax.legend(frameon=False,fontsize=6)

n='Figure_1';f,axs=canvas(n,'Aromatase inhibition withdraws import transcription','Untreated-biopsy comparison | paired patients | six series, 462 paired observations')
for ax,letter,oldname,scale,q in zip(axs[:2],'ab',['Figure_1_e','Figure_1_f'],['source','rank'],[.0410,.132]):
 z=data(oldname);panel(ax,n,letter,'POETIC | '+scale+' scale',z,oldname);rain(ax,{'AI\n(n=157)':z[z.contrast.str.startswith('AI_')].delta,'Untreated\n(n=56)':z[z.contrast.str.startswith('Control_')].delta},'Paired import-score change');note(ax,f'Controlled comparison; BH q={q:.3f}')
for ax,letter,oldname,title in zip(axs[2:6],'cdef',['Figure_1_a','Figure_1_b','Figure_1_c','Figure_1_d'],['Late letrozole | n=56','Extended letrozole | n=37','Anastrozole | n=81','NeoPalAna AI run-in | n=22']):
 z=data(oldname);panel(ax,n,letter,title,z,oldname);slope(ax,z)
z=data('Figure_1_g');panel(axs[6],n,'g','Within-sample rank changes',z,'Figure_1_g');forest(axs[6],z,[v.cohort+' ('+str(v.n)+')' for v in z.itertuples()],x='mean_delta');axs[6].set_xlabel('Mean change; 95% patient-bootstrap CI');axs[6].tick_params(axis='y',labelsize=5.4)
z=data('Figure_1_h');panel(axs[7],n,'h','Letrozole | complete patient trajectories',z,'Figure_1_h')
for _,v in z.groupby('patient'):axs[7].plot(v.time,v.score,c=BLUE,alpha=.14,lw=.6)
v=z.groupby('time').score.mean();axs[7].plot(v.index,v,'o-',c=INK,ms=3,lw=1.4);axs[7].set_xticks([0,1,2],['Baseline','10-14 d','90 d']);axs[7].set_ylabel('Import change from own baseline');zero(axs[7]);note(axs[7],'56 complete deposited triplets');save(f,n)

e=pd.read_csv(root/'revision8/analysis/GSE78167_module_trajectories.csv');es=pd.read_csv(root/'revision8/analysis/GSE89888_sample_scores.csv');est=pd.read_csv(root/'revision8/analysis/GSE89888_primary_tests.csv');sp=pd.read_csv(root/'revision10/analysis/specificity_primary_tests.csv')
n='Figure_2';f,axs=canvas(n,'Import transcription follows ER activation and blockade','Replicated cultures | activating ESR1 mutants | reversal | exploratory xenografts')
z=e[(e.scale=='source')&(e.module=='Import_TIM_TOM')];panel(axs[0],n,'a','Estradiol | complete import trajectories',z,'Figure_3_a')
for _,v in z.groupby('replicate'):axs[0].plot(v.hour,v.delta,c=BLUE,alpha=.3,lw=.7)
v=z.groupby('hour').delta.mean();axs[0].plot(v.index,v,'-o',c=BLUE,ms=2,lw=1.4);axs[0].set_xticks([0,2,6,12,24]);axs[0].set_xlabel('Hours after E2');axs[0].set_ylabel('Mean log2(TPM + 1) change');zero(axs[0]);note(axs[0],'Three biological trajectories; 23/24 import genes')
z=pd.read_csv(root/'revision8/analysis/GSE78167_replicate_window_changes.csv');z=z[(z.scale=='source')&(z.module=='Import_TIM_TOM')&(z.window=='late')];panel(axs[1],n,'b','Fixed 6-24 h contrast | each replicate',z);axs[1].scatter(np.arange(3),z.delta,c=BLUE,s=25);axs[1].axhline(z.delta.mean(),c=INK,lw=1);axs[1].set_xticks(range(3),[f'Rep {v}' for v in z.replicate]);axs[1].set_ylabel('Late minus baseline score');zero(axs[1]);note(axs[1],'Mean +0.117; prespecified P=0.00531')
for ax,letter,line in zip(axs[2:4],'cd',['MCF7','T47D']):
 z=es[(es.scale=='source')&(es.cell_line==line)];panel(ax,n,letter,line+' | activating ESR1 mutants',z,'Figure_3_'+('e' if line=='MCF7' else 'f'));cultures(ax,z);note(ax,'Four cultures/group; source clones pooled')
z=data('Figure_4_c');panel(axs[4],n,'e','Fulvestrant | reversal of E2 induction',z,'Figure_4_c');dots(axs[4],z,['Veh','E2','Fulv','Fulv+E2']);note(axs[4],'n=3/group; reversal source q=0.023')
z=data('Figure_4_h');panel(axs[5],n,'f','Giredestrant | all PDX dose contrasts',z,'Figure_4_h');forest(axs[5],z,[v.model+' '+str(v.dose_mg_kg)+' mg/kg' for v in z.itertuples()]);axs[5].set_xlabel('Import: treated minus vehicle');note(axs[5],'Exploratory; four animals/group; shared controls')
z=sp[(sp.scale=='source')&(sp.comparator=='Nonimport_mito')&(~sp.context.isin(['GSE20181','GSE111563','GSE105777','GSE152323']))];panel(axs[6],n,'g','Relative response | other mitochondrial genes',z);labs=[]
for c in z.context:
 labs.append('E2 trajectories' if c.startswith('GSE78167') else 'Fulv reversal' if c.startswith('GSE231397') else c.replace('GSE89888_','').replace('GSE246365_','').replace('_',' '))
forest(axs[6],z,labs);axs[6].tick_params(axis='y',labelsize=6.3);axs[6].set_xlabel('Import minus other-mito response');note(axs[6],'Original 48-test family per scale retained')
z=data('Figure_3_g');panel(axs[7],n,'h','Ligand-independent retention | effect estimates',z,'Figure_3_g');forest(axs[7],z,[v.cell_line+' '+v.contrast for v in z.itertuples()]);axs[7].set_xlabel('Mutant minus WT under vehicle');note(axs[7],'One T47D Y537S estimate remains uncertain');save(f,n)

n='Figure_3';f,axs=canvas(n,'Temporal and regulatory context of the import response','Observed kinetics | nascent transcription | NRF1 binding | patient coupling | growth control')
z=e[(e.scale=='source')&e.module.isin(['Import_TIM_TOM','ER_programme','Proliferation'])].copy();mx=z.groupby(['module','hour']).delta.mean().groupby('module').max();z['display_scaled_delta']=z.delta/z.module.map(mx);panel(axs[0],n,'a','Complete trajectories | timing display',z)
for module,col,label in zip(['Import_TIM_TOM','ER_programme','Proliferation'],[BLUE,TEAL,ORANGE],['Import','ER targets','Proliferation']):
 v=z[z.module==module]
 for _,sub in v.groupby('replicate'):axs[0].plot(sub.hour,sub.display_scaled_delta,c=col,lw=.5,alpha=.18)
 mean=v.groupby('hour').display_scaled_delta.mean();axs[0].plot(mean.index,mean,'-o',c=col,ms=2,lw=1.2,label=label)
axs[0].set_xticks([0,2,4,8,12,24]);axs[0].set_xlabel('Hours after E2');axs[0].set_ylabel('Fraction of observed maximum');axs[0].legend(frameon=False,fontsize=5.6);zero(axs[0]);note(axs[0],'Descriptive scaling; not response-magnitude comparison')
mods=json.loads((root/'revision6/analysis/canonical_modules.json').read_text())['Import_TIM_TOM'];a=pd.read_csv(root/'revision8/analysis/GSE78167_TPM.csv.gz',index_col=0);cross=pd.read_csv(root/'revision8/analysis/GSE78167_sample_crosswalk.csv');a=np.log2(a+1);rows=[]
for row in cross.itertuples():
 baseline=cross[(cross.replicate==row.replicate)&(cross.hour==0)].iloc[0]['sample']
 for gene in mods:rows.append(dict(gene=gene,replicate=row.replicate,hour=row.hour,delta=a.at[gene,row.sample]-a.at[gene,baseline] if gene in a.index else np.nan))
z=pd.DataFrame(rows);panel(axs[1],n,'b','Fixed import genes | temporal coherence',z);hh=z.pivot_table(index='gene',columns='hour',values='delta',dropna=False).reindex(mods);im=axs[1].imshow(hh,cmap='RdBu_r',aspect='auto',vmin=-.6,vmax=.6);axs[1].set_yticks(range(len(hh)),hh.index,fontsize=4.2);axs[1].set_xticks(range(len(hh.columns)),hh.columns,fontsize=4.6);axs[1].set_xlabel('Hours after E2');plt.colorbar(im,ax=axs[1],fraction=.04,pad=.03).set_label('Mean log2(TPM + 1) change',size=5.5);note(axs[1],'All fixed members; missing gene blank')
pro=pd.read_csv(root/'revision10/analysis/PROseq_all_tests.csv');z=pro[(pro.scale=='source')&(pro.module=='Import_TIM_TOM')&(pro.contrast=='E2_minus_Vehicle')];panel(axs[2],n,'c','Nascent import | all fixed E2 times',z,'Figure_8_a')
for line,col in [('MCF7',BLUE),('CAMA1',TEAL)]:
 v=z[z.line==line].sort_values('hour');axs[2].plot(range(4),v.estimate,'-o',c=col,ms=3,lw=1.1,label=line)
axs[2].set_xticks(range(4),['15 min','1 h','6 h','24 h']);axs[2].set_ylabel('E2 minus vehicle log2(CPM + 1)');axs[2].legend(frameon=False,fontsize=6);zero(axs[2]);note(axs[2],'Two source libraries/condition; descriptive')
z=pro[(pro.scale=='source')&(pro.hour==1)&(pro.contrast=='E2_minus_Vehicle')&pro.module.isin(['Import_TIM_TOM','ER_programme','Proliferation'])];panel(axs[3],n,'d','One-hour ER and import response',z,'Figure_8_d')
for i,line in enumerate(['MCF7','CAMA1']):
 v=z[z.line==line].set_index('module').loc[['Import_TIM_TOM','ER_programme','Proliferation']];axs[3].scatter(np.arange(3)+(i-.5)*.18,v.estimate,c=palette[i],s=22,label=line)
axs[3].set_xticks(range(3),['Import','ER targets','Proliferation'],rotation=15);axs[3].set_ylabel('E2 minus vehicle source score');zero(axs[3]);note(axs[3],'Culture independence not verified; no biological P')
z=data('Supplementary_Figure_S8_a');panel(axs[4],n,'e','NRF1 promoter binding | expressed genes',z,'Supplementary_Figure_S8_a');axs[4].bar(range(3),[z.K.iloc[0]/z.N.iloc[0],z.n_bound.iloc[0]/24,z.n_bound.iloc[1]/23],color=[GREY,BLUE,ORANGE],width=.55);axs[4].set_xticks(range(3),['Background','Import\n14/24','Nuclear IV\n12/23']);axs[4].set_ylim(0,.8);axs[4].set_ylabel('Fraction with nearby NRF1 peak');note(axs[4],'Import expression-matched q=0.0428')
corr=pd.read_csv(root/'revision8/analysis/Letrozole_control_correlations.csv')
for ax,letter,acc,title in zip(axs[5:7],'fg',['GSE20181','GSE111563'],['Late letrozole | patient-change coupling','Extended letrozole | patient-change coupling']):
 z=corr[(corr.scale=='source')&(corr.accession==acc)];panel(ax,n,letter,title,z,'Supplementary_Figure_S10_h');ax.barh(z.control,z.rho,color=[BLUE if q<.05 else '#B8C6CF' for q in z.q],height=.55);ax.axvline(0,c=GREY,lw=.65);ax.invert_yaxis();ax.set_xlim(-.2,.8);ax.set_xlabel('Spearman rho | paired import changes');ax.tick_params(axis='y',labelsize=5.5);note(ax,'All six controls; twelve-test family per scale')
nrtests=pd.read_csv(out/'NRF1_import_tests.csv');nrsamples=pd.read_csv(out/'NRF1_sample_scores.csv')
z=nrtests[(nrtests.accession=='GSE292988')&(nrtests.scale=='source')&(nrtests.family=='primary')];panel(axs[7],n,'h','NRF1 loss and rescue | HEK293',z,'revision12_NRF1_primary_tests');forest(axs[7],z,['NRF1 KD - control','FL rescue interaction','Deleted - FL rescue']);axs[7].set_xlabel('Import RNA contrast; 95% HC3 interval');note(axs[7],'Three samples/group; KD q<0.001; rescue unresolved');save(f,n)

n='Figure_4';f,axs=canvas(n,'Protein and clinical boundaries of the transcriptional finding','CPTAC | prespecified Ki67 models | exploratory metastatic outcomes')
for j,oldname,title,x,y in [(0,'Figure_5_a','CPTAC | matched RNA and protein','RNA','protein'),(1,'Figure_5_f','CPTAC ER+ | shared mitochondrial abundance','Nonimport_mito_protein','Import_protein'),(2,'Figure_5_c','Within ER+ | adjusted ER-programme RNA','x_partial','y_partial'),(3,'Figure_5_d','Within ER+ | adjusted ESR1 protein','x_partial','y_partial')]:
 z=data(oldname);panel(axs[j],n,'abcd'[j],title,z,oldname);axs[j].scatter(z[x],z[y],s=10,c=BLUE,alpha=.6,edgecolors='none');axs[j].set_xlabel({'RNA':'Import RNA (same 21 genes)','Nonimport_mito_protein':'Other mitochondrial protein','x_partial':'Residual predictor (SD units)'}[x]);axs[j].set_ylabel('Residual import protein' if j>1 else 'Import protein score')
 if j>1:
  beta=-.11893496798088923 if j==2 else -.095;xx=np.array([z[x].min(),z[x].max()]);axs[j].plot(xx,beta*xx,c=RED,lw=1);zero(axs[j]);note(axs[j],'n=81; adjusted q='+('0.111' if j==2 else '0.119'))
 elif j==0:note(axs[j],'n=122; rho=0.236; q=0.0264')
 else:note(axs[j],'n=81; cross-sectional abundance context')
baseline=data('Figure_6_a').assign(metric='baseline_model_coefficient')
k=pd.read_csv(root/'revision8/analysis/Ki67_patient_scores.csv');kd=pd.read_csv(root/'revision8/analysis/Ki67_delta_correlations.csv');corrrows=[];brng=np.random.default_rng(20261025)
for acc in ['GSE87411','GSE105777']:
 z=k[(k.scale=='source')&(k.accession==acc)&(k.arm=='AI')].dropna(subset=['ki67_pre','ki67_on','Import_TIM_TOM_pre','ER_programme_pre','Proliferation_pre'])
 x=z.Import_TIM_TOM_delta.to_numpy();y=np.log2((z.ki67_on.to_numpy()+.1)/(z.ki67_pre.to_numpy()+.1));rho=spearmanr(x,y).statistic;priorrow=kd[(kd.scale=='source')&(kd.accession==acc)].iloc[0];assert np.isclose(rho,priorrow.rho)
 vals=[]
 for _ in range(5000):
  idx=brng.integers(0,len(x),len(x));vals.append(spearmanr(x[idx],y[idx]).statistic)
 lo,hi=np.quantile(vals,[.025,.975]);corrrows.append(dict(accession=acc,scale='source',metric='delta_Spearman_rho',estimate=rho,low=lo,high=hi,n=len(x),q=priorrow.q,p=priorrow.p,interval_method='5000_patient_percentile_bootstrap'))
correlations=pd.DataFrame(corrrows);correlations.to_csv(out/'Ki67_display_bootstrap_intervals.csv',index=False);z=pd.concat([baseline,correlations],ignore_index=True);panel(axs[4],n,'e','Ki67 response | distinct clinical estimands',z)
pos=axs[4].get_position();axs[4].set_axis_off()
top=f.add_axes([pos.x0,pos.y0+pos.height*.61,pos.width,pos.height*.36]);bottom=f.add_axes([pos.x0,pos.y0+.003,pos.width,pos.height*.36])
forest(top,baseline,['Z1031 (98)','POETIC (144)']);top.set_xlabel('Adjusted on-Ki67 coefficient / SD import',fontsize=6);top.tick_params(labelsize=5.8)
forest(bottom,correlations,['Z1031 (98)','POETIC (144)']);bottom.set_xlabel('Import/Ki67 change | Spearman rho',fontsize=6);bottom.tick_params(labelsize=5.8)
z=data('Figure_6_h').rename(columns={'HR':'estimate'});panel(axs[5],n,'f','PFS | identical SET-complete patients',z,'Figure_6_h');forest(axs[5],z,['Unadjusted','Proliferation','Prolif. + SET'],ref=1);axs[5].set_xlabel('Import HR per SD | 95% Cox interval');note(axs[5],'47 patients / 44 events; final original q=0.092');save(f,n)

# Preserve established supplementary artwork. Relabel moved figures in vector PDF/SVG,
# then render their PNGs from the relabelled vector PDFs using Poppler.
prior=pd.read_csv(old/'panel_source_data_manifest.csv')
mapping={f'Supplementary_Figure_S{i}':f'Supplementary_Figure_S{i}' for i in range(1,13)}
mapping.update({'Figure_2':'Supplementary_Figure_S13','Figure_7':'Supplementary_Figure_S14','Figure_8':'Supplementary_Figure_S15'})
pop=Path(shutil.which('pdftoppm') or 'pdftoppm')
header_audit=[]
for original,new in mapping.items():
 svg=ET.fromstring((old/'figures'/f'{original}.svg').read_text(encoding='utf-8'));parents={child:parent for parent in svg.iter() for child in parent};headers=[]
 for element in list(svg.iter()):
  if element.tag.endswith('}text') and float(element.attrib.get('y','1000'))<65:
   headers.append(''.join(element.itertext()).strip());parent=parents[element]
   if parent.tag.endswith('}g') and parent.attrib.get('id','').startswith('text_'):parents[parent].remove(parent)
   else:parent.remove(element)
 view=list(map(float,svg.attrib['viewBox'].split()));cut=65.;svg.set('viewBox',f'{view[0]} {cut} {view[2]} {view[3]-cut}');svg.set('height',str(view[3]-cut)+'pt');ET.register_namespace('','http://www.w3.org/2000/svg');ET.ElementTree(svg).write(figdir/f'{new}.svg',encoding='utf-8',xml_declaration=True)
 pdf=PdfReader(old/'figures'/f'{original}.pdf');page=pdf.pages[0];width,height=float(page.mediabox.width),float(page.mediabox.height);stream=ContentStream(page.get_contents(),pdf);ops=[];block=[];inside=False;removed=[]
 for args,op in stream.operations:
  if op==b'BT':inside=True;block=[(args,op)]
  elif inside:
   block.append((args,op))
   if op==b'ET':
    text=''.join(str(value) for operands,operator in block if operator in [b'Tj',b'TJ'] for value in (operands[0] if operator==b'TJ' else operands) if isinstance(value,str)).strip()
    if text in headers:removed.append(text)
    else:ops.extend(block)
    inside=False
  else:ops.append((args,op))
 assert len(removed)==len(headers) and len(headers)>=2,(original,headers,removed)
 stream.operations=ops;page[NameObject('/Contents')]=stream;page.mediabox.upper_right=(width,height-cut);page.cropbox.upper_right=(width,height-cut);w=PdfWriter();w.add_page(page)
 with (figdir/f'{new}.pdf').open('wb') as output: w.write(output)
 subprocess.run([str(pop),'-singlefile','-r','350','-png',str(figdir/f'{new}.pdf'),str(figdir/new)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 header_audit.append(dict(figure=new,header_text_blocks_removed=len(removed),crop_points=cut,removed_text=' | '.join(headers)))
 for rec in prior[prior.figure==original].to_dict('records'):
  oldfile=rec['source_data'];newfile=oldfile.replace(original,new,1);shutil.copy2(old/'source_data'/oldfile,src/newfile);manifest.append(dict(figure=new,panel=rec['panel'],title=rec['title'],source_data=newfile));provenance.append(dict(figure=new,panel=rec['panel'],revision10_source=original+'_'+rec['panel']))

n='Supplementary_Figure_S16';f,axs=canvas(n,'Additional ER, protein and clinical controls','Complete primary families and contextual sensitivities are retained in result tables')
z=data('Figure_4_a');panel(axs[0],n,'a','ESR1 knockdown | target engagement',z,'Figure_4_a')
for i,key in enumerate(['ESR1','ER_targets_no_ESR1','Proliferation']):
 baseline=z.loc[z.group=='CTL',key].mean()
 for j,g in enumerate(['CTL','ESR1']):
  vals=z.loc[z.group==g,key]-baseline;axs[0].scatter(i+(j-.5)*.22+rng.uniform(-.03,.03,len(vals)),vals,s=12,c=[GREY,BLUE][j])
axs[0].set_xticks(range(3),['ESR1','ER targets','Proliferation'],rotation=15);axs[0].set_ylabel('Score minus control mean');zero(axs[0]);note(axs[0],'Three cultures/group')
z=data('Figure_4_d');panel(axs[1],n,'b','All three fixed ER-loss import contrasts',z,'Figure_4_d');forest(axs[1],z,['siESR1 - control','E2 + Fulv - E2','E2 x Fulv']);axs[1].set_xlabel('Source-scale contrast; 95% interval');note(axs[1],'Original source q: 0.060 / 0.023 / 0.060')
z=data('Figure_4_e');panel(axs[2],n,'c','Small-n ER contexts | descriptive',z,'Figure_4_e');axs[2].scatter(z.estimate,np.arange(len(z)),s=16,c=TEAL);axs[2].axvline(0,c=GREY,lw=.6);axs[2].set_yticks(range(len(z)),['ER KO 4 d','Giredestrant 4 d','siESR1 48 h','siESR1 72 h','Giredestrant mock','GDC-0927 mock'],fontsize=5.1);axs[2].invert_yaxis();axs[2].set_xlabel('Treated minus source control');note(axs[2],'Two source cultures/group; no P values')
z=data('Figure_6_f');panel(axs[3],n,'d','Fixed Ki67 models | prediction error',z,'Figure_6_f')
for i,(acc,v) in enumerate(z.groupby('accession',sort=False)):
 vals=v.set_index('model').loc[['covariates_only','plus_import'],'LOOCV_RMSE'];axs[3].plot([0,1],vals,'-o',c=palette[i],ms=3,lw=1.2,label='Z1031' if acc=='GSE87411' else 'POETIC')
axs[3].set_xticks([0,1],['Covariates','+ import']);axs[3].set_ylabel('Leave-one-out RMSE');axs[3].legend(frameon=False,fontsize=6);note(axs[3],'Fixed models; lower is better')
z=data('Figure_6_e');panel(axs[4],n,'e','POETIC | treatment interaction',z,'Figure_6_e');forest(axs[4],z,z.scale.tolist());axs[4].set_xlabel('Baseline import x AI coefficient');note(axs[4],'195 cases; 144 AI and 51 untreated')
z=data('Figure_6_g');panel(axs[5],n,'f','Clinical endpoint eligibility',z,'Figure_6_g');axs[5].barh(z.cohort,z.paired,color='#DCE5EA',height=.55,label='Paired');axs[5].barh(z.cohort,z.complete,color=BLUE,height=.55,label='Complete IHC/model');axs[5].invert_yaxis();axs[5].set_xlabel('Patients');axs[5].legend(frameon=False,fontsize=5.2)
z=est[(est.scale=='source')&(est.module=='Import_TIM_TOM')];panel(axs[6],n,'g','All ten fixed ESR1 import contrasts',z);forest(axs[6],z,[v.cell_line+' '+v.contrast for v in z.itertuples()]);axs[6].tick_params(axis='y',labelsize=4.6);axs[6].set_xlabel('Source-scale factorial contrast');note(axs[6],'Same ten-test family; cultures not new derivations')
z=data('Figure_5_g');panel(axs[7],n,'h','CPTAC | IHC ER status',z,'Figure_5_g');rain(axs[7],{'ER- (39)':z.loc[~z.ER_positive,'Import_protein'],'ER+ (81)':z.loc[z.ER_positive,'Import_protein']},'Import protein source score');save(f,n)

n='Supplementary_Figure_S17';f,axs=canvas(n,'Protein-synthesis blockade','GSE8597 | fixed factorial contrasts | source/rank sensitivity')
samples=pd.read_csv(out/'CHX_sample_scores.csv');tests=pd.read_csv(out/'CHX_all_module_tests.csv');specific=pd.read_csv(out/'CHX_specificity_tests.csv');genes=pd.read_csv(out/'CHX_gene_descriptive_responses.csv');coverage=pd.read_csv(out/'CHX_gene_coverage.csv')
z=samples[samples.scale=='source'];panel(axs[0],n,'a','E2 and cycloheximide | all RNA samples',z);dots(axs[0],z,['Vehicle','E2','CHX','E2+CHX']);axs[0].set_ylabel('Mean deposited RMA log2');note(axs[0],'Four RNA samples/condition; 22/24 genes')
z=tests[(tests.scale=='source')&(tests.module=='Import_TIM_TOM')];panel(axs[1],n,'b','All fixed import contrasts | source scale',z);forest(axs[1],z,['E2 without CHX','E2 during CHX','Interaction']);axs[1].set_xlabel('Source contrast; 95% HC3 interval');note(axs[1],'Three-test q: 0.0018 / 0.444 / 0.101')
z=tests[(tests.scale=='rank')&(tests.module=='Import_TIM_TOM')];panel(axs[2],n,'c','All fixed import contrasts | rank scale',z);forest(axs[2],z,['E2 without CHX','E2 during CHX','Interaction']);axs[2].set_xlabel('Rank contrast; 95% HC3 interval');note(axs[2],'Three-test q: 0.0052 / 0.561 / 0.062')
z=genes[genes.gene.isin(mods)];panel(axs[3],n,'d','All fixed genes | descriptive E2 responses',z);hh=z.pivot(index='gene',columns='CHX',values='E2_mean_response').reindex(mods);im=axs[3].imshow(hh,aspect='auto',cmap='RdBu_r',vmin=-.5,vmax=.5);axs[3].set_xticks([0,1],['No CHX','With CHX']);axs[3].set_yticks(range(len(hh)),hh.index,fontsize=4.5);plt.colorbar(im,ax=axs[3],fraction=.04,pad=.03).set_label('E2 minus matched vehicle',size=6)
z=tests[(tests.scale=='source')&(tests.contrast=='E2_by_CHX')];panel(axs[4],n,'e','Interaction | all fixed programmes',z);forest(axs[4],z,[v.replace('Import_TIM_TOM','Import').replace('Nonimport_mito','Other mito').replace('Cytosolic_ribosome','Ribosome').replace('Translation_factors','Translation').replace('ER_programme','ER targets') for v in z.module]);axs[4].set_xlabel('E2-by-CHX interaction; 95% interval');note(axs[4],'Separate three-test family per programme')
z=specific[(specific.scale=='source')&(specific.contrast=='E2_by_CHX')];panel(axs[5],n,'f','Direct interaction specificity | all comparators',z);forest(axs[5],z,['Other mito','Ribosome','Translation','Proliferation']);axs[5].set_xlabel('Import minus comparator interaction');note(axs[5],'Twelve-test family per scale');save(f,n)
n='Supplementary_Figure_S18';f,axs=canvas(n,'NRF1 perturbation and rescue','All prespecified conditions, contrasts and orthogonal platforms')
tests=pd.read_csv(out/'NRF1_import_tests.csv');samples=pd.read_csv(out/'NRF1_sample_scores.csv');genes=pd.read_csv(out/'NRF1_gene_descriptive_responses.csv');spec=pd.read_csv(out/'NRF1_specificity_tests.csv')
levels=['siNT','siNRF1','siNT+NRF1-FL','siNRF1+NRF1-FL','siNT+NRF1-dE7','siNRF1+NRF1-dE7'];labels=['Control','KD','Ctrl+FL','KD+FL','Ctrl+dE7','KD+dE7']
z=samples[(samples.accession=='GSE292988')&(samples.scale=='source')];panel(axs[0],n,'a','Import RNA | all six rescue conditions',z,'revision12_NRF1_sample_scores');dots(axs[0],z,levels);axs[0].set_xticks(range(6),labels,rotation=35,ha='right',fontsize=5.5);axs[0].set_ylabel('Mean log2(scaled abundance + 1)');note(axs[0],'HEK293 Flp-In; three RNA samples/group')
mat=pd.read_csv(out/'GSE292988_gene_log2_normalized_counts.csv.gz',index_col=0);z=z.copy();z['NRF1_log2_abundance']=z.column.map(mat.loc['NRF1']);z['NRF1_delta_control']=z.NRF1_log2_abundance-z.loc[z.group=='siNT','NRF1_log2_abundance'].mean();panel(axs[1],n,'b','NRF1 RNA | construct expression restored',z,'revision12_NRF1_target_check')
for i,g in enumerate(levels):
 v=z.loc[z.group==g,'NRF1_delta_control'];axs[1].scatter(i+rng.uniform(-.07,.07,len(v)),v,c=palette[i%4],s=14);axs[1].plot([i-.13,i+.13],[v.mean()]*2,c=INK,lw=1)
axs[1].set_xticks(range(6),labels,rotation=35,ha='right',fontsize=5.5);axs[1].set_ylabel('NRF1 RNA minus control mean');zero(axs[1]);note(axs[1],'Total RNA includes the rescue construct')
for ax,letter,scale in zip(axs[2:4],'cd',['source','rank']):
 z=tests[(tests.accession=='GSE292988')&(tests.scale==scale)&(tests.family=='primary')];panel(ax,n,letter,'All primary contrasts | '+scale,z,'revision12_NRF1_primary_tests');forest(ax,z,['KD - control','FL rescue interaction','Deleted - FL rescue']);ax.set_xlabel(scale.capitalize()+' score contrast; 95% HC3 interval');note(ax,'Three-test BH family; no module rescue')
z=genes[(genes.accession=='GSE292988')&genes.gene.isin(mods)].copy();base=z[z.group=='siNT'].set_index('gene').mean_log2_normalized_count;z['delta_control']=z.mean_log2_normalized_count-z.gene.map(base);panel(axs[4],n,'e','All fixed import genes | mean responses',z,'revision12_NRF1_gene_responses');hh=z.pivot(index='gene',columns='group',values='delta_control').reindex(index=mods,columns=levels);im=axs[4].imshow(hh,aspect='auto',cmap='RdBu_r',vmin=-.8,vmax=.8);axs[4].set_yticks(range(len(hh)),hh.index,fontsize=4.3);axs[4].set_xticks(range(6),labels,rotation=45,ha='right',fontsize=4.5);plt.colorbar(im,ax=axs[4],fraction=.04,pad=.03).set_label('Log2 RNA change vs control',size=5.5)
def dual_forest(ax,z,labels,splitkey,xlabel):
 pos=ax.get_position();ax.set_axis_off();children=[]
 for i,scale in enumerate(['source','rank']):
  child=f.add_axes([pos.x0,pos.y0+pos.height*(.59 if i==0 else .01),pos.width,pos.height*.36]);sub=z[z.scale==scale];forest(child,sub,labels);child.tick_params(labelsize=4.8);child.set_xlabel(scale.capitalize()+' '+xlabel,fontsize=5.3,labelpad=1);children.append(child)
 return children
z=tests[tests.accession=='GSE292986'];panel(axs[5],n,'f','Endogenous exon editing | all platforms',z,'revision12_NRF1_editing');dual_forest(axs[5],z,['SaCas9 ABE','SpCas9 ABE','SpCas9 CBE'],'contrast','score change')
z=spec[(spec.accession=='GSE292988')&(spec.contrast=='KD_vs_control')];panel(axs[6],n,'g','Depletion | direct programme specificity',z,'revision12_NRF1_specificity');dual_forest(axs[6],z,['Other mito','Ribosome','Translation','Proliferation'],'comparator','import minus comparator')
z=pd.read_csv(root/'revision10/analysis/CDK_control_primary_tests.csv');z=z[(z.scale=='source')&z.module.isin(['Import_TIM_TOM','ER_programme','Proliferation'])];panel(axs[7],n,'h','ER blockade versus CDK4/6 inhibition',z,'revision11_polished_Figure_3_h');forest(axs[7],z,[v.experiment+' | '+{'Import_TIM_TOM':'Import','ER_programme':'ER targets','Proliferation':'Proliferation'}[v.module] for v in z.itertuples()]);axs[7].set_xlabel('OP1 minus CDK25 source score');note(axs[7],'P/R: n=4/3 per group; exposures not equivalent');save(f,n)

pd.DataFrame(header_audit).to_csv(out/'figure_header_removal_audit.csv',index=False)
pd.DataFrame(manifest).to_csv(dest/'panel_source_data_manifest.csv',index=False);pd.DataFrame(provenance).to_csv(dest/'panel_revision_map.csv',index=False)
mean=e[e.scale=='source'].groupby(['module','hour']).delta.agg(['mean','min','max']).reset_index();mean.to_csv(out/'existing_E2_descriptive_kinetics.csv',index=False)
print('Built',len(set(v['figure'] for v in manifest)),'figures;',len(manifest),'linked panels')
