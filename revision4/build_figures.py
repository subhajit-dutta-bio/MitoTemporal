from pathlib import Path
import sys,json,re,shutil,textwrap
r=Path(__file__).parent;root=r.parent;old=root/'revision2';rep=root/'replication/analysis';out=r/'deliverables';figdir=out/'figures';src=out/'source_data';src.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'python_packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.font_manager import fontManager
from scipy.stats import spearmanr
font=Path('C:/Windows/Fonts/arial.ttf')
if font.exists():fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Arial' if font.exists() else 'DejaVu Sans','font.size':7,'axes.titlesize':7.8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#2374AB';ORANGE='#D77935';TEAL='#168C83';PURPLE='#8260A5';GREY='#73808C';PINK='#B05C84';INK='#172B3A'
pal=[BLUE,ORANGE,TEAL,PURPLE,PINK,'#C39C33','#63A8AD','#A78465','#67759D','#AA909E','#829B57','#6697CE','#BF6D6D']
mods=json.loads((root/'source/hummel2/gene_sets/mito_modules.json').read_text());ml=list(mods)
short={'OXPHOS_Complex_I':'C I','OXPHOS_Complex_II':'C II','OXPHOS_Complex_III':'C III','OXPHOS_Complex_IV':'C IV','OXPHOS_Complex_V':'C V','TCA_cycle':'TCA','Membrane_potential_proxy':'Transporters','Import_TIM_TOM':'TIM/TOM','mtDNA_encoded':'mtDNA','Biogenesis':'Biogenesis','Fission':'Fission','Fusion':'Fusion','Mitophagy':'Mitophagy','mtUPR':'mtUPR'}
ctshort={'B_cell':'B','Cancer_Epithelial':'Epithelial','Cycling':'Cycling','DC':'DC','Endothelial':'Endothelial','Fibroblast':'Fibroblast','Macrophage':'Macrophage','Mast':'Mast','Monocyte':'Monocyte','NK':'NK','Plasma':'Plasma','T_CD4':'CD4','T_CD8':'CD8'}
atlas=pd.read_csv(r/'analysis/atlas_display_and_full_summary.csv');full=atlas[atlas.record_type=='full_atlas_group_mean'].copy();display=atlas[atlas.record_type=='display_subsample'].copy()
traj=pd.read_csv(root/'progression_individual_trajectories.csv');prog=pd.read_csv(root/'progression_patient_compartment.csv');pat=pd.read_csv(root/'hummel_patient_scores.csv');treat=pd.read_csv(root/'treatment_patient_compartment.csv')
external=pd.read_csv(rep/'replication_primary_results.csv');extall=pd.read_csv(rep/'replication_all_results.csv');extd=pd.read_csv(rep/'replication_patient_deltas.csv');oldpairs=pd.read_csv(old/'analysis/additional_patient_pairs.csv');oldtests=pd.read_csv(old/'analysis/additional_all_tests.csv');newtests=pd.read_csv(r/'analysis/expansion_all_module_tests.csv');newd=pd.read_csv(r/'analysis/clinical_patient_deltas.csv');miss=pd.read_csv(r/'analysis/GSE87411_paired_observed_tests.csv');missd=pd.read_csv(r/'analysis/GSE87411_paired_observed_deltas.csv');cov=pd.read_csv(r/'analysis/gene_coverage.csv');manifest=[];figtext={};limits=[]

def main_canvas(n,title,subtitle):
 fig=plt.figure(figsize=(7.2,10.7));g=fig.add_gridspec(4,6,left=.115,right=.95,bottom=.07,top=.865,hspace=.83,wspace=1.65)
 spans=[(0,slice(0,3)),(0,slice(3,6)),(1,slice(0,3)),(1,slice(3,6)),(2,slice(0,2)),(2,slice(2,4)),(2,slice(4,6)),(3,slice(0,2)),(3,slice(2,4)),(3,slice(4,6))]
 axes=[fig.add_subplot(g[i,j]) for i,j in spans]
 fig.text(.115,.975,f'FIGURE {n}',size=8,weight='bold',color=BLUE);fig.text(.115,.947,title,size=11.1,weight='bold',color=INK);fig.text(.115,.923,textwrap.fill(subtitle,105),size=6.8,color=GREY)
 return fig,axes
def sup_canvas(n,title):
 fig=plt.figure(figsize=(7.2,9.2));g=fig.add_gridspec(3,6,left=.115,right=.93,bottom=.09,top=.865,hspace=.94,wspace=1.85)
 spans=[(0,slice(0,3)),(0,slice(3,6)),(1,slice(0,2)),(1,slice(2,4)),(1,slice(4,6)),(2,slice(0,2)),(2,slice(2,4)),(2,slice(4,6))]
 axes=[fig.add_subplot(g[i,j]) for i,j in spans];fig.text(.09,.965,f'SUPPLEMENTARY FIGURE {n}',size=8,weight='bold',color=BLUE);fig.text(.09,.93,title,size=11.5,weight='bold',color=INK);return fig,axes
def panel(ax,figure,letter,title,data,note=''):
 width=ax.get_position().width*ax.figure.get_figwidth();wrapped=textwrap.fill(title,max(16,min(46,int(width*18))))
 ax.set_title(wrapped,loc='left',pad=20,fontweight='bold',fontsize=7.4)
 pos=ax.get_position();ax.figure.text(pos.x0-.032,pos.y1+.026,letter,size=10.5,weight='bold',color=INK,ha='right')
 name=f'{figure}_{letter}.csv.gz';data.to_csv(src/name,index=False,compression='gzip');manifest.append({'figure':figure,'panel':letter,'title':title,'source_data':name,'note':note});return ax
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(figdir/(name+'.'+ext),dpi=360)
 plt.close(fig)
def heat(ax,a,cmap='RdBu_r',vmax=None,cbar=True,label='',annotations=None):
 v=a.values.astype(float);vmax=vmax or max(.01,np.nanmax(np.abs(v)));im=ax.imshow(v,aspect='auto',cmap=cmap,vmin=-vmax if cmap=='RdBu_r' else 0,vmax=vmax,interpolation='nearest')
 ax.set_xticks(range(a.shape[1]),[short.get(str(x),str(x)) for x in a.columns],rotation=65,ha='right');ax.set_yticks(range(a.shape[0]),[ctshort.get(str(x),str(x)) for x in a.index]);ax.tick_params(length=0,pad=2);ax.spines[['left','bottom']].set_visible(False)
 if cbar:
  cax=ax.inset_axes([1.025,.05,.025,.9]);cb=plt.colorbar(im,cax=cax);cb.ax.tick_params(labelsize=4.7,pad=1,length=2)
 return im
def includezero(ax):
 a,b=ax.get_ylim();ax.set_ylim(min(0,a),max(0,b));return ax
def paired(ax,a,b,color,label0='Baseline',label1='On treatment',zero=True):
 for x,y in zip(a,b):ax.plot([0,1],[x,y],color=color,alpha=.24,lw=.5,zorder=1)
 ax.scatter(np.zeros(len(a)),a,s=5,color=color,alpha=.45);ax.scatter(np.ones(len(b)),b,s=5,color=color,alpha=.45)
 ax.plot([0,1],[np.mean(a),np.mean(b)],'-o',color=INK,lw=1.5,ms=3,zorder=4);ax.set_xticks([0,1],[label0,label1]);ax.set_xlim(-.2,1.2)
 if zero:includezero(ax)
def note(ax,s):ax.text(.005,1.012,s,transform=ax.transAxes,va='bottom',size=5.4,color=GREY,clip_on=False)
def strip(ax,groups,colors=None,ylabel='Centered module expression'):
 for i,(name,v) in enumerate(groups.items()):
  v=np.asarray(v);color=(colors or pal)[i%len(colors or pal)];ax.scatter(i+np.linspace(-.065,.065,len(v)),v,s=18,c=color,edgecolors='white',linewidths=.35,zorder=3);ax.plot([i-.18,i+.18],[v.mean()]*2,c=INK,lw=1.4)
 ax.set_xticks(range(len(groups)),list(groups),rotation=30,ha='right');ax.axhline(0,color='#CBD1D5',lw=.7);ax.set_ylabel(ylabel);ax.margins(x=.22);includezero(ax)
def forest(ax,df,labels,color=BLUE):
 y=np.arange(len(df));ax.errorbar(df.mean_delta,y,xerr=np.array([df.mean_delta-df.ci_low,df.ci_high-df.mean_delta]),fmt='o',color=color,ms=3.4,capsize=2,lw=1);ax.axvline(0,color=GREY,lw=.7);ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Mean paired change (95% bootstrap CI)');lo,hi=ax.get_xlim();ax.set_xlim(min(lo,0),max(hi,0))
def clinical_delta(acc,module,scale='source'):
 if acc in ['GSE20181','GSE111563','GSE147271']:
  contrast={'GSE20181':'C_minus_A','GSE111563':'4_minus_1','GSE147271':'Tamoxifen_minus_Pretreatment'}[acc]
  if scale=='source':return oldpairs[(oldpairs.accession==acc)&(oldpairs.contrast==contrast)&(oldpairs.module==module)].set_index('patient').delta
  scores=pd.read_csv(old/'analysis'/f'{acc}_rank_scores.csv',index_col=0);md=pd.read_csv(old/'analysis/additional_clinical_sample_manifest.csv',dtype=str);md=md[md.accession==acc].set_index('sample');z=scores.join(md[['patient','timepoint']]);a=z.groupby(['patient','timepoint'])[module].mean().unstack();pre,post={'GSE20181':('A','C'),'GSE111563':('1','4'),'GSE147271':('Pretreatment','Tamoxifen')}[acc];a=a[[pre,post]].dropna();return a[post]-a[pre]
 a=newd if acc=='GSE153470' else missd;return a[(a.accession==acc)&(a.module==module)&(a.scale==scale)].set_index('patient').delta
def clinical_test(acc,module,scale='source'):
 if acc in ['GSE20181','GSE111563','GSE147271']:
  contrast={'GSE20181':'C_minus_A','GSE111563':'4_minus_1','GSE147271':'Tamoxifen_minus_Pretreatment'}[acc];return oldtests[(oldtests.accession==acc)&(oldtests.contrast==contrast)&(oldtests.module==module)&(oldtests.scale==('mean_log' if scale=='source' else 'rank'))].iloc[0]
 t=newtests if acc=='GSE153470' else miss;return t[(t.accession==acc)&(t.module==module)&(t.scale==scale)].iloc[0]

# Figure 1: original coordinates, full-atlas summaries and donor-level progression.
f,axs=main_canvas(1,'An atlas resolved by compartment and patient','708,338 cells | 14 predefined programmes | progression remains exploratory')
cts=sorted(display.cell_type.unique());ctpal=dict(zip(cts,pal))
ax=panel(axs[0],'Figure_1','a','Original integrated UMAP',display,'80,000 cells sampled uniformly with seed 20261024; original coordinates, no re-embedding')
for ct in cts:
 d=display[display.cell_type==ct];ax.scatter(d.UMAP1,d.UMAP2,s=.35,c=ctpal[ct],rasterized=True,alpha=.65,linewidths=0);ax.text(d.UMAP1.median(),d.UMAP2.median(),ctshort[ct],size=5.8,ha='center',bbox={'facecolor':'white','edgecolor':'none','alpha':.6,'pad':.5})
ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
ax=panel(axs[1],'Figure_1','b','Cohort provenance on the same UMAP',display[['UMAP1','UMAP2','dataset']])
for i,acc in enumerate(['GSE167036','GSE169246','GSE176078']):
 d=display[display.dataset==acc];ax.scatter(d.UMAP1,d.UMAP2,s=.32,c=pal[i],label=acc,rasterized=True,linewidths=0,alpha=.5)
ax.legend(fontsize=5.8,frameon=False,markerscale=6,loc='upper right',handletextpad=.2);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
w=full.copy();cols=['mito_'+v for v in ml];z=pd.DataFrame({ct:np.average(g[cols],axis=0,weights=g.n_cells) for ct,g in w.groupby('cell_type')},index=ml).T;z=(z-z.mean())/z.std(ddof=0)
ax=panel(axs[2],'Figure_1','c','Mitochondrial programmes across annotated cell types',z.reset_index(),'Full-atlas cell-weighted means, column z scores, no inferential tests');heat(ax,z,vmax=2.5,label='Column z score')
counts=full.groupby('dataset').n_cells.sum().reindex(['GSE167036','GSE169246','GSE176078']);ax=panel(axs[3],'Figure_1','d','Atlas size',counts.reset_index());ax.barh(range(3),counts/1000,color=pal[:3]);ax.set_yticks(range(3),[x.replace('GSE','') for x in counts.index]);ax.set_xlabel('Cells (thousands)');ax.invert_yaxis()
for i,v in enumerate(counts):ax.text(v/1000,i,f' {int(v):,}',va='center',size=5.8)
ax.set_xlim(0,570)
cc=full.pivot(index='cell_type',columns='dataset',values='n_cells').fillna(0);cc=cc.div(cc.sum())*100;ax=panel(axs[4],'Figure_1','e','Cell-type composition (%)',cc.reset_index());heat(ax,cc,cmap='Blues',vmax=60,label='% cells');ax.set_xticklabels(['167036','169246','176078'],rotation=45,ha='right')
for ax,ct,letter in zip(axs[5:7],['T_CD4','T_CD8'],['f','g']):
 d=traj[(traj.compartment=='blood')&(traj.cell_type==ct)&(traj.module=='OXPHOS_Complex_IV')];panel(ax,'Figure_1',letter,ctshort[ct]+' blood progression',d)
 for i,(_,z) in enumerate(d.iterrows()):ax.plot([0,1,2],z[['Pre','Post','Prog']].astype(float),'-o',color=pal[i],ms=2.1,lw=.7,alpha=.65)
 ax.plot([0,1,2],d[['Pre','Post','Prog']].mean(),'-o',color=INK,lw=1.6,ms=3);ax.set_xticks([0,1,2],['Pre','Post','Prog']);ax.set_ylabel('Complex IV score');ax.set_ylim(0,1.05);note(ax,'n=7; Prog-Pre q=0.410')
coverage=traj.groupby(['compartment','cell_type']).patient.nunique().unstack(0).fillna(0);ax=panel(axs[7],'Figure_1','h','Complete trajectory coverage',coverage.reset_index());y=np.arange(len(coverage));ax.barh(y-.15,coverage.blood,.3,label='Blood',color=BLUE);ax.barh(y+.15,coverage.tissue,.3,label='Tissue',color=ORANGE);ax.set_yticks(y,[ctshort.get(k,k) for k in coverage.index]);ax.invert_yaxis();ax.set_xlabel('Patients with all 3 timepoints');ax.legend(fontsize=5.5,frameon=False)
d=prog[(prog.compartment=='blood')&(prog.contrast=='Prog_minus_Pre')&(prog.n_pairs>=3)];a=d.pivot(index='cell_type',columns='module',values='mean_delta').reindex(columns=ml);ax=panel(axs[8],'Figure_1','i','Blood progression across modules',d);heat(ax,a,vmax=.18,label='Mean score change');ax.tick_params(axis='x',labelsize=4.8);ax.tick_params(axis='y',labelsize=5.3)
d=d[d.module=='OXPHOS_Complex_IV'];ax=panel(axs[9],'Figure_1','j','Complex IV across blood lineages',d);forest(ax,d,[ctshort.get(v,v) for v in d.cell_type]);ax.set_xlabel('Prog-Pre mean change (95% CI)',fontsize=6);note(ax,'Exploratory; all q>0.05')
save(f,'Figure_1')

# Figure 2: all external longitudinal blood results and sensitivity.
f,axs=main_canvas(2,'External blood cohorts do not confirm progression','Scheduled treatment endpoints | source-defined patient pairs | RNA-marker-defined T-cell subsets')
for ax,acc,lineage,letter in zip(axs[:4],['GSE205589','GSE205589','GSE300475','GSE300475'],['CD4','CD8','CD4','CD8'],list('abcd')):
 a=pd.read_csv(rep/f'{acc}_patient_timepoint_scores.csv');d=a[a.lineage==lineage];panel(ax,'Figure_2',letter,f'{acc} | {lineage}',d)
 times=['Pre','Post'] if acc=='GSE205589' else ['Week1','Week3','Week7','Week15'];pv=d.pivot(index='patient',columns='timepoint',values='mean').reindex(columns=times)
 for i,(_,z) in enumerate(pv.iterrows()):ax.plot(range(len(times)),z,'-o',color=pal[i],ms=2.2,alpha=.7)
 ax.plot(range(len(times)),pv.mean(),'-o',color=INK,lw=1.5,ms=3);ax.set_xticks(range(len(times)),times,rotation=0 if len(times)==2 else 30);ax.set_ylabel('Complex IV module score');ax.set_ylim(0,1);rr=external[(external.accession==acc)&(external.lineage==lineage)].iloc[0];note(ax,f'n={len(pv)}; endpoint q={rr.q:.3f}')
ax=panel(axs[4],'Figure_2','e','Endpoint effects and uncertainty',external);forest(ax,external,[a.replace('GSE','')+' '+b for a,b in zip(external.accession,external.lineage)]);ax.set_xlabel('Mean change (95% CI)',fontsize=6)
cc=pd.read_csv(rep/'replication_marker_subset_coverage.csv');cc=cc.groupby(['accession','patient','timepoint','lineage']).n_cells.sum().reset_index();cc['sample']=cc.accession.str.replace('GSE','')+' '+cc.patient+' '+cc.timepoint;hh=cc.pivot(index='sample',columns='lineage',values='n_cells').fillna(0);ax=panel(axs[5],'Figure_2','f','Cells contributing per sample',cc);heat(ax,np.log10(hh+1),cmap='Blues',vmax=4.2,label='log10(cells+1)');ax.set_yticks([]);ax.set_ylabel('Patient/timepoint samples');ax.set_xticklabels(hh.columns,rotation=0)
z=extall[(extall.analysis!='primary')&extall.lineage.isin(['CD4','CD8'])].copy();z['cohort_lineage']=z.accession.str.replace('GSE','')+' '+z.lineage;hh=z.pivot(index='analysis',columns='cohort_lineage',values='q');ax=panel(axs[6],'Figure_2','g','Sensitivity: corrected support',z);heat(ax,-np.log10(hh),cmap='Purples',vmax=1.3,label='-log10(q)');qlabels={'demux_dom10':'Demux 10','demux_dom5':'Demux 5','exclude_PT15_rep_S8':'Omit S8','mean_log':'Mean log','min50':'50 cells','mt10':'mt<10%','mt40':'mt<40%','nuclear':'Nuclear','nuclear_mean_log':'Nuc log'};ax.set_yticklabels([qlabels.get(x,x) for x in hh.index],size=5);ax.tick_params(axis='x',labelsize=4.8);note(ax,'All sensitivity q>0.05')
nu=extall[(extall.analysis.isin(['primary','nuclear']))&extall.lineage.isin(['CD4','CD8'])].copy();ax=panel(axs[7],'Figure_2','h','Full versus nuclear-only scores',nu);forest(ax,nu,[f'{a[-3:]} {b} '+('Full' if c=='primary' else 'Nuc') for a,b,c in zip(nu.accession,nu.lineage,nu.analysis)]);ax.set_xlabel('Mean change (95% CI)',fontsize=6);ax.tick_params(axis='y',labelsize=5)
dd=extd.pivot(index=['accession','patient'],columns='lineage',values='delta').reset_index();ax=panel(axs[8],'Figure_2','i','Paired CD4/CD8 response variation',dd)
for i,(acc,z) in enumerate(dd.groupby('accession')):ax.scatter(z.CD4,z.CD8,s=20,color=pal[i],label=acc.replace('GSE',''))
ax.axhline(0,color=GREY,lw=.6);ax.axvline(0,color=GREY,lw=.6);ax.set_xlabel('CD4 patient change');ax.set_ylabel('CD8 patient change');ax.legend(fontsize=5.5,frameon=False)
qc=[]
for threshold,file in [('Original','progression_patient_compartment.csv'),('<10% mt','sensitivity_10_progression.csv'),('<20% mt','sensitivity_20_progression.csv'),('<25% mt','sensitivity_25_progression.csv'),('<40% mt','sensitivity_40_progression.csv'),('50 cells','sensitivity_50_cells_progression.csv')]:
 a=pd.read_csv(root/file);a=a[(a.compartment=='blood')&(a.module=='OXPHOS_Complex_IV')&(a.contrast=='Prog_minus_Pre')&a.cell_type.isin(['T_CD4','T_CD8'])].copy();a['criterion']=threshold;qc.append(a)
qc=pd.concat(qc);ax=panel(axs[9],'Figure_2','j','Discovery QC sensitivity',qc)
for i,(ct,z) in enumerate(qc.groupby('cell_type')):ax.plot(range(len(z)),z.mean_delta,'-o',ms=3,color=pal[i],label=ctshort[ct])
ax.set_xticks(range(6),['Orig','10%','20%','25%','40%','50 cells'],rotation=55,ha='right');ax.axhline(0,color=GREY,lw=.7);ax.set_ylim(bottom=0);ax.set_ylabel('Prog-Pre Complex IV change');ax.legend(fontsize=5.5,frameon=False);note(ax,'No corrected progression result')
save(f,'Figure_2')

# Figure 3: clinical validation and gene-level support.
f,axs=main_canvas(3,'Endocrine therapy lowers import-associated transcription','Patient-level effects across distinct studies | source-scale and rank scoring | bulk tumour context')
for ax,acc,letter in zip(axs[:4],['GSE20181','GSE111563','GSE153470','GSE87411'],list('abcd')):
 d=clinical_delta(acc,'Import_TIM_TOM');rr=clinical_test(acc,'Import_TIM_TOM');panel(ax,'Figure_3',letter,acc+' | TIM/TOM',d.rename('delta').reset_index(),'GSE87411 uses paired-observed-gene sensitivity; other panels use fixed complete-gene scores')
 ax.scatter(np.arange(len(d))+1,d,s=9,color=BLUE,alpha=.65,linewidths=0);ax.axhline(0,color=GREY,lw=.7);ax.axhline(d.mean(),color=ORANGE,lw=1.2);ax.set_xlabel('Patients ordered by deposited identifier');ax.set_ylabel('Paired change (source units)');note(ax,f'n={len(d)}; q={rr.q:.2g}'+('\nPaired-gene sensitivity' if acc=='GSE87411' else ''));lo,hi=ax.get_ylim();ax.set_ylim(lo,hi+(hi-lo)*.18)
cohorts=['GSE20181','GSE111563','GSE153470','GSE87411'];dz=pd.DataFrame(index=ml,columns=cohorts,dtype=float)
for acc in cohorts:
 for module in ml:
  try:d=clinical_delta(acc,module);dz.loc[module,acc]=d.mean()/d.std(ddof=1)
  except (KeyError,IndexError):pass
ax=panel(axs[4],'Figure_3','e','Programme-wide paired shifts',dz.reset_index(),'Mean paired change divided by SD of patient changes; missing/ineligible modules blank');heat(ax,dz,vmax=1.3,label='Standardized paired change');ax.set_yticklabels([short[v] for v in dz.index]);ax.set_xticklabels([v.replace('GSE','') for v in dz.columns],rotation=50,ha='right')
rr=[]
for acc in cohorts:
 for scale in ['source','rank']:
  v=clinical_test(acc,'Import_TIM_TOM',scale);rr.append({'cohort':acc,'scale':scale,'q':v.q,'p':v.p})
rr=pd.DataFrame(rr);ax=panel(axs[5],'Figure_3','f','Import support across score scales',rr)
for k,scale in enumerate(['source','rank']):
 z=rr[rr.scale==scale];ax.scatter(-np.log10(z.q),np.arange(4)+(k-.5)*.18,s=24,color=pal[k],label=scale)
ax.axvline(-np.log10(.05),ls='--',color=GREY,lw=.7);ax.set_yticks(range(4),[x.replace('GSE','') for x in cohorts]);ax.invert_yaxis();ax.set_xlim(left=0);ax.set_xlabel('-log10(corrected q)');ax.legend(fontsize=5.5,frameon=False)
d=oldpairs[(oldpairs.accession=='GSE20181')&(oldpairs.module=='Import_TIM_TOM')&oldpairs.contrast.isin(['B_minus_A','C_minus_A'])].copy();pv=d.pivot(index='patient',columns='contrast',values='delta').dropna();ax=panel(axs[6],'Figure_3','g','Import change over letrozole time',pv.reset_index())
for _,z in pv.iterrows():ax.plot([0,1,2],[0,z.B_minus_A,z.C_minus_A],color=BLUE,alpha=.15,lw=.5)
ax.plot([0,1,2],[0,pv.B_minus_A.mean(),pv.C_minus_A.mean()],'-o',color=INK,ms=3,lw=1.5);ax.set_xticks([0,1,2],['Baseline','Early','Late'],rotation=30);ax.set_ylabel('TIM/TOM change from baseline');ax.axhline(0,color=GREY,lw=.6);note(ax,f'n={len(pv)} complete trajectories')
d=clinical_delta('GSE147271','Mitophagy');rr=clinical_test('GSE147271','Mitophagy');ax=panel(axs[7],'Figure_3','h','Tamoxifen: mitophagy transcripts',d.rename('delta').reset_index());ax.scatter(range(1,len(d)+1),d,s=13,color=TEAL);ax.axhline(0,c=GREY,lw=.7);ax.axhline(d.mean(),c=ORANGE,lw=1);ax.set_xlabel('Patient identifier order');ax.set_ylabel('Paired log-expression change');note(ax,f'n=28; q={rr.q:.3f}')
a=pd.read_pickle(r/'analysis/GSE153470_gene_expression.pkl');pairs=pd.read_csv(r/'analysis/GSE153470_pair_crosswalk.csv',dtype=str).set_index('patient');gs=[g for g in mods['Import_TIM_TOM'] if g in a.index];gd=a.loc[gs,pairs.on].set_axis(pairs.index,axis=1)-a.loc[gs,pairs.baseline].set_axis(pairs.index,axis=1);ax=panel(axs[8],'Figure_3','i','Anastrozole import-gene coherence',gd.reset_index(),'Genes in predefined order, patients in identifier order; source VST/RSN differences');heat(ax,gd,vmax=.5,label='Paired source-unit change');ax.set_xticks([]);ax.tick_params(axis='y',labelsize=4.8);ax.set_xlabel('81 patient pairs')
poetic=pd.read_csv(r/'analysis/POETIC_source_patient_deltas.csv',index_col=0);ax=panel(axs[9],'Figure_3','j','POETIC untreated-biopsy control',poetic[['Import_TIM_TOM','arm']].reset_index(),'Post-primary follow-up; one trial, 157 treated and 56 control patient pairs');strip(ax,{'AI':poetic.loc[poetic.arm=='AI','Import_TIM_TOM'],'Untreated':poetic.loc[poetic.arm=='Control','Import_TIM_TOM']},ylabel='TIM/TOM paired change');note(ax,'n=157/56; source q=0.041\nRank q=0.254')
save(f,'Figure_3')

# Figure 4: perturbations and resistance models, without claiming transcription measures function.
f,axs=main_canvas(4,'Perturbation context separates transcription from function','Source-defined experimental groups | all tested modules retained | culture replication has limited resolution')
ex=newtests[(newtests.family=='experimental_mito')&(newtests.scale=='source')].copy();rows=[]
for _,z in ex.iterrows():
 sc=pd.read_csv(r/'analysis'/f'{z.accession}_source_scores.csv',index_col=0);gm=pd.read_csv(r/'analysis'/f'{z.accession}_group_manifest.csv').set_index('sample').group;tx,co=z.contrast.split('_minus_');aa=sc.loc[gm[gm==co].index,z.module];bb=sc.loc[gm[gm==tx].index,z.module];sd=np.sqrt((aa.var(ddof=1)+bb.var(ddof=1))/2)
 label=z.accession.replace('GSE','')+' '+z.contrast
 if z.accession=='GSE67916':label='MCF7 '+tx
 elif z.accession=='GSE125738':label='T47D TamR'
 elif z.accession=='GSE118713':label='MCF7 '+tx
 elif 'PGC1aOE' in z.accession:label='TIL PGC1a OE'
 elif 'Ndufs4KO' in z.accession:label='TIL: tumour Nu4'
 elif 'ChronicStim' in z.accession:label=tx.replace('Acute','Ac').replace('Chronic','Chr').replace('Hypoxia','Hyp').replace('_','')+'/'+co.replace('Acute','Ac').replace('Chronic','Chr').replace('Hypoxia','Hyp').replace('_','')
 rows.append({'context':label,'accession':z.accession,'contrast':z.contrast,'module':z.module,'standardized_difference':z.mean_delta/sd if sd>0 else np.nan,'q':z.q,'n0':z.n0,'n1':z.n1})
ez=pd.DataFrame(rows);hh=ez.pivot(index='context',columns='module',values='standardized_difference');ax=panel(axs[0],'Figure_4','a','Mitochondrial programmes across perturbations',ez,'Mean difference / pooled sample SD; no pooling across models; color clipped at +/-3');heat(ax,hh,vmax=3,label='Standardized difference');ax.set_yticklabels([v.replace('GSE','').replace('ChronicStimHypoxia ','Stim ').replace('_',' ') for v in hh.index],size=4.9);ax.tick_params(axis='x',labelsize=5.3)
acc='GSE155192_PGC1aOE';ss=pd.read_csv(r/'analysis'/f'{acc}_source_scores.csv',index_col=0);gm=pd.read_csv(r/'analysis'/f'{acc}_group_manifest.csv').set_index('sample').group
ax=panel(axs[1],'Figure_4','b','PGC1a-overexpressing CD8 TIL',ss.reset_index());strip(ax,{name:ss.loc[gm[gm==name].index,'OXPHOS_Complex_IV'].values for name in ['EV','PGC1a_OE']},ylabel='Complex IV expression');note(ax,'3 vs 3 libraries; exact P=0.10\nSource genetic perturbation')
for ax,module,letter in zip(axs[2:4],['Mitophagy','Biogenesis'],['c','d']):
 panel(ax,'Figure_4',letter,'PGC1a: '+short[module],ss[[module]].reset_index());strip(ax,{name:ss.loc[gm[gm==name].index,module].values for name in ['EV','PGC1a_OE']});note(ax,'3 vs 3; no corrected support')
acc='GSE155192_ChronicStimHypoxia';ss=pd.read_csv(r/'analysis'/f'{acc}_source_scores.csv',index_col=0);gm=pd.read_csv(r/'analysis'/f'{acc}_group_manifest.csv').set_index('sample').group;ax=panel(axs[4],'Figure_4','e','Stimulation and hypoxia',ss.reset_index());names=['Acute_Norm','Acute_Hypoxia','Chronic_Norm','Chronic_Hypoxia'];strip(ax,{name.replace('_','\n'):ss.loc[gm[gm==name].index,'OXPHOS_Complex_IV'].values for name in names},ylabel='Complex IV expression');ax.tick_params(axis='x',labelsize=5.2);note(ax,'3 libraries/group')
acc='GSE155192_Ndufs4KO';ss=pd.read_csv(r/'analysis'/f'{acc}_source_scores.csv',index_col=0);gm=pd.read_csv(r/'analysis'/f'{acc}_group_manifest.csv').set_index('sample').group;ax=panel(axs[5],'Figure_4','f','Tumour Ndufs4 deficiency',ss.reset_index());strip(ax,{name.replace('_tumour',''):ss.loc[gm[gm==name].index,'OXPHOS_Complex_IV'].values for name in ['WT_tumour','Ndufs4_def_tumour']},ylabel='CD8 TIL Complex IV expression');note(ax,'Terminally exhausted CD8 TIL\nTumour perturbation; n=3/3')
exr=ex[ex.accession.isin(['GSE67916','GSE125738','GSE118713'])].copy()
for ax,module,letter in zip(axs[6:8],['Import_TIM_TOM','Mitophagy'],['g','h']):
 z=exr[exr.module==module];panel(ax,'Figure_4',letter,'Resistant models: '+short[module],z);y=np.arange(len(z));ax.scatter(z.mean_delta,y,s=22,c=[GREY if pd.isna(q) else TEAL for q in z.q]);ax.axvline(0,c=GREY,lw=.7);labels=[('T47D TR' if a=='GSE125738' else 'MCF7 '+b.split('_minus_')[0]) for a,b in zip(z.accession,z.contrast)];ax.set_yticks(y,labels,size=5.3);ax.invert_yaxis();ax.set_xlabel('Resistant-parental change');note(ax,'Grey: n=2, descriptive' if letter=='g' else 'Shared parental lineages')
d=pd.read_csv(r/'analysis/GSE153470_source_paired_deltas.csv',index_col=0);assoc=pd.read_csv(r/'analysis/clinical_control_associations.csv')
for ax,control,letter in zip(axs[8:10],['Proliferation','Leukocyte'],['i','j']):
 z=assoc[(assoc.target=='Import_TIM_TOM')&(assoc.control==control)].iloc[0];panel(ax,'Figure_4',letter,'Clinical import vs '+control.lower(),d[['Import_TIM_TOM',control]].reset_index());ax.scatter(d[control],d.Import_TIM_TOM,s=12,c=BLUE,alpha=.6,linewidths=0);ax.axhline(0,c=GREY,lw=.6);ax.axvline(0,c=GREY,lw=.6);ax.set_xlabel(control+' marker change');ax.set_ylabel('TIM/TOM change');note(ax,f'n=81; rho={z.rho:.2f}; q={z.q:.3g}')
save(f,'Figure_4')

# Supplementary 1: atlas QC and discovery robustness.
f,axs=sup_canvas('S1','Atlas structure and patient-level robustness')
a=pd.DataFrame([[len(set(mods[x])&set(mods[y]))/len(set(mods[x])|set(mods[y])) for y in ml] for x in ml],index=ml,columns=ml);ax=panel(axs[0],'Supplementary_Figure_S1','a','Predefined gene-set overlap',a.reset_index());heat(ax,a,cmap='Blues',vmax=1,label='Jaccard overlap');ax.set_yticklabels([short[x] for x in a.index],size=5.3)
ax=panel(axs[1],'Supplementary_Figure_S1','b','Mitochondrial RNA by cohort',display[['dataset','pct_counts_mt']]);groups={acc:display.loc[display.dataset==acc,'pct_counts_mt'].values for acc in ['GSE167036','GSE169246','GSE176078']};vp=ax.violinplot(groups.values(),showmedians=True,showextrema=False)
for body,color in zip(vp['bodies'],pal):body.set_facecolor(color);body.set_alpha(.65)
ax.set_xticks([1,2,3],['167036','169246','176078']);ax.set_ylim(bottom=0);ax.set_ylabel('Mitochondrial RNA (%)');note(ax,'80,000-cell display sample')
lo=[]
for ct in ['T_CD4','T_CD8']:
 d=traj[(traj.compartment=='blood')&(traj.cell_type==ct)&(traj.module=='OXPHOS_Complex_IV')].copy();dd=d.set_index('patient').Prog-d.set_index('patient').Pre
 for patient in dd.index:lo.append({'cell_type':ct,'omitted_patient':patient,'mean_delta':dd.drop(patient).mean()})
lo=pd.DataFrame(lo)
for ax,ct,letter in zip(axs[2:4],['T_CD4','T_CD8'],['c','d']):
 z=lo[lo.cell_type==ct];panel(ax,'Supplementary_Figure_S1',letter,ctshort[ct]+' patient omission',z);ax.scatter(z.mean_delta,range(len(z)),s=16,c=BLUE);ax.axvline(0,c=GREY,lw=.7);ax.set_yticks(range(len(z)),z.omitted_patient);ax.set_xlabel('Prog-Pre mean (n=6)');ax.set_xlim(left=0)
ax=panel(axs[4],'Supplementary_Figure_S1','e','Treatment effects across modules',treat);a=treat[(treat.compartment=='blood')&(treat.n_pairs>=3)].pivot(index='cell_type',columns='module',values='mean_delta').reindex(columns=ml);heat(ax,a,vmax=.2,label='Post-Pre mean score change');note(ax,'No patient-level q<0.05')
st=pd.read_csv(root/'source/SUBMISSION_PACKAGE/Results_CSV/module3/mito_state_profiles.csv').set_index('mito_state_id');st=st.reindex(columns=ml);st=(st-st.mean())/st.std(ddof=0);ax=panel(axs[5],'Supplementary_Figure_S1','f','Operational score-cluster profiles',st.reset_index());heat(ax,st,vmax=2.5,label='Column z score');ax.set_ylabel('Gaussian-mixture score cluster')
d=traj[(traj.compartment=='blood')&(traj.module=='OXPHOS_Complex_IV')&traj.cell_type.isin(['T_CD4','T_CD8'])].copy();d['delta']=d.Prog-d.Pre;v=d.pivot(index='patient',columns='cell_type',values='delta');ax=panel(axs[6],'Supplementary_Figure_S1','g','Patient concordance',v.reset_index());ax.scatter(v.T_CD4,v.T_CD8,s=20,c=TEAL);ax.set_xlabel('CD4 Prog-Pre');ax.set_ylabel('CD8 Prog-Pre');ax.set_xlim(left=0);ax.set_ylim(bottom=0)
cc=pat.groupby(['compartment','timepoint']).n_cells.sum().unstack();ax=panel(axs[7],'Supplementary_Figure_S1','h','Sampling by compartment',cc.reset_index());cc[['Pre','Post','Prog']].T.plot.bar(ax=ax,color=[BLUE,ORANGE],width=.8);ax.set_ylabel('Cells');ax.set_xlabel('');ax.legend(fontsize=5.2,frameon=False);ax.tick_params(axis='x',rotation=30)
save(f,'Supplementary_Figure_S1')

# Supplementary 2: gene, scoring and provenance controls.
f,axs=sup_canvas('S2','Gene-level, marker and provenance sensitivity')
gene_rows=[]
for acc in ['GSE20181','GSE111563','GSE153470']:
 a=pd.read_pickle(r/'analysis/GSE153470_gene_expression.pkl') if acc=='GSE153470' else pd.read_csv(old/'analysis/gene_expression'/f'{acc}_module_genes.csv.gz',index_col=0)
 if acc=='GSE153470':pp=pd.read_csv(r/'analysis/GSE153470_pair_crosswalk.csv',dtype=str).set_index('patient');pre=pp.baseline;post=pp.on
 else:
  md=pd.read_csv(old/'analysis/additional_clinical_sample_manifest.csv',dtype=str);md=md[md.accession==acc];pp=md.pivot(index='patient',columns='timepoint',values='sample');x,y=('A','C') if acc=='GSE20181' else ('1','4');pp=pp[[x,y]].dropna();pre=pp[x];post=pp[y]
 gs=[g for g in mods['Import_TIM_TOM'] if g in a.index];gd=a.loc[gs,post].set_axis(pp.index,axis=1)-a.loc[gs,pre].set_axis(pp.index,axis=1)
 for gene,z in gd.iterrows():gene_rows.append({'accession':acc,'gene':gene,'mean_delta':z.mean(),'standardized_change':z.mean()/z.std(ddof=1),'n_pairs':len(z)})
gg=pd.DataFrame(gene_rows);hh=gg.pivot(index='gene',columns='accession',values='standardized_change');ax=panel(axs[0],'Supplementary_Figure_S2','a','Import genes across clinical series',gg);heat(ax,hh,vmax=1,label='Standardized paired gene change');ax.set_xticklabels([v.replace('GSE','') for v in hh.columns],rotation=35,ha='right');ax.tick_params(axis='y',labelsize=5.1)
assoc=pd.read_csv(r/'analysis/clinical_control_associations.csv');hh=assoc[~assoc.control.str.contains('partial')].pivot(index='target',columns='control',values='rho');ax=panel(axs[1],'Supplementary_Figure_S2','b','Paired-change marker associations',assoc);heat(ax,hh,vmax=.6,label='Spearman rho');ax.set_yticklabels([short[v] for v in hh.index]);ax.set_xticklabels(hh.columns,rotation=50,ha='right');note(ax,'GSE153470; n=81')
gd=pd.read_csv(src/'Figure_3_i.csv.gz',index_col=0);om=[]
for gene in gd.index:om.append({'omitted_gene':gene,'mean_delta':gd.drop(gene).mean().mean()})
om=pd.DataFrame(om);ax=panel(axs[2],'Supplementary_Figure_S2','c','Gene omission',om);ax.scatter(om.mean_delta,range(len(om)),s=13,c=BLUE);ax.axvline(0,c=GREY,lw=.7);ax.set_yticks([]);ax.set_xlabel('Import mean change');ax.set_ylabel('Omit one measured import gene');ax.set_xlim(right=0)
loo=pd.read_csv(r/'analysis/leave_one_patient_out.csv');z=loo[(loo.accession=='GSE153470')&(loo.module=='Import_TIM_TOM')&(loo.scale=='source')];ax=panel(axs[3],'Supplementary_Figure_S2','d','Patient omission',z);ax.scatter(z.mean_delta,range(len(z)),s=5,c=TEAL);ax.axvline(0,c=GREY,lw=.7);ax.set_yticks([]);ax.set_ylabel('Omit one of 81 patients');ax.set_xlabel('Import mean change');ax.set_xlim(right=0)
sc=pd.read_csv(root/'revision3/analysis/overlap_sensitivity_key_results.csv');z=sc[(sc.module=='Import_TIM_TOM')&(sc.scale=='mean_log')];ax=panel(axs[4],'Supplementary_Figure_S2','e','Related letrozole-label exclusions',sc);ax.scatter(z.mean_delta,range(len(z)),s=24,c=[BLUE,BLUE,ORANGE]);ax.axvline(0,c=GREY,lw=.7);ax.set_yticks(range(len(z)),['50','32, 50','Broad']);ax.invert_yaxis();ax.set_xlabel('Import paired change')
for i,(_,v) in enumerate(z.iterrows()):ax.text(v.mean_delta,i+.19,f'n={int(v.n1)}; q={v.q:.3f}',size=5.2,ha='center')
ax.set_xlim(-.23,.01);ax.set_ylim(2.5,-.5)
pc=pd.read_csv(r/'analysis/GSE87411_pair_gene_coverage.csv');z=pc[pc.module.isin(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV'])];ax=panel(axs[5],'Supplementary_Figure_S2','f','Z1031 paired-observed gene coverage',z)
for i,module in enumerate(['Import_TIM_TOM','Mitophagy','OXPHOS_Complex_IV']):
 v=z[z.module==module];ax.scatter(i+np.linspace(-.12,.12,len(v)),v.n_paired_observed/v.n_original,s=6,c=pal[i],alpha=.35)
ax.axhline(.5,c=GREY,ls='--',lw=.7);ax.set_xticks([0,1,2],['TIM/TOM','Mitophagy','Complex IV'],rotation=35,ha='right');ax.set_ylabel('Fraction of original genes per pair');ax.set_ylim(0,1.05)
allgenes=pd.read_csv(r/'analysis/GSE153470_all_gene_paired_tests.csv');ax=panel(axs[6],'Supplementary_Figure_S2','g','Genome-wide paired expression',allgenes,'All 25,012 mapped complete genes; paired t tests with genome-wide BH');ax.scatter(allgenes.mean_delta,-np.log10(allgenes.q_genomewide.clip(lower=1e-30)),s=2,c='#CDD4D9',alpha=.5,rasterized=True);z=allgenes[allgenes.gene.isin(mods['Import_TIM_TOM'])];ax.scatter(z.mean_delta,-np.log10(z.q_genomewide.clip(lower=1e-30)),s=11,c=BLUE);ax.axhline(-np.log10(.05),ls='--',c=GREY,lw=.7);ax.axvline(0,c=GREY,lw=.6);ax.set_xlabel('Mean source-unit change');ax.set_ylabel('-log10(genome-wide q)');note(ax,'Blue: predefined import genes')
pp=pd.read_csv(r/'analysis/POETIC_between_arm_all_tests.csv');z=pp[pp.family=='mito'];hh=z.pivot(index='module',columns='scale',values='q');ax=panel(axs[7],'Supplementary_Figure_S2','h','POETIC: controlled support',z,'Two-sided Mann-Whitney tests of paired patient changes, BH over twelve eligible modules per scale');heat(ax,-np.log10(hh),cmap='Purples',vmax=2,label='-log10(q)');ax.set_yticklabels([short[v] for v in hh.index],size=4.8);ax.set_xticklabels(hh.columns,rotation=30,ha='right');note(ax,'157 AI / 56 untreated')
save(f,'Supplementary_Figure_S2')
pd.DataFrame(manifest).to_csv(out/'panel_source_data_manifest.csv',index=False)
print('GENERATED',len(manifest),'panels; 4 main x 10, 2 supplementary x 8',flush=True)
