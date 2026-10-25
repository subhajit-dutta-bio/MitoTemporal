from common import *
import shutil,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
dest=r/'deliverables';figdir=dest/'figures';src=dest/'source_data'
figdir.mkdir(parents=True,exist_ok=True);src.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';RED='#AF5464';INK='#20343F';GREY='#77868E';palette=[BLUE,TEAL,ORANGE,RED];rng=np.random.default_rng(20261024);manifest=[]
old=root/'revision7/deliverables';prior=pd.read_csv(old/'panel_source_data_manifest.csv')
for original,new in [('Figure_1','Figure_1'),('Figure_2','Figure_2')]+[(f'Supplementary_Figure_S{i}',f'Supplementary_Figure_S{i}') for i in range(1,8)]+[('Figure_3','Supplementary_Figure_S8'),('Figure_4','Supplementary_Figure_S9')]:
 for ext in ['pdf','svg','png']:shutil.copy2(old/'figures'/f'{original}.{ext}',figdir/f'{new}.{ext}')
 for rec in prior[prior.figure==original].to_dict('records'):
  oldname=rec['source_data'];newname=oldname.replace(original,new,1);shutil.copy2(old/'source_data'/oldname,src/newname)
  manifest.append(dict(figure=new,panel=rec['panel'],title=rec['title'],source_data=newname))
def canvas(name,title,subtitle):
 f,ax=plt.subplots(4,2,figsize=(7.2,10));f.subplots_adjust(left=.15,right=.96,bottom=.065,top=.87,hspace=1.02,wspace=.68)
 f.text(.1,.966,name.replace('_',' ').upper(),fontsize=8,weight='bold',color=BLUE);f.text(.1,.937,title,fontsize=12,weight='bold',color=INK);f.text(.1,.911,subtitle,fontsize=7,color=GREY);return f,ax.ravel()
def panel(ax,name,letter,title,data):
 ax.set_title(title,loc='left',pad=13,weight='bold',color=INK);pos=ax.get_position();ax.figure.text(pos.x0-.048,pos.y1+.017,letter,fontsize=11,weight='bold')
 file=f'{name}_{letter}.csv.gz';data.to_csv(src/file,index=False,compression='gzip');manifest.append(dict(figure=name,panel=letter,title=title,source_data=file))
def note(ax,text):ax.text(0,1.015,text,transform=ax.transAxes,fontsize=6.2,color=GREY,va='bottom')
def zero(ax):ax.axhline(0,color='#C9D3DA',lw=.65,zorder=0)
def forest(ax,z,labels,x='estimate',ref=0,col=BLUE):
 yy=np.arange(len(z));ax.errorbar(z[x],yy,xerr=np.asarray([z[x]-z.low,z.high-z[x]]),fmt='o',color=col,ms=3,capsize=2,lw=.9);ax.axvline(ref,color=GREY,lw=.7);ax.set_yticks(yy,labels);ax.invert_yaxis()
def save(f,name):
 for ext in ['pdf','svg','png']:f.savefig(figdir/f'{name}.{ext}',dpi=400,bbox_inches='tight',pad_inches=.13)
 plt.close(f)
def trajectories(ax,data,modules,colors):
 for module,col in zip(modules,colors):
  z=data[data.module==module]
  for _,v in z.groupby('replicate'):ax.plot(v.hour,v.delta,color=col,lw=.65,alpha=.25)
  v=z.groupby('hour').delta.mean();ax.plot(v.index,v,'-o',color=col,lw=1.4,ms=2,label=module.replace('Import_TIM_TOM','TIM/TOM').replace('ER_programme','ER programme').replace('CIV_nuclear','Nuclear IV'))
 ax.set_xticks([0,2,6,12,24]);ax.set_xlabel('Hours after E2');ax.set_ylabel('Change in mean log2(TPM + 1)');zero(ax)
 if len(modules)>1:ax.legend(frameon=False,fontsize=6,loc='best')
def cultures(ax,z):
 for j,(treatment,col) in enumerate([('veh',BLUE),('E2',ORANGE)]):
  means=[]
  for i,gen in enumerate(['WT','Y537S','D538G']):
   vals=z[(z.genotype==gen)&(z.treatment==treatment)].Import_TIM_TOM
   ax.scatter(i+(j-.5)*.25+rng.uniform(-.04,.04,len(vals)),vals,color=col,s=12,alpha=.8);means.append(vals.mean())
  ax.plot(np.arange(3)+(j-.5)*.25,means,'-_',lw=1,label='Vehicle' if treatment=='veh' else 'E2 24 h',color=col,ms=9)
 ax.set_xticks(range(3),['WT','Y537S','D538G']);ax.set_ylabel('Mean log2(coverage CPM + 1)');ax.legend(frameon=False,fontsize=6)
e=pd.read_csv(out/'GSE78167_module_trajectories.csv');et=pd.read_csv(out/'GSE78167_tests.csv');es=pd.read_csv(out/'GSE89888_sample_scores.csv');est=pd.read_csv(out/'GSE89888_primary_tests.csv');k=pd.read_csv(out/'Ki67_patient_scores.csv');kt=pd.read_csv(out/'Ki67_model_tests.csv');kd=pd.read_csv(out/'Ki67_delta_correlations.csv')
name='Figure_3';f,axs=canvas(name,'Import transcription responds to ER activation','Fixed module | independent culture trajectories | ligand-independent ESR1 mutants')
z=e[(e.scale=='source')&(e.module=='Import_TIM_TOM')];panel(axs[0],name,'a','Replicated estradiol time course',z);trajectories(axs[0],z,['Import_TIM_TOM'],[BLUE]);note(axs[0],'Three biological trajectories; 23/24 genes')
z=et[(et.scale=='source')&(et.window=='late')&et.module.isin(['Import_TIM_TOM','CIV_nuclear','NRF1','TFAM'])];panel(axs[1],name,'b','Fixed late-window contrasts',z);forest(axs[1],z,['TIM/TOM','Nuclear IV','NRF1','TFAM']);axs[1].set_xlabel('Mean 6-24 h minus baseline');note(axs[1],'95% t intervals; n=3 trajectories')
z=e[(e.scale=='source')&e.module.isin(['ER_programme','Proliferation'])];panel(axs[2],name,'c','ER and proliferation context',z);trajectories(axs[2],z,['ER_programme','Proliferation'],[TEAL,ORANGE])
z=e[(e.scale=='source')&e.module.isin(['NRF1','TFAM'])];panel(axs[3],name,'d','Regulator trajectories',z);trajectories(axs[3],z,['NRF1','TFAM'],[TEAL,RED]);note(axs[3],'Late-window regulator q=0.158 for both')
for ax,letter,line in zip(axs[4:6],'ef',['MCF7','T47D']):
 z=es[(es.scale=='source')&(es.cell_line==line)];panel(ax,name,letter,line+' | culture-level import',z);cultures(ax,z);note(ax,'Four cultures/group; source clones pooled')
z=est[(est.scale=='source')&(est.module=='Import_TIM_TOM')&est.contrast.isin(['Y537S','D538G'])];panel(axs[6],name,'g','Ligand-independent import expression',z);forest(axs[6],z,[f'{v.cell_line} {v.contrast}' for v in z.itertuples()]);axs[6].set_xlabel('Mutant minus WT under vehicle');note(axs[6],'T47D Y537S uncertain; other q<5e-8')
z=est[(est.scale=='source')&(est.module=='Import_TIM_TOM')&est.contrast.str.endswith('_E2')];panel(axs[7],name,'h','Attenuated E2 response in mutants',z);forest(axs[7],z,[f'{v.cell_line} {v.contrast[:5]}' for v in z.itertuples()],col=TEAL);axs[7].set_xlabel('Mutant-by-E2 interaction');note(axs[7],'95% HC3 intervals; ten-test import family')
save(f,name)
name='Figure_4';f,axs=canvas(name,'Clinical response tests limit biomarker claims','Ki67 protein endpoints | fixed baseline model | patient is the statistical unit')
for ax,letter,scale in zip(axs[:2],'ab',['source','rank']):
 z=kt[(kt.scale==scale)&(kt.model=='primary_full')];panel(ax,name,letter,'Baseline import | '+scale+' scoring',z);forest(ax,z,['Z1031 (n=98)' if v=='GSE87411' else 'POETIC (n=144)' for v in z.accession]);ax.set_xlabel('Adjusted log2(Ki67 + 0.1) per SD import');note(ax,'Baseline Ki67, ER and proliferation adjusted')
for ax,letter,acc,label in zip(axs[2:4],'cd',['GSE87411','GSE105777'],['Z1031','POETIC']):
 z=k[(k.scale=='source')&(k.accession==acc)&(k.arm=='AI')].dropna(subset=['ki67_pre','ki67_on','Import_TIM_TOM_pre','ER_programme_pre','Proliferation_pre']);z=z.assign(ki67_log_ratio=np.log2((z.ki67_on+.1)/(z.ki67_pre+.1)))
 panel(ax,name,letter,label+' | within-patient changes',z);ax.scatter(z.Import_TIM_TOM_delta,z.ki67_log_ratio,s=10,c=BLUE if acc=='GSE87411' else TEAL,alpha=.55,edgecolors='none');zero(ax);ax.axvline(0,color='#CBD4D9',lw=.6);ax.set_xlabel('Paired import-score change');ax.set_ylabel('log2 on/baseline Ki67 ratio');v=kd[(kd.accession==acc)&(kd.scale=='source')].iloc[0];note(ax,f'n={len(z)}; Spearman rho={v.rho:.2f}; q={v.q:.3f}')
z=kt[kt.model.str.startswith('POETIC_')];panel(axs[4],name,'e','POETIC | treatment interaction',z);forest(axs[4],z,z.scale.tolist(),col=TEAL);axs[4].set_xlabel('Baseline import x AI interaction');note(axs[4],'195 mapped cases, including 51 untreated')
z=pd.read_csv(out/'Ki67_LOOCV_prediction.csv').query('scale=="source"');panel(axs[5],name,'f','Incremental prediction error',z)
for i,(acc,v) in enumerate(z.groupby('accession',sort=False)):
 vals=v.set_index('model').loc[['covariates_only','plus_import'],'LOOCV_RMSE'];axs[5].plot([0,1],vals,'-o',c=palette[i],ms=3,lw=1.2,label='Z1031' if acc=='GSE87411' else 'POETIC')
axs[5].set_xticks([0,1],['Covariates','+ import']);axs[5].set_xlim(-.15,1.15);axs[5].set_ylabel('Leave-one-out RMSE');axs[5].legend(frameon=False,fontsize=6);note(axs[5],'Lower is better; no tuned cutpoint')
z=pd.DataFrame([dict(cohort='Z1031',paired=109,complete=98),dict(cohort='POETIC AI',paired=157,complete=144),dict(cohort='POETIC control',paired=56,complete=51),dict(cohort='FAIMoS',paired=81,complete=0)])
panel(axs[6],name,'g','Clinical endpoint eligibility',z);xx=np.arange(len(z));axs[6].bar(xx,z.complete,color=BLUE,width=.6);axs[6].bar(xx,z.paired-z.complete,bottom=z.complete,color='#DCE3E8',width=.6);axs[6].set_xticks(xx,z.cohort,rotation=18,ha='right');axs[6].set_ylabel('Deposited paired patients');note(axs[6],'Blue: complete; grey: unavailable / incomplete')
z=pd.read_csv(root/'revision7/assessment/PFS_same_population_diagnostic.csv').query('population=="SET_complete"');panel(axs[7],name,'h','PFS | identical 47-patient population',z);forest(axs[7],z,['Unadjusted','+ proliferation','+ proliferation + SET'],x='HR',ref=1,col=ORANGE);axs[7].set_xlabel('Import hazard ratio per SD');note(axs[7],'Exploratory; adjusted family q=0.092')
save(f,name)
name='Supplementary_Figure_S10';f,axs=canvas(name,'Scoring sensitivity and programme context','Secondary and exploratory analyses | all planned primary tests retained')
z=e[(e.scale=='rank')&(e.module=='Import_TIM_TOM')];panel(axs[0],name,'a','E2 import | rank-scale sensitivity',z);trajectories(axs[0],z,['Import_TIM_TOM'],[BLUE]);axs[0].set_ylabel('Mean percentile-rank change');note(axs[0],'Late-window P=0.0504; n=3')
z=pd.read_csv(out/'GSE78167_replicate_window_changes.csv').query('scale=="source" and window=="late"');panel(axs[1],name,'b','Every biological late-window change',z)
for i,module in enumerate(['Import_TIM_TOM','CIV_nuclear','NRF1','TFAM']):
 v=z[z.module==module];axs[1].scatter(i+rng.uniform(-.06,.06,len(v)),v.delta,s=15,c=palette[i]);axs[1].plot([i-.15,i+.15],[v.delta.mean()]*2,c=INK,lw=1)
axs[1].set_xticks(range(4),['TIM/TOM','Nuclear IV','NRF1','TFAM'],rotation=20);axs[1].set_ylabel('Mean log2-expression change');zero(axs[1])
z=est[est.module=='Import_TIM_TOM'].pivot(index=['cell_line','contrast'],columns='scale',values=['estimate','q']).reset_index();z.columns=['_'.join(c).strip('_') if isinstance(c,tuple) else c for c in z.columns];panel(axs[2],name,'c','ESR1 effects across scoring scales',z)
axs[2].scatter(z.estimate_source,z.estimate_rank,c=[BLUE if v=='MCF7' else TEAL for v in z.cell_line],s=18);zero(axs[2]);axs[2].axvline(0,c=GREY,lw=.6);axs[2].set_xlabel('Source-scale contrast');axs[2].set_ylabel('Rank-scale contrast');note(axs[2],'Ten fixed contrasts; one point per contrast')
a=pd.read_csv(out/'GSE89888_log2_coverage_CPM.csv.gz',index_col=0);md=pd.read_csv(out/'GSE89888_sample_crosswalk.csv').set_index('external_id');rows=[]
for line in ['MCF7','T47D']:
 for contrast in ['E2','Y537S','D538G']:
  p=md[(md.cell_line==line)&(md.genotype==('WT' if contrast=='E2' else contrast))&(md.treatment==('E2' if contrast=='E2' else 'veh'))].index;n=md[(md.cell_line==line)&(md.genotype=='WT')&(md.treatment=='veh')].index
  for gene in mods['Import_TIM_TOM']:rows.append(dict(gene=gene,contrast=line+' '+contrast,delta=a.loc[gene,p].mean()-a.loc[gene,n].mean()))
z=pd.DataFrame(rows);panel(axs[3],name,'d','Gene-level expression coherence',z);heat=z.pivot(index='gene',columns='contrast',values='delta');im=axs[3].imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-.7,vmax=.7);axs[3].set_xticks(range(6),heat.columns,rotation=35,ha='right',fontsize=5);axs[3].set_yticks(range(24),heat.index,fontsize=3.9);plt.colorbar(im,ax=axs[3],fraction=.035,pad=.03).set_label('Mean log2 change',size=6)
for ax,letter,acc,label in zip(axs[4:6],'ef',['GSE87411','GSE105777'],['Z1031','POETIC']):
 z=k[(k.scale=='source')&(k.accession==acc)&(k.arm=='AI')].dropna(subset=['ki67_pre','Import_TIM_TOM_pre']);panel(ax,name,letter,label+' | baseline covariation',z);ax.scatter(z.Import_TIM_TOM_pre,np.log2(z.ki67_pre+.1),s=10,alpha=.5,c=BLUE if acc=='GSE87411' else TEAL);ax.set_xlabel('Baseline import score');ax.set_ylabel('Baseline log2(Ki67 + 0.1)');note(ax,f'n={len(z)}; descriptive, not response prediction')
z=kt[(kt.scale=='source')&kt.model.isin(['primary_full','baseline_Ki67_only','pseudocount_1_sensitivity'])];panel(axs[6],name,'g','Ki67 model sensitivity',z);forest(axs[6],z,[('Z1031' if v.accession=='GSE87411' else 'POETIC')+' '+{'primary_full':'full','baseline_Ki67_only':'Ki67 adj.','pseudocount_1_sensitivity':'+1 offset'}[v.model] for v in z.itertuples()]);axs[6].set_xlabel('Baseline import coefficient');note(axs[6],'HC3 intervals; all comparisons retained')
z=pd.read_csv(out/'Letrozole_control_correlations.csv').query('scale=="source"');panel(axs[7],name,'h','Import change versus fixed controls',z);heat=z.pivot(index='control',columns='accession',values='rho');im=axs[7].imshow(heat,aspect='auto',vmin=-.7,vmax=.7,cmap='RdBu_r');axs[7].set_yticks(range(len(heat)),heat.index,fontsize=6);axs[7].set_xticks([0,1],['Extended letro','Late letro'],rotation=15)
for i in range(len(heat)):
 for j in range(2):axs[7].text(j,i,f'{heat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=6)
note(axs[7],'Exploratory; full twelve-test BH family')
save(f,name)
pd.DataFrame(manifest).to_csv(dest/'panel_source_data_manifest.csv',index=False)
print('Figures complete:',len(manifest),'panels')
import subprocess
subprocess.run([sys.executable,str(r/'rebuild_inherited_context.py')],check=True)
