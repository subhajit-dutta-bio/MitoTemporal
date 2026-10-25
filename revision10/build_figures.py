from analysis_utils import *
import shutil,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
dest=r/'deliverables';figdir=dest/'figures';src=dest/'source_data';figdir.mkdir(parents=True,exist_ok=True);src.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';RED='#AF5464';INK='#20343F';GREY='#77868E';palette=[BLUE,TEAL,ORANGE,RED];rng=np.random.default_rng(20261025);manifest=[]
code=(root/'revision8/build_figures.py').read_text();exec(code[code.index('def canvas('):code.index('e=pd.read_csv')],globals())
old=root/'revision9/deliverables';prior=pd.read_csv(old/'panel_source_data_manifest.csv')
for name in prior.figure.unique():
 for ext in ['pdf','svg','png']:shutil.copy2(old/'figures'/f'{name}.{ext}',figdir/f'{name}.{ext}')
 for rec in prior[prior.figure==name].to_dict('records'):
  shutil.copy2(old/'source_data'/rec['source_data'],src/rec['source_data']);manifest.append(rec)
sp=pd.read_csv(out/'specificity_primary_tests.csv');bench=pd.read_csv(out/'matched_gene_benchmark_tests.csv');null=pd.read_csv(out/'matched_gene_null_distributions.csv.gz');diag=pd.read_csv(out/'matched_gene_matching_diagnostics.csv');cdk=pd.read_csv(out/'CDK_control_primary_tests.csv')
contexts=sp[sp.scale=='source'].context.drop_duplicates().tolist();labels=['Letrozole 56','Letrozole 37','POETIC AI - control','E2 late, n=3','MCF7 WT E2','T47D WT E2','siESR1 48 h','E2 + Fulv - E2','HCI003 1 mg/kg','HCI003 80 mg/kg','HCI011 1 mg/kg','HCI011 80 mg/kg'];comps=list(comparators)
def matrix(ax,z,value='estimate',limit=.5):
 heat=z.pivot(index='context',columns='comparator',values=value).reindex(index=contexts,columns=comps);im=ax.imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit);ax.set_yticks(range(len(heat)),labels,fontsize=5.1);ax.set_xticks(range(4),['Other mito','Ribosome','Translation','Cell cycle'],rotation=25,ha='right',fontsize=5.5)
 for i,key in enumerate(heat.index):
  for j,comp in enumerate(heat.columns):
   v=z[(z.context==key)&(z.comparator==comp)]
   if not v.empty and v.iloc[0].q<.05:ax.text(j,i,'•',ha='center',va='center',color='white' if abs(heat.iloc[i,j])>limit*.5 else INK,fontsize=6)
 plt.colorbar(im,ax=ax,fraction=.035,pad=.03).set_label('Import minus comparator change',size=6)
def distributions(ax,universe):
 for i,acc in enumerate(['GSE20181','GSE111563','GSE105777']):
  z=null.query('scale=="source" and universe==@universe and accession==@acc');vals=z.contrast.to_numpy();grid=np.linspace(vals.min()-.015,vals.max()+.015,150);density=gaussian_kde(vals)(grid);density=density/density.max()*.55;ax.fill_between(grid,i,i+density,alpha=.2,color=palette[i]);ax.plot(grid,i+density,c=palette[i],lw=.8);obs=bench.query('scale=="source" and universe==@universe and accession==@acc').iloc[0].observed_import_contrast;ax.plot([obs,obs],[i-.04,i+.63],c=palette[i],lw=1.2)
 ax.set_yticks([.2,1.2,2.2],['Letrozole 56','Letrozole 37','POETIC']);ax.set_xlabel('Clinical mean change / arm difference');ax.set_ylim(-.15,2.85)
name='Figure_7';f,axs=canvas(name,'Import distinctiveness depends on comparator and assay','Direct programme contrasts | baseline-matched gene sets | ER blockade versus CDK4/6 inhibition')
for ax,letter,scale,limit in zip(axs[:2],'ab',['source','rank'],[.5,.035]):
 z=sp.query('scale==@scale');panel(ax,name,letter,'Direct specificity | '+scale,z);matrix(ax,z,limit=limit);note(ax,'48-test family; dots mark q<0.05')
z=sp.query('scale=="source" and comparator=="Nonimport_mito"').set_index('context').loc[contexts].reset_index();panel(axs[2],name,'c','Import versus other mitochondrial genes',z);forest(axs[2],z,labels);axs[2].tick_params(axis='y',labelsize=5);axs[2].set_xlabel('Source-scale differential response');note(axs[2],'Clinical intervals cross zero; perturbations differ')
for ax,letter,universe,title in zip(axs[3:5],'de',['all_coding','nonimport_mito'],['Coding-gene benchmark','Mitochondrial-gene benchmark']):
 z=null.query('scale=="source" and universe==@universe').merge(bench.query('scale=="source" and universe==@universe')[['accession','observed_import_contrast','q','n_import']],on='accession');panel(ax,name,letter,title,z);distributions(ax,universe);note(ax,'1,000 selected sets/cohort; lines: import subset')
z=cdk.query('scale=="source" and module=="Import_TIM_TOM"');panel(axs[5],name,'f','OP-1250 versus CDK4/6 inhibitor',z);forest(axs[5],z,['Versus palbociclib (4/group)','Versus ribociclib (3/group)']);axs[5].set_xlabel('OP1 minus CDK25 import score');note(axs[5],'Primary q=0.066 for both; 28-day xenografts')
z=cdk.query('scale=="source" and module in ["Import_TIM_TOM","Proliferation","ER_programme"]');panel(axs[6],name,'g','ER and proliferation drug context',z);forest(axs[6],z,[v.experiment+' | '+{'Import_TIM_TOM':'Import','Proliferation':'Cell cycle','ER_programme':'ER targets'}[v.module] for v in z.itertuples()],col=TEAL);axs[6].set_xlabel('OP1 minus CDK25 programme score');note(axs[6],'P/R experiments separate; cytostasis not equivalent')
z=diag.copy();panel(axs[7],name,'h','Baseline coexpression matching quality',z)
for i,acc in enumerate(['GSE20181','GSE111563','GSE105777']):
 for universe,marker in [('all_coding','o'),('nonimport_mito','^')]:
  v=z[(z.accession==acc)&(z.universe==universe)].iloc[0];axs[7].scatter(v.target_baseline_correlation,v.median_selected_correlation,c=palette[i],marker=marker,s=22)
axs[7].plot([0,.21],[0,.21],c=GREY,lw=.7);axs[7].set_xlabel('Import mean baseline correlation');axs[7].set_ylabel('Selected sets: median correlation');note(axs[7],'Circles: coding; triangles: mitochondrial; imperfect')
save(f,name)

pro=pd.read_csv(out/'PROseq_all_tests.csv');pros=pd.read_csv(out/'PROseq_sample_scores.csv');prsp=pd.read_csv(out/'PROseq_specificity_tests.csv');pc=pd.read_csv(out/'perturbational_proteomics_coverage.csv');pg=pd.read_csv(out/'proteomics_import_gene_ratios.csv');ps=pd.read_csv(out/'proteomics_detected_subset_context.csv')
name='Figure_8';f,axs=canvas(name,'Early transcription and proteins limit a uniform mechanism','Prespecified one-hour window | relative PRO-seq | coverage-gated perturbational proteomics')
z=pro.query('scale=="source" and module=="Import_TIM_TOM" and contrast=="E2_minus_Vehicle"');panel(axs[0],name,'a','E2 import | all fixed times',z)
for line,col in [('MCF7',BLUE),('CAMA1',TEAL)]:
 v=z[z.line==line].sort_values('hour');axs[0].plot(np.arange(4),v.estimate,'-o',c=col,lw=1.2,ms=3,label=line)
axs[0].set_xticks(range(4),['15 min','1 h','6 h','24 h']);axs[0].set_ylabel('Mean log2(CPM + 1) contrast');zero(axs[0]);axs[0].legend(frameon=False,fontsize=6);note(axs[0],'Two source libraries/condition; descriptive')
z=pro.query('scale=="source" and module=="Import_TIM_TOM" and hour==1');panel(axs[1],name,'b','One-hour import | fixed comparisons',z);forest(axs[1],z,[v.line+' '+v.contrast.split('_')[0] for v in z.itertuples()]);axs[1].set_xlabel('Treated minus time-matched vehicle');note(axs[1],'Conditional Welch intervals; independence unverified')
z=prsp.query('scale=="source" and hour==1 and contrast=="E2_minus_Vehicle"');panel(axs[2],name,'c','One-hour E2 | comparator context',z);heat=z.pivot(index='comparator',columns='line',values='estimate').reindex(index=comps,columns=['MCF7','CAMA1']);im=axs[2].imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-.3,vmax=.3);axs[2].set_xticks([0,1],['MCF7','CAMA1']);axs[2].set_yticks(range(4),['Other mito','Ribosome','Translation','Cell cycle']);plt.colorbar(im,ax=axs[2],fraction=.04,pad=.03).set_label('Import minus comparator',size=6);note(axs[2],'Source-library means; no biological P values')
z=pro.query('scale=="source" and hour==1 and contrast=="E2_minus_Vehicle" and module in ["Import_TIM_TOM","ER_programme","Proliferation"]');panel(axs[3],name,'d','ER targets respond more strongly at 1 h',z)
for i,line in enumerate(['MCF7','CAMA1']):
 v=z[z.line==line].set_index('module').loc[['Import_TIM_TOM','ER_programme','Proliferation']];axs[3].scatter(np.arange(3)+(i-.5)*.18,v.estimate,c=palette[i],s=21,label=line)
axs[3].set_xticks(range(3),['Import','ER targets','Cell cycle'],rotation=15);axs[3].set_ylabel('Source-scale E2 minus vehicle');zero(axs[3]);axs[3].legend(frameon=False,fontsize=6);note(axs[3],'Different programmes; no temporal-precedence claim')
z=pg.copy();panel(axs[4],name,'e','Fixed import-panel protein coverage',z);z['measured']=z.log2_E2_vehicle.notna().astype(float);heat=z.pivot_table(index='gene',columns=['study','condition'],values='measured',aggfunc='mean').reindex(index=mods['Import_TIM_TOM'],columns=[('PXD002797','6h'),('PXD002797','36h'),('PXD049343','Cyto'),('PXD049343','Nuclear'),('PXD049343','membrane')]);im=axs[4].imshow(heat,aspect='auto',cmap='Blues',vmin=0,vmax=1);axs[4].set_yticks(range(24),heat.index,fontsize=4.3);axs[4].set_xticks(range(5),['6 h','36 h','Cyto','Nuclear','Membrane'],rotation=30,ha='right',fontsize=5);plt.colorbar(im,ax=axs[4],fraction=.04,pad=.03).set_label('Fraction of source ratios measured',size=6);note(axs[4],'Replicated assay: 11/24 common; required >=12')
z=pg.query('study=="PXD002797" and eligible_common_core');panel(axs[5],name,'f','Detected proteins | time-dependent direction',z);heat=z.pivot_table(index='gene',columns='condition',values='log2_E2_vehicle',aggfunc='mean').reindex(columns=['6h','36h']);im=axs[5].imshow(heat,aspect='auto',cmap='RdBu_r',vmin=-.4,vmax=.4);axs[5].set_yticks(range(len(heat)),heat.index,fontsize=5);axs[5].set_xticks([0,1],['6 h','36 h']);plt.colorbar(im,ax=axs[5],fraction=.04,pad=.03).set_label('Mean log2(E2 / vehicle)',size=6);note(axs[5],'Eleven measured components; not a whole-module test')
z=ps.copy();panel(axs[6],name,'g','Observable 11-protein subset | no P',z)
for i,hr in enumerate([6,36]):
 v=z[z.hour==hr];axs[6].scatter(i+rng.uniform(-.06,.06,len(v)),v.detected_subset_mean,s=18,c=palette[i]);axs[6].plot([i-.16,i+.16],[v.detected_subset_mean.mean()]*2,c=INK,lw=1.1)
axs[6].set_xticks([0,1],['6 h','36 h']);axs[6].set_ylabel('Detected-subset mean log2 ratio');zero(axs[6]);note(axs[6],'Three biological ratios/time; below coverage gate')
z=pg.query('study=="PXD049343" and condition=="membrane"').dropna(subset=['log2_E2_vehicle']).sort_values('log2_E2_vehicle');panel(axs[7],name,'h','Independent membrane/organelle summary',z);axs[7].scatter(z.log2_E2_vehicle,np.arange(len(z)),s=13,c=ORANGE);axs[7].axvline(0,c=GREY,lw=.7);axs[7].set_yticks(range(len(z)),z.gene,fontsize=4.8);axs[7].set_xlabel('log2(E2 / vehicle), 18 h');note(axs[7],'14 genes; one fraction summary; no biological n')
save(f,name)

name='Supplementary_Figure_S12';f,axs=canvas(name,'Specificity, source design and availability audit','Complete planned contrasts | low replication and missing coverage retained explicitly')
z=sp.pivot(index=['context','comparator'],columns='scale',values=['estimate','q']).reset_index();z.columns=['_'.join(c).strip('_') for c in z.columns];panel(axs[0],name,'a','Direct specificity across scales',z);axs[0].scatter(z.estimate_source,z.estimate_rank,c=BLUE,s=11,alpha=.6);zero(axs[0]);axs[0].axvline(0,c=GREY,lw=.7);axs[0].set_xlabel('Source differential response');axs[0].set_ylabel('Rank differential response');note(axs[0],'48 contrasts; significance is scale-dependent')
z=bench.query('scale=="source" and universe=="nonimport_mito"');panel(axs[1],name,'b','Baseline-eligible import genes for matching',z);axs[1].bar(range(3),z.n_import,color=[BLUE,TEAL,ORANGE],width=.6);axs[1].axhline(24,c=GREY,lw=.7);axs[1].set_ylim(0,27);axs[1].set_xticks(range(3),['Late letro','Extended letro','POETIC'],rotation=15);axs[1].set_ylabel('Genes / 24 fixed panel');note(axs[1],'Observed and null subsets have identical sizes')
z=diag.copy();panel(axs[2],name,'c','Baseline matching residuals',z);axs[2].scatter(z.median_expression_SD_distance,z.median_absolute_correlation_mismatch,c=[BLUE if v=='all_coding' else TEAL for v in z.universe],s=20);axs[2].set_xlabel('Median baseline mean/SD distance');axs[2].set_ylabel('Median correlation mismatch');note(axs[2],'Coding: blue; mitochondrial: teal')
z=pro.query('scale=="source" and hour==1 and module=="Import_TIM_TOM"');panel(axs[3],name,'d','One-hour conditional inference sensitivity',z);axs[3].scatter(z.conditional_Welch_q,z.conditional_permutation_q,c=BLUE,s=19);axs[3].axvline(.05,c=GREY,lw=.7);axs[3].axhline(.05,c=GREY,lw=.7);axs[3].set_xlim(0,1.03);axs[3].set_ylim(0,1.03);axs[3].set_xlabel('Conditional Welch q');axs[3].set_ylabel('Conditional exact-permutation q');note(axs[3],'Not confirmed biological tests; n=2 libraries/group')
for ax,letter,line in zip(axs[4:6],'ef',['MCF7','CAMA1']):
 z=pro.query('scale=="source" and module=="Import_TIM_TOM" and line==@line');panel(ax,name,letter,line+' | all fixed drug time courses',z)
 for treatment,col in [('E2_minus_Vehicle',BLUE),('Fulv_minus_Vehicle',ORANGE),('OP1250_minus_Vehicle',TEAL)]:
  v=z[z.contrast==treatment].sort_values('hour');ax.plot(range(4),v.estimate,'-o',lw=1,ms=2,c=col,label=treatment.split('_')[0])
 ax.set_xticks(range(4),['15 min','1 h','6 h','24 h']);ax.set_ylabel('Mean import source contrast');zero(ax);ax.legend(frameon=False,fontsize=5.5)
z=pc.query('study=="PXD002797" and condition=="36h" and module in ["Import_TIM_TOM","Nonimport_mito","Cytosolic_ribosome","Translation_factors","Proliferation"]');panel(axs[6],name,'g','Replicated proteome | all control coverage',z);axs[6].bar(range(len(z)),z.n_measured_all_replicates/z.n_fixed,color=TEAL,width=.6);axs[6].axhline(.5,c=GREY,lw=.7);axs[6].set_xticks(range(len(z)),['Import','Other mito','Ribosome','Translation','Cell cycle'],rotation=25,ha='right',fontsize=5);axs[6].set_ylabel('Measured gene fraction');note(axs[6],'Coverage gate retained for every programme')
z=pd.read_csv(out/'perturbational_proteomics_input_audit.csv');panel(axs[7],name,'h','Complete proteome-group input audit',z);axs[7].bar(range(5),z.retained_rows,color=BLUE,width=.65);axs[7].bar(range(5),z.source_rows-z.retained_rows,bottom=z.retained_rows,color='#DCE3E8',width=.65);axs[7].set_xticks(range(5),['6 h','36 h','Cyto','Nuclear','Membrane'],rotation=25,ha='right',fontsize=5);axs[7].set_ylabel('Deposited protein groups');note(axs[7],'Grey: flagged, unmapped or ambiguous groups')
save(f,name)
pd.DataFrame(manifest).to_csv(dest/'panel_source_data_manifest.csv',index=False);print('Revision 10:',len(manifest),'source-linked panels,',len(set(v['figure'] for v in manifest)),'figures')
