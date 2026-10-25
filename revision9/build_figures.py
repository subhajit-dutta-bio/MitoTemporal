from common import *
import shutil
# Reuse the established figure grammar; do not execute prior analyses or copy stale numbering.
code=(root/'revision8/build_figures.py').read_text()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
dest=r/'deliverables';figdir=dest/'figures';src=dest/'source_data'
figdir.mkdir(parents=True,exist_ok=True);src.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';RED='#AF5464';INK='#20343F';GREY='#77868E';palette=[BLUE,TEAL,ORANGE,RED];rng=np.random.default_rng(20261024);manifest=[]
exec(code[code.index('def canvas('):code.index('e=pd.read_csv')],globals())
old=root/'revision8/deliverables';prior=pd.read_csv(old/'panel_source_data_manifest.csv')
for name in ['Figure_1','Figure_2','Figure_3']+[f'Supplementary_Figure_S{i}' for i in range(1,11)]:
 for ext in ['pdf','svg','png']:shutil.copy2(old/'figures'/f'{name}.{ext}',figdir/f'{name}.{ext}')
 for rec in prior[prior.figure==name].to_dict('records'):
  shutil.copy2(old/'source_data'/rec['source_data'],src/rec['source_data']);manifest.append(rec)
# Regenerate the clinical figure with its new visible number and matching source filenames.
newout=out;out=root/'revision8/analysis'
k=pd.read_csv(out/'Ki67_patient_scores.csv');kt=pd.read_csv(out/'Ki67_model_tests.csv');kd=pd.read_csv(out/'Ki67_delta_correlations.csv')
exec(code[code.index("name='Figure_4'"):code.index("name='Supplementary_Figure_S10'")].replace("name='Figure_4'","name='Figure_6'"),globals());out=newout

def dots(ax,z,groups,column='Import_TIM_TOM',labels=None):
 for i,g in enumerate(groups):
  v=z.loc[z.group==g,column].dropna();ax.scatter(i+rng.uniform(-.07,.07,len(v)),v,s=18,color=palette[i%4],alpha=.85);ax.plot([i-.17,i+.17],[v.mean()]*2,lw=1.3,color=INK)
 ax.set_xticks(range(len(groups)),labels or groups,rotation=17,ha='right');ax.set_ylabel('Mean log2(CPM + 1)')

er=pd.read_csv(out/'ER_loss_sample_scores.csv');ep=pd.read_csv(out/'ER_loss_primary_tests.csv');small=pd.read_csv(out/'ER_loss_small_n_context.csv');ps=pd.read_csv(out/'PDX_sample_scores.csv');pt=pd.read_csv(out/'PDX_exploratory_tests.csv')
name='Figure_4';f,axs=canvas(name,'ER loss and pharmacological suppression','Primary culture tests | small-n context | separately specified exploratory xenografts')
z=er.query('accession=="GSE152323" and scale=="source" and group in ["CTL","ESR1"]');panel(axs[0],name,'a','ESR1 knockdown | target engagement',z)
for i,c in enumerate(['ESR1','ER_targets_no_ESR1','Proliferation']):
 baseline=z.loc[z.group=='CTL',c].mean()
 for j,g in enumerate(['CTL','ESR1']):
  v=z.loc[z.group==g,c]-baseline;axs[0].scatter(i+(j-.5)*.22+rng.uniform(-.03,.03,len(v)),v,s=13,color=[GREY,BLUE][j],label=g if i==0 else None);axs[0].plot([i+(j-.5)*.22-.06,i+(j-.5)*.22+.06],[v.mean()]*2,c=INK,lw=.8)
axs[0].set_xticks(range(3),['ESR1 RNA','ER targets','Proliferation'],rotation=17);axs[0].set_ylabel('Score minus control mean');axs[0].legend(frameon=False,fontsize=6);zero(axs[0]);note(axs[0],'n=3 cultures/group; no pairing assumed')
panel(axs[1],name,'b','siESR1 | import remains uncertain',z);dots(axs[1],z,['CTL','ESR1'],labels=['Control','siESR1']);note(axs[1],'Welch source q=0.060; rank q=0.058')
z=er.query('accession=="GSE231397" and scale=="source"');panel(axs[2],name,'c','Fulvestrant | E2 factorial design',z);dots(axs[2],z,['Veh','E2','Fulv','Fulv+E2'],labels=['Vehicle','E2','Fulv','E2 + Fulv']);note(axs[2],'24 h; three cultures per condition')
z=ep.query('module=="Import_TIM_TOM" and scale=="source"');panel(axs[3],name,'d','Three fixed primary import contrasts',z);forest(axs[3],z,['siESR1 - control','E2 + Fulv - E2','E2 x Fulv']);axs[3].set_xlabel('Source-scale contrast; 95% interval');note(axs[3],'BH q: 0.060 / 0.023 / 0.060')
z=small.query('module=="Import_TIM_TOM" and scale=="source"');panel(axs[4],name,'e','Small-n context | descriptive only',z);yy=np.arange(len(z));axs[4].scatter(z.estimate,yy,s=17,c=TEAL);axs[4].axvline(0,c=GREY,lw=.6);axs[4].set_yticks(yy,['ER KO 4 d','Giredestrant 4 d','siESR1 48 h','siESR1 72 h','Giredestrant mock','GDC-0927 mock'],fontsize=5.8);axs[4].invert_yaxis();axs[4].set_xlabel('Treated minus source control');note(axs[4],'n=2/group; shared study and controls; no P')
for ax,letter,model in zip(axs[5:7],'fg',['HCI-003','HCI-011']):
 z=ps.query('scale=="source" and model==@model').copy();z['group']=z.dose.astype(str);panel(ax,name,letter,model+' | exploratory public PDX',z);dots(ax,z,['0','1','80'],labels=['Vehicle','1 mg/kg','80 mg/kg']);note(ax,'Four animals/group; five daily doses')
z=pt.query('module=="Import_TIM_TOM" and scale=="source"');panel(axs[7],name,'h','PDX | all four fixed dose contrasts',z);forest(axs[7],z,[v.model+' '+str(v.dose_mg_kg)+' mg/kg' for v in z.itertuples()],col=TEAL);axs[7].set_xlabel('Giredestrant minus vehicle');note(axs[7],'Exploratory; q=0.006-0.037; two models')
save(f,name)

cp=pd.read_csv(out/'CPTAC_sample_scores.csv');ct=pd.read_csv(out/'CPTAC_primary_tests.csv');cc=pd.read_csv(out/'CPTAC_matched_module_concordance.csv');loo=pd.read_csv(out/'CPTAC_leave_one_import_gene_out.csv')
name='Figure_5';f,axs=canvas(name,'Protein concordance does not establish ER coupling','CPTAC breast 2020 | 122 unique tumours | 81 ER-positive primary adjusted cases')
for ax,letter,scale in zip(axs[:2],'ab',['source','rank']):
 z=cc.query('scale==@scale and module=="Import_TIM_TOM"');panel(ax,name,letter,'Matched RNA-protein | '+scale,z);ax.scatter(z.RNA,z.protein,c=np.where(z.ER_positive,BLUE,ORANGE),s=12,alpha=.65,edgecolors='none');ax.set_xlabel('RNA score | same 21 import genes');ax.set_ylabel('Protein score | same 21 genes');v=ct.query('scale==@scale and module=="Import_TIM_TOM" and test=="matched_RNA_protein"').iloc[0];note(ax,f'n=122; rho={v.estimate:.3f}; q={v.q:.3g}')

def partial(z,xname):
 z=z.dropna(subset=['Import_protein',xname,'Proliferation_RNA','Nonimport_mito_protein','Stromal','Immune']).copy();cols=[xname,'Proliferation_RNA','Nonimport_mito_protein','Stromal','Immune'];x=z[cols].apply(lambda v:(v-v.mean())/v.std(ddof=1));y=(z.Import_protein-z.Import_protein.mean())/z.Import_protein.std(ddof=1);nuis=sm.add_constant(x.drop(columns=xname));z['x_partial']=sm.OLS(x[xname],nuis).fit().resid;z['y_partial']=sm.OLS(y,nuis).fit().resid;return z
for ax,letter,xname,test,title in zip(axs[2:4],'cd',['ER_programme_RNA','ESR1_protein'],['ER_RNA_adjusted_ERpos','ESR1_protein_adjusted_ERpos'],['ER-programme RNA | within ER+','ESR1 protein | within ER+']):
 z=partial(cp.query('scale=="source" and ER_positive'),xname);panel(ax,name,letter,title,z);v=ct.query('scale=="source" and module=="Import_TIM_TOM" and test==@test').iloc[0];ax.scatter(z.x_partial,z.y_partial,s=11,c=BLUE,alpha=.6,edgecolors='none');xx=np.array([z.x_partial.min(),z.x_partial.max()]);ax.plot(xx,v.estimate*xx,c=RED,lw=1);zero(ax);ax.set_xlabel('Residual predictor (standardized units)');ax.set_ylabel('Residual import protein');note(ax,f'n=81; beta={v.estimate:.3f}; q={v.q:.3f}')
z=ct.query('module=="Import_TIM_TOM" and test!="matched_RNA_protein"');panel(axs[4],name,'e','Adjusted ER slopes | scoring sensitivity',z);forest(axs[4],z,[{'source':'Source','rank':'Rank','gene_zscore_sensitivity':'Gene z'}[v.scale]+(' | ER RNA' if v.test.startswith('ER_RNA') else ' | ESR1 protein') for v in z.itertuples()],col=RED);axs[4].set_xlabel('Standardized beta; HC3 95% interval');note(axs[4],'Separate families; no predicted positive slope')
z=cp.query('scale=="source" and ER_positive');panel(axs[5],name,'f','Broader mitochondrial abundance context',z);axs[5].scatter(z.Nonimport_mito_protein,z.Import_protein,s=12,c=TEAL,alpha=.6,edgecolors='none');axs[5].set_xlabel('Non-import mitochondrial protein score');axs[5].set_ylabel('Import protein score');note(axs[5],'129 fixed comparator proteins; descriptive')
z=cp.query('scale=="source" and ER_known');panel(axs[6],name,'g','IHC ER status | protein distributions',z)
for i,v in enumerate([z.loc[~z.ER_positive,'Import_protein'],z.loc[z.ER_positive,'Import_protein']]):
 viol=axs[6].violinplot(v,positions=[i],widths=.7,showextrema=False);viol['bodies'][0].set_facecolor([ORANGE,BLUE][i]);viol['bodies'][0].set_alpha(.22);axs[6].scatter(i+rng.uniform(-.14,.14,len(v)),v,s=7,c=[ORANGE,BLUE][i],alpha=.5);axs[6].plot([i-.2,i+.2],[v.median()]*2,c=INK)
axs[6].set_xticks([0,1],['ER- (39)','ER+ (81)']);axs[6].set_ylabel('Import protein log2-ratio score');note(axs[6],'Two unknown ER statuses excluded; descriptive')
z=loo.copy();panel(axs[7],name,'h','Leave one import protein out',z);axs[7].scatter(z.ER_RNA_coefficient,np.arange(len(z)),s=12,c=TEAL);axs[7].axvline(0,c=GREY,lw=.7);axs[7].set_yticks([0,10,20],['1','11','21']);axs[7].invert_yaxis();axs[7].set_xlabel('Within-ER+ adjusted ER-RNA beta');axs[7].set_ylabel('Omitted-protein iteration');note(axs[7],'Range -0.140 to -0.086; not new patients')
save(f,name)

name='Supplementary_Figure_S11';f,axs=canvas(name,'ER-loss and proteomic robustness','All fixed primary results retained | sensitivity analyses do not replace primary tests')
z=ep.query('module=="Import_TIM_TOM" and scale=="rank"');panel(axs[0],name,'a','ER loss | rank-scale primary tests',z);forest(axs[0],z,['siESR1 - control','E2 + Fulv - E2','E2 x Fulv']);axs[0].set_xlabel('Mean percentile-rank contrast');note(axs[0],'All three corrected tests remain uncertain')
z=ep.query('module=="CIV_nuclear" and scale=="source"');panel(axs[1],name,'b','Nuclear Complex IV comparator',z);forest(axs[1],z,['siESR1 - control','E2 + Fulv - E2','E2 x Fulv'],col=ORANGE);axs[1].set_xlabel('Source-scale contrast');note(axs[1],'Separate three-test family; no specificity test')
z=small.query('module=="Import_TIM_TOM" and scale=="rank"');panel(axs[2],name,'c','Small-n contexts | rank direction',z);axs[2].scatter(z.estimate,np.arange(6),s=17,c=TEAL);axs[2].axvline(0,c=GREY,lw=.6);axs[2].set_yticks(range(6),['ER KO','Giredestrant 4 d','siESR1 48 h','siESR1 72 h','Giredestrant mock','GDC-0927 mock'],fontsize=5.8);axs[2].invert_yaxis();axs[2].set_xlabel('Rank-score mean difference');note(axs[2],'n=2/group; descriptive; no P values')
z=pt.query('module=="Import_TIM_TOM" and scale=="rank"');panel(axs[3],name,'d','Exploratory PDX | rank sensitivity',z);forest(axs[3],z,[v.model+' '+str(v.dose_mg_kg)+' mg/kg' for v in z.itertuples()],col=TEAL);axs[3].set_xlabel('Rank-score difference');note(axs[3],'Four-test family; shared vehicle per model')
z=pd.read_csv(out/'CPTAC_gene_coverage.csv').query('module=="Import_TIM_TOM"');panel(axs[4],name,'e','All 24 fixed import genes | coverage',z);heat=z.pivot(index='gene',columns='modality',values='n_measured').reindex(mods['Import_TIM_TOM']);im=axs[4].imshow(heat,aspect='auto',cmap='Blues',vmin=0,vmax=122);axs[4].set_xticks([0,1],heat.columns);axs[4].set_yticks(range(24),heat.index,fontsize=4.5);plt.colorbar(im,ax=axs[4],fraction=.04,pad=.04).set_label('Tumours measured',size=6);note(axs[4],'Reliability threshold: at least 80% of 122')
z=pd.read_csv(out/'CPTAC_import_gene_RNA_protein_correlations.csv');panel(axs[5],name,'f','Gene-level RNA-protein concordance',z);v=z.dropna(subset=['rho']).sort_values('rho');axs[5].scatter(v.rho,np.arange(len(v)),s=10,c=BLUE);axs[5].axvline(0,c=GREY,lw=.6);axs[5].set_yticks(range(len(v)),v.gene,fontsize=4.4);axs[5].set_xlabel('Spearman rho across tumours');note(axs[5],'Genes are components, not biological replicates')
z=ct.query('module=="CIV_nuclear" and scale=="source"');panel(axs[6],name,'g','CPTAC nuclear-IV full comparator family',z);forest(axs[6],z,['RNA-protein rho','Adjusted ER-RNA beta','Adjusted ESR1 beta'],col=ORANGE);axs[6].set_xlabel('Rho / standardized beta, as labelled');note(axs[6],'None corrected; not import-specificity evidence')
z=pd.read_csv(out/'CPTAC_context_sensitivity.csv').query('scale=="source"');z=z[z.test.isin(['complete_import_protein_core','all_tumor_ER_RNA_adjusted_context'])];panel(axs[7],name,'h','ER-RNA model contextual sensitivities',z);forest(axs[7],z,['Complete protein core (19)','All-tumour model'],col=RED);axs[7].set_xlabel('Adjusted standardized beta');note(axs[7],'Context only; unadjusted 95% HC3 intervals')
save(f,name)
pd.DataFrame(manifest).to_csv(dest/'panel_source_data_manifest.csv',index=False)
print('Revision 9:',len(manifest),'panels;',len(set(v['figure'] for v in manifest)),'figures')
