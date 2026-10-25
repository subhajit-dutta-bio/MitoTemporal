from v14_common import *
import shutil,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
dest=r/'deliverables';src=dest/'source_data';figdir=dest/'figures';prior=root/'revision13/deliverables'
dest.mkdir(exist_ok=True);src.mkdir(exist_ok=True);figdir.mkdir(exist_ok=True)
for folder in ['figures','source_data']:
 for p in (prior/folder).iterdir():
  if p.is_file():shutil.copy2(p,dest/folder/p.name)
for p in prior.glob('Supplementary_Table_*.csv'):shutil.copy2(p,dest/p.name)
manifest=pd.read_csv(prior/'panel_source_data_manifest.csv').to_dict('records');provenance=pd.read_csv(prior/'panel_revision_map.csv').to_dict('records')
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';RED='#AF5464';INK='#20343F';GREY='#77868E';rng=np.random.default_rng(20261025)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def load(name):return pd.read_csv(out/(name+'.csv'))
def canvas():
 f,aa=plt.subplots(3,2,figsize=(7.2,10.1));f.subplots_adjust(left=.17,right=.96,bottom=.075,top=.95,hspace=.83,wspace=.78);return f,aa.ravel()
def panel(ax,n,letter,title,z,origin):
 ax.set_title(title,loc='left',pad=17,weight='bold',color=INK);pos=ax.get_position();ax.figure.text(pos.x0-.065,pos.y1+.022,letter,fontsize=11,weight='bold')
 name=f'{n}_{letter}.csv.gz';z.to_csv(src/name,index=False,compression='gzip');manifest.append(dict(figure=n,panel=letter,title=title,source_data=name));provenance.append(dict(figure=n,panel=letter,revision10_source=origin))
def note(ax,t):ax.text(0,1.025,t,transform=ax.transAxes,fontsize=5.7,color=GREY,va='bottom')
def save(f,n):
 for ext in ['pdf','svg','png']:f.savefig(figdir/f'{n}.{ext}',dpi=350,bbox_inches='tight',pad_inches=.13)
 plt.close(f)
def forest(ax,z,labels,xlabel,fontsize=6.3):
 y=np.arange(len(z));ax.errorbar(z.estimate,y,xerr=np.array([z.estimate-z.low,z.high-z.estimate]),fmt='o',c=BLUE,ms=3,capsize=2,lw=.85);ax.axvline(0,c=GREY,lw=.65);ax.set_yticks(y,labels,fontsize=fontsize);ax.invert_yaxis();ax.set_xlabel(xlabel)
def points(ax,x,y,xlabel,ylabel,color=BLUE):
 ax.scatter(x,y,s=12,c=color,alpha=.7,edgecolors='none');ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
def residual(a,y,x,cov):
 z=a[[x,y,*cov]].dropna();X=sm.add_constant(z[cov]);return z.assign(x_res=sm.OLS(z[x],X).fit().resid,y_res=sm.OLS(z[y],X).fit().resid)
def heat(ax,array,rows,cols,binary=False):
 im=ax.imshow(array,aspect='auto',cmap=ListedColormap(['#EDF1F3',BLUE]) if binary else 'viridis_r',vmin=0,vmax=1 if binary else 100);ax.set_yticks(range(len(rows)),rows,fontsize=5.7);ax.set_xticks(range(len(cols)),cols,rotation=50,ha='right',fontsize=5.7);return im

n='Supplementary_Figure_S20';f,axs=canvas();d=load('DepMap_breast_analysis_data');tests=load('DepMap_primary_tests');z=d.dropna(subset=['ER_status','import_minus_other_GE']);panel(axs[0],n,'a','Breast models | differential dependency',z,'DepMap_breast_analysis_data')
for j,label in enumerate([0,1]):
 v=z.loc[z.ER_status==label,'import_minus_other_GE'];axs[0].scatter(j+rng.uniform(-.09,.09,len(v)),v,c=[GREY,BLUE][j],s=17);axs[0].plot([j-.16,j+.16],[v.mean()]*2,c=INK,lw=1)
axs[0].set_xticks([0,1],['ER- (15)','ER+ (7)']);axs[0].set_ylabel('Import minus other-mito gene effect');note(axs[0],'More negative = greater relative dependency; D1 q=0.057')
z=residual(d,'import_minus_other_GE','ER_programme_RNA',['global_GE_median','Proliferation_RNA']);panel(axs[1],n,'b','Breast ER activity | adjusted association',z,'DepMap_breast_analysis_data');points(axs[1],z.x_res,z.y_res,'Residual ER-programme RNA','Residual differential gene effect');xx=np.array([z.x_res.min(),z.x_res.max()]);slope=np.linalg.lstsq(z.x_res.to_numpy()[:,None],z.y_res,rcond=None)[0][0];axs[1].plot(xx,slope*xx,c=RED,lw=1);note(axs[1],'n=45; positive slope is opposite to predicted selectivity')
z=tests[tests.test.isin(['D3','D4','D5','D6'])];panel(axs[2],n,'c','Pan-cancer | ESR1 co-dependency',z,'DepMap_primary_tests');forest(axs[2],z,['Import','TOMM20','TOMM40','TIMM23'],'Gene-effect coefficient / SD ESR1 gene effect');note(axs[2],'Lineage/global-effect adjusted; n=1,073-1,080')
z=tests[tests.test.isin(['D7','D8'])];panel(axs[3],n,'d','Pan-cancer | expression and dependency',z,'DepMap_primary_tests');forest(axs[3],z,['ESR1 RNA','ER-programme RNA'],'Differential gene effect / SD predictor');note(axs[3],'n=1,008; full eight-test family has no q<0.05')
z=load('DepMap_gene_specific_ER_tests');panel(axs[4],n,'e','Every fixed import gene | breast ER contrast',z,'DepMap_gene_specific_ER_tests');forest(axs[4],z,z.gene.tolist(),'ER+ minus ER- differential gene effect',fontsize=5.4);note(axs[4],'Gene-specific family; genes are not biological replicates')
panel(axs[5],n,'f','Import and other mitochondrial dependency',d,'DepMap_breast_analysis_data');colors=d.ER_status.map({0:GREY,1:BLUE}).fillna('#CFD7DC');points(axs[5],d.Nonimport_mito_GE,d.Import_TIM_TOM_GE,'Other-mito mean gene effect','Import mean gene effect',colors);note(axs[5],'48 models; dark blue ER+, grey ER-, pale unannotated');save(f,n)

n='Supplementary_Figure_S21';f,axs=canvas();g=pd.read_csv(out/'ER_cistrome_gene_annotation.csv.gz',index_col=0);imp=pd.read_csv(out/'ER_import_gene_annotation.csv',index_col=0);rows=[]
for endpoint in ['NRF1_promoter','E2_promoter','E2_promoter_or_distal']:
 for group,z in [('Import',g[g.is_import]),('Background',g[~g.is_import])]:rows.append(dict(endpoint=endpoint,group=group,n=len(z),bound=int(z[endpoint].sum()),fraction=z[endpoint].mean()))
z=pd.DataFrame(rows);panel(axs[0],n,'a','Promoter binding and linked ER regions',z,'ER_cistrome_gene_annotation')
for j,label in enumerate(['Background','Import']):
 a=z[z.group==label];axs[0].bar(np.arange(3)+(j-.5)*.27,a.fraction,width=.26,color=[GREY,BLUE][j],label=label)
axs[0].set_xticks(range(3),['NRF1\npromoter','ER E2\npromoter','ER E2\n+ distal']);axs[0].set_ylabel('Fraction of expressed genes');axs[0].legend(frameon=False,fontsize=5.7);note(axs[0],'23 import genes; distinct assays, not a binding competition')
cols=['NRF1_promoter','E2_promoter','E2_distal','Vehicle_distal'];panel(axs[1],n,'b','All expressed import genes | occupancy',imp.reset_index(),'ER_import_gene_annotation');heat(axs[1],imp[cols].to_numpy(int),imp.index,['NRF1 prom','ER prom E2','ER distal E2','ER distal veh'],True);note(axs[1],'Blue: replicate-supported; promoter window +/-1 kb')
null=pd.read_csv(out/'ER_expression_TSS_matched_null_counts.csv.gz');z=null[null.test.eq('C4')];panel(axs[2],n,'c','ER-linked import genes | matched null',z,'ER_expression_TSS_matched_null_counts');axs[2].hist(z['count'],bins=np.arange(z['count'].max()+2)-.5,color='#C6D9E4',edgecolor='white');axs[2].axvline(5,c=RED,lw=1.5,label='Observed: 5');axs[2].set_xlabel('ER promoter-or-distal links / 23 genes');axs[2].set_ylabel('Matched panels');axs[2].legend(frameon=False,fontsize=6);note(axs[2],'10,000 panels; enrichment P=0.998, seven-test q=1')
z=load('ER_cistrome_positive_control_genes');panel(axs[3],n,'d','Known ER programme | mapping controls',z,'ER_cistrome_positive_control_genes');heat(axs[3],z[['E2_promoter','E2_distal','NRF1_promoter']].to_numpy(int),z.gene,['ER promoter','ER distal','NRF1 promoter'],True);note(axs[3],'All seven genes have reproducible ER promoter/distal support')
z=load('ER_bound_import_early_response');panel(axs[4],n,'e','ER-linked genes | early induction difference',z,'ER_bound_import_early_response')
for rep,a in z.groupby('replicate'):axs[4].plot(a.hour,a.difference,'-o',c=BLUE,alpha=.65,ms=3,lw=.8)
axs[4].axhline(0,c=GREY,lw=.65);axs[4].set_xticks([1,2],['1 h','2 h']);axs[4].set_ylabel('ER-linked minus unlinked import response');note(axs[4],'5 linked / 18 unlinked genes; three trajectories; 1-h q=1')
z=load('ER_cistrome_fixed_sensitivity_counts');panel(axs[5],n,'f','ER promoter support | fixed sensitivities',z,'ER_cistrome_fixed_sensitivity_counts')
for j,condition in enumerate(['Vehicle','E2']):
 a=z[z.condition.eq(condition)];axs[5].plot(np.arange(4)+j*.06,a.import_promoter_bound,'-o',c=[GREY,BLUE][j],ms=4,label=condition)
axs[5].set_xticks(range(4),['1 kb\nboth','1 kb\neither','2 kb\nboth','2 kb\neither']);axs[5].set_ylabel('Bound import promoters / 23');axs[5].set_ylim(-.2,3);axs[5].legend(frameon=False,fontsize=6);note(axs[5],'Both replicates: 0 at either promoter width');save(f,n)

n='Supplementary_Figure_S22';f,axs=canvas();c=load('CCLE_protein_analysis_data');c=c[c.scale.eq('source')];a=load('Asleh_protein_analysis_data');a=a[a.scale.eq('source')]
panel(axs[0],n,'a','CCLE | matched import RNA and protein',c,'CCLE_protein_analysis_data');points(axs[0],c.import_RNA,c.matched_import_protein,'Import RNA (same 20 genes)','Matched import protein');note(axs[0],'28 breast models; rho=-0.090; three-slot q=1')
panel(axs[1],n,'b','CCLE | independently annotated ER groups',c,'CCLE_protein_analysis_data')
for j,(label,mask) in enumerate([('ER- (12)',c.ER_status.eq(0)),('ER+ (4)',c.ER_status.eq(1)),('Unknown (12)',c.ER_status.isna())]):
 v=c.loc[mask,'Import_TIM_TOM'];axs[1].scatter(j+rng.uniform(-.10,.10,len(v)),v,s=16,c=[GREY,BLUE,'#C6D1D8'][j]);axs[1].plot([j-.16,j+.16],[v.mean()]*2,c=INK,lw=1)
axs[1].set_xticks(range(3),['ER- (12)','ER+ (4)','Unknown\n(12)']);axs[1].set_ylabel('Import protein source score');note(axs[1],'ER contrast below n=5/group; ESR1 protein unmeasured')
panel(axs[2],n,'c','Asleh | complete protein subset by subtype',a,'Asleh_protein_analysis_data')
for j,label in enumerate(['Luminal A','Luminal B','Her2-enriched','Basal-like']):
 z=a[a['PAM50 Subtype'].str.lower().eq(label.lower())];v=z.Import_TIM_TOM;axs[2].scatter(j+rng.uniform(-.16,.16,len(v)),v,s=5,c=[BLUE,TEAL,ORANGE,GREY][j],alpha=.55);axs[2].plot([j-.23,j+.23],[v.mean()]*2,c=INK,lw=1)
axs[2].set_xticks(range(4),['Luminal A','Luminal B','HER2','Basal'],rotation=15);axs[2].set_ylabel('18-gene import protein score');note(axs[2],'284 public QC patients; observational subtype comparison')
z=load('Asleh_protein_primary_tests');z=z[z.scale.eq('source')&z.test.isin(['A1','A2'])];panel(axs[3],n,'d','Asleh | conditional primary protein tests',z,'Asleh_protein_primary_tests');forest(axs[3],z,['Luminal vs basal (209)','ER+ vs ER- (274)'],'Import protein coefficient; 95% HC3 interval',5.8);note(axs[3],'Adjusted other-mito protein and era; q=0.644 / 1')
panel(axs[4],n,'e','Asleh | shared mitochondrial abundance',a,'Asleh_protein_analysis_data');points(axs[4],a.Nonimport_mito,a.Import_TIM_TOM,'Other-mito protein source score','Import protein source score',a.era.map({0:BLUE,1:ORANGE}));note(axs[4],'Blue: 2008-13; orange: 1986-92; no public survival times')
z=load('CCLE_protein_dependency_matched_models');panel(axs[5],n,'f','CCLE | protein and dependency intersection',z,'CCLE_protein_dependency_matched_models');points(axs[5],z.Import_TIM_TOM,z.import_minus_other_GE,'Import protein source score','Import minus other-mito gene effect');note(axs[5],'22 models; descriptive rho=0.067; not independent replication');save(f,n)

summary=load('Pathway_import_ranking_summary');order=summary[summary.benchmark.eq('all')].context.tolist();primary=order[:-1];labels=['E2 late','MCF7 Y537S','MCF7 D538G','T47D Y537S','T47D D538G','Fulvestrant','Late letrozole','Extended letrozole','Anastrozole','Z1031','NeoPalAna','POETIC AI'];z=summary[summary.benchmark.eq('all')&summary.primary_context].set_index('context').loc[primary].reset_index()
n='Supplementary_Figure_S23';f,axs=canvas();panel(axs[0],n,'a','Import | all primary genome-wide ranks',z,'Pathway_import_ranking_summary');axs[0].scatter(z.rank_effect_top_percent,np.arange(12),s=22,c=[TEAL]*6+[BLUE]*6);axs[0].axvline(5,c=GREY,ls=':',lw=.8);axs[0].set_yticks(range(12),labels,fontsize=6);axs[0].invert_yaxis();axs[0].set_xlim(0,30);axs[0].set_xlabel('Directional top-tail percentile (lower = stronger)');note(axs[0],'5,085-5,405 eligible sets/context; 6/12 in strongest 5%')
ss=summary[summary.primary_context];panel(axs[1],n,'b','Membership and overlap sensitivities',ss,'Pathway_import_ranking_summary')
for b,col,label in [('size15to40',TEAL,'15-40 genes'),('exclude_half_import_overlap',ORANGE,'Remove high overlap')]:
 yy=ss[ss.benchmark.eq(b)].set_index('context').loc[primary].rank_effect_top_percent;axs[1].scatter(z.rank_effect_top_percent,yy,c=col,s=16,label=label)
axs[1].plot([0,40],[0,40],c=GREY,lw=.7);axs[1].set_xlabel('All-set import percentile');axs[1].set_ylabel('Restricted-background percentile');axs[1].legend(frameon=False,fontsize=5.6);note(axs[1],'Fixed sensitivity backgrounds; no best-rank selection')
panel(axs[2],n,'c','Ranking scale | source and gene ranks',z,'Pathway_import_ranking_summary');points(axs[2],z.rank_effect_top_percent,z.source_effect_top_percent,'Within-sample rank-effect percentile','Source-effect percentile');axs[2].plot([0,40],[0,40],c=GREY,lw=.7);note(axs[2],'Source effects ranked within each assay, never pooled')
panel(axs[3],n,'d','Import | complete module-wide corrections',z,'Pathway_import_ranking_summary');values=-np.log10(z[['source_q','rank_q']].to_numpy());im=axs[3].imshow(values,aspect='auto',cmap='Blues',vmin=0,vmax=10);axs[3].set_yticks(range(12),labels,fontsize=5.7);axs[3].set_xticks([0,1],['Source q','Rank q']);plt.colorbar(im,ax=axs[3],fraction=.05,pad=.04).set_label('-log10(q)',size=6);note(axs[3],'All eligible sets + import, separate family per context/scale')
dd=pd.read_csv(out/'POETIC_controlled_all_module_ranking.csv.gz');eligible=dd[dd.eligible];panel(axs[4],n,'e','POETIC | treated-minus-untreated ranking',eligible[['module','context','direction','eligible','rank_effect','source_effect','rank_p','rank_q','source_p','source_q']],'POETIC_controlled_all_module_ranking');v=-eligible.rank_effect;axs[4].hist(v,bins=70,color='#C6D9E4',edgecolor='white');value=-eligible.loc[eligible.module.eq('Import_TIM_TOM'),'rank_effect'].iloc[0];axs[4].axvline(value,c=RED,lw=1.4,label='Import');axs[4].set_xlabel('Directional controlled rank-score effect');axs[4].set_ylabel('Number of eligible modules');axs[4].legend(frameon=False,fontsize=6);note(axs[4],'157 AI + 56 untreated; import top-tail percentile 9.93')
chosen=['HALLMARK_ESTROGEN_RESPONSE_EARLY','HALLMARK_ESTROGEN_RESPONSE_LATE','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION'];allr=pd.read_csv(out/'Pathway_all_contexts_complete_results.csv.gz');rows=[]
for context in primary:
 q=allr[allr.context.eq(context)&allr.eligible];bg=q[~q.module.eq('Import_TIM_TOM')];direction=q.direction.iloc[0]
 for name in chosen:
  hit=q[q.aliases.str.split(';').map(lambda names:name in names)]
  if len(hit):
   row=hit.iloc[0];values=direction*bg.rank_effect;position=1+(values>direction*row.rank_effect).sum()+.5*(values==direction*row.rank_effect).sum();rows.append(dict(context=context,module=name,percentile=100*position/(len(bg)+1)))
rows=pd.DataFrame(rows);panel(axs[5],n,'f','Named programme context | complete contexts',rows,'retrospective fixed-name Hallmark display');table=rows.pivot(index='module',columns='context',values='percentile').reindex(index=chosen,columns=primary);im=heat(axs[5],table.to_numpy(),['ER early','ER late','E2F','G2M','MYC V1','OXPHOS'],['E2','M YS','M DG','T YS','T DG','Fulv','Let L','Let E','Ana','Z1031','Neo','POET']);plt.colorbar(im,ax=axs[5],fraction=.05,pad=.04).set_label('Top-tail percentile',size=6);note(axs[5],'Retrospective context display; complete screen retained');save(f,n)
pd.DataFrame(manifest).to_csv(dest/'panel_source_data_manifest.csv',index=False);pd.DataFrame(provenance).to_csv(dest/'panel_revision_map.csv',index=False)
for name in ['DepMap_primary_tests','DepMap_gene_specific_ER_tests','DepMap_prespecified_sensitivities','ER_cistrome_primary_tests','ER_cistrome_fixed_sensitivity_counts','ER_import_gene_annotation','CCLE_protein_primary_tests','Asleh_protein_primary_tests','Asleh_fixed_era_sensitivities','Independent_proteome_gene_coverage','Pathway_import_ranking_summary']:
 shutil.copy2(out/(name+'.csv'),dest/('Supplementary_Table_'+name+'.csv'))
print('Added four supplementary figures, six panels each;',len(manifest),'total source-linked panels; all four main figures unchanged')
