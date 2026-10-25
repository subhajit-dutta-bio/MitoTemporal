from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde,spearmanr
import textwrap,shutil
figdir=r/'deliverables/figures';src=r/'deliverables/source_data';figdir.mkdir(parents=True,exist_ok=True);src.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.7,'ytick.labelsize':6.7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';INK='#20343F';GREY='#7C8991';RED='#AF5464';pal=[BLUE,TEAL,ORANGE,RED];rng=np.random.default_rng(20261024);manifest=[]
def read(n):return pd.read_csv(out/(n+'.csv'),float_precision='round_trip')
def canvas(n,title,sub):
 f,axs=plt.subplots(4,2,figsize=(7.2,10.0));f.subplots_adjust(left=.14,right=.96,bottom=.07,top=.87,hspace=.95,wspace=.65)
 f.text(.1,.966,f'FIGURE {n}',color=BLUE,fontsize=8,weight='bold');f.text(.1,.937,title,fontsize=12,weight='bold',color=INK);f.text(.1,.911,sub,fontsize=7,color=GREY);return f,axs.ravel()
def panel(ax,n,letter,title,z):
 ax.set_title(title,loc='left',pad=13,weight='bold',color=INK);pos=ax.get_position();ax.figure.text(pos.x0-.045,pos.y1+.017,letter,fontsize=11,weight='bold');name=f'Figure_{n}_{letter}.csv.gz';z.to_csv(src/name,index=False,compression='gzip');manifest.append(dict(figure=f'Figure_{n}',panel=letter,title=title,source_data=name))
def note(ax,t):
 for old in list(ax.texts):
  if old.get_gid()=='panel_note':old.remove()
 v=ax.text(0,1.015,t,transform=ax.transAxes,fontsize=6.5,color=GREY,va='bottom');v.set_gid('panel_note')
def zero(ax):ax.axhline(0,color='#CAD3D8',lw=.7,zorder=0)
def slope(ax,z,color=BLUE,ylabel='Mean log-expression'):
 for row in z.itertuples():ax.plot([0,1],[row.pre,row.post],color=color,alpha=.22,lw=.65)
 ax.plot([0,1],[z.pre.mean(),z.post.mean()],'-o',c=INK,lw=1.5,ms=3);ax.set_xticks([0,1],['Baseline','On treatment']);ax.set_xlim(-.2,1.2);ax.set_ylabel(ylabel);note(ax,f'{len(z)} paired patients')
def rain(ax,groups,ylabel='Paired change',colors=None):
 for i,(lab,arr) in enumerate(groups.items()):
  a=np.asarray(arr,dtype=float);a=a[np.isfinite(a)];c=(colors or pal)[i%len(colors or pal)]
  if len(a)>2 and a.std()>0:
   yy=np.linspace(a.min(),a.max(),120);dd=gaussian_kde(a)(yy);dd=.25*dd/dd.max();ax.fill_betweenx(yy,i-.08-dd,i-.08,color=c,alpha=.23,lw=0)
  ax.scatter(i+.05+rng.uniform(0,.12,len(a)),a,s=8,c=c,alpha=.6,edgecolors='none');ax.plot([i-.08,i+.23],[np.mean(a)]*2,c=INK,lw=1.5)
 ax.set_xticks(range(len(groups)),list(groups));ax.set_ylabel(ylabel);ax.set_xlim(-.45,len(groups)-.5+.35);zero(ax)
def heat(ax,z,vmax=1,label=''):
 im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax);ax.set_xticks(range(z.shape[1]),z.columns,rotation=30,ha='right');ax.set_yticks(range(z.shape[0]),z.index);ax.tick_params(length=0);cb=plt.colorbar(im,ax=ax,fraction=.035,pad=.03);cb.ax.tick_params(labelsize=6);cb.set_label(label,size=6)
def forest(ax,z,labels,x='mean_delta',lo='low',hi='high',ref=0,color=BLUE):
 y=np.arange(len(z));ax.errorbar(z[x],y,xerr=np.array([z[x]-z[lo],z[hi]-z[x]]),fmt='o',color=color,ms=3,capsize=2,lw=.9);ax.axvline(ref,color=GREY,lw=.7);ax.set_yticks(y,labels);ax.invert_yaxis()
def save(f,n):
 for ext in ['pdf','svg','png']:f.savefig(figdir/f'Figure_{n}.{ext}',dpi=400,bbox_inches='tight',pad_inches=.13)
 plt.close(f)
def boot(v):
 v=np.asarray(v);b=rng.choice(v,(5000,len(v)),replace=True).mean(axis=1);return np.quantile(b,[.025,.975])
mp=read('mechanism_patient_scores');can=pd.read_csv(root/'revision6/analysis/canonical_annotation_patient_deltas.csv',float_precision='round_trip');ct=pd.read_csv(root/'revision6/analysis/canonical_annotation_all_tests.csv',float_precision='round_trip');neo=read('NeoPalAna_patient_deltas')
f,axs=canvas(1,'Recurring import-gene suppression during endocrine therapy','Patient pairing across studies | canonical annotation | untreated-biopsy context')
for ax,letter,acc,label in zip(axs[:4],'abcd',['GSE20181','GSE111563','GSE153470','GSE93204'],['Late letrozole','Extended letrozole','Two-week anastrozole','NeoPalAna AI run-in']):
 z=neo[(neo.module=='Import_TIM_TOM')&(neo.scale=='source')&(neo.window=='anastrozole_run_in')] if acc=='GSE93204' else mp[(mp.accession==acc)&(mp.module=='Import_TIM_TOM')&(mp.scale=='source')]
 panel(ax,1,letter,label+' | '+acc,z);slope(ax,z,ylabel='TIM/TOM score\n(source platform units)')
for ax,letter,scale in zip(axs[4:6],'ef',['source','rank']):
 z=can[(can.accession=='GSE105777')&(can.module=='Import_TIM_TOM')&(can.scale==scale)];panel(ax,1,letter,'POETIC | '+scale+' scoring',z);rain(ax,{'AI\n(n=157)':z[z.contrast.str.startswith('AI_')].delta,'Untreated\n(n=56)':z[z.contrast.str.startswith('Control_')].delta},ylabel='TIM/TOM paired change');q=ct[(ct.accession=='GSE105777')&(ct.module=='Import_TIM_TOM')&(ct.scale==scale)].q.iloc[0];note(ax,f'Mann-Whitney BH q={q:.3f}')
rows=[]
for acc,lab in [('GSE20181','Letrozole late'),('GSE111563','Letrozole long'),('GSE153470','Anastrozole'),('GSE87411','Z1031 paired genes')]:
 z=can[(can.accession==acc)&(can.module=='Import_TIM_TOM')&(can.scale=='rank')&(~can.contrast.str.startswith('B_'))];v=z.delta.to_numpy();lo,hi=boot(v);rows.append(dict(cohort=lab,n=len(v),mean_delta=v.mean(),low=lo,high=hi))
z=pd.DataFrame(rows);panel(axs[6],1,'g','Import changes on within-sample ranks',z);forest(axs[6],z,[f'{x.cohort} ({x.n})' if x.n!=109 else 'Z1031\npaired genes (109)' for x in z.itertuples()]);axs[6].set_xlabel('Mean rank change (95% bootstrap CI)')
a=oldmatrix('GSE20181').dropna(axis=0,how='any');s=scores(a);p=pd.read_csv(out/'GSE20181_all_identifier_audit.csv',dtype={'patient':str}).set_index('patient')[['A','B','C']].dropna();zz=[]
for pt,row in p.iterrows():
 for j,c in enumerate(['A','B','C']):zz.append(dict(patient=pt,time=j,score=s.loc[row[c],'Import_TIM_TOM']-s.loc[row.A,'Import_TIM_TOM']))
z=pd.DataFrame(zz);panel(axs[7],1,'h','Letrozole trajectory | complete triplets',z)
for _,v in z.groupby('patient'):axs[7].plot(v.time,v.score,c=BLUE,alpha=.12,lw=.6)
m=z.groupby('time').score.mean();axs[7].plot(m.index,m.values,'o-',c=INK,lw=1.5,ms=3);axs[7].set_xticks([0,1,2],['Pre','10-14 days','90 days']);axs[7].set_ylabel('TIM/TOM change from baseline');note(axs[7],'56 complete deposited trajectories');zero(axs[7]);save(f,1)

ft=read('FELINE_patient_time_scores');fd=read('FELINE_patient_deltas');tests=read('FELINE_all_tests');ps=read('FELINE_pseudobulk_tests');response=read('FELINE_response_tests')
f,axs=canvas(2,'Malignant-cell longitudinal data define the cellular limit','FELINE | source CNA-labelled cancer nuclei | 10x analysed separately from ICELL8')
z=ft[ft.platform=='10x'];hh=z.pivot(index='patient',columns='time',values='n_cells').reindex(columns=['S','M','E']);panel(axs[0],2,'a','Cancer-cell coverage by patient',z);im=axs[0].imshow(np.log10(hh+1),aspect='auto',cmap='Blues',vmin=0,vmax=4);axs[0].set_xticks([0,1,2],['Day 0','Day 14','Late']);axs[0].set_yticks([0,10,20,30],[hh.index[i] for i in [0,10,20,30]]);axs[0].set_ylabel('34 profiled 10x patients');plt.colorbar(im,ax=axs[0],fraction=.04).set_label('log10(nuclei + 1)',size=6)
for ax,letter,arm,title in zip(axs[1:4],'bcd',['A','B','C'],['Letrozole alone','Intermittent ribociclib + AI','Continuous ribociclib + AI']):
 z=fd[(fd.platform=='10x')&(fd.arm==arm)&(fd.module=='Import_TIM_TOM')&(fd.contrast=='E_minus_S')];panel(ax,2,letter,title,z);slope(ax,z,color=pal[['A','B','C'].index(arm)],ylabel='TIM/TOM mean ln(1 + CPM)');note(ax,f'{len(z)} late pairs; >=10 nuclei/time')
z=tests[tests.family=='primary_10x_letrozole'].copy();ci=[]
for row in z.itertuples():v=fd[(fd.platform=='10x')&(fd.arm=='A')&(fd.module==row.module)&(fd.contrast==row.contrast)].delta;ci.append(boot(v))
z['low']=[v[0] for v in ci];z['high']=[v[1] for v in ci];panel(axs[4],2,'e','Letrozole alone | primary tests',z);forest(axs[4],z,[('Import' if v.module=='Import_TIM_TOM' else 'Nuclear IV')+(' day 14' if v.contrast.startswith('M') else ' late') for v in z.itertuples()]);axs[4].set_xlabel('Mean change (95% bootstrap CI)');note(axs[4],'Four-test family: q=0.863-1.000')
z=fd[(fd.platform=='10x')&(fd.arm=='A')&(fd.module=='Import_TIM_TOM')&(fd.contrast=='E_minus_S')];panel(axs[5],2,'f','Source dynamic response | late import',z);rain(axs[5],{'Response\n(n=3)':z[z.response=='Response'].delta,'Non-response\n(n=4)':z[z.response=='Non-response'].delta},ylabel='TIM/TOM paired change');note(axs[5],'Between-response test: q=1.000')
z=ps[ps.arm=='A'].copy();panel(axs[6],2,'g','Pseudobulk sensitivity | AI alone',z);labels=[('Import' if v.module=='Import_TIM_TOM' else 'Nuclear IV')+(' d14' if v.contrast.startswith('M') else ' late') for v in z.itertuples()];axs[6].scatter(z.mean_delta,np.arange(len(z)),c=TEAL,s=23);axs[6].axvline(0,c=GREY,lw=.7);axs[6].set_yticks(range(len(z)),labels);axs[6].invert_yaxis();axs[6].set_xlabel('Mean log(1 + pseudobulk CPM) change');note(axs[6],'Separate four-test sensitivity: q>=0.168')
cov=read('FELINE_gene_coverage');panel(axs[7],2,'h','Measured genes and biopsy attrition',cov);x=cov[cov.platform=='10x'];axs[7].barh(['Import','Nuclear IV','Mitophagy'],x.n_measured,color=[BLUE,TEAL,ORANGE],height=.55);axs[7].set_xlabel('Measured canonical genes');axs[7].invert_yaxis();note(axs[7],'11 AI patients; 10 early / 7 late pairs');save(f,2)

en=read('NRF1_binding_enrichment');bind=read('NRF1_module_gene_binding');gro=read('GROseq_descriptive_timecourse');reg=read('mechanism_regulator_correlations');adj=read('mechanism_adjusted_tests');cox=read('GSE124647_Cox_tests')
f,axs=canvas(3,'Regulatory plausibility and tests of pathway specificity','Independent public regulatory assays | covariate adjustment | clinical outcomes')
panel(axs[0],3,'a','NRF1 promoter binding | MCF7',en);base=en.K.iloc[0]/en.N.iloc[0];axs[0].bar([0,1,2],[base,en.n_bound.iloc[0]/en.n_module_expressed.iloc[0],en.n_bound.iloc[1]/en.n_module_expressed.iloc[1]],color=[GREY,BLUE,ORANGE],width=.55);axs[0].set_xticks([0,1,2],['Expressed\nbackground','TIM/TOM\n14/24','Nuclear IV\n12/23']);axs[0].set_ylabel('Fraction with peak near a TSS');axs[0].set_ylim(0,.8);note(axs[0],'Import: matched-background q=0.043')
z=bind[bind.expressed].copy();panel(axs[1],3,'b','Binding across expression levels',z)
for j,(mod,c) in enumerate([('Import_TIM_TOM',BLUE),('CIV_nuclear',ORANGE)]):
 v=z[z.module==mod];axs[1].scatter(np.log10(v.baseline_GRO_value+1),v.NRF1_bound.astype(int)+rng.uniform(-.06,.06,len(v)),s=16,c=c,alpha=.7,label='Import' if j==0 else 'Nuclear IV')
axs[1].set_yticks([0,1],['No nearby peak','Nearby peak']);axs[1].set_xlabel('log10(baseline GRO value + 1)');axs[1].legend(frameon=False,fontsize=6,loc='center right');note(axs[1],'+/-1 kb of any GENCODE v47 TSS')
panel(axs[2],3,'c','Estradiol nascent-transcription profiles',gro)
for name,c in [('Import_TIM_TOM',BLUE),('NRF1',RED),('TFAM',TEAL)]:
 v=gro[gro.programme==name];axs[2].plot([0,10,40,160],v.mean_log2_response,'o-',c=c,ms=3,label='Import' if name=='Import_TIM_TOM' else name)
axs[2].set_xticks([0,40,160]);axs[2].set_xlabel('Estradiol exposure (min)');axs[2].set_ylabel('Mean log2((value+1)/(pre+1))');axs[2].legend(frameon=False,fontsize=6);zero(axs[2]);note(axs[2],'Deposited aggregate profiles: descriptive')
z=reg[(reg.scale=='source')&(reg.kind=='unadjusted')].copy();z['label']=z.accession.str.replace('GSE','')+' '+z.module.map({'Import_TIM_TOM':'Import','CIV_nuclear':'Nuc IV'});hh=z.pivot(index='label',columns='gene',values='rho').reindex(columns=['NRF1','TFAM','PPARGC1A']);panel(axs[3],3,'d','Patient-change regulator correlations',z);heat(axs[3],hh,vmax=1,label='Spearman rho');note(axs[3],'12 planned tests; unavailable TFAM blank')
z=adj[adj.scale=='source'];panel(axs[4],3,'e','Adjusted module-change intercepts',z);forest(axs[4],z,[v.accession.replace('GSE','')+' '+('Import' if v.module=='Import_TIM_TOM' else 'Nuc IV') for v in z.itertuples()],x='estimate');axs[4].set_xlabel('Conditional change (95% HC3 CI)');note(axs[4],'Proliferation + stromal/immune changes')
zs=[]
for acc in ['GSE20181','GSE111563']:
 a=oldmatrix(acc).dropna(axis=0,how='any');s=scores(a);p=pd.read_csv(out/f'{acc}_mechanism_pair_crosswalk.csv',dtype=str).set_index('patient');c0,c1=p.columns;d=s.loc[p[c1]].set_axis(p.index)-s.loc[p[c0]].set_axis(p.index);d['accession']=acc;d['patient']=d.index;zs.append(d)
z=pd.concat(zs);panel(axs[5],3,'f','Import and proliferation changes',z)
for acc,c in [('GSE20181',BLUE),('GSE111563',ORANGE)]:
 v=z[z.accession==acc];axs[5].scatter(v.Proliferation,v.Import_TIM_TOM,s=10,c=c,alpha=.6,label=acc.replace('GSE',''))
axs[5].set_xlabel('Proliferation change');axs[5].set_ylabel('TIM/TOM change');axs[5].legend(frameon=False,fontsize=6);zero(axs[5])
zz=[]
for acc in ['GSE20181','GSE111563']:
 es=read(acc+'_ESTIMATE_scores').set_index('Unnamed: 0');p=pd.read_csv(out/f'{acc}_mechanism_pair_crosswalk.csv',dtype=str).set_index('patient');c0,c1=p.columns;d=es.loc[p[c1]].set_axis(p.index)-es.loc[p[c0]].set_axis(p.index)
 for name in ['StromalSignature','ImmuneSignature']:zz.extend(dict(accession=acc,component=name,delta=v) for v in d[name])
z=pd.DataFrame(zz);panel(axs[6],3,'g','ESTIMATE component changes',z);rain(axs[6],{'20181\nstromal':z[(z.accession=='GSE20181')&(z.component=='StromalSignature')].delta,'20181\nimmune':z[(z.accession=='GSE20181')&(z.component=='ImmuneSignature')].delta,'111563\nstromal':z[(z.accession=='GSE111563')&(z.component=='StromalSignature')].delta,'111563\nimmune':z[(z.accession=='GSE111563')&(z.component=='ImmuneSignature')].delta},ylabel='Signature score change')
z=cox[(cox.endpoint=='PFS')&(cox.model.isin(['proliferation_adjusted','SET_RNAseq_complete_case']))];panel(axs[7],3,'h','Metastatic endocrine-treatment PFS',z);forest(axs[7],z,[('Import' if v.module=='Import_TIM_TOM' else 'Nuc IV')+(' | Prolif, n=97' if v.n==97 else ' | SET+Prolif, n=47') for v in z.itertuples()],x='hazard_ratio',ref=1);axs[7].set_xlabel('Hazard ratio per score SD (95% CI)');note(axs[7],'All two-module PFS families q>=0.092');save(f,3)

f,axs=canvas(4,'Mitophagy-associated transcription depends on context','Tamoxifen-treated tumours | other endocrine windows | resistant-culture contrasts')
z=mp[(mp.accession=='GSE147271')&(mp.module=='Mitophagy')&(mp.scale=='source')];panel(axs[0],4,'a','Tamoxifen | GSE147271',z);slope(axs[0],z,color=ORANGE,ylabel='Mitophagy-associated expression');note(axs[0],'28 pairs; canonical module 10/10 genes')
a=oldmatrix('GSE147271').dropna(axis=0,how='any');p=pd.read_csv(out/'GSE147271_mechanism_pair_crosswalk.csv',dtype=str).set_index('patient');c0,c1=p.columns;ds=a.loc[mods['Mitophagy'],p[c1]].set_axis(p.index,axis=1)-a.loc[mods['Mitophagy'],p[c0]].set_axis(p.index,axis=1);rows=[]
for g in ds.index:lo,hi=boot(ds.loc[g]);rows.append(dict(gene=g,mean_delta=ds.loc[g].mean(),low=lo,high=hi))
z=pd.DataFrame(rows);panel(axs[1],4,'b','Tamoxifen | all module genes',z);forest(axs[1],z,z.gene,color=ORANGE);axs[1].set_xlabel('Gene change (95% bootstrap CI)');axs[1].tick_params(axis='y',labelsize=6)
z=ct[ct.module=='Mitophagy'].copy();z['label']=z.accession.str.replace('GSE','')+' '+z.contrast.map({'B_minus_A':'early','C_minus_A':'late','4_minus_1':'long','Tamoxifen_minus_Pretreatment':'Tam','on_minus_baseline':'AI','paired_observed_on_minus_baseline':'AI','AI_minus_control_patient_changes':'AI-control'});hh=z.pivot(index='label',columns='scale',values='q');panel(axs[2],4,'c','Mitophagy | full clinical correction families',z);im=axs[2].imshow(-np.log10(hh),aspect='auto',cmap='Purples',vmin=0,vmax=3);axs[2].set_yticks(range(len(hh)),hh.index,fontsize=6);axs[2].set_xticks([0,1],hh.columns);plt.colorbar(im,ax=axs[2],fraction=.035).set_label('-log10(q)',size=6)
z=can[(can.module=='Mitophagy')&(can.scale=='rank')&(can.accession.isin(['GSE153470','GSE87411']))];panel(axs[3],4,'d','AI mitophagy changes are not uniformly null',z);rain(axs[3],{'Anastrozole\n(n=81)':z[z.accession=='GSE153470'].delta,'Z1031\n(n=109)':z[z.accession=='GSE87411'].delta},ylabel='Mitophagy mean-rank change');note(axs[3],'Z1031 uses paired-observed genes')
z=neo[(neo.module=='Mitophagy')&(neo.scale=='source')&(neo.window.isin(['anastrozole_run_in','added_palbociclib']))];panel(axs[4],4,'e','NeoPalAna | consecutive drug windows',z);rain(axs[4],{'AI run-in\n(n=22)':z[z.window=='anastrozole_run_in'].delta,'Added palbo\n(n=23)':z[z.window=='added_palbociclib'].delta},ylabel='Mitophagy source-score change');note(axs[4],'No corrected support in either window')
z=read('GSE75971_LTED_tests');hh=z.pivot(index='model',columns='module',values='standardized_difference').reindex(columns=['Import_TIM_TOM','CIV_nuclear','Mitophagy']);hh.columns=['Import','Nuclear IV','Mitophagy'];panel(axs[5],4,'f','Long-term estrogen-deprived cultures',z);heat(axs[5],hh,vmax=4,label='Difference / pooled within-group SD');note(axs[5],'3 vs 3 libraries; exact P>=0.10; colour clipped')
z=pd.concat([pd.read_csv(root/'revision6/deliverables/source_data'/f'Figure_2_{l}.csv.gz') for l in ['g','h']]);panel(axs[6],4,'g','Resistant derivatives | shared parental context',z)
labels=list(z.context.unique())
for mod,c,off in [('Import_TIM_TOM',BLUE,-.12),('Mitophagy',ORANGE,.12)]:
 v=z[z.module==mod].set_index('context').reindex(labels);axs[6].scatter(v.standardized_difference,np.arange(len(v))+off,c=c,s=17,label='Import' if c==BLUE else 'Mitophagy')
axs[6].set_yticks(range(len(labels)),labels,fontsize=5.8);axs[6].axvline(0,c=GREY,lw=.7);axs[6].invert_yaxis();axs[6].set_xlabel('Within-model standardized difference');axs[6].legend(frameon=False,fontsize=6)
z=pd.read_csv(root/'revision6/deliverables/source_data/Figure_2_b.csv.gz');panel(axs[7],4,'h','PGC1a-overexpressing CD8 TIL context',z)
for j,(mod,c,lab) in enumerate([('OXPHOS_Complex_IV',BLUE,'Complex IV'),('Mitophagy',ORANGE,'Mitophagy')]):
 ev=z[z['sample'].str.startswith(('A1_','C1_','D1_'))][mod];oe=z[z['sample'].str.startswith(('A2_','G1_','H1_'))][mod]
 for k,vals in enumerate([ev,oe]):axs[7].scatter(j*3+k+rng.uniform(-.1,.1,len(vals)),vals,s=17,c=c,alpha=.8);axs[7].plot([j*3+k-.2,j*3+k+.2],[vals.mean()]*2,c=INK,lw=1)
axs[7].set_xticks([0,1,3,4],['EV','PGC1a','EV','PGC1a']);axs[7].set_ylabel('Centered source module expression');note(axs[7],'Complex IV (blue) / mitophagy (orange); n=3');zero(axs[7]);save(f,4)
pd.DataFrame(manifest).to_csv(r/'deliverables/main_panel_source_manifest.csv',index=False);print('Created four main figures, eight panels each')
