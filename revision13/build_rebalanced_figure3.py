from pathlib import Path
import sys,json,shutil
r=Path(__file__).resolve().parent;root=r.parent;prior=root/'revision12';dest=r/'deliverables';src=dest/'source_data';figdir=dest/'figures'
sys.path.insert(0,str(root/'python_packages'))
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BLUE='#246A91';TEAL='#218C83';ORANGE='#D48237';INK='#20343F';GREY='#77868E';palette=[BLUE,TEAL,ORANGE,'#AF5464']
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.3,'ytick.labelsize':6.3,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
manifest=pd.read_csv(prior/'deliverables/panel_source_data_manifest.csv');manifest=manifest[manifest.figure!='Figure_3'].to_dict('records')
prov=pd.read_csv(prior/'deliverables/panel_revision_map.csv');prov=prov[prov.figure!='Figure_3'].to_dict('records')
def read(name):return pd.read_csv(prior/'deliverables/source_data'/f'{name}.csv.gz')
def panel(ax,figure,letter,title,z,origin):
 ax.set_title(title,loc='left',pad=13,weight='bold',color=INK);pos=ax.get_position();ax.figure.text(pos.x0-.048,pos.y1+.017,letter,fontsize=11,weight='bold')
 file=f'{figure}_{letter}.csv.gz';z.to_csv(src/file,index=False,compression='gzip');manifest.append(dict(figure=figure,panel=letter,title=title,source_data=file));prov.append(dict(figure=figure,panel=letter,revision10_source=origin))
def note(ax,t):ax.text(0,1.015,t,transform=ax.transAxes,fontsize=5.9,color=GREY,va='bottom')
def zero(ax):ax.axhline(0,c='#CAD3D8',lw=.65,zorder=0)
def save(f,name):
 for ext in ['pdf','svg','png']:f.savefig(figdir/f'{name}.{ext}',dpi=350,bbox_inches='tight',pad_inches=.13)
 plt.close(f)
n='Figure_3';f,aa=plt.subplots(4,2,figsize=(7.2,9.3));f.subplots_adjust(left=.15,right=.96,bottom=.07,top=.95,hspace=1.03,wspace=.68);axs=aa.ravel()
z=read('Figure_3_a');panel(axs[0],n,'a','Complete trajectories | timing display',z,'revision12_Figure_3_a')
for module,col,label in zip(['Import_TIM_TOM','ER_programme','Proliferation'],[BLUE,TEAL,ORANGE],['Import','ER programme','Proliferation']):
 v=z[z.module==module]
 for _,sub in v.groupby('replicate'):axs[0].plot(sub.hour,sub.display_scaled_delta,c=col,lw=.5,alpha=.18)
 mean=v.groupby('hour').display_scaled_delta.mean();axs[0].plot(mean.index,mean,'-o',c=col,ms=2,lw=1.2,label=label)
axs[0].set_xticks([0,2,4,8,12,24]);axs[0].set_xlabel('Hours after E2');axs[0].set_ylabel('Fraction of observed maximum');axs[0].legend(frameon=False,fontsize=5.6);zero(axs[0]);note(axs[0],'Descriptive scaling; not response-magnitude comparison')
z=read('Supplementary_Figure_S17_b');panel(axs[1],n,'b','E2-by-CHX | fixed import contrasts',z,'revision12_Supplementary_Figure_S17_b');y=np.arange(3);axs[1].errorbar(z.estimate,y,xerr=np.array([z.estimate-z.low,z.high-z.estimate]),fmt='o',c=BLUE,ms=3,capsize=2,lw=.85);axs[1].axvline(0,c=GREY,lw=.65);axs[1].set_yticks(y,['E2 without CHX','E2 during CHX','Interaction']);axs[1].invert_yaxis();axs[1].set_xlabel('RMA log2 contrast; 95% HC3 interval');note(axs[1],'Three-test q: 0.0018 / 0.444 / 0.101')
z=read('Figure_3_c');panel(axs[2],n,'c','Nascent import | all fixed E2 times',z,'revision12_Figure_3_c')
for line,col in [('MCF7',BLUE),('CAMA1',TEAL)]:
 v=z[z.line==line].sort_values('hour');axs[2].plot(range(4),v.estimate,'-o',c=col,ms=3,lw=1.1,label=line)
axs[2].set_xticks(range(4),['15 min','1 h','6 h','24 h']);axs[2].set_ylabel('E2 minus vehicle log2(CPM + 1)');axs[2].legend(frameon=False,fontsize=6);zero(axs[2]);note(axs[2],'Two source libraries/condition; descriptive')
z=read('Figure_3_d');panel(axs[3],n,'d','One-hour ER programme and import',z,'revision12_Figure_3_d')
for i,line in enumerate(['MCF7','CAMA1']):
 v=z[z.line==line].set_index('module').loc[['Import_TIM_TOM','ER_programme','Proliferation']];axs[3].scatter(np.arange(3)+(i-.5)*.18,v.estimate,c=palette[i],s=22,label=line)
axs[3].set_xticks(range(3),['Import','ER programme','Proliferation'],rotation=15);axs[3].set_ylabel('E2 minus vehicle source score');zero(axs[3]);note(axs[3],'Culture independence not verified; no biological P')
z=read('Figure_3_e');panel(axs[4],n,'e','NRF1 promoter binding | expressed genes',z,'revision12_Figure_3_e');axs[4].bar(range(3),[z.K.iloc[0]/z.N.iloc[0],z.n_bound.iloc[0]/24,z.n_bound.iloc[1]/23],color=[GREY,BLUE,ORANGE],width=.55);axs[4].set_xticks(range(3),['Background','Import\n14/24','Nuclear IV\n12/23']);axs[4].set_ylim(0,.8);axs[4].set_ylabel('Fraction with nearby NRF1 peak');note(axs[4],'Import expression-matched q=0.0428')
for ax,letter,title in zip(axs[5:7],'fg',['Late letrozole | patient-change coupling','Extended letrozole | patient-change coupling']):
 z=read('Figure_3_'+letter);panel(ax,n,letter,title,z,'revision12_Figure_3_'+letter);labels=z.control.str.replace('ER_programme','ER programme',regex=False)
 ax.barh(labels,z.rho,color=[BLUE if q<.05 else 'white' for q in z.q],edgecolor=BLUE,linewidth=.7,height=.55);ax.axvline(0,c=GREY,lw=.65);ax.invert_yaxis();ax.set_xlim(-.2,.8);ax.set_xlabel('Spearman rho | paired import changes');ax.tick_params(axis='y',labelsize=5.5);note(ax,'Filled: q<0.05; open: q>=0.05')
z=read('Supplementary_Figure_S17_a');panel(axs[7],n,'h','MCF7 | E2 and cycloheximide',z,'revision12_Supplementary_Figure_S17_a');rng=np.random.default_rng(20261025)
for i,g in enumerate(['Vehicle','E2','CHX','E2+CHX']):
 v=z.loc[z.group==g,'Import_TIM_TOM'];axs[7].scatter(i+rng.uniform(-.07,.07,len(v)),v,s=16,c=palette[i]);axs[7].plot([i-.15,i+.15],[v.mean()]*2,c=INK,lw=1)
axs[7].set_xticks(range(4),['Vehicle','E2','CHX','E2+CHX'],rotation=15);axs[7].set_ylabel('Mean deposited RMA log2');note(axs[7],'Interaction -0.090; source q=0.101');save(f,n)

# Enlarge the unchanged descriptive heatmap in a standalone supplementary panel.
n='Supplementary_Figure_S19';f,ax=plt.subplots(figsize=(7.2,7.6));f.subplots_adjust(left=.19,right=.88,bottom=.11,top=.93)
z=read('Figure_3_b');panel(ax,n,'a','Complete E2 time course | fixed import genes',z,'revision12_Figure_3_b');mods=json.loads((root/'revision6/analysis/canonical_modules.json').read_text())['Import_TIM_TOM'];hh=z.pivot_table(index='gene',columns='hour',values='delta',dropna=False).reindex(mods);im=ax.imshow(hh,cmap='RdBu_r',aspect='auto',vmin=-.6,vmax=.6);ax.set_yticks(range(len(hh)),hh.index,fontsize=8);ax.set_xticks(range(len(hh.columns)),hh.columns,fontsize=7);ax.set_xlabel('Hours after E2');plt.colorbar(im,ax=ax,fraction=.04,pad=.04).set_label('Mean log2(TPM + 1) change',size=7);save(f,n)
manifest=pd.DataFrame(manifest).sort_values(['figure','panel']);manifest.to_csv(dest/'panel_source_data_manifest.csv',index=False);pd.DataFrame(prov).sort_values(['figure','panel']).to_csv(dest/'panel_revision_map.csv',index=False)
print('Updated Figure 3; enlarged unchanged gene heatmap in S19; 173 source-linked panels')
