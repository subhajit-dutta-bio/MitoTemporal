from pathlib import Path
r=Path(__file__).parent;root=r.parent
for folder in ['deliverables/figures','deliverables/source_data']:(r/folder).mkdir(parents=True,exist_ok=True)
s=(root/'revision4/build_figures.py').read_text(encoding='utf-8')
s=s.replace("r=Path(__file__).parent;root=r.parent;", "r=Path(__file__).parent;root=r.parent;data_r=root/'revision4';")
s=s.replace("r/'analysis", "data_r/'analysis")
s=s.replace("def main_canvas(n,title,subtitle):\n", "def main_canvas(n,title,subtitle):\n n={1:3,2:4,3:1,4:2}[n]\n")
s=s.replace("def panel(ax,figure,letter,title,data,note=''):\n", "def panel(ax,figure,letter,title,data,note=''):\n figure={'Figure_1':'Figure_3','Figure_2':'Figure_4','Figure_3':'Figure_1','Figure_4':'Figure_2'}.get(figure,figure)\n")
s=s.replace("def save(fig,name):\n", "def save(fig,name):\n name={'Figure_1':'Figure_3','Figure_2':'Figure_4','Figure_3':'Figure_1','Figure_4':'Figure_2'}.get(name,name)\n")
s=s.replace("src/'Figure_3_i.csv.gz'", "src/'Figure_1_i.csv.gz'")
s=s.replace("'Figure_4':'Figure_2'}.get(figure,figure)","'Figure_4':'Figure_2','Supplementary_Figure_S1':'Supplementary_Figure_S4','Supplementary_Figure_S2':'Supplementary_Figure_S1','Supplementary_Figure_S3':'Supplementary_Figure_S2','Supplementary_Figure_S4':'Supplementary_Figure_S3'}.get(figure,figure)")
s=s.replace("'Figure_4':'Figure_2'}.get(name,name)","'Figure_4':'Figure_2','Supplementary_Figure_S1':'Supplementary_Figure_S4','Supplementary_Figure_S2':'Supplementary_Figure_S1'}.get(name,name)")
s=s.replace("def sup_canvas(n,title):\n", "def sup_canvas(n,title):\n n={1:4,2:1,'S1':'S4','S2':'S1'}[n]\n")
s=s.replace("sd=np.sqrt((aa.var(ddof=1)+bb.var(ddof=1))/2)","sd=np.sqrt(((len(aa)-1)*aa.var(ddof=1)+(len(bb)-1)*bb.var(ddof=1))/(len(aa)+len(bb)-2))")
s=s.replace("exr=ex[ex.accession.isin(['GSE67916','GSE125738','GSE118713'])].copy()", "exr=ez[ez.accession.isin(['GSE67916','GSE125738','GSE118713'])].copy()")
s=s.replace("ax.scatter(z.mean_delta,y,s=22", "ax.scatter(z.standardized_difference,y,s=22")
s=s.replace("ax.set_xlabel('Resistant-parental change')", "ax.set_xlabel('Within-model standardized difference',fontsize=5.5)")
s=s.replace("'Clinical validation and gene-level support'", "'Clinical treatment associations and gene-level checks'")
s=s.replace("Endocrine therapy lowers import-associated transcription", "Lower import-associated expression during aromatase inhibition")
s=s.replace("External blood cohorts do not confirm progression", "External blood profiles during scheduled treatment")
s=s.replace("'Perturbation context separates transcription from function'", "'Experimental context and clinical marker associations'")
s=s.replace("'New experimental module contrasts'", "'Experimental module differences'")
s=s.replace("'Clinical import vs '", "'Clinical import vs '")
oldline=next(l for l in s.splitlines() if l.startswith("st=pd.read_csv(root/'source/SUBMISSION_PACKAGE"))
newline="ccp=pat[pat.compartment=='blood'].groupby(['patient','timepoint']).n_cells.sum().unstack().reindex(columns=['Pre','Post','Prog']);ax=panel(axs[5],'Supplementary_Figure_S1','f','Blood cell coverage by patient',ccp.reset_index());heat(ax,np.log10(ccp+1),cmap='Blues',vmax=5,label='log10(cells+1)');ax.set_ylabel('Patient identifier');ax.set_xticklabels(ccp.columns,rotation=30);ax.tick_params(axis='y',labelsize=5)"
s=s.replace(oldline,newline)
exec(compile(s,str(__file__),'exec'))

# Additional post-primary annotation sensitivity: no reuse of old result labels.
tests=pd.read_csv(r/'analysis/canonical_annotation_all_tests.csv');delta=pd.read_csv(r/'analysis/canonical_annotation_patient_deltas.csv');coverage=pd.read_csv(r/'analysis/canonical_annotation_gene_coverage.csv');ocov=pd.read_csv(data_r/'analysis/gene_coverage.csv');oldcov=pd.read_csv(root/'revision2/analysis/gene_coverage.csv');oldcov=pd.concat([oldcov,ocov,pd.read_csv(data_r/'analysis/POETIC_gene_coverage.csv')],ignore_index=True)
f,axs=plt.subplots(3,2,figsize=(7.2,10.1));f.subplots_adjust(left=.16,right=.94,bottom=.09,top=.86,hspace=.95,wspace=.8);axs=axs.ravel()
f.text(.115,.975,'SUPPLEMENTARY FIGURE S2',size=8,weight='bold',color=BLUE);f.text(.115,.944,'Post-primary gene-annotation sensitivity',size=11,weight='bold',color=INK);f.text(.115,.919,'HGNC canonical symbols | duplicate aliases removed | original patient populations retained',size=6.5,color=GREY)
for ax,module,letter in zip(axs[:2],['Import_TIM_TOM','OXPHOS_Complex_V'],['a','b']):
 z=coverage[coverage.module==module].copy();old=oldcov[oldcov.module==module].drop_duplicates('accession').set_index('accession');z['old_measured']=z.accession.map(old.n_mapped);panel(ax,'Supplementary_Figure_S3',letter,'Measured '+short[module]+' genes',z)
 for j,v in enumerate(z.itertuples()):ax.plot([v.old_measured,v.n_measured],[j,j],c=GREY,lw=1);ax.scatter(v.old_measured,j,c=GREY,s=17);ax.scatter(v.n_measured,j,c=BLUE,s=17)
 ax.set_yticks(range(len(z)),z.accession.str.replace('GSE',''),size=5.5);ax.set_xlabel('Gene count');ax.invert_yaxis();ax.set_xlim(0,25 if module=='Import_TIM_TOM' else 18);note(ax,'Grey: exact-symbol / blue: canonical')
for ax,module,letter in zip(axs[2:4],['Import_TIM_TOM','Mitophagy'],['c','d']):
 z=tests[tests.module==module].copy();z['label']=z.accession.str.replace('GSE','')+' '+z.contrast.map({'B_minus_A':'early','C_minus_A':'late','4_minus_1':'long','Tamoxifen_minus_Pretreatment':'Tam','on_minus_baseline':'AI','paired_observed_on_minus_baseline':'paired genes','AI_minus_control_patient_changes':'AI vs control'});hh=z.pivot(index='label',columns='scale',values='q');panel(ax,'Supplementary_Figure_S3',letter,short[module]+': canonical-test q values',z);heat(ax,-np.log10(hh),cmap='Purples',vmax=5,label='-log10(q)');ax.set_yticklabels(hh.index,size=5);ax.set_xticklabels(hh.columns,rotation=20);note(ax,'Correction families differ; see legend')
for ax,scale,letter in zip(axs[4:6],['source','rank'],['e','f']):
 z=delta[(delta.accession=='GSE105777')&(delta.module=='Import_TIM_TOM')&(delta.scale==scale)];panel(ax,'Supplementary_Figure_S3',letter,'POETIC canonical '+scale+' changes',z);strip(ax,{name:z.loc[z.contrast.str.startswith(name+'_'),'delta'] for name in ['AI','Control']},ylabel='TIM/TOM paired change');q=tests[(tests.accession=='GSE105777')&(tests.module=='Import_TIM_TOM')&(tests.scale==scale)].q.iloc[0];note(ax,f'n=157/56; q={q:.3f}')
for ext in ['pdf','svg','png']:f.savefig(figdir/f'Supplementary_Figure_S2.{ext}',dpi=360)
plt.close(f)

f,axs=plt.subplots(2,2,figsize=(7.2,7.8));f.subplots_adjust(left=.15,right=.94,bottom=.11,top=.84,hspace=.85,wspace=.65);axs=axs.ravel();f.text(.115,.97,'SUPPLEMENTARY FIGURE S3',size=8,weight='bold',color=BLUE);f.text(.115,.937,'Additional Complex IV perturbation context',size=11,weight='bold',color=INK)
sc=pd.read_csv(root/'revision2/analysis/GSE138459_module_scores.csv',index_col=0)
for ax,letter,groups,title in zip(axs[:3],list('abc'),[['polyT_ac','polyT_chr'],['OTI_ac','OTI_chr'],['veh_ac','veh_chr','nac_chr']],['Polyclonal CD8 stimulation','OT-I CD8 stimulation','Antioxidant context']):
 z=sc.loc[sc.index.str.startswith(tuple(v+'_' for v in groups)),['OXPHOS_Complex_IV']].reset_index();panel(ax,'Supplementary_Figure_S4',letter,title,z);labels={'polyT_ac':'Acute','polyT_chr':'Chronic','OTI_ac':'Acute','OTI_chr':'Chronic','veh_ac':'Acute','veh_chr':'Chronic','nac_chr':'Chronic + NAC'};strip(ax,{labels[g]:sc.loc[sc.index.str.startswith(g+'_'),'OXPHOS_Complex_IV'] for g in groups},ylabel='Complex IV mean log expression');note(ax,'3 samples/group; exact tests')
sc=pd.read_csv(root/'revision2/analysis/GSE296966_module_scores.csv',index_col=0);panel(axs[3],'Supplementary_Figure_S4','d','ER-Y537S MCF7: enzalutamide',sc.reset_index());strip(axs[3],{g:sc.loc[sc.index.str.startswith(g+'_'),'OXPHOS_Complex_IV'] for g in ['DMSO','Enza']},ylabel='Complex IV mean log expression');note(axs[3],'3 samples/group; no corrected support')
for ext in ['pdf','svg','png']:f.savefig(figdir/f'Supplementary_Figure_S3.{ext}',dpi=360)
plt.close(f)
pd.DataFrame(manifest).to_csv(out/'panel_source_data_manifest.csv',index=False)
print('Current figures:',len(manifest),'panels; 4 main x 10, supplementary 8/6/4/8',flush=True)
