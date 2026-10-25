from pathlib import Path
import shutil
r=Path(__file__).parent;root=r.parent
# Reproduce only the inherited context figures, relabelled as supplementary.
s=(root/'revision4/build_figures.py').read_text(encoding='utf-8');s=s.replace("r=Path(__file__).parent;root=r.parent;", "r=Path(__file__).parent;root=r.parent;data_r=root/'revision4';")
s=s.replace("r/'analysis", "data_r/'analysis")
s=s.replace("sd=np.sqrt((aa.var(ddof=1)+bb.var(ddof=1))/2)","sd=np.sqrt(((len(aa)-1)*aa.var(ddof=1)+(len(bb)-1)*bb.var(ddof=1))/(len(aa)+len(bb)-2))")
s=s.replace("exr=ex[ex.accession.isin(['GSE67916','GSE125738','GSE118713'])].copy()", "exr=ez[ez.accession.isin(['GSE67916','GSE125738','GSE118713'])].copy()")
s=s.replace("ax.scatter(z.mean_delta,y,s=22","ax.scatter(z.standardized_difference,y,s=22").replace("ax.set_xlabel('Resistant-parental change')","ax.set_xlabel('Within-model standardized difference',fontsize=5.5)")
start=s.index('# Figure 3:');end=s.index('# Figure 4:');s=s[:start]+s[end:];s=s[:s.index('# Supplementary')]
s=s.replace("n,title,subtitle):\n", "n,title,subtitle):\n n={1:5,2:6,4:7}[n]\n")
s=s.replace("f'FIGURE {n}'", "f'SUPPLEMENTARY FIGURE S{n}'")
s=s.replace("def panel(ax,figure,letter,title,data,note=''):\n", "def panel(ax,figure,letter,title,data,note=''):\n figure={'Figure_1':'Supplementary_Figure_S5','Figure_2':'Supplementary_Figure_S6','Figure_4':'Supplementary_Figure_S7'}.get(figure,figure)\n")
s=s.replace("def save(fig,name):\n", "def save(fig,name):\n name={'Figure_1':'Supplementary_Figure_S5','Figure_2':'Supplementary_Figure_S6','Figure_4':'Supplementary_Figure_S7'}.get(name,name)\n")
s=s.replace('External blood cohorts do not confirm progression','External blood profiles during scheduled treatment')
exec(compile(s,str(__file__),'exec'))
pd.DataFrame(manifest).to_csv(r/'deliverables/supplementary_context_panel_manifest.csv',index=False)
for n in range(1,5):
 for ext in ['pdf','svg','png']:shutil.copy2(root/'revision6/deliverables/figures'/f'Supplementary_Figure_S{n}.{ext}',r/'deliverables/figures'/f'Supplementary_Figure_S{n}.{ext}')
print('Context moved to S5-S7; previous S1-S4 preserved')
