"""Create PNG/PDF/SVG scientific figures with actual Arial (no silent fallback)."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ros-mpl-cache')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/ros-font-cache')
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from analyze_ros import LABELS

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FIG = HERE / 'figures'
FIG.mkdir(exist_ok=True)
RED, BLUE, GREY = '#C94745', '#3479AB', '#E5E7EB'
GROUPS = ['Canonical scavenger enzyme','Other antioxidant / removal','Glutathione / NADPH support',
    'ROS biosynthesis','Other ROS metabolism','ROS regulation','Oxidative stress response']
SHORT = {'Canonical scavenger enzyme':'Scavenger enzyme','Other antioxidant / removal':'Other antioxidant / removal',
    'Glutathione / NADPH support':'Glutathione / NADPH','ROS biosynthesis':'ROS biosynthesis',
    'Other ROS metabolism':'Other ROS metabolism','ROS regulation':'ROS regulation','Oxidative stress response':'Stress response'}

def read(name): return pd.read_csv(RES / f'{name}.csv')

def save(fig, name):
    for ext in ['png','pdf','svg']:
        fig.savefig(FIG / f'{name}.{ext}', dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def dotplot(acc, tissue, name):
    d = read(f'{acc}_ROS_DEG')
    d['order'] = d.display_group.map({v:i for i,v in enumerate(GROUPS)})
    d = d.sort_values(['order','log2FC']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.8, max(6.8, len(d)*.27+2.2)))
    for group in GROUPS:
        ids = d.index[d.display_group == group]
        if len(ids) and GROUPS.index(group)%2 == 0:
            ax.axhspan(ids.min()-.5,ids.max()+.5,color='#F3F5F7',zorder=0)
    color = np.where(d.log2FC > 0, RED, BLUE)
    sizes = 24 + 20*np.minimum(-np.log10(d.FDR_padj),6)
    ax.hlines(d.index, 0, d.log2FC, color=color, lw=1.1, alpha=.6)
    ax.scatter(d.log2FC, d.index, s=sizes, c=color, edgecolor='white', lw=.6, zorder=3)
    ax.axvline(0, color='#777777', lw=.8)
    ax.set_yticks(d.index, labels=[f'{g}   |   {SHORT[c]}' for g,c in zip(d.gene_symbol,d.display_group)])
    ax.invert_yaxis()
    ax.set_xlabel('log2 fold change (PBM / disease-model control)')
    ax.set_title(f'{acc} | {tissue}\nROS-related DEGs from the supplied FDR < 0.05 table',loc='left',weight='bold',pad=35)
    ax.text(0,1.012,f'{len(d)} genes: {(d.direction=="up").sum()} up, {(d.direction=="down").sum()} down',transform=ax.transAxes,color='#4B5563')
    ax.grid(axis='x',alpha=.17)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y', length=0, labelsize=9)
    handles=[ax.scatter([],[],s=24+20*k,color='#64748B',label=f'FDR = {10**(-k):g}') for k in [2,4,6]]
    ax.legend(handles=handles,title='Point size',loc='lower right',frameon=False,fontsize=8)
    fig.text(.02,.012,'GO/KEGG membership defines ROS relevance; membership does not establish ROS removal or pathway activation.',size=8,color='#4B5563')
    fig.tight_layout(rect=(0,.035,1,1))
    save(fig,name)

def overlap_plot():
    d=read('ROS_DEG_both_studies').sort_values(['concordant','gene_symbol'],ascending=[False,True])
    fig, ax=plt.subplots(figsize=(9,4.3))
    y=np.arange(len(d))
    ax.barh(y-.17,d.log2FC_A,height=.3,color='#6554A4',label='GSE245658 | hippocampus')
    ax.barh(y+.17,d.log2FC_B,height=.3,color='#249C9A',label='GSE325653 | prefrontal cortex')
    ax.axvline(0,color='#777777',lw=.9)
    ax.set_yticks(y,labels=[f'{r.gene_symbol}  ({"same: down" if r.concordant else "opposite"})' for r in d.itertuples()])
    ax.invert_yaxis();ax.set_xlabel('log2 fold change (PBM / disease-model control)')
    ax.set_title('ROS-related genes significant in both supplied DEG tables',loc='left',weight='bold',pad=20)
    ax.legend(loc='lower left',bbox_to_anchor=(0,-.42),frameon=False,ncol=1)
    ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',alpha=.15)
    fig.text(.02,.01,'Penk is an oxidative-stress-response gene; no common-direction antioxidant/removal DEG was identified.',size=9)
    fig.tight_layout(rect=(0,.06,1,1));save(fig,'fig03_shared_ROS_DEG')

def heatmap(d, title, name, source_note, cols, fdr_cols, labels):
    d=d.copy().reset_index(drop=True)
    matrix=d[cols].to_numpy(dtype=float)
    vmax=max(.5,float(np.nanquantile(abs(matrix),.95))) if name=='fig05_GEO_scavenger_panel' else max(.5,float(np.nanmax(abs(matrix))))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad(GREY)
    fig,ax=plt.subplots(figsize=(7.8,max(5.6,.27*len(d)+2.1)))
    im=ax.imshow(np.ma.masked_invalid(matrix),aspect='auto',cmap=cmap,vmin=-vmax,vmax=vmax)
    for i in range(len(d)):
        for j in range(2):
            v=matrix[i,j]
            q=d.iloc[i][fdr_cols[j]]
            text='NA' if not np.isfinite(v) else f'{v:+.2f}'+('*' if np.isfinite(q) and q<.05 else '')
            ax.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if np.isfinite(v) and abs(v)>.7*vmax else '#222222')
    ax.set_yticks(range(len(d)),labels=d.gene_symbol,fontsize=9)
    ax.set_xticks([0,1],labels=labels,fontsize=9)
    ax.tick_params(length=0);ax.spines[:].set_visible(False)
    ax.set_title(title,loc='left',weight='bold',pad=18)
    fig.colorbar(im,ax=ax,pad=.05,shrink=.55,label='log2 fold change')
    fig.text(.02,.015,source_note,fontsize=8,wrap=True)
    fig.tight_layout(rect=(0,.07,1,1));save(fig,name)

def scavenger_figures():
    d=read('ROS_DEG_union')
    hit=d.antioxidant_activity | d.hydrogen_peroxide_catabolism | d.superoxide_removal
    heatmap(d[hit].sort_values('gene_symbol'), 'Antioxidant / ROS-removal DEGs', 'fig04_antioxidant_DEG',
        'Supplied DEG tables only. Grey/NA: absent from the FDR < 0.05 list; effect unknown.\n* Gene FDR < 0.05. GO annotation includes enzymes and other removal-related proteins.',
        ['GSE245658_log2FC','GSE325653_log2FC'],['GSE245658_FDR','GSE325653_FDR'],['Hippocampus\n808 nm','Prefrontal cortex\n1064 nm'])
    p=read('canonical_and_support_reference_panel')
    heatmap(p, 'Canonical scavenger enzymes and selected support genes\nComplete GEO statistics, including non-DEGs', 'fig05_GEO_scavenger_panel',
        '* Gene FDR < 0.05. A: exploratory two-group PyDESeq2 refit of rounded estimated counts.\nB: deposited DESeq2 statistics. NA: no unambiguous finite GEO statistic. Colours clipped at the 95th percentile.',
        ['GSE245658_GEO_log2FC','GSE325653_GEO_log2FC'],['GSE245658_GEO_FDR','GSE325653_GEO_FDR'],['Hippocampus\nGEO refit (exploratory)','Prefrontal cortex\nGEO deposited statistics'])

def pathways():
    d=read('focused_GSEA')
    pivot=d.pivot(index='Term',columns='dataset',values='NES').reindex(list(LABELS))
    fig,ax=plt.subplots(figsize=(11,7))
    y=np.arange(len(pivot))
    for j,acc in enumerate(['GSE245658','GSE325653']):
        vals=pivot[acc];offset=(-.15 if j==0 else .15)
        ax.scatter(vals,y+offset,s=80,marker='o' if j==0 else 's',color='#6554A4' if j==0 else '#249C9A',label=acc,zorder=3)
        for i,n in enumerate(pivot.index):
            r=d[(d.Term==n)&(d.dataset==acc)]
            if len(r) and r.iloc[0].FDR_BH_across_studies<.05:
                ax.text(vals.iloc[i]+.07,i+offset,'*',va='center',fontsize=14)
    ax.axvline(0,color='#777777',lw=.8)
    ax.set_yticks(y,labels=[LABELS[n] for n in pivot.index]);ax.invert_yaxis()
    ax.set_xlabel('GSEA normalized enrichment score (positive: PBM higher)')
    ax.set_title('Exploratory gene-label GSEA | complete GEO statistics',loc='left',weight='bold',pad=18)
    ax.legend(frameon=False,loc='lower right');ax.grid(axis='x',alpha=.15);ax.spines[['top','right','left']].set_visible(False)
    fig.text(.02,.015,'* BH FDR < 0.05 across both studies. Gene-label permutations are exploratory; sample-level exact tests are reported separately.',size=8)
    fig.tight_layout(rect=(0,.04,1,1));save(fig,'fig06_ROS_GSEA')

    legacy=read('legacy_ROS_ORA_global_FDR')
    focus=['GO:0072593','GO:0016209','GO:0042744','GO:0019430','GO:0006749','mmu00480','GO:0006979']
    l=legacy[legacy.set_id.isin(focus)].copy()
    l['column']=l.dataset.str.split('_').str[0]+' | '+l.direction
    mat=l.pivot(index='set_id',columns='column',values='padj_BH_full_catalog').reindex(focus)
    fig,ax=plt.subplots(figsize=(10,5.5))
    im=ax.imshow(-np.log10(np.maximum(mat.to_numpy(),1e-12)),aspect='auto',cmap='Blues',vmin=0,vmax=1.30103)
    for i in range(len(mat)):
        for j in range(len(mat.columns)):
            ax.text(j,i,f'{mat.iloc[i,j]:.3f}',ha='center',va='center',fontsize=9)
    ax.set_yticks(range(len(mat)),labels=[l[l.set_id==s].set_name.iloc[0] for s in mat.index],fontsize=9)
    ax.set_xticks(range(len(mat.columns)),labels=mat.columns,rotation=25,ha='right',fontsize=8)
    ax.set_title('Existing ORA | original full-catalog BH FDR',loc='left',weight='bold',pad=18)
    fig.colorbar(im,ax=ax,label='-log10(FDR); FDR = 0.05 at upper limit')
    fig.text(.02,.01,'Numbers are the original globally adjusted pathway FDR, not gene FDR. None of these ROS terms passes FDR < 0.05.',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1));save(fig,'fig07_existing_ORA_FDR')

def sample_score_plot():
    scores = read('sample_level_pathway_scores')
    tests = read('sample_level_exact_permutation')
    focus=['ros_metabolism','antioxidant_activity','hydrogen_peroxide_catabolism']
    displayed=scores[scores.set_name.isin(focus)].score
    fig,axes=plt.subplots(2,3,figsize=(12.5,7),sharey=True)
    for i,acc in enumerate(['GSE245658','GSE325653']):
        for j,name in enumerate(focus):
            ax=axes[i,j]
            d=scores[(scores.dataset==acc)&(scores.set_name==name)]
            row=tests[(tests.dataset==acc)&(tests.set_name==name)].iloc[0]
            for k,condition in enumerate(['control','PBM']):
                values=d[d.condition==condition].score.to_numpy()
                ax.scatter(k+np.linspace(-.08,.08,len(values)),values,s=42,color=BLUE if k==0 else RED,alpha=.85)
                ax.plot([k-.16,k+.16],[values.mean()]*2,color='#222222',lw=1.7)
            ax.set_xticks([0,1],labels=['Control','PBM'])
            ax.set_title(f'{LABELS[name]}\n{acc}',fontsize=10)
            ax.text(.5,.98,f'exact p = {row.exact_two_sided_p:.3f}\nBH FDR = {row.FDR_BH_across_studies:.3f}',transform=ax.transAxes,ha='center',va='top',fontsize=9)
            ax.spines[['top','right']].set_visible(False)
            ax.set_ylim(displayed.min()-.15,displayed.max()+.45)
            if j==0:ax.set_ylabel('Mean standardized gene expression')
    fig.suptitle('Sample-level pathway scores | exact two-sided label permutations',x=.03,ha='left',weight='bold')
    fig.text(.02,.015,'Each point is one GEO sample. A: log2 normalized estimated counts; B: log2(FPKM + 1).\nBH correction covers 20 set-by-study tests. These displayed sets were selected to illustrate the GSEA signals.',fontsize=8)
    fig.tight_layout(rect=(0,.06,1,.95));save(fig,'fig08_sample_level_pathway_scores')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--font-dir',type=Path,default=Path('/tmp/ros_arial'));args=parser.parse_args()
    fonts=list(args.font_dir.glob('*.TTF'))+list(args.font_dir.glob('*.ttf'))
    if not fonts: raise RuntimeError('Arial files missing: run setup_arial.py or pass --font-dir')
    for path in fonts:font_manager.fontManager.addfont(str(path))
    actual=font_manager.findfont('Arial',fallback_to_default=False)
    mpl.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path','axes.unicode_minus':False})
    dotplot('GSE245658','hippocampus | 808 nm','fig01_hippocampus_ROS_DEG')
    dotplot('GSE325653','prefrontal cortex | 1064 nm','fig02_cortex_ROS_DEG')
    overlap_plot();scavenger_figures();pathways();sample_score_plot()
    (RES/'font_provenance.json').write_text(json.dumps(dict(family='Arial',actual_font=Path(actual).name,
        sha256=hashlib.sha256(Path(actual).read_bytes()).hexdigest(),source_url='https://downloads.sourceforge.net/corefonts/arial32.exe',
        PDF='Embedded TrueType Arial',SVG='Arial glyphs converted to paths to preserve appearance',PNG='Rendered using Arial'),indent=2)+'\n')
    print('Created 8 figures in PNG / PDF / SVG with Arial:',actual)

if __name__=='__main__':main()
