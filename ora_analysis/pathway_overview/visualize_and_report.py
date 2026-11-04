"""Arial figures and a Korean report for per-dataset KEGG/ROS analyses."""
from pathlib import Path
import os
import json
import hashlib
import textwrap
import sys
os.environ.setdefault('MPLCONFIGDIR','/tmp/pathway-mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/pathway-font-cache')
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from analyze_pathways import HERE,ROOT,OLD,OUT,DATASETS,ROS_KEGG

FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
RED,BLUE,GREY='#BF4945','#347CAB','#B9BEC5'
Q='FDR_KEGG_plus_ROS_per_dataset'
SQ='sample_FDR_KEGG_plus_ROS_per_dataset'
OQ='FDR_all_sets_directions_per_dataset'
SHORT={'GO:0016209':'Antioxidant activity','GO:0072593':'ROS metabolism',
 'GO:0042744':'Hydrogen peroxide removal','GO:0019430':'Superoxide removal',
 'GO:1903409':'ROS biosynthesis','GO:0006979':'Oxidative stress response',
 'GO:2000377':'ROS metabolism regulation','GO:0006740':'NADPH regeneration',
 'GO:0006749':'Glutathione metabolism (GO)',
 'mmu05208':'Chemical carcinogenesis - ROS',
 'mmu00480':'Glutathione metabolism (KEGG)'}

def read(acc,name):return pd.read_csv(OUT/acc/f'{name}.csv')
def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',dpi=220,facecolor='white')
    plt.close(fig)
def source(acc):return 'Hippocampus, 808 nm | exploratory GEO count refit' if acc==DATASETS[0] else 'Prefrontal cortex, 1064 nm | deposited GEO DESeq2 statistics'
def short(s,n=42):return '\n'.join(textwrap.wrap(str(s),n))

def top_gsea(acc,g):
    fig,axes=plt.subplots(1,2,figsize=(15,7))
    for j,(direction,title) in enumerate([(1,'PBM-higher enrichment'),(-1,'PBM-lower enrichment')]):
        ax=axes[j];d=g[(g.source=='KEGG')&(np.sign(g.NES)==direction)].sort_values(['nominal_p','set_id']).head(12).iloc[::-1]
        color=[(RED if direction>0 else BLUE) if q<.05 else GREY for q in d[Q]]
        ax.barh(range(len(d)),abs(d.NES),color=color,height=.65)
        ax.set_yticks(range(len(d)),labels=[short(x,36) for x in d.set_name],fontsize=8)
        xmax=max(2.5,float(abs(d.NES).max())+.75);ax.set_xlim(0,xmax)
        for i,r in enumerate(d.itertuples()):ax.text(abs(r.NES)+.04,i,f'q={getattr(r,Q):.3g}',va='center',fontsize=8)
        ax.set_title(title,loc='left',weight='bold',fontsize=11)
        ax.set_xlabel('|NES| (sign indicated by panel)')
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.13)
    fig.suptitle(f'{acc} | whole-KEGG preranked GSEA\n{source(acc)}',ha='left',x=.015,weight='bold',fontsize=13)
    fig.text(.015,.012,'Top 12 by nominal P within each direction; grey = q >= 0.05. q: BH across all eligible KEGG + 9 pre-specified ROS GO sets in this dataset.\nCompetitive gene-label permutations; pathway names referring to diseases do not establish those diseases or their improvement.',fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.92));save(fig,f'{acc}_01_KEGG_GSEA_overview')

def ros_panel(acc,g):
    order=['GO:0016209','GO:0042744','GO:0019430','GO:0072593','GO:1903409','GO:2000377','GO:0006979','GO:0006740','GO:0006749']+list(ROS_KEGG)
    d=g.set_index('set_id').reindex(order).dropna(subset=['NES']).reset_index()
    fig,(ax,tableax)=plt.subplots(1,2,figsize=(13,8),gridspec_kw={'width_ratios':[1.6,1]})
    colors=np.where(d[Q]<.05,np.where(d.NES>0,RED,BLUE),GREY)
    ax.hlines(range(len(d)),0,d.NES,color=colors,alpha=.55,lw=1.2)
    ax.scatter(d.NES,range(len(d)),c=colors,s=68,zorder=3,edgecolor='white',lw=.6)
    ax.axvline(0,color='#777777',lw=.8);ax.set_xlim(-2.6,1.7);ax.set_ylim(len(d)-.5,-1.6)
    ax.set_yticks(range(len(d)),labels=[short(SHORT.get(r.set_id,r.set_name),37) for r in d.itertuples()],fontsize=9)
    ax.set_xlabel('NES: negative = PBM lower, positive = PBM higher')
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.15)
    tableax.set_xlim(0,1);tableax.set_ylim(len(d)-.5,-1.6);tableax.axis('off')
    headings=[('Gene-label q',.05),('Genes',.41),('P-rank',.61),('Sample q',.86)]
    for text,x in headings:tableax.text(x,-1,text,ha='center',weight='bold',fontsize=9)
    for i,r in enumerate(d.itertuples()):
        for value,x in [(f'{getattr(r,Q):.3g}',.05),(str(r.matched_genes),.41),(str(r.rank_by_p),.61),(f'{getattr(r,SQ):.3g}',.86)]:
            tableax.text(x,i,value,ha='center',va='center',fontsize=9,weight='bold' if getattr(r,Q)<.05 and x==.05 else 'normal')
    fig.suptitle(f'{acc} | ROS / antioxidant pathway context\n{source(acc)}',x=.015,ha='left',weight='bold',fontsize=13)
    fig.text(.015,.012,'Coloured points: gene-label q < 0.05; grey: not significant. Both q columns are corrected over the full eligible set collection in this dataset.\nP-rank is within all sets. KEGG context terms include broad redox signaling, not exclusively ROS scavenging. Sample-score and GSEA tests have different nulls.',fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.91));save(fig,f'{acc}_02_ROS_pathway_significance')

def curves(acc,g):
    ranks=pd.read_csv(OLD/f'{acc}_GSEA_ranked_genes.csv')
    cat=pd.read_csv(OUT/'gene_set_catalog.csv').set_index('set_id')
    symbols=ranks.gene_symbol.to_numpy();weights=abs(ranks.signed_z.to_numpy());N=len(ranks)
    fig=plt.figure(figsize=(15,6.6));grid=fig.add_gridspec(3,3,height_ratios=[3,.35,1],hspace=.14,wspace=.28)
    for j,sid in enumerate(['GO:0016209','mmu00480','GO:0042744']):
        r=g[g.set_id==sid].iloc[0];hit=np.isin(symbols,cat.loc[sid,'genes'].split(';'))
        increments=np.full(N,-1/(N-hit.sum()));increments[hit]=weights[hit]/weights[hit].sum()
        walk=np.r_[0,np.cumsum(increments)];xs=np.arange(N+1)/N*100
        ax=fig.add_subplot(grid[0,j]);rug=fig.add_subplot(grid[1,j],sharex=ax);metric=fig.add_subplot(grid[2,j],sharex=ax)
        ax.plot(xs,walk,color='#428B56',lw=1.5);ax.axhline(0,color='#888888',lw=.7)
        ax.set_title(f'{SHORT[sid]}\nNES={r.NES:+.2f}, q={r[Q]:.3g}',fontsize=11,weight='bold',pad=12)
        ax.set_ylabel('Running enrichment score' if j==0 else '')
        ax.tick_params(labelbottom=False);ax.spines[['top','right']].set_visible(False)
        rug.vlines((np.flatnonzero(hit)+1)/N*100,0,1,color='#222222',lw=.7);rug.set_axis_off()
        metric.plot(np.arange(N)/N*100,ranks.signed_z,color='#727A83',lw=.7)
        metric.axhline(0,color='#888888',lw=.5);metric.set_xlim(0,100);metric.set_xticks([0,25,50,75,100])
        metric.set_xlabel('Rank percentile: PBM higher -> PBM lower',fontsize=9)
        if j==0:metric.set_ylabel('Signed z',fontsize=9)
        metric.spines[['top','right']].set_visible(False)
    fig.suptitle(f'{acc} | enrichment curves for antioxidant / scavenger functions\n{source(acc)}',x=.015,ha='left',weight='bold',fontsize=13)
    fig.text(.015,.015,'Black ticks: gene-set members. Entire ranked transcriptome is used, including non-DEGs.\nq = BH over eligible KEGG + 9 pre-specified ROS GO sets. Enrichment reflects transcript ranks, not measured enzyme activity or ROS flux.',fontsize=8)
    fig.subplots_adjust(top=.76,bottom=.18,left=.055,right=.985);save(fig,f'{acc}_03_ROS_GSEA_curves')

def ora_plot(acc,o):
    fig,axes=plt.subplots(1,2,figsize=(15,7))
    for j,direction in enumerate(['up','down']):
        ax=axes[j];d=o[(o.source=='KEGG')&(o.direction==direction)&(o.k>0)].sort_values('pvalue').head(12).iloc[::-1]
        vals=-np.log10(np.maximum(d[OQ],1e-16));colors=[(RED if direction=='up' else BLUE) if q<.05 else GREY for q in d[OQ]]
        ax.scatter(vals,range(len(d)),s=25+5*np.minimum(d.k,40),c=colors,edgecolor='white',lw=.6,zorder=3)
        ax.set_yticks(range(len(d)),labels=[short(x,37) for x in d.set_name],fontsize=8)
        ax.axvline(-np.log10(.05),color='#555555',ls='--',lw=.8)
        ax.set_xlabel('-log10(BH FDR)');ax.set_title(f'{direction.upper()} DEGs | top KEGG terms',loc='left',weight='bold',fontsize=11)
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',alpha=.13)
        ax.set_xlim(-.1,max(2,float(vals.max())+.3))
    fig.suptitle(f'{acc} | KEGG over-representation analysis\n{source(acc)}',x=.015,ha='left',weight='bold',fontsize=13)
    fig.text(.015,.012,'FDR family: all eligible KEGG + 9 ROS GO sets x all/up/down in this dataset. Dashed line: FDR = 0.05. Point size: DEG hits (capped at 40).\nORA uses gene FDR < 0.05; a pathway may be significant in GSEA but not ORA, because the selection rules and nulls differ.',fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.92));save(fig,f'{acc}_04_KEGG_ORA_overview')

def score_plot(acc,g):
    scores=read(acc,'sample_pathway_scores');fig,axes=plt.subplots(1,3,figsize=(12,4.8))
    focus=['GO:0016209','mmu00480','GO:0042744']
    for j,sid in enumerate(focus):
        ax=axes[j];d=scores[scores.set_id==sid];r=g[g.set_id==sid].iloc[0]
        values_all=d.score.to_numpy();lo,hi=min(values_all),max(values_all);span=max(.4,hi-lo)
        for k,group in enumerate(['control','PBM']):
            values=d[d.condition==group].score.to_numpy()
            ax.scatter(k+np.linspace(-.08,.08,len(values)),values,c=BLUE if k==0 else RED,s=46)
            ax.plot([k-.16,k+.16],[values.mean()]*2,c='#222222',lw=1.6)
        ax.set_xticks([0,1],labels=['Control','PBM']);ax.set_ylim(lo-.12*span,hi+.5*span)
        ax.set_title(SHORT[sid],fontsize=10,weight='bold')
        ax.text(.5,.98,f'exact P={r.sample_exact_p:.3g}\nBH FDR={r[SQ]:.3g}',transform=ax.transAxes,ha='center',va='top',fontsize=9)
        if j==0:ax.set_ylabel('Mean standardized gene expression')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle(f'{acc} | complementary sample-level pathway scores',x=.015,ha='left',weight='bold',fontsize=13)
    fig.text(.015,.015,'Each point is one GEO sample. Exact two-sided sample-label permutations preserve inter-gene correlation.\nFDR corrected over the full eligible set collection in this dataset. This score test has a different null from competitive GSEA.',fontsize=8)
    fig.tight_layout(rect=(0,.09,1,.9));save(fig,f'{acc}_05_ROS_sample_scores')

def mdtable(df):
    def fmt(v):return f'{v:.4g}' if isinstance(v,float) else str(v).replace('|',' / ')
    return '\n'.join(['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)])

def report():
    comparison=pd.read_csv(OUT/'ROS_pathway_comparison.csv')
    highlights=['GO:0016209','GO:0042744','GO:0019430','GO:0072593','mmu00480','mmu05208']
    c=comparison[comparison.set_id.isin(highlights)][['dataset','set_id','set_name','NES',Q,'rank_by_p','matched_genes','leading_edge_count',SQ]].copy()
    c['set_name']=[SHORT.get(i,n) for i,n in zip(c.set_id,c.set_name)]
    c.columns=['Dataset','ID','기능군','NES','GSEA BH FDR','P-rank','Genes','Leading edge','Sample BH FDR']
    provenance=json.loads((OUT/'provenance.json').read_text())
    parts=['''# 데이터셋별 전체 KEGG 및 ROS·항산화 경로 분석

기존 저장소의 `ROS_scavenger_20261103` 브랜치에는 ROS 기능군 10개의 [GSEA 요약 그림](../ros_scavenger/figures/fig06_ROS_GSEA.png), [기존 ORA 그림](../ros_scavenger/figures/fig07_existing_ORA_FDR.png), [표본 단위 score 그림](../ros_scavenger/figures/fig08_sample_level_pathway_scores.png)이 있었다. **전체 KEGG 경로를 비교하는 데이터셋별 그림은 없어 이번에 추가했다.** 그림의 실제 폰트는 Arial이며 PNG·PDF·SVG로 저장했다.

## 핵심 해석

- **GSE245658(해마):** 전체 경로를 포함한 검정에서 antioxidant activity는 감소 방향이며 GSEA FDR<0.05다. ROS metabolism과 그 regulation도 감소 방향의 탐색적 신호다. 다만 직접적인 hydrogen peroxide/superoxide removal 기능군과 KEGG glutathione metabolism은 유의하지 않다.
- **GSE325653(전전두피질):** antioxidant activity와 hydrogen peroxide removal의 NES는 음수지만, 이번 전체 경로 검정에서는 FDR≥0.05다. Oxidative stress response는 감소 방향의 GSEA 신호가 있고, KEGG glutathione metabolism은 유의하지 않다.
- **두 연구 모두 KEGG chemical carcinogenesis–ROS는 감소 방향 GSEA 신호가 있다.** 이 집합은 신호 전달·대사·호흡사슬 등을 함께 포함하는 넓은 경로이며 직접적인 scavenger 경로가 아니다. 암 발생·치료 효과나 실제 ROS 감소로 해석하지 않는다.
- 표본 단위 평균 경로 score 검정은 전체 경로 보정 후 유의한 항목이 없다. 따라서 gene-label GSEA 결과를 생체 항산화 능력 저하 또는 ROS 생성량의 변화로 확정할 수 없다.

여기서 ‘영향력’은 **전사체 순위에서의 상대적 enrichment(NES), 전체 경로 중 위치, FDR, 기여 유전자 범위**로 제시한다. 인과적 기여율·ROS 제거 능력·생체 내 flux의 크기를 계산한 것은 아니다.

## 비교 조건과 통계 출처

| Dataset | 비교 | 시료 | 전체 통계 |
|---|---|---|---|
| GSE245658 | 해마 808 nm, 5xFAD+light / 5xFAD | 4 / 4 | 기존 ROS 분석에서 저장한 GEO 추정 count 반올림 기반의 새 두 군 PyDESeq2 fit |
| GSE325653 | 전전두피질 1064 nm, 저산소+PBM / 저산소 | 5 / 5 | 원 GEO 제출 DESeq2 통계; FPKM을 count로 재계산하지 않음 |

원 GEO 파일·시료 확인, A의 소수 count·길이 보정 부재·새 fit의 한계, B의 기존 DEG 누락 21개는 [기존 상세 보고서](../ros_scavenger/ROS_scavenger_report_ko.md)에 설명돼 있다. 이번 분석은 이미 확보한 전체 통계와 순위 파일을 재사용했다. 조직·파장·질환 배경이 함께 다르므로 두 연구의 차이를 순수한 피질/해마 또는 파장 효과로 분해하지 않는다.

## 전체 경로에서 ROS 기능군의 위치

아래 FDR는 데이터셋별 **분석 가능한 전체 KEGG + 사전 지정 ROS GO 9개**에 대한 BH 값이다. P-rank는 해당 전체 집합에서 nominal P 순위이며 작은 값이 더 강한 통계적 enrichment다. 두 연구를 함께 보정한 민감도 값과 KEGG만 보정한 값도 Excel/CSV에 저장했다.
''',mdtable(c)]
    parts.append('''
기존 ROS 한정 분석의 FDR는 10기능군×2연구=20검정에 대한 값이다. **이번 FDR와 같은 값이 아니다.** GO의 antioxidant activity는 A에서 이번에도 유의하지만, B에서는 기존 ROS 한정 분석의 유의성이 전체 KEGG를 포함한 검정으로 유지되지 않는다. BH는 검정 수뿐 아니라 전체 P값 분포에도 영향을 받으며, 데이터셋별/두 연구 합산 보정과 Monte Carlo 표본 수·seed도 달라 단순히 ‘검정이 많아져 FDR가 항상 커진다’고 해석하지 않는다.

GO antioxidant activity는 peroxidase·환원효소 등도 포함한다. Ptgs1/2 같은 prostaglandin 대사 유전자까지 들어갈 수 있으므로 전체를 순수한 직접 scavenger 집합이라고 부르지 않는다. 제거 기능 자체는 별도의 H2O2 catabolism 및 superoxide removal로 구분했다.
''')
    for acc in DATASETS:
        g=read(acc,'GSEA_all_pathways');o=read(acc,'ORA_all_pathways');audit=json.loads((OUT/acc/'audit.json').read_text())
        a=g[g.set_id=='GO:0016209'].iloc[0];gh=g[g.set_id=='mmu00480'].iloc[0]
        ox=g[g.set_id=='mmu00190'].iloc[0];rib=g[g.set_id=='mmu03010'].iloc[0]
        stronger=int((g.nominal_p<a.nominal_p).sum())
        parts.append(f'''\n## {acc}: 데이터셋별 결과

분석된 KEGG는 **{audit['eligible_KEGG']}개**, ROS GO는 9개로 GSEA family는 {audit['GSEA_family_size']}개다. 순위 유전자 수는 {audit['ranked_genes']:,}개다. 전체 확장 family에서 gene-label GSEA FDR<0.05인 검정은 {audit['GSEA_significant']}개이며, KEGG만 보정했을 때 유의한 KEGG는 {audit['GSEA_KEGG_significant']}개다. 서로 겹치는 경로를 독립적인 현상 수로 세지 않는다.

**Antioxidant activity:** NES {a.NES:+.3f}, FDR {a[Q]:.4g}, P 순위 {int(a.rank_by_p)}/{len(g)}. |NES|는 전체 집합의 {a.abs_NES_percentile:.1f} percentile(절대 NES가 더 작은 집합의 비율)에 해당한다. {int(a.matched_genes)}개 matched gene 중 leading edge는 {int(a.leading_edge_count)}개다. 이 크기·순위는 relative enrichment이며 설명 분산이나 항산화 능력의 백분율이 아니다.

비교를 위해 oxidative phosphorylation은 NES {ox.NES:+.3f}/FDR {ox[Q]:.4g}, ribosome은 NES {rib.NES:+.3f}/FDR {rib[Q]:.4g}다. 항산화 집합보다 nominal P가 작은 경로는 {stronger}개 있다. 따라서 항산화 축만이 유일하거나 가장 강한 전사체 반응이라고 보지 않는다.

**KEGG glutathione metabolism:** NES {gh.NES:+.3f}, nominal P {gh.nominal_p:.4g}, 전체 family FDR {gh[Q]:.4g}; KEGG-only FDR {gh.FDR_KEGG_only_per_dataset:.4g}. 이 경로는 FDR<0.05가 아니다.

![전체 KEGG GSEA](figures/{acc}_01_KEGG_GSEA_overview.png)

![ROS 기능군의 방향과 유의성](figures/{acc}_02_ROS_pathway_significance.png)

![ROS 기능군 GSEA 곡선](figures/{acc}_03_ROS_GSEA_curves.png)

### 기여 유전자

Antioxidant activity의 leading edge: **{a.leading_edge_genes.replace(';',', ')}**.

전체 leading edge 목록과 각 gene log2FC/FDR는 [`{acc}/ROS_leading_edge_gene_statistics.csv`](results/{acc}/ROS_leading_edge_gene_statistics.csv)에 있다. Leading edge 포함이 gene-level FDR<0.05나 인과적인 driver임을 뜻하지 않는다. GSEA는 많은 작은 변화를 합칠 수 있어, 유의한 DEG가 적어도 pathway enrichment가 나타날 수 있다.

### KEGG ORA와 표본 단위 검정

ORA는 gene FDR<0.05의 전체/up/down DEG를 사용했고, {audit['ORA_family_size']}개 경로×방향 검정 전체에 BH를 적용했다. {audit['ORA_significant']}개 검정이 유의하다. 기능군 10개만 분석한 기존 focused ORA와 보정 family가 달라 기존 수치를 덮어쓰지 않았다.

![전체 KEGG ORA](figures/{acc}_04_KEGG_ORA_overview.png)

![표본 단위 ROS 경로 score](figures/{acc}_05_ROS_sample_scores.png)

표본 score 검정에서 FDR<0.05인 항목은 {audit['sample_significant']}개다. Antioxidant activity raw exact P={a.sample_exact_p:.4g}, BH FDR={a[SQ]:.4g}다. 이 검정은 유전자 간 상관을 보존하지만 평균 표준화 발현 차이를 검정하는 별도 self-contained 가설이며 competitive GSEA와 귀무가설이 다르다. 둘을 같은 검정의 통과/실패로 부르지 않는다.
''')
        ros_ora=o[(o.ROS_related)&(o[OQ]<.05)]
        if len(ros_ora):
            parts.append('ROS context에 지정한 항목 중 ORA 유의 결과:\n\n'+mdtable(ros_ora[['set_name','direction','k','K','fold_enrichment',OQ]])+'\n\nFoxO 같은 넓은 신호 경로의 ORA를 직접적인 scavenger 활성의 증거로 해석하지 않는다.\n')
        else:parts.append('사전 지정 ROS/항산화 context 항목의 ORA는 모두 다중검정 후 유의하지 않다.\n')
    parts.append('''
## 분석 방법과 해석 범위

1. KEGG의 최신 마우스 pathway 목록과 gene–pathway link를 직접 내려받았다. 368개 경로를 catalog에 보관하고, 발현/통계 universe와 교차한 크기가 15–500인 KEGG만 검정했다. 전체 대사 map이나 작은 경로가 제외된 이유는 `set_coverage.csv`에 남겼다. ROS GO 9개는 기존 QuickGO snapshot을 그대로 사용하고 크기 5–500으로 검정했다. 이는 전체 GO catalog 분석이 아니라 사전에 정의한 ROS 기능군을 KEGG에 추가한 분석이다.
2. KEGG Entrez를 저장된 공식 mouse gene-name 목록의 primary symbol로 매핑했다. 유효하고 고유한 심볼만 있는 기존 순위 universe를 유지한다. 모든 alias를 임의로 합치지 않았으며 매핑되지 않은 ID와 matched size는 audit에 남겼다. 명칭이 바뀐 유전자 등의 매핑 한계가 있다.
3. 전체 순위는 `sign(log2FC) × normal_inverse_survival(P/2)`이며 DEG로 자르지 않았다. NES<0은 PBM에서 감소 쪽 순위에 상대적으로 몰린다는 뜻이다. 제거와 생성 유전자가 섞인 경로의 음의 NES는 순 ROS 농도 감소나 항산화 기능 향상을 뜻하지 않는다.
4. Weighted GSEA 지수 1, streaming gene-label null을 사용했다. 전체 catalog에는 크기별 10,000회, 사전 지정 ROS 관련 set의 matched size에는 50,000회 순열을 계산했다. 동일한 크기의 집합은 같은 gene-label null 분포를 가지므로 null을 재사용한다. seed는 `20261104 + crc32(dataset:K)`다. 같은 부호 null의 평균 |ES|로 NES를 계산하고, tail P에는 +1을 적용했다. nominal P의 Monte Carlo Wilson 구간도 보관했다. 이는 생물학적 신뢰구간이 아니다.
5. 1차 GSEA 보정은 데이터셋별 전체 eligible KEGG+ROS GO의 nominal P에 대한 BH다. KEGG-only와 두 데이터셋 전체 합산 보정도 민감도로 저장했다. GSEA의 pooled-NES FDR q와는 다른, 명시적인 BH 보정이다. p순위 동률이 있을 수 있으며 permutation 해상도 이하로 순위를 세밀하게 단정하지 않는다.
6. KEGG ORA는 finite gene FDR의 전체 universe에서 FDR<0.05 전체/up/down을 각각 hypergeometric 우측검정했고 데이터셋별 전체 set×3방향에 BH를 적용했다. DEG 수나 pathway 이름만으로 GSEA 방향을 판정하지 않았다.
7. 보완적인 표본 score는 A log2(normalized estimated count+1), B log2(FPKM+1)에서 각 유전자를 전체 시료 기준 z-standardization한 후 유전자 평균으로 정의했다. A 70개/B 252개의 모든 균형 시료 라벨 배정을 열거한 양측 exact P에 데이터셋별 전체 경로 BH를 적용했다. 원 시료의 batch/짝지음 정보가 확인되지 않아 독립 시료를 가정한다. 작은 n에서 비유의 결과가 효과 부재나 충분한 검출력을 뜻하지는 않는다.

**판단:** 해마에는 항산화 주석 집합의 상대적 감소라는 탐색적 전사체 신호가 있다. 전전두피질은 ROS/스트레스의 넓은 경로 신호는 있지만, 직접적인 scavenger·antioxidant 경로가 전체 경로 보정에서 명확히 유의하다고 하기 어렵다. 두 연구 모두 glutathione metabolism 증가·감소가 통계적으로 확증됐다는 결론은 지지하지 않는다. 단백질 활성·GSH/GSSG·ROS 측정과 연결하기 전에는 기능적 영향의 크기로 번역하지 않는다.

## 파일

- [통합 Excel](results/Dataset_pathway_GSEA_KEGG.xlsx): 모든 경로, ROS 비교, ORA, leading edge.
- [`results/ROS_pathway_comparison.csv`](results/ROS_pathway_comparison.csv): NES·각 FDR family·전체 순위·leading edge 수·표본 score.
- `results/GSE245658/`, `results/GSE325653/`: 각 데이터셋 전체 결과와 set coverage.
- `results/gene_set_catalog.csv`, `kegg_gene_mapping.csv`: 경로 정의와 매핑.
- `results/provenance.json`, `../input/pathway_overview/source_manifest.json`: 입력 해시와 다운로드 근거.
- `figures/`: 각 데이터셋 5개씩, 총 10개 그림을 PNG·PDF·SVG로 저장.
- [재실행 방법](README.md).

이번 전체 경로 분석은 기존 ROS 한정 분석과 구분해 보관한다. 코드·입력·결과 및 Arial 그림은 GitHub의 `ROS_scavenger_20261103` 브랜치에 함께 저장한다.
''')
    (HERE/'Dataset_pathway_report_ko.md').write_text('\n\n'.join(parts))

def main():
    for p in Path('/tmp/ros_arial').glob('*.TTF'):font_manager.fontManager.addfont(str(p))
    font=Path(font_manager.findfont('Arial',fallback_to_default=False))
    mpl.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path','axes.unicode_minus':False})
    for acc in DATASETS:
        g=read(acc,'GSEA_all_pathways');o=read(acc,'ORA_all_pathways')
        top_gsea(acc,g);ros_panel(acc,g);curves(acc,g);ora_plot(acc,o);score_plot(acc,g)
    report()
    (OUT/'font_provenance.json').write_text(json.dumps(dict(family='Arial',font=font.name,sha256=hashlib.sha256(font.read_bytes()).hexdigest(),PDF='Embedded TrueType',SVG='Arial glyphs as paths',figure_count=10),indent=2)+'\n')
    print('Saved 10 Arial figures (PNG/PDF/SVG) and Korean report',flush=True)

if __name__=='__main__':main()
