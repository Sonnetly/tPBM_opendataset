"""Dataset-specific whole-KEGG GSEA/ORA with pre-specified ROS GO comparators.

Uses saved complete GEO statistics/ranks from the ROS analysis. Gene-label nulls
are streamed and reused by matched set size: their distribution depends on the
ranking weights, universe size and K, not the identity of the observed pathway.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import os
import sys
import zlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/pathway-mpl')
import numpy as np
import pandas as pd
from scipy import stats

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'ora_analysis/ros_scavenger/results'
INP=ROOT/'ora_analysis/input/pathway_overview'
OUT=HERE/'results'
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'ora_analysis/ros_scavenger'))
from analyze_ros import enrichment_score, LABELS

DATASETS=['GSE245658','GSE325653']
ROS_KEGG={
 'mmu00480':'Glutathione metabolism: detoxification and GSH support',
 'mmu04146':'Peroxisome: ROS generation and removal',
 'mmu04216':'Ferroptosis: lipid peroxidation / GSH-GPX4 context',
 'mmu05208':'Chemical carcinogenesis - ROS: broad ROS signaling context',
 'mmu00030':'Pentose phosphate pathway: NADPH supply',
 'mmu00450':'Selenocompound metabolism: selenoprotein context',
 'mmu04068':'FoxO signaling: broad stress-response context',
 'mmu04066':'HIF-1 signaling: broad hypoxia / redox context',
}

def save(df,path):df.to_csv(path,index=False)
def bh(p):return stats.false_discovery_control(np.asarray(p,float),method='bh')

def prepare():
    names={}
    for line in (INP/'kegg_mouse_pathways.tsv').read_text().splitlines():
        key,value=line.split('\t')
        names[key.removeprefix('path:')]=value.removesuffix(' - Mus musculus (house mouse)')
    gene_names={}
    for line in (ROOT/'ora_analysis/input/ros_geo/kegg_mouse_gene_names.tsv').read_text().splitlines():
        cols=line.split('\t')
        if ';' in cols[-1]:gene_names[cols[0]]=cols[-1].split(';')[0].split(', ')[0]
    links=[]
    for line in (INP/'kegg_mouse_gene_pathways.tsv').read_text().splitlines():
        gene,path=line.split('\t');path=path.removeprefix('path:')
        links.append(dict(set_id=path,entrez=gene.removeprefix('mmu:'),gene_symbol=gene_names.get(gene)))
    link=pd.DataFrame(links)
    save(link,OUT/'kegg_gene_mapping.csv')
    rows=[]
    for sid,s in link.dropna(subset=['gene_symbol']).groupby('set_id'):
        rows.append(dict(set_id=sid,set_name=names[sid],source='KEGG',
            min_size=15,max_size=500,ROS_related=sid in ROS_KEGG,
            ROS_context=ROS_KEGG.get(sid,''),genes=';'.join(sorted(set(s.gene_symbol)))))
    catalog=pd.read_csv(OLD/'ros_gene_catalog.csv')
    definitions=pd.read_csv(OLD/'gene_set_definitions.csv')
    for r in definitions.itertuples():
        if r.set_name=='kegg_glutathione':continue
        genes=catalog.loc[catalog[r.set_name], 'gene_symbol']
        rows.append(dict(set_id=r.root_id,set_name=LABELS[r.set_name],source='GO_ROS',
            min_size=5,max_size=500,ROS_related=True,ROS_context=r.set_name,genes=';'.join(sorted(genes))))
    sets=pd.DataFrame(rows).sort_values('set_id')
    save(sets,OUT/'gene_set_catalog.csv')
    assert set(ROS_KEGG)<=set(sets.set_id)
    assert sets.set_id.is_unique
    print('Prepared',len(sets.query('source=="KEGG"')),'KEGG pathways +',len(sets.query('source=="GO_ROS"')),'ROS GO sets',flush=True)

def sample_scores(acc,ranks,eligible):
    if acc==DATASETS[0]:
        expr=pd.read_csv(OLD/f'{acc}_log2_normalized_counts.csv.gz',index_col=0)
        mapping=pd.read_csv(OLD/f'{acc}_GEO_full_statistics.csv').set_index('gene_id').gene_symbol
        expr.index=mapping.loc[expr.index].values
    else:expr=pd.read_csv(OLD/f'{acc}_log2_FPKM.csv.gz',index_col=0)
    # Duplicated symbols are excluded, consistent with complete ranking universe.
    expr=expr[~expr.index.duplicated(keep=False)]
    expr=expr.loc[expr.index.intersection(ranks.gene_symbol)]
    std=expr.std(axis=1,ddof=1)
    expr=expr.loc[std>1e-8]
    z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1,ddof=1),axis=0)
    n=expr.shape[1]//2
    masks=np.zeros((len(list(itertools.combinations(range(2*n),n))),2*n))
    for i,indices in enumerate(itertools.combinations(range(2*n),n)):masks[i,list(indices)]=1
    contrasts=(2*masks-1)/n
    rows,scores=[],[]
    for r in eligible.itertuples():
        members=z.index.intersection(r.genes.split(';'))
        score=z.loc[members].mean(axis=0).to_numpy()
        diff=float(score[n:].mean()-score[:n].mean())
        null=contrasts@score
        p=float(np.mean(abs(null)>=abs(diff)-1e-12))
        pooled=np.sqrt((score[:n].var(ddof=1)+score[n:].var(ddof=1))/2)
        rows.append(dict(set_id=r.set_id,score_genes=len(members),score_difference=diff,
            score_standardized_difference=diff/pooled if pooled else np.nan,sample_exact_p=p,
            exact_label_assignments=len(masks)))
        for col,value in zip(expr.columns,score):
            scores.append(dict(set_id=r.set_id,sample=col,condition='control' if col in expr.columns[:n] else 'PBM',score=value))
    return pd.DataFrame(rows),pd.DataFrame(scores)

def analyze(acc):
    dest=OUT/acc;dest.mkdir(exist_ok=True)
    catalog=pd.read_csv(OUT/'gene_set_catalog.csv').fillna('')
    ranks=pd.read_csv(OLD/f'{acc}_GSEA_ranked_genes.csv')
    assert ranks.gene_symbol.is_unique and ranks.signed_z.is_monotonic_decreasing
    symbols=ranks.gene_symbol.to_numpy();weights=abs(ranks.signed_z.to_numpy());N=len(ranks)
    catalog['matched_genes']=catalog.genes.map(lambda s:len(set(s.split(';'))&set(symbols)))
    catalog['eligible']=catalog.matched_genes.between(catalog.min_size,catalog.max_size)
    save(catalog,dest/'set_coverage.csv')
    eligible=catalog[catalog.eligible].copy()
    # 10k for the catalog; 50k nulls at ROS-related set sizes for greater precision.
    permutations={int(k):10000 for k in eligible.matched_genes}
    for k in eligible.loc[eligible.ROS_related,'matched_genes']:permutations[int(k)]=50000
    cache=dest/'null_cache';cache.mkdir(exist_ok=True)
    signature=hashlib.sha256(ranks.to_csv(index=False).encode()).hexdigest()[:12]
    nulls={}
    for i,(K,B) in enumerate(sorted(permutations.items())):
        file=cache/f'{signature}_K{K}_B{B}.npz'
        if file.exists():null=np.load(file)['ES']
        else:
            rng=np.random.default_rng(20261104+zlib.crc32(f'{acc}:{K}'.encode()))
            null=np.empty(B)
            for j in range(B):null[j]=enrichment_score(np.sort(rng.choice(N,K,replace=False)),weights)[0]
            np.savez_compressed(file,ES=null)
        nulls[K]=null
        if i%20==0 or i+1==len(permutations):print(acc,'null sizes',i+1,'/',len(permutations),flush=True)
    results=[]
    for r in eligible.itertuples():
        pos=np.flatnonzero(np.isin(symbols,r.genes.split(';')))
        es,boundary=enrichment_score(pos,weights);null=nulls[len(pos)]
        same=null[null>=0] if es>=0 else null[null<0]
        exceed=int(np.sum(abs(same)>=abs(es)))
        p=(exceed+1)/(len(same)+1)
        # Wilson interval describes Monte Carlo uncertainty, not biological confidence.
        lo,hi=stats.binomtest(exceed,len(same)).proportion_ci(confidence_level=.95,method='wilson')
        lead=pos[pos<=boundary] if es>=0 else pos[pos>=boundary]
        results.append(dict(set_id=r.set_id,set_name=r.set_name,source=r.source,ROS_related=r.ROS_related,
            ROS_context=r.ROS_context,matched_genes=len(pos),universe=N,ES=es,NES=es/np.mean(abs(same)),
            nominal_p=p,MC_p_lower=lo,MC_p_upper=hi,permutations=len(null),same_sign_null=len(same),
            null_exceedances=exceed,leading_edge_count=len(lead),leading_edge_genes=';'.join(symbols[lead])))
    gsea=pd.DataFrame(results)
    gsea['FDR_KEGG_plus_ROS_per_dataset']=bh(gsea.nominal_p)
    gsea['FDR_KEGG_only_per_dataset']=np.nan
    mask=gsea.source=='KEGG'
    gsea.loc[mask,'FDR_KEGG_only_per_dataset']=bh(gsea.loc[mask,'nominal_p'])
    gsea['rank_by_p']=gsea.nominal_p.rank(method='min').astype(int)
    gsea['abs_NES_percentile']=gsea.NES.abs().rank(pct=True)*100
    sample,scores=sample_scores(acc,ranks,eligible)
    sample['sample_FDR_KEGG_plus_ROS_per_dataset']=bh(sample.sample_exact_p)
    gsea=gsea.merge(sample,on='set_id',validate='one_to_one')
    save(gsea.sort_values('nominal_p'),dest/'GSEA_all_pathways.csv')
    save(gsea[gsea.ROS_related].sort_values('nominal_p'),dest/'ROS_pathway_summary.csv')
    save(scores,dest/'sample_pathway_scores.csv')
    # GEO ORA, independent of the GSEA null. Family = all eligible sets x 3 directions.
    full=pd.read_csv(OLD/f'{acc}_GEO_full_statistics.csv')
    full=full[full.gene_symbol.isin(symbols)&np.isfinite(full.FDR_padj)]
    assert full.gene_symbol.is_unique
    bg=set(full.gene_symbol);deg=full[full.FDR_padj<.05];ora=[]
    for direction in ['all','up','down']:
        hit=set(deg.gene_symbol if direction=='all' else deg.loc[(deg.log2FC>0) if direction=='up' else (deg.log2FC<0),'gene_symbol'])
        for r in eligible.itertuples():
            genes=set(r.genes.split(';'))&bg;hits=genes&hit
            M,K,n,k=len(bg),len(genes),len(hit),len(hits)
            ora.append(dict(set_id=r.set_id,set_name=r.set_name,source=r.source,ROS_related=r.ROS_related,
                direction=direction,N=M,K=K,n=n,k=k,expected=n*K/M,
                fold_enrichment=k*M/(n*K) if n*K else np.nan,pvalue=stats.hypergeom.sf(k-1,M,K,n),hits=';'.join(sorted(hits))))
    ora=pd.DataFrame(ora);ora['FDR_all_sets_directions_per_dataset']=bh(ora.pvalue)
    save(ora.sort_values('pvalue'),dest/'ORA_all_pathways.csv')
    save(ora[ora.ROS_related],dest/'ROS_pathway_ORA.csv')
    # Leading-edge statistics help assess whether ROS signals are broadly distributed.
    details=[]
    for r in gsea[gsea.ROS_related].itertuples():
        for gene in r.leading_edge_genes.split(';'):
            d=full[full.gene_symbol==gene]
            if len(d):details.append(dict(set_id=r.set_id,gene_symbol=gene,log2FC=float(d.iloc[0].log2FC),gene_FDR=float(d.iloc[0].FDR_padj)))
    save(pd.DataFrame(details),dest/'ROS_leading_edge_gene_statistics.csv')
    audit=dict(dataset=acc,ranked_genes=N,eligible_KEGG=int(mask.sum()),eligible_ROS_GO=int((~mask).sum()),
        GSEA_family_size=len(gsea),GSEA_significant=int((gsea.FDR_KEGG_plus_ROS_per_dataset<.05).sum()),
        GSEA_KEGG_significant=int(((gsea.FDR_KEGG_only_per_dataset<.05)&mask).sum()),
        sample_significant=int((gsea.sample_FDR_KEGG_plus_ROS_per_dataset<.05).sum()),
        ORA_family_size=len(ora),ORA_significant=int((ora.FDR_all_sets_directions_per_dataset<.05).sum()),
        random_seed='20261104 + crc32(dataset:K)',null_cache='Recomputable cache excluded from version control',
        interpretation='Exploratory competitive gene-label GSEA; sample-score exact test is a different self-contained null')
    (dest/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)

def summarize():
    allg=[];allo=[];rows=[]
    for acc in DATASETS:
        g=pd.read_csv(OUT/acc/'GSEA_all_pathways.csv');g.insert(0,'dataset',acc);allg.append(g)
        o=pd.read_csv(OUT/acc/'ORA_all_pathways.csv');o.insert(0,'dataset',acc);allo.append(o)
        rows.append(json.loads((OUT/acc/'audit.json').read_text()))
    g=pd.concat(allg,ignore_index=True);o=pd.concat(allo,ignore_index=True)
    g['FDR_all_sets_both_datasets_sensitivity']=bh(g.nominal_p)
    o['FDR_all_sets_directions_both_datasets_sensitivity']=bh(o.pvalue)
    save(g,OUT/'GSEA_both_datasets.csv');save(o,OUT/'ORA_both_datasets.csv')
    save(g[g.ROS_related],OUT/'ROS_pathway_comparison.csv')
    with pd.ExcelWriter(OUT/'Dataset_pathway_GSEA_KEGG.xlsx',engine='openpyxl') as writer:
        g.to_excel(writer,sheet_name='GSEA_all',index=False)
        g[g.ROS_related].to_excel(writer,sheet_name='ROS_GSEA',index=False)
        o.to_excel(writer,sheet_name='ORA_all',index=False)
        for acc in DATASETS:
            pd.read_csv(OUT/acc/'ROS_leading_edge_gene_statistics.csv').to_excel(writer,sheet_name=acc+'_leading_edge',index=False)
        for s in writer.book:s.freeze_panes='C2';s.auto_filter.ref=s.dimensions
    paths=[OUT/'gene_set_catalog.csv',ROOT/'ora_analysis/input/ros_geo/kegg_mouse_gene_names.tsv']
    paths+=list(INP.glob('*'))
    for acc in DATASETS:paths += [OLD/f'{acc}_GEO_full_statistics.csv',OLD/f'{acc}_GSEA_ranked_genes.csv']
    (OUT/'provenance.json').write_text(json.dumps(dict(datasets=rows,
        SHA256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        correction='Primary: each dataset, all eligible KEGG plus 9 pre-specified ROS GO sets. KEGG-only and both-dataset sensitivity also saved.'),indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','analyze','summarize']);parser.add_argument('--dataset',choices=DATASETS);args=parser.parse_args()
    if args.action=='prepare':prepare()
    elif args.action=='summarize':summarize()
    else:
        assert args.dataset;analyze(args.dataset)
