"""ROS DEG annotation, cross-study comparison, and focused pathway statistics."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ros-mpl-cache')
import hashlib
import itertools
import json
import numpy as np
import pandas as pd
from scipy import stats, optimize, integrate
from fetch_sources import TERMS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INP = ROOT / 'ora_analysis/input/ros_geo'
OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)
SETS = {**TERMS, 'kegg_glutathione': 'mmu00480'}
CANONICAL = ['Sod1', 'Sod2', 'Sod3', 'Cat'] + [f'Gpx{i}' for i in range(1, 9)] + [f'Prdx{i}' for i in range(1, 7)]
LABELS = {
    'ros_metabolism': 'ROS metabolism', 'antioxidant_activity': 'Antioxidant activity',
    'oxidative_stress_response': 'Oxidative stress response', 'glutathione_metabolism': 'Glutathione metabolism (GO)',
    'ros_biosynthesis': 'ROS biosynthesis', 'hydrogen_peroxide_catabolism': 'Hydrogen peroxide removal',
    'superoxide_removal': 'Superoxide removal', 'ros_metabolism_regulation': 'ROS metabolism regulation',
    'nadph_regeneration': 'NADPH regeneration', 'kegg_glutathione': 'Glutathione metabolism (KEGG)',
}

def bh(p):
    p = np.asarray(p, dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    return stats.false_discovery_control(p, method='bh')

def save(df, name):
    df.to_csv(OUT / f'{name}.csv', index=False)
    return df

def valid_symbols(df):
    s = df.gene_symbol
    return s.notna() & ~s.isin(['--', '-', '', 'NA'])

def annotations():
    annotations, sets = [], {}
    for name in TERMS:
        a = pd.read_csv(INP / f'quickgo_{name}.tsv', sep='\t')
        assert set(a['TAXON ID']) == {10090}
        a = a[~a.QUALIFIER.fillna('').str.contains('NOT')].copy()
        a['set_name'] = name
        a['root_id'] = TERMS[name]
        annotations.append(a)
        sets[name] = set(a.SYMBOL.dropna())
    ann = pd.concat(annotations, ignore_index=True)
    save(ann, 'quickgo_annotation_evidence')
    # KEGG list/mmu has four columns. Only explicit symbols before ';' are used.
    name_map = {}
    for line in (INP / 'kegg_mouse_gene_names.tsv').read_text().splitlines():
        cols = line.split('\t')
        if ';' in cols[-1]:
            name_map[cols[0]] = cols[-1].split(';')[0].split(', ')[0]
    kegg_rows = []
    for line in (INP / 'kegg_glutathione_mouse.tsv').read_text().splitlines():
        _, gid = line.split('\t')
        kegg_rows.append({'set_name': 'kegg_glutathione', 'root_id': 'mmu00480',
            'entrez': gid.split(':')[1], 'gene_symbol': name_map.get(gid)})
    kegg = save(pd.DataFrame(kegg_rows), 'kegg_gene_mapping')
    sets['kegg_glutathione'] = set(kegg.gene_symbol.dropna())
    catalog = pd.DataFrame({'gene_symbol': sorted(set.union(*sets.values()))})
    for name, genes in sets.items():
        catalog[name] = catalog.gene_symbol.isin(genes)
    catalog['canonical_scavenger_enzyme'] = catalog.gene_symbol.isin(CANONICAL)
    catalog['set_memberships'] = catalog.apply(lambda r: ';'.join(n for n in sets if r[n]), axis=1)
    def group(r):
        if r['canonical_scavenger_enzyme']: return 'Canonical scavenger enzyme'
        if r['antioxidant_activity'] or r['hydrogen_peroxide_catabolism'] or r['superoxide_removal']:
            return 'Other antioxidant / removal'
        if r['ros_biosynthesis']: return 'ROS biosynthesis'
        if r['glutathione_metabolism'] or r['kegg_glutathione'] or r['nadph_regeneration']:
            return 'Glutathione / NADPH support'
        if r['ros_metabolism']: return 'Other ROS metabolism'
        if r['ros_metabolism_regulation']: return 'ROS regulation'
        return 'Oxidative stress response'
    catalog['display_group'] = catalog.apply(group, axis=1)
    save(catalog, 'ros_gene_catalog')
    save(pd.DataFrame([dict(set_name=n, label=LABELS[n], root_id=SETS[n], mouse_symbols=len(g)) for n, g in sets.items()]), 'gene_set_definitions')
    return sets, catalog

def supplied_deg(catalog):
    selected = {}
    for acc, filename in [('GSE245658', 'selected_hippocampus_PBM'), ('GSE325653', 'selected_cortex_PBM')]:
        d = pd.read_csv(ROOT / f'ora_analysis/deg_validation/{filename}.csv')
        d = d[valid_symbols(d)].copy()
        duplicated = d.gene_symbol.duplicated(keep=False)
        save(d[duplicated], f'{acc}_duplicate_symbol_audit')
        d = d[~duplicated]
        assert (d.FDR_padj < .05).all()
        assert (np.where(d.log2FC > 0, 'up', 'down') == d.direction).all()
        selected[acc] = d
        ros = d.merge(catalog, on='gene_symbol', validate='one_to_one')
        save(ros.sort_values(['display_group', 'FDR_padj']), f'{acc}_ROS_DEG')
    a, b = selected.values()
    both = a.merge(b, on='gene_symbol', suffixes=('_A', '_B'), validate='one_to_one')
    both = both.merge(catalog, on='gene_symbol', validate='one_to_one')
    both['concordant'] = both.direction_A == both.direction_B
    save(both, 'ROS_DEG_both_studies')
    save(both[both.concordant], 'ROS_same_direction_DEG')
    save(both[~both.concordant], 'ROS_opposite_direction_DEG')
    # Missing from selected list is an explicit unknown, never an imputed zero.
    union = sorted(set(a.gene_symbol) | set(b.gene_symbol))
    heat = catalog[catalog.gene_symbol.isin(union)].copy()
    for acc, d in selected.items():
        heat = heat.merge(d[['gene_symbol', 'log2FC', 'FDR_padj']].rename(columns={
            'log2FC': f'{acc}_log2FC', 'FDR_padj': f'{acc}_FDR'}), on='gene_symbol', how='left', validate='one_to_one')
    save(heat, 'ROS_DEG_union')
    return selected

def legacy_ora():
    d = pd.read_csv(ROOT / 'ora_analysis/input/ORA_fullcatalog_KEGG_GO.csv')
    mask = d.set_name.str.contains('reactive oxygen|superoxide|hydrogen peroxide|antioxidant|glutathione|thioredoxin|peroxidase|oxidative stress', case=False)
    chosen = d[mask].copy()
    # Verify original global-catalog statistics before reporting the ROS subset.
    for _, group in d.groupby(['dataset', 'direction']):
        p = stats.hypergeom.sf(group.k.to_numpy() - 1, group.n_background_genes.to_numpy(),
            group.set_size_tested_in_background.to_numpy(), group.n_DEG_in_contrast.to_numpy())
        assert np.allclose(p, group.pval, atol=1e-12, rtol=1e-7)
        assert np.allclose(bh(p), group.padj_BH_full_catalog, atol=1e-12, rtol=1e-7)
    save(chosen, 'legacy_ROS_ORA_global_FDR')
    return chosen

def full_statistics():
    full, audit = {}, []
    for acc in ['GSE245658', 'GSE325653']:
        d = pd.read_csv(OUT / f'{acc}_GEO_full_statistics.csv')
        valid = valid_symbols(d)
        dup = valid & d.gene_symbol.duplicated(keep=False)
        save(d[dup], f'{acc}_full_stats_duplicate_symbols')
        finite = np.isfinite(d.pvalue) & np.isfinite(d.log2FC)
        keep = valid & ~dup & finite
        s = d[keep].copy()
        assert s.gene_symbol.is_unique
        assert ((s.pvalue >= 0) & (s.pvalue <= 1)).all()
        full[acc] = s.set_index('gene_symbol', drop=False).rename_axis(None)
        audit.append(dict(dataset=acc, source_rows=len(d), ambiguous_symbol_rows=int(dup.sum()),
            invalid_symbol_rows=int((~valid).sum()), finite_unique_ranked_genes=int(keep.sum()),
            finite_adjusted_p_genes=int(np.isfinite(s.FDR_padj).sum()),
            DEG_FDR05=int((s.FDR_padj < .05).sum())))
    save(pd.DataFrame(audit), 'full_statistics_audit')
    return full

def gene_panels(full, selected, catalog):
    support = ['Gclc', 'Gclm', 'Gss', 'Gsr', 'Txn1', 'Txn2', 'Txnrd1', 'Txnrd2',
        'G6pdx', 'Pgd', 'Idh1', 'Idh2', 'Me1', 'Nnt', 'Nfe2l2', 'Keap1', 'Ppargc1a', 'Txnip']
    panel = pd.DataFrame({'gene_symbol': CANONICAL + support})
    panel['panel_class'] = np.where(panel.gene_symbol.isin(CANONICAL), 'Canonical scavenger enzyme', 'Context / support (selected reference panel)')
    # Reference panel does not expand the tested ontology gene sets.
    for acc, d in full.items():
        panel = panel.merge(d[['gene_symbol', 'log2FC', 'FDR_padj']], on='gene_symbol', how='left').rename(columns={
            'log2FC': f'{acc}_GEO_log2FC', 'FDR_padj': f'{acc}_GEO_FDR'})
        panel[f'{acc}_provided_DEG'] = panel.gene_symbol.isin(selected[acc].gene_symbol)
    save(panel, 'canonical_and_support_reference_panel')
    # Assess compatibility, without merging independently refitted results into tasks 1/2.
    rows = []
    for acc, d in selected.items():
        geo_ros = full[acc].loc[full[acc].FDR_padj < .05].reset_index(drop=True).merge(catalog, on='gene_symbol', validate='one_to_one')
        geo_ros['in_provided_DEG_table'] = geo_ros.gene_symbol.isin(d.gene_symbol)
        save(geo_ros, f'{acc}_GEO_ROS_DEG')
        extra = full[acc].loc[(full[acc].FDR_padj < .05) & ~full[acc].gene_symbol.isin(d.gene_symbol)].reset_index(drop=True)
        extra['ROS_catalog_member'] = extra.gene_symbol.isin(catalog.gene_symbol)
        save(extra, f'{acc}_GEO_extra_DEG_vs_provided')
        m = d.merge(full[acc][['log2FC', 'FDR_padj']], left_on='gene_symbol', right_index=True, suffixes=('_provided', '_GEO'))
        save(m, f'{acc}_provided_vs_GEO_statistics')
        rows.append(dict(dataset=acc, supplied_unique_genes=len(d), matched=len(m),
            log2FC_spearman=stats.spearmanr(m.log2FC_provided, m.log2FC_GEO).statistic,
            sign_agreement=float((np.sign(m.log2FC_provided) == np.sign(m.log2FC_GEO)).mean()),
            GEO_FDR05=int((m.FDR_padj_GEO < .05).sum()),
            max_abs_log2FC_difference=float(abs(m.log2FC_provided - m.log2FC_GEO).max()),
            max_abs_FDR_difference=float(abs(m.FDR_padj_provided - m.FDR_padj_GEO).max())))
    save(pd.DataFrame(rows), 'provided_vs_GEO_audit')
    ga = pd.read_csv(OUT / 'GSE245658_GEO_ROS_DEG.csv')
    gb = pd.read_csv(OUT / 'GSE325653_GEO_ROS_DEG.csv')
    shared = ga[['gene_symbol', 'log2FC', 'FDR_padj']].merge(gb[['gene_symbol','log2FC','FDR_padj']], on='gene_symbol', suffixes=('_A_refit','_B_deposited'))
    shared['same_direction'] = np.sign(shared.log2FC_A_refit) == np.sign(shared.log2FC_B_deposited)
    save(shared, 'GEO_ROS_shared_DEG_sensitivity')

def focused_ora(full, sets):
    rows = []
    for acc, d in full.items():
        # The ORA universe includes finite q values because selection uses q<.05;
        # genes removed by independent filtering cannot be selected as DEGs.
        eligible = d[np.isfinite(d.FDR_padj)]
        bg = set(eligible.gene_symbol)
        deg = eligible[eligible.FDR_padj < .05]
        for direction in ['all', 'up', 'down']:
            hit = set(deg.gene_symbol if direction == 'all' else deg.loc[(deg.log2FC > 0) if direction == 'up' else (deg.log2FC < 0), 'gene_symbol'])
            for name, genes in sets.items():
                kset, hits = bg & genes, hit & genes
                n, N, K, k = len(hit), len(bg), len(kset), len(hits)
                p = stats.hypergeom.sf(k - 1, N, K, n)
                rows.append(dict(dataset=acc, direction=direction, set_name=name, root_id=SETS[name],
                    N=N, K=K, n=n, k=k, expected=n*K/N, fold_enrichment=(k*N/(n*K) if n*K else np.nan),
                    pvalue=p, hits=';'.join(sorted(hits))))
    ora = pd.DataFrame(rows)
    ora['FDR_BH_60_tests'] = bh(ora.pvalue)
    save(ora, 'focused_GEO_ORA')
    return ora

def ranked_enrichment(full, sets):
    results = []
    for acc, d in full.items():
        # A common signed-z transform puts A and B on a stated ranking definition;
        # cortex has no deposited Wald statistic. No q threshold or FC cutoff.
        z = np.sign(d.log2FC) * stats.norm.isf(np.maximum(d.pvalue.to_numpy(), 1e-300) / 2)
        ranks = pd.DataFrame({'gene_symbol': d.gene_symbol.to_numpy(), 'signed_z': z})
        ranks = ranks.sort_values(['signed_z', 'gene_symbol'], ascending=[False, True])
        save(ranks, f'{acc}_GSEA_ranked_genes')
        # Streaming gene-label null: O(N + K) memory, rather than storing a
        # transcriptome x 10,000-permutation matrix. Weighted GSEA, exponent=1.
        symbols = ranks.gene_symbol.to_numpy()
        weights = abs(ranks.signed_z.to_numpy())
        rng = np.random.default_rng(20261102)
        N = len(symbols)
        for name, genes in sets.items():
            positions = np.flatnonzero(np.isin(symbols, list(genes)))
            K = len(positions)
            if not 5 <= K <= min(2000, N-1): continue
            es, boundary = enrichment_score(positions, weights)
            null = np.empty(10000)
            for j in range(len(null)):
                pos = np.sort(rng.choice(N, K, replace=False))
                null[j] = enrichment_score(pos, weights)[0]
            sign_null = null[null >= 0] if es >= 0 else null[null < 0]
            p = (1 + np.sum(abs(sign_null) >= abs(es))) / (1 + len(sign_null))
            nes = es / np.mean(abs(sign_null))
            lead = positions[positions <= boundary] if es >= 0 else positions[positions >= boundary]
            results.append(dict(dataset=acc, Term=name, ES=es, NES=nes,
                **{'NOM p-val': p}, matched_genes=K, background_genes=N,
                null_permutations=10000, same_sign_null_permutations=len(sign_null),
                Lead_genes=';'.join(symbols[lead])))
    gsea = pd.DataFrame(results)
    # Primary correction covers both studies and all pre-specified sets that passed size filtering.
    # Conservative +1 floor avoids zero permutation p values being interpreted as exact zero.
    gsea['permutation_p_plus1'] = gsea['NOM p-val']
    gsea['FDR_BH_across_studies'] = bh(gsea.permutation_p_plus1)
    save(gsea, 'focused_GSEA')
    return gsea

def enrichment_score(positions, weights):
    """Weighted KS running-sum extrema, evaluated just before and after hits."""
    N, K = len(weights), len(positions)
    hit = weights[positions]
    hit = hit / hit.sum() if hit.sum() else np.full(K, 1/K)
    misses = (positions - np.arange(K)) / (N-K)
    after = np.cumsum(hit) - misses
    before = after - hit
    hi, lo = int(np.argmax(after)), int(np.argmin(before))
    # Resolve theoretical equal-magnitude extrema toward the positive tail,
    # rather than allowing roundoff at machine precision to choose the sign.
    return (float(after[hi]), int(positions[hi])) if after[hi] >= abs(before[lo])-1e-12 else (float(before[lo]), int(positions[lo]))

def sample_tests(full, sets):
    rows, scores = [], []
    for acc in full:
        fname = 'GSE245658_log2_normalized_counts.csv.gz' if acc == 'GSE245658' else 'GSE325653_log2_FPKM.csv.gz'
        expr = pd.read_csv(OUT / fname, index_col=0)
        if acc == 'GSE245658':
            mapping = full[acc].set_index('gene_id').gene_symbol
            expr = expr.loc[expr.index.intersection(mapping.index)]
            expr.index = mapping.loc[expr.index].values
        else:
            expr = expr.loc[expr.index.intersection(full[acc].index)]
        assert expr.index.is_unique
        n = expr.shape[1] // 2
        assert expr.shape[1] == n * 2
        std = expr.std(axis=1, ddof=1)
        expr = expr.loc[std > 1e-8]
        z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1), axis=0)
        assignments = list(itertools.combinations(range(n*2), n))
        for name, genes in sets.items():
            members = z.index.intersection(sorted(genes))
            if len(members) < 5: continue
            score = z.loc[members].mean(axis=0).to_numpy()
            obs = score[n:].mean() - score[:n].mean()
            differences = []
            for treated in assignments:
                mask = np.zeros(n*2, dtype=bool); mask[list(treated)] = True
                differences.append(score[mask].mean() - score[~mask].mean())
            p = np.mean(np.abs(differences) >= abs(obs) - 1e-12)
            pooled = np.sqrt((score[n:].var(ddof=1) + score[:n].var(ddof=1)) / 2)
            rows.append(dict(dataset=acc, set_name=name, genes=len(members), n_control=n, n_PBM=n,
                score_difference=obs, standardized_difference=obs/pooled if pooled else np.nan,
                exact_two_sided_p=p, permutations=len(assignments), members=';'.join(members)))
            for s, value in zip(expr.columns, score):
                scores.append(dict(dataset=acc, set_name=name, sample=s, condition='PBM' if s in expr.columns[n:] else 'control', score=value))
    tests = pd.DataFrame(rows)
    tests['FDR_BH_across_studies'] = bh(tests.exact_two_sided_p)
    save(tests, 'sample_level_exact_permutation')
    save(pd.DataFrame(scores), 'sample_level_pathway_scores')
    return tests

def sensitivity():
    rows = []
    for n in [4, 5]:
        for alpha in [.05, .05/20]:
            df = 2*n - 2
            critical = stats.t.ppf(1-alpha/2, df)
            def power(d):
                nc = d * np.sqrt(n/2)
                # Integrate the normal/chi-square definition directly; avoids
                # noncentral-t tail NaNs in some installed SciPy versions.
                def integrand(v):
                    bound = critical * np.sqrt(v/df)
                    return (stats.norm.sf(bound-nc) + stats.norm.cdf(-bound-nc)) * stats.chi2.pdf(v, df)
                return integrate.quad(integrand, 0, np.inf, epsabs=1e-9)[0]
            d80 = optimize.brentq(lambda d: power(d)-.8, .01, 8)
            rows.append(dict(n_per_group=n, alpha=alpha, d_for_80pct_power=d80,
                assumed_test='Two-sided equal-variance normal two-sample t-test on a single pathway score',
                note='Design sensitivity only; not achieved power, not exact permutation power, not a gene-count effect'))
    save(pd.DataFrame(rows), 'design_sensitivity_assumptions')

def main():
    sets, catalog = annotations()
    selected = supplied_deg(catalog)
    legacy_ora()
    full = full_statistics()
    gene_panels(full, selected, catalog)
    ora = focused_ora(full, sets)
    gsea = ranked_enrichment(full, sets)
    perm = sample_tests(full, sets)
    sensitivity()
    # Convenient workbook, raw full transcriptome statistics remain CSVs.
    names = ['gene_set_definitions', 'GSE245658_ROS_DEG', 'GSE325653_ROS_DEG',
        'ROS_same_direction_DEG', 'ROS_opposite_direction_DEG', 'ROS_DEG_union',
        'GSE245658_GEO_ROS_DEG', 'GSE325653_GEO_ROS_DEG', 'GEO_ROS_shared_DEG_sensitivity',
        'GSE325653_GEO_extra_DEG_vs_provided',
        'canonical_and_support_reference_panel', 'legacy_ROS_ORA_global_FDR',
        'focused_GEO_ORA', 'focused_GSEA', 'sample_level_exact_permutation',
        'full_statistics_audit', 'provided_vs_GEO_audit', 'design_sensitivity_assumptions']
    with pd.ExcelWriter(OUT / 'ROS_scavenger_analysis.xlsx', engine='openpyxl') as writer:
        for name in names:
            pd.read_csv(OUT / f'{name}.csv').to_excel(writer, sheet_name=name[:31], index=False)
        for sheet in writer.book:
            sheet.freeze_panes = 'B2'
            sheet.auto_filter.ref = sheet.dimensions
            for col in sheet.columns:
                sheet.column_dimensions[col[0].column_letter].width = min(48, max(14, max(len(str(c.value or '')) for c in col[:100]) + 2))
    import importlib.metadata
    versions = {name: importlib.metadata.version(name) for name in ['numpy','pandas','scipy','matplotlib','pydeseq2','openpyxl','libarchive-c']}
    (HERE / 'requirements.txt').write_text('\n'.join(f'{n}=={v}' for n,v in versions.items()) + '\n')
    inputs = list(INP.glob('*')) + [ROOT / 'ora_analysis/input/DEG_FDR05_GSE245658_GSE325653.csv', ROOT / 'ora_analysis/input/ORA_fullcatalog_KEGG_GO.csv']
    summary = dict(branch='ROS_scavenger', gene_set_count=len(sets), versions=versions,
        broad_ROS_DEG={a:len(pd.read_csv(OUT/f'{a}_ROS_DEG.csv')) for a in full},
        same_direction=len(pd.read_csv(OUT/'ROS_same_direction_DEG.csv')),
        opposite_direction=len(pd.read_csv(OUT/'ROS_opposite_direction_DEG.csv')),
        significant_new_ORA=int((ora.FDR_BH_60_tests<.05).sum()),
        significant_GSEA_BH=int((gsea.FDR_BH_across_studies<.05).sum()),
        significant_sample_permutation_BH=int((perm.FDR_BH_across_studies<.05).sum()),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs if p.is_file()})
    (OUT / 'validation_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'input_sha256'}, indent=2), flush=True)

if __name__ == '__main__': main()
