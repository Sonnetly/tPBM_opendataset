"""Recover complete statistics and sample expression from the original GEO files.

The hippocampal refit is a NEW two-group analysis, kept separate from supplied DEG.
The cortex uses deposited DESeq2 statistics, never DESeq2 on FPKM.
"""
from pathlib import Path
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ros-mpl-cache')
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[2]
INP = ROOT / 'ora_analysis/input/ros_geo'
OUT = Path(__file__).resolve().parent / 'results'
OUT.mkdir(exist_ok=True)

def main():
    samples = []
    for acc in ['GSE245658', 'GSE325653']:
        for block in (INP / f'{acc}_samples.txt').read_text().split('^SAMPLE = ')[1:]:
            lines = block.splitlines()
            def values(key):
                return '; '.join(l.split(' = ', 1)[1] for l in lines if l.startswith(key + ' = '))
            samples.append(dict(dataset=acc, gsm=lines[0], title=values('!Sample_title'),
                source=values('!Sample_source_name_ch1'), characteristics=values('!Sample_characteristics_ch1'),
                processing=values('!Sample_data_processing')))
    pd.DataFrame(samples).to_csv(OUT / 'geo_sample_metadata.csv', index=False)

    a = pd.read_csv(INP / 'GSE245658_all.genes.expression.annot.LT.txt.gz', sep='\t')
    controls = [f'5xFAD{i}_HP_count' for i in range(1, 5)]
    treated = [f'5xFAD+L{i}_HP_count' for i in range(1, 5)]
    columns = controls + treated
    counts = a.set_index('id')[columns].T
    assert not counts.columns.duplicated().any()
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy() >= 0).all()
    fractional_count_fraction = float(((counts.to_numpy() % 1) != 0).mean())
    # Submitter count columns contain estimated fractional counts. Explicit rounding
    # is required by PyDESeq2; this approximation is audited, not called raw counts.
    counts = counts.round().astype(int)
    # Reproducible, treatment-blind expression prefilter: >=10 reads in >=4/8 samples.
    keep = (counts >= 10).sum(axis=0) >= 4
    counts = counts.loc[:, keep].astype(int)
    meta = pd.DataFrame({'condition': ['control'] * 4 + ['PBM'] * 4}, index=columns)
    dds = DeseqDataSet(counts=counts, metadata=meta, design='~condition', n_cpus=2, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=['condition', 'PBM', 'control'], alpha=0.05, n_cpus=2, quiet=True)
    ds.summary()
    result = ds.results_df.rename_axis('gene_id').reset_index().rename(columns={
        'log2FoldChange': 'log2FC', 'stat': 'wald_stat', 'padj': 'FDR_padj'})
    result = result.merge(a[['id', 'Symbol']].rename(columns={'id': 'gene_id', 'Symbol': 'gene_symbol'}), on='gene_id', validate='one_to_one')
    result.to_csv(OUT / 'GSE245658_GEO_full_statistics.csv', index=False)
    expr = pd.DataFrame(np.log2(dds.layers['normed_counts'].T + 1), index=counts.columns, columns=counts.index)
    expr.index.name = 'gene_id'
    expr.to_csv(OUT / 'GSE245658_log2_normalized_counts.csv.gz')
    meta.to_csv(OUT / 'GSE245658_analysis_samples.csv', index_label='sample')

    b = pd.read_excel(INP / 'GSE325653_raw_counts.xlsx')
    result_b = b[['gene_id', 'BaseMean_DESeq2', 'log2FoldChange', 'p-value', 'q-value']].rename(columns={
        'gene_id': 'gene_symbol', 'BaseMean_DESeq2': 'baseMean', 'log2FoldChange': 'log2FC',
        'p-value': 'pvalue', 'q-value': 'FDR_padj'})
    result_b.to_csv(OUT / 'GSE325653_GEO_full_statistics.csv', index=False)
    columns_b = [f'FPKM_Hypoxia{i:02}' for i in range(1, 6)] + [f'FPKM_Laser{i:02}' for i in range(1, 6)]
    expr_b = b.set_index('gene_id')[columns_b]
    assert not expr_b.index.duplicated().any()
    assert np.isfinite(expr_b.to_numpy()).all() and (expr_b.to_numpy() >= 0).all()
    np.log2(expr_b + 1).to_csv(OUT / 'GSE325653_log2_FPKM.csv.gz', index_label='gene_symbol')
    pd.DataFrame({'sample': columns_b, 'condition': ['control'] * 5 + ['PBM'] * 5}).to_csv(OUT / 'GSE325653_analysis_samples.csv', index=False)
    fractional = ((expr_b.to_numpy() % 1) != 0).mean()
    audit = dict(hippocampus=dict(original_rows=len(a), prefilter_genes=int(keep.sum()),
        control_n=4, PBM_n=4, design='~condition; PBM vs control; WT excluded',
        prefilter='rounded estimated count >=10 in >=4 of 8 samples', independent_filtering=True,
        fractional_count_fraction=fractional_count_fraction,
        count_handling='Nearest-integer rounding (numpy/pandas half-to-even); no length offsets available; exploratory refit',
        note='New PyDESeq2 fit, not the supplied FDR05 table; do not substitute in tasks 1/2'),
        cortex=dict(original_rows=len(b), control_n=5, PBM_n=5,
        source='Deposited full DESeq2 statistics; log2(Laser/Hypoxia)',
        expression_type='FPKM, despite raw_counts.xlsx filename', fractional_FPKM_fraction=float(fractional),
        raw_counts_available_in_this_file=False))
    (OUT / 'geo_data_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2), flush=True)

if __name__ == '__main__':
    main()
