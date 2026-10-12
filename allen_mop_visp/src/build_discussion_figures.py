"""Summarize existing descriptive CSVs; do not fit regional or PBM effects."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

os.environ.setdefault('MPLCONFIGDIR', '/tmp/tpbm-discussion-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures' / 'discussion'
TABLES = ROOT / 'results' / 'discussion'
OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)
COLORS = {'MOp': '#235a9f', 'VISp': '#cc6228'}
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.facecolor': 'white',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
sources = [
    'results/merfish/composition_region_major.csv',
    'results/merfish/composition_layer_major.csv',
    'results/rnaseq/panel22_P58_cohort_cluster_sensitivity.csv',
    'results/connectivity/projection_by_experiment_target_hemisphere.csv',
    'results/connectivity/WT_target_projection_descriptive_summary.csv',
]
outputs = []


def save(fig, stem):
    for ext in ('png', 'pdf', 'svg'):
        path = OUT / f'{stem}.{ext}'
        fig.savefig(path, dpi=190, bbox_inches='tight')
        outputs.append(str(path.relative_to(ROOT)))
    plt.close(fig)


# Percentages retain each region/layer's all-QC-cell denominator.
major = pd.read_csv(ROOT / sources[0])
layer = pd.read_csv(ROOT / sources[1])
assert major.groupby('region').cell_count.sum().to_dict() == {'MOp': 83976, 'VISp': 61884}
assert np.allclose(major.groupby('region').percentage.sum(), 100)
assert set(major.independent_animals) == {1}
categories = ['Glutamatergic neuron', 'GABAergic neuron', 'Astrocyte',
              'Oligodendrocyte', 'OPC', 'Microglia', 'Endothelial',
              'Pericyte/SMC', 'VLMC', 'Other non-neuronal']
fig = plt.figure(figsize=(12.8, 8))
grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.45, wspace=0.40)
ax = fig.add_subplot(grid[:, 0])
y = np.arange(len(categories))
for region, shift in [('MOp', -0.18), ('VISp', 0.18)]:
    rows = major[major.region.eq(region)].set_index('major_group').loc[categories]
    values = rows.percentage.to_numpy()
    ax.barh(y + shift, values, height=.32, color=COLORS[region], label=region)
    for pos, value in zip(y + shift, values):
        ax.text(value + .45, pos, f'{value:.2f}', va='center', fontsize=8)
ax.set(yticks=y, yticklabels=categories, xlim=(0, 79),
       xlabel='% of all QC-annotated cells in region', title='A. Whole-region composition')
ax.invert_yaxis()
ax.legend(frameon=False, loc='lower right')
ax.grid(axis='x', alpha=.18)
layers = ['1', '2/3', '5', '6a', '6b']
layer_rows = []
for row_i, (group, letter) in enumerate([('GABAergic neuron', 'B'), ('Astrocyte', 'C')]):
    ax = fig.add_subplot(grid[row_i, 1])
    for region in ['MOp', 'VISp']:
        rows = layer[layer.region.eq(region) & layer.major_group.eq(group)].set_index('anatomical_layer').loc[layers]
        layer_rows.append(rows.reset_index())
        values = rows.percentage.to_numpy()
        ax.plot(['L' + s for s in layers], values, 'o-', color=COLORS[region], label=region)
    ax.set(title=f'{letter}. {group}: shared anatomical layers',
           ylabel='% of all QC cells within layer', ylim=(0, 16 if row_i == 0 else 42))
    ax.grid(alpha=.18)
    ax.legend(frameon=False, fontsize=9)
fig.suptitle('MOp (M1 atlas region) vs VISp (V1 atlas region): composition', fontsize=15, y=.985)
fig.text(.08, -.015, 'One WT C57BL/6J male brain. MOp n=83,976; VISp n=61,884 cells. Descriptive proportions, not density.\n'
         'Shared CCF layers shown; MOp has no separate CCF L4 label. Whole-region results do not define the illuminated ROI.', fontsize=10)
save(fig, '01_composition_and_depth')
major.to_csv(TABLES / 'composition_whole_region.csv', index=False)
pd.concat(layer_rows).to_csv(TABLES / 'composition_shared_layers.csv', index=False)


# Use the same eight genes and ten adequately sampled neuronal subclasses
# in both panels. Every numeric value is copied from existing P58 summaries.
expression = pd.read_csv(ROOT / sources[2])
genes = ['Entpd1', 'Nt5e', 'Slc29a1', 'Slc29a2', 'Adk', 'Ada', 'Adora1', 'Adora2a']
subclasses = ['Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Vip', 'L2/3 IT CTX',
              'L4/5 IT CTX', 'L5 IT CTX', 'L6 CT CTX', 'L6 IT CTX']
sel = expression[expression.gene.isin(genes) & expression.subclass.isin(subclasses)].copy()
assert len(sel) == 80 and not sel.duplicated(['subclass', 'gene']).any()
assert (sel[['n_cells_MOp', 'n_cells_VISp']] >= 20).all().all()
metrics = ['P58_cell_weighted_MOp_minus_VISp',
           'P58_shared_cluster_equal_weight_MOp_minus_VISp']
matrices = [sel.pivot(index='subclass', columns='gene', values=m).loc[subclasses, genes].to_numpy() for m in metrics]
assert all(np.isfinite(a).all() for a in matrices)
scale = max(float(np.abs(a).max()) for a in matrices)
fig, axes = plt.subplots(1, 2, figsize=(15.5, 7), gridspec_kw={'wspace': .48})
counts = sel.drop_duplicates('subclass').set_index('subclass')
rowlabels = [f"{s.replace(' CTX', '')}  ({int(counts.loc[s, 'n_cells_MOp'])}/{int(counts.loc[s, 'n_cells_VISp'])})" for s in subclasses]
for ax, values, title in zip(axes, matrices, ['A. Cell-weighted subclass mean difference', 'B. Equal weighting of shared clusters']):
    im = ax.imshow(values, cmap='RdBu_r', vmin=-scale, vmax=scale, aspect='auto')
    ax.set(xticks=np.arange(8), xticklabels=genes,
           yticks=np.arange(10), yticklabels=rowlabels, title=title)
    ax.tick_params(axis='x', rotation=45)
    for i in range(10):
        for j in range(8):
            ax.text(j, i, f'{values[i, j]:+.3f}', ha='center', va='center',
                    color='white' if abs(values[i, j]) > scale*.62 else '#111', fontsize=8)
fig.subplots_adjust(right=.90, bottom=.19, top=.83)
cb_ax = fig.add_axes([.92, .25, .014, .50])
fig.colorbar(im, cax=cb_ax, label='MOp - VISp mean log1p(CP10k); not log2FC')
fig.suptitle('Baseline adenosine-related mRNA: direction depends on cell type and weighting', fontsize=15, y=.96)
fig.text(.08, .035, 'P58 male Snap25-IRES2-Cre;Ai14 reporter. Row n = full subclass cells (MOp/VISp), before shared-cluster restriction; donors: 1/2.\n'
         'Red = MOp higher; blue = VISp higher. Preparation date is confounded with region; no regional P/FDR or PBM effect.\n'
         'Ten selected neuronal subclasses; sparse non-neuronal groups excluded from this figure. Low mRNA detection is not biological absence.', fontsize=10)
save(fig, '02_adenosine_gene_sensitivity')
sel.to_csv(TABLES / 'adenosine_neuronal_plot_values.csv', index=False)
expression[expression.gene.isin(genes) & expression.subclass.eq('Astro')].to_csv(TABLES / 'adenosine_astro_sparse_audit.csv', index=False)


# Show every injection; medians are descriptive, not confidence intervals.
projection = pd.read_csv(ROOT / sources[3])
wt = projection[projection.cohort.eq('WT_C57BL6J') & projection.hemisphere_id.eq(3)
                & projection.target_acronym.isin(['ZI', 'RT'])].copy()
assert len(wt) == 82
scopes = [wt, wt[wt.injection_age_P56_only & (wt.injection_source_parent_fraction >= .8)]]
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), sharey='row')
validation_medians = []
for col, data in enumerate(scopes):
    for row, target in enumerate(['ZI', 'RT']):
        ax = axes[row, col]
        labels = []
        for x, region in enumerate(['MOp', 'VISp']):
            selected = data[data.target_acronym.eq(target) & data.source_region.eq(region)].sort_values('experiment_id')
            values = selected.normalized_projection_volume.to_numpy()
            assert (values > 0).all()
            jitter = np.linspace(-.15, .15, len(values))
            ax.scatter(x + jitter, values, s=35, alpha=.70, color=COLORS[region])
            median = float(np.median(values))
            ax.plot([x-.24, x+.24], [median, median], color=COLORS[region], lw=3)
            labels.append(f'{region}\nn={len(values)}; median={median:.4f}')
            validation_medians.append({'scope': ['all', 'P56_purity_0.8'][col], 'target': target,
                                       'region': region, 'n': len(values), 'median_NPV': median})
        ax.set(xticks=[0, 1], xticklabels=labels, xlim=(-.5, 1.5),
               title=f"{target}: {['all primary-source injections', 'P56 + source fraction >=0.8'][col]}")
        ax.set_yscale('log')
        ax.grid(axis='y', alpha=.18)
        if col == 0:
            ax.set_ylabel('Bilateral normalized projection volume\n(log scale; target signal / injection signal)')
fig.suptitle('ZI is an anatomical candidate; RT shows a less distinct source pattern', fontsize=15, y=.98)
fig.subplots_adjust(hspace=.48, bottom=.17, top=.90, wspace=.20)
fig.text(.08, .035, 'WT C57BL/6J male tracer experiments. Each dot = one injection; horizontal line = median; no regional significance test.\n'
         'Source fraction >=0.8 is an exploratory spread filter. Viral construct/protocol differences remain.\n'
         'Axonal signal does not identify postsynaptic Lhx6 cells, synapses, functional sign, NREM effects or PBM responses.', fontsize=10)
save(fig, '03_projection_individual_experiments')
pd.concat([d.assign(plot_scope=['all', 'P56_purity_0.8'][i]) for i, d in enumerate(scopes)]).to_csv(TABLES / 'projection_plot_values.csv', index=False)
medians = pd.DataFrame(validation_medians)
medians.to_csv(TABLES / 'projection_plot_medians.csv', index=False)
summary = pd.read_csv(ROOT / sources[4])
for record in validation_medians:
    scope = ('all_primary_source_injections' if record['scope'] == 'all'
             else 'exact_P56_and_source_parent_signal_fraction_at_least_0.8_exploratory')
    expected = summary[summary.source_region.eq(record['region']) & summary.target_acronym.eq(record['target']) & summary.injection_scope.eq(scope)]
    assert len(expected) == 1, (scope, expected.injection_scope.tolist())
    assert int(expected.iloc[0].n_experiments) == record['n']
    assert np.isclose(expected.iloc[0].median_NPV, record['median_NPV'])

manifest = {
    'source_repository_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(),
    'scope': 'Replot existing descriptive result CSVs; no new PBM analysis, DE, module score, animal test or causal inference.',
    'source_files': [{'path': p, 'sha256': hashlib.sha256((ROOT / p).read_bytes()).hexdigest()} for p in sources],
    'figures': outputs,
    'checks': {'composition_counts_and_denominators': 'passed',
               'expression_unique_80_gene_subclass_rows_and_finite_values': 'passed',
               'projection_medians_and_n_match_original_summary': 'passed'},
    'expression_selection': {'genes': genes, 'subclasses': subclasses,
                             'sparse_astro_in_separate_audit_table': True},
    'projection_medians': validation_medians,
}
(TABLES / 'figure_provenance_and_validation.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps({'figures': outputs, 'validation': manifest['checks']}, indent=2))
