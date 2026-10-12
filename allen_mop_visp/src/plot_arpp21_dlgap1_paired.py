"""Rearrange two existing exploratory dotplot columns into adjacent region rows."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/tpbm-paired-dotplot-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
source = root / 'results/rnaseq/exploratory_cohort_region_subclass_gene_summary.csv'
frame = pd.read_csv(source)
genes = ['Arpp21', 'Dlgap1']
subclasses = ['Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Vip', 'L2/3 IT CTX',
              'L4/5 IT CTX', 'L5 IT CTX', 'L6 CT CTX', 'L6 IT CTX']
data = frame[frame.gene.isin(genes) & frame.subclass.isin(subclasses)].copy()
assert len(data) == 40 and not data.duplicated(['region', 'subclass', 'gene']).any()
assert data.n_cells.ge(20).all()
ordered = [(s, r) for s in subclasses for r in ['MOp', 'VISp']]
index = data.set_index(['subclass', 'region', 'gene'])
figure_dir = root / 'figures/discussion'
result_dir = root / 'results/discussion'
figure_dir.mkdir(parents=True, exist_ok=True)
result_dir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(8.4, 9.3))
norm = Normalize(vmin=0, vmax=float(data.mean_log1p_CP10k.max()))
for i in range(len(subclasses)):
    if i % 2 == 0:
        ax.axhspan(2*i-.5, 2*i+1.5, color='#f1f3f5', zorder=0)
    ax.axhline(2*i+1.5, color='#d5d8dd', lw=.6, zorder=0)
for x, gene in enumerate(genes):
    rows = [index.loc[(s, r, gene)] for s, r in ordered]
    means = np.array([row.mean_log1p_CP10k for row in rows])
    fractions = np.array([row.detection_fraction for row in rows])
    points = ax.scatter(np.full(20, x), np.arange(20), s=4+fractions*160,
                        c=means, cmap='viridis', norm=norm, edgecolors='#333', lw=.35)
    for y, (mean, fraction) in enumerate(zip(means, fractions)):
        ax.text(x+.16, y, f'{mean:.3f} | {fraction:.1%}', va='center', fontsize=9)
labels = [f"{s.replace(' CTX','')} | {r} (n={int(index.loc[(s,r,genes[0])].n_cells):,})" for s, r in ordered]
ax.set(xticks=[0, 1], xticklabels=genes, yticks=np.arange(20), yticklabels=labels,
       xlim=(-.3, 1.95), ylim=(19.7, -.7))
ax.tick_params(axis='x', length=0)
ax.set_title('Adjacent MOp/VISp rows within each neuronal subclass\n'
             'Color = mean expression; size = detection fraction; labels = mean | detected', fontsize=11, pad=15)
cb = fig.colorbar(points, ax=ax, fraction=.045, pad=.035)
cb.set_label('Mean log1p(CP10k); common scale for these two genes')
fig.subplots_adjust(left=.40, right=.85, bottom=.15, top=.89)
fig.text(.04, .045, 'Existing exploratory P57-P59 cohort; age not exactly matched. Male Snap25;Ai14, 10x v2.\n'
         'n = cells, not animals. Region/preparation date confounded; no regional P/FDR or PBM treatment contrast.\n'
         'Arpp21: original PBM common-down candidate. Dlgap1: original PBM common-up candidate.\n'
         'This figure shows baseline atlas expression, copied unchanged from the original result CSV.', fontsize=9)
exports=[]
for ext in ['png', 'pdf', 'svg']:
    path = figure_dir / ('04_arpp21_dlgap1_paired_dotplot.'+ext)
    fig.savefig(path, dpi=190, bbox_inches='tight', facecolor='white')
    exports.append(str(path.relative_to(root)))
plt.close(fig)
data.to_csv(result_dir/'arpp21_dlgap1_paired_dotplot_values.csv', index=False)
manifest = {'source': str(source.relative_to(root)), 'source_SHA256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'scope': 'Replot existing means and detection fractions; no new expression analysis or hypothesis test.',
            'genes': genes, 'subclasses': subclasses, 'numeric_rows_copied': len(data), 'exports': exports,
            'color_scale': 'Shared across both genes in this figure; independent of original up34/down18 scales.',
            'validation': 'Unique region/subclass/gene rows; all selected cell counts >=20; exact existing values copied.'}
(result_dir/'arpp21_dlgap1_paired_dotplot_provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps(manifest, indent=2))
