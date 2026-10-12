"""Replot all six existing 10x dotplots from result CSVs with adjacent region rows.

Does not recalculate expression, combine age strata, impute unavailable groups,
or perform regional hypothesis tests.
"""
from pathlib import Path
import hashlib
import json
import math
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/tpbm-all-paired-dotplots-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / 'figures/rnaseq/paired'
RESULTS = ROOT / 'results/rnaseq/paired_dotplots'
FIGURES.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)
REGIONS = ['MOp', 'VISp']
CLASSES = ['GABAergic', 'Glutamatergic', 'Non-Neuronal']
SOURCES = ['results/rnaseq/exploratory_cohort_region_subclass_gene_summary.csv',
           'results/rnaseq/region_subclass_gene_summary.csv',
           'figures/rnaseq/dotplot_strata.csv']
exploratory = pd.read_csv(ROOT / SOURCES[0])
strict = pd.read_csv(ROOT / SOURCES[1])
strata = pd.read_csv(ROOT / SOURCES[2])
exploratory['source_row_index'] = exploratory.index
strict['source_row_index'] = strict.index
plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
maximum = max(exploratory.mean_log1p_CP10k.max(), strict.mean_expression.max())
COLOR_MAX = math.ceil(float(maximum) * 100) / 100
NORM = Normalize(vmin=0, vmax=COLOR_MAX)
jobs = []
for name, selection in [('panel22', exploratory.panel22),
                        ('common_up34', exploratory.common52_direction.eq('up')),
                        ('common_down18', exploratory.common52_direction.eq('down'))]:
    frame = exploratory.loc[selection].copy()
    assert len(frame.gene.unique()) == {'panel22': 22, 'common_up34': 34, 'common_down18': 18}[name]
    jobs.append({'stem': f'exploratory_cohort_{name}', 'frame': frame,
                 'source': SOURCES[0], 'original': f'exploratory_cohort_{name}_dotplot.png',
                 'title': f'{name}: baseline atlas | P57-P59 exploratory cohort',
                 'condition': 'Male Snap25;Ai14, 10x v2; ages not exactly matched.'})
for row in strata.itertuples(index=False):
    selected = pd.Series(True, index=strict.index)
    for field in ['chemistry', 'genotype', 'sex', 'age', 'facs']:
        selected &= strict[field].eq(getattr(row, field))
    frame = strict.loc[selected & strict.in_panel22].copy()
    frame = frame.rename(columns={'mean_expression': 'mean_log1p_CP10k'})
    assert frame.gene.nunique() == 22
    jobs.append({'stem': Path(row.filename).stem + '_' + row.age, 'frame': frame,
                 'source': SOURCES[1], 'original': row.filename,
                 'title': f'panel22: baseline atlas | exact-age {row.age} stratum',
                 'condition': f'Male Snap25;Ai14, 10x v2, RFP-positive; age {row.age}.'})
assert len(jobs) == 6


def prepare(job):
    data = job['frame']
    assert not data.duplicated(['region', 'class', 'subclass', 'gene']).any()
    genes = sorted(data.gene.unique())
    groups = sorted(set(zip(data['class'], data.subclass)),
                    key=lambda x: (CLASSES.index(x[0]), x[1]))
    rows = [(cls, sub, region) for cls, sub in groups for region in REGIONS]
    index = data.set_index(['class', 'subclass', 'region', 'gene'])
    counts = data.groupby(['class', 'subclass', 'region']).n_cells.agg(['first', 'nunique'])
    assert counts['nunique'].eq(1).all()
    assert data.groupby(['class', 'subclass', 'region']).gene.nunique().eq(len(genes)).all()
    table = []
    for pos, (cls, sub, region) in enumerate(rows):
        for gene in genes:
            key = (cls, sub, region, gene)
            observed = key in index.index
            record = index.loc[key] if observed else None
            table.append({'plot_name': job['stem'], 'row_order': pos,
                          'class': cls, 'subclass': sub, 'region': region, 'gene': gene,
                          'observed': observed, 'n_cells': int(record.n_cells) if observed else None,
                          'mean_log1p_CP10k': float(record.mean_log1p_CP10k) if observed else None,
                          'detection_fraction': float(record.detection_fraction) if observed else None,
                          'source_file': job['source'],
                          'source_row_index': int(record.source_row_index) if observed else None})
    values = pd.DataFrame(table)
    assert values.loc[values.observed, 'detection_fraction'].between(0, 1).all()
    assert np.isfinite(values.loc[values.observed, ['mean_log1p_CP10k', 'detection_fraction']]).all().all()
    assert values.loc[~values.observed, ['mean_log1p_CP10k', 'detection_fraction', 'n_cells']].isna().all().all()
    return genes, rows, counts, values


def draw(job, genes, rows, counts, values, numeric=False, page_label=''):
    cell_width = .88 if numeric else .43
    left, right = 4.65, 1.95
    plot_width = max(3.6, len(genes) * cell_width)
    plot_height = len(rows) * (.34 if numeric else .28)
    bottom, top = 1.95, 1.10
    width, height = left + plot_width + right, bottom + plot_height + top
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([left/width, bottom/height, plot_width/width, plot_height/height])
    for pair in range(len(rows)//2):
        if pair % 2 == 0:
            ax.axhspan(pair*2-.5, pair*2+1.5, color='#f1f3f5', zorder=0)
        ax.axhline(pair*2+1.5, color='#d4d8dd', lw=.5, zorder=0)
    subset = values[values.gene.isin(genes) & values.observed]
    xmap = {g: i for i, g in enumerate(genes)}
    x = subset.gene.map(xmap).to_numpy(dtype=float)
    if numeric:
        x -= .20
    points = ax.scatter(x, subset.row_order, s=4 + 120*subset.detection_fraction,
                        c=subset.mean_log1p_CP10k, cmap='viridis', norm=NORM,
                        edgecolors='#444' if numeric else 'none', linewidths=.25)
    if numeric:
        for record in subset.itertuples(index=False):
            ax.text(xmap[record.gene]+.03, record.row_order,
                    f'{record.mean_log1p_CP10k:.3f}\n{record.detection_fraction:.1%}',
                    ha='left', va='center', fontsize=6.5, linespacing=1.1)
    labels = []
    missing_rows = []
    for pos, (cls, sub, region) in enumerate(rows):
        group_key = (cls, sub, region)
        if group_key in counts.index:
            n = int(counts.loc[group_key, 'first'])
            labels.append(f'{cls} | {sub} | {region} (n={n:,})' + (' [small n]' if n < 10 else ''))
        else:
            labels.append(f'{cls} | {sub} | {region} (N/A)')
            missing_rows.append(pos)
            ax.axhspan(pos-.45, pos+.45, facecolor='#ececec', hatch='///',
                       edgecolor='#d6d6d6', linewidth=0, zorder=.5)
            ax.text((len(genes)-1)/2, pos, 'N/A: no cells in this age/subclass',
                    ha='center', va='center', fontsize=7.5, color='#555')
    ax.set(xticks=range(len(genes)), xticklabels=genes,
           yticks=range(len(rows)), yticklabels=labels,
           xlim=(-.55, len(genes)-.35), ylim=(len(rows)-.5, -.5))
    ax.tick_params(axis='x', rotation=65, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_ha('right')
    ax.tick_params(axis='y', length=0, labelsize=7.7)
    for label, key in zip(ax.get_yticklabels(), rows):
        region = key[2]
        label.set_color('#235a9f' if region == 'MOp' else '#b55b24')
        if key not in counts.index or int(counts.loc[key, 'first']) < 10:
            label.set_color('#777')
    color_height = min(5.0, plot_height*.60)
    cax = fig.add_axes([(left+plot_width+.25)/width,
                       (bottom+plot_height-color_height)/height, .16/width, color_height/height])
    cb = fig.colorbar(points, cax=cax)
    cb.set_label('Mean log1p(CP10k); shared 0-4.25 scale', fontsize=8)
    cb.ax.tick_params(labelsize=8)
    legend_x = left+plot_width+.28
    legend_y = max(bottom+.25, bottom+plot_height-color_height-1.8)
    legend_ax = fig.add_axes([legend_x/width, legend_y/height, 1.45/width, 1.45/height])
    legend_ax.axis('off')
    for fraction in [.25, .5, 1.0]:
        legend_ax.scatter([], [], s=4+120*fraction, c='#777', label=f'{fraction:.0%} detected')
    legend_ax.legend(loc='upper left', frameon=False, title='Dot area', fontsize=7.5, title_fontsize=8)
    fig.text(.025, 1-.27/height, job['title'] + (' | ' + page_label if page_label else ''),
             fontsize=12, va='top')
    fig.text(.025, 1-.57/height,
             'Adjacent MOp/VISp rows within each subclass; color = mean expression; dot area = detection fraction.'
             + (' Numbers: mean / detected %.' if numeric else ''), fontsize=9, va='top')
    fig.text(.025, .12/height,
             job['condition']+'\n'
             'n = cells, not animals. [small n] = fewer than 10 cells, retained for completeness. N/A is unobserved, not zero expression.\n'
             'Region/preparation date confounded; descriptive baseline expression, no regional P/FDR or PBM effect.\n'
             'up34/down18 are candidate directions in the original PBM studies; values and detection rates are copied from existing CSVs.',
             fontsize=8.5, va='bottom')
    return fig


manifest = {'status': 'all_six_existing_dotplots_replotted',
            'analysis_scope': 'Replot existing result tables only; no expression recalculation or regional tests.',
            'shared_color_limits': [0, COLOR_MAX],
            'dot_area_formula': '4 + 120 * detection_fraction',
            'numeric_gene_columns_per_page': 8,
            'missing_group_policy': 'Retain adjacent placeholder row with N/A and missing numeric values; never impute zero.',
            'small_group_policy': 'Retain all original groups; n<10 explicitly labeled, no new exclusion cutoff.',
            'sources': [{'path': p, 'SHA256': hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in SOURCES],
            'plots': []}
all_values = []
overview_bundle = FIGURES / 'all_six_paired_dotplots.pdf'
numeric_bundle = FIGURES / 'all_six_numeric_detail_pages.pdf'
with PdfPages(overview_bundle) as overview_pdf, PdfPages(numeric_bundle) as numeric_pdf:
    for job in jobs:
        genes, rows, counts, values = prepare(job)
        all_values.append(values)
        values.to_csv(RESULTS/(job['stem']+'_plot_values.csv'), index=False)
        overview = draw(job, genes, rows, counts, values)
        exports = []
        for ext in ['png', 'pdf', 'svg']:
            path = FIGURES/(job['stem']+'_paired_dotplot.'+ext)
            overview.savefig(path, dpi=180, bbox_inches='tight', facecolor='white')
            exports.append(str(path.relative_to(ROOT)))
        overview_pdf.savefig(overview, bbox_inches='tight', facecolor='white')
        plt.close(overview)
        numeric_path = FIGURES/(job['stem']+'_numeric_detail.pdf')
        gene_pages = [genes[i:i+8] for i in range(0, len(genes), 8)]
        with PdfPages(numeric_path) as panel_pdf:
            for page, page_genes in enumerate(gene_pages, 1):
                detail = draw(job, page_genes, rows, counts, values, numeric=True,
                              page_label=f'genes {page}/{len(gene_pages)}')
                panel_pdf.savefig(detail, bbox_inches='tight', facecolor='white')
                numeric_pdf.savefig(detail, bbox_inches='tight', facecolor='white')
                if job['stem'] == 'panel22_stratum_001_P58' and page == 1:
                    detail.savefig(FIGURES/'P58_panel22_numeric_page01_preview.png', dpi=180,
                                   bbox_inches='tight', facecolor='white')
                plt.close(detail)
        manifest['plots'].append({'name': job['stem'], 'original_figure': job['original'],
                                  'gene_count': len(genes), 'row_count_including_NA': len(rows),
                                  'observed_group_count': len(counts),
                                  'missing_group_count': len(rows)-len(counts),
                                  'observed_value_count': int(values.observed.sum()),
                                  'missing_value_count': int((~values.observed).sum()),
                                  'region_cell_totals': job['frame'].drop_duplicates(['region','class','subclass']).groupby('region').n_cells.sum().to_dict(),
                                  'overview_exports': exports,
                                  'numeric_detail_pdf': str(numeric_path.relative_to(ROOT)),
                                  'numeric_detail_page_count': len(gene_pages),
                                  'plot_values_csv': str((RESULTS/(job['stem']+'_plot_values.csv')).relative_to(ROOT))})
pd.concat(all_values, ignore_index=True).to_csv(RESULTS/'all_plot_values_and_availability.csv', index=False)
manifest['bundles'] = [str(overview_bundle.relative_to(ROOT)), str(numeric_bundle.relative_to(ROOT))]
(RESULTS/'replot_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
lines = ['# RNA-seq: MOp–VISp를 인접 행으로 비교하는 dotplot', '',
         '기존 6개 dotplot의 수치를 그대로 재배치했다. 같은 class/subclass의 MOp 행 바로 아래에 VISp 행을 놓는다.', '',
         '- 색: 평균 log1p(CP10k). 모든 그림에 공통 0–4.25 범위.',
         '- 점 면적: 검출률(raw count > 0). 세포 수는 각 행의 n.',
         '- 숫자판: 한 페이지에 유전자 최대 8개, 각 점 옆에 평균·검출률을 표시.',
         '- N/A: 해당 연령/세포 유형에 관측 세포가 없음. 0 발현과 구분.',
         '- `[small n]`: 표집 세포가 10개 미만인 행. 기존 값을 보존하고 불안정성을 표시.', '',
         '**P58만 두 영역이 모두 표집되었다. P57은 MOp만, P59는 VISp만 있어 해당 연령에서 직접 지역 비교가 불가능하다.**',
         'P57–P59 탐색 묶음은 지역별 연령 구성이 다르고 모든 자료에 region–library prep 교란이 남는다.',
         '`up34/down18`은 이전 PBM 연구의 후보 방향이며 현재 atlas의 Region 차이나 PBM 처리 효과가 아니다.', '',
         '[전체 6개 요약 그림 PDF](all_six_paired_dotplots.pdf) · [전체 숫자판 PDF](all_six_numeric_detail_pages.pdf)', '',
         '| 묶음 | 요약 PNG | 요약 PDF | 숫자판 PDF |', '|---|---|---|---|']
for item in manifest['plots']:
    name = item['name']
    lines.append(f'| {name} | [PNG]({name}_paired_dotplot.png) | [PDF]({name}_paired_dotplot.pdf) | [숫자판]({name}_numeric_detail.pdf) |')
lines += ['', '[원수치·관측 여부 CSV](../../../results/rnaseq/paired_dotplots/all_plot_values_and_availability.csv)',
          '[출처·패널 구성 manifest](../../../results/rnaseq/paired_dotplots/replot_manifest.json)',
          '[재현 코드](../../../src/replot_all_rnaseq_dotplots_paired.py)', '',
          '저장소 루트에서 `python3 allen_mop_visp/src/replot_all_rnaseq_dotplots_paired.py`로 재생성한다.', '']
(FIGURES/'README_ko.md').write_text('\n'.join(lines), encoding='utf-8')
print(json.dumps({'plots': len(jobs), 'color_max': COLOR_MAX,
                  'overview_pages': len(jobs),
                  'numeric_pages': sum(p['numeric_detail_page_count'] for p in manifest['plots']),
                  'manifest': str((RESULTS/'replot_manifest.json').relative_to(ROOT))}, indent=2))
