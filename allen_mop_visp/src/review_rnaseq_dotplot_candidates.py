"""Review existing descriptive gene/subclass contrasts without DE or P/FDR."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/rnaseq/paired_dotplots'
OUT.mkdir(parents=True, exist_ok=True)
sources = ['results/rnaseq/all62_P58_cohort_cluster_sensitivity.csv',
           'results/rnaseq/shared_cluster_region_difference.csv']
all_values = pd.read_csv(ROOT / sources[0])
clusters = pd.read_csv(ROOT / sources[1])
clusters = clusters[clusters.age.eq('P58')].copy()
A = 'P58_cell_weighted_MOp_minus_VISp'
B = 'P58_shared_cluster_equal_weight_MOp_minus_VISp'
C = 'P57_59_shared_cluster_equal_weight_MOp_minus_VISp'
assert not all_values.duplicated(['gene', 'class', 'subclass']).any()
assert not clusters.duplicated(['gene', 'class', 'subclass', 'cluster']).any()
summaries = []
for key, group in clusters.groupby(['class', 'subclass', 'gene']):
    eligible = group[group.n_cells_MOp.ge(10) & group.n_cells_VISp.ge(10)]
    summaries.append(dict(zip(['class', 'subclass', 'gene'], key),
                          equal_cluster_mean_recomputed=group.MOp_minus_VISp.mean(),
                          shared_clusters_ge10_both_regions=len(eligible),
                          exploratory_cluster_ge10_difference=eligible.MOp_minus_VISp.mean(),
                          positive_shared_cluster_count=int(group.MOp_minus_VISp.gt(0).sum()),
                          negative_shared_cluster_count=int(group.MOp_minus_VISp.lt(0).sum())))
review = all_values.merge(pd.DataFrame(summaries), on=['class','subclass','gene'], how='left', validate='one_to_one')
available = review[B].notna()
assert np.allclose(review.loc[available, B], review.loc[available, 'equal_cluster_mean_recomputed'], atol=1e-12)
review['detection_difference_percentage_points'] = 100*(review.detection_fraction_MOp-review.detection_fraction_VISp)
review['same_direction_A_B_cohort'] = np.sign(review[A]).eq(np.sign(review[B])) & np.sign(review[B]).eq(np.sign(review[C]))
review['same_direction_B_ge10'] = np.sign(review[B]).eq(np.sign(review.exploratory_cluster_ge10_difference))
review.to_csv(OUT/'all62_descriptive_candidate_review.csv', index=False)
selections = [
    ('Adenosine hypothesis', 'Adora1', 'Pvalb'),
    ('Adenosine hypothesis', 'Adora1', 'L4/5 IT CTX'),
    ('Adenosine hypothesis', 'Adk', 'L4/5 IT CTX'),
    ('Metabolic hypothesis', 'Ppargc1a', 'L4/5 IT CTX'),
    ('Larger observed regional contrasts', 'Zeb2', 'L4/5 IT CTX'),
    ('Larger observed regional contrasts', 'Grid1', 'L4/5 IT CTX'),
    ('Larger observed regional contrasts', 'Tnr', 'L4/5 IT CTX'),
    ('Larger observed regional contrasts', 'Sphkap', 'L6b CTX'),
    ('Larger observed regional contrasts', 'Stxbp6', 'L5/6 NP CTX'),
    ('Physiology-linked candidates', 'Cacna1c', 'L6 CT CTX'),
    ('Physiology-linked candidates', 'Grin2a', 'Pvalb'),
    ('Physiology-linked candidates', 'Penk', 'Vip'),
    ('Previously discussed dotplot columns', 'Dlgap1', 'L5 IT CTX'),
    ('Previously discussed dotplot columns', 'Arpp21', 'L2/3 IT CTX'),
    ('Weight/cohort-sensitive examples', 'Adora1', 'L6 CT CTX'),
    ('Weight/cohort-sensitive examples', 'Dlgap1', 'L6 CT CTX'),
    ('Weight/cohort-sensitive examples', 'Kcnj6', 'L6 CT CTX'),
    ('Weight/cohort-sensitive examples', 'Grm2', 'L4/5 IT CTX'),
    ('Weight/cohort-sensitive examples', 'Adk', 'L6 IT CTX'),
]
index = review.set_index(['gene','subclass'])
selected = index.loc[[(gene,sub) for _,gene,sub in selections]].reset_index()
selected.insert(0, 'review_category', [category for category,_,_ in selections])
selected.to_csv(OUT/'selected_gene_subclass_observations.csv', index=False)
neuronal = review[review['class'].isin(['GABAergic','Glutamatergic']) & review.n_cells_MOp.ge(100) & review.n_cells_VISp.ge(100)]
overview = {'neuronal_subclasses_both_regions_n_ge100': neuronal[['class','subclass']].drop_duplicates().shape[0],
            'gene_subclass_comparisons': len(neuronal),
            'median_absolute_P58_mean_difference': float(neuronal[A].abs().median())}
lines = ['# RNA-seq dotplot에서 우선 확인할 gene × 세포 유형', '',
         '**모든 유전자가 같은 것은 아니다. 아래는 기존 P58 결과에서 고른 기술적 관측 후보이며 통계적으로 유의한 지역 DEG나 최종 qPCR 순위가 아니다.**', '',
         'P58는 두 영역의 연령이 같지만 region–prep date 교란은 남는다. P57–P59 탐색 PDF와 P58 평균은 표본 범위가 달라 수치가 다를 수 있다.', '',
         '## 선별 방식', '',
         '기존 62개 후보 전부의 평균·검출률과 공유 cluster 민감도를 확인했다. Adenosine 가설 관련성과 상대적으로 큰 지역 대비를 함께 고려해 예시를 선택했다. 모든 칸의 단위는 평균 log1p(CP10k)이며 차이는 MOp−VISp다.', '',
         '- A: P58의 전체 subclass 세포 평균 차이.',
         '- B: P58 공통 cluster에 동일 가중치를 준 차이.',
         '- C: P57–P59 탐색 cohort의 공통 cluster 동일 가중치 차이. P58 세포를 포함하므로 독립 재현 자료가 아니다.',
         '- 추가 민감도: P58 공통 cluster 중 양쪽 세포 수가 각각 10 이상인 cluster만 다시 같은 비중으로 집계. 탐색적으로 선택한 범위이며 official QC나 검정 기준이 아니다.', '',
         '## 우선 볼 관측', '',
         '| Gene | Subclass | n MOp/VISp | P58 평균 MOp/VISp | 검출률 MOp/VISp | A | B | C |',
         '|---|---|---:|---:|---:|---:|---:|---:|']
for row in selected.iloc[:14].itertuples(index=False):
    lines.append(f'| {row.gene} | {row.subclass} | {int(row.n_cells_MOp)}/{int(row.n_cells_VISp)} | '
                 f'{row.mean_expression_MOp:.3f}/{row.mean_expression_VISp:.3f} | '
                 f'{100*row.detection_fraction_MOp:.1f}%/{100*row.detection_fraction_VISp:.1f}% | '
                 f'{getattr(row,A):+.3f} | {getattr(row,B):+.3f} | {getattr(row,C):+.3f} |')
lines += ['',
          '연구 질문의 우선 후보는 Adora1/Pvalb와 Adora1/IT의 반대 방향, Adk/L4/5 IT, Ppargc1a/IT다. ADK는 adenosine을 AMP로 전환하는 대사 효소이며 ADK 발현이 높다고 adenosine 순축적이 크다고 할 수 없다. Ppargc1a의 VISp 우세도 실제 ATP 생산이나 PBM 처리 반응을 측정한 값이 아니다.',
          'Zeb2/Grid1/Tnr/L4/5 IT, Sphkap/L6b, Stxbp6/L5/6 NP는 62개 후보 중 비교적 큰 관측 대비를 보여주는 추가 후보다. Zeb2는 세포 정체·지역 전사 프로그램, Grid1/Tnr는 시냅스·세포외 환경, Sphkap/Stxbp6는 신호·시냅스 조절 축의 후보로 살핀다. 기능적 활성이나 수면 인과성은 미확인이다.',
          'Adora1/Pvalb는 공유 cluster 13개 중 12개에서 MOp 쪽 차이가 양수지만, 이를 독립 동물 반복으로 세지 않는다. Dlgap1/L5 IT와 Arpp21/L2/3 IT는 비교적 작은 차이가 여러 요약 방식에서 같은 방향으로 남는다.', '',
          '## 단순 지역 후보로 해석하기 어려운 예', '',
          '| Gene | Subclass | A | B | C |', '|---|---|---:|---:|---:|']
for row in selected.iloc[14:].itertuples(index=False):
    lines.append(f'| {row.gene} | {row.subclass} | {getattr(row,A):+.3f} | {getattr(row,B):+.3f} | {getattr(row,C):+.3f} |')
lines += ['',
          '방향 반전이나 차이 축소는 공유 cluster 선택·가중치·연령 범위에 민감하다는 뜻이다. C가 거의 0이면 반대 지역의 안정한 우세로 확대 해석하지 않는다. 희소 cluster 제외 민감도도 batch 교란을 해결하지 않는다.', '',
          '## 왜 그림은 비슷해 보이는가', '',
          f"두 영역에 각각 세포 100개 이상인 neuronal subclass {overview['neuronal_subclasses_both_regions_n_ge100']}개에서 62개 후보의 {overview['gene_subclass_comparisons']}개 대비를 보면, 평균 차이 절댓값 중앙값은 {overview['median_absolute_P58_mean_difference']:.3f}다. 이는 보기 위한 탐색적 표본 제한에서의 기술 요약이며 생물학적 동등성을 증명하지 않는다.",
          '대부분의 공유 세포 유형에서 기본 발현 패턴이 유사하고, 공통 색 범위 0–4.25에서는 작은 차이가 잘 드러나지 않는다. 검출률이 양쪽 약100%이면 점 크기도 비슷하다. 숫자판과 실제 대비를 함께 확인한다.', '',
          '[전체 후보 검토 CSV](all62_descriptive_candidate_review.csv) · [선택 관측 CSV](selected_gene_subclass_observations.csv) · [출처·검증](candidate_review_provenance.json)',
          '[인접 행 dotplot 목록](../../../figures/rnaseq/paired/README_ko.md)',
          '[재현 코드](../../../src/review_rnaseq_dotplot_candidates.py)', '']
(OUT/'gene_observation_review_ko.md').write_text('\n'.join(lines), encoding='utf-8')
manifest = {'scope':'Descriptive review of existing contrasts; no DE, animal inference, P/FDR or definitive qPCR ranking.',
            'sources':[{'path':p,'SHA256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in sources],
            'selected_observations':len(selected), 'overview':overview,
            'exploratory_extra_sensitivity':'Equal P58 shared-cluster weights after requiring n>=10 in each region; not official QC.',
            'checks':{'duplicate_keys_absent':True,'P58_cluster_means_match_existing_B':True}}
(OUT/'candidate_review_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(manifest,indent=2))
