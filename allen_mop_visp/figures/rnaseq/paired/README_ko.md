# RNA-seq: MOp–VISp를 인접 행으로 비교하는 dotplot

기존 6개 dotplot의 수치를 그대로 재배치했다. 같은 class/subclass의 MOp 행 바로 아래에 VISp 행을 놓는다.

- 색: 평균 log1p(CP10k). 모든 그림에 공통 0–4.25 범위.
- 점 면적: 검출률(raw count > 0). 세포 수는 각 행의 n.
- 숫자판: 한 페이지에 유전자 최대 8개, 각 점 옆에 평균·검출률을 표시.
- N/A: 해당 연령/세포 유형에 관측 세포가 없음. 0 발현과 구분.
- `[small n]`: 표집 세포가 10개 미만인 행. 기존 값을 보존하고 불안정성을 표시.

**P58만 두 영역이 모두 표집되었다. P57은 MOp만, P59는 VISp만 있어 해당 연령에서 직접 지역 비교가 불가능하다.**
P57–P59 탐색 묶음은 지역별 연령 구성이 다르고 모든 자료에 region–library prep 교란이 남는다.
`up34/down18`은 이전 PBM 연구의 후보 방향이며 현재 atlas의 Region 차이나 PBM 처리 효과가 아니다.

[전체 6개 요약 그림 PDF](all_six_paired_dotplots.pdf) · [전체 숫자판 PDF](all_six_numeric_detail_pages.pdf)

| 묶음 | 요약 PNG | 요약 PDF | 숫자판 PDF |
|---|---|---|---|
| exploratory_cohort_panel22 | [PNG](exploratory_cohort_panel22_paired_dotplot.png) | [PDF](exploratory_cohort_panel22_paired_dotplot.pdf) | [숫자판](exploratory_cohort_panel22_numeric_detail.pdf) |
| exploratory_cohort_common_up34 | [PNG](exploratory_cohort_common_up34_paired_dotplot.png) | [PDF](exploratory_cohort_common_up34_paired_dotplot.pdf) | [숫자판](exploratory_cohort_common_up34_numeric_detail.pdf) |
| exploratory_cohort_common_down18 | [PNG](exploratory_cohort_common_down18_paired_dotplot.png) | [PDF](exploratory_cohort_common_down18_paired_dotplot.pdf) | [숫자판](exploratory_cohort_common_down18_numeric_detail.pdf) |
| panel22_stratum_000_P57 | [PNG](panel22_stratum_000_P57_paired_dotplot.png) | [PDF](panel22_stratum_000_P57_paired_dotplot.pdf) | [숫자판](panel22_stratum_000_P57_numeric_detail.pdf) |
| panel22_stratum_001_P58 | [PNG](panel22_stratum_001_P58_paired_dotplot.png) | [PDF](panel22_stratum_001_P58_paired_dotplot.pdf) | [숫자판](panel22_stratum_001_P58_numeric_detail.pdf) |
| panel22_stratum_002_P59 | [PNG](panel22_stratum_002_P59_paired_dotplot.png) | [PDF](panel22_stratum_002_P59_paired_dotplot.pdf) | [숫자판](panel22_stratum_002_P59_numeric_detail.pdf) |

[원수치·관측 여부 CSV](../../../results/rnaseq/paired_dotplots/all_plot_values_and_availability.csv)
[출처·패널 구성 manifest](../../../results/rnaseq/paired_dotplots/replot_manifest.json)
[재현 코드](../../../src/replot_all_rnaseq_dotplots_paired.py)

저장소 루트에서 `python3 allen_mop_visp/src/replot_all_rnaseq_dotplots_paired.py`로 재생성한다.
