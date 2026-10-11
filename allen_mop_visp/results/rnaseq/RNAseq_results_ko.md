# CTX-HPF 10x MOp–VISp 실제 후보 발현 분석

공식 GEO GSE185862 원본에서 **MOp 123,727개와 VISp 29,763개, 총 153,490개 세포의 62개 후보 유전자**를 실제로 분석했다. 처음의 접속 실패 이후 Allen 설명 페이지와 GEO HTTPS archive가 열렸다. Allen 페이지가 연결하는 공식 GSE185862 archive의 자료를 사용했으며 후보 22개와 공통 52개의 합집합 62개가 모두 원본 matrix에 있다.

31,053개 원본 유전자의 세포별 raw-count 합을 분모로 CP10k를 계산한 뒤 `log1p`했다. 후보 62개만 더한 값을 정규화 분모로 사용하지 않았다. 결과는 관찰된 세포 평균의 기술적 비교이다. **영역별 library preparation date가 완전히 분리되어 생물학적 지역 효과의 P/FDR 또는 PBM 반응을 추정하지 않았다.**

## 표본과 직접 비교 범위

| 조건 | MOp | VISp |
|---|---:|---:|
| QC/taxonomy에 있는 선택 세포 | 123,727 | 29,763 |
| GABAergic 표집 세포 | 27,373 | 4,905 |
| Glutamatergic 표집 세포 | 95,154 | 24,486 |
| Non-Neuronal 표집 세포 | 1,200 | 372 |
| 기록된 donor ID | 3 | 3 |
| 나이 | P57, P58 | P58, P59 |
| P58 세포 | 43,829 | 16,822 |
| P58 기록된 donor ID | 1 | 2 |
| Library preparation date label | 171026 | 171120, 180115 |

양쪽 모두 `10x_v2`, 수컷, `Snap25-IRES2-Cre/wt;Ai14(RCL-tdT)/wt`, `RFP-positive` 선별이다. GEO protocol에서 C57BL/6J background와 Chromium Single Cell 3′ Reagent Kit v2를 확인했다. **WT B6J 자료가 아니라 B6J 배경 reporter 자료**이다. RNA-seq 표집 세포 수나 glia 비율을 조직의 실제 cellular composition으로 해석하지 않는다. 모든 선택 세포의 cortical layer field는 원본에서 `NA`이며 transcriptomic subclass를 해부학적 layer로 바꾸지 않는다.

원본 XLSX에서 donor ID마다 고유한 age를 확인해 메타데이터에 연결했다. MOp와 VISp donor는 겹치지 않으며 원본 protocol은 같은 ROI의 sample pooling을 기술한다. 독립 biological pool 관계는 확인하지 않았다. 특히 P58 MOp는 기록된 donor가 한 명이므로 많은 세포를 독립 동물로 세지 않는다.

## 나이 일치 비교와 전체 cohort 탐색

**P58 나이 일치 결과**는 `within_subclass_region_difference.csv`의 `both_regions_observed=True` 행이다. chemistry, genotype, sex, age, FACS가 같은 strata끼리만 대응한다. 직접 대응되는 22개 subclass × 62 genes, 1,364개 기술적 대조가 있다. P57/P59는 한 영역만 있어 같은 유형의 지역 비교로 연결하지 않는다. `shared_cluster_standardized_difference.csv`는 알려진 공유 cluster만 남겨 같은 cluster 가중치를 준다.

**P57–P59 전체 cohort 탐색**은 `exploratory_cohort_*` 파일이다. 원본 age를 유지하면서 같은 chemistry/genotype/sex/FACS 안에서 age가 완전히 일치하지 않는 전체 표집을 기술적으로 요약한다. P58 결과의 보조 민감도 화면이며 나이 일치 비교로 주장하지 않는다. 전체 cohort도 region–preparation date 교란 때문에 동물 검정을 하지 않는다.

공유 cluster로 좁혀 같은 가중치를 주었을 때 방향이 바뀌면, 세부 유형 선택·가중치가 gross subclass 대비에 기여할 수 있다. 이것만으로 composition의 인과 효과나 batch 제거가 증명되지는 않는다. 공유 cluster가 없으면 미평가로 남기며 방향 반전으로 세지 않는다.

## 핵심 관찰

차이는 **MOp minus VISp 평균 `log1p(CP10k)`**이며 log2FC가 아니다. 아래는 P58 수치이다.

| Gene / subclass | 표집 n MOp / VISp | 검출률 MOp / VISp | subclass 평균 차이 | 공유 cluster 동일 가중치 차이 |
|---|---:|---:|---:|---:|
| Adora1 / L2/3 IT CTX | 4,014 / 2,104 | 55.5% / 72.3% | −0.1370 | −0.0890 |
| Adora1 / Pvalb | 1,576 / 422 | 45.8% / 35.3% | +0.1478 | +0.1392 |
| Adk / L4/5 IT CTX | 13,373 / 6,567 | 65.4% / 54.1% | +0.2211 | +0.0956 |
| Grm2 / L2/3 IT CTX | 4,014 / 2,104 | 21.5% / 36.4% | −0.0555 | −0.0530 |
| Grm2 / L4/5 IT CTX | 13,373 / 6,567 | 48.6% / 52.4% | +0.0282 | −0.0568 |
| Penk / Vip | 2,025 / 882 | 46.7% / 59.4% | −0.1821 | −0.0956 |
| Penk / L5 IT CTX | 4,109 / 620 | 9.3% / 4.2% | +0.0538 | −0.0206 |
| Kcnj6 / L6 CT CTX | 7,790 / 1,948 | 73.4% / 70.0% | +0.1982 | −0.0792 |
| Nt5e / Pvalb | 1,576 / 422 | 10.5% / 19.7% | −0.0457 | −0.0207 |
| Adora2a / L6 CT CTX | 7,790 / 1,948 | 1.16% / 0.62% | +0.0047 | +0.0037 |

**Adora1은 세포 유형에 따라 방향이 다르다.** L2/3 IT에서는 VISp 쪽, Pvalb에서는 MOp 쪽 관찰 평균이 높고, 이 두 방향은 공유 cluster와 전체 P57–P59 탐색에서도 유지된다. A1 receptor를 영역 전체 하나의 수치로만 비교하는 것보다 IT·Pvalb를 분리해 위치·발현을 확인할 근거다. mRNA 차이를 adenosine의 수면 유도 threshold로 환산하지 않는다.

**Adk/L4/5 IT, Grm2/L2/3 IT, Penk/Vip**도 P58과 전체 cohort의 공유 cluster 화면에서 같은 방향을 보이는 후속 위치 확인 후보이다. 반면 `Grm2/L4/5 IT`, `Kcnj6/L6 CT`, `Penk/L5 IT`는 공유 cluster 선택·같은 가중치 적용 후 방향이 바뀐다. gross subclass 평균만으로 지역 차이를 주장하기 어렵다.

**Nt5e는 나이 범위에 민감한 subtype이 있고, Adora2a·Ada·Slc29a2 등은 검출률이 낮다.** Nt5e/L2/3 IT의 subclass 차이는 P58 −0.0124에서 전체 cohort +0.0027로 바뀐다. Adora2a/L6 CT 검출률은 약 1%이고 많은 다른 뉴런 subclass는 더 낮다. 미검출은 수용체 부재나 생물학적 0을 뜻하지 않는다. 기존 PBM의 Nt5e·Adora2a 감소 정보도 정상 atlas의 MOp–VISp 기저 방향을 결정하지 않는다.

최종 qPCR 순위는 배정하지 않았다. 현재 화면은 Adora1의 세포 유형별 위치, Adk의 대사 축, Grm2와 Penk의 해당 subtype 발현을 독립 표본·실제 자극 위치에서 확인할 근거를 제공한다. Nt5e·Adora2a는 SMART-seq 및 RNAscope 등 검출 가능한 측정으로 보완해야 한다. M1 자극 좌표는 확정되었고 V1 AP는 잠정값이다. Bregma 좌표를 CCF로 등록하는 방법과 국소 광분포가 아직 확정되지 않아 실제 자극 ROI 특이 결론은 보류했다.

## 원본과 검증 기록

공식 Allen CTX-HPF 다운로드 entry가 연결하는 `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/` 원본을 받았다. 정상 HTML GEO sample 페이지의 HTTP 200은 reCAPTCHA라 reject하고 공식 SOFT export endpoint에서 accession 본문을 확인했다. TLS 검증은 유지했다.

- 10x matrix: 5,716,426,072 bytes; SHA256 `1ec6e788edfde84a875121762b651931a39f731dc35da1db3089fa82211c2163`.
- 실제 `data/counts`: 31,053 genes × 1,169,320 cells; int32; gzip chunks 1,000 × 10,000. `data/shape`는 반대인 cell × gene이므로 gene/sample identity 배열 길이로 축을 확인했다.
- 원본 taxonomy 메타데이터: 1,169,213 cells. matrix에만 있는 107개는 audit에 기록했으며 선택된 MOp/VISp는 모두 원본 matrix에 정확히 대응했다.
- `ctx_hpf_mop_visp_candidate_counts.h5ad`: 153,490 × 62. `obs.full_library_counts`는 원본 31,053개 유전자 합계다.
- `download_manifest.json`, matrix `.download.json`, `GEO_sample_source_manifest.json`, `subset_manifest.json`: URL·크기·SHA256 및 schema/정규화 근거.

큰 원본은 `input/rnaseq/raw/`에서 보관하고 portable ZIP/Git 산출물에서 제외한다. subset, 선택 메타데이터, donor age 원본 레코드, SOFT, 작은 readme와 provenance는 보존한다. 62개 후보 결과를 재현할 때 큰 원본을 다시 내려받을 필요가 없다.

## 재실행

보존한 subset으로 P58 strict 비교와 전체 cohort 탐색을 다시 계산한다. 프로젝트 위치를 코드 경로로 찾으므로 Windows로 프로젝트 폴더를 옮겨 자신의 Python 환경에서 실행할 수 있다.

```bash
python src/rnaseq/analyze_expression.py \
  --metadata input/rnaseq/selected_cell_metadata.csv.gz \
  --h5ad input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad \
  --value-kind raw_counts \
  --library-total-column full_library_counts
python src/rnaseq/cohort_screening.py
python src/rnaseq/summarize_findings.py
```

큰 원본부터 subset을 재현하려면:

```bash
python src/rnaseq/download_geo.py GSE185862_expression_matrix_10x.hdf5
python src/rnaseq/download_geo.py GSE185862_metadata_10x.csv.gz
python src/rnaseq/download_geo.py GSE185862_CTX-HPF-73538sample_metadata.xlsx
python src/rnaseq/extract_geo_subset.py --workers 2
```

extractor는 이미 보존한 `GSM5629156.soft.txt`·`GSM5629158.soft.txt`에서 chemistry를 확인한다. 원본이 이미 완전한 크기로 저장되어 있으면 downloader는 다시 받지 않고 SHA256을 표시한다. 62개만 있는 subset에 `--full-transcriptome-counts-confirmed`를 사용하면 안 된다. 반드시 전체 원본 library 합계 열을 지정한다.

## 결과 파일

- `panel22_P58_cohort_cluster_sensitivity.csv`: 모든 22개 × 대응 subtype의 exact-age 평균·전체 cohort·공유 cluster 민감도.
- `all62_P58_cohort_cluster_sensitivity.csv`: 공통 52개까지 포함한 전체 탐색 표.
- `panel22_screening_coverage_and_sensitivity.csv`: 검출률·방향 민감도이며 유의성/최종 qPCR 순위가 아님.
- `core_gene_observation_examples.csv`: 위 대표 표의 검증 가능한 원 수치.
- `recorded_donor_subclass_candidate_summaries.csv`: 기록된 donor별 기술 통계; 독립 pool pseudobulk 검정 결과가 아님.
- `figures/rnaseq/panel22_stratum_001.png`: P58 dot plot. 다른 strata의 나이는 `dotplot_strata.csv`에 기록. 그림 행마다 표집 n 표시.
- `figures/rnaseq/exploratory_cohort_panel22_dotplot.{png,pdf}` 및 `common_up34`, `common_down18`: 전체 cohort 그림에 age-unmatched·교란 문구 표시.

`results/rnaseq/status.json`은 `descriptive_expression_complete`이다. SMART-seq는 별도 보조 결과를 사용하며 10x와 정규화 발현값을 합치지 않는다.
