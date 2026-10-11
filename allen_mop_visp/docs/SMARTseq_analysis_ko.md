# SMART-seq MOp–VISp 보조 분석

GEO GSE185862의 SMART-seq v4 실제 발현 자료에서 정확한 `MOp`와 `VISp`를 선택해, PBM 패널 22개와 같은 방향 공통 DEG 52개의 합집합 62개를 분석했다. 62개 모두 원본 행렬에 존재한다. 이 자료는 **C57BL/6J 배경의 transgenic/reporter 자료**이며, 정상 야생형 B6의 직접 비교를 제공하지 않는다.

## 분석 대상과 메타데이터

| 단계 | MOp | VISp | 비교 층화 수 |
|---|---:|---:|---:|
| 원본에서 정확한 지역 선택 | 6,588 | 15,572 | — |
| 알려진 genotype·sex·정확한 age·FACS·subclass, 두 지역 모두 존재 | 2,403 | 2,206 | 67 |
| 위 조건에서 양쪽 각 10 cells 이상 표시 | 2,110 | 1,871 | 38 |
| 같은 실제 dissection-layer 표기까지 일치, 각 10 cells 이상 | 결과 CSV에 기록 | 결과 CSV에 기록 | 22 |

원본 GEO XLSX에서 donor ID별 유일한 보고 연령을 연결했다. 연령이 확인되지 않은 2,077 cells는 엄격한 지역 비교에서 제외했으며, 검출률은 `missing_age_candidate_detection_audit.csv`에 남겼다. Injection의 material/ROI/method가 기록된 1,641 cells도 엄격한 비교에서 제외했다. 빈 injection 필드는 **기록된 injection 없음**으로 취급하며, 실험 이력 부재를 독립적으로 증명하지 않는다. 이 제외 집단은 일부 겹치므로 두 숫자를 단순 합산하면 안 된다.

이 RNA-seq cell 수는 FACS·해부학적 절취·library 선택을 거친 표집 수이며, 조직의 실제 세포 구성 비율로 해석하지 않는다.

`n ≥ 10`은 매우 작은 표집 층화를 그림·요약에서 제외하기 위한 표시 필터이며, 통계적 검정·검출력의 기준이 아니다. 한 cell 이상으로 양쪽이 존재하는 모든 공유 층화 결과를 별도 파일에 보존했다. 38개 표시 층화의 연령은 P53, P54, P55, P56, P59이고 reporter line은 Slc32a1, Snap25, Pvalb, Slc17a7, Nr5a1, Rbp4 계열이다. 각각의 정확한 genotype은 `displayed_stratum_key.csv`에 있다. 각 지역에서 13 donors가 포함되지만 지역 사이에 같은 donor ID나 같은 library batch ID가 없다. 이 자료에서는 지역·batch·dissection의 효과를 분리한 동물 수준 유의성 검정을 수행하지 않았다.

## 원본 계수와 정규화 검증

원본 `GSE185862_expression_matrix_SSv4.hdf5`는 73,363 cells × 45,768 genes의 sparse 행렬이다. `data/exon`은 cell × gene의 CSC이고 `data/t_exon`은 그 transpose이다. 유전자·cell ID의 길이와 sparse pointer를 통해 방향을 확인했다. 메타데이터는 73,347 cells이며, 선택한 22,160 cells는 모두 행렬 ID와 정확하게 한 번씩 대응한다.

GEO의 MOp_L2-3 (`GSM5629189`)와 VISp_L2-3 (`GSM5629176`) 처리 기록은 **uniquely mapped exon read counts**, STAR PCR duplicate 제거, SMART-seq v4 전장 cDNA 증폭을 설명한다. 따라서 해당 값은 UMI 또는 TPM이 아니다. 분석에서 사용하는 정규화는 다음과 같다.

```text
CPM = gene exon reads / cell의 전체 exon reads × 1,000,000
mean_expression = cell별 log1p(CPM)의 평균
```

22,160 cells 각각에 대해 **45,768개 전체 유전자 read counts를 다시 합산**했고, 모두 원본 `data/total_exon_counts`와 정확히 일치했다. 62개 후보만 합산한 값을 분모로 쓰지 않았다. 후보 정수·비음수 검증 및 10개 cells에서 두 sparse 방향 간 후보 계수 일치도 확인했다. 같은 유전자 안에서 지역을 비교하는 상대 read abundance 지표이며, gene length와 assay가 다른 유전자·기술 사이의 절대 발현량을 직접 비교하는 척도가 아니다.

## 실제 결과를 읽는 방법

`shared_strata_region_difference_min10.csv`에서 양수는 해당 메타데이터 층화에서 MOp의 평균 log1p(read CPM)이 높고, 음수는 VISp가 높다는 뜻이다. 이 평균의 차이는 log2 fold change가 아니다. 세포 검출률과 각 지역 cell 수를 함께 확인한다.

38개 표시 층화에서의 방향 패턴은 다음과 같다. 서로 다른 genotype·연령 층화를 독립 반복으로 취급하거나, 아래 숫자를 전체 뇌의 통합 지역 효과로 해석하면 안 된다.

| Gene | MOp가 높은 층화 | VISp가 높은 층화 | 두 지역 같은 값 | 보조 해석 |
|---|---:|---:|---:|---|
| Adora1 | 10 | 28 | 0 | 여러 뉴런 층화에서 VISp가 높으나 조건별 결과를 확인해야 함 |
| Adora2a | 5 | 10 | 23 | 23개 층화는 양쪽 모두 0 read; 최대 검출률 MOp 5.94%, VISp 8.82% |
| Nt5e | 14 | 24 | 0 | 원본 PBM 후보와 연결해 세포 유형별 검토 가능 |
| Grm2 | 18 | 19 | 1 | 방향이 조건에 따라 나뉨 |
| Penk | 21 | 17 | 0 | 방향이 조건에 따라 나뉨 |
| Cacna1c | 5 | 33 | 0 | 이 보조 자료의 표시 층화에서는 VISp 방향이 많음 |
| Mrps6 | 30 | 8 | 0 | 이 보조 자료의 표시 층화에서는 MOp 방향이 많음 |

`Adora2a`의 0 read는 발현이 없다는 증거가 아니며, 감소·증가의 확정적 근거로 삼을 수 없다. `panel22_direction_sensitivity_audit.csv`는 22개 모두의 방향·검출률 및 공유 cluster 표준화 후 방향을 제공한다.

공유 cluster만 남겨 두 지역에 동일한 pooled cell-count 가중치를 적용한 표준화 결과도 계산했다. 이렇게 해도 공유 cluster 안의 genotype·sex·age·FACS는 유지한다. 각 cluster의 작은 표집은 원자료에 보존되며, 표준화가 batch나 동물 설계를 보완하지는 않는다. `cortical_layer_label`은 원본 해부학적 dissection 표기이고 transcriptomic subclass와 별개다. 동일 dissection-layer 결과 역시 별도의 민감도 분석으로 제공한다.

## 주 10x 분석과의 관계

주 10x 분석과 **같은 genotype·male·P58·FACS·실제 subclass가 양쪽 지역에 모두 있는 SMART-seq 층화는 0개**다. 이는 각 지역 1 cell 이상만 요구해도 마찬가지다. 따라서 정확한 조건에서의 기술 간 재현을 확인했다고 결론 내릴 수 없다. `exact_metadata_10x_direction_crosscheck.csv`는 결과 행 없이 header를 남기며 `status.json`은 0개로 기록한다.

SMART-seq 결과는 다른 연령·reporter 층화에서의 별도 기술통계 자료로 사용한다. 10x UMI와 SMART-seq read counts를 합산하거나 평균·표준오차·P 값을 통합하지 않았다. 재현 코드는 동일 조건이 실제로 존재할 때만 별도 CP10k read/UMI 열의 방향을 비교하도록 되어 있다.

## 저장 위치와 재실행

- `input/rnaseq/smartseq/source_manifest.json`: 공식 다운로드 URL, 크기, SHA-256, TLS 검증 기록.
- `input/rnaseq/smartseq/GSE185862_CTX_HPF_2021_SS_readme.txt`: 대용량 raw cache 밖의 공식 readme 사본.
- `input/rnaseq/smartseq/GSM*_protocol_excerpt.soft.txt`와 `official_protocol_manifest.json`: 공식 GEO title·전체 protocol·processing 발췌 및 전체 응답/발췌 hash.
- `input/rnaseq/smartseq/mop_visp_candidate_exon_counts.h5ad`: 62개 후보 raw integer read counts, 원본 annotation, 전체 read-library 분모.
- `results/rnaseq/smartseq/SMARTseq_MOp_VISp_candidate_crosscheck.xlsx`: 표시 층화, 62개 결과, 동일 층 민감도, 10x 교차 비교 상태.
- `results/rnaseq/smartseq/*.csv.gz`: 용량을 줄인 전체 층화·cluster 원자료 요약. pandas `read_csv`로 바로 읽을 수 있음.
- `figures/rnaseq/smartseq/panel22_strict_metadata_difference_heatmap.png` 및 `.svg`: 38개 표시 층화의 22개 후보 차이.
- `results/validation/smartseq_actual_data_validation.json`: 실제 자료 기반 검증 기록.

프로젝트 루트에서 현재 환경의 Python을 사용해 실행한다. 전체 3.586 GB 원본 matrix는 `input/rnaseq/raw/`의 무시된 cache에 두고, 위 작은 재현용 입력과 출처·결과는 저장소에 보존한다. Windows에서도 프로젝트 폴더와 같은 Python 의존성 환경으로 실행할 수 있으며, root에 대한 절대 경로는 코드에 저장하지 않았다.

```bash
python src/rnaseq/smartseq_download.py GSE185862_expression_matrix_SSv4.hdf5
python src/rnaseq/smartseq_download.py GSE185862_metadata_ssv4.csv.gz
python src/rnaseq/smartseq_download.py GSE185862_CTX_HPF_2021_SS_readme.txt
# 공통 원본 XLSX는 주 RNA-seq 입력 준비 단계에서 확보됨
python src/rnaseq/smartseq_extract.py
python src/rnaseq/smartseq_analyze.py
python src/rnaseq/smartseq_validate.py
```

기존 cache가 있으면 downloader는 원본을 다시 받지 않고 공식 Content-Length와 비교한 뒤 hash를 기록한다. 실제 분석을 다시 실행하려면 주 10x 결과와 PBM 후보 CSV가 먼저 존재해야 한다. Git commit 또는 push는 수행하지 않았다.
