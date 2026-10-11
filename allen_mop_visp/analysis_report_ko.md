# MOp–VISp 실제 Allen 분석 결과

작성일: 2026-10-11. **공식 원본 다운로드에 성공했고, MOp–VISp 세포 구성 및 22개 핵심·62개 탐색 후보 발현의 실제 계산을 완료했다. NREM 관련 영역으로의 해부학적 투사도 실제 Connectivity 자료로 계산했다.** MERFISH의 공간 구성, RNA-seq의 세포 유형별 기저 발현, tracer 실험의 축삭 신호는 서로 다른 표본에서 얻은 기술 통계로 제시한다.

## 자료와 분석 범위

| 목적 | 실제 사용 자료 | 분석 표본 | 해석 범위 |
|---|---|---|---|
| 공간 세포 구성·층 분포 | ABC B6J MERFISH-C57BL6J-638850 | MOp 83,976 · VISp 61,884개; WT 수컷 한 마리 | 이 동물의 표집·annotation된 세포 구성; 동물 수준 검정·조직 density 없음 |
| 후보 발현 | GEO GSE185862 CTX-HPF 10x v2 | MOp 123,727 · VISp 29,763개; 지역별 donor 3개 | B6J 배경 Snap25-IRES2-Cre;Ai14 수컷 reporter, RFP-positive 선별; WT 공간 자료와 분리 |
| 별도 발현 보조 자료 | GEO GSE185862 SMART-seq v4 | MOp 6,588 · VISp 15,572개 | reporter·sex·age·FACS를 맞춘 자체 strata에서 비교; 10x와 합산하지 않음 |
| NREM 관련 영역으로의 투사 | Allen Mouse Brain Connectivity Atlas 공식 API | WT MOp 8 · VISp 33개 tracer 실험; Cre line 127개 별도 | 영역별 축삭 신호·연령·주입 순도 민감도; 기능 연결·NREM 유발 효과 아님 |
| 기존 PBM 연결 | 저장소의 원본 DEG·ORA | 공통 DEG 52개, 핵심 패널 22개; 합집합 62개 | 저산소 PFC·5xFAD 해마 처리 효과와 정상 atlas 기저 차이를 구분 |

ABC 공식 manifest에서 필요한 metadata 5개 **1,571,075,446 bytes**를 받아 크기·공개 MD5를 검증했다. RNA-seq은 기존 Allen S3 링크가 계속 차단되어 [같은 연구의 공식 GEO archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/)에서 10x 발현 행렬 **5,716,426,072 bytes** 및 metadata·원본 설명·donor 기록을 받았다. GEO 파일은 TLS·공식 응답 크기를 확인하고 로컬 SHA256을 기록했다. 공개 비교 checksum이 없는 GEO 파일의 로컬 SHA256을 publisher checksum으로 표현하지 않는다.

## 세포 구성

분모는 각 영역의 QC·annotation을 갖춘 전체 선택 세포이다. 유효 좌표와 annotation이 없는 선택 세포는 없었고, raw class/subclass·CCF 영역·해부학적 층·좌표를 보존했다.

| 그룹 | MOp (%) | VISp (%) | MOp − VISp (percentage points) |
|---|---:|---:|---:|
| 흥분성 뉴런 | 58.28 | 69.48 | −11.20 |
| 억제성 뉴런 | 9.21 | 8.07 | +1.14 |
| Astrocyte | 10.86 | 8.58 | +2.29 |
| Oligodendrocyte | 5.79 | 3.06 | +2.73 |
| OPC | 1.52 | 1.24 | +0.27 |
| Microglia | 1.86 | 1.21 | +0.65 |
| Endothelial | 8.81 | 6.06 | +2.74 |
| Pericyte/SMC | 2.11 | 1.09 | +1.02 |
| VLMC | 1.09 | 0.71 | +0.37 |
| 기타 비뉴런 | 0.48 | 0.50 | −0.02 |

GABA 세포 내부에서 Pvalb 비율은 MOp 35.83%, VISp 42.74%이며 Lamp5는 16.19%, 11.60%다. 전체 세포 분모와 GABA 내부 분모는 구분한다. 뉴런 내부 L4/5 IT 비율은 MOp 18.05%, VISp 34.49%다.

**전체 영역에서 관찰한 차이를 모든 층에 적용할 수 없다.** 전체 세포 중 GABA 비율은 MOp가 높지만, L5는 MOp 10.77%·VISp 12.56%, L6b는 3.55%·4.20%로 방향이 다르다. CCF가 MOp에 별도 L4 영역 label을 제공하지 않는 것은 해당 층의 세포가 생물학적으로 없다는 뜻이 아니다. 해부학적 층과 transcriptomic `L4/5 IT` annotation은 별도 정보다.

[전체 구성 보고서](results/merfish/composition_report_ko.md), [구성 그림](figures/merfish/regional_major_composition.png), [층별 그림](figures/merfish/layer_major_composition.png)에 세부 결과가 있다. 32개 구성 CSV의 모든 층화 분모·합계·지역 차이를 독립 검증했다. 단일 동물의 절편·반구·세포를 동물 반복으로 취급하지 않았다.

## 같은 세포 유형의 후보 발현

원본 10x 행렬은 31,053 genes × 1,169,320 cells이다. metadata의 exact MOp/VISp 153,490개 cell ID를 행렬과 전부 정확히 대응시켰고, 후보 62개를 모두 확보했다. 후보 raw counts와 **원본 31,053개 유전자의 세포별 전체 count 합계**를 함께 보존했다. 정규화는 `log1p(10,000 × gene count / full-library count)`이고 검출률은 raw count > 0인 세포 비율이다. 후보 62개의 합계를 분모로 사용하지 않는다. 평균 차이는 log1p(CP10k)의 차이이며 gene log2FC가 아니다.

sex·genotype·chemistry·FACS·정확한 연령·subclass를 맞춘 비교는 공통 연령 **P58**에서 가능하다. 다만 P58도 MOp donor 1개와 VISp donor 2개이며 library-preparation date가 완전히 분리되어 있다. 같은 ROI pooling의 독립 단위도 확정하지 못했으므로 **지역 효과와 batch 효과를 분리하는 동물 수준 p-value/FDR를 생성하지 않았다.** P57–59 전체 표본의 탐색 결과는 연령까지 맞춘 P58 결과와 별도로 표시한다.

P58의 실제 예시는 다음과 같다. 모두 같은 subclass 안의 평균이며, 검출률은 괄호에 표시했다.

| Gene / subclass | MOp 평균 (검출률 %) | VISp 평균 (검출률 %) | MOp / VISp 세포 수 |
|---|---:|---:|---:|
| Adora1 / Pvalb | 0.3545 (45.81) | 0.2067 (35.31) | 1,576 / 422 |
| Adora1 / L2/3 IT CTX | 0.4410 (55.51) | 0.5779 (72.29) | 4,014 / 2,104 |
| Adora1 / L4/5 IT CTX | 0.8984 (78.91) | 1.0006 (86.78) | 13,373 / 6,567 |
| Adora1 / L6 CT CTX | 1.0793 (87.89) | 1.1805 (93.99) | 7,790 / 1,948 |
| Nt5e / Pvalb | 0.0655 (10.47) | 0.1112 (19.67) | 1,576 / 422 |
| Penk / Vip | 0.6840 (46.67) | 0.8661 (59.41) | 2,025 / 882 |

Adora1은 Pvalb에서 MOp가 높고, 위 IT·CT 그룹에서는 VISp가 높았다. 양쪽에 공유되는 cluster를 같은 가중치로 비교해도 이 예시의 방향은 유지됐다. 따라서 adenosine 수용체 발현의 지역 차이를 뉴런 전체에 하나의 방향으로 요약하는 것은 적절하지 않다.

반면 **subclass 평균과 공유 cluster 표준화 결과가 다를 수 있다.** 예를 들어 Nt5e/Lamp5의 MOp−VISp 차이는 일반 평균 −0.0139에서 공유 cluster 동일 가중치 +0.0298로, Penk/L2/3 IT는 +0.0230에서 −0.0358로 바뀐다. 비교에 포함되는 cluster와 구성비가 달라지므로 이 결과만으로 순수한 지역 발현 차이를 확정하지 않는다. 동일 가중치 표에는 작은 cluster도 포함되며 희소 표본의 과도한 해석을 피한다.

Adora2a는 주요 IT·Pvalb·Sst·Vip 등에서 10x 검출이 매우 낮다. 관측 0을 생물학적 미발현으로 판정하거나 낮은 검출의 작은 차이를 우선 후보 선정 근거로 삼지 않는다. 비뉴런도 RFP-positive 선별에서 적게 표집되므로 glia·혈관 발현의 대표성에 제한이 있다. SMART-seq 보조 자료의 조건별 확인은 별도 결과로 관리한다.

실제 발현 표는 `results/rnaseq/region_subclass_gene_summary.csv`, `within_subclass_region_difference.csv`, `shared_cluster_standardized_difference.csv`에 있다. [RNA-seq 상세 보고서](results/rnaseq/RNAseq_results_ko.md)와 [실행·해석 안내](docs/RNAseq_analysis_ko.md)를 함께 참고한다.

## SMART-seq 보조 자료

공식 GEO SMART-seq v4 exon read-count 원본 **3,585,764,160 bytes**도 확보해 exact MOp 6,588개·VISp 15,572개 세포와 후보 62개를 별도로 분석했다. 정규화 분모는 원본 전체 **45,768개 feature의 exon read 합계**이며 결과 단위는 `log1p(exon read CPM)`이다. UMI 자료와 같은 값으로 합치지 않는다.

연령이 없거나 기록된 injection이 있는 세포를 직접 비교에서 제외하고 genotype·sex·정확한 age·FACS·subclass를 맞췄다. 각 영역 n≥10으로 표시한 비교 strata는 38개이며 해당 표집은 MOp 2,110개·VISp 1,871개다. 같은 해부학적 dissection-layer label도 맞춘 민감도 비교는 22개 strata다. n≥10은 희소 화면을 줄이는 표시 기준이며 통계적 유의성 기준이 아니다.

**10x와 모든 metadata 조건이 같은 지역 대조는 0개다.** 따라서 SMART-seq은 다른 reporter·연령·선별 조건에서 후보 발현을 확인하는 보조 자료이며, 위 10x의 지역 차이 방향이 재현됐다는 근거로 사용하지 않는다. WT 반복 자료도 아니다. [SMART-seq 보고서·재현 방법](docs/SMARTseq_analysis_ko.md)과 `results/rnaseq/smartseq/`의 CSV·Excel에 원본 조건과 결과를 보존했다.

## NREM 관련 영역으로의 실제 투사

공식 ontology와 기능 논문을 대조해 VLPO, MnPO(Allen `MEPO`), ZI, SI, MA, NDB, 시상 망상핵(`RT`)을 각각 조회했다. PZ는 exact atlas label이 없어 **미평가**로 남겼다. `TRN`은 Allen에서는 다른 뇌간 영역이므로 시상 망상핵으로 사용하지 않았다. [표적 정의와 문헌](docs/NREM_projection_target_definition_ko.md)에 NREM·REM·각성의 근거와 세포 유형을 구분했다.

MOp·VISp primary injection 168개 중 WT 41개를 주 분석으로 사용하고, Cre 127개는 line별로 분리했다. 일곱 mapped target의 좌·우·양측 3,528개 기록을 모두 확보했다. 아래는 WT 실험의 **양측 normalized projection volume(NPV) 중앙값**이다. NPV는 표적의 segmented projection signal volume을 전체 injection signal volume으로 나눈 값이며, 투사 확률·뉴런 비율·시냅스 수가 아니다.

| 표적 | MOp (8개 실험) | VISp (33개 실험) | 해석 |
|---|---:|---:|---|
| ZI | 0.17054 | 0.02560 | 후속 해부학 후보; Lhx6 세포에 도달하는지는 미확인 |
| VLPO | 4.16×10⁻⁸ | 8.23×10⁻⁶ | 신호가 작고 주입 조건 제한에 따라 비교 방향이 바뀜 |
| MEPO (MnPO) | 2.16×10⁻⁷ | 6.70×10⁻⁹ | 중앙값이 매우 작음; Vgat 세포 targeting 미확인 |
| RT | 0.10245 | 0.09633 | 양쪽 source에서 신호; NREM spindle·delta 관련 영역 |
| SI | 0.00302 | 0.00010 | mixed-function BF 영역; SOM 세포 targeting 미확인 |

MOp WT 실험은 모두 P56이고 VISp는 P49–80이다. 전체 injection signal 중 source parent 영역에 포함되는 비율이 0.8 이상인 실험은 MOp 4개·VISp 24개다. **P56이면서 이 비율이 0.8 이상인 실험**으로 제한하면 MOp 4개·VISp 10개이며, ZI 중앙값은 0.26841·0.04803으로 MOp가 큰 경향이 남는다. 0.8은 탐색적 주입 확산 민감도 기준이며 공식 QC 또는 기능 threshold가 아니다. Source만 annotation된 실험은 별도 기준으로 MOp 2개·VISp 16개다.

연령 제한과 signal-volume 정규화만으로 tracer 주입 부피·주변 영역 확산·프로토콜·실험 배치의 차이가 제거되는 것은 아니다. P value/FDR를 계산하지 않았으며, **MOp가 더 강한 NREM을 유발한다는 결론은 없다.** Segmented signal의 0은 이 실험 조건의 관측값이며 축삭이나 기능의 완전한 부재가 아니다. PZ 미평가는 0으로 바꾸지 않았다.

우선 회로 검증 후보는 **MOp→ZI**다. Tracing에서 ventral ZI·Lhx6 세포 targeting을 확인하고, 시냅스/생리 및 projection별 개입과 EEG/EMG로 NREM·REM·wake를 측정해야 한다. 사용자 자극 좌표와 atlas injection의 위치 대응은 아직 등록하지 않았다. 새 분자 pathway 점수는 실행하지 않았으며, 전장 유전자/경로 구성 유전자 추출과 아데노신 하위 신호 분석 방법은 [후속 분석 제안](docs/downstream_pathway_analysis_options_ko.md)에 있다. [투사 보고서](results/connectivity/NREM_projection_report_ko.md)와 [별도 Excel](results/connectivity/MOp_VISp_NREM_projection.xlsx)에 개별 실험·민감도·해석을 정리했으며 그림은 `figures/connectivity/`에 있다.

## PBM 연결과 미완료 범위

기존 PBM 공통 DEG는 52개(up 34/down 18), 반대 방향 23개다. 동시에 유의한 GO는 38개(BP 29/MF 9), 실제 공통 hit 합집합은 19개다. FDR05 표의 부재를 미발현·0효과로 바꾸지 않았다. [PBM 연결 보고서](results/pbm/PBM_bridge_report_ko.md)에 원본·검증·후보 연결이 있다.

Nt5e의 피질 PBM log2FC는 −0.5367(FDR 0.00744556), Adora2a는 −1.2929(FDR 0.00932146)였다. 이는 **저산소 PFC의 처리 효과**이고 정상 MOp–VISp의 기저 차이·정상 B6 PBM 반응 방향과 같은 값이 아니다. 현재 atlas 결과로 adenosine 농도·효소 flux·단백질 기능·수면 threshold를 산출하지 않는다.

B6J MERFISH 500-gene 패널에서는 핵심 22개 중 Penk만 측정하며 탐색 62개 중 7개를 측정한다. 미측정 유전자는 0 발현이 아니다. 이번 MERFISH는 metadata 구성 분석을 완료했고 별도 대형 raw matrix의 Penk 공간 발현은 아직 계산하지 않았다. B6N Zhuang-1 패널의 12개 탐색 후보를 B6J 동물 반복으로 합치지 않았다.

국소 자극 ROI는 여전히 보류했다. 사용자가 확인한 mm·bregma AP/ML, 뇌 표면 DV, core 직경 400 µm는 [기록 파일](docs/stimulation_coordinates_user.json)에 있다. M1은 `(1.34, 2.46, 0)`, V1은 `(약 −2.3, 2.46, 0)`이며 V1 AP가 잠정값이다. bregma→CCF 등록과 광학적 범위를 확정하기 전에는 해당 좌표를 CCF에 직접 대입하지 않는다. 최종 qPCR 3–5개 패널은 이번 기술 통계·측정 가능성과 실제 실험 조건을 함께 평가해 좁혀야 한다.

## 검증·보관·재현

실제 MERFISH 145,860개 선택 cell과 RNA-seq 153,490개 선택 cell의 annotation·행렬 대응을 검증했다. RNA-seq 원본에서 donor별 표본 cell을 다시 읽어 전체 유전자 정규화 분모와 후보 raw counts를 대조했고, 저장 subset과 정확히 일치했다. subclass·cluster의 평균·검출률·세포 수·공유 cluster 가중치 결과도 독립 재계산했다. 상세 기록은 `results/merfish/actual_data_validation.json`과 `results/validation/`에 있다.

Connectivity는 공식 API 원본 응답·URL·hash를 보존하고, 선택 실험 조건·NPV 정규화 분모·전체 표적 기록·다섯 조건별 WT 요약·Cre line별 요약을 독립 재계산했다. 원본/코드는 `input/connectivity/`, `src/connectivity/`에 있으며 검증 기록은 `results/validation/connectivity_actual_data_review.json`에 있다.

정본은 `/workspace/tPBM_opendataset/allen_mop_visp`이다. Windows `G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`은 클라우드에 연결되어 있지 않다. 대형 원본 캐시는 Git·이동용 ZIP에서 제외하고, 선택 metadata·후보 counts·전체 library totals·코드·URL·hash·결과를 보존한다. [재현 방법](README_ko.md), [파일 관리 안내](docs/file_management_ko.md), 통합 `results/MOp_VISp_Allen_comparison.xlsx`, 최신 `exports/Allen_Dataset_M1_V1.zip`을 이용한다.
