# Allen B6J MERFISH: MOp–VISp 실제 세포 구성 비교

2026-10-11

**공식 Allen 원본을 확보하여 MOp 83,976개, VISp 61,884개, 총 145,860개의 QC annotated 세포를 비교했다.** 이 C57BL/6J WT 수컷 한 동물에서 VISp는 전체 표집 세포 중 glutamatergic neuron 비율이 더 높았고, MOp는 astrocyte·oligodendrocyte 및 vascular 범주 비율이 더 높았다. 뉴런 내부에서도 원본 L4/5 IT 및 L6 IT subclass의 구성 차이가 확인되었다. 아래 값은 **한 동물의 공개 자료에서 표집·등록된 세포의 기술 통계**이며 동물 간 재현성 또는 PBM 반응 차이를 입증하지 않는다.

## 원본과 분석 범위

- 데이터: Allen ABC `MERFISH-C57BL6J-638850`, CCF mapping `MERFISH-C57BL6J-638850-CCF`.
- 실제 donor metadata: `C57BL6J-638850`, `wt/wt`, 수컷 `M`. 독립 동물 **1개**.
- 입력: live release manifest `20250531`; 세포 annotation `20241115`, 좌표 `20231215`, CCF parcellation `20230630`.
- 공식 HTTPS에서 필수 입력 5개, 총 **1,571,075,446 bytes**를 확보했고 manifest의 원본 MD5와 파일 크기를 모두 검증했다. TLS 검증은 유지했다.
- 공개 annotation 전체 3,938,808개, reconstructed-coordinate 전체 3,739,961개. 두 전체 파일의 수가 다르므로 아래 분석 분모는 **선택된 MOp/VISp의 QC annotated 세포**로 명시한다.
- CCF `structure`가 정확히 MOp 또는 VISp인 145,860개를 선택했다. annotation join 결측 0개, CCF 좌표 결측 0개, QC annotation 제외 0개, 중복 cell label 0개.
- MOp 17개 절편, VISp 14개 절편을 포함했다. 절편·반구·세포를 독립 동물로 세지 않았다.

원본 Allen class, subclass, supertype, cluster, neurotransmitter와 correlation score를 보존했다. 표시용 그룹은 별도 매핑 열이며 단일 marker 발현으로 세포를 재분류하지 않았다. 새 correlation-score cutoff를 만들지 않았다.

## 영역 전체 구성

분모는 각 영역의 **전체 QC annotated 세포**다.

| 표시 그룹 | MOp n | MOp % | VISp n | VISp % |
|---|---:|---:|---:|---:|
| Glutamatergic neuron | 48,941 | 58.28 | 42,995 | 69.48 |
| GABAergic neuron | 7,731 | 9.21 | 4,993 | 8.07 |
| Astrocyte | 9,123 | 10.86 | 5,307 | 8.58 |
| Oligodendrocyte | 4,865 | 5.79 | 1,896 | 3.06 |
| OPC | 1,273 | 1.52 | 768 | 1.24 |
| Microglia | 1,560 | 1.86 | 748 | 1.21 |
| Endothelial | 7,395 | 8.81 | 3,751 | 6.06 |
| Pericyte/SMC | 1,776 | 2.11 | 677 | 1.09 |
| VLMC | 912 | 1.09 | 442 | 0.71 |
| Other non-neuronal | 400 | 0.48 | 307 | 0.50 |
| 분모 합계 | **83,976** | **100** | **61,884** | **100** |

뉴런 분모는 MOp 56,672개, VISp 47,988개다. **뉴런 중 GABAergic neuron** 비율은 MOp **13.64%**, VISp **10.40%**다. 전체 세포 대비 GABA 비율 9.21%/8.07%와 분모가 다르다.

Other non-neuronal 범주는 원본 label을 확인할 수 있도록 유지했으며 다른 glial/vascular 범주로 강제 합치지 않았다. 대부분은 원본 `329 ABC NN`(MOp 368개, VISp 299개)이며 나머지는 BAM·DC·lymphoid 등 원본 immune 범주와 소수 Ependymal/CHOR다. 전체 세포 대비 비율은 density나 조직 전체 절대 세포 수가 아니다. 현재 공간 면적 마스크가 없어 cells/mm² 또는 cells/mm³를 계산하지 않았다.

## GABAergic neuron 내부 구성

분모는 MOp **7,731개**, VISp **4,993개**의 GABAergic neuron이다. 원본 Pvalb chandelier는 표시상 Pvalb에, Sst Chodl은 Sst에 포함했으며 원본 subclass는 별도로 보존했다.

| 표시 family | MOp n | MOp % of GABA | VISp n | VISp % of GABA |
|---|---:|---:|---:|---:|
| Pvalb | 2,770 | 35.83 | 2,134 | 42.74 |
| Sst | 2,314 | 29.93 | 1,485 | 29.74 |
| Vip | 1,046 | 13.53 | 713 | 14.28 |
| Lamp5 | 1,252 | 16.19 | 579 | 11.60 |
| Sncg | 317 | 4.10 | 76 | 1.52 |
| Other GABA | 32 | 0.41 | 6 | 0.12 |

이 동물에서 GABA 내부 Pvalb 비율은 VISp가 약 **6.91 percentage points** 높았다. Sst의 GABA 내부 비율은 29.93%와 29.74%로 비슷한 기술적 값이다. Sst를 전체 QC annotated 세포로 나누면 MOp **2.76%**, VISp **2.40%**가 된다. 이는 분모를 바꾸면 같은 세포군의 지역 차이가 다르게 보일 수 있음을 보여준다.

## 원본 뉴런 subclass 구성

분모는 각 영역의 전체 뉴런이다. 큰 차이를 보이는 원본 subclass 일부를 표시했다. 이 표의 transcriptomic layer 명칭은 CCF 해부학적 layer가 아니다.

| 원본 subclass | MOp % of neurons | VISp % of neurons | MOp−VISp percentage points |
|---|---:|---:|---:|
| `006 L4/5 IT CTX Glut` | 18.05 | 34.49 | −16.44 |
| `004 L6 IT CTX Glut` | 12.68 | 6.53 | +6.15 |
| `005 L5 IT CTX Glut` | 6.46 | 4.00 | +2.46 |
| `029 L6b CTX Glut` | 2.82 | 0.88 | +1.94 |
| `030 L6 CT CTX Glut` | 17.32 | 15.64 | +1.67 |
| `007 L2/3 IT CTX Glut` | 19.66 | 20.56 | −0.90 |

`L4/5 IT`의 구성 차이는 후보 유전자의 전체 영역 평균이 세포 유형 비율에 영향을 받을 가능성을 제시한다. 실제 gene expression 비교에서는 원본 subclass 내부 및 가능한 공유 하위 유형에서 지역 차이를 따로 확인해야 한다. 큰 subclass에서 지역 평균 발현이 다르더라도 within-type 조절 변화로 바로 해석하지 않는다.

## 해부학적 층·반구·절편

CCF `substructure`에서 가져온 해부학적 layer별 분모다. MOp에는 이 parcellation에서 해부학적 L4 label이 없으며 이를 0 발현 또는 transcriptomic L4/5 IT 부재로 해석하지 않는다. 두 영역 모두 anatomical-layer unknown 세포는 0개였다.

| CCF anatomical layer | MOp QC annotated n | VISp QC annotated n |
|---|---:|---:|
| L1 | 6,739 | 4,311 |
| L2/3 | 33,777 | 21,574 |
| L4 | 이 영역의 해당 CCF label 없음 | 12,977 |
| L5 | 19,352 | 10,991 |
| L6a | 21,882 | 10,128 |
| L6b | 2,226 | 1,903 |

공유하는 해부학적 층에서 전체 QC annotated 세포 대비 비율을 비교하면 다음과 같다.

| CCF layer | Glut MOp % | Glut VISp % | GABA MOp % | GABA VISp % |
|---|---:|---:|---:|---:|
| L1 | 8.07 | 13.38 | 10.86 | 9.12 |
| L2/3 | 62.78 | 76.02 | 9.66 | 6.50 |
| L5 | 57.52 | 63.34 | 10.77 | 12.56 |
| L6a | 67.27 | 80.78 | 7.19 | 5.41 |
| L6b | 60.15 | 75.25 | 3.55 | 4.20 |

이 동물에서 glutamatergic neuron 비율은 모든 공유 CCF layer에서 VISp가 더 높았다. GABA 비율은 층에 따라 방향이 달랐다. 특히 L5와 L6b에서는 영역 전체와 달리 VISp가 더 높으므로 **전체 영역 GABA 차이를 모든 자극 깊이에 동일하게 적용하지 않는다.** L1의 astrocyte 비율은 MOp 37.05%, VISp 37.14%로 유사한 기술적 값이었다. 해당 패턴은 `composition_layer_major.csv`에서 원본 세포 수와 분모까지 확인할 수 있다.

CCF z=5.7 mm midline 기준의 left/right 표집 수는 MOp **39,938/44,038**, VISp **32,293/29,591**이다. 반구 unknown 및 midline 세포는 0개였다. 원본 mm CCF 좌표와 section label을 저장했으며, 각 영역에서 layer, hemisphere, section 및 그 조합별 구성 CSV를 산출했다. 절편별 비교는 동일 동물 내 공간 표집 차이를 확인하는 자료다.

CCF layer 구성의 차이와 자극 주변 국소 ROI가 아직 같다고 가정하지 않는다. 사용자 AP/ML은 bregma, DV는 뇌 표면 기준이므로 검증된 bregma–CCF 변환과 광분포·깊이 정보가 있어야 국소 ROI를 계산할 수 있다.

## 후보 유전자 실측 여부와 다음 분석

공식 B6J gene metadata를 실제로 읽어 22개 atlas screening 후보의 coverage를 재확인했다. **Penk만 직접 측정되고 나머지 21개는 이 500-gene 패널에서 미측정**이다. 이 21개를 0 발현으로 표시하지 않는다.

현재 필수 세포 metadata 분석을 우선 완료했고 선택 사항인 7.63 GB raw H5AD는 다운로드하지 않았다. **Penk의 실제 발현 평균·검출률은 아직 계산하지 않았다.** 폭넓은 22개 및 PBM 공통 후보의 발현 비교는 별도의 CTX-HPF RNA-seq 분석에서 진행한다. 세포 구성만으로 adenosine 생성·순축적·수면 threshold 또는 PBM 반응 방향을 판정하지 않는다.

## 검증과 파일

실제 세포의 고유 ID, exact region, 단일 WT 수컷 donor, CCF 좌표 완전성을 확인했다. **32개 구성 CSV**의 모든 층화 단위에서 count 합계가 분모와 일치하고 fraction 합계가 1임을 검증했다. 공식 원본 checksum 5개와 coverage Penk-only 조건도 확인했다. 단일 동물이므로 animal-level 검정, p value 또는 독립 동물 confidence interval은 계산하지 않았다.

- `readiness_report.json`: 실제 실행 `exit 0`, `composition_completed`, `actual_analysis_computed=true`와 원본 검증.
- `actual_data_validation.json`: 실제 결과 분모·합계·좌표 검증.
- `selected_cells_raw_labels_and_coordinates.csv.gz`: 모든 선택 세포의 원본 annotation·좌표·표시 그룹.
- `raw_label_display_group_mapping.csv`: 버전별 원본 label → 표시 그룹.
- `composition_region_major.csv`, `composition_region_gaba_within_gaba.csv`: 위 주요 표의 정확한 수치.
- `regional_descriptive_difference_neuron_subclass.csv`: 원본 뉴런 subclass 전체 차이.
- `sample_coverage.csv` 및 `composition_layer*`, `composition_hemisphere*`, `composition_section*`: 공간 표집 분모와 구성.
- `../../figures/merfish/regional_major_composition.png`, `layer_major_composition.png`, `regional_neuron_subclasses.png`: 실제 결과 그림 3개.

재현 실행은 프로젝트 최상위에서 `python src/merfish/run_merfish.py`다. 원본이 없을 때는 `--download`로 먼저 확보한다. 원본 cache는 commit 및 일반 산출물 패키지에서 제외하고 공식 URL·크기·MD5와 재현 코드를 함께 보존한다.
