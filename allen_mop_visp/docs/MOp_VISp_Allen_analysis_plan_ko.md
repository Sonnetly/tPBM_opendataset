# MOp와 VISp의 세포 구성 및 유전자 발현 비교 계획

2026년 10월 11일

목표는 현재 조건에서 관찰된 M1과 V1의 tPBM 수면 반응 차이에 연결할 세포 유형과 분자 후보를 선정하는 것이다. 분석은 **세포 구성, 같은 세포 유형 안의 기저 발현, PBM 처리 후 변화**를 구분한다. 기존 tPBM 공개 전사체에서 관찰된 변화와 Allen의 정상 상태 atlas를 연결하되, 정상 atlas를 PBM 반응 데이터로 해석하지 않는다.

이번 작업에서는 레포의 최신 DEG 검증 자료를 읽고, 공통 유전자 및 GO 수치를 점검하며, Allen 공간 데이터의 후보 유전자 실측 여부를 확인하였다. **Allen MOp–VISp의 실제 발현 차이와 세포 비율은 아직 계산하지 않았다.** 아래 순위는 다음 분석에 사용할 후보 순위이며 지역 차이가 확인된 결과가 아니다.

## 비교할 연구 질문

1. MOp와 VISp의 세포 유형 비율과 층별 분포가 다른가.
2. 서로 대응되는 세포 유형 안에서 adenosine 생성·수송·대사·수용체 관련 기저 발현이 다른가.
3. 기존 PBM 데이터에서 선별한 Ca²⁺ 처리·시냅스·막전위 관련 후보의 기저 발현이 다른가.
4. 이 차이가 어떤 국소 측정과 기능 검증을 우선하게 하는가.

랩미팅의 ‘adenosine production 차이’는 먼저 **세포외 adenosine의 순축적 차이**로 정의한다. 측정 농도의 차이는 생성뿐 아니라 수송·대사·제거에도 영향을 받기 때문이다. ‘수면 유도 threshold 차이’는 유전자 발현만으로 판정할 수 없다. 수용체와 세포 유형별 발현은 동일한 adenosine 입력에 대한 기능 비교의 후보를 제공한다.

## 사용할 데이터와 역할

| 데이터 | 역할 | 적용 조건 및 한계 |
|---|---|---|
| Allen ABC MERFISH-C57BL6J-638850 및 CCF 좌표 | C57BL/6J에서 MOp와 VISp의 세포 구성·공간 분포 비교 | 공개된 단일 뇌 데이터이다. 절편과 반구를 독립 동물로 세지 않는다. 실측 패널은 500 genes이다. |
| Allen Mouse Whole Cortex and Hippocampus 10x 및 2021 Cell taxonomy | MOp와 VISp가 별도로 채취된 자료에서 세포 유형별 후보 발현 비교 | 약 8주령, 암수, B6J 배경의 pan-neuronal transgenic 자료가 중심이다. WT B6J와 구분하고 donor/pooling 및 선별 조건을 먼저 확인한다. |
| Whole Cortex and Hippocampus SMART-seq | 낮은 검출률을 보이는 후보의 보완 확인 | 같은 자료 내 MOp–VISp 비교를 따로 수행한다. 10x와 발현값을 직접 합치지 않는다. |
| Allen WMB-10X | 전체 세포 유형의 후보 발현 분포 참고 | 공개 메타데이터에서 MOp는 있으나 visual ROI는 VIS/VIS-PTLp 등이다. 이를 VISp로 바꾸어 직접 비교하지 않는다. |
| Zhuang-ABCA-1–4, Zhang et al. Nature 2023 | 별도 공간 자료에서 구성·층 분포 및 실측 후보의 재현 방향 확인 | C57BL/6NCrl이다. B6J 주 분석과 합산하지 않는다. B6J에 대한 독립 동물 반복을 대신하지 않는다. |
| Allen Mouse Brain Atlas ISH | 상위 후보의 영역·층 분포 보완 | 같은 유전자의 MOp–VISp 패턴을 확인한다. probe가 다른 유전자끼리 영상 강도를 비교하지 않으며 세포 유형별 단백질 정량을 대신하지 않는다. |

‘B6 기준’은 우선 C57BL/6J를 주 분석으로 설정한다. B6J 배경의 reporter/Cre line과 WT는 별도 표시한다. WT-only 조건에서 세포 유형별 RNA-seq 비교가 성립하지 않으면 그 제약을 명시하고, B6J 배경 자료를 보조 분석으로 표시한다.

현재 확인된 B6J 500-gene 패널에는 adenosine 후보 8개가 모두 없고, 아래 22개 중 Penk만 직접 측정된다. Zhuang-ABCA-1의 1122-gene 패널에는 Adora1, Adora2a, Trpc5, Penk가 있다. 따라서 세포 구성은 MERFISH, 후보의 폭넓은 발현 비교는 RNA-seq를 사용한다. (2026-10-16 진행) Zhuang-ABCA-1~4 네 동물의 subclass 안 MOp 대 VISp 비교와 10x 전사체 전체의 뉴런 유형별 분석을 수행했다. [Zhuang 보고서](../results/zhuang/zhuang_report_ko.md), [뉴런 유형별 발현 보고서](../results/neuron_type_report/M1_V1_neuron_type_expression_report_ko.md). **패널에 없는 유전자는 미측정이며 0 발현이 아니다.** Imputed expression은 예측값으로 별도 표시하고 실제 지역 차이 검증의 주 결과로 사용하지 않는다.

## 기억하는 Nature 논문과 연결

정확한 논문은 현재 설명만으로 하나로 특정하지 않는다. 다음 논문들이 관련성이 높다.

- **Tasic et al., Nature 2018, Shared and distinct transcriptomic cell types across neocortical areas**: VISp와 ALM의 비교이다. ALM을 MOp로 대체하지 않는다. shared/region-specific cell type의 해석에 참고한다.
- **Zhang et al., Nature 2021, Spatially resolved cell atlas of the mouse primary motor cortex by MERFISH**: MOp의 뉴런·비뉴런 구성 및 cortical depth 분포를 직접 보여준다. Figure 1과 Extended Data Figure 2f, 3가 현재 구성·층 분석 설계와 연결된다. VISp가 없는 MOp 중심 자료이므로 직접 지역 비교의 주 데이터로 사용하지 않는다.
- **Yao et al., Nature 2023, A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain**: 현재 Allen ABC의 전체 뇌 taxonomy 및 B6J 공간 분석을 뒷받침한다.
- **Yao et al., Cell 2021, A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation**: MOp와 VISp를 모두 포함하는 발현 비교의 우선 후보 자료이다. 저널은 Nature가 아니라 Cell이다.

Allen과 Nature 논문은 일부가 같은 데이터의 포털과 출판물 관계이다. 출처 이름이 다르다는 이유로 독립적 재현 자료로 세지 않는다.

## 세포 유형 구분

기존 Allen annotation을 유지하고, 분석용으로 대응 가능한 큰 그룹을 먼저 만든다. 유전자 하나의 발현 여부로 모든 세포를 새로 분류하지 않는다.

| 단계 | 구분 | 주 목적 |
|---|---|---|
| 전체 조직 | Glutamatergic neuron, GABAergic neuron, Astrocyte, Oligodendrocyte, OPC, Microglia, Endothelial, Pericyte/SMC, VLMC 및 기타 | 신경세포와 glia/vascular 세포의 구성 비교 |
| 흥분성 뉴런 | L2/3 IT, L4/5 IT, L5 IT, L6 IT/Car3, L5 ET, L5/6 NP, L6 CT, L6b 등 실제 annotation에서 확보되는 그룹 | 층·투사 유형과 후보 발현의 연결 |
| 억제성 뉴런 | Pvalb, Sst, Vip, Lamp5, Sncg를 출발점으로 최신 taxonomy와 대응 | 억제성 회로 구성과 후보 발현 비교 |
| 세부 유형 | 두 영역에서 충분히 확보되는 supertype/cluster | 큰 그룹에서 보인 차이가 하위 유형 비율 때문에 생겼는지 확인 |

Allen 2021과 2023 taxonomy의 명칭을 그대로 동일 유형으로 간주하지 않는다. 분석 시작 시 raw label → 표시용 그룹 매핑을 버전별로 저장한다. 두 영역에 한쪽만 있는 유형은 ‘유형 구성의 차이’로 보고, 같은 유형의 지역별 발현 차이로 검정하지 않는다.

전체 QC·분류 세포 대비 비율과 뉴런 안에서의 비율을 구분한다. 예를 들어 Sst 비율은 전체 세포 대비와 GABAergic neuron 대비를 각각 표시한다. 공간 면적 마스크가 확보될 때만 cells/mm²를 추가하며, 3D 표본 체적 보정 없이 cells/mm³나 전체 뇌 세포 수를 보고하지 않는다. 해리·FACS로 얻은 RNA-seq cell count를 조직의 실제 세포 비율로 사용하지 않는다.

## 층과 실제 자극 위치

MOp 전체와 VISp 전체를 먼저 비교한 뒤 실제 자극 AP/ML/DV 및 깊이에 대응되는 국소 ROI로 좁힌다. 영역 전체의 차이가 자극 부위 주변에서도 유지되는지가 중요하다.

해부학적 layer label과 transcriptomic subclass를 별도 열로 유지한다. MOp의 L4/5 IT-like transcriptomic population을 VISp의 해부학적 L4와 동일시하지 않는다. 두 영역의 층 구조를 맞추기 위해 superficial/deep grouping 또는 정규화된 cortical depth를 보조적으로 사용하되 CCF layer 분석도 보존한다. 국소 ROI의 크기는 fiber core 직경만으로 정하지 않고 실제 광분포·자극 깊이 정보로 설정한다. 영역 경계 포함과 경계 제외 분석에서 결론이 유지되는지도 확인한다.

## 기존 tPBM 분석에서 확인한 내용

레포 commit e430e8d2d2827ccdd222e8828d77227e937b6f1b의 최신 DEG 검증 자료를 기준으로 한다.

- 같은 방향 overlap은 **52개: up 34개, down 18개**이다.
- 두 up ORA에서 모두 유의한 GO는 **38개**이며, 그 항목들의 실제 공통 hit의 합집합은 **19개**이다.
- 기존 1차 qPCR 후보 8개는 Grm2, Ryr2, Trpc5, Atp2b3, Cacna1c, Grin2a, Mrps6, Penk이다.
- 피질 PBM DEG 표에서 **Nt5e는 down, log2FC −0.5367, gene FDR 0.00745**, **Adora2a는 down, log2FC −1.2929, gene FDR 0.00932**이다.
- 나머지 adenosine 후보 6개가 제공 FDR05 표에 없다는 것은 미발현·변화 없음·PBM 비관련을 뜻하지 않는다. 전체 tested-gene 통계가 있어야 유의성 외의 정보를 확인할 수 있다.

해마 비교는 5xFAD + light vs 5xFAD, 피질 비교는 hypobaric hypoxia + PBM vs hypobaric hypoxia이다. 이를 정상 B6의 MOp–VISp 비교로 해석하지 않는다. 공개 PBM 반응 유전자의 기저 발현이 높은 지역이 반드시 PBM에 더 크게 반응하는 것도 아니다.

## Allen 1차 확인 패널 22개

A 축 8개는 랩미팅 가설에 근거한 후보, B 축 8개는 레포의 기존 1차 qPCR 패널, C 축 4개는 전기생리 연결 후보, D 축 2개는 미토콘드리아 확장 후보이다. **Atlas screening panel이며 22개를 모두 qPCR로 바로 수행한다는 의미는 아니다.**

| Gene | 분석 축 | 해마 방향 / log2FC / FDR | 피질 방향 / log2FC / FDR | B6J 500-gene 실측 | B6N 1122-gene 실측 |
|---|---|---|---|---|---|
| Slc29a2 | A 아데노신 수송 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 없음 |
| Slc29a1 | A 아데노신 수송 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 없음 |
| Entpd1 | A 아데노신 생성 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 없음 |
| Nt5e | A 아데노신 생성 | FDR05 표에 없음 | down / -0.5367 / 0.00745 | 없음 | 없음 |
| Adk | A 아데노신 대사 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 없음 |
| Ada | A 아데노신 대사 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 없음 |
| Adora1 | A 아데노신 반응 | FDR05 표에 없음 | FDR05 표에 없음 | 없음 | 있음 |
| Adora2a | A 아데노신 반응 | FDR05 표에 없음 | down / -1.2929 / 0.00932 | 없음 | 있음 |
| Grm2 | B 기존 PBM 1차 | up / +1.1062 / 0.0214 | up / +0.4419 / 1.95e-07 | 없음 | 없음 |
| Ryr2 | B 기존 PBM 1차 | up / +0.6854 / 0.0239 | up / +0.3068 / 0.00237 | 없음 | 없음 |
| Trpc5 | B 기존 PBM 1차 | up / +0.8097 / 0.023 | up / +0.3639 / 0.00272 | 없음 | 있음 |
| Atp2b3 | B 기존 PBM 1차 | up / +0.6907 / 0.0222 | up / +0.2607 / 0.003 | 없음 | 없음 |
| Cacna1c | B 기존 PBM 1차 | up / +1.0542 / 4.28e-07 | up / +0.3065 / 0.0353 | 없음 | 없음 |
| Grin2a | B 기존 PBM 1차 | up / +0.8364 / 0.00078 | up / +0.3480 / 0.0498 | 없음 | 없음 |
| Mrps6 | B 기존 PBM 1차 | down / -0.8467 / 0.0309 | down / -0.3386 / 0.000431 | 없음 | 없음 |
| Penk | B 기존 PBM 1차 | down / -1.8638 / 0.00671 | down / -1.2446 / 7.26e-05 | 있음 | 있음 |
| Kcnj6 | C 전기생리 연결 | up / +0.7272 / 0.015 | up / +0.2514 / 0.0487 | 없음 | 없음 |
| Shisa7 | C 전기생리 연결 | up / +0.7042 / 0.0155 | up / +0.2232 / 0.0359 | 없음 | 없음 |
| Slitrk3 | C 전기생리 연결 | up / +0.5693 / 0.0388 | up / +0.2359 / 0.0142 | 없음 | 없음 |
| Dlgap1 | C 전기생리 연결 | up / +0.3575 / 0.00976 | up / +0.1583 / 0.00321 | 없음 | 없음 |
| Bcs1l | D 미토콘드리아 확장 | FDR05 표에 없음 | up / +0.2525 / 0.00915 | 없음 | 없음 |
| Ppargc1a | D 미토콘드리아 확장 | FDR05 표에 없음 | up / +0.2328 / 0.0388 | 없음 | 없음 |

FDR05 표에 없음은 제공된 유의 DEG 목록에 없음이다. 해마·피질 간 log2FC 또는 baseMean을 동일한 절대 발현 척도로 비교하지 않는다. B6N 패널 실측 유무는 Zhuang-ABCA-1의 공식 gene list에서 확인한 것이다.

### 후보의 기능과 해석 조건

- **Slc29a2 (ENT2 핵심 후보)**: 세포 안팎의 adenosine 수송. 높은 발현만으로 방출 방향·flux를 추정하지 않는다.
- **Slc29a1 (ENT1 비교 축)**: 세포 안팎의 adenosine 수송. ENT2의 상대적 발현과 함께 비교한다.
- **Entpd1 (CD39)**: 세포외 ATP/ADP → AMP. 전사량은 효소 활성·기질 공급을 대신하지 않는다.
- **Nt5e (CD73)**: 세포외 AMP → adenosine. 피질 PBM에서 down; MOp–VISp 기저 차이는 아직 모른다.
- **Adk (adenosine kinase)**: 세포내 adenosine → AMP. 세포내 대사를 extracellular clearance와 구분한다.
- **Ada (adenosine deaminase)**: adenosine → inosine. 관련 compartment와 세포 유형을 확인한다.
- **Adora1 (A1 receptor)**: adenosine 수용체 반응. 발현량을 수면 유도 threshold로 환산하지 않는다.
- **Adora2a (A2A receptor)**: adenosine 수용체 반응. 피질 PBM에서 down; receptor mRNA와 기능은 다르다.
- **Grm2 (mGluR2)**: 글루탐산 조절. up를 일괄 흥분성 증가로 해석하지 않는다.
- **Ryr2 (ryanodine receptor 2)**: ER Ca2+ 방출. 기저 발현과 acute 자극 반응을 구분한다.
- **Trpc5 (TRPC5)**: 양이온 전도. 전사량은 채널 전류가 아니다.
- **Atp2b3 (PMCA3)**: ATP 의존성 Ca2+ 배출. ATP synthase 또는 ATP 생성 유전자가 아니다.
- **Cacna1c (Cav1.2)**: 전압 의존성 Ca2+ 유입. 기저 발현과 PBM 유도변화를 구분한다.
- **Grin2a (GluN2A)**: NMDA 수용체. 수용체 구성·단백질 활성 확인이 별도로 필요하다.
- **Mrps6 (mitochondrial ribosomal protein S6)**: 미토콘드리아 번역. down을 ATP 생산 감소로 해석하지 않는다.
- **Penk (proenkephalin)**: 신경펩타이드. 전체 변화에 세포 유형 비율이 기여하는지 확인한다.
- **Kcnj6 (GIRK2)**: GPCR 관련 K+ 전도. OFF-state를 직접 측정하는 marker가 아니다.
- **Shisa7 (SHISA7)**: GABA 수용체 조절. 흥분성 뉴런의 억제성 입력 관련 발현일 수도 있다.
- **Slitrk3 (SLITRK3)**: 억제성 시냅스 조직화. GABA neuron 수를 대신하는 marker로 쓰지 않는다.
- **Dlgap1 (DLGAP1)**: postsynaptic 조직화. 시냅스 관련 발현을 세포 유형 내에서 비교한다.
- **Bcs1l (BCS1L)**: 호흡사슬 complex III 조립. 공통 DEG가 아니며 단일 유전자에서 호흡을 추정하지 않는다.
- **Ppargc1a (PGC1alpha)**: 미토콘드리아·대사 조절. 세포 유형 구성 차이와 전사량 차이를 구분한다.

Grid1, Btbd9, Gjb6 등은 공통 52개 전체의 탐색 화면에 유지한다. 기존 qPCR 후보 순위는 유지하되 Allen 비교에서는 receptor–GIRK 및 세포 유형의 연결을 보기 위해 Kcnj6, Shisa7, Slitrk3, Dlgap1을 포함한다. 이 재배치는 qPCR 수치 근거가 새로 더 강해졌다는 뜻이 아니다.

## Gene set 비교 방법

Gene-level 결과를 먼저 제시하고 adenosine 관련 8개, Ca²⁺ handling, synaptic/membrane signaling, mitochondrial regulation을 설명용으로 묶는다. 수동으로 묶은 패널을 새로 통계적으로 유의한 pathway라고 부르지 않는다.

기존의 common-up 34와 common-down 18은 방향이 있으므로 별도 집합으로 표시한다. 38개 GO는 부모·자식 항목이 중복되므로 독립된 38개 현상으로 세지 않는다. 기존 GO hit 19개는 기전 후보의 보조 근거로 연결한다. 최초 분석은 dot plot과 세포 유형별 발현 차이로 충분하다. Gene-set score는 탐색 단계에서만 보조적으로 사용하며, 검출된 유전자 수와 각 유전자의 기여를 같이 표시한다. 적은 수의 mitochondrial genes로 ATP 생성이나 호흡량을 판정하지 않는다.

## 발현 분석과 통계 단위

1. 각 후보의 평균 정규화 발현 및 검출 세포 비율을 region × cell type별로 표시한다. Dot plot은 dot 크기 = 검출률, 색 = 정규화 발현으로 만든다. 검출률은 기술·read depth에도 영향을 받는다.
2. 동일 dataset/chemistry, 대응 가능한 genotype·성별·연령·선별 조건에서 비교한다. 전체 세포 평균과 같은 유형 내 평균을 함께 보아 composition과 within-type 차이를 구분한다.
3. 큰 subclass 내부의 세부 cluster 구성이 다른지 확인한다. 필요하면 공유 하위 유형별 비교나 공통 가중치로 표준화한 탐색 결과를 추가한다.
4. 독립 donor 또는 독립 biological pool 정보를 확보할 때 raw counts를 donor/pool × region × cell type으로 합친 pseudobulk 비교를 사용한다. 같은 동물의 두 영역이면 짝지음/동물 효과를 반영한다.
5. pooling을 되돌릴 수 없거나 영역과 batch가 완전히 겹쳐 지역 효과를 분리할 수 없으면 개체 수준 통계 검정을 하지 않는다. 유전자별 효과 방향과 데이터의 제한을 보고한다. 세포·절편·양 반구를 independent animal n으로 세지 않는다.
6. 비교 수가 성립할 때만 사전에 정한 패널·세포 유형 비교 범위에서 다중검정을 보정한다. 최소 표본 수를 임의로 정해 검출력이 충분하다고 가정하지 않는다.

B6J 공간 자료는 한 동물이므로 지역 내 공간 패턴의 기술적 비교가 1차 결과이다. 절편 반복이나 bootstrap은 그 동물 내 공간적 불확실성을 보여줄 수 있지만 동물 간 재현성을 증명하지 않는다. B6N 공간 자료는 별도 결과로 대조한다.

## 실행 순서와 산출물

| 순서 | 실행 | 완료 기준 |
|---|---|---|
| 1 메타데이터 확인 | MOp/VISp label, genotype, donor/pool, sex, age, chemistry, FACS 조건, cell type별 수, gene coverage 확인 | 어떤 세포군에서 직접 비교·통계 검정이 가능한지 표로 확정 |
| 2 세포 구성 | B6J MERFISH에서 영역 전체 및 층별 큰 그룹·뉴런 subclass 비교 | 구성 bar plot, 층별 heatmap, 분모 및 sampled-cell count 표시 |
| 3 후보 발현 | CTX-HPF 자료에서 22개 패널, 공통 52개는 탐색 범위로 확인 | cell type × gene dot plot 및 대응 유형 내 MOp–VISp 차이 표 |
| 4 공간 및 자료 간 확인 | 실제 자극 ROI, 경계 제외, measured ISH/MERFISH, 별도 자료에서 방향 확인 | 상위 후보가 전체 영역 평균에만 의존하는지 판정 |
| 5 실험 우선순위 | 후보·세포 유형·후속 readout 연결 | 우선 검증할 후보 3–5개와 각 후보가 판별할 가설 제시 |

초기 스크리닝에서 세포 구성과 후보 발현을 병행한다. Glia 후보가 해당 RNA-seq 자료에서 충분히 표집되지 않으면 더 넓은 VIS 자료로 바꾸어 VISp라고 보고하지 않는다. 그 세포군은 보조 발현 참고 및 실제 동물의 위치·세포 유형을 유지하는 검증 대상으로 남긴다.

## 결과가 다음 실험으로 이어지는 기준

| 확인되는 패턴 | 해석할 수 있는 범위 | 우선 다음 단계 |
|---|---|---|
| Slc29a2, Nt5e, Adk 등에서 같은 유형 내 지역 차이가 유지됨 | adenosine 순축적 차이를 검증할 분자 후보 | M1–V1 국소 extracellular adenosine 시간 경과 비교 및 후보 위치 확인 |
| Adora1/Adora2a 및 Kcnj6 등에서 지역 차이가 유지됨 | adenosine 입력에 대한 반응 차이를 검증할 후보 | 같은 adenosine 입력에 대한 LFP 및 sleep 반응 비교 |
| 특정 subtype 비율이 다르고 해당 subtype에 후보가 집중됨 | composition이 전체 평균에 기여할 가능성 | 자극 위치 주변의 세포 구성과 후보 발현을 함께 검증 |
| 후보 발현 차이가 작거나 다른 자료에서 일관되지 않음 | 특정 기저 전사체 차이의 근거가 약함 | PBM dose·광전달과 자극 유도 dynamics/adenosine을 직접 비교 |

GRAB-Ado는 extracellular adenosine의 상대적 변화, LFP는 국소 전기생리 변화를 각각 측정한다. 서로 다른 부위·동물의 dF/F를 동일한 절대 adenosine 농도로 간주하지 않는다. 단일 LFP의 slow-wave 지표를 neuronal OFF occupancy와 동일시하지 않는다. 따라서 현재 2채널+EMG 조건에서 가능한 local LFP·상태 분석과 연결하되 atlas 결과가 OFF-state나 threshold를 직접 입증한다고 표현하지 않는다.

실제 검증에서 sham/PBM과 MOp/VISp를 비교한다면 region × treatment 차이를 평가한다. 한 영역에서 유의하고 다른 영역에서 비유의하다는 것만으로 지역 간 처리 효과가 다르다고 판정하지 않는다. 채취 ZT, 자극 후 시간, 수면/각성 상태 및 실제 국소 광량을 반영한다.

## 공통 52개 전체 탐색 목록

공통 up 34개: Arsb, Atp2b3, Bcl7a, Btbd3, Btbd9, Cacna1c, Ccbe1, Cd276, Cdk14, Dlgap1, Epha4, Faah, Grid1, Grin2a, Grm2, Hsd11b1, Irs1, Kcnj6, Map2k6, Npnt, Pcdh19, Proz, Ptpn3, Ryr2, Samd12, Shisa7, Slitrk3, Sphkap, Stra6, Stxbp6, Tnr, Trpc5, Wnt9a, Zeb2

공통 down 18개: Adipor2, Arid5b, Arpp21, Cdkn1c, Chrm4, Ecel1, Gimap6, Gjb6, Ier5, Ifitm2, Mrps6, Nupr1, Penk, S100a10, Sap30, Slc10a4, Sorbs3, Tmem40

## 주요 출처

- [레포 최신 DEG 검증 보고서](https://github.com/Sonnetly/tPBM_opendataset/blob/e430e8d2d2827ccdd222e8828d77227e937b6f1b/ora_analysis/deg_validation/DEG_aware_qPCR_report_ko.md)
- [레포 전체 선택 DEG 합집합](https://github.com/Sonnetly/tPBM_opendataset/blob/e430e8d2d2827ccdd222e8828d77227e937b6f1b/ora_analysis/deg_validation/all_selected_symbol_union.csv)
- [Allen B6J MERFISH 자료](https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850.html)
- [Allen B6J 실측 gene list](https://alleninstitute.github.io/abc_atlas_access/_downloads/54f8c94c59bb7c322c5f933e0b151681/gene_list.html)
- [Zhuang 실측 gene list](https://alleninstitute.github.io/abc_atlas_access/_downloads/94cc4adeccbb444b45046cd9ea170138/gene_list.html)
- [Allen WMB 10x 메타데이터·접근 안내](https://alleninstitute.github.io/abc_atlas_access/notebooks/10x_snRNASeq_tutorial_part_1.html)
- [Allen CTX-HPF 10x 자료](https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data/mouse-whole-cortex-and-hippocampus-10x)
- [MOp GEO sample GSM5629156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5629156)
- [VISp GEO sample GSM5629158](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5629158)
- [Yao et al. Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8195859/)
- [Yao et al. Nature 2023](https://www.nature.com/articles/s41586-023-06812-z)
- [Zhang et al. Nature 2023](https://www.nature.com/articles/s41586-023-06808-9)
- [Zhang et al. Nature 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8494645/)
- [Tasic et al. Nature 2018](https://www.nature.com/articles/s41586-018-0654-5)
- [Zhou et al. Cell Research 2024 40 Hz light flickering promotes sleep through cortical adenosine signaling](https://pmc.ncbi.nlm.nih.gov/articles/PMC10907382/)

