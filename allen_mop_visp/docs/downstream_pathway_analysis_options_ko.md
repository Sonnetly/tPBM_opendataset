# MOp–VISp downstream pathway 후속 분석 제안

사용자가 선택한 후속 방향은 **MOp·VISp에서 NREM 관련 영역으로 향하는 투사(projection)를 먼저 확인**하고, 아데노신 관련 분자 경로를 별도 단계로 연결하는 것이다. 이 문서는 분석 방법과 해석 기준이다. 새 분자 pathway 점수·ORA·GSEA·TF 분석은 실행하지 않았다. 실제 투사 계산과 결과는 별도 connectivity 산출물에서 확인한다.

## 우선 실행: NREM 관련 영역으로의 투사

Allen Mouse Brain Connectivity Atlas의 anterograde tracer 자료에서 MOp와 VISp injection 실험을 추출하고, injection 중심·영역 순도·부피·세포 유형·hemisphere·QC를 확인한다. 주 비교는 조건이 가능한 한 유사한 실험 안에서 수행하며, WT와 Cre line 실험을 구분한다. 실험 수와 조건 차이를 숨기지 않고 개별 injection 결과를 보존한다.

표적은 NREM을 촉진한 **실험적 세포 집단과 해당 해부학적 영역을 분리하여** 정의한다. 우선 VLPO·MnPO를 평가하고, ZI Lhx6 집단과 basal forebrain SOM 집단은 세포 유형 특이 근거를 붙인 추가 표적으로 다룬다. Parafacial zone(PZ)은 Allen ontology에 정확한 영역이 없으면 미평가로 남기며 다른 뇌간 영역으로 치환하지 않는다. Thalamic reticular nucleus는 NREM 관련 추가 목록으로 구분한다. 기능 문헌과 실제 ontology ID는 `NREM_projection_target_definition_ko.md` 및 `nrem_target_definitions` 파일에 기록한다.

각 target의 ipsilateral·contralateral projection density/volume 및 공식 정규화 지표를 확인하고, injection volume·source purity·threshold·target definition에 따른 민감도를 함께 제시한다. 전체 target label로 들어오는 신호는 그 영역의 수면촉진 세포에 도달했다는 증거가 아니다. 투사 신호는 시냅스 수·기능적 연결 강도·NREM 유발·PBM 유발 전파를 직접 측정하지 않는다. atlas의 source injection 위치와 사용자 자극 ROI의 일치 여부도 별도 확인해야 한다.

권장 산출물은 source injection/QC 목록, experiment별 target 지표, MOp·VISp 조건별 요약, 해부학적 커버리지와 미평가 표적, projection heatmap·개별 injection 민감도 그림이다. 기능 해석은 **관찰된 해부학적 투사 → 표적 세포 유형 → 실제 NREM 효과**의 각 단계에서 근거와 미확인 사항을 구분한다. 현재 RNA atlas만으로 downstream 회로를 확정할 수 없다.

## 현재 자료로 가능한 범위

10x 원본에는 31,053개 유전자, SMART-seq 원본에는 45,768개 exon-count feature가 있다. 그러나 보관한 두 H5AD는 **선택된 후보 62개만 포함**한다. 전체 유전자 library 합계를 알고 있어도 발현 행렬이 전장 자료가 되는 것은 아니다. 이 62개를 배경으로 GO·GSEA·GSVA·TF 분석을 하면 후보 선별에 따른 편향이 생긴다.

원본 캐시에서 필요한 세포의 전체 유전자 요약 또는 사전에 정한 경로의 **모든 구성 유전자**를 추가 추출해야 한다. 대형 행렬을 메모리에 한꺼번에 펼치기보다 chunk 단위로 정규화·집계한다. 정규화 분모는 기존처럼 전체 원본 유전자 합계를 유지한다. 10x UMI와 SMART-seq read counts는 별도로 분석한다.

10x의 정확한 조건 비교는 P58에서 가능하지만 지역별 library-preparation date가 완전히 분리되어 지역 효과와 날짜 효과를 구분할 수 없다. donor별 요약이나 batch 보정만으로 이 문제를 해결할 수 없다. 현재 결과는 **관측한 표본의 서술적 차이**이며, 지역의 생물학적 효과에 대한 P value/FDR를 제시할 설계가 아니다. SMART-seq도 donor·batch·해부학적 표집의 제약이 있고, 표시한 비교 조건은 P53/54/55/56/59여서 주 10x P58 조건의 직접 재현 자료가 아니다.

## 우선 살펴볼 분자 축

| 분석 축 | 질문과 포함할 범위 | 별도로 확인해야 할 기능 |
|---|---|---|
| Adenosine 생성·수송·대사 | ENTPD1–NT5E, SLC29A1/2, ADK, ADA와 관련 반응의 전체 구성 유전자 | 세포외 adenosine의 생성·이동·제거 및 순축적. 수송체 발현만으로 수송 방향을 정하지 않는다. |
| A1 receptor 하위 신호 | ADORA1–Gi/o–adenylyl cyclase–cAMP/PKA, Gβγ–GIRK 및 Ca²⁺ 채널 조절 | cAMP/PKA 변화, GIRK 전류와 Ca²⁺ 반응. mRNA 점수를 수용체 활성으로 환산하지 않는다. |
| A2A receptor 하위 신호 | ADORA2A–Gs 계열–cAMP/PKA/CREB | 수용체·cAMP·인산화 반응. 현재 10x의 낮은 Adora2a 검출을 신호 부재로 판정하지 않는다. |
| 미토콘드리아·대사 | OXPHOS, ATP 관련 반응, 산화 스트레스 대응의 충분한 유전자 집합 | ATP·호흡·ROS 측정. 전사량과 실제 ATP 생성량 또는 세포외 adenosine 농도는 다른 값이다. |
| Ca²⁺·시냅스·막전위 | Ca²⁺ 수송·저장·채널, glutamate/GABA 신호와 시냅스 조절 | Ca²⁺ dynamics, 시냅스 전류, 국소 LFP. baseline atlas를 PBM 처리 후 변화로 해석하지 않는다. |

경로는 Reactome/MSigDB 등의 버전이 명시된 유전자 집합을 우선 사용하고, 위 가설 축의 수동 모듈은 별도로 표시한다. 중첩되는 경로를 독립된 현상 여러 개로 세지 않는다. 원본 PBM의 공통 52개·GO hit 19개는 연결 근거로 주석에 유지하며, 새 MOp–VISp pathway 유의성의 근거로 사용하지 않는다.

## 권장 실행 순서

1. **유전자 범위와 경로 정의를 고정한다.** Mouse gene ID/symbol, 자료별 측정 여부, 경로 출처·버전·구성 유전자·중복 매핑을 저장한다. ORA/GSEA의 배경은 해당 assay와 세포 유형에서 측정·분석 가능한 전체 유전자다. 발현 필터를 쓴다면 두 영역에 동일한 사전 기준을 적용한다.
2. **Gene-level 결과와 모듈 점수를 함께 본다.** P58 10x의 동일 subclass에서 각 경로 구성 유전자의 평균 정규화 발현·검출률을 먼저 확인한다. 이어 기준 집단과 control-gene 선택이 명시된 탐색적 module score를 계산한다. 경로별 coverage와 유전자별 기여를 함께 표시하여 몇 개 유전자가 점수를 지배하는지 확인한다. 모든 모듈은 같은 세포 집합과 정규화 기준에서 계산한다.
3. **공유 subtype/cluster의 구성 효과를 점검한다.** 두 영역에 공유되는 cluster만 사용한 공통 가중치 결과를 추가하고, cell-weighted 결과와 방향·크기가 유지되는지 비교한다. 가중치와 제외된 세포·cluster 수를 저장한다. donor/library별 관측값과 검출률도 확인하되 이를 batch 효과를 제거한 지역 효과라고 부르지 않는다. 희소 glia/vascular 결과는 reporter/FACS 표집 범위 안에서만 해석한다.
4. **필요하면 전장 유전자 순위로 탐색을 넓힌다.** 동일 조건·세포 유형의 전체 유전자에 대해 `MOp − VISp 평균 log1p 정규화 발현` 같은 기술적 효과 순위를 만든다. 그 순위로 preranked GSEA의 enrichment score와 leading-edge 유전자를 살핀다. 현재 1차 보고서는 유의성 판정 없이 점수·기여 유전자·민감도를 제시한다. 유전자집합 permutation의 nominal P/FDR가 산출되더라도 이는 유전자 순위 귀무모형에 대한 값이며 독립 동물에 근거한 지역 차이의 증거가 아니다.
5. **별도 자료와 기능 검증을 연결한다.** SMART-seq은 같은 자료 안에서 reporter·sex·정확한 age·FACS·subclass를 맞춘 보조 결과로 유지한다. 두 assay의 발현값이나 점수를 합치지 않는다. 표본 조건이 다른 방향 비교는 제한을 붙여 제시한다. 우선 실험은 경로별 분자 후보와 판별할 기능 readout을 연결해 정한다.

GSVA/ssGSEA는 전장 발현 범위를 확보한 뒤 보조 점수로 사용할 수 있다. 점수가 높다는 사실은 경로의 실제 활성이나 flux를 뜻하지 않는다. TF/regulon 분석도 decoupler 등의 타깃 네트워크와 충분한 전장 유전자 coverage를 확보한 뒤 탐색적으로 수행하며, TF mRNA와 TF 활성 추정을 구분한다. CellChat 같은 발현 기반 상호작용 추정은 인과적 통신 증거가 아니고, 대사물질인 adenosine의 생성·방출·농도를 직접 모델링하는 분석을 대신하지 않는다.

## 계획한 산출물과 검증

| 산출물 | 반드시 포함할 내용 |
|---|---|
| `gene_universe_and_coverage.csv`, `module_definitions.tsv` | assay별 전체 배경, 경로 출처·버전, 구성 유전자, ID 대응, 미측정·중복 처리 |
| `module_celltype_summary.csv` | 영역·조건·세포 유형·모듈 점수, gene-level 발현·검출률·기여, 세포 수 |
| `module_cluster_sensitivity.csv` | 공유 cluster, 동일 가중치, 제외 수, donor/library별 기술 요약, 방향 변화·미평가 구분 |
| `exploratory_preranked_gene_sets.csv` | 전장 순위 정의, enrichment score, leading edge, coverage, 생물학적 유의성 검정 미실행 표시 |
| 모듈 heatmap·주요 유전자 dot plot·해석 보고서 | 분모와 조건, 핵심 유전자, 표집·batch 한계, 다음 기능 실험 |

선택 세포 ID의 원본 대응, 전체 library 분모, 측정 유전자 범위, 동일 세포 유형 매칭을 검증한다. 일부 세포의 원본 counts와 점수 기여를 다시 계산하고, 전체 평균과 공유 cluster 결과를 비교한다. 미측정 유전자나 공유 cluster 부재를 0 또는 방향 반전으로 바꾸지 않는다. 원본·경로 목록·코드·산출물의 버전과 checksum을 남긴다.

## 지역별 PBM 기전을 검증하려면

후속 자료는 WT B6J에서 **region × sham/PBM**을 포함하고, 독립 동물 반복과 연령·성별·ZT·수면 상태·처리 후 시간·prep batch를 균형 있게 배치하는 것이 우선이다. 양쪽 영역을 같은 동물에서 확보한다면 짝지음 정보를 보존한다. 지역별 처리 효과의 차이는 region × treatment 항으로 평가하며, 한 영역만 유의하다는 이유로 지역 효과가 다르다고 판정하지 않는다.

핵심 기능 검증은 국소 extracellular adenosine dynamics, A1/A2A 관련 약리·cAMP/인산화 반응, GIRK/Ca²⁺ 또는 LFP readout을 조합하는 방향이다. 실제 자극 ROI와 광량·깊이가 확인되어야 전체 영역 atlas와 자극 주변 결과를 연결할 수 있다.

## 참고 자료

- Reactome 경로·분석 안내: <https://reactome.org/> / <https://reactome.org/userguide/analysis>
- MSigDB 유전자 집합: <https://www.gsea-msigdb.org/gsea/msigdb/>
- GSEA 공식 문서: <https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html>
- decoupler 공식 문서: <https://decoupler.readthedocs.io/>
- Allen Mouse Brain Connectivity Atlas: <https://connectivity.brain-map.org/>

경로 ID와 유전자 목록은 실제 분석 시 공식 자료에서 확인하여 고정한다. 이 제안서에는 확인하지 않은 pathway ID나 새로운 분석 수치를 넣지 않았다.
