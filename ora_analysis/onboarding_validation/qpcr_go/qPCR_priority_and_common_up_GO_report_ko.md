# qPCR 우선 후보와 공통 up 유의 GO

> 이 보고서는 gene-level 통계를 연결하기 전의 ORA 기반 후보 제안이다. DEG FDR·log2FC를 반영한 최신 방향 검증 패널은 [DEG-aware 보고서](../deg_validation/DEG_aware_qPCR_report_ko.md)를 참고한다.

## 사용자 가설과 우선 제안

가설: **Regional differences in the sleep response to tPBM arise from differences in local cortical responsiveness to stimulation.** qPCR 시료는 cortex 중심이며, 영역별 전기생리 반응과 이온 신호/미토콘드리아 기전을 탐색한다.

**1차 패널 8개: Cacna1c, Grin2a, Atp2b3, Kcnj6, Ryr2, Grm2, Ppargc1a, Mrps6.**

앞의 6개는 두 데이터셋 공통 up인 이온/글루탐산 관련 후보다. Ppargc1a는 전전두피질 up hits에서만 관찰되고, Mrps6은 두 연구 공통 down인 미토콘드리아 번역 관련 후보다. 이 8개 전부가 공통 up이라는 뜻은 아니다.

6개로 줄이면 Cacna1c, Grin2a, Atp2b3, Kcnj6, Ppargc1a, Mrps6을 제안한다. 12개로 확장하면 Shisa7, Slitrk3, Dlgap1, Bcs1l을 추가한다. Per2/Cry2는 일주기 가설을 별도로 검증할 때 고려할 추가 후보이고, Chrm4/Penk는 추가 down 검증 대상이다. 이들은 내부 정규화 reference gene이 아니다.

현재 표는 ORA 결과로, log2FC/유전자별 FDR/평균 발현량/분산이 없다. 우선순위는 통계적으로 가장 강한 DEG나 실제 qPCR 검출 용이성을 뜻하지 않는다. 공통 발현 방향과 기능 대표성을 바탕으로 한 잠정적인 실험 순서다. 시료의 종·조직·처리 조건과 원본 DEG를 확인하면 조정해야 한다.

| 순서 | 후보 | 방향 | 패널 | 선택 이유 | 대표 근거 |
|---:|---|---|---|---|---|
| 1 | Cacna1c | up | A_cortex_priority_8 | 공통 유의 gated channel GO와 calcium signaling를 근거로 피질의 칼슘/전기생리 반응 관련 발현을 우선 확인 | [GO:0022836](https://www.ebi.ac.uk/QuickGO/term/GO:0022836) (양쪽 ORA 유의); [GO:0005244](https://www.ebi.ac.uk/QuickGO/term/GO:0005244) (양쪽 ORA 유의); [mmu04020](https://www.kegg.jp/pathway/mmu04020) (양쪽 ORA 유의) |
| 2 | Grin2a | up | A_cortex_priority_8 | 공통 glutamatergic synapse, calcium signaling 및 시냅스 가소성 GO를 함께 대표 | [mmu04724](https://www.kegg.jp/pathway/mmu04724) (양쪽 ORA 유의); [mmu04020](https://www.kegg.jp/pathway/mmu04020) (양쪽 ORA 유의); [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) (양쪽 ORA 유의) |
| 3 | Atp2b3 | up | A_cortex_priority_8 | 공통 calcium signaling와 금속 이온 수송체 GO를 근거로 채널 외 칼슘 수송/항상성 관련 발현을 확인 | [mmu04020](https://www.kegg.jp/pathway/mmu04020) (양쪽 ORA 유의); [GO:0046873](https://www.ebi.ac.uk/QuickGO/term/GO:0046873) (양쪽 ORA 유의); [GO:0006874](https://www.ebi.ac.uk/QuickGO/term/GO:0006874) (기능 주석; 양쪽 유의 조건 미충족) |
| 4 | Kcnj6 | up | A_cortex_priority_8 | 공통 막전위 조절 GO와 endocannabinoid KEGG를 근거로 칼륨 관련 신호를 포함 | [GO:0099505](https://www.ebi.ac.uk/QuickGO/term/GO:0099505) (양쪽 ORA 유의); [mmu04723](https://www.kegg.jp/pathway/mmu04723) (양쪽 ORA 유의); [GO:0006813](https://www.ebi.ac.uk/QuickGO/term/GO:0006813) (기능 주석; 양쪽 유의 조건 미충족) |
| 5 | Ryr2 | up | A_cortex_priority_8 | 공통 calcium signaling와 이온채널 GO에 연결되는 추가 칼슘 관련 후보 | [mmu04020](https://www.kegg.jp/pathway/mmu04020) (양쪽 ORA 유의); [GO:0022836](https://www.ebi.ac.uk/QuickGO/term/GO:0022836) (양쪽 ORA 유의) |
| 6 | Grm2 | up | A_cortex_priority_8 | 공통 glutamatergic synapse와 channel regulator GO를 근거로 단순 흥분성 증가 해석을 보완하는 조절형 수용체 후보 | [mmu04724](https://www.kegg.jp/pathway/mmu04724) (양쪽 ORA 유의); [GO:0016247](https://www.ebi.ac.uk/QuickGO/term/GO:0016247) (양쪽 ORA 유의); [GO:0007193](https://www.ebi.ac.uk/QuickGO/term/GO:0007193) (기능 주석; 양쪽 유의 조건 미충족) |
| 7 | Ppargc1a | A:not_observed;B:up | A_cortex_priority_8 | 전전두피질 up hits에서 관찰된 미토콘드리아 관련 조절 후보; 대표 미토 GO는 양쪽 유의하지 않아 가설 검증 목적으로 포함 | [GO:0007005](https://www.ebi.ac.uk/QuickGO/term/GO:0007005) (기능 주석; 양쪽 유의 조건 미충족); [GO:0090140](https://www.ebi.ac.uk/QuickGO/term/GO:0090140) (기능 주석; 양쪽 유의 조건 미충족) |
| 8 | Mrps6 | down | A_cortex_priority_8 | 두 데이터셋 공통 down이며 미토콘드리아 번역에 연결됨; ATP 생산 자체의 지표나 유의한 공통 미토 경로로 해석하지 않음 | [GO:0032543](https://www.ebi.ac.uk/QuickGO/term/GO:0032543) (기능 주석; 양쪽 유의 조건 미충족); [GO:0140053](https://www.ebi.ac.uk/QuickGO/term/GO:0140053) (기능 주석; 양쪽 유의 조건 미충족) |
| 9 | Shisa7 | up | B_extension_4 | 공통 시냅스 조직화/가소성과 GABA 기능 주석을 근거로 전기생리적 조절 가설을 확장 | [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) (양쪽 ORA 유의); [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) (양쪽 ORA 유의); [GO:0007214](https://www.ebi.ac.uk/QuickGO/term/GO:0007214) (기능 주석; 양쪽 유의 조건 미충족) |
| 10 | Slitrk3 | up | B_extension_4 | 공통 시냅스 조직화와 신경 분화 관련 GO 및 GABA성 전달 주석을 보완 | [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) (양쪽 ORA 유의); [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) (양쪽 ORA 유의); [GO:0051932](https://www.ebi.ac.uk/QuickGO/term/GO:0051932) (기능 주석; 양쪽 유의 조건 미충족) |
| 11 | Dlgap1 | up | B_extension_4 | 공통 synapse organization과 glutamatergic synapse에 연결되어 시냅스 구조 관련 발현을 확인 | [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) (양쪽 ORA 유의); [GO:0099173](https://www.ebi.ac.uk/QuickGO/term/GO:0099173) (양쪽 ORA 유의); [mmu04724](https://www.kegg.jp/pathway/mmu04724) (양쪽 ORA 유의) |
| 12 | Bcs1l | A:not_observed;B:up | B_extension_4 | 전전두피질 up hits에서 관찰된 직접적인 호흡사슬 조립 관련 후보; 해당 GO FDR는 유의하지 않음 | [GO:0034551](https://www.ebi.ac.uk/QuickGO/term/GO:0034551) (기능 주석; 양쪽 유의 조건 미충족) |
| 13 | Epha4 | up | C_other_shared_GO_10 | 공통 유의 axon guidance GO의 공통 hit가 Epha4 하나여서 신경 연결 변화 확인에 유용 | [GO:0007411](https://www.ebi.ac.uk/QuickGO/term/GO:0007411) (양쪽 ORA 유의); [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) (양쪽 ORA 유의) |
| 14 | Tnr | up | C_other_shared_GO_10 | 공통 시냅스 조직화/가소성과 ECM 기능 주석을 연결하는 후보 | [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) (양쪽 ORA 유의); [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) (양쪽 ORA 유의); [mmu04512](https://www.kegg.jp/pathway/mmu04512) (기능 주석; 양쪽 유의 조건 미충족) |
| 15 | Grid1 | up | C_other_shared_GO_10 | 공통 glutamate receptor activity와 synapse structure/activity 조절 GO를 보완 | [GO:0008066](https://www.ebi.ac.uk/QuickGO/term/GO:0008066) (양쪽 ORA 유의); [GO:0050803](https://www.ebi.ac.uk/QuickGO/term/GO:0050803) (양쪽 ORA 유의) |
| 16 | Trpc5 | up | C_other_shared_GO_10 | 공통 이온채널 GO와 axonogenesis 조절 GO를 보완 | [GO:0005216](https://www.ebi.ac.uk/QuickGO/term/GO:0005216) (양쪽 ORA 유의); [GO:0050770](https://www.ebi.ac.uk/QuickGO/term/GO:0050770) (양쪽 ORA 유의) |
| 17 | Zeb2 | up | C_other_shared_GO_10 | 공통 신경 분화 관련 형태 형성과 developmental growth GO의 대표 후보 | [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) (양쪽 ORA 유의); [GO:0060560](https://www.ebi.ac.uk/QuickGO/term/GO:0060560) (양쪽 ORA 유의) |
| 18 | Bcl7a | up | C_other_shared_GO_10 | 공통 시냅스 가소성 GO와 염색질 재구성 주석을 연결 | [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) (양쪽 ORA 유의); [mmu03082](https://www.kegg.jp/pathway/mmu03082) (기능 주석; 양쪽 유의 조건 미충족) |
| 19 | Btbd3 | up | C_other_shared_GO_10 | 공통 신경 분화 관련 형태 형성 GO와 수상돌기 형태 형성 주석을 보완 | [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) (양쪽 ORA 유의); [GO:0048813](https://www.ebi.ac.uk/QuickGO/term/GO:0048813) (기능 주석; 양쪽 유의 조건 미충족) |
| 20 | Btbd9 | up | C_other_shared_GO_10 | 공통 cognition/learning or memory GO와 일주기/수면 주석을 보완; 급성 수면 반응의 특이적 지표로 확정하지 않음 | [GO:0050890](https://www.ebi.ac.uk/QuickGO/term/GO:0050890) (양쪽 ORA 유의); [GO:0007611](https://www.ebi.ac.uk/QuickGO/term/GO:0007611) (양쪽 ORA 유의); [GO:0007623](https://www.ebi.ac.uk/QuickGO/term/GO:0007623) (기능 주석; 양쪽 유의 조건 미충족) |
| 21 | Ptpn3 | up | C_other_shared_GO_10 | 공통 이온 막수송 조절 및 channel regulator GO를 보완 | [GO:0034765](https://www.ebi.ac.uk/QuickGO/term/GO:0034765) (양쪽 ORA 유의); [GO:0016247](https://www.ebi.ac.uk/QuickGO/term/GO:0016247) (양쪽 ORA 유의) |
| 22 | Stra6 | up | C_other_shared_GO_10 | 공통 cognition/learning or memory GO에 포함되며 기능 주석은 비타민 수송에 연결됨 | [GO:0050890](https://www.ebi.ac.uk/QuickGO/term/GO:0050890) (양쪽 ORA 유의); [GO:0007611](https://www.ebi.ac.uk/QuickGO/term/GO:0007611) (양쪽 ORA 유의); [GO:0051180](https://www.ebi.ac.uk/QuickGO/term/GO:0051180) (기능 주석; 양쪽 유의 조건 미충족) |
| 23 | Faah | up | D_KEGG_specific_up_1 | 공통 유의 endocannabinoid KEGG에는 포함되지만 공통 유의 GO의 공통 hit에는 포함되지 않는 별도 기전 후보 | [mmu04723](https://www.kegg.jp/pathway/mmu04723) (양쪽 ORA 유의) |
| 24 | Per2 | A:not_observed;B:down | E_optional_clock_2 | 전전두피질 down hits와 유의 circadian rhythm GO에서 관찰; 채취 시간/ZT를 통제하는 별도 일주기 가설용 후보 | [GO:0007623](https://www.ebi.ac.uk/QuickGO/term/GO:0007623) (기능 주석; 양쪽 유의 조건 미충족) |
| 25 | Cry2 | A:up;B:down | E_optional_clock_2 | 전전두피질 down·해마 up으로 관찰되나 연구/파장이 달라 영역 특이성을 입증하지 못함; 시간 통제하에 탐색 | [GO:0007623](https://www.ebi.ac.uk/QuickGO/term/GO:0007623) (기능 주석; 양쪽 유의 조건 미충족) |
| 26 | Chrm4 | down | F_optional_down_2 | 공통 down 방향을 확인하기 위한 신경 관련 추가 후보; 양쪽 모두 유의한 동일 GO/KEGG 근거는 없음 | [mmu04725](https://www.kegg.jp/pathway/mmu04725) (기능 주석; 양쪽 유의 조건 미충족); [GO:0007193](https://www.ebi.ac.uk/QuickGO/term/GO:0007193) (기능 주석; 양쪽 유의 조건 미충족) |
| 27 | Penk | down | F_optional_down_2 | 공통 down 방향을 확인하기 위한 신경펩타이드 추가 후보; 양쪽 모두 유의한 동일 GO/KEGG 근거는 없음 | [GO:0007218](https://www.ebi.ac.uk/QuickGO/term/GO:0007218) (기능 주석; 양쪽 유의 조건 미충족); [mmu04080](https://www.kegg.jp/pathway/mmu04080) (기능 주석; 양쪽 유의 조건 미충족) |

## 미토콘드리아 가설의 현재 근거

현재 데이터에서는 이온채널/시냅스 관련 ORA 근거가 강하고 미토콘드리아 경로의 유의한 근거는 약하다. 전전두피질의 KEGG Oxidative phosphorylation(mmu00190)은 up/down 모두 k=0, FDR=1이다. 이름에 mitochond/oxidative phosphorylation/respiratory chain/cellular respiration/electron transport가 들어가는 전전두피질 ORA 항목을 추출한 부분집합에서도 FDR<0.05 항목은 없다. 이는 정의한 항목명 검색의 결과이며 모든 미토콘드리아 기능을 완전히 포괄하는 검정은 아니다.

Ppargc1a는 전전두피질 up이며 GO:0090140(regulation of mitochondrial fission)의 FDR는 약 0.604다. Bcs1l도 전전두피질 up이며 GO:0034551(mitochondrial respiratory chain complex III assembly)의 FDR는 약 0.491이다. Mrps6은 공통 down이지만 mitochondrial translation GO의 FDR는 두 연구 모두 1이다. 이 수치는 유전자별 유의확률이 아니다. 세 후보는 데이터에서 관찰되는 가설 관련 표적이며, 유의한 미토 경로 활성/억제의 증거로 제안한 것이 아니다.

Cox4i1/Tfam/Sod2 등 추가 기전 표적은 현재 ORA hits에서 관찰되지 않는다. 이것이 비-DEG이거나 발현하지 않음을 뜻하지는 않는다. 후보를 추가하려면 원본 DEG/발현량과 별도 주석 근거를 확인해야 한다.

## 공통 up GO 추출 기준

두 데이터셋의 up ORA에 같은 GO ID가 존재하고, 각 기존 full-catalog BH FDR가 모두 0.05 미만인 항목을 선택했다. 유전자 연결은 해당 GO의 두 hits 목록에 실제로 같은 유전자가 있는지로 확인했다. 최대 FDR값은 정렬 편의를 위한 값이며 새로 결합한 P/FDR가 아니다.

**38개 항목: Biological Process 29개, Molecular Function 9개.** 이 중 공통 hit가 있는 항목은 37개, 그 유전자 합집합은 19개다. GO:1990138(neuron projection extension)은 두 데이터셋 모두 유의하지만 공통 hit는 0개다.

공통 유의 GO의 공통 hit 19개: Atp2b3, Bcl7a, Btbd3, Btbd9, Cacna1c, Dlgap1, Epha4, Grid1, Grin2a, Grm2, Kcnj6, Ptpn3, Ryr2, Shisa7, Slitrk3, Stra6, Tnr, Trpc5, Zeb2.

## 기능군별 요약

아래 묶음은 설명용 수동 분류다. GO ontology를 재구성하거나 독립적 경로로 검정한 것이 아니다.

| 기능군 | GO 항목 수 | 공통 유전자 합집합 | 우선 패널에서 포함한 유전자 |
|---|---:|---|---|
| 시냅스 조직화·구조 | 6 | Dlgap1, Epha4, Grid1, Grin2a, Shisa7, Slitrk3, Tnr | Grin2a |
| 시냅스 신호·가소성 | 6 | Bcl7a, Cacna1c, Dlgap1, Epha4, Grid1, Grin2a, Grm2, Shisa7, Slitrk3, Tnr | Cacna1c, Grin2a, Grm2 |
| 이온채널·수송·막전위 | 14 | Atp2b3, Cacna1c, Grid1, Grin2a, Grm2, Kcnj6, Ptpn3, Ryr2, Shisa7, Trpc5 | Atp2b3, Cacna1c, Grin2a, Grm2, Kcnj6, Ryr2 |
| 신경 형태·축삭 발달 | 9 | Btbd3, Epha4, Slitrk3, Tnr, Trpc5, Zeb2 |  |
| 인지·학습·행동 | 3 | Bcl7a, Btbd9, Cacna1c, Epha4, Grin2a, Grm2, Shisa7, Stra6, Tnr | Cacna1c, Grin2a, Grm2 |

## 공통 up 유의 GO 38개 전체

| GO | 항목 | 분류 | A ORA FDR | B ORA FDR | 공통 up 유전자 |
|---|---|---|---:|---:|---|
| [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) | synapse organization | biological_process | 8.121e-08 | 7.96e-08 | Dlgap1, Epha4, Grin2a, Shisa7, Slitrk3, Tnr |
| [GO:0022836](https://www.ebi.ac.uk/QuickGO/term/GO:0022836) | gated channel activity | molecular_function | 0.0001232 | 6.639e-05 | Cacna1c, Grid1, Grin2a, Kcnj6, Ryr2 |
| [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) | cell morphogenesis involved in neuron differentiation | biological_process | 1.25e-05 | 0.0004004 | Btbd3, Epha4, Slitrk3, Zeb2 |
| [GO:0061564](https://www.ebi.ac.uk/QuickGO/term/GO:0061564) | axon development | biological_process | 5.998e-05 | 0.0007145 | Epha4, Slitrk3, Zeb2 |
| [GO:0007409](https://www.ebi.ac.uk/QuickGO/term/GO:0007409) | axonogenesis | biological_process | 1.232e-05 | 0.0008082 | Epha4, Slitrk3, Zeb2 |
| [GO:0050803](https://www.ebi.ac.uk/QuickGO/term/GO:0050803) | regulation of synapse structure or activity | biological_process | 6.815e-06 | 0.0008082 | Epha4, Grid1, Slitrk3 |
| [GO:0022832](https://www.ebi.ac.uk/QuickGO/term/GO:0022832) | voltage-gated channel activity | molecular_function | 0.0004623 | 0.0008602 | Cacna1c, Grin2a, Kcnj6 |
| [GO:0005216](https://www.ebi.ac.uk/QuickGO/term/GO:0005216) | monoatomic ion channel activity | molecular_function | 1.232e-05 | 0.001443 | Cacna1c, Grid1, Grin2a, Kcnj6, Ryr2, Trpc5 |
| [GO:0046873](https://www.ebi.ac.uk/QuickGO/term/GO:0046873) | metal ion transmembrane transporter activity | molecular_function | 5.87e-06 | 0.001448 | Atp2b3, Cacna1c, Grin2a, Kcnj6, Ryr2, Trpc5 |
| [GO:0050890](https://www.ebi.ac.uk/QuickGO/term/GO:0050890) | cognition | biological_process | 0.0006284 | 0.002457 | Bcl7a, Btbd9, Cacna1c, Grin2a, Shisa7, Stra6, Tnr |
| [GO:0050807](https://www.ebi.ac.uk/QuickGO/term/GO:0050807) | regulation of synapse organization | biological_process | 1.232e-05 | 0.002519 | Epha4, Grid1, Slitrk3 |
| [GO:0099537](https://www.ebi.ac.uk/QuickGO/term/GO:0099537) | trans-synaptic signaling | biological_process | 0.004284 | 2.129e-05 | Cacna1c, Dlgap1, Grid1, Grin2a, Slitrk3 |
| [GO:0005244](https://www.ebi.ac.uk/QuickGO/term/GO:0005244) | voltage-gated monoatomic ion channel activity | molecular_function | 0.0003679 | 0.005835 | Cacna1c, Grin2a, Kcnj6 |
| [GO:0007611](https://www.ebi.ac.uk/QuickGO/term/GO:0007611) | learning or memory | biological_process | 0.007485 | 0.001269 | Btbd9, Cacna1c, Grin2a, Shisa7, Stra6, Tnr |
| [GO:0007268](https://www.ebi.ac.uk/QuickGO/term/GO:0007268) | chemical synaptic transmission | biological_process | 0.007838 | 6.643e-05 | Cacna1c, Dlgap1, Grid1, Grin2a, Slitrk3 |
| [GO:0098916](https://www.ebi.ac.uk/QuickGO/term/GO:0098916) | anterograde trans-synaptic signaling | biological_process | 0.007838 | 6.643e-05 | Cacna1c, Dlgap1, Grid1, Grin2a, Slitrk3 |
| [GO:0022843](https://www.ebi.ac.uk/QuickGO/term/GO:0022843) | voltage-gated monoatomic cation channel activity | molecular_function | 0.0002156 | 0.00992 | Cacna1c, Grin2a, Kcnj6 |
| [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) | regulation of synaptic plasticity | biological_process | 0.01343 | 6.643e-05 | Bcl7a, Epha4, Grid1, Grin2a, Grm2, Shisa7, Tnr |
| [GO:0007411](https://www.ebi.ac.uk/QuickGO/term/GO:0007411) | axon guidance | biological_process | 0.007574 | 0.01352 | Epha4 |
| [GO:0097485](https://www.ebi.ac.uk/QuickGO/term/GO:0097485) | neuron projection guidance | biological_process | 0.007838 | 0.01419 | Epha4 |
| [GO:0043270](https://www.ebi.ac.uk/QuickGO/term/GO:0043270) | positive regulation of monoatomic ion transport | biological_process | 0.01635 | 0.01189 | Cacna1c |
| [GO:0099536](https://www.ebi.ac.uk/QuickGO/term/GO:0099536) | synaptic signaling | biological_process | 0.01635 | 4.799e-06 | Cacna1c, Dlgap1, Grid1, Grin2a, Slitrk3 |
| [GO:0099173](https://www.ebi.ac.uk/QuickGO/term/GO:0099173) | postsynapse organization | biological_process | 8.023e-05 | 0.01835 | Dlgap1, Grin2a, Shisa7, Slitrk3 |
| [GO:0099084](https://www.ebi.ac.uk/QuickGO/term/GO:0099084) | postsynaptic specialization organization | biological_process | 0.01867 | 0.01846 | Dlgap1, Slitrk3 |
| [GO:0099505](https://www.ebi.ac.uk/QuickGO/term/GO:0099505) | regulation of presynaptic membrane potential | biological_process | 0.02386 | 0.0008082 | Grin2a, Kcnj6 |
| [GO:0005261](https://www.ebi.ac.uk/QuickGO/term/GO:0005261) | monoatomic cation channel activity | molecular_function | 3.157e-06 | 0.02452 | Cacna1c, Grin2a, Kcnj6, Ryr2, Trpc5 |
| [GO:0060078](https://www.ebi.ac.uk/QuickGO/term/GO:0060078) | regulation of postsynaptic membrane potential | biological_process | 0.0008507 | 0.02883 | Grid1, Grin2a |
| [GO:0060560](https://www.ebi.ac.uk/QuickGO/term/GO:0060560) | developmental growth involved in morphogenesis | biological_process | 0.01777 | 0.02883 | Zeb2 |
| [GO:0022898](https://www.ebi.ac.uk/QuickGO/term/GO:0022898) | regulation of transmembrane transporter activity | biological_process | 0.04038 | 0.01035 | Grin2a, Shisa7 |
| [GO:0034765](https://www.ebi.ac.uk/QuickGO/term/GO:0034765) | regulation of monoatomic ion transmembrane transport | biological_process | 0.04038 | 0.01352 | Cacna1c, Grin2a, Ptpn3, Ryr2 |
| [GO:0051966](https://www.ebi.ac.uk/QuickGO/term/GO:0051966) | regulation of synaptic transmission, glutamatergic | biological_process | 0.04189 | 0.001516 | Grin2a, Grm2, Tnr |
| [GO:1990138](https://www.ebi.ac.uk/QuickGO/term/GO:1990138) | neuron projection extension | biological_process | 0.04189 | 0.01352 | 없음: 경로만 공통 |
| [GO:0050770](https://www.ebi.ac.uk/QuickGO/term/GO:0050770) | regulation of axonogenesis | biological_process | 0.04353 | 3.445e-05 | Epha4, Tnr, Trpc5, Zeb2 |
| [GO:0008066](https://www.ebi.ac.uk/QuickGO/term/GO:0008066) | glutamate receptor activity | molecular_function | 0.0206 | 0.04409 | Grid1, Grin2a, Grm2 |
| [GO:0099174](https://www.ebi.ac.uk/QuickGO/term/GO:0099174) | regulation of presynapse organization | biological_process | 0.02139 | 0.04585 | Slitrk3 |
| [GO:0030534](https://www.ebi.ac.uk/QuickGO/term/GO:0030534) | adult behavior | biological_process | 0.04701 | 0.001491 | Btbd9, Cacna1c, Epha4, Grm2 |
| [GO:0048588](https://www.ebi.ac.uk/QuickGO/term/GO:0048588) | developmental cell growth | biological_process | 0.01343 | 0.04781 | Zeb2 |
| [GO:0016247](https://www.ebi.ac.uk/QuickGO/term/GO:0016247) | channel regulator activity | molecular_function | 0.03189 | 0.04787 | Grm2, Ptpn3 |

## 실험 해석과 사전 확인

- qPCR은 전사체 발현 방향의 재현을 확인한다. 신경 흥분성, 채널 전류, 단백질 활성, 시냅스 기능이나 PBM 인과 효과를 직접 입증하지 않는다.
- 새 qPCR 시료에서 원 연구와 비교 방향·조직·조건이 대응하는지 확인한다. 다른 시료에서 같은 방향이 반드시 나올 것이라고 예측할 수는 없다.
- 영역별 반응성을 보려면 같은 자극 조건에서 각 cortex 영역의 sham/control과 tPBM을 비교한다. 영역 간 baseline 차이만으로 반응성 차이를 판정하지 않는다. region×treatment 상호작용을 평가하고, 한 동물에서 여러 영역을 측정했다면 동물 내 짝지음을 분석에 반영한다.
- 두 원 연구는 조직(해마/전전두피질)과 파장(808/1064 nm)이 함께 달라, 연구 간 유전자 차이로 cortex 영역별 반응성을 입증할 수 없다.
- 채취 시점, 자극 후 경과 시간, 수면/각성 상태와 circadian time/ZT를 통제한다. Per2/Cry2 발현은 특히 시간 맥락과 함께 해석하며, 두 시점 또는 영역의 차이만으로 수면 반응을 단정하지 않는다.
- 가능하면 원본 gene-level log2FC/FDR와 평균 발현량을 확인해 실험 표적을 조정한다. 현재 파일에서는 실제 DEG 기준을 재판정할 수 없다.
- 내부 정규화 유전자는 실제 시료에서 별도로 안정성을 검증한다. GAPDH/ACTB 등의 안정성은 이 ORA 표나 hits 부재만으로 판단하지 않는다. 적절한 생물학적 반복과 assay 특이성/효율 확인이 필요하다.
- 표적별 primer/assay는 아직 선정하거나 검증하지 않았다. 종과 transcript/isoform 정보를 확보한 뒤 선택해야 한다.
- GO 부모·자식 항목과 동일 hit 목록은 중복될 수 있어 항목 수를 독립적 근거 수로 세지 않는다. A/B 수치는 각 연구의 기존 DEG 집합에 대한 ORA 값이며 이번 qPCR 패널만으로 새로 검정하지 않았다.
- 대표 기능 주석이 있어도 그 특정 항목이 두 연구 모두에서 유의한 것은 아닐 수 있다. `primary_candidate_evidence.csv`의 유의 여부와 FDR를 함께 확인한다.
- 현재 근거는 첨부 CSV의 KEGG/GO 주석이다. 공식 링크는 참조용이며 최신 외부 데이터/문헌의 직접 검증은 수행하지 못했다.

## 재실행과 파일

저장소 루트에서 `python3 ora_analysis/prioritize_qpcr_and_go.py`를 실행한다. 기존 ORA 분석 파일은 그대로 두고 `ora_analysis/qpcr_go/`에 출력한다.

- `qPCR_priority_and_common_up_GO.xlsx`: 후보, 전체 GO/BP/MF, 19개 공통 유전자, 기능군, 유전자별 근거.
- `qpcr_candidates.csv`: cortex 가설 기준 8개 우선 + 4개 확장 + 나머지 공통 GO 후보 10개 + Faah + 일주기 후보 2개 + 추가 down 2개.
- `common_up_GO_significant.csv`: 38개 전체와 두 연구 FDR 및 실제 hits.
- `common_up_GO_BP.csv`, `common_up_GO_MF.csv`: BP/MF 분리 결과.
- `candidate_and_concordant_evidence.csv`: 기존 같은 방향 51개와 추가 cortex/일주기 후보의 공통 유의 항목 연결 여부 및 각 연구 관찰 방향.
- `cortex_mito_related_ORA.csv`: 정의한 미토콘드리아 관련 항목명 검색으로 추출한 전전두피질 ORA 부분집합.
- `primary_candidate_evidence.csv`: 수동 선정의 각 대표 근거와 원본 행 번호.
- 후보 설정은 `ora_analysis/qpcr_candidate_config.tsv`에 저장했다. 다른 데이터셋에는 설정/검증 조건을 수정해야 한다.
