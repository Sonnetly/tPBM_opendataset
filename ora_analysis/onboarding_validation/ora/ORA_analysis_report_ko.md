# ORA 전체 catalog 기반 유전자 방향 일치 분석

> 이 보고서는 ORA hits에서 복원한 범위의 분석이다. 원본 DEG FDR·log2FC를 반영한 최신 overlap 및 qPCR 후보는 [최신 결과 안내](README.md)의 DEG 검증 보고서를 참고한다.

**확인 가능한 같은 방향 overlap은 51개(공통 up 34개, 공통 down 17개)입니다.** 이전에 언급된 7개는 이온채널·칼슘 수송 기능군의 부분집합일 가능성이 있습니다. 기존 7개 목록 자체는 제공되지 않았습니다.

이 결과는 첨부 ORA 표에서 복원되는 유전자 기준입니다. 완전한 DEG 표가 없어 전체 DEG 기준 최종 overlap을 확정할 수 없습니다. 두 데이터셋의 비교 방향은 파일의 up/down을 그대로 사용했습니다.

## 사용 범위와 누락

| 데이터셋 | 방향 | 표에 기록된 DEG 수 | 이름 복원 | 미복원 |
|---|---|---:|---:|---:|
| GSE245658_hippocampus_808nm | down | 292 | 259 | 33 |
| GSE245658_hippocampus_808nm | up | 354 | 278 | 76 |
| GSE325653_prefrontal_1064nm | down | 527 | 484 | 43 |
| GSE325653_prefrontal_1064nm | up | 471 | 419 | 52 |

- CSV의 28,092행 전체를 사용하고, 각 hits의 세미콜론 구분 유전자 이름을 모아 중복 제거했습니다. 경로 FDR로 필터하지 않았습니다.
- 복원된 유전자 중 두 데이터셋에 공통으로 있는 것은 74개: 같은 방향 51개, 반대 방향 23개입니다. 한 데이터셋에만 보이는 경우도 비-DEG라고 판정하지 않았습니다.
- A에서 up/B에서 down은 12개, A에서 down/B에서 up은 11개입니다. 전체 복원 유전자의 합집합 1,366개도 엑셀에 포함했습니다.
- gene-level log2FC, adjusted P, 원래 DEG 선정 기준, contrast 정의 및 전체 배경 유전자 이름은 이 CSV에 없습니다.

## 과거 7개와 전체 51개의 관계

**Atp2b3, Cacna1c, Grid1, Grin2a, Kcnj6, Ryr2, Trpc5**는 모두 공통 up입니다. 현재 표의 GO:0005216(ion channel activity)에 나오는 공통 유전자는 Atp2b3를 제외한 6개이며, 칼슘 수송/항상성과 KEGG calcium signaling에 연결되는 Atp2b3를 함께 묶으면 7개입니다. 이 분류가 과거의 7개와 같다는 확정은 아닙니다.

FDR<0.05인 GO/KEGG 항목에 속한 유전자만 각각 모아 교집합을 구하는 정의에서는 26개(22 up, 4 down)가 나옵니다. 따라서 ‘FDR를 적용하면 항상 7개’라는 설명은 이 파일에서 재현되지 않습니다. FDR 조건과 특정 기능군 조건은 별개입니다.

## 51개 유전자별 기능과 근거

아래 근거는 업로드된 ORA 표의 gene–term 매핑입니다. KEGG는 경로, GO BP는 생물학적 과정, GO MF는 분자 기능입니다. 링크는 공식 항목 참조용이며 이번 환경에서는 외부 페이지를 직접 검증하지 못했습니다.

| 유전자 | 공통 방향 | 기능군 | 주석에 근거한 해석 | KEGG/GO 근거 |
|---|---|---|---|---|
| Ccbe1 | up | ECM·접착·혈관·면역 | 림프관 형성 및 혈관 형성 | [GO:0001946](https://www.ebi.ac.uk/QuickGO/term/GO:0001946) — lymphangiogenesis; [GO:0001525](https://www.ebi.ac.uk/QuickGO/term/GO:0001525) — angiogenesis |
| Cd276 | up | ECM·접착·혈관·면역 | 세포접착분자 경로 및 면역반응 조절 관련 GO | [mmu04514](https://www.kegg.jp/pathway/mmu04514) — Cell adhesion molecule (CAM) interaction; [GO:0002764](https://www.ebi.ac.uk/QuickGO/term/GO:0002764) — immune response-regulating signaling pathway |
| Npnt | up | ECM·접착·혈관·면역 | ECM-수용체 및 integrin 신호; 세포-기질 접착 | [mmu04512](https://www.kegg.jp/pathway/mmu04512) — ECM-receptor interaction; [mmu04518](https://www.kegg.jp/pathway/mmu04518) — Integrin signaling; [GO:0007160](https://www.ebi.ac.uk/QuickGO/term/GO:0007160) — cell-matrix adhesion |
| Pcdh19 | up | ECM·접착·혈관·면역 | cadherin 신호 및 동종 세포 간 접착 | [mmu04519](https://www.kegg.jp/pathway/mmu04519) — Cadherin signaling; [GO:0007156](https://www.ebi.ac.uk/QuickGO/term/GO:0007156) — homophilic cell-cell adhesion |
| Proz | up | ECM·접착·혈관·면역 | 혈액 응고 및 지혈 | [GO:0007596](https://www.ebi.ac.uk/QuickGO/term/GO:0007596) — blood coagulation; [GO:0007599](https://www.ebi.ac.uk/QuickGO/term/GO:0007599) — hemostasis |
| Hsd11b1 | up | 대사·호르몬·MAPK | 스테로이드 호르몬 생합성 및 glucocorticoid 대사 | [mmu00140](https://www.kegg.jp/pathway/mmu00140) — Steroid hormone biosynthesis; [GO:0008211](https://www.ebi.ac.uk/QuickGO/term/GO:0008211) — glucocorticoid metabolic process |
| Irs1 | up | 대사·호르몬·MAPK | 인슐린, PI3K-Akt 및 AMPK 신호 경로 | [mmu04910](https://www.kegg.jp/pathway/mmu04910) — Insulin signaling pathway; [mmu04151](https://www.kegg.jp/pathway/mmu04151) — PI3K-Akt signaling pathway; [mmu04152](https://www.kegg.jp/pathway/mmu04152) — AMPK signaling pathway |
| Map2k6 | up | 대사·호르몬·MAPK | MAPK 및 TNF 신호 경로 | [mmu04010](https://www.kegg.jp/pathway/mmu04010) — MAPK signaling pathway; [mmu04668](https://www.kegg.jp/pathway/mmu04668) — TNF signaling pathway |
| Arsb | up | 분해·영양소 수송 | glycosaminoglycan 분해 및 lysosome 관련 경로 | [mmu00531](https://www.kegg.jp/pathway/mmu00531) — Glycosaminoglycan degradation; [mmu04142](https://www.kegg.jp/pathway/mmu04142) — Lysosome biogenesis |
| Stra6 | up | 분해·영양소 수송 | 비타민 수송 및 막 수송체 활성 | [GO:0051180](https://www.ebi.ac.uk/QuickGO/term/GO:0051180) — vitamin transport; [GO:0090482](https://www.ebi.ac.uk/QuickGO/term/GO:0090482) — vitamin transmembrane transporter activity |
| Btbd9 | up | 수면·일주기·기억 | 수면, 일주기 리듬 및 기억 관련 GO | [GO:0030431](https://www.ebi.ac.uk/QuickGO/term/GO:0030431) — sleep; [GO:0007623](https://www.ebi.ac.uk/QuickGO/term/GO:0007623) — circadian rhythm; [GO:0007613](https://www.ebi.ac.uk/QuickGO/term/GO:0007613) — memory |
| Dlgap1 | up | 시냅스 전달·조직화 | 글루탐산성 시냅스 및 postsynaptic density 조직화 | [mmu04724](https://www.kegg.jp/pathway/mmu04724) — Glutamatergic synapse; [GO:0097106](https://www.ebi.ac.uk/QuickGO/term/GO:0097106) — postsynaptic density organization |
| Faah | up | 시냅스 전달·조직화 | 역행성 endocannabinoid 신호 경로 및 지방산 분해 | [mmu04723](https://www.kegg.jp/pathway/mmu04723) — Retrograde endocannabinoid signaling; [GO:0009062](https://www.ebi.ac.uk/QuickGO/term/GO:0009062) — fatty acid catabolic process |
| Grm2 | up | 시냅스 전달·조직화 | G protein 결합 글루탐산 수용체 신호 및 adenylate cyclase 억제 관련 GO | [mmu04724](https://www.kegg.jp/pathway/mmu04724) — Glutamatergic synapse; [GO:0007216](https://www.ebi.ac.uk/QuickGO/term/GO:0007216) — G protein-coupled glutamate receptor signaling pathway; [GO:0007193](https://www.ebi.ac.uk/QuickGO/term/GO:0007193) — adenylate cyclase-inhibiting G protein-coupled receptor signaling pathway |
| Shisa7 | up | 시냅스 전달·조직화 | GABA 신호, 수용체 군집화 및 시냅스 조직화 | [GO:0007214](https://www.ebi.ac.uk/QuickGO/term/GO:0007214) — gamma-aminobutyric acid signaling pathway; [GO:0043113](https://www.ebi.ac.uk/QuickGO/term/GO:0043113) — receptor clustering; [GO:0050808](https://www.ebi.ac.uk/QuickGO/term/GO:0050808) — synapse organization |
| Slitrk3 | up | 시냅스 전달·조직화 | GABA성 시냅스 전달 및 시냅스 조립 | [GO:0051932](https://www.ebi.ac.uk/QuickGO/term/GO:0051932) — synaptic transmission, GABAergic; [GO:0007416](https://www.ebi.ac.uk/QuickGO/term/GO:0007416) — synapse assembly; [mmu04514](https://www.kegg.jp/pathway/mmu04514) — Cell adhesion molecule (CAM) interaction |
| Stxbp6 | up | 시냅스 전달·조직화 | 소포 융합 및 SNARE 복합체 조립 | [GO:0006906](https://www.ebi.ac.uk/QuickGO/term/GO:0006906) — vesicle fusion; [GO:0035493](https://www.ebi.ac.uk/QuickGO/term/GO:0035493) — SNARE complex assembly |
| Btbd3 | up | 신경 연결·발달 | 수상돌기 형태 형성 및 신경 분화 관련 세포 형태 형성 | [GO:0048813](https://www.ebi.ac.uk/QuickGO/term/GO:0048813) — dendrite morphogenesis; [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) — cell morphogenesis involved in neuron differentiation |
| Epha4 | up | 신경 연결·발달 | 축삭 유도 및 시냅스 구조/활성 조절 | [mmu04360](https://www.kegg.jp/pathway/mmu04360) — Axon guidance; [GO:0050803](https://www.ebi.ac.uk/QuickGO/term/GO:0050803) — regulation of synapse structure or activity |
| Tnr | up | 신경 연결·발달 | ECM-수용체 상호작용, 신경세포 접착 및 시냅스 가소성 조절 | [mmu04512](https://www.kegg.jp/pathway/mmu04512) — ECM-receptor interaction; [GO:0007158](https://www.ebi.ac.uk/QuickGO/term/GO:0007158) — neuron cell-cell adhesion; [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) — regulation of synaptic plasticity |
| Zeb2 | up | 신경 연결·발달 | 축삭 형성 및 신경 분화 관련 세포 형태 형성 | [GO:0007409](https://www.ebi.ac.uk/QuickGO/term/GO:0007409) — axonogenesis; [GO:0048667](https://www.ebi.ac.uk/QuickGO/term/GO:0048667) — cell morphogenesis involved in neuron differentiation |
| Ptpn3 | up | 신호 조절 단백질 | 단백질 tyrosine phosphatase 및 이온채널 조절 기능 | [GO:0004725](https://www.ebi.ac.uk/QuickGO/term/GO:0004725) — protein tyrosine phosphatase activity; [GO:0099106](https://www.ebi.ac.uk/QuickGO/term/GO:0099106) — ion channel regulator activity |
| Samd12 | up | 신호 조절 단백질 | 수용체 tyrosine kinase 신호 관련 GO; 이 파일의 기능 근거가 제한적 | [GO:0007169](https://www.ebi.ac.uk/QuickGO/term/GO:0007169) — cell surface receptor protein tyrosine kinase signaling pathway |
| Sphkap | up | 신호 조절 단백질 | protein kinase A 결합; 이 파일에서는 MF 주석만 확인됨 | [GO:0051018](https://www.ebi.ac.uk/QuickGO/term/GO:0051018) — protein kinase A binding |
| Bcl7a | up | 염색질·Wnt·세포주기 | ATP 의존성 염색질 재구성과 시냅스 가소성 조절 | [mmu03082](https://www.kegg.jp/pathway/mmu03082) — ATP-dependent chromatin remodeling; [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) — regulation of synaptic plasticity |
| Cdk14 | up | 염색질·Wnt·세포주기 | canonical Wnt 신호 조절 및 세포주기 G2/M 전환 | [GO:0060828](https://www.ebi.ac.uk/QuickGO/term/GO:0060828) — regulation of canonical Wnt signaling pathway; [GO:0000086](https://www.ebi.ac.uk/QuickGO/term/GO:0000086) — G2/M transition of mitotic cell cycle |
| Wnt9a | up | 염색질·Wnt·세포주기 | Wnt 신호 경로 | [mmu04310](https://www.kegg.jp/pathway/mmu04310) — Wnt signaling pathway; [GO:0016055](https://www.ebi.ac.uk/QuickGO/term/GO:0016055) — Wnt signaling pathway |
| Atp2b3 | up | 이온채널·칼슘 수송 | 칼슘 수송 및 세포 내 칼슘 항상성; 채널과 구분되는 칼슘 펌프 | [mmu04020](https://www.kegg.jp/pathway/mmu04020) — Calcium signaling pathway; [GO:0006874](https://www.ebi.ac.uk/QuickGO/term/GO:0006874) — intracellular calcium ion homeostasis |
| Cacna1c | up | 이온채널·칼슘 수송 | 전압 의존성 이온채널; 칼슘 신호 및 glutamatergic synapse | [GO:0005244](https://www.ebi.ac.uk/QuickGO/term/GO:0005244) — voltage-gated monoatomic ion channel activity; [mmu04020](https://www.kegg.jp/pathway/mmu04020) — Calcium signaling pathway; [mmu04724](https://www.kegg.jp/pathway/mmu04724) — Glutamatergic synapse |
| Grid1 | up | 이온채널·칼슘 수송 | 글루탐산 수용체/채널 관련 주석 및 시냅스 가소성 조절 | [GO:0008066](https://www.ebi.ac.uk/QuickGO/term/GO:0008066) — glutamate receptor activity; [GO:0005216](https://www.ebi.ac.uk/QuickGO/term/GO:0005216) — monoatomic ion channel activity; [GO:0048167](https://www.ebi.ac.uk/QuickGO/term/GO:0048167) — regulation of synaptic plasticity |
| Grin2a | up | 이온채널·칼슘 수송 | 글루탐산 수용체; 칼슘 신호 및 glutamatergic synapse | [GO:0008066](https://www.ebi.ac.uk/QuickGO/term/GO:0008066) — glutamate receptor activity; [mmu04020](https://www.kegg.jp/pathway/mmu04020) — Calcium signaling pathway; [mmu04724](https://www.kegg.jp/pathway/mmu04724) — Glutamatergic synapse |
| Kcnj6 | up | 이온채널·칼슘 수송 | 칼륨 수송; 시냅스 막전위 및 endocannabinoid 관련 경로 | [GO:0006813](https://www.ebi.ac.uk/QuickGO/term/GO:0006813) — potassium ion transport; [GO:0099505](https://www.ebi.ac.uk/QuickGO/term/GO:0099505) — regulation of presynaptic membrane potential; [mmu04723](https://www.kegg.jp/pathway/mmu04723) — Retrograde endocannabinoid signaling |
| Ryr2 | up | 이온채널·칼슘 수송 | 칼슘 신호 및 circadian entrainment 관련 이온채널 | [GO:0005216](https://www.ebi.ac.uk/QuickGO/term/GO:0005216) — monoatomic ion channel activity; [mmu04020](https://www.kegg.jp/pathway/mmu04020) — Calcium signaling pathway; [mmu04713](https://www.kegg.jp/pathway/mmu04713) — Circadian entrainment |
| Trpc5 | up | 이온채널·칼슘 수송 | 칼슘 수송 및 축삭 유도 관련 이온채널 | [GO:0006816](https://www.ebi.ac.uk/QuickGO/term/GO:0006816) — calcium ion transport; [GO:0005216](https://www.ebi.ac.uk/QuickGO/term/GO:0005216) — monoatomic ion channel activity; [mmu04360](https://www.kegg.jp/pathway/mmu04360) — Axon guidance |
| Adipor2 | down | 대사·미토콘드리아 | adipocytokine 및 AMPK 신호 경로 | [mmu04920](https://www.kegg.jp/pathway/mmu04920) — Adipocytokine signaling pathway; [mmu04152](https://www.kegg.jp/pathway/mmu04152) — AMPK signaling pathway |
| Mrps6 | down | 대사·미토콘드리아 | 리보솜 및 미토콘드리아 번역 | [mmu03010](https://www.kegg.jp/pathway/mmu03010) — Ribosome; [GO:0032543](https://www.ebi.ac.uk/QuickGO/term/GO:0032543) — mitochondrial translation |
| Gimap6 | down | 스트레스·면역·자가포식 | 자가포식 및 macroautophagy | [GO:0006914](https://www.ebi.ac.uk/QuickGO/term/GO:0006914) — autophagy; [GO:0016236](https://www.ebi.ac.uk/QuickGO/term/GO:0016236) — macroautophagy |
| Ier5 | down | 스트레스·면역·자가포식 | 열에 대한 세포 반응 | [GO:0034605](https://www.ebi.ac.uk/QuickGO/term/GO:0034605) — cellular response to heat |
| Ifitm2 | down | 스트레스·면역·자가포식 | I형 및 II형 인터페론 반응 | [GO:0034340](https://www.ebi.ac.uk/QuickGO/term/GO:0034340) — response to type I interferon; [GO:0034341](https://www.ebi.ac.uk/QuickGO/term/GO:0034341) — response to type II interferon |
| Nupr1 | down | 스트레스·면역·자가포식 | 산화적 인산화 조절 및 DNA 손상 반응성 intrinsic apoptosis 관련 GO | [GO:0002082](https://www.ebi.ac.uk/QuickGO/term/GO:0002082) — regulation of oxidative phosphorylation; [GO:0008630](https://www.ebi.ac.uk/QuickGO/term/GO:0008630) — intrinsic apoptotic signaling pathway in response to DNA damage |
| Chrm4 | down | 신경전달·세포 간 연결 | 콜린성 시냅스 및 adenylate cyclase 억제형 GPCR 신호 | [mmu04725](https://www.kegg.jp/pathway/mmu04725) — Cholinergic synapse; [GO:0007193](https://www.ebi.ac.uk/QuickGO/term/GO:0007193) — adenylate cyclase-inhibiting G protein-coupled receptor signaling pathway |
| Ecel1 | down | 신경전달·세포 간 연결 | 신경펩타이드 신호 및 단백질 processing | [GO:0007218](https://www.ebi.ac.uk/QuickGO/term/GO:0007218) — neuropeptide signaling pathway; [GO:0016485](https://www.ebi.ac.uk/QuickGO/term/GO:0016485) — protein processing |
| Gjb6 | down | 신경전달·세포 간 연결 | 전기적 coupling에 의한 세포 간 통신 및 세포 접합 | [GO:0010644](https://www.ebi.ac.uk/QuickGO/term/GO:0010644) — cell communication by electrical coupling; [GO:0007043](https://www.ebi.ac.uk/QuickGO/term/GO:0007043) — cell-cell junction assembly |
| Penk | down | 신경전달·세포 간 연결 | 신경펩타이드 신호 및 neuroactive ligand-receptor interaction | [GO:0007218](https://www.ebi.ac.uk/QuickGO/term/GO:0007218) — neuropeptide signaling pathway; [mmu04080](https://www.kegg.jp/pathway/mmu04080) — Neuroactive ligand-receptor interaction |
| Slc10a4 | down | 신경전달·세포 간 연결 | 신경전달물질 수송 조절 및 나트륨 수송 | [GO:0051588](https://www.ebi.ac.uk/QuickGO/term/GO:0051588) — regulation of neurotransmitter transport; [GO:0006814](https://www.ebi.ac.uk/QuickGO/term/GO:0006814) — sodium ion transport |
| Arid5b | down | 전사·RNA·세포주기 | DNA 기반 전사 및 RNA polymerase II 전사 개시 | [GO:0006351](https://www.ebi.ac.uk/QuickGO/term/GO:0006351) — DNA-templated transcription; [GO:0006367](https://www.ebi.ac.uk/QuickGO/term/GO:0006367) — transcription initiation at RNA polymerase II promoter |
| Arpp21 | down | 전사·RNA·세포주기 | mRNA 및 mRNA 3'-UTR 결합 | [GO:0003730](https://www.ebi.ac.uk/QuickGO/term/GO:0003730) — mRNA 3'-UTR binding; [GO:0003729](https://www.ebi.ac.uk/QuickGO/term/GO:0003729) — mRNA binding |
| Cdkn1c | down | 전사·RNA·세포주기 | 세포주기 및 세포주기 과정의 음성 조절 | [mmu04110](https://www.kegg.jp/pathway/mmu04110) — Cell cycle; [GO:0010948](https://www.ebi.ac.uk/QuickGO/term/GO:0010948) — negative regulation of cell cycle process |
| Sap30 | down | 전사·RNA·세포주기 | 전사 corepressor 활성 | [GO:0003714](https://www.ebi.ac.uk/QuickGO/term/GO:0003714) — transcription corepressor activity |
| S100a10 | down | 접착·세포골격 | 세포-기질 접착 및 plasminogen 활성화 조절 | [GO:0001954](https://www.ebi.ac.uk/QuickGO/term/GO:0001954) — positive regulation of cell-matrix adhesion; [GO:0010755](https://www.ebi.ac.uk/QuickGO/term/GO:0010755) — regulation of plasminogen activation |
| Sorbs3 | down | 접착·세포골격 | 세포-기질 접착 및 actin stress fiber 조립 조절 | [GO:0031589](https://www.ebi.ac.uk/QuickGO/term/GO:0031589) — cell-substrate adhesion; [GO:0051496](https://www.ebi.ac.uk/QuickGO/term/GO:0051496) — positive regulation of stress fiber assembly |

## 경로 수준에서 상대적으로 강한 공통 근거

아래 4개 KEGG 경로는 각 데이터셋의 up DEG 전체로 수행된 기존 ORA에서 모두 FDR<0.05입니다. **51개 overlap만으로 새로 검정한 enrichment가 아닙니다.** FDR는 유전자 자체의 값이 아닙니다.

| KEGG 경로 | 같은 방향 공통 유전자 | A ORA FDR | B ORA FDR |
|---|---|---:|---:|
| [Circadian entrainment (mmu04713)](https://www.kegg.jp/pathway/mmu04713) | Cacna1c, Grin2a, Kcnj6, Ryr2 | 1.23e-05 | 4.33e-06 |
| [Calcium signaling pathway (mmu04020)](https://www.kegg.jp/pathway/mmu04020) | Atp2b3, Cacna1c, Grin2a, Ryr2 | 0.000648 | 0.0066 |
| [Glutamatergic synapse (mmu04724)](https://www.kegg.jp/pathway/mmu04724) | Cacna1c, Dlgap1, Grin2a, Grm2 | 0.0404 | 0.00086 |
| [Retrograde endocannabinoid signaling (mmu04723)](https://www.kegg.jp/pathway/mmu04723) | Cacna1c, Faah, Kcnj6 | 0.0422 | 0.0104 |

공통 up에서 GO 유의 항목도 38개가 겹칩니다. 예: synapse organization(GO:0050808), regulation of synaptic plasticity(GO:0048167), axon guidance(GO:0007411). GO는 부모·자식 항목이 중복되므로 38개의 독립적인 현상으로 해석하지 않습니다. 공통 유의 GO 한 항목(GO:1990138)은 두 데이터셋의 hits 교집합이 0개입니다. 경로가 공통인 것과 그 경로의 유전자가 공통인 것은 구분해야 합니다.

공통 down에서는 두 데이터셋 모두 FDR<0.05인 동일 KEGG/GO 항목이 없습니다. 따라서 down 17개에 대해서는 기능적 연관성을 제시할 수 있지만, 공통적으로 유의한 경로 감소라고 표현할 수 없습니다.

## 해석 시 주의할 점

- 발현 up과 신경 흥분성 증가를 동일시하지 않습니다. 이온 수송 방향, 수용체 조절, 소포/시냅스 조직화 역할이 섞여 있습니다. 예를 들어 Grm2에는 adenylate cyclase-inhibiting GPCR GO가, Shisa7/Slitrk3에는 GABA 관련 GO가 함께 있습니다.
- 입력 라벨상 해마 808 nm와 전전두피질 1064 nm 데이터입니다. 두 연구의 contrast, 처리 조건과 표본 구성은 확인되지 않았으므로 공통 PBM 기전이나 효과의 재현을 단정하지 않습니다.
- KEGG의 심장질환·암·감염 경로에 연결되는 공유 유전자가 있어도 해당 조직에서 그 질환이 발생하거나 개선된다는 증거는 아닙니다. 본 요약은 기능적으로 직접적인 경로/GO를 우선 제시했습니다.
- 51개 집합 자체의 enrichment/overlap 유의성 검정, 효과크기 비교, gene-level FDR 재판정은 수행하지 않았습니다.

## 수치 검증과 산출물

28,092개 행에서 k와 hits 수가 일치했고 dataset×direction별 hypergeometric P를 입력 집계값으로 재계산해 일치함을 확인했습니다. GO+KEGG 전체 catalog를 함께 BH 보정한 값도 최대 절대오차 2.3e-16 이내에서 재현되었습니다. 이는 입력된 집계값의 내부 일관성 검증이며 원래 DEG 선정이나 annotation 정확성을 외부 검증한 것은 아닙니다.

- `ORA_direction_overlap_analysis.xlsx`: 전체 복원 유전자, 공통 51개, 반대방향 23개, 각 유전자에 대한 모든 주석/원본 행 번호/기존 ORA 통계, 유의 경로 비교 및 원본 전체 표.
- `concordant_51.csv`: 같은 방향 51개 요약. `all_recovered_genes.csv`: 복원 가능한 모든 유전자.
- `primary_evidence.csv`: 요약에 사용한 유전자별 KEGG/GO 근거와 두 데이터셋의 기존 ORA 수치.
- `all_concordant_annotations.csv`: 51개에 연결된 전체 KEGG/GO 주석.
- `analyze_ora.py`, `gene_function_map.tsv`: 재현용 분석 코드와 설명용 기능 분류.

완전한 DEG 기준으로 마무리하려면 두 데이터셋의 원본 DEG 표(gene ID/symbol, log2FC, adjusted P와 비교 조건)가 필요합니다. 현재 결과의 미복원 유전자 이름은 ORA 표에서 추정하여 채우지 않았습니다.
