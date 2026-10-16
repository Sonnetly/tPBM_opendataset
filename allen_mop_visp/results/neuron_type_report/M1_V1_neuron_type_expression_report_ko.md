# M1(MOp)과 V1(VISp) 뉴런 유형별 유전자 발현 보고서

작성일 2026-10-15. 목적은 PBM 다음 실험을 고르기 위해 어떤 뉴런 유형이 어떤 유전자를 더 발현하는지, 그리고 같은 뉴런 유형이 M1과 V1에서 어떻게 다른지를 정상 상태 공개 atlas로 정리하는 것이다. 이 보고서의 어떤 결과도 PBM 처리 효과를 측정한 것이 아니다.

그림은 `figures/report_m1v1/`에 PNG(180 dpi)와 PDF로 있다. 본문 그림은 Fig1~Fig7(Fig2는 a, b, c, Fig4는 a, b), 보조 그림은 Fig S1~S6이다.

## 1. 두 데이터셋과 각각이 답하는 질문

| | 10x RNA-seq (GSE185862, Allen CTX-HPF) | MERFISH (Zhuang-ABCA-1~4) |
|---|---|---|
| 유전자 | 전사체 전체 31,053개 | 1,122개 panel(세포 유형 구분용) |
| 세포 | Snap25-IRES2-Cre;Ai14 RFP+ 정렬 세포(뉴런 99%), 수컷 P57~P59. MOp 123,727개, VISp 29,763개 | 조직 절편의 모든 세포. 네 동물 MOp 98,334개, VISp 80,373개 |
| 동물 구조 | MOp donor 3마리와 VISp donor 3마리가 서로 다르고 library batch도 다름 | 한 동물의 같은 절편에 MOp와 VISp가 함께 있음 |
| 쓸 수 있는 비교 | 같은 영역 안에서 뉴런 유형 간 비교 | 같은 뉴런 유형 안에서 MOp 대 VISp 비교 |
| 쓸 수 없는 비교 | MOp 대 VISp (지역과 donor, batch가 완전히 겹침) | 미토콘드리아와 adenosine 경로 대부분 (panel에 OXPHOS 서브유닛 없음, adenosine은 Adora1과 Adora2a만 있음) |

따라서 질문 A(어떤 뉴런이 어떤 유전자를 더 발현하나)는 10x로, 질문 B(같은 뉴런이 M1과 V1에서 다른가)는 MERFISH로 답했다.

## 2. 방법 요약

### A. 10x 뉴런 유형 marker (`src/rnaseq/pseudobulk_10x.py`, `src/rnaseq/subclass_markers_10x.py`)

- GEO 원본 행렬 전체에서 영역×donor×subclass pseudobulk(원 count 합과 검출 세포 수)를 만들었다. 131개 group의 count 합이 저장된 전체 library 합과 일치함을 확인했다.
- 각 영역 안에서 subclass와 나머지 뉴런을 donor별로 짝지어 log2 CPM 비(fold change)를 구하고 donor 간 t 검정과 BH 보정을 했다.
- MOp marker 기준은 log2FC 1 이상, donor 3마리 모두 같은 방향, subclass 내 검출률 25% 이상, q 0.05 미만이다.
- 재현 기준은 VISp의 다른 donor 3마리(다른 batch)에서 log2FC 1 이상, 3마리 모두 같은 방향, 검출률 25% 이상이다. donor 3마리 t 검정은 자유도가 2라 재현 판정에는 q 기준을 쓰지 않았다.
- 주제 유전자 세트는 결과를 보기 전에 `src/rnaseq/focus_gene_sets.py`에 고정했다. OXPHOS 복합체 I~V는 공식 명명 규칙으로 정했고 핵 유전자와 mt 유전자를 나눴다. adenosine은 수용체, 생성, 제거, 수송, ATP 방출로 나눴다.
- 세트 점수는 해당 세트가 library에서 차지하는 비율을 같은 donor의 전체 뉴런 비율로 나눈 log2 값이다.

### B. MERFISH 동물 짝지은 지역 DEG (`src/merfish/zhuang_paired_deg.py`)

- 동물마다 subclass 안 MOp 평균 빼기 VISp 평균(공식 log2 단위)을 구하고 동물 간 one-sample t 검정(짝지은 검정)을 했다.
- 기준은 3마리 이상, 모든 동물 같은 부호, 차이 0.10 이상, subclass 내 BH q 0.10 미만이다.
- supertype 표준화 결과로 subclass 내부 구성 차이에 따른 것인지도 확인했다.

## 3. 결과

### 3.1 뉴런 유형 marker는 교과서적 패턴을 재현한다 (Fig1)

![Fig1](../../figures/report_m1v1/fig1_10x_subclass_markers_heatmap.png)

- MOp에서 marker 4,745개(subclass×gene)를 찾았고 subclass마다 50~87%가 별도 batch의 VISp donor에서 재현됐다. Sst Chodl은 VISp 세포가 34개뿐이라 재현 판정에서 제외했다.
- 상위 marker는 알려진 정체성 유전자와 일치한다. L4/5 IT의 Rorb, L5 PT의 Pou3f1, Npr3, Nefh, L6 CT의 Foxp2, Syt6, Pvalb의 Pvalb, Syt2, Vip의 Vip, Crh, Lamp5의 Ndnf 등이다. 분석 파이프라인의 타당성 근거로 사용한다.
- 전체 표는 `results/rnaseq/whole_transcriptome/subclass_markers/MOp_subclass_markers.csv`이다.

### 3.2 에너지 대사 유전자는 Pvalb와 L5 PT에 집중된다 (Fig3, Fig S1~S5)

![Fig3](../../figures/report_m1v1/fig3_10x_energy_module_heatmap.png)

세트 점수의 정의. 각 행은 기능 세트 하나이다. 세트에 속한 유전자 count의 합이 subclass library에서 차지하는 비율을, 같은 donor의 같은 영역 전체 뉴런에서의 비율로 나눈 log2 값을 donor 3마리에서 평균했다. 0은 그 영역 평균 뉴런과 같다는 뜻이다. OXPHOS 행은 핵 유전자 서브유닛만 쓰고, mtDNA가 만드는 13개 서브유닛은 mtDNA 복제 수와 세포 상태를 함께 반영하므로 한 행으로 따로 두었다. 각 영역은 자기 영역 뉴런을 기준으로 따로 계산했기 때문에 MOp와 VISp 패널은 같은 패턴의 독립 재현으로 읽고, 패널 사이 칸끼리의 크기 비교는 하지 않는다.

- OXPHOS 복합체 I~V가 가장 높은 유형은 Pvalb와 L5 PT이다. MOp에서 Pvalb +0.62~+0.77, L5 PT +0.47~+0.78이고, 다른 donor와 batch인 VISp에서도 Pvalb +0.69~+0.81, L5 PT +0.35~+0.53이다. 두 영역 세트 점수의 Pearson r은 0.92이다.
- 810 nm PBM의 1차 수용체로 거론되는 복합체 IV(cytochrome c oxidase)도 같은 순서이다.
- 미토콘드리아 생합성 세트는 Pvalb(MOp +1.56, VISp +1.28)와 L5 PT(+0.80, +0.64)에서 높다. Ppargc1a(PGC-1α)는 기존 PBM 22개 가설 panel에 포함된 유전자이다.
- L5 PT는 creatine 완충(MOp +1.11)과 해당 및 젖산 세트(+0.77)도 함께 높다. creatine 차이는 VISp에서 +0.39로 작아진다.
- Car3과 IT 계열은 OXPHOS가 평균 이하이다.

특정 유전자 하나에 끌려간 점수인지 확인 (Fig S1~S5, leave-one-out)

![Fig S3](../../figures/report_m1v1/figS3_energy_genes_complexIV_cycs.png)

- 보조 그림은 세트마다 굵은 세트 점수 행과 그 아래 구성 유전자의 값(유전자 CPM의 전체 뉴런 대비 log2 비)을 MOp와 VISp에 나란히 보이고, 오른쪽 막대와 숫자로 각 유전자가 세트 count에서 차지하는 비율을 보인다. 유전자 이름은 각 행 바로 왼쪽에 있고 세트 안에서 기여 비율 순으로 정렬했다.
- OXPHOS 복합체 I~V와 mtDNA 행에서 Pvalb와 L5 PT의 높은 점수는 대부분의 구성 유전자가 같은 방향으로 올라간 결과이다. 유전자를 하나씩 빼고 다시 계산해도 Pvalb는 복합체 I 0.69~0.73, 복합체 IV 0.57~0.82, 복합체 V 0.74~0.79 범위에 머물고 L5 PT도 0.45~0.83 범위에서 유지된다.
- 복합체 IV에서 Cox6a2는 Pvalb와 Lamp5에 특이적인 isoform이라 다른 유형에서는 강한 음의 값으로 보인다. 세트 count의 0.4%라 점수에는 거의 기여하지 않는다.
- 반대로 몇 개 유전자에 의존하는 세트가 있다. 생합성 세트는 Ppargc1a를 빼면 Pvalb가 MOp +1.56에서 +0.71로, L5 PT가 +0.80에서 +0.15로 떨어진다. 따라서 이 행은 사실상 Ppargc1a 한 유전자의 신호로 읽어야 한다. creatine 세트는 Ckb가 count의 76%를 차지해서 GABA 유형의 음의 값(Sst, Vip, Lamp5)은 Ckb를 빼면 0 근처로 바뀐다.
- 구성 유전자가 둘 이상인 261칸 중 유전자 하나를 빼서 부호가 바뀌는 칸은 16칸이다. 15칸이 생합성(Ppargc1a)과 creatine(Ckb), 해당(Pkm) 행이고, OXPHOS 행에서는 0 근처 값인 MOp L6b 복합체 II(+0.12, Sdha를 빼면 −0.04) 1칸뿐이다.
- Cox4i1과 Ndufa13은 이 GEO release에서 모든 세포의 count가 0이다. 생물학적 부재가 아니라 주석상 인접 유전자와 겹쳐 read가 집계되지 않은 것으로 보인다. 세트 count 기여가 0이라 점수에는 영향이 없다.
- 표는 `results/rnaseq/whole_transcriptome/energy_detail/`에 있다(`energy_set_summary.csv`에 세트 점수와 leave-one-out 범위, `energy_gene_contribution.csv`에 유전자별 기여 비율).

### 3.3 adenosine 경로의 뉴런 유형별 발현과 MOp 대 VISp 차이 (Fig2a~2c, Fig S6)

![Fig2a](../../figures/report_m1v1/fig2a_10x_adenosine_dotplot.png)

발현 수준 (Fig2a)

- Adora1은 Glut 뉴런에서 59~98% 세포에 검출되고 GABA 뉴런에서는 낮다(Lamp5, Vip, Sncg 검출 7~13%). Pvalb와 Sst는 중간이다(검출 41~49%).
- Glut 안에서 Adora1이 가장 높은 유형은 L6b(log2FC +1.13), L6 CT(+1.03), L5/6 NP(+0.96)이다. 가장 낮은 Glut는 L2/3 IT(−1.17)와 L6 IT(−0.69)이다.
- Adora2a, Adora2b, Adora3는 뉴런에서 거의 검출되지 않는다.

10x 지역 차이를 믿을 수 있는가 (Fig S6)

10x에서 MOp와 VISp는 donor와 library batch가 다르다. 그래서 MOp 빼기 VISp 값에는 지역 차이와 batch 차이가 섞인다. 이를 점검하기 위해 두 자료에 모두 있는 유전자(subclass마다 576~697개)에서 10x 지역 차이와 MERFISH 네 동물 짝지은 차이를 비교했다.

![Fig S6](../../figures/report_m1v1/figS6_10x_vs_merfish_region_concordance.png)

- Glut 7개 subclass에서는 MERFISH 차이가 0.10 이상인 유전자의 86~93%가 10x에서도 같은 부호였고, MERFISH 기준을 통과한 지역 DEG로 좁히면 92~100%였다. Glut 뉴런의 10x 지역 차이는 batch가 섞였어도 대부분 실제 지역 차이를 반영한다.
- GABA 4개 subclass는 일치도가 69~84%로 낮고 상관도 r 0.20~0.31로 약하다. MERFISH의 GABA 지역 차이 자체가 작아서 비교할 유전자도 적다(46~75개).

MOp 대 VISp 차이 (Fig2b)

![Fig2b](../../figures/report_m1v1/fig2b_adenosine_region_difference.png)

- 색은 subclass 안 MOp 빼기 VISp(10x, donor 평균 log2(CPM+1))이고 빨강이 MOp 쪽이다. 각 subclass 이름 아래 숫자는 위의 MERFISH 일치도이며, 85% 미만이거나 검증할 수 없는 subclass는 주황 글씨로 표시했다. 칼표는 MOp donor 3마리가 모두 VISp donor 3마리보다 높거나 낮은 칸이다.
- 검증 일치도가 높은 Glut IT 계열에서 두 가지 방향이 함께 보인다. MOp에서 Adora1은 낮고(L2/3 IT −0.47, L4/5 IT −0.32, L5 IT −0.34) adenosine 제거 효소 Adk는 높다(L2/3 IT +0.57, L4/5 IT +1.00, L5 IT +0.36). 모두 donor 범위가 겹치지 않는다. 생성 효소 Nt5e(CD73)도 Glut 다수에서 MOp가 높다(+0.59~+1.01). 다만 Adk와 Nt5e는 MERFISH panel에 없어 유전자 단위의 직접 검증은 없다.
- 수송체(Slc29a1, Slc29a2)와 ATP 방출(Panx1, Calhm2)은 subclass마다 방향이 다르고 크기가 작다.

Adora1을 두 자료에서 직접 비교 (Fig2c)

![Fig2c](../../figures/report_m1v1/fig2c_adora1_10x_vs_merfish.png)

- Glut에서는 두 자료가 같은 방향이다. L2/3 IT, L4/5 IT, L5 IT에서 MOp가 낮고 크기는 10x가 더 크다(두 자료는 단위가 달라 크기를 직접 비교하지 않는다).
- GABA에서는 10x가 MOp +0.7~+1.75로 크게 높지만 MERFISH 네 동물은 Pvalb +0.07, 나머지는 0 근처이다. 10x GABA의 Adora1 지역 차이는 batch 차이일 가능성이 높아 해석하지 않는다.

정리하면 MOp의 표층과 중층 IT 뉴런은 VISp보다 A1R 전사체가 적고 adenosine 제거 효소 전사체가 많다. 같은 adenosine 농도에서 MOp IT 뉴런이 덜 억제되는 쪽의 차이이다. 이 해석은 mRNA 수준이고 Adk의 지역 차이는 10x 단독 관찰이다.

### 3.4 PBM 공통 DEG의 기저 발현 분포 (Fig4a, Fig4b)

![Fig4a](../../figures/report_m1v1/fig4a_10x_pbm_candidates_MOp_VISp.png)

- 기존 PBM에서 해마와 피질 모두 같은 방향으로 변한 52개 중 10x에 있는 유전자를 뉴런 유형별 기저 발현으로 정렬했다. 행 순서는 PBM 방향과 MOp에서 발현이 가장 높은 유형 순이다.
- Fig4a는 MOp와 VISp를 나란히 놓았다. 행 z-score를 두 영역 열 전체로 계산해 두 패널의 같은 색이 같은 log2(CPM+1)를 뜻한다. 두 영역의 뉴런 유형별 패턴은 대체로 같다. 예를 들어 Arpp21은 두 영역 모두 IT에서 높고 Pvalb·Sst에서 낮으며, Nupr1과 Sorbs3는 두 영역 모두 Pvalb·Sst에서 높다.
- PBM으로 증가한 34개 중 8개가 MOp L2/3 IT에서 기저 발현이 가장 높았다. tPBM의 빛은 표층에 더 많이 도달하므로 표층 IT 뉴런이 증가 유전자를 주로 담고 있을 가능성이 있다. 다만 기존 PBM DEG는 glia가 섞인 조직 bulk에서 나왔고 이 atlas는 뉴런만 담고 있어 이 분포는 가설 생성용이다.

![Fig4b](../../figures/report_m1v1/fig4b_10x_pbm_candidates_region_difference.png)

- Fig4b는 같은 뉴런 유형 안의 MOp 빼기 VISp(donor 평균 log2(CPM+1))이다. Fig2b와 같은 규칙으로 subclass 이름 아래에 MERFISH 대조 부호 일치도를 붙였고, 85% 미만은 주황으로 표시했다. 칼표는 MOp donor 3마리가 모두 VISp donor 3마리 범위 밖에 있는 칸이며, 두 영역 모두 5 CPM 미만인 칸은 low로 가렸다.
- 신뢰도가 확인된 Glut 유형에서 칼표가 붙고 차이가 큰 예로는 MOp 쪽 Tnr(L4/5 IT +1.06, L5 PT +1.16), Pcdh19(L4/5 IT +1.59, L6 CT +1.09), Zeb2(L4/5 IT +1.30), Ccbe1(L6 IT +2.34), Sphkap(L6 CT +1.79), Trpc5(L6 IT +1.05)가 있다. VISp 쪽으로는 Stxbp6(L5/6 NP −2.39, L5 PT −1.06), Hsd11b1(L5 IT −1.77), Npnt(L2/3 IT −1.07)가 있다. 유형마다 방향이 바뀌는 유전자(Stxbp6, Sphkap, Hsd11b1 등)가 많아, 한 유전자의 지역 차이를 조직 전체 방향으로 일반화할 수 없다.
- 이 차이들은 PBM 반응이 아닌 정상 뇌의 기저 발현이며, 10x 지역 차이는 donor와 batch가 섞여 있다. MERFISH panel에 있는 후보(Penk, Trpc5, Zeb2 등)는 Zhuang 짝지은 DEG 표에서 따로 확인할 수 있다.
- 표는 `results/rnaseq/whole_transcriptome/subclass_markers/pbm_candidates_peak_subclass_MOp.csv`와 `results/rnaseq/whole_transcriptome/pbm_candidates_region/`에 있다.

### 3.5 M1과 V1의 같은 뉴런 유형 비교는 주로 영역 정체성 유전자 차이이다 (Fig5, Fig6)

![Fig5](../../figures/report_m1v1/fig5_merfish_paired_deg_counts.png)

![Fig6](../../figures/report_m1v1/fig6_merfish_paired_deg_top_heatmap.png)

- 네 동물 짝지은 기준을 통과한 subclass×gene은 496개(유전자 260개)이다. 거의 모두 Glut 뉴런에서 나왔고 L5 ET가 가장 많다(122개). GABA 뉴런은 Pvalb와 Sst에서 1개씩뿐이다.
- 여러 subclass에서 반복된 유전자는 영역 정체성과 투사 관련 유전자이다. MOp 쪽은 Adgrl2(7개 subclass), Scn4b(6개), Nefh, Ppp1r1b, Epha4, Epha7이고 VISp 쪽은 Nr2f1(COUP-TFI, 5개), Adcy8(5개), Satb1, Chn2이다. Nr2f1은 후방 피질 정체성 인자로 알려져 있어 방법의 타당성 근거가 된다.
- 표는 `results/zhuang/paired_deg/`에 있다.

### 3.6 Adora1과 Penk의 지역 차이 (Fig7)

![Fig7](../../figures/report_m1v1/fig7_merfish_adora1_penk_animals.png)

- Adora1은 L2/3 IT에서 MOp가 낮다(네 동물 평균 −0.16, q 0.098, supertype 표준화 후 −0.16). L4/5 IT도 네 동물 모두 음이지만(−0.22) q가 0.16이다. L6 CT의 음의 차이는 이전 점검에서 대부분 supertype 구성 효과였다.
- Penk(PBM 공통 DEG, 감소 방향)은 L2/3 IT에서 MOp가 높다(+0.30, q 0.054, 표준화 후 +0.31). L6 IT도 +0.28이지만 q 0.19이다.
- 두 차이 모두 log2 단위 0.2~0.3 이하로 작다.

### 3.7 MOp와 VISp의 뉴런 구성 차이 (MERFISH 네 동물)

| subclass | MOp 뉴런 중 비율 | VISp 뉴런 중 비율 |
|---|---|---|
| Pvalb | 5.2~6.4% | 3.7~4.3% |
| L5 ET | 4.0~4.2% | 3.7~5.2% |
| L4/5 IT | 11.9~18.0% | 26.4~28.9% |
| L5 NP | 2.4~2.8% | 1.5~2.0% |

에너지 대사가 가장 높은 두 유형 가운데 Pvalb는 네 동물 모두 MOp에 더 많고 L5 ET는 차이가 없다.

## 4. 종합 해석

1. 뉴런 고유 수준에서 PBM의 1차 표적(복합체 IV)과 미토콘드리아 생합성 능력이 가장 큰 유형은 Pvalb와 L5 PT(ET)이다.
2. Adora1은 심층 Glut(L6 CT, L6b, L5/6 NP)에 가장 많고 표층 L2/3 IT와 GABA 뉴런에 적다. L6 CT는 시상으로 투사하므로 피질 adenosine이 시상피질 수면 진동에 관여하는 경로 후보가 된다.
3. 같은 뉴런 유형의 M1 대 V1 차이는 주로 영역 정체성 유전자이다. adenosine 쪽에서는 MOp IT 뉴런이 VISp보다 Adora1이 적고 Adk가 많은 방향이 보인다. Adora1은 두 자료에서 확인됐고 Adk는 MERFISH 검증 일치도가 높은 Glut에서의 10x 단독 관찰이다. 크기는 log2 단위 0.3~1.0 수준이다. M1과 V1의 PBM 반응이 다르다면 이 adenosine 감수성 차이와 함께 구성 차이(MOp의 Pvalb 비율이 높고 L4/5 IT 비율이 낮음), 깊이에 따른 빛 도달, glia와 혈관 같은 비뉴런 요인을 같이 고려해야 한다.
4. 미토콘드리아와 adenosine 경로 대부분은 MERFISH panel에 없다. 이 유전자들의 지역 차이는 10x로만 볼 수 있고, 10x 지역 차이는 Glut에서만 MERFISH로 신뢰도가 확인됐다. GABA 뉴런의 지역 차이는 공개 atlas로 판단할 수 없어 실험으로 확인해야 한다.

## 5. 다음 실험 후보

| 우선순위 | 질문 | 실험 | 이 보고서의 근거 |
|---|---|---|---|
| 1 | PBM에 먼저 반응하는 세포는 Pvalb와 L5 PT인가 | tPBM 후 M1에서 c-Fos 또는 즉시 초기 유전자와 세포 유형 marker의 동시 FISH(RNAscope). Pvalb, Pou3f1 또는 Fezf2(L5 PT), Cux2(L2/3 IT), Foxp2(L6 CT) | Fig3에서 OXPHOS와 Ppargc1a가 가장 높은 유형. Fig1 marker가 probe 선택 근거 |
| 2 | PBM의 수면 효과는 A1R을 거치는가 | tPBM 전 A1R 길항제(DPCPX) 투여 후 NREM과 SWA 변화가 사라지는지 확인 | Adora1만 뉴런에 높게 발현되고 Adora2a, 2b, 3은 거의 없음(Fig2a) |
| 3 | 심층 L6 CT의 A1R이 경로인가 | 1번 FISH에 Adora1 probe를 더해 PBM 후 L6 CT 반응 확인. 이후 L6 CT 선택적 조작(Ntsr1-Cre 등) 고려 | Adora1이 L6 CT와 L6b에서 가장 높음(Fig2a) |
| 4 | M1 대 V1 차이가 Pvalb 구성 때문인가 | 같은 자극 조건에서 M1과 V1의 Pvalb c-Fos 비율 비교 | Pvalb 비율이 MOp에서 네 동물 모두 높음(3.7) |
| 5 | 표층 IT가 PBM 유전자 변화를 주도하는가 | 기존 PBM 조직의 표층과 심층 분리 qPCR 또는 공간 전사체 | PBM 증가 유전자가 L2/3 IT 기저 발현에 몰림(Fig4) |

1번이 가장 정보량이 크다. 1번 결과가 어떤 세포 유형을 가리키느냐에 따라 2번과 3번의 조작 대상이 정해진다.

## 6. 한계

- 10x는 Snap25-Cre;Ai14 reporter 수컷의 정렬 뉴런이며 glia와 혈관 세포가 없다. donor 수가 영역당 3마리라 t 검정 검정력이 낮다.
- MERFISH는 B6N 네 동물이고 panel이 세포 유형 구분용이다. sagittal 동물(ABCA-3, -4)의 차이 크기가 coronal 동물보다 2~3배 크다.
- 두 자료 모두 mRNA이다. 복합체 IV 효소 활성이나 A1R 단백질 양과 다를 수 있다.
- 어떤 결과도 PBM 처리 조건을 포함하지 않는다.

## 7. 재현

```bash
.venv/bin/python src/rnaseq/download_geo.py GSE185862_expression_matrix_10x.hdf5
.venv/bin/python src/rnaseq/pseudobulk_10x.py --workers 2
.venv/bin/python src/rnaseq/subclass_markers_10x.py
.venv/bin/python src/merfish/zhuang_paired_deg.py
.venv/bin/python src/report_figures_m1v1.py
.venv/bin/python src/rnaseq/energy_gene_detail.py
.venv/bin/python src/rnaseq/adenosine_region_compare.py
.venv/bin/python src/rnaseq/pbm_candidates_region.py
```

`pbm_candidates_region.py`는 Fig4a와 Fig4b를 만든다.

`report_figures_m1v1.py`는 Fig1, Fig2a, Fig5~Fig7(그리고 이전 MOp 단독 Fig4)을, `energy_gene_detail.py`는 Fig3과 Fig S1~S5를, `adenosine_region_compare.py`는 Fig2b, Fig2c, Fig S6을 만든다.

MERFISH 입력은 `src/merfish/run_zhuang.py`와 `src/merfish/zhuang_supertype_check.py` 실행 결과를 사용한다([Zhuang 보고서](../zhuang/zhuang_report_ko.md)).
