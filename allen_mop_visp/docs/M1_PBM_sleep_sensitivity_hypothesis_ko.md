# M1의 PBM 수면 반응이 V1보다 클 수 있다는 가설

2026-10-12. 현재 Allen 구성·발현·투사 분석과 확인한 원 논문을 연결한 **검증할 가설**이다. M1–V1 PBM 처리 실험이나 새로운 유의성 분석을 수행한 결과가 아니다. 기존 그림은 변경하지 않았다.

## 제안하는 가설 문장

> **810 nm 연속광 PBM을 180분 동안 적용할 때, M1의 sham 대비 NREM 증가량은 V1보다 클 것이다. 이 지역 차이는 세포 유형별 adenosine 반응과 수면 관련 하위 회로로의 신호 전달 차이에서 발생할 수 있으며, MOp→ZI 회로를 우선 검증할 후보로 둔다.**

이 가설의 1차 검증 대상은 **특정 자극 조건에서의 반응 크기**다. ‘M1이 더 민감하다’를 **더 낮은 광량에서 같은 수면 효과가 나타난다**는 의미로 사용하려면 여러 광량의 용량–반응 비교가 추가로 필요하다. 현재 atlas 또는 한 자극 조건만으로 수면 유도 역치가 낮다고 결론내리지 않는다.

논문 서술용으로는 다음과 같이 쓸 수 있다.

> PBM은 신경세포의 에너지 대사와 세포외 adenosine 신호를 변화시켜 수면에 영향을 줄 가능성이 있다. 본 분석에서 MOp와 VISp는 세포 구성, 일부 GABAergic subclass의 Adora1 기저 발현, ZI로의 해부학적 투사 신호에서 차이를 보였다. 이에 동일한 유효 광량과 생리적 조건에서 M1 자극이 V1 자극보다 큰 NREM 반응을 유도할 수 있다는 가설을 제안한다. 특히 국소 adenosine 변화가 세포 유형별 반응을 거쳐 수면 관련 ZI 집단에 전달되는 과정이 지역 차이를 매개할 수 있다. 이 연결은 현재의 해부학적·전사체 관찰에서 도출한 중개가설이며, 표적 세포와 기능 방향의 검증이 필요하다.

## 실제 자극 조건과 용량의 범위

사용자가 확인한 값은 **810 nm, 연속광, fiber 끝 기준 약 20 mW/cm², 180분, core 직경 400 µm**다. 일정한 조사세기가 180분 내내 유지된다고 가정한 nominal fiber-tip fluence는 다음과 같다.

`0.020 W/cm² × 10,800 s = 216 J/cm²`

이는 **fiber 끝에서의 명목값**이며 뇌 조직이 흡수한 광량이나 각 층의 도달광량은 아니다. 현재 두 영역의 실제 조사 부피·깊이별 광량·온도 및 CCF 등록은 검증하지 않았다. 동일 fiber-tip 값만으로 두 영역의 유효 조직 선량이 같다고 가정하지 않는다. 한 번의 자극은 3시간이므로 ZT·자극 시작 상태·기저 수면압력 및 자극 중/후 분석 시간을 맞춰야 한다.

M1 좌표는 bregma AP +1.34, ML 2.46 mm, 뇌 표면 DV 0이며 V1 AP 약 −2.3 mm는 아직 잠정값이다. 전체 MOp/VISp atlas의 값은 해당 광섬유 주변 ROI의 값이 아니다. 조건 정본은 [사용자 자극 기록](stimulation_coordinates_user.json)에 있다.

## 현재 분석에서 가설을 제안하게 하는 관측

모든 숫자 쌍은 **MOp / VISp** 순서다. 세포 구성은 각 영역의 전체 QC 세포를 분모로 한다. 발현 평균은 log1p(CP10k), 공유 cluster 차이는 MOp−VISp 평균 log1p(CP10k)이며 log2FC가 아니다.

| 관측 근거 | 실제 결과 | 가설에 기여하는 범위 |
|---|---|---|
| ZI로의 regional axonal signal | WT MOp 8개·VISp 33개 실험의 양측 NPV 중앙값 **0.17054 / 0.02560** | MOp→ZI를 후속 tracing·생리 후보로 정할 근거 |
| ZI 주입 조건 민감도 | P56 및 source signal fraction ≥0.8에서 n=4/10, NPV **0.26841 / 0.04803**; 다섯 선택 조건에서 MOp가 큰 방향 유지 | 다른 주입 선택에서도 남는 해부학적 관측. 지역의 생리적 효과 크기 추정은 아님 |
| 전체 세포 구성 | Astro **10.86% / 8.58%**, GABA **9.21% / 8.07%** | 국소 대사와 E/I 회로를 세포 유형별로 조사할 근거 |
| GABAergic Adora1 | P58 Pvalb 평균 log1p(CP10k) **0.3545 / 0.2067**, 검출률 **45.81% / 35.31%**, n=1,576/422; 공유 cluster 차이 +0.1392 | 일부 interneuron의 A1 반응성이 다를 가능성 |
| 흥분성 Adora1의 반대 방향 | L2/3 IT 평균 **0.4410 / 0.5779**, L4/5 IT **0.8984 / 1.0006** | 수용체 차이가 모든 뉴런에서 MOp 우세라는 설명을 배제함 |

앞의 세포 구성은 **WT B6J 수컷 한 동물**의 분모이며 세포 밀도·기능 활성·독립 동물 반복이 아니다. RNA-seq은 WT가 아닌 male Snap25-IRES2-Cre;Ai14 reporter 자료이고, P58 기록 donor는 MOp 1개·VISp 2개다. 지역과 library-preparation date가 분리되어 있다. 공유 cluster 비교도 이 교란을 제거하지 않으므로 생물학적 지역 효과의 P value/FDR를 제시하지 않는다.

NPV는 `표적 projection signal volume / 전체 injection signal volume`이다. 연결 확률·시냅스 수·Lhx6 세포에 대한 입력·PBM 전달 효율이 아니다. WT 내에서도 주입 spread·viral construct·관찰 기간 등의 차이가 남는다. Cre 실험은 별도 요약했다. 따라서 현재 가설의 우선 회로 후보는 ZI이지만, 이 표는 **MOp가 ZI의 수면촉진 세포를 더 강하게 구동한다는 증거가 아니다**.

### 2026-10-16 추가 관측

이후 B6N Zhuang MERFISH 네 동물과 10x 전사체 전체를 분석했다. [뉴런 유형별 발현 보고서](../results/neuron_type_report/M1_V1_neuron_type_expression_report_ko.md), [Zhuang 보고서](../results/zhuang/zhuang_report_ko.md)에 근거가 있다.

| 관측 근거 | 실제 결과 | 가설에 기여하는 범위 |
|---|---|---|
| 흥분성 IT의 Adora1 | MERFISH 네 동물에서 MOp가 낮음(L2/3 IT −0.16, L4/5 IT −0.22, supertype 표준화 후 유지). 10x도 같은 방향 | 같은 adenosine에서 MOp IT 뉴런이 덜 억제될 가능성 |
| 흥분성 IT의 Adk | 10x 전사체 전체에서 MOp가 높음(L4/5 IT log2(CPM+1) +1.00, L2/3 IT +0.57). MERFISH panel에는 없음 | MOp IT 주변의 adenosine 제거가 더 클 가능성. 순축적이 MOp에서 크다는 가정과 반대 방향일 수 있음 |
| GABAergic Adora1 | MERFISH에서 Pvalb +0.07만 작게 남고 Sst, Vip, Lamp5는 0 근처. 10x의 큰 MOp 우세는 MERFISH와 불일치 | 위 표의 interneuron A1 차이 근거를 약화함 |
| 에너지 대사 | OXPHOS 복합체 I~V와 Ppargc1a가 Pvalb와 L5 PT에서 가장 높음. Pvalb는 뉴런 중 비율이 MOp에서 높음(5.2~6.4% 대 3.7~4.3%) | PBM의 1차 표적(복합체 IV)이 많은 세포가 MOp에 상대적으로 많을 가능성 |

이 관측은 모두 정상 뇌의 mRNA이며 PBM 처리 효과나 adenosine 농도를 측정한 것이 아니다.

## 가장 중요한 반대 근거와 미확인 연결

1. **Astro가 많다고 adenosine이 더 축적되는 것은 아니다.** 생성뿐 아니라 수송·대사·제거도 영향을 준다. 표면 L1의 Astro 비율은 MOp 37.05%·VISp 37.14%로 거의 같고, L5·L6b GABA 비율은 전체 영역과 달리 VISp가 높다. 현재 표면 자극에 전체 영역 비율을 그대로 적용할 수 없다.
2. **생성·수송 유전자가 일관되게 MOp 우세가 아니다.** Nt5e/Pvalb는 VISp가 높고, Nt5e·Entpd1·Slc29a1의 여러 subclass는 cluster 가중치에 따라 방향이 바뀐다. Adk/L4/5 IT는 MOp가 높지만 ADK는 adenosine 대사·제거와 연결되는 후보이므로 순축적 우세를 의미하지 않는다. 이 Adk 차이는 이후 10x 전사체 전체에서도 확인됐다(+1.00). B6J MERFISH에서는 adenosine 후보 8개가 모두 미측정이고, B6N MERFISH에서는 Adora1과 Adora2a만 측정된다.
3. **Astro RNA 결과는 희소 표집에 민감하다.** P58 Astro n=125/59이며 공유 cluster는 n=123/57과 n=2/2 두 개뿐이다. 동일 cluster 가중치에서는 n=2/2 집단도 절반을 차지하므로 큰 Nt5e·Entpd1·Adk 차이를 안정한 지역 기전으로 제시할 수 없다.
4. **A1의 기능 방향은 세포 유형에 따라 다를 수 있다.** Interneuron의 A1 신호가 억제를 유발하면 pyramidal neuron의 disinhibition으로 이어질 수도 있다. 다른 세포의 A1은 흥분성 출력 감소와 연결될 수 있다. `GABA Adora1↑ → 피질 억제↑ → NREM↑`라는 일률적 화살표를 사용하지 않는다.
5. **MOp→ZI의 postsynaptic identity와 부호가 미확인이다.** Tracer는 Lhx6 세포의 시냅스를 분리하지 않는다. PBM·adenosine이 MOp 출력을 억제한다면, 더 많은 ZI 축삭이 반드시 수면촉진 세포를 더 활성화하는 것도 아니다. ZI의 다른 집단, 국소 억제와 disinhibition, REM 중심 효과가 경쟁 설명이다. Lhx6 문헌도 양측 chemogenetic/ZT 조건과 단측 optogenetic 조건에서 수면 효과가 다르다 [5,6].
6. **다른 수면 표적으로 일반화할 수 없다.** VLPO·MEPO 신호는 작고 선택 조건에 민감하며, 큰 신호의 실험에는 주변 영역 주입이 있다. PZ는 미평가다. 모든 NREM 관련 표적에 MOp 투사가 더 크다고 할 수 없다.
7. **기존 PBM DEG는 현재 처리 비교가 아니다.** Nt5e와 Adora2a 감소는 1064 nm hypobaric-hypoxia PFC 자료의 결과다. 정상 WT MOp–VISp, 현재 810 nm 연속광 또는 급성 extracellular adenosine 반응과 동일시하지 않는다.

## 원 논문과 직접 지지하는 주장

| 문헌 | 확인한 결과 | 이 가설에 적용할 때의 한계·확인 범위 |
|---|---|---|
| **Jung et al., 2026 [1]** | Mouse에서 sham 대비 PBM 처리 시 NREM·NREM delta 증가, frontal cortex adenosine 증가; 배양세포 ATP 증가 | 직접 PBM–수면/adenosine 근거. 공식 초록·공개 Figure 3 caption·공개 보충자료 확인. 보충자료의 배양 조건에 810 nm·3시간 명시. 구독 본문 전체를 확인하지 못했으며 MOp–VISp 우열, in vivo 연속광 여부·20 mW/cm² 동등성, adenosine/ENT 차단의 매개 인과성은 확정하지 않음 |
| **Sharma et al., 2011 [2]** | 배양 mouse cortical neuron에 810 nm 연속광·25 mW/cm²; ATP는 3 J/cm²에서 최고이고 10/30 J/cm²에서는 control 수준으로 회귀 | 전문 확인. ATP·미토콘드리아 반응과 용량 의존성을 지지. 수면·세포외 adenosine·지역 차이를 측정하지 않음. 배양세포 fluence를 현재 216 J/cm² fiber-tip 값에 외삽하지 않음 |
| **Porkka-Heiskanen et al., 1997 [3]** | Cat basal forebrain에서 prolonged waking의 extracellular adenosine 증가와 회복수면의 감소; 국소 transport inhibitor로 수면박탈 후 양상 재현 | 공식 초록 확인. Adenosine–sleep pressure의 관찰·조작 근거이며 MOp–VISp 비교가 아님 |
| **Benington et al., 1995 [4]** | Rat A1 agonist로 NREM slow-wave activity의 용량 의존 증가 | 공식 초록 확인. A1과 수면압력의 EEG 특징을 연결하지만 SWA 증가와 NREM 시간 증가를 같은 endpoint로 쓰지 않음 |
| **Liu et al., 2017 [5]** | ZI Lhx6-positive GABA 세포의 수면 조절; NREM 효과는 ZT 의존적이며 REM도 변화 | 전문 및 기존 표적 문서 확인. 해당 세포를 실제 MOp input이 조절하는지, PBM이 같은 효과를 내는지 미확인 |
| **Vidal-Ortiz et al., 2024 (online 2023) [6]** | Unilateral ZI Lhx6 optogenetic stimulation은 REM 증가, NREM/wake 양 불변 | 전문 및 기존 표적 문서 확인. ZI/Lhx6을 보편적 NREM 촉진 경로로 부르는 것을 제한하는 근거 |
| **Zhou et al., 2024 [7]** | 40 Hz 백색 시각 flicker에서 V1 adenosine·이후 수면 증가; ENT2/V1 neuron 조작으로 관련 경로 검증 | 전문 확인. **810 nm 연속광 PBM와 다른 자극**이며 M1 비교 없음. V1도 광자극을 adenosine·수면 반응으로 변환할 수 있다는 근거. V1 A1 매개성은 후속 검증 사항이며, CD73 KO에서도 adenosine 증가가 유지됨 |

Jung 공개 보충자료의 배양 조사강도는 brain depth의 추정 도달량에 대응하는 값이다. Fiber-tip irradiance와 같은 값이 아니다. 파장·3시간이 유사하더라도 현재 프로토콜과 생체 유효 선량의 동등성은 미확인이다. 전문에서 검증되지 않은 receptor blockade나 ENT2 의존성을 확인된 PBM 결과로 쓰지 않는다.

Zhou의 ENT2 결과는 **세포내 대사에서 생성한 adenosine의 수송** 경로를 지지한다. 이를 `PBM → extracellular ATP 방출 → CD39/CD73 분해`의 증거로 쓰지 않는다. 현재 PBM에서 ENT2 경로인지 extracellular ATP 분해인지, 두 경로가 함께 작동하는지는 별도로 검증할 중개가설이다. 단순 CCO/NO 설명도 현재 분석에서 직접 검증하지 않았다.

이번 원문·메타데이터 확인 범위에서는 **현재 조건의 PBM을 MOp와 VISp에 직접 비교해 MOp 우세를 확정한 문헌은 확인하지 못했다**. 이는 문헌의 완전한 부재를 증명한다는 뜻은 아니다.

## 가설을 구분할 검증과 예측

| 경쟁 설명 | 판별할 관측 | 지지·약화 기준 |
|---|---|---|
| M1의 adenosine 순축적이 더 큼 | Region×sham/PBM에서 상승·AUC·감쇠 비교 | M1 상승/AUC가 크면 생성·방출/제거 모델을 지지. 비슷하면 생성량 우세보다 반응·회로 모델에 무게 |
| 같은 adenosine 입력에 M1 회로가 더 크게 반응 | 같은 국소 입력에서 Pvalb/Sst·IT의 생리, LFP와 EEG/EMG 반응 비교 | PBM 없이도 같은 입력에 지역 차이가 남으면 수용체·회로 반응 모델을 지지 |
| MOp→ZI 전달이 M1의 수면 이득을 매개 | Target Lhx6 identity·시냅스/생리 확인 후 projection별 개입 | 국소 adenosine 반응은 남고 sham 대비 M1의 추가 NREM 효과가 줄면 전달 모델을 지지. 개입 자체의 baseline 수면 효과와 구분 |
| 광전달·층·혈관/온도의 차이가 원인 | 실제 ROI·깊이별 dose·temperature 등을 맞춘 비교 | 유효 광량을 맞추자 지역 차이가 줄면 명목 입력보다 광전달이 중요한 설명 |
| M1 motor 억제가 수면처럼 보임 | EEG/EMG로 NREM·REM·quiet wake 분리 | 움직임/EMG 감소만 있고 EEG 수면이 늘지 않으면 수면촉진 가설을 지지하지 않음 |

GRAB-Ado를 사용한다면 서로 다른 영역의 dF/F를 같은 절대 adenosine 농도로 간주하지 않는다. Sensor expression·광학 수집·artifact/control과 각 영역의 sham 대비 변화를 확인한다. ATP 생성, extracellular adenosine, 단백질·전류, 상태 분류를 서로 다른 측정으로 유지한다.

주 endpoint는 사전에 정한 창의 **sham 대비 NREM 시간 변화**로 두고, 수면 시작 latency·NREM delta power·REM·wake를 별도로 보고한다. 예를 들어 자극 중 180분과 자극 후의 동일하게 고정한 창을 분리한다. SWA 또는 국소 slow wave만으로 전체 NREM 시간 증가를 판정하지 않는다.

독립 WT 동물 반복, 연령·성별·ZT·기저 수면/각성·조건 및 prep/실험 순서를 균형 있게 배치한다. 같은 동물에서 여러 조건을 사용하면 순서·이월 효과와 동물 짝지음을 보존한다. **Region×PBM 상호작용** 또는 sham 대비 변화량의 직접 차이를 평가한다. ‘M1만 유의하고 V1은 비유의’라는 비교는 지역 차이의 검정이 아니다. 민감도/역치 주장에는 여러 선량의 실제 용량–반응 자료가 필요하고, PBM의 이상성 반응을 고려해 단순 단조 증가 곡선을 가정하지 않는다.

현재 근거를 바탕으로 검증할 가설은 **‘M1의 adenosine 생성량 우세’ 하나로 고정하기보다, 국소 adenosine dynamics와 세포 유형별 반응 및 MOp→ZI 전달을 구분한 지역 반응 가설**이다.

## 참고 문헌

1. Jung J et al. (2026). *Transcranial photobiomodulation to enhance sleep by elevating adenosine signaling: a translational study.* Molecular Psychiatry. DOI [10.1038/s41380-026-03875-y](https://doi.org/10.1038/s41380-026-03875-y); PMID [42744976](https://pubmed.ncbi.nlm.nih.gov/42744976/). [공개 보충자료](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41380-026-03875-y/MediaObjects/41380_2026_3875_MOESM1_ESM.docx), [공개 Figure 3](https://www.nature.com/articles/s41380-026-03875-y/figures/3).
2. Sharma SK et al. (2011). *Dose response effects of 810 nm laser light on mouse primary cortical neurons.* Lasers in Surgery and Medicine. DOI [10.1002/lsm.21100](https://doi.org/10.1002/lsm.21100); [PMC3199299](https://pmc.ncbi.nlm.nih.gov/articles/PMC3199299/).
3. Porkka-Heiskanen T et al. (1997). *Adenosine: a mediator of the sleep-inducing effects of prolonged wakefulness.* Science. DOI [10.1126/science.276.5316.1265](https://doi.org/10.1126/science.276.5316.1265); PMID [9157887](https://pubmed.ncbi.nlm.nih.gov/9157887/).
4. Benington JH et al. (1995). *Stimulation of A1 adenosine receptors mimics the electroencephalographic effects of sleep deprivation.* Brain Research. DOI [10.1016/0006-8993(95)00590-M](https://doi.org/10.1016/0006-8993(95)00590-M); PMID [8548323](https://pubmed.ncbi.nlm.nih.gov/8548323/).
5. Liu K et al. (2017). *Lhx6-positive GABA-releasing neurons of the zona incerta promote sleep.* Nature. DOI [10.1038/nature23663](https://doi.org/10.1038/nature23663); [PMC5958617](https://pmc.ncbi.nlm.nih.gov/articles/PMC5958617/).
6. Vidal-Ortiz A et al. (2024; online 2023). *Unilateral optogenetic stimulation of Lhx6 neurons in the zona incerta increases REM sleep.* Sleep. DOI [10.1093/sleep/zsad217](https://doi.org/10.1093/sleep/zsad217); [PMC11502959](https://pmc.ncbi.nlm.nih.gov/articles/PMC11502959/).
7. Zhou X et al. (2024). *40 Hz light flickering promotes sleep through cortical adenosine signaling.* Cell Research. DOI [10.1038/s41422-023-00920-1](https://doi.org/10.1038/s41422-023-00920-1); [PMC10907382](https://pmc.ncbi.nlm.nih.gov/articles/PMC10907382/).

서지·초록 정본은 [공식 EuropePMC 응답과 검증 registry](../input/literature/PBM_sleep_hypothesis_reference_registry.json)에 보존했다. 원 논문의 전체 접근 가능 범위는 위 표에 구분했다.

## 분석 수치의 출처

- [세포 구성 보고서](../results/merfish/composition_report_ko.md), [영역 비율](../results/merfish/composition_region_major.csv), [층별 비율](../results/merfish/composition_layer_major.csv)
- [P58·cluster 민감도](../results/rnaseq/panel22_P58_cohort_cluster_sensitivity.csv), [RNA-seq 보고서](../results/rnaseq/RNAseq_results_ko.md)
- [Astro cluster별 원수치](../results/rnaseq/region_cluster_gene_summary.csv), [공유 cluster 차이](../results/rnaseq/shared_cluster_region_difference.csv)
- [실제 투사 보고서](../results/connectivity/NREM_projection_report_ko.md), [WT 주입 조건별 요약](../results/connectivity/WT_target_projection_descriptive_summary.csv), [NREM 표적 문헌](NREM_projection_target_definition_ko.md)
- [기존 PBM 연결 보고서](../results/pbm/PBM_bridge_report_ko.md)
- [L4·그림·분모 표현 검토](composition_visualization_review_ko.md)
