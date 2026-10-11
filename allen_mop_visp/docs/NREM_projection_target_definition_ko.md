# NREM 관련 projection target의 기능·해부학적 정의

**MOp·VISp에서 NREM 관련 영역으로 투사하는지**를 먼저 평가한다. 아래 정의는 기능 문헌과 Allen Structure graph 1을 연결한 표적 목록이며, 투사가 수면을 유발한다는 결과가 아니다. 기계 판독용 정본은 `input/connectivity/nrem_target_definitions.json`과 같은 내용의 CSV다.

2026-10-11에 공식 Allen ontology 1,327개 항목의 ID·acronym·명칭을 확인하고, 아래 원 논문 10편의 PMC 원문과 DOI·PMID 메타데이터를 직접 검증했다. 해부학적 label의 정확성과 기능 세포군의 대응 정확성을 별도 필드로 저장했다. MA와 논문의 MCPO 사이 경계 대응은 **provisional**이며, PZ는 **미대응**이다.

## 분석에 사용할 구분

| 문헌 표적 | Allen exact label / ID | 역할과 근거 | 현재 projection으로 확인할 수 없는 것 |
|---|---|---|---|
| VLPO | VLPO — Ventrolateral preoptic nucleus / **689** | **Core:** galanin/GABA 세포군의 인과적 NREM 촉진 [1] | Gal/GABA 세포에 직접 접촉하는지, synapse·효능 |
| MnPO / MnPN | **MEPO** — Median preoptic nucleus / **452** | **Core:** Vgat 세포의 인과적 NREM 촉진 및 sleep homeostasis 관련 [2,3] | Vgat와 Vglut2 세포를 구분한 postsynaptic targeting |
| Ventral ZI, Lhx6 | ZI — Zona incerta / **797** | **세포 유형·조작 조건 특이:** NREM과 REM 모두 관련 [4,5] | Lhx6 세포에 도달하는지, ventral ZI 세부 위치, NREM 특이성 |
| Caudal BF, SOM/Sst | SI — Substantia innominata / **342** | **세포 유형 특이·혼합 근거:** NREM 촉진 연구와 급성 효과가 없는 연구가 함께 있음 [6,7] | SOM 세포에 도달하는지; SI 전체의 기능 |
| BF의 magnocellular preoptic 맥락 | MA — Magnocellular nucleus / **298** | **선택적 해부학 맥락:** 논문 MCPO와의 경계 대응 미검증 [6,7] | MA=논문 MCPO라는 정확한 대응, MA 특이 NREM 효과 |
| BF의 horizontal diagonal band 맥락 | NDB — Diagonal band nucleus / **596** | **선택적 해부학 맥락:** 논문의 HDB보다 atlas label이 넓음 [6,7] | horizontal/vertical limb 및 SOM 세포 targeting |
| Thalamic reticular nucleus | **RT** — Reticular nucleus of the thalamus / **262** | **NREM 관련 추가 표적:** local spindle·delta 등 수면 리듬 조절 [8] | 전체 NREM 양 증가·수면 개시를 유발한다는 결론 |
| Parafacial zone, PZ | **exact label 없음 / ID 없음** | **기능상 core, atlas 미평가:** inhibitory/GABA 세포군의 SWS/NREM·SWA 촉진 [9,10] | 정규 unionize로 PZ projection을 계산하는 것 |

**Allen `TRN` / 574는 Tegmental reticular nucleus라는 뇌간 영역이다.** 문헌에서 약칭 TRN으로 부르는 시상 망상핵은 Allen **RT / 262**와 대응하며, 574는 분석 표적에서 제외한다. PZ를 PGRN이나 인접 reticular 영역으로 치환하지 않는다. PZ의 추후 분석은 원문 조직학적 위치를 CCF에 등록한 별도 voxel ROI가 있어야 가능하다.

Core 두 영역, ZI·BF, RT, 선택적 BF 해부학 맥락은 각각의 tier로 보고한다. 이들의 신호를 하나의 ‘NREM-promoting projection score’로 합하지 않는다. MA·NDB는 unionize를 조회할 수 있는 실제 ontology label이지만, 기능 표적 대응이 확정된 core는 아니다.

## NREM·REM·각성의 근거를 구분한다

- **VLPO Gal:** 저빈도 optogenetic activation은 NREM을 늘리고 photoinhibition은 줄였다. 그 조작에서는 REM 양이 변하지 않았지만, chemogenetic activation은 NREM 증가와 REM 감소·저체온을 함께 유발했다. 고빈도 광자극의 depolarization block·각성 효과도 있어 ‘모든 활성화는 NREM만 증가’로 일반화하지 않는다 [1].
- **MnPO Vgat:** activation은 NREM latency를 줄이고 NREM을 늘리며 REM을 줄였다. ablation은 낮 NREM을 감소시키고 정상 수면량에 도달하는 sleep pressure 요구를 바꿨다. 수면박탈 후 회복량의 군 간 차이는 유의하지 않았으므로 ‘homeostatic rebound 소실’로 해석하지 않는다 [2]. Rat MnPN GABA/Fos와 sleep pressure의 관계는 활동 상관 근거이고, 그 자체로 인과성이 아니다 [3].
- **ZI Lhx6:** 양측 chemogenetic 조작에서 NREM 효과는 ZT에 의존했고 REM도 변화했다 [4]. 별도 단측 저빈도 optogenetic 연구는 REM을 늘렸으며 NREM과 wake 양은 바뀌지 않았다 [5]. 전체 ZI와 모든 Lhx6 세포를 NREM 촉진 세포로 묶지 않는다.
- **BF SOM:** Xu 연구의 광자극은 휴식기 wake→NREM transition을 높였고, 일부 SOM 세포만 NREM-active였다 [6]. Anaclet 연구의 급성 activation/inhibition은 뚜렷한 NREM 또는 wake 촉진을 보이지 않았고, 제한적·시간 의존적인 arousal 영향이 관찰됐다 [7]. 따라서 확정 core가 아닌 **mixed evidence**다. BF에는 각성을 촉진하는 cholinergic·glutamatergic·PV 집단도 있다.
- **RT:** 기능 근거는 sector별 inhibitory bursting과 CaV3.3에 따른 local NREM spindle/delta 조절이다 [8]. Local slow wave나 spindle 변화는 전체 수면량 증가와 같은 endpoint가 아니며, quiet wake와 EEG/EMG로 판정한 NREM도 구분해야 한다.
- **PZ:** sleep-active inhibitory 집단의 필요성과 GABA 세포 activation에 의한 SWS/NREM·EEG SWA 촉진 근거가 있다 [9,10]. Lesion과 VGAT transmission disruption의 REM 효과는 다르며, 2014 저자들도 REM-specific 역할 부여를 경계했다. VGAT 표지를 GABA-only로 바꾸지 않는다.

BF 원문 조작은 caudal BF의 horizontal diagonal band·magnocellular preoptic nucleus·SI를 포함하며, rostral medial septum·vertical diagonal band와 구분된다 [6]. 현재 SI만으로 ‘전체 basal forebrain’을 대표한다고 하지 않으며, MA/NDB는 경계 한계를 붙인 별도 행으로 유지한다.

## Projection 결과의 해석 기준

`structure_id`가 있는 일곱 label은 개별 target으로 조회한다. PZ는 `include_in_unionize_query=false`, `structure_id=null`로 보존하며 미평가를 0으로 바꾸지 않는다. source injection의 순도·부피·Cre 조건과 target의 hemisphere·volume·측정 threshold를 기록한다. 서로 겹치는 parent/child label을 합산하여 신호를 중복 계산하지 않는다.

전체 영역 projection signal은 **영역으로 들어오는 해부학적 축삭 신호**다. 해당 영역의 수면촉진 세포에 대한 synaptic connectivity, excitation/inhibition의 기능 방향, 수면 유도, PBM 신호 전파는 이 값만으로 판정할 수 없다. 이 단계에서 가능한 결론은 ‘관측한 injection 조건에서 해당 atlas 영역에 어떤 projection signal이 검출되는가’이다.

예를 들어 ZI signal이 다른 target보다 크거나 MOp와 VISp에서 차이가 있더라도, 그 차이가 **Lhx6 세포로의 input 차이**인지는 미확인이다. Cortex→ZI input이 어떤 세포의 흥분·억제와 downstream disinhibition을 유발하는지, NREM·REM·wake에 어떤 영향을 주는지는 별도 회로 생리와 EEG/EMG가 필요하다. 후속 실험은 source→ZI projection을 구분한 tracing·projection-specific manipulation과 target Lhx6 확인을 조합하고, 시간대·PBM/sham 조건별 상태 변화를 측정한다. 단순 축삭 신호 순위를 NREM 촉진 능력 순위로 바꾸지 않는다.

현재 mapped label들의 unionize 행을 빠짐없이 확보하는 것과 **기능상 중요한 PZ가 ontology에 없어 미평가인 것**은 별개다. 전자는 API/측정 커버리지 검증이고 후자는 해부학적 표적 정의의 한계다. 검출 0도 이 조건의 assay·threshold 아래 신호가 없거나 작다는 관측이며, 투사 또는 수면 기능의 완전한 부재를 뜻하지 않는다.

기능 검증에는 target의 Gal/Vgat/Lhx6/SOM 등 세포 유형을 확인하는 tracing·조직학, synapse/회로 생리, EEG/EMG NREM·REM·wake 및 시간대별 개입이 추가로 필요하다. 사용자 자극 위치와 atlas source injection의 일치도 역시 별도 검증 항목이다.

## 직접 확인한 원 논문

1. Kroeger et al. (2018). *Galanin neurons in the ventrolateral preoptic area promote sleep and heat loss in mice.* Nature Communications. DOI [10.1038/s41467-018-06590-7](https://doi.org/10.1038/s41467-018-06590-7); [PMC6175893](https://pmc.ncbi.nlm.nih.gov/articles/PMC6175893/); PMID [30297727](https://pubmed.ncbi.nlm.nih.gov/30297727/).
2. Machado et al. (2022). *Median preoptic GABA and glutamate neurons exert differential control over sleep behavior.* Current Biology. DOI [10.1016/j.cub.2022.03.039](https://doi.org/10.1016/j.cub.2022.03.039); [PMC9090993](https://pmc.ncbi.nlm.nih.gov/articles/PMC9090993/); PMID [35385692](https://pubmed.ncbi.nlm.nih.gov/35385692/).
3. Gvilia et al. (2006). *Homeostatic Regulation of Sleep: A Role for Preoptic Area Neurons.* Journal of Neuroscience. DOI [10.1523/JNEUROSCI.2012-06.2006](https://doi.org/10.1523/JNEUROSCI.2012-06.2006); [PMC6674606](https://pmc.ncbi.nlm.nih.gov/articles/PMC6674606/); PMID [16971526](https://pubmed.ncbi.nlm.nih.gov/16971526/).
4. Liu et al. (2017). *Lhx6-positive GABA-releasing neurons of the zona incerta promote sleep.* Nature. DOI [10.1038/nature23663](https://doi.org/10.1038/nature23663); [PMC5958617](https://pmc.ncbi.nlm.nih.gov/articles/PMC5958617/); PMID [28847002](https://pubmed.ncbi.nlm.nih.gov/28847002/).
5. Vidal-Ortiz et al. (2023). *Unilateral optogenetic stimulation of Lhx6 neurons in the zona incerta increases REM sleep.* Sleep. DOI [10.1093/sleep/zsad217](https://doi.org/10.1093/sleep/zsad217); [PMC11502959](https://pmc.ncbi.nlm.nih.gov/articles/PMC11502959/); PMID [37599437](https://pubmed.ncbi.nlm.nih.gov/37599437/).
6. Xu et al. (2015). *Basal Forebrain Circuit for Sleep-Wake Control.* Nature Neuroscience. DOI [10.1038/nn.4143](https://doi.org/10.1038/nn.4143); [PMC5776144](https://pmc.ncbi.nlm.nih.gov/articles/PMC5776144/); PMID [26457552](https://pubmed.ncbi.nlm.nih.gov/26457552/).
7. Anaclet et al. (2018). *Genetic Activation, Inactivation, and Deletion Reveal a Limited And Nuanced Role for Somatostatin-Containing Basal Forebrain Neurons in Behavioral State Control.* Journal of Neuroscience. DOI [10.1523/JNEUROSCI.2955-17.2018](https://doi.org/10.1523/JNEUROSCI.2955-17.2018); [PMC5977448](https://pmc.ncbi.nlm.nih.gov/articles/PMC5977448/); PMID [29735555](https://pubmed.ncbi.nlm.nih.gov/29735555/).
8. Fernandez et al. (2018). *Thalamic reticular control of local sleep in mouse sensory cortex.* eLife. DOI [10.7554/eLife.39111](https://doi.org/10.7554/eLife.39111); [PMC6342525](https://pmc.ncbi.nlm.nih.gov/articles/PMC6342525/); PMID [30583750](https://pubmed.ncbi.nlm.nih.gov/30583750/).
9. Anaclet et al. (2012). *Identification and Characterization of a Sleep-Active Cell Group in the Rostral Medullary Brainstem.* Journal of Neuroscience. DOI [10.1523/JNEUROSCI.0620-12.2012](https://doi.org/10.1523/JNEUROSCI.0620-12.2012); [PMC3564016](https://pmc.ncbi.nlm.nih.gov/articles/PMC3564016/); PMID [23238713](https://pubmed.ncbi.nlm.nih.gov/23238713/).
10. Anaclet et al. (2014). *The GABAergic parafacial zone is a medullary slow wave sleep-promoting center.* Nature Neuroscience. DOI [10.1038/nn.3789](https://doi.org/10.1038/nn.3789); [PMC4214681](https://pmc.ncbi.nlm.nih.gov/articles/PMC4214681/); PMID [25129078](https://pubmed.ncbi.nlm.nih.gov/25129078/).

Ontology 출처: [Allen Structure graph 1 공식 API](https://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[graph_id$eq1],rma::options[num_rows$eqall]). 원본 응답과 checksum은 `input/connectivity/ontology_structure_graph1.json` 및 target definitions의 provenance 필드에 저장했다.
