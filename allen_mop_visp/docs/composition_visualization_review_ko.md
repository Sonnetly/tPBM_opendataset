# MOp L4와 세포 구성 그림의 표현 검토

2026-10-12. 사용자 요청에 따른 **검토 문서**이며 기존 데이터, plotting code, heatmap·bar 그림은 변경하지 않았다.

## 1. MOp L4가 비어 있는 이유

이번 Allen CCF ontology/parcellation에는 MOp1·MOp2/3·MOp5·MOp6a·MOp6b가 있고 별도 **MOp4 label은 없다**. VISp4는 존재한다. 공식 structure graph 1에서 MOp ID는 985, 그 하위 영역 ID는 320·943·648·844·882이며 VISp4 ID는 721이다. 따라서 MOp의 L4 열은 이 atlas의 해부학적 분류로 계산할 수 없는 항목이다. 분석 중 L4 세포를 삭제하거나 join에서 잃은 결과가 아니다. 선택 세포 145,860개에 annotation/CCF coordinate join 누락과 anatomical-layer unknown은 없었다.

현재 그림은 이 칸을 흰색으로 mask하고 축에는 `n=0` 및 `blank = not sampled`를 표시한다. 이러한 표현은 단순 표집 실패 또는 생물학적 부재로 오해될 수 있다. 추후 그림에서는 **`N/A: 해당 CCF label 없음`**으로 표시하는 편이 정확하다. 숫자 0인 세포 유형 비율과 같은 색·기호를 쓰지 않는다.

Motor cortex에서 L4와 유사한 세포·회로가 없다는 뜻도 아니다. 실제 transcriptomic `006 L4/5 IT CTX Glut`은 MOp **10,232개**, 전체 뉴런의 **18.05%**이며 VISp는 16,551개·34.49%다. MOp의 해당 subclass는 CCF L2/3와 L5 등에 등록되어 있다. 이 transcriptomic 이름을 해부학적 L4로 바꾸거나 이 세포를 모두 새 MOp4로 재분류하지 않는다.

Yamawaki et al. (2014), *A genuine layer 4 in motor cortex with prototypical synaptic circuit connectivity*, eLife, DOI [10.7554/eLife.05422](https://doi.org/10.7554/eLife.05422)는 mouse forelimb M1의 얇은 L3/5A 경계 영역에서 L4형 회로를 보고했다. 이는 M1을 단순히 ‘L4 관련 세포가 없는 영역’으로 해석하지 않게 하는 문헌이며, 현재 MERFISH subclass와 논문의 세포를 일대일로 대응하는 근거는 아니다.

## 2. 현재 heatmap의 색상 압축

`src/merfish/run_merfish.py`의 `plot_composition`은 각 층의 **전체 QC annotated 세포**를 분모로 계산하고, 두 영역에 공통 `viridis`, 선형 `vmin=0`, `vmax=100`을 적용한다. 대부분 층의 Glut은 57.52–80.78%이고 작은 세포군은 낮은 비율이므로 비슷한 보라색으로 보인다. Glut 비율 때문에 다른 세포군의 값이 잘못 계산된 것은 아니며, 낮은 값에서 색상 구별이 부족한 문제다. Glut 행만 지우고 0–100% 색 범위를 유지하면 이 문제는 그대로다.

Glut 제외 그룹 중 전체 최대는 VISp L1 astrocyte 37.14%이고, L1을 제외하면 MOp L6b oligodendrocyte 17.25%다. 별도 nonGlut heatmap에서 **원래 전체 세포 분모를 유지하고 양쪽에 같은 0–40% 색 범위·숫자 label**을 쓰는 방안이 합리적이다. 희소 그룹의 작은 차이를 더 보려면 해당 그룹별 bar/dot를 쓰거나 별도 범위를 명시한 패널로 나눈다. 영역별 자동 색 범위는 같은 색의 의미가 달라지므로 피한다.

## 3. Glut을 빼는 두 방법은 다른 질문이다

표시에서 Glut을 빼되 분모를 유지하면 각 그룹의 값은 그대로 `100 × 해당 세포 수 / 전체 층 세포 수`다. 이 동물에서 표집된 QC 세포 구성의 지역 차이를 계속 비교할 수 있다.

반면 nonGlut만으로 다시 100%를 만들면 `100 × 해당 세포 수 / (전체 세포 수 − Glut 수)`가 된다. 이것은 **nonGlut 내부 구성**이라는 다른 질문이다. 분모가 바뀌므로 지역 비교 방향도 바뀔 수 있다.

| 영역 전체 GABA 비율의 분모 | MOp | VISp |
|---|---:|---:|
| 전체 QC 세포: 83,976 / 61,884 | 9.21% | 8.07% |
| nonGlut: 35,035 / 18,889 | 22.07% | 26.43% |

둘 다 계산상 맞지만, 뒤의 값을 ‘전체 세포에서 VISp에 GABA가 더 많다’고 해석하면 틀린다. 전체 분모를 유지하면서 Glut만 뺀 stack의 영역 전체 합은 MOp 41.72%·VISp 30.52%이며, 이 그림을 100% stacked bar라고 부르면 안 된다. 비율은 세포 density나 절대 조직 내 세포 수가 아니다.

## 4. 권장 표현

| 목적 | 권장 그림 | 분모·표시 원칙 |
|---|---|---|
| 작은 세포군의 MOp–VISp 직접 비교 | **층별 grouped bar 또는 dot plot** | 전체 층 세포 분모 유지; 세포군별 패널, 같은 층의 두 영역을 나란히 표시 |
| 전체 part–whole 구성 | 층별 100% stacked bar, Glut 포함 | 전체 분모, 고정 범주 색, 원본 n 표시; 작은 차이의 주 비교에는 별도 bar 필요 |
| Glut을 제외한 원래 비율을 한 번에 확인 | 별도 nonGlut heatmap | 전체 분모 유지, 공통 0–40% 범위와 수치 label |
| nonGlut 내부 상대 구성 | 별도 nonGlut 100% stacked bar | 재정규화했다는 명시, nonGlut n 표시; 원래 전체 비율과 구분 |
| 공유 층에서 지역 차이만 확인 | MOp−VISp percentage-point heatmap | 0 중심 공통 발산 색 범위; 단일 동물의 기술적 차이이며 유의성 아님 |
| 한 층의 대략적 구성 소개 | Pie 가능, 주 비교에는 비권장 | 작은 조각·여러 층의 각도 비교가 어려움 |

주 비교는 **전체 분모를 유지한 nonGlut grouped bar**를 권한다. 전체 구성은 Glut 포함 stacked bar로 보조하고, nonGlut 재정규화는 내부 구성이라는 별도 목적이 있을 때만 추가한다. 공유하는 L1·L2/3·L5·L6a·L6b를 직접 비교하고 VISp L4는 단독 표시한다. MOp L4에 0을 대입하거나 VISp L4를 L5로 임의 합산하지 않는다. 해부학적 층과 transcriptomic subclass 비교는 계속 분리한다.

현재 자료는 WT B6J 수컷 한 동물이다. Bar의 각 세포·절편을 독립 동물로 취급한 error bar나 P value는 만들지 않는다.

## 검토에 사용한 정본

- [실제 세포 구성 보고서](../results/merfish/composition_report_ko.md)
- [영역 구성 CSV](../results/merfish/composition_region_major.csv), [층 구성 CSV](../results/merfish/composition_layer_major.csv)
- [뉴런 subclass CSV](../results/merfish/composition_region_subclass_neurons.csv)
- [공식 ontology snapshot](../input/connectivity/ontology_structure_graph1.json)
- [현재 layer heatmap](../figures/merfish/layer_major_composition.png), [plot 코드](../src/merfish/run_merfish.py)
