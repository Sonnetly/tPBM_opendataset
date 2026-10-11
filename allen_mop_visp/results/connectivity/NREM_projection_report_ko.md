# MOp–VISp에서 수면 관련 영역으로의 Allen 투사 분석

2026-10-11 UTC에 확보한 **실제 Allen Mouse Connectivity Atlas 원본 API 자료**를 분석했다. WT C57BL/6J의 MOp 주입 8개, VISp 주입 33개를 비교하면 ZI·RT에 segmented axonal signal이 관측되고, SI에도 대부분의 실험에서 신호가 있다. ZI의 주입량 정규화 신호 중앙값은 MOp에서 높으며 아래의 순도·연령 민감도 분석에서도 같은 방향이다.

**이 결과만으로 MOp나 VISp가 NREM을 촉진하는 세포에 기능적으로 연결된다고 결론낼 수는 없다.** 핵심 preoptic 표적인 VLPO·MEPO는 중앙값이 작고 관측된 0이 많으며, VLPO 지역 비교 방향은 주입 선택 조건에 따라 바뀐다. 가장 큰 preoptic 신호를 보이는 두 실험에는 주변 영역 주입 spread가 상당하다. PZ는 exact atlas label이 없어 미평가다. ZI와 BF의 수면 기능은 특정 세포 유형과 실험 조건에 의존하며 상반된 문헌도 있어, region-level signal을 NREM 촉진 능력으로 환산하지 않았다.

## 실제 원본과 실험 선택

공식 HTTPS RMA API의 structure graph 1, `ApiConnectivity`, `SectionDataSet`, `ProjectionStructureUnionize`를 사용했다. 공개 SDK와 지표 문서를 원문으로 보존했다. 제품 5(`MouseConnectivityProjection`)에서 primary injection structure가 MOp(985)·VISp(385) 또는 각각의 ontology 자손인 실험을 선택했다. 현재 snapshot에는 자손을 primary로 갖는 추가 실험이 없어 아래의 168개 모두 parent primary와 일치한다. retinal projection 제품 31과 기타 제품은 제외했다.

| Source | WT C57BL/6J | Cre line: 별도 탐색 | WT donor | WT 주입 연령 |
|---|---:|---:|---:|---|
| MOp | 8 | 23 | 8마리, 모두 male | 모두 P56 |
| VISp | 33 | 104 | 33마리, 모두 male | P49–P80 |

WT는 API와 donor strain이 모두 C57BL/6J이고 Cre 정보가 없는 실험이다. WT 41개의 donor ID는 서로 다르다. 선택된 168개 실험의 manual `failed` flag는 모두 false다. 주입 quality 원본 코드, viral construct, 주입 방법, age, days post injection을 보존했으나 공식 quality 코드의 점수 의미를 별도로 검증하지 않았으므로 코드 값에 따른 cutoff를 적용하지 않았다. WT 사이에도 construct·연령·주입 프로토콜·관찰 기간이 달라 **연령과 주입 조건이 일치하는 지역 효과 추정**으로 해석할 수 없다. Cre 실험은 각 line별로 따로 요약했다.

주입 spread는 전체 injection signal volume 중 source parent에 포함된 signal fraction으로 계산했다. Parent와 layer child를 합산하지 않았다.

| WT source | Source signal fraction: 최소 / 중앙값 / 최대 | Source 밖 annotation이 있는 실험 | Source-only annotation 실험 | 전체 injection signal volume 범위, mm³ |
|---|---|---:|---:|---|
| MOp | 0.5279 / 0.8211 / 0.9999 | 6/8 | 2/8 | 0.1726–1.0151 |
| VISp | 0.3093 / 0.9690 / 1.0000 | 17/33 | 16/33 | 0.03872–1.0830 |

`source-only annotation`은 원본 `injection_structures` 목록이 source 안의 ID만 포함한다는 뜻이다. 수치상 source fraction이 완전히 1이라는 뜻은 아니다. 전체 MOp/VISp 주입 실험을 다루며, 사용자의 광섬유 자극 위치나 400 μm core에 일치하는 injection ROI로 간주하지 않았다.

## 표적의 해부학적 범위와 기능 근거

상세 문헌·DOI와 상반된 결과는 [표적 정의 문서](../../docs/NREM_projection_target_definition_ko.md)와 [원본 정의 JSON](../../input/connectivity/nrem_target_definitions.json)에 보존했다. Atlas region mapping과 수면 기능을 보고한 세포 유형의 대응은 별도 문제다.

| Target | Allen ID | 분석 구분 | 해석 범위 |
|---|---:|---|---|
| VLPO | 689 | Core | Gal/GABA 세포의 NREM 근거가 있으나 tracer가 그 세포에 도달하는지 미측정 |
| MEPO, MnPO alias | 452 | Core | Preoptic 수면 관련 세포 기능과 regional signal을 구분 |
| ZI | 797 | 세포 유형·조건 의존 | Lhx6 세포의 NREM+REM 근거와 후속 REM 중심 결과가 함께 존재 |
| SI | 342 | BF의 세포 유형별 혼합 근거 | SOM 관련 NREM 근거와 급성 효과가 없었던 결과가 공존; 전체 BF를 대표하지 않음 |
| MA | 298 | 선택적 지역 맥락 | Allen Magnocellular nucleus; 문헌 MCPO와의 정확한 대응은 미검증 |
| NDB | 596 | 선택적 지역 맥락 | Horizontal/vertical diagonal band가 혼재하여 문헌 caudal BF 범위와 일치하지 않음 |
| RT | 262 | NREM 관련, core와 구분 | Thalamic reticular nucleus; canonical NREM switch로 합산하지 않음 |
| PZ | 없음 | Core, atlas 미매핑 | Exact PZ label 없음; PGRN으로 치환하지 않고 미평가 유지 |

Thalamic reticular nucleus는 **RT(262)**다. TRN(574)은 tegmental reticular nucleus이므로 사용하지 않았다. 매핑된 7개 표적 사이에 ancestor/descendant 중복은 없으며, 서로 다른 표적을 합산한 수면 점수나 NREM 능력 순위를 계산하지 않았다.

## 지표와 관측 여부

주 지표는 양측(`hemisphere_id=3`)의 `normalized_projection_volume`, 이하 NPV다.

\[
\mathrm{NPV} = \frac{\text{target의 projection signal volume}}{\text{injection site 전체의 signal volume}}
\]

분자 `projection_volume`과 분모 `injection_volume`의 단위는 mm³이며, NPV는 무차원이다. Target 전체 부피를 분모로 쓰는 비율이나 연결 확률·synapse 수·수면 촉진 효과가 아니다. NPV가 1보다 클 수도 있다. Target `projection_density`는 valid target voxel 중 segmented signal 비율이라는 별도 지표다. `is_injection=false` unionize는 injection-site voxel 밖의 signal만 포함한다.

Hemisphere 1은 left, 2는 right, 3은 양측 aggregate다. 3에 1과 2를 더하면 중복이므로 요약에는 3만 사용했다. 실험별 관측된 injection root signal로 주입 방향을 확인했으며 현재 선택 실험은 모두 right다. Root의 left injection unionize record가 없는 상태를 숫자 0으로 바꾸지 않고 원본 결측으로 유지했다.

매핑된 표적은 168개 실험 × 7개 표적 × 3개 반구 = **3,528개 예상 record 모두 관측**되었고 중복 key는 없다. API가 반환한 실제 0과 미반환 record를 구분했다. PZ의 504개 placeholder 행은 결측·미평가다. 아래의 `positive`는 오직 원본 segmented projection volume이 0보다 크다는 뜻이며, 잡음·경계·등록 오차를 제외한 검증된 projection detection 기준은 아니다. 별도의 신호 cutoff를 만들지 않았다.

## WT 전체 주입의 양측 결과

각 값은 실험별 NPV의 중앙값이다. Positive는 양측 segmented signal이 0보다 큰 실험 수다. 작은 값의 지역 간 배수는 분모에 민감하므로 아래에서는 절대값을 표시한다.

| Target | MOp NPV 중앙값, n=8 | VISp NPV 중앙값, n=33 | Positive MOp / VISp | 관측된 0 MOp / VISp |
|---|---:|---:|---|---|
| VLPO | 4.158e-08 | 8.231e-06 | 6/8 · 26/33 | 2 · 7 |
| MEPO | 2.158e-07 | 6.701e-09 | 4/8 · 17/33 | 4 · 16 |
| ZI | 0.1705 | 0.02560 | 8/8 · 33/33 | 0 · 0 |
| SI | 0.003025 | 9.752e-05 | 8/8 · 31/33 | 0 · 2 |
| RT | 0.1025 | 0.09633 | 8/8 · 33/33 | 0 · 0 |
| MA | 2.771e-05 | 2.641e-08 | 6/8 · 17/33 | 2 · 16 |
| NDB | 4.294e-05 | 1.628e-04 | 8/8 · 33/33 | 0 · 0 |
| PZ | 미평가 | 미평가 | 해당 없음 | 해당 없음 |

ZI의 중앙값은 MOp에서 높고 RT는 두 source에서 비슷한 수준이다. 이는 공개 tracer 주입 표본의 기술 통계다. ZI의 postsynaptic Lhx6 세포, SI의 SOM 세포, 특정 transmitter identity를 분리하지 않았고, synapse나 기능적 회로 전달을 확인하지 않았다. 따라서 해당 값으로 지역의 NREM 촉진 능력이나 PBM 반응을 비교할 수 없다.

## 주입 spread와 P56 민감도

Source fraction ≥0.8 조건은 분석에서 정한 **탐색적 주입 spread 민감도 조건**이다. Official quality 기준, 신호 검출 기준 또는 기능적 근거 기준이 아니다. Source-only annotation, P56, fraction 조건은 서로 다른 선택이므로 혼합해서 하나의 matched cohort라고 부르지 않는다.

| WT 선택 조건 | MOp / VISp n | ZI NPV MOp / VISp | VLPO NPV MOp / VISp | MEPO NPV MOp / VISp |
|---|---|---|---|---|
| 전체 primary source 주입 | 8 / 33 | 0.1705 / 0.02560 | 4.158e-08 / 8.231e-06 | 2.158e-07 / 6.701e-09 |
| Source-only annotation | 2 / 16 | 0.4745 / 0.03270 | 2.480e-05 / 6.042e-07 | 5.460e-06 / 0 |
| 정확 P56 | 8 / 17 | 0.1705 / 0.01938 | 4.158e-08 / 4.704e-05 | 2.158e-07 / 4.227e-07 |
| Source fraction ≥0.8 | 4 / 24 | 0.2684 / 0.02782 | 1.005e-06 / 6.042e-07 | 2.158e-07 / 0 |
| 정확 P56 및 source fraction ≥0.8 | 4 / 10 | 0.2684 / 0.04803 | 1.005e-06 / 3.720e-05 | 2.158e-07 / 6.306e-08 |

ZI에서 MOp 중앙값이 높은 방향은 다섯 조건 모두에서 유지된다. P56와 source fraction으로 제한해도 viral construct·프로토콜 등 다른 차이는 남으며 표본 수도 줄어든다. 유의성 검정이나 matched regional effect는 추정하지 않았다. VLPO는 전체 WT에서 VISp가 높지만 source-only 및 fraction 조건에서는 MOp가 높아져, 작은 preoptic signal로 일관된 source 순위를 만들기 어렵다.

두 source의 VLPO·MEPO 최대 신호 실험에는 다음과 같은 spread가 있다.

| Source | Experiment | Source fraction | Source 밖 annotation | VLPO / MEPO 최대 NPV |
|---|---:|---:|---|---|
| MOp | [584903636](https://connectivity.brain-map.org/projection/experiment/584903636) | 0.5279 | MOs | 9.964e-05 / 1.395e-04 |
| VISp | [114008926](https://connectivity.brain-map.org/projection/experiment/114008926) | 0.3093 | CA1, SUB, DG, ProS | 0.02511 / 0.01470 |

이 실험의 preoptic signal을 순수 MOp 또는 VISp 투사에 귀속할 수 없다. 미량 신호의 segmentation·registration 검토와 더 좁은 source tracer validation이 후속 해부학 검증의 대상이다. Allen region signal이 없거나 작은 경우에도 실험의 검출 한계와 표본 조건 때문에 생물학적 연결 부재로 단정하지 않았다.

## 재현과 검증

Raw 응답 11개(ontology·실험 metadata 2개와 3개 batch의 세 API 종류)를 `input/connectivity/`에 보관했다. 각 원본의 URL·확보 또는 기존 파일 확인 UTC 시각·응답 byte 수·SHA256은 sidecar와 `acquisition_manifest.json`에 기록했다. 최초 ontology·ApiConnectivity 두 응답은 기존 snapshot을 확인하면서 provenance를 작성했으므로 최초 HTTP status와 최초 조회 시각은 기록되지 않았다. 나머지 아홉 응답은 실제 HTTPS GET의 status 200과 조회 시각을 기록했다. Cached 재사용 시 provenance를 유지하고, raw hash 불일치나 다른 URL이면 오류로 종료한다. 원본을 자동 교체하지 않는다.

실제 raw arithmetic과 요약 CSV를 별도 코드로 재계산한 [독립 검증](../validation/connectivity_actual_data_review.json)에서 1,101개 검증이 통과했다. NPV의 전체 injection signal volume 분모, 원본 지표의 보존, 실제 WT 41개·Cre 분리·다섯 민감도 요약·Excel row 수를 확인했다. Cached provenance 보존/변조 거부는 별도의 synthetic fixture로 검증했으며 실제 atlas 결과의 근거와 구분했다.

```bash
python src/connectivity/run_connectivity.py          # 보존된 실제 원본에서 재현
python src/connectivity/build_connectivity_workbook.py
python src/connectivity/check_provenance_cache.py    # synthetic integrity test
python src/connectivity/validate_connectivity.py     # 실제 자료 독립 재계산
```

공식 metric 정의는 [보존된 문서](../../input/connectivity/official_unionizes.rst)와 [AllenSDK 원문](https://raw.githubusercontent.com/AllenInstitute/AllenSDK/master/doc_template/unionizes.rst)에 있다. API는 [Allen RMA endpoint](https://api.brain-map.org/api/v2/data/query.json)를 사용했다. 실행 방법과 입력 검증 동작은 [코드 README](../../src/connectivity/README.md)에 정리했다.

검토용 결과:

- [Excel: 실제 실험·표적·반구와 민감도 요약](MOp_VISp_NREM_projection.xlsx)
- [실험별 target×hemisphere 지표 CSV](projection_by_experiment_target_hemisphere.csv)
- [WT 전체 및 민감도 기술 통계 CSV](WT_target_projection_descriptive_summary.csv)
- [실험별 주입 metadata·spread audit CSV](selected_experiment_metadata_and_injection_audit.csv)
- [Cre line별 별도 탐색 CSV](Cre_line_separate_exploratory_summary.csv)
- [WT 실험별 heatmap](../../figures/connectivity/WT_individual_target_projection_heatmap.png)
- [WT 주입 source fraction 그림](../../figures/connectivity/WT_injection_source_purity.png)

이 분석은 공개 anterograde tracer의 해부학적 신호를 측정한다. 기존 MERFISH·RNA-seq 표본 수나 PBM 처리 효과와 합치지 않았다. 수면 상태·처리 효과·postsynaptic identity를 측정하지 않으므로 실제 NREM 기능의 확인에는 별도의 세포 유형과 회로 실험이 필요하다.
