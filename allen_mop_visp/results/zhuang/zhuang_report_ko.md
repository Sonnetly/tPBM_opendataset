# Zhuang-ABCA MERFISH MOp와 VISp 비교 보고서

작성일 2026-10-15. 실행 대상은 Zhuang-ABCA-1과 Zhuang-ABCA-2이며 Zhuang-ABCA-3과 Zhuang-ABCA-4는 아직 포함하지 않았다.

## 데이터셋

- 출처는 Allen Brain Cell (ABC) Atlas에 공개된 Zhuang lab whole mouse brain MERFISH 자료(Zhang et al., Nature 2023)이다. 공식 S3 release manifest `20250531`을 사용했고 manifest SHA256은 `d1c3f8d8c4e6657dbb8da53e6403f11e6116a125a3682dba9180cc6e998f24cd`이다.
- 동물은 C57BL/6NCrl(B6N) 배경이다. 기존 주 공간 자료인 B6J MERFISH-C57BL6J-638850과 계통이 다르므로 두 자료를 합산하지 않고 B6J의 독립 동물 반복으로 간주하지 않는다.
- 측정 panel은 1,122개 유전자이다. 세포 label은 Allen WMB taxonomy(class, subclass, supertype, cluster)이고 위치는 Allen CCF 2020 parcellation으로 주어진다.
- 사용한 파일은 동물별 `cell_metadata_with_cluster_annotation`, `ccf_coordinates`, `log2` 발현 행렬(h5ad)과 CCF parcellation 대응표이다. 모든 파일은 크기와 공개 MD5를 검증한 뒤 사용했다. 원본(약 4.2 GB)은 `input/zhuang/cache/`에 두고 Git에서 제외했으며 live manifest 사본은 `input/zhuang/live_manifest_20250531.json`에 보존했다.
- 이 atlas에는 PBM 처리 조건이 없다. 정상 상태 지역 차이만 기술하며 PBM 반응을 추정하지 않는다.

## 수행한 작업

1. `src/merfish/run_zhuang.py --list`로 manifest를 확인했다. ABCA-1과 ABCA-2, parcellation 대응표를 합한 다운로드 크기는 4.24 GB였다.
2. `--download --animals Zhuang-ABCA-1 Zhuang-ABCA-2`로 다운로드와 분석을 실행했다.
3. 첫 실행은 ABCA-1 발현 행렬 로드 단계에서 메모리 초과로 종료됐다(exit 137). 이 log2 행렬은 sparse가 아닌 dense float32 4,167,870 × 1,122 배열(gzip chunk, 압축 해제 시 약 18.7 GB)이며 anndata backed 모드의 fancy indexing이 7 GB RAM을 넘었다.
4. `load_expression`을 수정해 h5py로 저장 chunk 경계에 맞춘 행 블록을 순서대로 읽고 MOp와 VISp 세포 행과 필요한 유전자 열만 남기도록 했다. 반환 값과 세포 순서는 이전 코드와 같고 최대 메모리만 달라진다. 기존 `src/merfish/test_zhuang_logic.py` 9개 테스트를 통과했다.
5. 수정 후 재실행은 정상 종료됐다(status `zhuang_within_subclass_analysis_completed`).

## 분석 방법 요약

- CCF parcellation에서 MOp와 VISp에 속한 세포만 선택했다.
- 두 영역에 모두 있고 영역별 50개 이상인 subclass 안에서만 MOp 평균 빼기 VISp 평균을 계산했다. subclass 구성 차이가 발현 차이로 보이는 것을 막기 위한 설계이다. 단위는 공식 log2 발현 평균의 차이이며 log2 fold change가 아니다.
- 유전자별 cell-class discordance 점수는 GABA subclass와 Glut subclass에서 차이 부호가 반대로 일관된 정도이다. panel 전체 1,122개 유전자에 대해 계산해 후보 유전자를 panel 전체 분포 안의 순위로 읽는다.
- 동물 수준 p value는 계산하지 않고 동물 간 부호 일치를 보고한다.

## 결과

### 규모

| 동물 | MOp 세포 | VISp 세포 | 비교 subclass | 제외 strata |
|---|---|---|---|---|
| Zhuang-ABCA-1 | 48,400 | 38,475 | 23 | 32 |
| Zhuang-ABCA-2 | 21,194 | 15,410 | 19 | 27 |

### 후보 유전자 측정 범위

탐색 후보 62개 중 12개가 panel에 있다. Adora1, Adora2a, Ccbe1, Cdkn1c, Ecel1, Epha4, Npnt, Penk, Slc10a4, Stxbp6, Trpc5, Zeb2이다. Adk, Nt5e, Entpd1, Slc29a1, Slc29a2, Ada 등 나머지 adenosine 경로 유전자는 측정되지 않았으며 0 발현을 뜻하지 않는다.

### Adora1

| subclass | ABCA-1 | ABCA-2 |
|---|---|---|
| L2/3 IT CTX Glut | −0.121 | −0.089 |
| L4/5 IT CTX Glut | −0.116 | −0.099 |
| L6 CT CTX Glut | −0.070 | −0.069 |
| Pvalb Gaba | +0.088 | +0.089 |
| Sst Gaba | +0.034 | +0.002 |
| Vip Gaba | +0.008 | −0.010 |
| Lamp5 Gaba | +0.017 | −0.018 |
| discordance 순위(1,122개 중) | 63 | 545 |

- IT와 CT 계열 Glut subclass에서 MOp가 VISp보다 낮은 방향은 두 동물에서 같다.
- GABA 쪽에서 두 동물이 일치하는 것은 Pvalb의 양의 방향뿐이다. Sst, Vip, Lamp5는 0 근처에서 부호가 갈린다. 10x에서 관찰한 GABA 전반의 양의 방향은 이 자료에서 Pvalb에 국한된다.
- ABCA-2에서 discordance 순위가 떨어진 것은 GABA subclass 4~5개의 부호 다수결 점수에서 Pvalb 신호가 희석된 결과로 해석된다.
- 차이 크기는 log2 단위 0.07~0.14로 작다.

### 동물 간 재현성

- 두 동물에서 모두 계산된 subclass×gene 21,318개 strata에서 MOp−VISp 차이의 Pearson r은 0.746이다. 차이 절댓값 평균이 0.05 이상인 2,764개 strata의 93.5%가 두 동물에서 같은 부호이다.
- 반면 panel 전체 discordance 점수의 동물 간 Spearman ρ는 0.095이다. 이 점수는 GABA subclass 수가 적어 불안정하므로 후보 순위의 근거로 단독 사용하지 않는다. subclass별 차이 자체를 우선 읽는다.

### 주의할 후보

Zeb2는 discordance 순위가 두 동물 모두 상위(6위, 1위)이지만 Glut 쪽 차이(+0.23~+0.26)가 크다. subclass 안의 supertype 구성 차이일 가능성을 먼저 확인해야 한다. MOp에는 고전적 L4가 거의 없으므로 L4/5 IT 결과에도 같은 우려가 적용된다.

## 산출물

- `zhuang_report.json` 실행 조건, manifest, 세포 수, 측정 후보
- `candidate_panel_measured_coverage.csv` 후보 62개의 측정 여부
- `within_subclass_region_difference.csv.gz` 동물×subclass×gene의 MOp, VISp 평균과 차이
- `panel_discordance_scores.csv.gz` panel 전체 discordance 점수와 순위
- `candidate62_discordance_ranks.csv` 측정 후보의 순위와 백분위
- `cross_animal_sign_agreement.csv.gz` subclass×gene별 동물 간 부호 일치
- `Zhuang-ABCA-*_dropped_strata.csv` 제외 strata와 사유
- `figures/zhuang/zhuang_candidate_by_subclass.png`, `figures/zhuang/zhuang_panel_discordance_null.png`

## 한계

- 동물 2마리의 기술 통계이며 효과 크기 추정이나 유의성 검정이 아니다.
- B6N 계통이며 B6J 주 분석과 합산하지 않는다.
- subclass 수준 비교이므로 subclass 내부 구성 차이는 아직 분리하지 않았다.
- PBM 처리 조건이 없어 PBM 반응을 측정하지 않는다.

## 재현

```bash
.venv/bin/python src/merfish/run_zhuang.py --list
.venv/bin/python src/merfish/run_zhuang.py --download --animals Zhuang-ABCA-1 Zhuang-ABCA-2
.venv/bin/python -m pytest -q src/merfish/test_zhuang_logic.py
```
