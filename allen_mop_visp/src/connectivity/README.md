# Allen Mouse Connectivity Atlas 분석

실제 공개 API 자료로 MOp/VISp 주입에서 수면 관련 atlas 영역으로 가는 segmented axonal signal을 비교한다. NREM 기능을 입증한 세포 유형과 region-level tracer signal을 구분한다. 기존 RNA-seq·MERFISH 결과를 이 tracer 실험의 표본 수나 처리 효과로 합치지 않는다.

프로젝트 최상위에서 실행한다.

```bash
# 최초: 없는 원본 API 응답만 공식 HTTPS에서 확보
python src/connectivity/run_connectivity.py --download

# 보존된 원본 응답·SHA256으로 분석 재현
python src/connectivity/run_connectivity.py

# 검토용 Excel 생성
python src/connectivity/build_connectivity_workbook.py

# cached provenance 보존과 변조된 응답 거부를 synthetic fixture로 검증
python src/connectivity/check_provenance_cache.py

# 분석 모듈을 import하지 않고 실제 원본·요약·Excel을 독립 대조
python src/connectivity/validate_connectivity.py
```

클라우드 Python 실행 파일은 `/workspace/Allen Dataset_M1_V1/.venv/bin/python`을 사용했다. Windows의 기존 프로젝트 폴더에서는 `setup_windows.ps1`로 만든 Python 환경을 사용한다. 기본 라이브러리는 numpy, pandas, matplotlib, openpyxl이며 AllenSDK 설치 없이 공식 SDK에 문서화된 RMA REST endpoint를 사용한다.

`acquire_connectivity.py`는 그래프1 ontology와 ApiConnectivity 실험을 확보하고 product5에서 primary structure가 MOp985/VISp385 또는 그 자손인 실험만 선택한다. 현재 snapshot에서 primary child로 추가되는 실험은 없다. product31은 retinal projection이며 다른 제품과 함께 제외한다. 168개 실험을 80개씩 나누어 manual metadata, injection unionize 전체 구조, curated target noninjection unionize를 확보했다. API 응답의 success 및 반환 행 수/total_rows가 일치해야 사용한다. TLS를 유지하고 URL·응답 byte 수·SHA256·조회 시각을 함께 저장한다.

이미 보관한 응답은 initial provenance의 SHA256과 URL이 일치해야 재사용한다. 초기 조회 timestamp와 HTTP status를 유지하며 마지막 cached-load 시각을 별도로 기록한다. 원본 hash가 바뀌면 오류로 종료하며 원본과 provenance를 덮어쓰지 않는다. `--download`는 기존 frozen 원본을 새 API 응답으로 자동 교체하지 않는다. 동적 API source는 `input/connectivity/acquisition_manifest.json`의 응답 snapshot으로 재현한다.

`analyze_connectivity.py`는 다음을 계산한다.

- 원본 manual failed flag, donor·strain·sex·Cre line·age·주입 quality raw code·viral materials·days post injection 보존.
- whole-brain injection signal volume 대비 MOp/VISp **parent**에 든 injection signal fraction. Parent와 layer child를 합산하지 않는다.
- WT B6J와 Cre line을 분리하고, WT 전체·source-only annotated labels·P56·순도≥0.8·P56와 순도≥0.8의 기술 통계. 0.8은 탐색적 주입 spread 민감도 조건이며 official quality cutoff가 아니다.
- target×실험×hemisphere별 원본 projection_volume, normalized_projection_volume, density·energy 및 관측 여부.
- `NPV = target projection_volume / total injection signal volume` 일치 확인. 지표는 무차원이며 비율%, 연결 확률, synapse 수, 기능적 NREM 효과가 아니다.
- hemisphere1left/2right/3both를 보존하며 3을 1/2에 더하지 않는다. Ipsilateral/contralateral 표시에는 실험의 관측된 injection-side root signal을 사용한다. 현재 모든 선택 실험에서 관측된 root injection hemisphere는2다. 존재하지 않는 root hemisphere1 injection record를 숫자0으로 바꾸지 않는다.
- exact mapped targets 7개에 대한 record coverage 확인. PZ에는 exact Allen term이 없어 projection 지표를 결측·미평가로 유지한다. target hierarchy 중복이 있으면 audit에 기록하고 표적들을 합산한 sleep score를 만들지 않는다.

공식 metric·API 해석 문서는 `input/connectivity/official_unionizes.rst`와 보관한 SDK source다. 문헌에 근거한 표적 정의는 scientific review에서 작성한 `input/connectivity/nrem_target_definitions.json` 및 `docs/NREM_projection_target_definition_ko.md`에 보존했다.

주요 결과는 `results/connectivity/NREM_projection_report_ko.md`, `MOp_VISp_NREM_projection.xlsx`, `projection_by_experiment_target_hemisphere.csv`, `WT_target_projection_descriptive_summary.csv`다. 그림은 `figures/connectivity/WT_individual_target_projection_heatmap.png`와 `WT_injection_source_purity.png`다. 수면 상태·PBM 처리·postsynaptic identity가 측정되지 않아 기능적 효과나 표적 우선순위는 별도 실험의 대상이다.
