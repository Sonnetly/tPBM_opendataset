# B6J MERFISH MOp–VISp 분석

현재 상태: **공식 입력 다운로드와 실제 MOp–VISp 구성 분석을 완료**했다. 초기에는 공식 Allen S3 호스트가 클라우드 outbound proxy의 `HTTP 403 CONNECT tunnel failed`로 차단되었으나, 사용자 요청 후 재확인한 현재 런타임에서는 HTTPS GET이 성공했다. 공식 필수 입력 5개(1,571,075,446 bytes)를 다운로드하고 원본 MD5와 크기를 검증했다. MOp 83,976개, VISp 61,884개, 합계 145,860개 QC annotated 세포를 분석했으며 annotation 및 CCF 좌표 결측은 0개였다. 결과는 `results/merfish/composition_report_ko.md`에서 확인한다. 선택 사항인 raw matrix와 measured Penk 발현은 아직 다운로드·계산하지 않았다.

## 입력과 출처

Allen ABC의 `MERFISH-C57BL6J-638850`은 C57BL/6J WT 수컷 **한 동물**의 공간 자료다. 원본 클래스와 subclass를 그대로 보존하며, `WMB-taxonomy 20231215`의 표시용 큰 그룹을 별도 열로 만든다.

- 공식 SDK: <https://github.com/AllenInstitute/abc_atlas_access>
- SDK 확인 commit: `b16ee223d3466a57b30066973413ccc948c2e983`
- 확보한 문서: SDK `descriptions/MERFISH-C57BL6J-638850.md`, `notebooks/merfish_ccf_registration_tutorial.ipynb`, `notebooks/ccf_and_parcellation_annotation_tutorial.ipynb`
- 최초 확보한 manifest: SDK의 `tests/abc_atlas_cache/resources/manifest.json`에 저장된 공식 URL·크기·MD5 목록. **이 manifest 자체에 실제 세포 자료는 포함되어 있지 않다.** 원문 사본은 `input/merfish/official_sdk_manifest_snapshot.json`, 초기 기록은 `input/merfish/provenance.json`의 `initial_SDK_snapshot_record`에 보존했다. 실제 분석에는 S3에서 확보한 `input/merfish/live_manifest_20250531.json`과 checksum을 통과한 5개 원본 파일을 사용했다.
- 자동 다운로드 시 같은 release의 실제 S3 manifest를 먼저 읽고, 각 파일의 공식 MD5 및 크기를 검증한다. TLS 검증을 유지한다.

| 필수 입력 | 버전 | byte |
|---|---:|---:|
| `cell_metadata_with_cluster_annotation` | 20241115 | 1,018,613,800 |
| `reconstructed_coordinates` | 20231215 | 256,489,776 |
| `ccf_coordinates` | 20231215 | 295,901,193 |
| `parcellation_to_parcellation_term_membership_acronym` | 20230630 | 22,261 |
| `gene` | 20241115 | 48,416 |
| 필수 합계 | | **1,571,075,446** |
| 선택 입력: measured raw H5AD | 20230830 | **7,627,589,574** |

Penk 실측을 포함하면 입력 합계 약 **9.20 GB**이며, 다운로드와 산출물에 추가 여유 공간이 필요하다. Imputed expression은 사용하지 않는다. 필요한 공개 호스트는 `allen-brain-cell-atlas.s3.us-west-2.amazonaws.com`이며 별도의 자격증명은 필요하지 않다.

## 실행

프로젝트 최상위에서 실행한다. 경로는 스크립트 위치에서 계산되므로 Windows의 `G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`에서도 동일한 프로젝트 구조를 사용할 수 있다.

```bash
# 현재 로컬 파일의 유무와 checksum을 검사; 입력 부재 시 exit code 2
.venv/bin/python src/merfish/run_merfish.py

# 환경에서 공식 S3 호스트가 허용된 후 실제 다운로드와 구성 분석
.venv/bin/python src/merfish/run_merfish.py --download

# 선택: 7.63 GB raw matrix를 추가 확보하고 measured Penk 기술 통계 산출
.venv/bin/python src/merfish/run_merfish.py --download --with-penk

# 독립된 임시 synthetic fixture로 코드 동작만 검증
.venv/bin/python src/merfish/test_merfish_synthetic.py
```

Windows에서는 위 Python 실행 파일을 `.venv\Scripts\python.exe`로 바꾸면 된다. `--manifest PATH`로 공식 release manifest를 명시할 수 있다. 원본 cache 위치는 manifest의 `relative_path` 구조를 유지한 `input/merfish/cache/`다. 수동으로 파일을 배치한 경우에도 동일한 공식 checksum을 통과해야 분석한다.

## 계산 규칙과 산출물

CCF `structure`가 **정확히** `MOp` 또는 `VISp`인 세포를 먼저 선택하고 원본 `substructure`에서 해부학적 layer를 가져온다. MOs, ALM, VIS, VIS-PTLp를 MOp 또는 VISp로 바꾸지 않는다. `L4/5 IT` 같은 transcriptomic label은 해부학적 L4로 바꾸지 않는다. layer 없는 세포는 `unknown`으로 보존한다.

원본 클래스에 분류되고 원본 subclass와 cluster alias가 있는 공개 세포 metadata를 QC annotated 분모로 사용한다. 새 correlation-score 임계치를 만들지 않는다. 원본 `average_correlation_score`도 산출물에 유지한다. annotation 또는 좌표 join 실패 수는 audit에 보고한다.

- `results/merfish/readiness_report.json`: 입력 유무·checksum·실제 분석 수행 여부·pending 이유·완료 시 cell audit.
- `results/merfish/network_acquisition_status.json`: 초기 proxy 차단 및 현재 HTTPS 다운로드 성공 기록.
- `results/merfish/actual_data_validation.json`: 실제 145,860개 세포와 32개 구성 CSV의 분모·합계·좌표·원본 검증 확인.
- `results/merfish/composition_report_ko.md`: 실제 세포 구성 결과와 단일 동물 해석 조건.
- `selected_cells_raw_labels_and_coordinates.csv.gz`: 영역·해부학적 층·원본 taxonomy·correlation score·CCF 좌표·표시 그룹·QC flag.
- `raw_label_display_group_mapping.csv`: raw label → 표시용 그룹과 taxonomy 버전.
- `sample_coverage.csv`: donor, genotype, sex, region, layer, hemisphere, section별 세포 수.
- `composition_*.csv`: 영역 전체, layer, hemisphere, section, layer×hemisphere, layer×section, layer×hemisphere×section 및 superficial/deep 보조 그룹. 전체 QC annotated 세포 대비, 뉴런 대비, GABAergic neuron 대비 **분모를 별도 열과 파일로 유지**한다. 존재하는 범주가 한 영역에서 관측되지 않으면 0 cell count를 표시한다.
- `regional_descriptive_difference_*.csv`: MOp−VISp percentage point 차이. 동물 수준 검정 또는 유의성 결과가 아니다.
- `figures/merfish/*.png`: 전체 큰 그룹 구성 bar plot, 해부학적 층별 heatmap, 원본 뉴런 subclass 비교.
- `candidate_panel_measured_coverage.csv`: 실제 확보한 gene list에서 22개 후보의 측정 여부. 미측정은 0 발현이 아니다.
- `measured_Penk_*.csv`: 선택 실행 시에만 생성. raw count>0 검출률, raw count 평균, `log1p(Penk / measured 500-gene panel total × 10000)` 평균. 500개 패널 내 정규화이며 RNA-seq의 전체 transcriptome 정규화와 직접 합치지 않는다. panel count가 0인 세포의 정규화 값은 결측으로 유지한다.

CCF 좌표는 mm이며 ASL 방향이다. x는 anterior→posterior, y는 dorsal→ventral, z는 left→right다. 원본 CCF 1140×10 μm 측방 범위의 중앙 `z=5.7 mm`를 반구 구분에 사용한다. 반구와 절편은 공간 표집 단위이며 독립 동물 반복이 아니다.

한 동물의 기술 통계만 보고한다. p value, 독립 동물 n으로 계산한 confidence interval, `cells/mm²` 또는 `cells/mm³`는 산출하지 않는다. 면적 마스크 없이 세포 밀도나 전체 영역 세포 수를 추정하지 않는다. 자극 AP/ML은 bregma 기준이고 DV는 뇌 표면 기준이므로 검증된 bregma–CCF 변환 및 광분포가 없으면 국소 ROI를 자동 생성하지 않는다.

## 검증 상태

2026-10-11에 독립적인 **synthetic fixture 4개**를 검증했다. 원본 taxonomy 표시 매핑, 세 종류의 분모와 영역별 0 범주, 정확한 영역 join·층·반구 및 그림 산출, raw Penk 검출·명시적 정규화와 zero-total 결측 처리를 확인했다. fixture 파일과 그림은 임시 디렉터리에서 제거되며 실제 `results/merfish`에 저장되지 않는다. 최초 입력 검사에서는 missing input 상태를 정확하게 `exit 2`, `actual_analysis_computed=false`로 보고했다. 현재 실제 입력을 확보한 실행은 `exit 0`, `status=composition_completed`, `actual_analysis_computed=true`이며, 실제 결과에 대한 독립적인 출력 검증도 통과했다. 총 32개 구성 CSV의 각 층화 단위에서 cell count의 합이 분모와 일치하고 fraction의 합이 1임을 확인했다.
