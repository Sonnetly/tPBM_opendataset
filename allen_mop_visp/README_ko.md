# Allen MOp–VISp 분석

목적은 정상 상태 atlas에서 MOp와 VISp의 세포 구성과 대응 세포 유형 안의 후보 유전자 발현을 비교하고, 기존 PBM 결과와 연결하는 것이다. 후속 회로 분석은 NREM 관련 영역으로의 투사를 우선 확인한다. 정상 상태 발현과 atlas의 축삭 신호를 PBM 처리 효과 또는 수면 유발 효과로 해석하지 않는다.

2026-10-12 후속 검토: [MOp L4·구성 그림의 분모와 표현 검토](docs/composition_visualization_review_ko.md), [실제 810 nm 연속광·180분 조건의 M1 PBM 수면 반응 가설](docs/M1_PBM_sleep_sensitivity_hypothesis_ko.md)을 추가했다. 기존 구성 그림은 유지했으며, 가설은 실제 지역별 PBM 효과를 계산한 결과가 아니다. 사용자 자극 조건과 새로 확인한 원 논문·공개 보충자료를 함께 기록했다.

## 저장 위치

현재 정본은 저장소의 `allen_mop_visp/` 폴더다. 클라우드 절대 경로는 `/workspace/tPBM_opendataset/allen_mop_visp`이다. 사용자의 커밋·푸시 요청에 따라 이전 `/workspace/Allen Dataset_M1_V1`에서 파일 68개와 inventory를 복사하고 SHA256으로 일치를 확인했다. 지정한 Windows 경로는 `G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`이며 클라우드에 연결되어 있지 않다. 로컬 저장이 완료됐다는 의미가 아니다.

`docs/`에 계획서와 자극 정보, `input/`에 출처와 입력, `src/`에 코드, `results/`에 계산 결과, `figures/`에 그림을 보관한다. 임시 네트워크 진단은 `logs/`에 두고 Git에서 제외한다. 검증 기록은 `results/validation/`과 각 분석의 결과 폴더에 보존한다. 세부 방침은 [파일 관리 안내](docs/file_management_ko.md)를 참고한다.

## 현재 완료와 미완료

- 기존 PBM 입력을 다시 확인하고 22개 패널·공통 DEG 52개·GO 38개·공통 GO hit 19개의 연결 표를 만들었다. 22개 패널과 공통 52개를 합친 탐색 목록은 62개이다. 원본 복사본·SHA256·검증 결과는 `input/pbm/`과 `results/pbm/`에 있다.
- 공식 유전자 목록으로 실측 패널을 확인했다. B6J 500-gene 패널은 핵심 후보 중 Penk만 측정하고, 탐색 62개 중 7개를 측정한다. B6N Zhuang-1 1122-gene 패널은 핵심 후보 중 Adora1·Adora2a·Penk·Trpc5를 측정한다. 미측정 유전자는 0 발현이 아니다. `results/coverage/`에 전체 표가 있다.
- 원본 다운로드 재시도가 성공했다. 공식 ABC manifest의 B6J MERFISH metadata 5개(1.57 GB)를 받아 크기와 공개 MD5를 검증했다. MOp 83,976개·VISp 61,884개 세포의 전체 영역·층·subclass 구성 비율을 실제 계산했다. [MERFISH 보고서](results/merfish/composition_report_ko.md)에 결과와 분모가 있다.
- RNA-seq은 공식 GEO GSE185862의 10x 원본(발현 행렬 5.72 GB 및 metadata)을 확보해 exact MOp 123,727개·VISp 29,763개 세포의 62개 후보 발현을 실제 계산했다. 원본 전체 31,053개 유전자의 세포별 count 합계를 정규화 분모로 보존했다. [실제 발현 보고서](results/rnaseq/RNAseq_results_ko.md), [RNA-seq 상태 파일](results/rnaseq/status.json), [재현 방법](docs/RNAseq_analysis_ko.md)에 결과와 조건이 있다.
- 공식 SMART-seq v4 원본도 받아 MOp 6,588개·VISp 15,572개와 62개 후보를 별도로 분석했다. 조건을 맞춘 자체 비교 38개 strata를 표시하며 10x와 직접 대응되는 strata는 없다. 다른 reporter·연령의 보조 결과를 10x 방향 재현으로 해석하지 않는다. `results/rnaseq/smartseq/`에 별도 결과가 있다.
- Allen Connectivity의 MOp·VISp primary injection 실험 168개를 실제 확보했다. WT MOp 8개·VISp 33개를 주 분석으로, Cre line MOp 23개·VISp 104개는 line별 별도 탐색으로 보존했다. NREM 관련 atlas label 7개의 좌·우·양측 3,528개 unionize 기록을 확인하고 연령·주입 순도 민감도를 계산했다. PZ는 ontology exact label이 없어 미평가로 남겼다. ZI에서 관찰된 차이는 후속 해부학 후보이며 Lhx6 세포 targeting·시냅스·NREM 효과는 미확인이다. `results/connectivity/`에 개별 실험·요약·상태·보고서·별도 Excel을 보관한다.
- 공간 주 자료는 B6J 단일 동물 MERFISH이며 동물 간 통계 검정을 수행하지 않는다. RNA-seq의 세포 수를 조직 구성 비율이나 독립 동물 수로 사용하지 않는다.
- WT B6J와 B6J 배경 transgenic은 분리한다. 정확한 MOp와 VISp가 있는 자료만 직접 비교하며 ALM이나 VIS 전체를 대신 사용하지 않는다.
- 국소 자극 ROI는 V1 AP가 잠정값이고 bregma→CCF 등록을 수행하지 않아 보류했다. mm·bregma AP/ML 및 뇌 표면 DV, core 직경 400 µm는 사용자 확인을 받았다. M1은 (1.34, 2.46, 0), V1은 (약 −2.3, 2.46, 0)이다. 사용자 좌표는 [기록 파일](docs/stimulation_coordinates_user.json)에 있다. 400 µm만으로 광학적 ROI를 정하지 않는다.

## 설치

Python 3.12를 사용한다. Linux/macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
```

Windows에서는 이 폴더에서 PowerShell로 `./setup_windows.ps1`를 실행하거나 아래 명령을 사용한다. PowerShell 설치 경로는 이 클라우드에서 실행 검증하지 않았다.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
```

## 분석 재개

모든 명령은 이 프로젝트 폴더에서 실행한다. Windows에서는 `.venv/bin/python` 대신 `.venv\Scripts\python.exe`를 사용한다.

이미 보존한 PBM 원본과 공식 gene list에서 완료된 연결 표·실측 audit를 다시 생성한다:

```bash
.venv/bin/python src/build_pbm_panels.py
.venv/bin/python src/audit_gene_coverage.py
```

B6J MERFISH의 공식 manifest와 체크섬으로 metadata를 받아 전체 영역·해부학적 층·subclass의 구성 비율을 계산한다:

```bash
.venv/bin/python src/merfish/run_merfish.py --download
```

공식 GEO에서 10x 입력을 내려받아 후보 subset과 전체 유전자 정규화 분모를 다시 생성한다. 대형 원본을 받지 않아도 아래의 발현 계산은 이동용 묶음에 포함된 subset으로 재현할 수 있다:

```bash
.venv/bin/python src/rnaseq/download_geo.py GSE185862_expression_matrix_10x.hdf5
.venv/bin/python src/rnaseq/download_geo.py GSE185862_metadata_10x.csv.gz
.venv/bin/python src/rnaseq/download_geo.py GSE185862_CTX-HPF-73538sample_metadata.xlsx
.venv/bin/python src/rnaseq/extract_geo_subset.py --workers 2
```

검증된 raw counts 후보 subset에서 같은 세포 유형 안의 기술 통계를 생성한다. 후보 62개의 합계를 정규화 분모로 사용하지 않는다:

```bash
.venv/bin/python src/rnaseq/analyze_expression.py --project-root . --h5ad input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad --metadata input/rnaseq/selected_cell_metadata.csv.gz --value-kind raw_counts --library-total-column full_library_counts
.venv/bin/python src/rnaseq/cohort_screening.py
.venv/bin/python src/rnaseq/summarize_findings.py
.venv/bin/python src/build_atlas_workbook.py
```

다른 입력의 메타데이터 스키마와 gene ID는 `--schema-map`, `--gene-map`으로 명시한다. 해당 10x 자료는 수컷 `Snap25-IRES2-Cre;Ai14` reporter이며 WT MERFISH와 구분한다. 지역별 donor·library date·일부 연령이 달라 지역 효과를 독립적으로 추정하는 p-value/FDR를 계산하지 않는다. 공유되는 P58 연령의 비교와 P57–59 탐색 비교의 범위를 구분한다. RNA-seq 세포 수를 조직 구성 비율로 쓰지 않는다.

NREM 관련 projection의 실제 공식 API 응답을 확보하고, 보존된 응답으로 표적별 기술 통계를 재현한다:

```bash
.venv/bin/python src/connectivity/run_connectivity.py --download
.venv/bin/python src/connectivity/build_connectivity_workbook.py
.venv/bin/python src/connectivity/validate_connectivity.py
```

이미 원본 응답이 있으면 `--download` 없이 같은 runner로 재현할 수 있다. 실제 관측 결과는 [투사 보고서](results/connectivity/NREM_projection_report_ko.md)와 [별도 Excel](results/connectivity/MOp_VISp_NREM_projection.xlsx)에 있다. 기능 문헌과 ontology 대응은 [NREM 표적 정의](docs/NREM_projection_target_definition_ko.md)에 있다. 시상 망상핵의 Allen label은 `RT`이며 뇌간의 `TRN`과 구분한다. 분자 경로의 전체 유전자 추출·모듈·GSEA 및 회로 기능 검증 방법은 [후속 분석 제안](docs/downstream_pathway_analysis_options_ko.md)에 정리했다. 새 분자 경로 점수는 아직 계산하지 않았다.

## 네트워크와 원본 보존

현재 런타임에서 ABC S3, 공식 GEO HTTPS archive와 Allen Connectivity HTTPS API에 실제 접근·다운로드했다. 과거 Allen CTX-HPF S3 경로는 여전히 막혀 있어 동일 연구의 공식 GEO 원본을 사용한다. 네트워크 설정 초안 저장은 런타임 적용을 뜻하지 않는다. 주소·응답·검증 정보는 `input/merfish/`, `input/rnaseq/`, `input/connectivity/`, [네트워크 기록](docs/network_requirements.json)에 보존한다.

대형 다운로드는 `input/merfish/cache/`와 `input/rnaseq/raw/`에 남기고 Git·이동용 ZIP에서 제외한다. 분석에 사용한 MOp/VISp metadata, 후보 raw-count subset, 전체 유전자 count 합계, 공식 출처·hash는 포함한다. 원본이 필요하면 같은 코드와 URL로 다시 받는다.

## 이동용 묶음

```bash
.venv/bin/python src/package_project.py
```

이 명령은 `exports/Allen_Dataset_M1_V1.zip`을 만들고, `file_inventory.json`에 각 파일의 SHA256을 기록한다. ZIP의 CRC 및 파일별 해시를 확인한다. OS별 가상환경·임시 로그·이전 ZIP·대형 원본 캐시는 제외한다. 재현용 subset·코드·결과·보고서·출처는 포함한다. 압축 안의 최상위 폴더 이름은 로컬 사용을 위해 `Allen Dataset_M1_V1`로 유지한다. Windows 지정 위치의 부모 폴더에 압축을 풀면 같은 프로젝트 구조를 사용할 수 있다.

이전 클라우드 작업 폴더에서 만든 최초 ZIP도 `exports/Allen_Dataset_M1_V1_initial.zip`으로 보존했다. 최신 재현 안내는 현재 README와 이름에 `initial`이 없는 ZIP을 사용한다.
