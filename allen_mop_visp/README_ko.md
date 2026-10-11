# Allen MOp–VISp 분석

목적은 정상 상태 atlas에서 MOp와 VISp의 세포 구성과 대응 세포 유형 안의 후보 유전자 발현을 비교하고, 기존 PBM 결과와 연결하는 것이다. 정상 상태 발현을 PBM 처리 효과 또는 수면 threshold로 해석하지 않는다.

## 저장 위치

현재 정본은 저장소의 `allen_mop_visp/` 폴더다. 클라우드 절대 경로는 `/workspace/tPBM_opendataset/allen_mop_visp`이다. 사용자의 커밋·푸시 요청에 따라 이전 `/workspace/Allen Dataset_M1_V1`에서 파일 68개와 inventory를 복사하고 SHA256으로 일치를 확인했다. 지정한 Windows 경로는 `G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`이며 클라우드에 연결되어 있지 않다. 로컬 저장이 완료됐다는 의미가 아니다.

`docs/`에 계획서와 자극 정보, `input/`에 출처와 입력, `src/`에 코드, `results/`에 계산 결과, `figures/`에 그림을 보관한다. 임시 네트워크 진단은 `logs/`에 두고 Git에서 제외한다. 검증 기록은 `results/validation/`과 각 분석의 결과 폴더에 보존한다. 세부 방침은 [파일 관리 안내](docs/file_management_ko.md)를 참고한다.

## 현재 완료와 미완료

- 기존 PBM 입력을 다시 확인하고 22개 패널·공통 DEG 52개·GO 38개·공통 GO hit 19개의 연결 표를 만들었다. 22개 패널과 공통 52개를 합친 탐색 목록은 62개이다. 원본 복사본·SHA256·검증 결과는 `input/pbm/`과 `results/pbm/`에 있다.
- 공식 유전자 목록으로 실측 패널을 확인했다. B6J 500-gene 패널은 핵심 후보 중 Penk만 측정하고, 탐색 62개 중 7개를 측정한다. B6N Zhuang-1 1122-gene 패널은 핵심 후보 중 Adora1·Adora2a·Penk·Trpc5를 측정한다. 미측정 유전자는 0 발현이 아니다. `results/coverage/`에 전체 표가 있다.
- Allen 원본 다운로드는 현재 클라우드 egress proxy의 CONNECT 403으로 막혔다. Allen 세포 비율·발현 차이·후보의 지역별 순위는 아직 계산하지 않았다.
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

RNA-seq 공식 다운로드 경로를 확인하고 MOp/VISp cell metadata의 생물학적 표본 조건부터 점검한다:

```bash
.venv/bin/python src/rnaseq/discover_downloads.py --project-root .
.venv/bin/python src/rnaseq/analyze_expression.py --project-root . --audit-only --metadata input/rnaseq/cell_metadata.csv
```

발현 행렬을 확보한 뒤 원본의 단위와 feature ID를 확인하고 명시적으로 입력한다. 아래 명령은 정규화에 필요한 전체 유전자 집합을 포함하는 `raw_counts` H5AD임을 원본 문서로 확인한 경우에만 적용한다. target-only 또는 HVG 행렬에서는 아래 확인 플래그를 쓰지 말고, 전체 행렬에서 계산한 세포별 library total을 metadata에 넣어 `--library-total-column`으로 지정한다. 다른 형식은 코드의 도움말과 원본 자료 설명에 맞춰 지정한다:

```bash
.venv/bin/python src/rnaseq/analyze_expression.py --project-root . --h5ad input/rnaseq/expression.h5ad --metadata input/rnaseq/cell_metadata.csv --value-kind raw_counts --full-transcriptome-counts-confirmed
```

메타데이터 스키마 매핑과 gene ID 매핑이 필요한 경우 각각 `--schema-map`, `--gene-map`을 사용한다. donor/pool과 batch를 확인하기 전에는 개체 수준 p-value를 생성하지 않는다. 실제 데이터 없이 synthetic fixture로 검증한 코드는 준비 상태만 의미하며 분석 완료를 의미하지 않는다.

## 현재 외부 선행조건

네트워크 설정 초안에 Allen 공식 S3·포털·문서와 NCBI 도메인 추가를 저장했다. 설정 저장만으로 현재 런타임 접근이 바뀌었다고 간주하지 않는다. 환경 설정에서 네트워크 변경을 검토·저장하고 제품 안내에 따라 게시한 뒤, 실제 자료 접근을 다시 확인해야 한다. Windows 로컬에서 다운로드하거나 이미 보유한 입력을 제공하는 방법도 가능하다. 공개 Allen 자료 다운로드에는 새 비밀키를 요구하지 않는다.

## 이동용 묶음

```bash
.venv/bin/python src/package_project.py
```

이 명령은 `exports/Allen_Dataset_M1_V1.zip`을 만들고, `file_inventory.json`에 각 파일의 SHA256을 기록한다. ZIP의 CRC 및 파일별 해시를 확인한다. OS별 가상환경·임시 로그·이전 ZIP은 제외한다. 원본 데이터·코드·보고서·출처는 포함한다. 압축 안의 최상위 폴더 이름은 로컬 사용을 위해 `Allen Dataset_M1_V1`로 유지한다. Windows 지정 위치의 부모 폴더에 압축을 풀면 같은 프로젝트 구조를 사용할 수 있다.

이전 클라우드 작업 폴더에서 만든 최초 ZIP도 `exports/Allen_Dataset_M1_V1_initial.zip`으로 보존했다. 최신 재현 안내는 현재 README와 이름에 `initial`이 없는 ZIP을 사용한다.
