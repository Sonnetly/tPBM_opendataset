# 파일 관리와 로컬 이동

사용자가 지정한 최종 Windows 작업 위치:

`G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`

현재 정본은 Linux 클라우드의 `/workspace/tPBM_opendataset/allen_mop_visp`다. 이전 `/workspace/Allen Dataset_M1_V1`의 파일을 사용자의 커밋·푸시 요청에 따라 복사하고 SHA256으로 대조했다. Windows G: 드라이브는 이 머신에 연결되지 않았다. 클라우드에 저장한 파일을 Windows 경로에 저장했다고 간주하지 않는다.

배포용 프로젝트 폴더를 같은 이름의 로컬 폴더로 옮기면 상대경로 기반 코드를 사용할 수 있다. 원본 계획서는 docs, 코드는 src, 출처·분석 입력은 input, 계산 결과는 results, 그림은 figures, 임시 진단은 logs로 분리한다. .venv는 OS별로 새로 만들며 복사하지 않는다. 원본 대용량 캐시는 일반 프로젝트 패키지에 포함되지 않으므로, 해당 캐시까지 필요한 로컬 재실행에서는 공식 원본을 별도로 재다운로드하거나 아래 캐시 위치를 그대로 유지해 별도 복사한다.

이번 요청으로 분석 코드·입력·매핑·표·그림·보고서를 저장소의 `allen_mop_visp/`로 옮겼다. 기존 ORA 결과와 구분해 보존하며, 이후 해당 폴더를 정본으로 사용한다. `docs/migration_provenance.json`에 최초 복사한 파일의 원본 해시를, `file_inventory.json`에 현재 패키지 파일의 해시를 기록한다. 이동하면서 저장 위치·재현 안내와 패키징 경로만 갱신했다.

현재 공식 Allen ABC와 GEO HTTPS 원본 다운로드를 실제로 완료했다. 초기 네트워크 차단 기록은 이력을 위해 남겼으며, 현재 주 분석 자료의 성공 상태는 `docs/network_requirements.json`과 각 acquisition/provenance 파일에 기록한다. 일부 legacy Allen 링크의 두 호스트는 아직 차단되어 있지만 공식 ABC 및 GEO를 통한 주 분석 원본 확보는 완료했다.

| 보관 자료 | 프로젝트 내 경로 | Git·일반 패키지 | 역할 |
|---|---|---|---|
| MERFISH 전체 metadata 5개, 합계 1,571,075,446 bytes | `input/merfish/cache/` | 제외 | 원본 CSV를 공식 relative path와 함께 보관; 각 파일의 공식 MD5·크기 검증 |
| GEO 10x 전체 원본 HDF5, 5,716,426,072 bytes 및 전체 metadata | `input/rnaseq/raw/` | 제외 | 전체 RNA-seq 재추출용 원본; URL·크기·SHA256 기록 |
| GEO SMART-seq v4 전체 원본 HDF5, 3,585,764,160 bytes 및 metadata | `input/rnaseq/raw/` | 제외 | 별도 exon read-count 보조 자료; 10x와 합산하지 않음 |
| MERFISH 선택 세포 annotation·좌표 및 구성 결과 | `results/merfish/` | 포함 | 실제 분석한 145,860개 세포와 분모·표·그림의 재현 근거 |
| RNA-seq 선택 원본 metadata·donor 기록 | `input/rnaseq/selected_*` | 포함 | 선택 세포의 원본 annotation과 biological pool·batch 조건 보존 |
| RNA-seq 후보 62개 count 부분집합·전체 library total | `input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad`, `input/rnaseq/selected_cell_metadata.csv.gz` | 포함 | 22개 screening panel과 PBM 공통 52개의 합집합; 전체 31,053개 원본 유전자의 count 합계를 정규화 분모로 별도 보존 |
| SMART-seq 후보 count 부분집합·metadata·전체 exon totals | `input/rnaseq/smartseq/` | 포함 | 원본 45,768개 feature의 합계를 분모로 보존; 별도 genotype·sex·age·FACS 비교 재현 |
| Connectivity 공식 API 원본·ontology·표적 정의·출처·hash | `input/connectivity/` | 포함 | 대형 voxel 자료 없이 실제 실험 metadata와 표적별 unionize 분석 재현 |
| NREM 관련 투사 개별 실험·조건별 요약·보고서·Excel·그림 | `results/connectivity/`, `figures/connectivity/` | 포함 | WT·Cre line, 연령·주입 순도 민감도, PZ 미평가와 기능 해석 한계 보존 |
| 원본 URL·해시·출처·재현 코드·보고서·표·그림 | `input/`, `src/`, `docs/`, `results/`, `figures/`의 캐시 밖 파일 | 포함 | 대형 전체 원본 없이도 선택한 분석 범위와 결과를 검토하고 재현 |

RNA-seq 전체 matrix의 SHA256은 다운로드된 파일의 로컬 해시이며, 해당 matrix에는 공개 publisher checksum이 없었다. 이를 원본 제공자가 발표한 checksum 검증으로 표시하지 않는다. TLS 검증을 유지하고 공식 GEO 출처, 응답 크기, SHA256을 함께 기록했다. MERFISH 5개 파일은 publisher MD5와 크기를 모두 대조했다.

대형 전체 원본은 Git ignore 및 `src/package_project.py`의 패키지 제외 규칙으로 관리한다. 작은 선택 자료와 원본 해시는 제외된 raw/cache 디렉터리 밖에 보존한다. RNA-seq 62개 유전자 부분집합으로 정규화할 때도 62개만 합산한 값을 library size로 쓰지 않고, 전체 원본 count 합계를 사용한다. 실제 실행·산출 여부는 `results/rnaseq/status.json` 및 `input/rnaseq/subset_manifest.json`에서 확인한다.

대형 원본을 새로 확보하거나 재실행하기 전에 현재 디스크 여유와 필요한 입력 크기를 확인한다. 기본 MERFISH 구성에는 약 1.57 GB, 10x 전체 matrix에는 약 5.72 GB, SMART-seq 보조 matrix에는 약 3.59 GB와 metadata·산출물 여유가 필요하다. 선택 사항인 Penk 실측 raw MERFISH matrix 약 7.63 GB는 현재 확보하지 않았다. 원본 캐시가 있는 클라우드에서 재추출은 가능하며, 로컬 패키지에는 선택 자료·출처·재다운로드 코드를 보존한다.

사용자의 커밋·푸시 요청에는 캐시·가상환경·실제 자격증명을 제외한 누적 분석 산출물 전체를 반영한다. 정확한 커밋과 원격 반영 상태는 Git 이력과 push 뒤 원격 검증으로 확인한다. 분석 상태 파일과 이동용 묶음은 과학 결과를 기록하며 특정 시점의 Git 상태를 대신하지 않는다.
