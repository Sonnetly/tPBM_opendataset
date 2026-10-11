# 파일 관리와 로컬 이동

사용자가 지정한 최종 Windows 작업 위치:

`G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`

현재 정본은 Linux 클라우드의 `/workspace/tPBM_opendataset/allen_mop_visp`다. 이전 `/workspace/Allen Dataset_M1_V1`의 파일을 사용자의 커밋·푸시 요청에 따라 복사하고 SHA256으로 대조했다. Windows G: 드라이브는 이 머신에 연결되지 않았다. 클라우드에 저장한 파일을 Windows 경로에 저장했다고 간주하지 않는다.

폴더 전체를 같은 이름의 로컬 폴더로 옮기면 상대경로 기반 코드를 사용할 수 있다. 원본 계획서는 docs, 코드는 src, 출처·분석 입력은 input, 계산 결과는 results, 그림은 figures, 임시 진단은 logs로 분리한다. .venv는 OS별로 새로 만들며 복사하지 않는다.

이번 요청으로 분석 코드·입력·매핑·표·그림·보고서를 저장소의 `allen_mop_visp/`로 옮겼다. 기존 ORA 결과와 구분해 보존하며, 이후 해당 폴더를 정본으로 사용한다. `docs/migration_provenance.json`에 최초 복사한 파일의 원본 해시를, `file_inventory.json`에 현재 패키지 파일의 해시를 기록한다. 이동하면서 저장 위치·재현 안내와 패키징 경로만 갱신했다.

대형 전체 atlas 다운로드가 필요하면 용량을 먼저 확인한다. GitHub 단일 파일 크기 제한을 넘는 원본은 분석에 사용한 부분집합과 원본 URL·해시·재다운로드 절차를 함께 보존하고, 전체 원본의 보관 위치는 추후 GitHub 반영 시 사용자와 정한다. 캐시, 가상환경, 실제 자격증명은 반영하지 않는다.
