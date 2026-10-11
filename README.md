# tPBM_opendataset

[ORA 분석 보고서](ORA_analysis_report_ko.md)와 [전체 분석 파일·재실행 방법](ora_analysis/README.md)을 함께 보관한다. 코드, 기능 매핑, 원본 입력, CSV/Excel 결과는 `ora_analysis/`에 있다.

**최신 DEG FDR·log2FC 반영 결과:** [cortex qPCR 후보 재평가](ora_analysis/deg_validation/DEG_aware_qPCR_report_ko.md). 기존 ORA에서 복원한 공통 51개에 원본 DEG의 Tmem40 공통 down이 추가되어, 제공된 두 PBM 비교의 같은 방향 overlap은 52개다.

[qPCR 후보 요약 Markdown](qPCR_candidate_summary_ko.md)에는 1차 8개·미토콘드리아 확장 10개 패널과 참고 CSV 링크를 모았다.

## Allen MOp–VISp 분석

[분석 파일과 재실행 방법](allen_mop_visp/README_ko.md), [진행 보고서](allen_mop_visp/analysis_report_ko.md), [PBM 후보 연결 Excel](allen_mop_visp/results/pbm/PBM_panel22_and_gene_sets.xlsx)은 `allen_mop_visp/`에 있다. 기존 PBM 후보 연결과 공식 실측 패널 확인을 완료했으며, 실제 Allen 지역별 세포 구성·발현 비교는 원본 다운로드 접근 문제로 보류 중이다.

## GitHub 반영

Codex에게 **“커밋하고 푸쉬해 줘”**라고 요청하면 `AGENTS.md`의 규칙에 따라 현재까지 누적된 프로젝트 변경 전체를 반영한다. 파일을 생성할 때부터 저장소 안에 저장하도록 같은 지침에 기록했다.

직접 실행할 때는 저장소 루트에서 다음 명령을 사용한다.

```bash
bash scripts/commit_and_push.sh "분석 결과 업데이트"
```

이 명령은 새 파일·수정·삭제를 모두 stage하고, 변경이 있으면 커밋한 뒤 `origin/main`에 push한다. 작업 파일 변경 없이 기존 로컬 커밋만 남아 있어도 push한다. 메시지, 원격, 브랜치를 각각 지정할 수 있다.

```bash
bash scripts/commit_and_push.sh "분석 결과 업데이트" origin main
```

원격 브랜치의 이력이 현재 로컬 커밋에 포함되어 있는지 먼저 확인한다. 원격이 앞서거나 이력이 갈라졌다면 파일을 stage하거나 커밋하기 전에 중단하므로, 차이를 해결한 뒤 다시 실행한다. 강제 push는 사용하지 않는다. Environment Publish는 GitHub 반영과 별도다.

자동화 명령은 외부 네트워크 없이 로컬 Git 원격을 만들어 검증할 수 있다.

```bash
python3 -m unittest discover -s tests -v
```
