# tPBM_opendataset

[ORA 분석 보고서](ORA_analysis_report_ko.md)와 [전체 분석 파일·재실행 방법](ora_analysis/README.md)을 함께 보관한다. 코드, 기능 매핑, 원본 입력, CSV/Excel 결과는 `ora_analysis/`에 있다.

**최신 DEG FDR·log2FC 반영 결과:** [cortex qPCR 후보 재평가](ora_analysis/deg_validation/DEG_aware_qPCR_report_ko.md). 기존 ORA에서 복원한 공통 51개에 원본 DEG의 Tmem40 공통 down이 추가되어, 제공된 두 PBM 비교의 같은 방향 overlap은 52개다.

[qPCR 후보 요약 Markdown](qPCR_candidate_summary_ko.md)에는 1차 8개·미토콘드리아 확장 10개 패널과 참고 CSV 링크를 모았다.

## Allen MOp–VISp 분석

[M1–V1 데이터 해석 안내서와 요약 그림](allen_mop_visp/docs/M1_V1_data_analysis_guide_ko.md)은 실제 수행된 분석, PBM 반응에 연결할 가설, 추가 분석·실험의 우선순위와 한계를 한곳에 정리한다.

[분석 파일과 재실행 방법](allen_mop_visp/README_ko.md), [실제 분석 보고서](allen_mop_visp/analysis_report_ko.md), [통합 비교 Excel](allen_mop_visp/results/MOp_VISp_Allen_comparison.xlsx)은 `allen_mop_visp/`에 있다. 공식 원본 다운로드에 성공해 B6J MERFISH의 MOp 83,976개·VISp 61,884개 세포 구성과 CTX-HPF RNA-seq 153,490개 세포의 62개 후보 발현을 계산했다. 단일 WT 동물의 구성 통계와 batch가 분리된 reporter RNA-seq의 기술 통계를 구분하며, 동물 수준 유의성 검정은 하지 않았다.

B6N Zhuang-ABCA MERFISH(Zhang et al., Nature 2023, 1,122-gene panel) 중 Zhuang-ABCA-1~4 네 동물로 subclass 내 MOp−VISp 발현 차이, 동물 간 부호 일치, supertype 구성 표준화 점검을 계산했다. Adora1은 L2/3 IT와 L4/5 IT에서 MOp가 낮고 Pvalb에서 MOp가 높은 방향이 네 동물에서 일치했고 supertype 표준화 후에도 유지됐다. L6 CT의 차이는 대부분 supertype 구성 효과였다. B6J 주 분석과 합산하지 않으며 PBM 반응을 측정한 것이 아니다. 데이터셋 설명, 수행 내용, 결과는 [Zhuang 보고서](allen_mop_visp/results/zhuang/zhuang_report_ko.md)에 있다.

[M1·V1 뉴런 유형별 발현 보고서](allen_mop_visp/results/neuron_type_report/M1_V1_neuron_type_expression_report_ko.md)는 PBM 다음 실험 선택을 위해 두 자료를 역할별로 나눠 썼다. 10x 전사체 전체로 영역 안 뉴런 유형 marker와 OXPHOS·adenosine 세트를, Zhuang MERFISH 네 동물로 같은 뉴런 유형의 MOp 대 VISp 짝지은 DEG를 계산했다. 보고서용 그림 7종은 `allen_mop_visp/figures/report_m1v1/`에 있다.

Allen Mouse Brain Connectivity Atlas의 실제 MOp·VISp 주입 실험에서 NREM 관련 영역으로의 투사도 계산했다. WT 실험 MOp 8개·VISp 33개의 개별 지표와 연령·주입 순도 민감도 결과는 [투사 보고서](allen_mop_visp/results/connectivity/NREM_projection_report_ko.md)와 [별도 Excel](allen_mop_visp/results/connectivity/MOp_VISp_NREM_projection.xlsx)에 있다. ZI는 후속 회로 검증 후보이며, 축삭 신호를 수면촉진 세포와의 연결이나 NREM 유발 효과로 해석하지 않는다. [표적 정의와 기능 문헌](allen_mop_visp/docs/NREM_projection_target_definition_ko.md), [투사·분자 경로 후속 방법](allen_mop_visp/docs/downstream_pathway_analysis_options_ko.md)을 함께 참고한다.

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
