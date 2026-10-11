# 기존 PBM 결과와 Allen 후보 패널 연결

기준 저장소 commit: `e430e8d2d2827ccdd222e8828d77227e937b6f1b`. 원본 입력 13개를 `input/pbm/repository/`에 보존하고 SHA-256과 Git commit 내용을 검증했다.

같은 방향 공통 DEG는 **52개(up 34/down 18)**, 반대 방향은 23개다. 두 up ORA에서 모두 유의한 GO는 **38개(BP29/MF9)**이고 실제 공통 hit 합집합은 **19개**다. GO 1개에는 실제 공통 hit가 없다. 전체 DEG를 사용하면 기존 ORA hit 기반 51개에 Tmem40 공통 down이 추가된다.

22개 screening panel은 아데노신 8개, 기존 PBM 1차 8개, 전기생리 4개, 미토콘드리아 2개다. 공통 52개와 합친 탐색 범위는 62개다. 이 패널은 Atlas 기저 발현을 확인하기 위한 후보 목록이며 아직 MOp–VISp 차이가 확인된 결과나 최종 qPCR 순위는 아니다.

| Gene | 피질 PBM log2FC | 원본 gene FDR |
|---|---:|---:|
| Nt5e | -0.5367 | 0.0074455603 |
| Adora2a | -1.2929 | 0.0093214593 |

Entpd1, Slc29a1, Slc29a2, Adk, Ada, Adora1은 두 비교의 제공 FDR05 표에 없다. 이는 미발현·변화 없음·유의하지 않음의 확정 결과가 아니다. `missing_from_FDR05_table`과 빈 수치로 보존하고 0으로 대체하지 않았다. 전체 tested-gene 통계가 있어야 비선택 결과를 확인할 수 있다.

해마는 5xFAD + light vs 5xFAD(808 nm), 피질은 hypobaric hypoxia + PBM vs hypobaric hypoxia(1064 nm)이다. 정상 B6의 MOp–VISp 기저 차이 또는 지역별 PBM 효과로 직접 해석하지 않는다. Atlas 발현이 높다는 사실만으로 PBM 반응이 크거나 adenosine 순축적이 많다고 판정할 수 없다.

원본의 gene FDR/log2FC와 GO 항목 FDR를 각각 보존했다. FDR05 subset에서 BH를 다시 계산하거나 새로운 ORA를 실행하지 않았다. 중복 심볼의 최소 FDR/최대 효과를 골라 쓰지 않고 개별 행을 audit에 보존했다. 공유 52개와 panel22에는 중복 심볼의 수치 모호성이 없다.

재계산한 52개와 GO38의 심볼·통계·hit가 저장소 기존 결과와 일치했고, ORA 입력 개수와 모든 hit의 DEG 방향을 확인했다. 계획서의 Slitrk3 피질 FDR `0.01…`은 반올림 표기이며 결과 CSV는 원본 정밀도를 사용한다.

재실행: `python src/build_pbm_panels.py` (pandas, numpy, openpyxl 필요). 첫 입력 보존만 `--source-repo /workspace/tPBM_opendataset`을 지정한다. 이 분석 명령은 보존된 원본을 읽고 결과 파일을 생성하며 Git 커밋·푸시는 별도 작업이다.

- `PBM_panel22_and_gene_sets.xlsx`: 주요 패널/공통 유전자/GO/검증 결과. 저장 후 모든 sheet 행 수와 gene 순서를 다시 읽어 검증했다.
- `panel22.csv`: gene별 가설 축·역할·해마/피질 방향/log2FC/FDR/status.
- `same_direction_52.csv`, `common_up34.csv`, `common_down18.csv`: 탐색 공통 목록.
- `common_up_GO_38.csv`, `common_GO_hit19.csv`, `common_GO_gene_pairs.csv`: GO와 실제 hit, 서로 분리한 항목/gene FDR.
- `exploration_genes62.csv`, `gene_group_manifest.json`: Allen 입력용 통합 유전자와 그룹.
- `duplicate_placeholder_audit.csv`, `ORA_DEG_alignment.csv`, `validation_summary.json`: 원본 audit 및 검증 기록.
