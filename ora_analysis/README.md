# ORA 유전자 방향 비교 분석

**최신 결과는 [DEG FDR·log2FC 반영 보고서](deg_validation/DEG_aware_qPCR_report_ko.md)를 먼저 참고한다.** 아래 ORA-only 결과의 51개는 hits에서 복원한 범위이며, 제공된 원본 DEG 표에서는 같은 방향 overlap이 52개(up34/down18)로 확인된다. 이전 qPCR 제안에 gene-level 수치를 연결한 최신 우선 패널도 해당 보고서에 있다.

[별도 qPCR 후보 요약](../qPCR_candidate_summary_ko.md)에서 패널과 참고 CSV를 바로 확인할 수 있다.

현재 목표는 B6 cortex tPBM 실험에서 방향성을 qPCR로 검증하는 것이다. 1차 공통 방향 8개는 Grm2/Ryr2/Trpc5/Atp2b3/Cacna1c/Grin2a/Mrps6/Penk이며, 피질에서만 확인되는 미토콘드리아 후보 Bcs1l/Ppargc1a를 추가하면 10개 패널이다. FDR·log2FC는 검증 순서의 참고값이며, 질환 맥락을 후보 제외 조건으로 쓰지 않는다.

첨부된 `ORA_fullcatalog_KEGG_GO.csv`의 전체 ORA 행에서 유전자별 up/down을 복원하고, 두 데이터셋의 같은 방향·반대 방향 overlap과 KEGG/GO 근거를 정리한다.

이 분석은 `GSE245658_hippocampus_808nm`와 `GSE325653_prefrontal_1064nm` 및 현재 기능 매핑에 맞춰 작성되었다. 다른 데이터에 적용하려면 분석의 검증 조건과 `gene_function_map.tsv`를 함께 검토해야 한다.

## 파일

- `input/ORA_fullcatalog_KEGG_GO.csv`: 제공받은 원본 입력. 원본 파일 SHA256은 `validation_summary.json`에 기록되어 있다.
- `analyze_ora.py`: 분석·내부 수치 검증·CSV/Excel/Markdown 출력 코드.
- `gene_function_map.tsv`: 51개 공통 유전자의 설명용 기능 분류와 대표 KEGG/GO ID. 데이터베이스 전체 annotation이나 통계적 enrichment 결과가 아니다.
- `requirements.txt`: 이번 실행에서 확인한 Python 패키지 버전.
- `ORA_analysis_report_ko.md`: 유전자별 근거를 포함한 한국어 보고서.
- `ORA_direction_overlap_analysis.xlsx`: 전체 복원 유전자, 공통·반대방향 유전자, 근거, 원본 ORA 행을 담은 Excel.
- 결과 CSV와 `validation_summary.json`: 표와 수치 검증 결과.

## 재실행

Python 3.12 환경에서 검증했다. 아래 명령은 **저장소 루트**에서 실행한다.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r ora_analysis/requirements.txt
.venv/bin/python ora_analysis/analyze_ora.py \
  ora_analysis/input/ORA_fullcatalog_KEGG_GO.csv ora_analysis
cp ora_analysis/ORA_analysis_report_ko.md ORA_analysis_report_ko.md
```

출력 폴더를 다른 곳으로 지정하면 기존 결과를 덮어쓰지 않고 재현성을 확인할 수 있다.

```bash
.venv/bin/python ora_analysis/analyze_ora.py \
  ora_analysis/input/ORA_fullcatalog_KEGG_GO.csv /tmp/ora_reproduction
```

현재 입력의 결과는 복원 유전자 1,366개, 같은 방향 51개(up 34/down 17), 반대 방향 23개다. 검증에 실패하면 프로그램은 성공 결과를 출력하지 않고 종료한다.

원본 DEG 수와 ORA hits에서 복원되는 수는 다르다. 완전한 DEG 기준 overlap, gene-level FDR와 효과 크기를 판단하려면 각 데이터셋의 원본 DEG 표가 필요하다. 기존 ORA FDR는 pathway/GO 항목의 값이며, overlap 51개를 대상으로 새로 검정한 값이 아니다. 공식 항목 링크는 제공하지만 외부 자료는 이번 환경에서 직접 조회하지 못했다.

## Cortex qPCR 후보와 공통 up GO

사용자의 가설은 피질 영역별 tPBM 수면 반응의 차이가 국소 피질 반응성 차이에서 생긴다는 것이다. 이온/칼슘 신호와 미토콘드리아 관련 표적을 중심으로 잠정적인 실험 후보를 선정했다.

```bash
python3 ora_analysis/prioritize_qpcr_and_go.py
```

같은 Python 의존성을 사용한다. 후보와 근거 설정은 `qpcr_candidate_config.tsv`, 결과는 `qpcr_go/`에 저장한다. 기존 ORA 분석을 변경하지 않는다.

- [qPCR 우선 후보와 공통 up GO 보고서](qpcr_go/qPCR_priority_and_common_up_GO_report_ko.md)
- [Excel 전체 표](qpcr_go/qPCR_priority_and_common_up_GO.xlsx)
- `qpcr_go/qpcr_candidates.csv`: cortex 가설 기준 후보 27개. 처음 8개가 우선 제안이며 통계적 gene-level 순위가 아니다.
- `qpcr_go/common_up_GO_significant.csv`: 두 데이터셋 공통 up에서 기존 ORA FDR가 모두 0.05 미만인 GO 38개. BP 29개, MF 9개를 별도 CSV로도 저장한다.
- `qpcr_go/primary_candidate_evidence.csv`: 후보별 원본 행, 관찰 방향, 주석과 기존 ORA 통계. cortex에서만 관찰되는 후보와 공통 방향 후보를 구분한다.
- `qpcr_go/cortex_mito_related_ORA.csv`: 지정한 항목명 패턴으로 추출한 전전두피질 미토콘드리아 관련 ORA 부분집합.

공통 up 유의 GO의 실제 공통 hit는 19개 유전자다. 38개 GO 중 하나는 두 연구 모두 유의하지만 hit 교집합이 비어 있다. 미토콘드리아 우선 후보는 가설 검증용이며 유의한 미토 경로 enrichment를 뜻하지 않는다. 연구 간 조직과 파장이 함께 달라, 영역별 반응성은 별도의 region×treatment 설계로 검증해야 한다.

기존 결과를 유지하며 재현성을 확인하려면 출력 위치를 따로 지정한다.

```bash
python3 ora_analysis/prioritize_qpcr_and_go.py --output-dir /tmp/qpcr_go_reproduction
```

## 원본 DEG FDR·log2FC 검증과 qPCR 우선순위 갱신

입력은 `input/DEG_FDR05_GSE245658_GSE325653.csv`다. 파일에 여러 조직·질병/처리 비교가 있으므로 해마의 `PBM effect: 5xFAD+light vs 5xFAD`와 전전두피질의 `PBM effect: hypobaric hypoxia + PBM vs hypobaric hypoxia`만 대응시킨다. 파일에 실린 분석 방법과 비교 조건은 외부 GEO에서 독립적으로 재검증한 것이 아니다.

```bash
python3 ora_analysis/analyze_deg_and_reprioritize.py
```

같은 Python 의존성을 사용하고 `deg_validation/`에 출력한다. 수동 우선순위와 이유는 `deg_aware_priority_config.tsv`에 기록했다. 원본 파일은 그대로 보존한다. 이미 FDR로 선택된 표에서 gene FDR를 다시 계산하거나 새로운 enrichment를 수행하지 않는다.

- [DEG-aware qPCR/overlap Excel](deg_validation/DEG_aware_qPCR_and_overlap.xlsx)
- [재평가 보고서](deg_validation/DEG_aware_qPCR_report_ko.md)
- `same_direction_52.csv`: 같은 방향 52개와 두 비교의 log2FC·gene FDR.
- `qPCR_candidates_with_DEG.csv`: 이전 후보와 추가 Tmem40, 원본 gene 통계 및 현재 우선순위.
- `common_up_GO_with_gene_stats.csv`: 기존 공통 up 유의 GO의 공통 hit 19개에 gene-level 수치를 연결한 표. GO 항목 FDR와 gene FDR를 별도 열로 보관한다.
- `ORA_DEG_alignment.csv`: 기존 ORA 입력 수와 visible hit 방향이 해당 PBM DEG와 일치하는지 검사한 결과.
- `duplicate_symbols_audit.csv`, `placeholder_symbols_audit.csv`: 중복 및 미지정 심볼을 보존한 audit. 유리한 FDR/log2FC를 골라 중복을 합치지 않는다.

해마 PBM 결과는 657행(유효 고유 심볼 651개)이며 pathway-analysis 사용 표시가 True인 646행이 기존 ORA 입력 수와 일치한다. 피질은 998행/998개 심볼이다. 미토콘드리아 후보의 gene FDR가 유의해도 미토 pathway/GO enrichment가 유의하거나 ATP 생산이 바뀌었다는 뜻은 아니다. 피질 lfcSE는 없고 baseMean의 연구 간 동일 척도도 확인되지 않아, 표본 수/검출력이나 절대 발현량 비교를 계산하지 않는다.

```bash
python3 ora_analysis/analyze_deg_and_reprioritize.py --output-dir /tmp/deg_qpcr_reproduction
```
