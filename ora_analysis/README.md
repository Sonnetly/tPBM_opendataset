# ORA 유전자 방향 비교 분석

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
