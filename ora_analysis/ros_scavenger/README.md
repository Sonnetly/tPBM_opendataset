# ROS scavenger 분석

**전체 KEGG와 비교하는 데이터셋별 후속 결과:** [보고서](../pathway_overview/Dataset_pathway_report_ko.md), [그림 10종](../pathway_overview/README.md). 이 폴더의 GSEA FDR는 ROS 한정 20검정의 값이며, 새 분석의 전체 경로 FDR와 다르다.

[한국어 분석 보고서](ROS_scavenger_report_ko.md) · [통합 Excel](results/ROS_scavenger_analysis.xlsx)

기존 제공 DEG 표를 이용한 ROS 유전자·공통 방향 분석과 원 GEO 전체 통계를 이용한 경로 검정을 구분한다. 분석은 `ROS_scavenger`에서 수행했고, GitHub 게시 브랜치는 `ROS_scavenger_20261103`이다.

## 그림

모든 그림은 `figures/`에 PNG·PDF·SVG로 있다. 실제 Arial을 사용하며 PDF에는 TrueType을 임베딩하고 SVG는 글리프를 path로 저장한다. 원 글꼴 파일은 저장소에 재배포하지 않는다.

| 그림 | 내용 |
|---|---|
| [01](figures/fig01_hippocampus_ROS_DEG.png) | 기존 해마 ROS DEG 23개 |
| [02](figures/fig02_cortex_ROS_DEG.png) | 기존 전전두피질 ROS DEG 46개 |
| [03](figures/fig03_shared_ROS_DEG.png) | 공통 방향 Penk 및 반대 방향 Col1a1·Foxo1 |
| [04](figures/fig04_antioxidant_DEG.png) | 기존 표의 antioxidant/removal DEG |
| [05](figures/fig05_GEO_scavenger_panel.png) | 원 GEO의 canonical scavenger·support 패널, 비유의 유전자 포함 |
| [06](figures/fig06_ROS_GSEA.png) | 탐색적 gene-label GSEA |
| [07](figures/fig07_existing_ORA_FDR.png) | 기존 full-catalog ORA FDR |
| [08](figures/fig08_sample_level_pathway_scores.png) | 표본 단위 경로 score와 정확 순열검정 |

## 재실행

저장소 루트에서 실행한다. 입력 snapshot은 `ora_analysis/input/ros_geo/`에, 산출물은 이 폴더 안에 저장한다. 기존 ORA/qPCR 분석을 덮어쓰지 않는다. 전체 source 파일 및 결과 CSV·Excel은 재현성을 위해 보관한다. 원 source를 최신으로 다시 받으려면 기존 snapshot을 별도 백업한 후 다운로드해야 한다. 기본 fetch는 저장된 파일을 재사용한다.

```bash
python3 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -r ora_analysis/ros_scavenger/requirements.txt
.venv/bin/python ora_analysis/ros_scavenger/fetch_sources.py
.venv/bin/python ora_analysis/ros_scavenger/setup_arial.py --directory /tmp/ros_arial
.venv/bin/python ora_analysis/ros_scavenger/geo_statistics.py
.venv/bin/python ora_analysis/ros_scavenger/analyze_ros.py
.venv/bin/python ora_analysis/ros_scavenger/visualize_ros.py --font-dir /tmp/ros_arial
.venv/bin/python ora_analysis/ros_scavenger/build_report.py
.venv/bin/python ora_analysis/ros_scavenger/validate_results.py
```

Arial이 이미 있으면 `setup_arial.py` 대신 해당 TTF 폴더를 `--font-dir`로 전달할 수 있다. Arial을 찾지 못하면 시각화가 중단되어 다른 폰트로 조용히 대체되지 않는다. Matplotlib·fontconfig 캐시는 `/tmp`를 사용한다.

## 통계 해석

- 1·2번의 기본 결과는 기존 FDR05 표다. 결측은 0이 아니다. 중복 심볼은 수치가 가장 좋은 행을 선택하지 않고 모두 제외한다.
- A GEO 결과는 두 군 4/4의 추정 count를 반올림한 별도 탐색적 PyDESeq2 fit이다. B는 5/5의 **제출 DESeq2 통계**를 쓰며 FPKM으로 DESeq2를 재계산하지 않는다.
- B 원 GEO에 기존 표에서 누락된 21개 유의 유전자가 있고, ROS catalog에는 Ptgs1·Ptgs2·Tacr1이 추가된다. 원 GEO 보완 결과를 별도 CSV와 Excel 시트에 저장한다.
- GO/KEGG 기능군 10개: focused ORA의 60검정, GSEA의 20검정, 표본 경로 score의 20검정을 각각 BH 보정한다. 기존 ORA는 원래 full-catalog FDR를 유지한다.
- GSEA는 signed-z 순위, weighted ES, 10,000 gene-label nulls, seed 20261102, 동일 부호 null 평균으로 NES 정규화를 사용한다. 메모리를 제한하는 streaming 구현이며 작은 universe의 모든 hit 패턴에서 dense running sum과 일치함을 검증한다.
- GSEA 6개 항목은 탐색적 신호다. 표본 score는 유전자 상관을 유지하는 보완적인 self-contained 검정이며 다중검정 후 유의하지 않다. 두 검정은 귀무가설이 다르다.
- power는 관측 효과로 사후 계산하지 않는다. 설계 민감도 표는 정규·등분산 t-test 가정의 참고값이다.

`validate_results.py`는 sparse ES와 dense KS, ORA와 Fisher exact, 모든 표본 순열 P와 SciPy exact enumeration, 입력 해시 및 출력 파일을 독립적으로 대조한다. 결과는 `results/independent_checks.json`에 기록한다.
