# 데이터셋별 KEGG·GSEA와 ROS 경로 비교

[한국어 보고서](Dataset_pathway_report_ko.md) · [Excel](results/Dataset_pathway_GSEA_KEGG.xlsx) · [ROS 비교 CSV](results/ROS_pathway_comparison.csv)

기존 ROS 한정 분석은 [별도 폴더](../ros_scavenger/README.md)에 보존했다. 이번 분석은 전체 KEGG catalog에서 각 데이터셋의 항산화·scavenger 기능군이 차지하는 위치를 평가한다. A/B의 기존 GEO 전체 통계·동일한 순위 파일을 재사용한다.

| Dataset | 전체 KEGG | 사전 정의 ROS GO | 전체 GSEA 검정 수 |
|---|---:|---:|---:|
| GSE245658 | 327 | 9 | 336 |
| GSE325653 | 338 | 9 | 347 |

KEGG는 matched size 15–500, ROS GO는 5–500으로 제한한다. 1차 FDR는 각 데이터셋에서 전체 KEGG+ROS GO nominal P를 BH 보정한 값이다. KEGG-only와 두 데이터셋 합산 보정도 제공한다. GO 9개는 전체 GO catalog가 아니다. RNA enrichment를 실제 ROS 제거 능력이나 인과적 기여율로 해석하지 않는다.

## Arial 그림

각 파일은 같은 이름의 PNG·PDF·SVG가 있다. PDF에 Arial TrueType을 임베딩하고 SVG는 Arial 글리프를 path로 변환한다.

| Dataset | KEGG GSEA | ROS 유의성 | GSEA 곡선 | KEGG ORA | 시료별 점수 |
|---|---|---|---|---|---|
| GSE245658 | [01](figures/GSE245658_01_KEGG_GSEA_overview.png) | [02](figures/GSE245658_02_ROS_pathway_significance.png) | [03](figures/GSE245658_03_ROS_GSEA_curves.png) | [04](figures/GSE245658_04_KEGG_ORA_overview.png) | [05](figures/GSE245658_05_ROS_sample_scores.png) |
| GSE325653 | [01](figures/GSE325653_01_KEGG_GSEA_overview.png) | [02](figures/GSE325653_02_ROS_pathway_significance.png) | [03](figures/GSE325653_03_ROS_GSEA_curves.png) | [04](figures/GSE325653_04_KEGG_ORA_overview.png) | [05](figures/GSE325653_05_ROS_sample_scores.png) |

## 재실행

저장소 루트에서 기존 ROS 분석의 Python 환경과 결과를 사용한다. 새 KEGG 원자료는 `ora_analysis/input/pathway_overview/`에 있다. 첫 다운로드 후에는 snapshot을 그대로 재사용한다.

```bash
.venv/bin/python ora_analysis/pathway_overview/fetch_kegg.py
.venv/bin/python ora_analysis/ros_scavenger/setup_arial.py --directory /tmp/ros_arial
.venv/bin/python ora_analysis/pathway_overview/analyze_pathways.py prepare
.venv/bin/python ora_analysis/pathway_overview/analyze_pathways.py analyze --dataset GSE245658
.venv/bin/python ora_analysis/pathway_overview/analyze_pathways.py analyze --dataset GSE325653
.venv/bin/python ora_analysis/pathway_overview/analyze_pathways.py summarize
.venv/bin/python ora_analysis/pathway_overview/visualize_and_report.py
.venv/bin/python ora_analysis/pathway_overview/validate_pathways.py
```

패키지 버전은 [기존 ROS requirements](../ros_scavenger/requirements.txt)와 같다. 원 GEO 재계산까지 처음부터 수행하려면 [기존 재실행 안내](../ros_scavenger/README.md)를 먼저 따른다. Arial이 없으면 시각화는 자동으로 다른 폰트로 대체되지 않고 중단한다.

GSEA null은 matched size별 10,000회, 사전 지정 ROS/context 경로의 크기에는 50,000회다. 동일 크기의 gene-label null은 재사용하며, cache는 Git에서 제외한다. cache를 삭제해도 seed로 재현할 수 있다. sample score는 모든 균형 시료 라벨을 열거해 정확 양측 P를 계산하며 gene-label GSEA와 귀무가설이 다르다.

`validate_pathways.py`는 모든 경로의 observed ES를 dense running sum으로, ORA P를 Fisher exact로, ROS 관련 sample score P를 SciPy exhaustive permutation으로 독립 대조한다. 기존 ROS 분석과 observed ES가 같은지도 확인한다. 입력 해시·보고서 링크·Arial PDF 임베딩까지 검사한다.
