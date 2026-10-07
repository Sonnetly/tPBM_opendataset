# DEG FDR·log2FC를 반영한 cortex qPCR 후보 재평가

## 핵심 변경

**제공된 DEG 전체에서 같은 방향 overlap은 52개(up 34/down 18)다. 이전 ORA hits 기반 51개에 Tmem40 공통 down이 추가된다. 반대 방향 overlap은 23개로 같다.**

**B6 cortex 방향 검증 1차 8개: Grm2, Ryr2, Trpc5, Atp2b3, Cacna1c, Grin2a, Mrps6, Penk.**

사용자는 업로드 자료의 공통 방향을 후보 선정에 활용하고, 실제 B6 cortex tPBM 시료에서 qPCR로 확인할 계획이라고 설명했다. 위 6개 공통 up은 글루탐산·이온·칼슘 축을, Mrps6/Penk 공통 down은 미토콘드리아 번역/신경펩타이드 축을 대표한다. 질환 모델 맥락을 후보 제외 조건으로 사용하지 않고, 방향 일치와 기능을 중심으로 선정했다.

FDR와 log2FC는 qPCR 검증 순서의 참고값이다. Grin2a/Cacna1c의 기능적 중요성과 수치 근거를 함께 표시하고, Penk는 공통 down이며 변화 폭이 커서 첫 패널에 포함했다. 수동 순위는 가장 큰 fold change 순이나 FDR 순과 같지 않으며, 별도의 FDR 정렬 표도 제공한다.

**미토콘드리아를 포함해 10개로 확장할 경우 Bcs1l/Ppargc1a를 추가한다.** 두 유전자는 피질 PBM FDR05 표에서 up이며 해마 해당 비교의 표에는 없어서 공통 방향 8개와 구분했다. Kcnj6/Slitrk3/Shisa7/Dlgap1은 추가 기전 후보로 유지한다.

Per2/Cry2는 일주기 가설을 검증할 때 별도로 고려한다. 특히 Per2의 cortex 감소는 강한 수치 근거가 있으나, ZT/채취 시점/자극 후 시간과 수면 상태를 통제하지 않으면 해석하기 어렵다.

## 비교 조건 확인

- A: GSE245658, hippocampus, 808 nm, **PBM effect: 5xFAD+light vs 5xFAD**.
- B: GSE325653, prefrontal cortex, 1064 nm, **PBM effect: hypobaric hypoxia + PBM vs hypobaric hypoxia**.
- 전체 파일의 질병 효과, PBM-treated 5xFAD vs WT, 수막 및 림프관 내피세포 결과를 위 PBM 비교에 합치지 않았다.
- 두 비교는 조직·파장뿐 아니라 5xFAD/저산소 모델 배경과 분석 방법도 다르다. 이 차이를 cortex 영역별 반응성의 직접 근거로 해석하지 않는다. 조건은 제공 파일의 metadata에서 확인한 것이며 원 GEO를 외부에서 재검증한 것은 아니다.
- 목표 실험은 B6 cortex tPBM이다. B6의 정확한 substrain, 성별, 연령, 처리/채취 시간 등은 미제공이며, 원 연구의 질병 맥락을 B6 실험 조건으로 그대로 대입하지 않는다.

## DEG와 ORA 일치 여부

- A 해당 PBM 비교는 657행이다. 그 중 pathway-analysis 사용 표시가 True인 646행의 방향(up354/down292)이 이전 ORA 입력 수와 일치한다.
- A에는 placeholder `--` 4행, Zfp809/Bpifb9a 중복 심볼이 있어, 정확한 심볼 overlap은 유효한 651개 고유 심볼로 계산했다. B는 998행/998개 고유 심볼이다. 중복의 수치 중 유리한 값을 선택하지 않았으며 audit 표에 개별 행을 보존했다.
- 기존 ORA hits 537개(A)/903개(B)는 모두 올바른 조직·PBM 비교의 pathway-eligible DEG에 포함되고 방향이 일치한다.
- 같은 방향 52개와 반대 방향 23개에는 ambiguous duplicate symbol이 없다. 두 표에 없는 심볼은 미발현이나 0 효과로 처리하지 않았다.

## 1차 후보 수치

| 후보 | 해마 log2FC | 해마 gene FDR | 피질 log2FC | 피질 gene FDR | 피질 처리/대조 fold | 선택 이유 |
|---|---:|---:|---:|---:|---:|---|
| Grm2 | +1.1062 | 0.0214 | +0.4419 | 1.95e-07 | 1.358 | 피질 log2FC와 FDR가 기존 기능 후보 중 강하고 공통 up 및 glutamatergic 근거가 있어 우선 유지 |
| Ryr2 | +0.6854 | 0.0239 | +0.3068 | 0.00237 | 1.237 | 공통 up이며 피질 FDR가 0.01 미만인 칼슘 신호 후보 |
| Trpc5 | +0.8097 | 0.023 | +0.3639 | 0.00272 | 1.287 | 공통 up이며 피질 효과 크기와 FDR가 Kcnj6보다 뒷받침되어 이온채널 검증 우선순위를 올림 |
| Atp2b3 | +0.6907 | 0.0222 | +0.2607 | 0.003 | 1.198 | 공통 up이며 피질 FDR가 0.01 미만으로 칼슘 수송·항상성 축을 대표 |
| Cacna1c | +1.0542 | 4.28e-07 | +0.3065 | 0.0353 | 1.237 | 피질 gene FDR 근거는 다른 1차 후보보다 약하지만 전압 의존성 칼슘채널이라는 핵심 기능을 대표하도록 유지 |
| Grin2a | +0.8364 | 0.00078 | +0.3480 | 0.0498 | 1.273 | 공통 up이며 NMDA/글루탐산·칼슘 축을 대표하므로 첫 기전 패널에 포함; FDR와 효과 크기는 검증 순서의 참고값으로 제공 |
| Mrps6 | -0.8467 | 0.0309 | -0.3386 | 0.000431 | 0.791 | 공통 down이며 피질 gene FDR가 낮아 미토콘드리아 번역 가설의 우선 후보로 유지; 미토 기능 감소를 뜻하지 않음 |
| Penk | -1.8638 | 0.00671 | -1.2446 | 7.26e-05 | 0.422 | 두 PBM 비교에서 방향이 일치하고 피질 감소 폭이 커서 qPCR 방향 검증의 우선 후보로 포함 |

## 확장/일주기/새 overlap 후보

| 후보 | 분류 | 해마 log2FC | 피질 log2FC | 피질 gene FDR | 판단 |
|---|---|---:|---:|---:|---|
| Bcs1l | B_cortex_mito_extra | 표에 없음 | +0.2525 | 0.00915 | 피질 DEG 표의 up과 gene FDR를 근거로 미토콘드리아 호흡사슬 조립 가설의 추가 후보로 제안; 해마 PBM FDR05 표에는 없음 |
| Ppargc1a | B_cortex_mito_extra | 표에 없음 | +0.2328 | 0.0388 | 피질 DEG 표의 up과 gene FDR를 근거로 미토콘드리아 조절 가설의 추가 후보로 제안; 해마 PBM FDR05 표에는 없음 |
| Kcnj6 | B_mechanism_extension | +0.7272 | +0.2514 | 0.0487 | 피질 FDR가 약 0.049라 효과 크기·유의성 근거를 반영해 첫 패널보다 확장 검증에 배치 |
| Slitrk3 | B_mechanism_extension | +0.5693 | +0.2359 | 0.0142 | 공통 up의 시냅스/GABA 기능을 추가로 검증할 후보 |
| Shisa7 | B_mechanism_extension | +0.7042 | +0.2232 | 0.0359 | 공통 up의 시냅스/GABA 조절 가설을 보완하되 피질 효과가 작고 FDR 약 0.036임을 함께 고려 |
| Dlgap1 | B_mechanism_extension | +0.3575 | +0.1583 | 0.00321 | 공통 up이고 피질 FDR는 낮지만 약 12% 증가인 작은 효과여서 시냅스 구조 검증 확장에 배치 |
| Per2 | C_clock_time_controlled | 표에 없음 | -0.5525 | 6.58e-18 | 피질에서 유의한 감소가 강하지만 일주기/ZT·자극 후 시간 통제하에 별도 가설로 검증 |
| Cry2 | C_clock_time_controlled | +0.3458 | -0.2560 | 4.06e-06 | 피질 감소는 유의하나 해마 PBM에서는 증가함; 모델·파장 차이로 영역 특이성을 단정하지 않음 |
| Tmem40 | D_new_unannotated_overlap | -2.2219 | -1.0335 | 0.00112 | 원본 DEG 표에서 새로 확인되는 공통 down; 효과는 크지만 제공된 KEGG/GO 기능 근거가 없어 탐색 후보로 분류 |

## GO/미토콘드리아 해석

기존 [공통 up 유의 GO 38개](common_up_GO_38.csv)(BP29/MF9)는 이번 DEG 표로 변경되지 않는다. 해당 항목들의 실제 공통 hit 19개 모두 두 PBM 비교의 gene FDR<0.05 및 up 방향을 확인했다. [GO–gene 통계 연결 CSV](common_up_GO_with_gene_stats.csv)에는 원래 GO 항목 FDR와 각 유전자 FDR/log2FC를 별개 열로 연결했다.

미토콘드리아 관련 Ppargc1a/Bcs1l/Mrps6은 이제 gene-level FDR<0.05라는 근거가 확인된다. 그렇더라도 미토콘드리아 GO/KEGG enrichment가 유의하다는 뜻은 아니다. 기존 cortex Oxidative phosphorylation은 hits 0, FDR1이고 대표 미토 GO 항목의 FDR도 유의하지 않다. 단일 유전자 발현으로 ATP 생산 증가/감소를 판정하지 않는다.

Tmem40은 해마 log2FC -2.2219/FDR0.0252, 피질 -1.0335/FDR0.00112이며 공통 down이다. 이 파일들의 기존 KEGG/GO hits에는 없고 제공된 기능 주석이 없으므로 기전 중심 1차 패널에 자동 추가하지 않고 탐색 대상으로 표시한다.

## 실험 및 수치 해석

- 피질 1차 후보 중 많은 증가 효과는 처리/대조 약 1.2–1.4배다. 통계적으로 유의한 것과 qPCR에서 큰 변화가 예상되는 것은 다르다. assay 효율·정규화·생물학적 변동을 확인하며, Cortex lfcSE가 없어 이 표만으로 필요한 표본 수나 검출력을 계산하지 않았다.
- 독립적인 cortex 영역별 control/tPBM 비교와 region×treatment 차이를 평가한다. 동일 동물의 여러 영역을 측정하면 짝지음/동물 효과를 반영한다. 채취 시점, ZT, 자극 후 시간과 수면/각성 상태를 통제한다.
- qPCR은 발현을 검증하며 전기생리 반응, 채널 전류, 단백질 활성이나 미토콘드리아 호흡을 직접 측정하지 않는다. 해당 기능 readout과 연결해서 가설을 평가한다.
- gene FDR는 제공 값이다. FDR<0.05로 이미 선택된 파일만으로 BH를 재계산하거나 전체 발현 배경으로 새로운 ORA를 수행하지 않았다.
- log2FC/통계는 해당 비교 안에서 해석한다. baseMean은 서로 다른 처리/정규화 맥락에서 온 제공 값으로, 두 연구 간 절대 발현량 척도로 비교하지 않는다.
- Cortex 통계는 submitter-provided DESeq2 statistics라고 기록되어 있고 FPKM 제출 파일에서 새로 계산한 값이 아니다. A는 count columns의 PyDESeq2 재계산으로 기록되어 있다. 원시 count 및 디자인은 이번 파일에 없어 그 분석을 다시 수행하지 않았다.

## 재실행 및 결과 파일

`python3 ora_analysis/analyze_deg_and_reprioritize.py`로 재현한다. 수동 선별 이유는 [우선순위 설정 TSV](../deg_aware_priority_config.tsv)에 저장했다. [별도 qPCR 요약 Markdown](../../qPCR_candidate_summary_ko.md)에는 패널과 주요 CSV 링크를 모았다.

- [전체 Excel](DEG_aware_qPCR_and_overlap.xlsx): 전체 후보, 정확한 비교별 DEG, 완전한 심볼 overlap, 원본·중복 audit, GO 연결.
- [공통 방향 52개 CSV](same_direction_52.csv), [반대 방향 23개 CSV](opposite_direction_23.csv): log2FC/gene FDR와 방향을 포함한 overlap.
- [qPCR 후보 CSV](qPCR_candidates_with_DEG.csv): 이전 모든 후보와 새로운 Tmem40의 수치 및 현재 우선순위.
- [피질 FDR 정렬 CSV](candidate_cortex_FDR_sort.csv): 수동 순위와 구분되는 단순 cortex FDR 정렬 참고표.
- [공통 up 유의 GO 38개 CSV](common_up_GO_38.csv): 항목별 ORA 통계·공통 hit·공식 항목 링크.
- [GO–gene 통계 연결 CSV](common_up_GO_with_gene_stats.csv): 공통 유의 GO 38개의 실제 공통 hit에 대한 gene-level 수치.
- [해마 PBM DEG CSV](selected_hippocampus_PBM.csv), [전전두피질 PBM DEG CSV](selected_cortex_PBM.csv): 이번 분석에 사용한 정확한 비교별 입력.
- [ORA–DEG 대조 CSV](ORA_DEG_alignment.csv), [전체 비교 요약 CSV](DEG_group_summary.csv): 원본 비교와 이전 ORA 입력·방향 확인.
- [새 overlap CSV](new_unannotated_overlap.csv), [전체 심볼 합집합 CSV](all_selected_symbol_union.csv): ORA hit에서 빠졌던 공통 유전자와 전체 비교 범위.
- [중복 심볼 CSV](duplicate_symbols_audit.csv), [미지정 심볼 CSV](placeholder_symbols_audit.csv): 원본 행을 보존한 audit.
- [원본 DEG CSV](../input/DEG_FDR05_GSE245658_GSE325653.csv), [원본 ORA CSV](../input/ORA_fullcatalog_KEGG_GO.csv).
