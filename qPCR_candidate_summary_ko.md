# B6 cortex tPBM qPCR 후보 요약

가설: **Regional differences in the sleep response to tPBM arise from differences in local cortical responsiveness to stimulation.**

두 데이터셋에서 **같은 방향으로 변하는 유전자**를 우선 선정하고, 제공된 FDR·log2FC를 검증 순서의 참고값으로 사용한다. 실제 방향은 B6 cortex tPBM 시료의 qPCR로 확인한다.

## 1차 qPCR 패널: 8개

아래 방향은 해마·전전두피질 PBM 비교에서 일치한다. 수치는 GSE325653 전전두피질 PBM 비교의 제공 값이며, 전체 후보와 두 비교의 수치는 [qPCR 후보 CSV](ora_analysis/deg_validation/qPCR_candidates_with_DEG.csv)에서 확인할 수 있다.

| Gene | 공통 방향 | 피질 log2FC | 피질 gene FDR | 검증할 축 |
|---|---|---:|---:|---|
| **Grm2** | up | +0.4419 | 1.95e−7 | 글루탐산 신호 조절 |
| **Ryr2** | up | +0.3068 | 0.00237 | 세포 내 칼슘 방출 |
| **Trpc5** | up | +0.3639 | 0.00272 | 양이온 채널 |
| **Atp2b3** | up | +0.2607 | 0.00300 | 칼슘 배출·항상성 |
| **Cacna1c** | up | +0.3065 | 0.0353 | 전압 의존성 칼슘채널 |
| **Grin2a** | up | +0.3480 | 0.0498 | NMDA 수용체 |
| **Mrps6** | down | −0.3386 | 0.000431 | 미토콘드리아 번역 |
| **Penk** | down | −1.2446 | 7.26e−5 | 신경펩타이드 신호 |

첫 6개는 글루탐산·이온·칼슘 축을, Mrps6와 Penk는 미토콘드리아 번역·신경펩타이드 축을 대표한다. 기능과 방향을 함께 고려한 패널이며, 단순 FDR 정렬은 [피질 FDR 정렬 CSV](ora_analysis/deg_validation/candidate_cortex_FDR_sort.csv)에 따로 저장했다.

## 미토콘드리아 확장: 10개 패널

위 8개에 **Bcs1l·Ppargc1a**를 추가한다. 두 유전자는 피질 FDR<0.05 표에서 up으로 확인되며, 해마 해당 PBM 비교의 FDR<0.05 표에는 없다.

| Gene | 피질 방향 | 피질 log2FC | 피질 gene FDR | 검증할 축 |
|---|---|---:|---:|---|
| **Bcs1l** | up | +0.2525 | 0.00915 | 미토콘드리아 호흡사슬 복합체 III 조립 |
| **Ppargc1a** | up | +0.2328 | 0.0388 | 미토콘드리아 조직화·조절 |

Kcnj6·Slitrk3·Shisa7·Dlgap1은 추가 기전 후보로 유지한다. 후보별 KEGG/GO ID는 [qPCR 후보 CSV](ora_analysis/deg_validation/qPCR_candidates_with_DEG.csv)의 `evidence_ids`, 해당 원본 ORA 행은 [후보·공통 유전자 근거 CSV](ora_analysis/qpcr_go/candidate_and_concordant_evidence.csv)에서 확인한다.

## 전체 overlap과 공통 up 유의 GO

제공된 FDR<0.05 DEG 표의 두 PBM 비교에서 **공통 방향은 52개(up 34/down 18)**다. 기존 ORA hits 기반 51개에 **Tmem40 공통 down**이 추가됐으며, 반대 방향 overlap은 23개다.

- [공통 방향 52개 CSV](ora_analysis/deg_validation/same_direction_52.csv): 두 비교의 방향·log2FC·gene FDR.
- [반대 방향 23개 CSV](ora_analysis/deg_validation/opposite_direction_23.csv): 방향이 다른 유전자와 제공 통계.
- [추가된 Tmem40 CSV](ora_analysis/deg_validation/new_unannotated_overlap.csv): 새로 확인된 공통 down 유전자.

두 데이터셋의 up ORA에서 모두 FDR<0.05인 **GO 항목은 38개(BP 29/MF 9)**다. 시냅스 조직화·신호·가소성, 이온채널·수송·막전위, 신경 형태·축삭 발달, 인지·학습·행동 관련 항목이 포함된다. 실제 공통 hit의 합집합은 19개이며, 모두 두 PBM 비교에서 up이고 gene FDR<0.05다.

- [공통 up 유의 GO 38개 CSV](ora_analysis/deg_validation/common_up_GO_38.csv): GO ID·이름·양쪽 ORA FDR·공통 hit·공식 항목 링크.
- [GO–gene 통계 연결 CSV](ora_analysis/deg_validation/common_up_GO_with_gene_stats.csv): 19개 유전자의 131개 GO–gene 연결과 두 비교의 log2FC·gene FDR. GO 항목 FDR와 gene FDR는 별도 열이다.
- [공통 GO hit 19개 CSV](ora_analysis/qpcr_go/common_up_GO_genes_19.csv): 유전자별 관련 GO 항목.
- [GO 주제별 요약 CSV](ora_analysis/qpcr_go/GO_theme_summary.csv): 38개 항목을 기능별로 묶은 설명용 요약. 이 파일의 `core_panel_genes_in_theme` 열은 DEG 수치를 받기 전 패널 기준이며, 현재 1차 패널은 위 8개다.

## 원본 및 상세 자료

- [원본 DEG CSV](ora_analysis/input/DEG_FDR05_GSE245658_GSE325653.csv)
- [원본 KEGG/GO ORA CSV](ora_analysis/input/ORA_fullcatalog_KEGG_GO.csv)
- [사용한 해마 PBM DEG CSV](ora_analysis/deg_validation/selected_hippocampus_PBM.csv)
- [사용한 전전두피질 PBM DEG CSV](ora_analysis/deg_validation/selected_cortex_PBM.csv)
- [ORA–DEG 입력·방향 대조 CSV](ora_analysis/deg_validation/ORA_DEG_alignment.csv)
- [상세 분석 보고서](ora_analysis/deg_validation/DEG_aware_qPCR_report_ko.md)
- [전체 Excel](ora_analysis/deg_validation/DEG_aware_qPCR_and_overlap.xlsx)
- [재실행 방법](ora_analysis/README.md)

원본 후보 수치를 대조했고, 분석을 재실행해 CSV·JSON·상세 Markdown의 동일성과 Excel 셀 값의 동일성을 확인했다.
