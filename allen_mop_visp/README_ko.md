# Allen MOp–VISp 분석

목적은 정상 상태 atlas에서 MOp와 VISp의 세포 구성과 대응 세포 유형 안의 후보 유전자 발현을 비교하고, 기존 PBM 결과와 연결하는 것이다. 후속 회로 분석은 NREM 관련 영역으로의 투사를 우선 확인한다. 정상 상태 발현과 atlas의 축삭 신호를 PBM 처리 효과 또는 수면 유발 효과로 해석하지 않는다.

처음 읽을 문서: [M1–V1 데이터 해석 안내서](docs/M1_V1_data_analysis_guide_ko.md)(0절에 2026-10-16 갱신 요약). 뉴런 유형별 발현과 다음 실험 후보는 [M1·V1 뉴런 유형별 발현 보고서](results/neuron_type_report/M1_V1_neuron_type_expression_report_ko.md)에 있다. 기존 결과 CSV의 구성·발현·개별 투사 요약 그림과, 완료/미완료 분석·가설·추가 실험·한계를 함께 정리했다.

[RNA-seq 전체 6개 dotplot의 MOp–VISp 인접 행 비교](figures/rnaseq/paired/README_ko.md): 공통 색 범위의 요약 그림과 평균·검출률 숫자판을 제공한다. [전체 요약 PDF](figures/rnaseq/paired/all_six_paired_dotplots.pdf), [전체 숫자판 PDF](figures/rnaseq/paired/all_six_numeric_detail_pages.pdf). P57/P59의 미표집 영역은 N/A로 유지했다.

[Dotplot에서 우선 확인할 gene × 세포 유형](results/rnaseq/paired_dotplots/gene_observation_review_ko.md): 전체 62개 후보를 검토하고 P58 평균·검출률·공유 cluster·cohort 민감도를 연결했다. 통계적으로 유의한 지역 DEG나 최종 qPCR 순위가 아닌 관측 후보이다.

2026-10-12 후속 검토: [MOp L4·구성 그림의 분모와 표현 검토](docs/composition_visualization_review_ko.md), [실제 810 nm 연속광·180분 조건의 M1 PBM 수면 반응 가설](docs/M1_PBM_sleep_sensitivity_hypothesis_ko.md)을 추가했다. 기존 구성 그림은 유지했으며, 가설은 실제 지역별 PBM 효과를 계산한 결과가 아니다. 사용자 자극 조건과 새로 확인한 원 논문·공개 보충자료를 함께 기록했다.

## 저장 위치

현재 정본은 저장소의 `allen_mop_visp/` 폴더다. 클라우드 절대 경로는 `/workspace/tPBM_opendataset/allen_mop_visp`이다. 사용자의 커밋·푸시 요청에 따라 이전 `/workspace/Allen Dataset_M1_V1`에서 파일 68개와 inventory를 복사하고 SHA256으로 일치를 확인했다. 지정한 Windows 경로는 `G:\내 드라이브\0. Research\2026\openRNAseq_PBM\Allen Dataset_M1_V1`이며 클라우드에 연결되어 있지 않다. 로컬 저장이 완료됐다는 의미가 아니다.

`docs/`에 계획서와 자극 정보, `input/`에 출처와 입력, `src/`에 코드, `results/`에 계산 결과, `figures/`에 그림을 보관한다. 임시 네트워크 진단은 `logs/`에 두고 Git에서 제외한다. 검증 기록은 `results/validation/`과 각 분석의 결과 폴더에 보존한다. 세부 방침은 [파일 관리 안내](docs/file_management_ko.md)를 참고한다.

## 현재 완료와 미완료

- 기존 PBM 입력을 다시 확인하고 22개 패널·공통 DEG 52개·GO 38개·공통 GO hit 19개의 연결 표를 만들었다. 22개 패널과 공통 52개를 합친 탐색 목록은 62개이다. 원본 복사본·SHA256·검증 결과는 `input/pbm/`과 `results/pbm/`에 있다.
- 공식 유전자 목록으로 실측 패널을 확인했다. B6J 500-gene 패널은 핵심 후보 중 Penk만 측정하고, 탐색 62개 중 7개를 측정한다. B6N Zhuang-1 1122-gene 패널은 핵심 후보 중 Adora1·Adora2a·Penk·Trpc5를 측정한다. 미측정 유전자는 0 발현이 아니다. `results/coverage/`에 전체 표가 있다.
- 원본 다운로드 재시도가 성공했다. 공식 ABC manifest의 B6J MERFISH metadata 5개(1.57 GB)를 받아 크기와 공개 MD5를 검증했다. MOp 83,976개·VISp 61,884개 세포의 전체 영역·층·subclass 구성 비율을 실제 계산했다. [MERFISH 보고서](results/merfish/composition_report_ko.md)에 결과와 분모가 있다.
- RNA-seq은 공식 GEO GSE185862의 10x 원본(발현 행렬 5.72 GB 및 metadata)을 확보해 exact MOp 123,727개·VISp 29,763개 세포의 62개 후보 발현을 실제 계산했다. 원본 전체 31,053개 유전자의 세포별 count 합계를 정규화 분모로 보존했다. [실제 발현 보고서](results/rnaseq/RNAseq_results_ko.md), [RNA-seq 상태 파일](results/rnaseq/status.json), [재현 방법](docs/RNAseq_analysis_ko.md)에 결과와 조건이 있다.
- 공식 SMART-seq v4 원본도 받아 MOp 6,588개·VISp 15,572개와 62개 후보를 별도로 분석했다. 조건을 맞춘 자체 비교 38개 strata를 표시하며 10x와 직접 대응되는 strata는 없다. 다른 reporter·연령의 보조 결과를 10x 방향 재현으로 해석하지 않는다. `results/rnaseq/smartseq/`에 별도 결과가 있다.
- Allen Connectivity의 MOp·VISp primary injection 실험 168개를 실제 확보했다. WT MOp 8개·VISp 33개를 주 분석으로, Cre line MOp 23개·VISp 104개는 line별 별도 탐색으로 보존했다. NREM 관련 atlas label 7개의 좌·우·양측 3,528개 unionize 기록을 확인하고 연령·주입 순도 민감도를 계산했다. PZ는 ontology exact label이 없어 미평가로 남겼다. ZI에서 관찰된 차이는 후속 해부학 후보이며 Lhx6 세포 targeting·시냅스·NREM 효과는 미확인이다. `results/connectivity/`에 개별 실험·요약·상태·보고서·별도 Excel을 보관한다.
- 2026-10-15 Zhuang-ABCA MERFISH(Zhang et al., Nature 2023, ABC Atlas release `20250531`) 중 Zhuang-ABCA-1, -2(coronal)와 -3, -4(sagittal) 네 동물을 분석했다. B6N(C57BL/6NCrl) 동물이고 1,122-gene panel이며 Allen WMB taxonomy와 CCF 2020 parcellation이 붙어 있다. 공개 MD5를 검증한 원본 약 6.1 GB로 MOp 98,334개와 VISp 80,373개 세포의 subclass 내 MOp−VISp log2 평균 차이, panel 전체 cell-class discordance, 동물 간 부호 일치를 계산했다. 후보 62개 중 12개가 측정된다. 이어서 subclass 안의 supertype 구성 차이를 분리하는 표준화 점검(`src/merfish/zhuang_supertype_check.py`)을 추가했다. Adora1은 L2/3 IT와 L4/5 IT의 MOp 저발현, Pvalb의 MOp 고발현(+0.07로 작음)이 네 동물에서 같고 supertype 표준화 후에도 유지됐다. 네 동물 짝지은 지역 DEG는 `src/merfish/zhuang_paired_deg.py`로 계산했다(subclass×gene 496개 통과, 대부분 Glut의 영역 정체성 유전자). L6 CT의 음의 차이는 대부분 supertype 구성 효과였다. 차이 크기는 log2 단위 0.1~0.3으로 작고 sagittal 동물에서 2~3배 크다. ABCA-1 발현 행렬은 dense 약 18.7 GB여서 기존 로더가 메모리 초과로 종료됐고, `src/merfish/run_zhuang.py`의 `load_expression`을 행 블록 읽기로 고쳤다(값 동일). B6J 주 분석과 합산하지 않는다. [Zhuang 보고서](results/zhuang/zhuang_report_ko.md)에 데이터셋 설명, 수행 내용, 결과와 한계가 있다.
- 2026-10-15~16 [M1·V1 뉴런 유형별 발현 보고서](results/neuron_type_report/M1_V1_neuron_type_expression_report_ko.md)를 추가했다. GEO 10x 원본 전체에서 영역×donor×subclass pseudobulk(`src/rnaseq/pseudobulk_10x.py`)를 만들고, 영역 안 뉴런 유형 marker를 MOp donor에서 찾아 별도 batch의 VISp donor에서 재현 여부를 확인했다(`src/rnaseq/subclass_markers_10x.py`). OXPHOS와 adenosine 세트는 결과 확인 전에 `src/rnaseq/focus_gene_sets.py`에 고정했다.
  - 에너지 대사(`src/rnaseq/energy_gene_detail.py`). OXPHOS 복합체 I~V와 Ppargc1a는 Pvalb와 L5 PT에서 가장 높고 MOp와 VISp에서 같은 패턴이다(세트 점수 r 0.92). 세트마다 유전자별 값, count 기여 비율, leave-one-out 결과를 보조 그림 S1~S5와 `results/rnaseq/whole_transcriptome/energy_detail/`에 저장했다. 생합성 세트는 사실상 Ppargc1a, creatine 세트는 Ckb가 좌우한다.
  - adenosine 경로(`src/rnaseq/adenosine_region_compare.py`). Adora1은 L6 CT, L6b, L5/6 NP에서 가장 높다. 10x 지역 차이는 donor와 batch에 혼입되므로 공통 유전자로 MERFISH 짝지은 차이와 대조했다. Glut는 부호 일치 86~93%, GABA는 69~84%였다. 이 범위 안에서 MOp IT 뉴런은 VISp보다 Adora1이 낮고 Adk가 높다. GABA의 10x Adora1 지역 차이(MOp +0.7 이상)는 MERFISH에서 재현되지 않아 해석하지 않는다. 결과는 `results/rnaseq/whole_transcriptome/adenosine_region/`에 있다.
  - Cox4i1과 Ndufa13은 이 GEO release에서 모든 세포의 count가 0이며 주석 문제로 보인다.
  - 그림은 `figures/report_m1v1/`의 본문 Fig1~Fig7(Fig2는 a, b, c)과 보조 Fig S1~S6이다.
- 공간 주 자료는 B6J 단일 동물 MERFISH이며 동물 간 통계 검정을 수행하지 않는다. RNA-seq의 세포 수를 조직 구성 비율이나 독립 동물 수로 사용하지 않는다.
- WT B6J와 B6J 배경 transgenic은 분리한다. 정확한 MOp와 VISp가 있는 자료만 직접 비교하며 ALM이나 VIS 전체를 대신 사용하지 않는다.
- 국소 자극 ROI는 V1 AP가 잠정값이고 bregma→CCF 등록을 수행하지 않아 보류했다. mm·bregma AP/ML 및 뇌 표면 DV, core 직경 400 µm는 사용자 확인을 받았다. M1은 (1.34, 2.46, 0), V1은 (약 −2.3, 2.46, 0)이다. 사용자 좌표는 [기록 파일](docs/stimulation_coordinates_user.json)에 있다. 400 µm만으로 광학적 ROI를 정하지 않는다.

## 설치

Python 3.12를 사용한다. Linux/macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
```

Windows에서는 이 폴더에서 PowerShell로 `./setup_windows.ps1`를 실행하거나 아래 명령을 사용한다. PowerShell 설치 경로는 이 클라우드에서 실행 검증하지 않았다.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip check
```

## 분석 재개

모든 명령은 이 프로젝트 폴더에서 실행한다. Windows에서는 `.venv/bin/python` 대신 `.venv\Scripts\python.exe`를 사용한다.

이미 보존한 PBM 원본과 공식 gene list에서 완료된 연결 표·실측 audit를 다시 생성한다:

```bash
.venv/bin/python src/build_pbm_panels.py
.venv/bin/python src/audit_gene_coverage.py
```

B6J MERFISH의 공식 manifest와 체크섬으로 metadata를 받아 전체 영역·해부학적 층·subclass의 구성 비율을 계산한다:

```bash
.venv/bin/python src/merfish/run_merfish.py --download
```

B6N Zhuang-ABCA MERFISH의 manifest를 확인하고 동물별 원본을 받아 subclass 내 MOp−VISp 비교를 계산한다. 원본은 `input/zhuang/cache/`에 저장되며 Git에서 제외된다. 메모리 7 GB 환경에서 실행을 확인했다.

```bash
.venv/bin/python src/merfish/run_zhuang.py --list --animals Zhuang-ABCA-1 Zhuang-ABCA-2 Zhuang-ABCA-3 Zhuang-ABCA-4
.venv/bin/python src/merfish/run_zhuang.py --download --animals Zhuang-ABCA-1 Zhuang-ABCA-2 Zhuang-ABCA-3 Zhuang-ABCA-4
.venv/bin/python src/merfish/zhuang_supertype_check.py --animals Zhuang-ABCA-1 Zhuang-ABCA-2 Zhuang-ABCA-3 Zhuang-ABCA-4
.venv/bin/python src/merfish/plot_zhuang_supertype.py
```

`--list`도 `results/zhuang/zhuang_report.json`을 manifest 확인 상태로 덮어쓰므로 `--download` 뒤에 따로 실행하지 않는다.

뉴런 유형별 발현 보고서의 계산과 그림을 다시 만든다. GEO 10x 원본(5.72 GB)을 `input/rnaseq/raw/`에 받은 뒤 실행하며, pseudobulk 추출은 2개 worker로 약 2분 걸렸다.

```bash
.venv/bin/python src/rnaseq/download_geo.py GSE185862_expression_matrix_10x.hdf5
.venv/bin/python src/rnaseq/pseudobulk_10x.py --workers 2
.venv/bin/python src/rnaseq/subclass_markers_10x.py
.venv/bin/python src/merfish/zhuang_paired_deg.py
.venv/bin/python src/report_figures_m1v1.py
.venv/bin/python src/rnaseq/energy_gene_detail.py
.venv/bin/python src/rnaseq/adenosine_region_compare.py
```

공식 GEO에서 10x 입력을 내려받아 후보 subset과 전체 유전자 정규화 분모를 다시 생성한다. 대형 원본을 받지 않아도 아래의 발현 계산은 이동용 묶음에 포함된 subset으로 재현할 수 있다:

```bash
.venv/bin/python src/rnaseq/download_geo.py GSE185862_expression_matrix_10x.hdf5
.venv/bin/python src/rnaseq/download_geo.py GSE185862_metadata_10x.csv.gz
.venv/bin/python src/rnaseq/download_geo.py GSE185862_CTX-HPF-73538sample_metadata.xlsx
.venv/bin/python src/rnaseq/extract_geo_subset.py --workers 2
```

검증된 raw counts 후보 subset에서 같은 세포 유형 안의 기술 통계를 생성한다. 후보 62개의 합계를 정규화 분모로 사용하지 않는다:

```bash
.venv/bin/python src/rnaseq/analyze_expression.py --project-root . --h5ad input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad --metadata input/rnaseq/selected_cell_metadata.csv.gz --value-kind raw_counts --library-total-column full_library_counts
.venv/bin/python src/rnaseq/cohort_screening.py
.venv/bin/python src/rnaseq/summarize_findings.py
.venv/bin/python src/build_atlas_workbook.py
```

다른 입력의 메타데이터 스키마와 gene ID는 `--schema-map`, `--gene-map`으로 명시한다. 해당 10x 자료는 수컷 `Snap25-IRES2-Cre;Ai14` reporter이며 WT MERFISH와 구분한다. 지역별 donor·library date·일부 연령이 달라 지역 효과를 독립적으로 추정하는 p-value/FDR를 계산하지 않는다. 공유되는 P58 연령의 비교와 P57–59 탐색 비교의 범위를 구분한다. RNA-seq 세포 수를 조직 구성 비율로 쓰지 않는다.

NREM 관련 projection의 실제 공식 API 응답을 확보하고, 보존된 응답으로 표적별 기술 통계를 재현한다:

```bash
.venv/bin/python src/connectivity/run_connectivity.py --download
.venv/bin/python src/connectivity/build_connectivity_workbook.py
.venv/bin/python src/connectivity/validate_connectivity.py
```

이미 원본 응답이 있으면 `--download` 없이 같은 runner로 재현할 수 있다. 실제 관측 결과는 [투사 보고서](results/connectivity/NREM_projection_report_ko.md)와 [별도 Excel](results/connectivity/MOp_VISp_NREM_projection.xlsx)에 있다. 기능 문헌과 ontology 대응은 [NREM 표적 정의](docs/NREM_projection_target_definition_ko.md)에 있다. 시상 망상핵의 Allen label은 `RT`이며 뇌간의 `TRN`과 구분한다. 분자 경로의 전체 유전자 추출·모듈·GSEA 및 회로 기능 검증 방법은 [후속 분석 제안](docs/downstream_pathway_analysis_options_ko.md)에 정리했다. 새 분자 경로 점수는 아직 계산하지 않았다.

## 네트워크와 원본 보존

현재 런타임에서 ABC S3, 공식 GEO HTTPS archive와 Allen Connectivity HTTPS API에 실제 접근·다운로드했다. 과거 Allen CTX-HPF S3 경로는 여전히 막혀 있어 동일 연구의 공식 GEO 원본을 사용한다. 네트워크 설정 초안 저장은 런타임 적용을 뜻하지 않는다. 주소·응답·검증 정보는 `input/merfish/`, `input/rnaseq/`, `input/connectivity/`, [네트워크 기록](docs/network_requirements.json)에 보존한다.

대형 다운로드는 `input/merfish/cache/`, `input/zhuang/cache/`와 `input/rnaseq/raw/`에 남기고 Git·이동용 ZIP에서 제외한다. 분석에 사용한 MOp/VISp metadata, 후보 raw-count subset, 전체 유전자 count 합계, 공식 출처·hash는 포함한다. 원본이 필요하면 같은 코드와 URL로 다시 받는다.

## 이동용 묶음

```bash
.venv/bin/python src/package_project.py
```

이 명령은 `exports/Allen_Dataset_M1_V1.zip`을 만들고, `file_inventory.json`에 각 파일의 SHA256을 기록한다. ZIP의 CRC 및 파일별 해시를 확인한다. OS별 가상환경·임시 로그·이전 ZIP·대형 원본 캐시는 제외한다. 재현용 subset·코드·결과·보고서·출처는 포함한다. 압축 안의 최상위 폴더 이름은 로컬 사용을 위해 `Allen Dataset_M1_V1`로 유지한다. Windows 지정 위치의 부모 폴더에 압축을 풀면 같은 프로젝트 구조를 사용할 수 있다.

이전 클라우드 작업 폴더에서 만든 최초 ZIP도 `exports/Allen_Dataset_M1_V1_initial.zip`으로 보존했다. 최신 재현 안내는 현재 README와 이름에 `initial`이 없는 ZIP을 사용한다.
