# CTX-HPF 후보 발현 비교: 실행 상태와 재개 방법

현재 상태는 **외부 자료 확보 대기**이다. Allen 공식 다운로드 페이지와 NCBI 메타데이터에 접속할 때 HTTPS proxy CONNECT 단계에서 `403 Forbidden`이 발생했다. 라이브러리 설치나 TLS 인증서 검증 오류가 아니다. 현재 세포 수, genotype, donor/pool 수, MOp–VISp 발현 차이를 실제로 계산한 결과는 없다. 계획서에 적힌 약 8주령·B6J 배경 등의 조건은 실제 메타데이터 확인 전에는 분석 표본의 확정 조건으로 사용하지 않는다.

완료한 준비는 후보 22개와 공통 52개를 읽는 발현 분석 코드, 원본 메타데이터 필드 확인, 정확한 `MOp`·`VISp`만 선택하는 규칙, 대응 subclass 및 공유 cluster 비교, dot plot, 합성 입력으로 수행하는 코드 검증이다. 합성 입력은 임시 폴더에서 삭제되며 실제 프로젝트 결과에 저장되지 않는다. 후보 합집합은 62개이다.

## 먼저 확보·확인할 자료

공식 출처는 `input/rnaseq/source_inventory.csv`에 있다. `src/rnaseq/discover_downloads.py`는 공식 Allen 페이지에서 실제 다운로드 링크를 수집한다. 아직 접근할 수 없는 큰 expression artifact의 파일명·크기·정규화 방법은 추정하지 않았다. 공식 링크가 확인되면 동일 dataset/chemistry의 전체 셀 메타데이터와 발현 자료를 확보하고, checksum 및 다운로드 출처를 기록해야 한다.

메타데이터 확인 항목은 cell ID, 정확한 ROI, 원본 class/subclass/cluster, donor, biological pool, genotype, sex, age, chemistry, FACS 조건, batch/library이다. `schema_map.example.json`은 입력 열 이름을 지정하는 형식 예시이며 실제 Allen CSV에서 검증한 스키마가 아니다. 없는 필드를 새로 만들어 채우지 말고 예시에서 삭제한다. 필요한 열만 실제 확인된 이름으로 지정한다.

발현 분석은 H5AD 형식을 지원한다. 공식 원본이 다른 형식이면 공식 포맷에 맞춘 변환이 필요하다. cell ID와 gene ID를 원본 annotation에 정확히 맞추고 원본 정규화 단위를 유지해야 한다. 입력 파일이 CTX-HPF인지, 10x인지 SMART-seq인지 확인하며 두 기술을 같은 평균값으로 합치지 않는다. WMB의 `VIS`·`VIS-PTLp` 또는 Tasic 2018의 `ALM`을 `VISp`·`MOp`로 바꿔 쓰지 않는다.

## 실행

프로젝트 폴더에서 가상환경 Python을 사용한다. 아래 명령은 Linux 예시이며 Windows에서는 같은 스크립트에 `.venv\Scripts\python.exe`를 사용할 수 있다. `--project-root`를 생략하면 스크립트 위치에서 프로젝트 루트를 찾으므로 프로젝트 폴더를 옮길 수 있다.

```bash
python src/rnaseq/discover_downloads.py
python src/rnaseq/analyze_expression.py --metadata input/rnaseq/cell_metadata.csv --audit-only
```

실제 메타데이터에 따라 필요하면 `--schema-map input/rnaseq/schema_map.json`을 추가한다. 두 정확한 영역과 class/subclass annotation을 확보한 뒤 진행한다. donor/pool 정보가 없다고 세포를 독립 동물로 세지 않는다.

전체 transcriptome raw-count matrix임이 원본 문서에서 확인될 때:

```bash
python src/rnaseq/analyze_expression.py \
  --metadata input/rnaseq/cell_metadata.csv \
  --h5ad input/rnaseq/ctx_hpf.h5ad \
  --value-kind raw_counts \
  --full-transcriptome-counts-confirmed
```

HVG 또는 후보 유전자만 남은 counts에는 해당 확인 옵션을 붙이면 안 된다. 대신 메타데이터에 원본 전체 transcriptome 기준 UMI 합계가 보존되어 있다면 `--library-total-column 실제_열이름`을 사용한다. 스크립트는 외부 합계가 양의 정수이며 제공 matrix 합보다 작지 않은지 확인한다. raw counts는 전체 transcriptome 기준 CP10k를 만든 뒤 `log1p`한다. 단위가 이미 확인된 정규화 입력에는 `--value-kind log1p_cp10k`, `log2_cpm_plus1`, `linear_cpm` 중 실제 값에 맞는 것을 사용한다. matrix 값만 보고 단위를 추정하지 않는다.

유전자 이름이 Ensembl ID이면 `--gene-map`으로 `gene_id,gene_symbol` CSV를 지정한다. 중복 gene ID 또는 한 후보의 복수 feature mapping은 임의로 첫 항목을 고르지 않고 중단한다.

```bash
python src/rnaseq/check_pipeline.py
```

검증은 정확한 ROI 선택, 전체 유전자 기준 library 분모, 외부 library totals, 검출률, 미측정 후보 처리, 공유 cluster 가중치, unknown subtype 제외, raw counts 검증, 입력 부재 상태를 확인한다. 생물학적 결론을 검증하는 시험은 아니다.

## 실제 입력 확보 후 생성할 산출물

- `metadata_fields.csv`, `metadata_audit.json`: 실제 열과 missing 조건을 명시한다.
- `selected_cell_counts.csv`, donor/pool/batch × region count 표: RNA-seq 표집량이다. 조직 세포 비율로 해석하지 않는다.
- `candidate_gene_coverage.csv`: 없는 후보는 미확보이며 0 발현으로 바꾸지 않는다.
- `region_subclass_gene_summary.csv`: chemistry, genotype, sex, age, FACS별로 분리한 평균 발현과 검출률이다.
- `region_all_selected_cells_gene_summary.csv`: FACS/RNA-seq에 표집된 세포 평균이다. 조직 평균을 추정하지 않는다.
- `within_subclass_region_difference.csv`: MOp minus VISp 평균 발현과 검출률 차이이다. 값의 단위가 log1p(CP10k)일 때 이 차이는 유전자 log2FC가 아니다.
- `shared_cluster_region_difference.csv`, `shared_cluster_standardized_difference.csv`: 알려진 공유 cluster만 비교하고 같은 cluster 가중치를 주는 탐색 결과이다. unknown cluster 제외 수는 `shared_cluster_audit.json`에 남긴다.
- `figures/rnaseq/panel22_stratum_*.png`: 점 크기는 검출률, 색은 기록된 발현 단위이다. 서로 다른 chemistry/genotype 등의 층화 조건은 별도 그림에 둔다.
- `selected_candidate_expression.npz`: 62개 후보 중 확보한 유전자만 저장한 부분 matrix이다.

`status.json`이 `descriptive_expression_complete`일 때만 해당 실행의 결과가 완성된 것으로 판단한다. 입력 없이 실행하면 `pending_external_data`와 종료 코드 2를 반환한다. 실패한 실행에서는 이전 결과 파일이 남아 있을 수 있으므로 상태와 입력 출처를 함께 확인해야 한다.

이 단계의 결과는 기술적 비교이다. 알려진 annotation이 같은 경우에도 genotype/chemistry/age/sex/FACS가 unknown이면 대응 조건을 확인한 비교로 주장하지 않는다. 동물 또는 독립 biological pool의 관계와 영역–batch confounding을 실제 원본에서 확인할 때까지 pseudobulk 동물 검정 및 P/FDR는 생성하지 않는다. 지역 간 기저 mRNA 차이를 정상 동물의 PBM 반응 차이, adenosine 농도, 효소 flux, 단백질 기능 또는 수면 threshold로 환산하지 않는다.
