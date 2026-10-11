# MOp–VISp 분석 진행 결과

작성일: 2026-10-11. 현재까지 실제 계산한 결과는 기존 PBM 자료의 후보 연결·재검증과 Allen 공식 실측 패널 audit이다. **MOp–VISp 세포 비율과 후보 발현 차이는 원본 자료 다운로드가 차단되어 아직 계산하지 않았다.**

## 완료한 분석

기존 저장소 commit `e430e8d2d2827ccdd222e8828d77227e937b6f1b`의 원본 13개를 보존하고 Git 원본·SHA256과 대조했다. 분석을 원본 저장소 없이 다시 실행할 수 있다. 공통 DEG는 52개(up34/down18), 반대 방향은 23개다. 두 연구 up ORA에서 동시에 유의한 GO는 38개(BP29/MF9)이며 실제 공통 hit의 합집합은 19개다. 공통 hit가 없는 GO 1개도 그대로 보존했다.

22개 screening panel과 공통 52개의 합집합은 62개다. 22개를 핵심 표시 패널, 62개를 탐색 범위로 고정했다. gene FDR와 GO 항목 FDR를 분리하고, 제공 FDR05 목록에서 빠진 유전자는 `missing_from_FDR05_table`로 기록했다. 빈 수치를 0으로 바꾸거나 selected subset에서 FDR를 재계산하지 않았다.

| 후보 | 피질 PBM log2FC | gene FDR | 현재 해석 |
|---|---:|---:|---|
| Nt5e | −0.5367 | 0.00744556 | 기존 저산소 PFC PBM 자료의 감소; 정상 MOp–VISp 차이는 미확인 |
| Adora2a | −1.2929 | 0.00932146 | 동일 조건의 감소; 지역별 기저 발현과 기능은 별도 검증 |

기존 PBM 연구의 해마는 5xFAD, 피질은 저산소 모델이다. 정상 B6 atlas의 기저 발현과 이 PBM 처리 효과를 같은 값으로 합치지 않는다. Entpd1·Slc29a1·Slc29a2·Adk·Ada·Adora1이 FDR05 표에 없는 것을 변화 없음으로 판정하지 않는다.

## 공식 실측 패널 확인

Allen 공식 SDK 저장소 commit `b16ee223d3466a57b30066973413ccc948c2e983`의 HTML gene list를 파싱해 패널의 전체 수·중복·후보 교집합을 검사했다. 이는 실측 가능성 확인이며 해당 유전자의 실제 세포 발현량을 계산한 결과가 아니다.

| 자료 | 패널 수 | 핵심 22개 중 실측 | 탐색 62개 중 실측 |
|---|---:|---|---:|
| B6J MERFISH-C57BL6J-638850 | 500 | Penk | 7 |
| B6N Zhuang-ABCA-1 | 1122 | Adora1, Adora2a, Penk, Trpc5 | 12 |

B6J 탐색 가능 7개는 Ccbe1·Ecel1·Epha4·Npnt·Penk·Stxbp6·Zeb2다. 핵심 adenosine 후보 8개는 B6J 패널에 모두 없다. 세포 구성은 B6J MERFISH에서, 폭넓은 후보 발현은 CTX-HPF RNA-seq에서 보는 설계를 유지한다. 패널에 없는 유전자는 미측정이며 0 발현이 아니다. B6N 결과를 B6J의 동물 반복으로 합산하지 않는다.

## 아직 계산하지 않은 결과와 준비한 코드

| 단계 | 현재 상태 | 필요한 실제 입력 |
|---|---|---|
| MOp–VISp 구성·해부학적 층 분포 | 코드 준비, 실제 계산 보류 | B6J annotated cell metadata·CCF 영역·좌표 |
| 같은 세포 유형의 22/62개 발현 | 코드 준비, 실제 계산 보류 | exact MOp/VISp CTX-HPF cell metadata·발현 행렬·단위·gene ID |
| donor/pool pseudobulk 검정 | 미실행 | 독립 동물/pool, genotype·chemistry·선별·batch 설계의 적합성 |
| 국소 ROI·경계 제외 분석 | 보류 | V1 수술 AP 확정, bregma–CCF 변환, 광분포·공간 자료 |
| 지역 차이에 따른 최종 후보 3–5개 | 미선정 | 대응 세포 유형 안의 실제 차이와 자료 간 확인 |

RNA-seq 코드에서는 exact MOp와 VISp만 비교하며 VIS 전체·VIS-PTLp·ALM으로 대체하지 않는다. 평균 정규화 발현·검출률과 같은 subclass 안의 차이를 기술적으로 제시하고, 공유된 세부 cluster의 공통 가중치 결과를 별도로 낸다. 독립 동물 정보가 확인되지 않은 세포 수로 p-value를 계산하지 않는다. 비뉴런 표집이 부족할 경우 그 범위를 제한한다.

B6J 공간 자료는 단일 동물이다. 전체 QC·분류 세포 대비와 뉴런·GABAergic 세포 대비의 분모를 구분하고, raw taxonomy label과 표시용 매핑을 보존한다. 면적 마스크가 없으므로 density나 전체 뇌 세포 수를 산출하지 않는다. 해부학적 layer와 transcriptomic subclass를 별도 열로 둔다.

## 자극 정보

사용자가 확인한 값은 mm 단위 bregma 기준 AP/ML, 뇌 표면 기준 DV와 광섬유 core 직경 400 µm이다. M1은 `(1.34, 2.46, 0)`, V1은 `(약 −2.3, 2.46, 0)`이며 V1 AP는 기억에 의존한 잠정값이다. bregma 좌표를 Allen CCF 좌표로 직접 대입하지 않는다. ML의 반구 부호도 실제 수술 기록과 atlas 규약을 확인할 때 해석한다. core 직경만으로 광학적 자극 범위를 정하지 않는다.

## 다운로드 차단과 재개 조건

현재 실행 환경에서 Allen S3·포털·NCBI의 HTTPS 요청이 TLS 연결 전 proxy CONNECT 403으로 거부됐다. 이를 자료 없음이나 인증키 누락으로 해석하지 않았다. 필요한 도메인을 클라우드 환경 설정 초안에 추가해 저장했다. 이 저장은 런타임 접근 허용 또는 게시 완료를 뜻하지 않는다.

환경 설정에서 네트워크 추가를 검토·저장하고 안내에 따라 게시한 뒤 실제 접근을 확인하면 다운로드와 미완료 분석을 재개할 수 있다. 로컬 Windows에서 공식 자료를 다운로드하거나 보유한 입력을 이 폴더 구조로 제공해도 된다. RNA-seq 원본의 정확한 다운로드 URL·단위·독립 표본 수는 아직 확인되지 않아 추정값으로 채우지 않았다.

## 결과와 재현

- [PBM 연결 보고서](results/pbm/PBM_bridge_report_ko.md), [주요 표 Excel](results/pbm/PBM_panel22_and_gene_sets.xlsx), `results/pbm/validation_summary.json`
- `results/coverage/panel22_measurement_coverage.csv`, `results/coverage/exploration62_measurement_coverage.csv`, `results/coverage/validation_summary.json`
- [실행·재개 방법](README_ko.md), [파일 관리·Windows 이동](docs/file_management_ko.md), [사용자 자극 좌표](docs/stimulation_coordinates_user.json)

PBM 재실행의 CSV·manifest·보고서·검증 결과 해시 일치와 Excel 재독해를 확인했다. 합성 입력을 사용하는 Allen 코드 테스트는 코드의 동작 확인이며 실제 Allen 데이터 분석 완료가 아니다. 전체 프로젝트 ZIP은 원본 입력·코드·보고서·파일별 해시를 포함하고 OS별 가상환경·임시 로그는 제외한다.
