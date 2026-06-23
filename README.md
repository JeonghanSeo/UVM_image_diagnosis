# UVM Image Diagnosis

Uveal Melanoma (포도막 흑색종) 멀티모달 예측 모델 개발 프로젝트

## Overview

TCGA-UVM 코호트 (80명) 의 병리 슬라이드 (WSI), 임상 데이터, 유전체 데이터를 활용한 멀티모달 딥러닝 모델.

**예측 목표:**
- 분자 서브타입 (SCNA Cluster 1~4)
- 생존율 / 재발 위험도 (OS, PFI)

## Molecular Subtypes (Robertson et al. Cancer Cell 2017)

| Cluster | Chr3 | Key Mutation | Old Class | Prognosis |
|---------|------|-------------|-----------|-----------|
| Cluster 1 | Disomy 3 | EIF1AX | Class 1A | Best |
| Cluster 2 | Disomy 3 | SF3B1  | Class 1B | Very Good |
| Cluster 3 | Monosomy 3 | BAP1 | Class 2 | Poor |
| Cluster 4 | Monosomy 3 | BAP1 + iso8q | Class 2 | Worst |

## Project Structure

```
UVM_image_diagnosis/
├── 01.Raw_Data/
│   ├── 01.WSI/          # 병리 슬라이드 (SVS) — TCIA (git 제외)
│   ├── 02.Clinical/     # 임상 데이터 CSV — GDC Portal
│   └── 03.Genomics/     # RNA-seq — GDC Portal (git 제외)
├── 02.Preprocessing/
│   ├── 01.Tiles/        # WSI 타일 패치 (git 제외)
│   └── 02.Features/     # 추출된 feature 벡터 (git 제외)
├── 03.EDA/              # 탐색적 데이터 분석 스크립트
├── 04.Models/           # 모델 아키텍처
├── 05.Training/
│   ├── 01.Checkpoints/  # 학습 가중치 (git 제외)
│   └── 02.Logs/         # 학습 로그
├── 06.Results/          # 평가 결과
├── 07.Figures/          # 논문/보고서용 그래프
└── 99.Scripts/          # 유틸리티 스크립트
```

## Data Sources

| Data | Source | Note |
|------|--------|------|
| WSI (SVS) | [TCIA](https://www.cancerimagingarchive.net) | TCGA-UVM collection |
| Clinical / RNA-seq | [GDC Portal](https://portal.gdc.cancer.gov) | TCGA-UVM project |
| Molecular subtypes | Robertson et al. 2017 Table S1 | SCNA Cluster 1~4 |
| Survival data | TCGA-CDR (Liu et al. Cell 2018) | OS, PFI, DSS, DFI |

## Environment

```bash
conda activate uvm
# Python 3.10 / PyTorch 2.12.1 (CPU) / lifelines 0.30.0
```

## Status

- [x] 임상 데이터 수집 및 마스터 테이블 구성
- [x] EDA (Cluster 분포, KM 생존 곡선, 돌연변이 landscape)
- [x] RNA-seq 다운로드 및 행렬 구성
- [ ] WSI 다운로드 (집 GPU PC)
- [ ] WSI 전처리 (타일링, feature 추출)
- [ ] 멀티모달 모델 구현

## References

1. **Robertson et al. (2017)**
   Integrative Analysis Identifies Four Molecular and Clinical Subsets in Uveal Melanoma.
   *Cancer Cell*, 32(2), 204–220.
   https://doi.org/10.1016/j.ccell.2017.07.003
   > TCGA-UVM 코호트 주 논문. SCNA Cluster 1~4 분자 서브타입 정의, Table S1 (환자별 서브타입 레이블) 출처.

2. **Liu et al. (2018)**
   An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics.
   *Cell*, 173(2), 400–416.
   https://doi.org/10.1016/j.cell.2018.02.052
   > TCGA-CDR. OS / DSS / DFI / PFI 표준화 생존 데이터 출처.

3. **TCGA Research Network (2013)**
   The Cancer Genome Atlas Pan-Cancer analysis project.
   *Nature Genetics*, 45, 1113–1120.
   https://doi.org/10.1038/ng.2764
   > TCGA 전체 프로젝트 기준 논문.

4. **GDC (Grossman et al., 2016)**
   Toward a Shared Vision for Cancer Genomic Data.
   *New England Journal of Medicine*, 375, 1109–1112.
   https://doi.org/10.1056/NEJMp1607591
   > GDC Portal (portal.gdc.cancer.gov) 데이터 다운로드 플랫폼.

5. **TCIA (Clark et al., 2013)**
   The Cancer Imaging Archive (TCIA): Maintaining and Operating a Public Information Repository.
   *Journal of Digital Imaging*, 26, 1045–1057.
   https://doi.org/10.1007/s10278-013-9622-7
   > TCIA (cancerimagingarchive.net) WSI 이미지 다운로드 플랫폼.
