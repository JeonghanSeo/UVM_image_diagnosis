# UVM Image Diagnosis - Project Structure

Uveal Melanoma (포도막 흑색종) 멀티모달 예측 모델 개발 프로젝트

## Folder Structure

```
UVM_image_diagnosis/
│
├── 01.Raw_Data/
│   ├── 01.WSI/          # 병리 슬라이드 원본 (SVS, TIFF) - TCIA 다운로드
│   ├── 02.Clinical/     # 임상 데이터 CSV (생존기간, 나이, 성별 등) - GDC 다운로드
│   └── 03.Genomics/     # RNA-seq, 돌연변이 데이터 - GDC 다운로드
│
├── 02.Preprocessing/
│   ├── 01.Tiles/        # WSI에서 잘라낸 타일 이미지 (패치)
│   └── 02.Features/     # 타일에서 추출한 feature 벡터
│
├── 03.EDA/              # 탐색적 데이터 분석 노트북 및 결과물
│
├── 04.Models/           # 모델 아키텍처 정의 (Python 파일)
│
├── 05.Training/
│   ├── 01.Checkpoints/  # 학습된 모델 가중치 (.pth)
│   └── 02.Logs/         # 학습 로그 (loss, accuracy 등)
│
├── 06.Results/          # 평가 결과, 예측값, 메트릭
│
├── 07.Figures/          # 논문/보고서용 그래프, 시각화
│
└── 99.Scripts/          # 실행 스크립트, 노트북, 이 문서
```

## Data Sources

| 데이터 | 출처 | 접근 방법 |
|---|---|---|
| WSI (병리 슬라이드) | TCIA | https://www.cancerimagingarchive.net |
| 임상 데이터 / RNA-seq | GDC Portal | https://portal.gdc.cancer.gov |

## Environment

- Python 3.10
- PyTorch 2.12.1 (CPU)
- conda env: `uvm`

## Prediction Targets

- **1단계:** 생존율/재발 위험도 + 분자 서브타입 (Class 1A / 1B / 2)
- **2단계:** 약물 반응성 예측
- **3단계:** 신약 후보 탐색
