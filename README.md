# Floor-Flatness-Analyzer_V3.2
# 🏢 건물 바닥 슬래브 평활도 측정 프로그램 v3.2

건물 바닥 슬래브의 평활도를 자동으로 측정하고 분석하는 Python 기반 프로그램입니다. 
LiDAR 포인트 클라우드 데이터(PLY 파일)를 입력받아 KCS 41 46 01 기준에 따른 등급 판정과 크랙 탐지를 수행합니다.

---

## 📋 주요 기능

### 1. **자동 좌표계 보정**
- 스캔 데이터의 기울기 자동 감지
- 회전 행렬 계산 및 좌표 자동 보정
- 다양한 각도로 스캔된 데이터에 대응

### 2. **강한 필터링 (Smart Floor Detection)**
- 스캔 최저 Z 기준 상대값 기반 필터링
- 명확한 바닥 영역만 추출
- 노이즈 및 불필요한 데이터 제거

### 3. **평활도 분석**
- 최소자승법(LSQ) 평면 피팅
- 편차 계산 및 자동 분류:
  - 🟢 평탄: -2mm ~ +2mm
  - 🟡 함몰: < -2mm
  - 🔴 요철: > +2mm

### 4. **KCS 41 46 01 기준 등급 판정**
- **A 우수**: 초과율 ≤ 5%
- **B 양호**: 초과율 ≤ 15%
- **C 보통**: 초과율 ≤ 30%
- **D 불량**: 초과율 > 30%

### 5. **다중신호 크랙 탐지**
- 곡률 기반 탐지
- 법선 벡터 불연속성 감지
- 이방성(Anisotropy) 분석
- 통합 크랙 점수 계산

### 6. **3m 직선자 기준 검증**
- 실제 측정 데이터 기반 직선자 시뮬레이션
- 3.0mm 기준 자동 판정
- 참고용 추가 검증

### 7. **다양한 시각화 및 리포트**
- 히트맵(실측 기반, 보간 없음)
- 3D 포인트 클라우드 시각화
- 편차 분포 그래프
- 상세 텍스트 리포트
- JSON 형식 데이터 리포트

---

## 🛠️ 설치 및 환경 설정

### 필수 패키지
```bash
pip install numpy scipy matplotlib pathlib
```

### 권장 Python 버전
- Python 3.7 이상

---

## 🚀 사용 방법

### 1. 기본 사용법

```python
from flatness_analyzer_v3_2_KCS_CrackDetection import main

# PLY 파일 경로와 출력 디렉토리 지정
ply_file = '/path/to/your/scan.ply'
output_dir = '/path/to/output'

# 분석 실행
result = main(ply_file, output_dir, z_threshold=0.5)
```

### 2. 스캔 미리보기 (ROI 좌표 확인)

특정 구역(예: 테스트보드)만 분석하려면, 먼저 전체 스캔의 위에서 내려다본 평면도를 확인하세요.

```python
from flatness_analyzer_v3_2_KCS_CrackDetection import preview_scan

# 미리보기 생성
roi_info = preview_scan(ply_file, output_dir, z_threshold=0.5)

# 출력에서 테스트보드가 차지하는 X/Y 범위 확인
# preview_topdown.png를 열어서 좌표 확인
print(f"X 범위: {roi_info['x_range']}")
print(f"Y 범위: {roi_info['y_range']}")
```

### 3. ROI(관심영역) 지정하여 분석

```python
# 테스트보드의 좌표 범위 지정
roi_bounds = (x_min, x_max, y_min, y_max)  # 예: (0.5, 2.5, 1.0, 3.0)

result = main(ply_file, output_dir, z_threshold=0.5, roi_bounds=roi_bounds)
```

### 4. 명령줄 실행

```bash
python flatness_analyzer_v3_2_KCS_CrackDetection.py
```
> ⚠️ 코드 내 `__main__` 섹션의 경로를 수정한 후 실행하세요.

---

## 📊 입출력 파일

### 입력 파일
- **PLY 파일** (포인트 클라우드)
  - 형식: `binary_little_endian` 또는 `binary_big_endian`
  - 필수 속성: X, Y, Z 좌표 (float)
  - 선택 속성: RGB, 강도 등

### 출력 파일

분석 완료 후 다음 파일들이 생성됩니다:

| 파일명 | 설명 |
|--------|------|
| `report.txt` | 상세 텍스트 리포트 (한글) |
| `report.json` | 구조화된 JSON 데이터 |
| `heatmap.png` | 편차 히트맵 (평탄/함몰/요철) |
| `additional_graphs.png` | 편차 분포, 통계 그래프 |
| `pointcloud_3d.pcd` | 3D 포인트 클라우드 (PCD 형식) |

---

## 🎯 파라미터 설명

### `main()` 함수 파라미터

| 파라미터 | 타입 | 기본값 | 설명 |
|---------|------|--------|------|
| `ply_filepath` | str | - | PLY 파일 경로 (필수) |
| `output_dir` | str | `~/Desktop/flatness_output` | 출력 디렉토리 |
| `z_threshold` | float | 0.5 | 필터링 Z 범위 (m) |
| `roi_bounds` | tuple | None | ROI 좌표 (x_min, x_max, y_min, y_max) |

### 상수 설정 (코드 수정)

```python
# KCS 기준: 초과율에 따른 등급
GRADE_BOUNDARIES = {
    'A': 0.05,   # 우수
    'B': 0.15,   # 양호
    'C': 0.30,   # 보통
    'D': 1.0     # 불량
}

# 크랙 탐지 설정
THRESHOLD_CRACK_SCORE = 0.5    # 크랙 점수 임계값 (0~1)
CURVATURE_K = 15               # k-nearest neighbor 개수
CRACK_THRESHOLD_MM = 3.0       # 크랙 판정 편차 (mm)
```

---

## 📈 분석 프로세스

```
┌─────────────────────────────────────────┐
│     Phase 1: 데이터 전처리              │
│  ├─ PLY 파일 읽기                      │
│  ├─ 강한 필터링 (바닥 추출)             │
│  ├─ ROI 크롭 (선택)                    │
│  ├─ 좌표계 기울기 감지                  │
│  └─ 좌표 자동 보정                      │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│     Phase 2: 평활도 계산                │
│  ├─ RANSAC 평면 검출                   │
│  ├─ LSQ 평면 피팅                      │
│  ├─ 편차 계산                          │
│  └─ 카테고리 분류                      │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│     Phase 3: 검증 및 탐지               │
│  ├─ 직선자 시뮬레이션                   │
│  ├─ 다중신호 크랙 탐지                  │
│  └─ KCS 등급 판정                      │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│     Phase 4: 시각화                    │
│  ├─ 히트맵 생성                        │
│  ├─ 3D 포인트 클라우드 렌더링           │
│  └─ 분포 그래프 생성                    │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│     Phase 5: 리포트 생성                │
│  ├─ 텍스트 리포트 작성                  │
│  └─ JSON 데이터 내보내기                │
└─────────────────────────────────────────┘
```

---

## 📊 출력 리포트 예시

### 텍스트 리포트 (report.txt)
```
================================================================================
🏢 건물 바닥 슬래브 평활도 측정 프로그램 v3.2
================================================================================

【 1. 스캔 정보 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
파일:              S8동_복도_바닥.ply
포인트 수:         125,436
스캔 시각:         2024-01-15 14:30:25

【 2. 데이터 필터링 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
원본:              125,436 → 필터 후: 95,234 (유지율: 75.9%)
Z 필터링 범위:     -0.500m ~ 0.000m

【 3. 자동 좌표계 보정 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
기울기:            2.45° (자동 보정됨)
법선 벡터:         [0.001, -0.043, 0.999]

【 4. 평활도 통계 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
평균 편차:         +0.12mm
표준편차:          1.84mm
RMS:               1.95mm
PV (P-V):          8.56mm

【 5. 카테고리별 분석 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【 초록색 (평탄) 】 -2mm ~ +2mm
  포인트: 87,234 (91.7%)
  평균 편차: -0.05mm
  판정: ✅ 합격

【 노란색 (함몰) 】 < -2mm
  포인트: 3,456 (3.6%)
  평균 편차: -3.24mm
  판정: ⚠️  주의

【 빨강색 (요철) 】 > +2mm
  포인트: 4,544 (4.8%)
  평균 편차: +3.78mm
  판정: ❌ 부적합

【 6. 크랙 탐지 (다중신호 기반) 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
크랙 포인트:       2,145 / 95,234 (2.25%)
크랙 점수 평균:    0.573
크랙 점수 최댓값:  0.892
탐지 신호:        곡률 + 법선불연속성 + 이방성

【 7. KCS 41 46 01 등급 판정 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
등급:              B 양호
기준값:            ±2.0mm
초과 포인트:       8,000 / 95,234 (8.40%)

【 8. 3m 직선자 기준 검증 (참고) 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
기준: 3.0mm 이하
테스트 직선: 48개
최대 틈새: 2.84mm
판정: ✅ 합격

【 9. 최종 판정 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
KCS 등급:          B 양호
3m 직선자 기준:    ✅ 만족 (참고용)
```

---

## 🔍 각 분석 단계의 상세 설명

### Phase 1: 데이터 전처리

#### 1.1 PLY 파일 읽기
- Binary 형식의 PLY 파일을 로드
- 헤더 파싱으로 포인트 개수 및 속성 확인
- X, Y, Z 좌표 추출

#### 1.2 강한 필터링
- 스캔 최저 Z 값을 기준으로 설정 범위 내의 포인트만 추출
- `z_threshold` 파라미터로 조정 가능
- 천장, 벽 등 불필요한 부분 제거

#### 1.3-1.5 자동 좌표계 보정
- PCA(주성분 분석)를 통한 평면 법선 벡터 추정
- 법선 벡터와 수직(Z축)의 각도 계산
- 기울기 >= 0.5°일 때 자동 보정
- Rodrigues 회전 공식으로 회전 행렬 생성

#### 1.6 RANSAC 평면 검출
- 노이즈에 강인한 평면 검출 알고리즘
- 이상치 포인트 식별 및 제거
- 다음 단계의 기준면 제공

### Phase 2: 평활도 계산

#### 2.1 최소자승법(LSQ) 평면 피팅
- RANSAC 결과를 초기값으로 사용
- 모든 인라이어 포인트로 최적 평면 계산
- 각 포인트의 평면으로부터의 거리(편차) 계산

#### 2.2-2.4 편차 계산 및 분류
- 각 포인트의 편차를 mm 단위로 변환
- 편차 범위에 따라 자동 분류:
  - **FLAT (평탄)**: -2 ~ +2 mm
  - **DEPRESSION (함몰)**: < -2 mm
  - **PROTRUSION (요철)**: > +2 mm

### Phase 3: 검증 및 탐지

#### 3.1 직선자 시뮬레이션
- 바닥 영역을 격자로 나누어 여러 직선 생성
- 각 직선에서 최대 편차(틈새) 측정
- 3.0mm 기준 만족 여부 판정

#### 3.2 다중신호 크랙 탐지
세 가지 신호를 조합하여 크랙 탐지:

1. **곡률 기반**: k-nearest neighbor의 곡률 계산
2. **법선 불연속성**: 인접 포인트의 법선 각도 변화
3. **이방성**: 주방향의 편차 분포 불균형

최종 크랙 점수 = 정규화된 세 신호의 가중 평균

#### 3.3 KCS 41 46 01 등급 판정
- 편차 허용 범위(±2.0mm)를 초과하는 포인트의 비율 계산
- 초과율에 따라 A~D 등급 자동 판정

### Phase 4-5: 시각화 및 리포트

#### 히트맵
- 보간 없는 실제 측정 포인트 기반
- 색상으로 편차 범주 표시 (초록/노랑/빨강)

#### 3D 포인트 클라우드
- PCD 형식으로 저장
- CloudCompare 등의 도구로 시각화 가능

#### 그래프
- 편차 분포 히스토그램
- 누적 분포 함수
- 카테고리별 통계 막대 그래프

---

## ⚠️ 주의사항 및 팁

### 데이터 전처리
- **z_threshold 설정**
  - 너무 작으면 바닥 포인트 부족 → 에러 발생
  - 너무 크면 천장/벽 포함 → 부정확한 결과
  - 권장: 0.3 ~ 1.0m

- **PLY 파일 형식**
  - Binary 형식만 지원 (ASCII 미지원)
  - float 형식의 X, Y, Z 필수
  - 인코딩: UTF-8 또는 latin1

### ROI 지정
- 먼저 `preview_scan()`으로 미리보기 생성
- PNG 이미지에서 테스트보드의 좌표 범위 확인
- `(x_min, x_max, y_min, y_max)` 형식으로 입력

### 크랙 탐지 임계값 조정
- `THRESHOLD_CRACK_SCORE`: 0 ~ 1 사이의 값
  - 낮을수록 더 많은 포인트를 크랙으로 탐지
  - 기본값(0.5): 중간 수준의 민감도
  
- `CRACK_THRESHOLD_MM`: 편차 기준 (mm)
  - 기본값(3.0mm)로 충분한 대부분의 경우

---

## 📖 API 사용 예시

### 예시 1: 기본 분석

```python
from flatness_analyzer_v3_2_KCS_CrackDetection import main

result = main(
    ply_filepath='/path/to/scan.ply',
    output_dir='/path/to/output',
    z_threshold=0.5
)

# 결과 접근
print(f"KCS 등급: {result['kcs_results']['grade_name']}")
print(f"평탄 비율: {result['stats']['FLAT']['percentage']:.1f}%")
print(f"크랙 비율: {result['crack_results']['crack_percentage']:.2f}%")
```

### 예시 2: 특정 영역만 분석

```python
from flatness_analyzer_v3_2_KCS_CrackDetection import preview_scan, main

# 1단계: 미리보기로 좌표 확인
roi_info = preview_scan('/path/to/scan.ply', '/path/to/output')

# 2단계: ROI 지정하여 분석
result = main(
    ply_filepath='/path/to/scan.ply',
    output_dir='/path/to/output',
    z_threshold=0.5,
    roi_bounds=(0.5, 2.5, 1.0, 3.0)  # 테스트보드 좌표
)
```

### 예시 3: 배치 처리

```python
import os
from pathlib import Path
from flatness_analyzer_v3_2_KCS_CrackDetection import main

scan_dir = '/path/to/ply/files'
output_base = '/path/to/output'

for ply_file in Path(scan_dir).glob('*.ply'):
    output_subdir = f"{output_base}/{ply_file.stem}"
    
    try:
        result = main(str(ply_file), output_subdir)
        print(f"✅ {ply_file.name}: {result['kcs_results']['grade_name']}")
    except Exception as e:
        print(f"❌ {ply_file.name}: {e}")
```

---

## 🐛 트러블슈팅

### 문제 1: "지원하지 않는 PLY 포맷입니다" 에러
**원인**: ASCII 형식의 PLY 파일 또는 지원하지 않는 바이너리 형식  
**해결**: Binary Little Endian 형식으로 변환 필요

### 문제 2: "Z < ... 조건을 만족하는 바닥 포인트가 없습니다" 에러
**원인**: `z_threshold` 값이 너무 작음  
**해결**: `z_threshold` 값을 증가시켜 다시 실행
```python
result = main(ply_file, output_dir, z_threshold=1.0)  # 기본값 0.5에서 1.0으로 증가
```

### 문제 3: 결과가 부정확하거나 이상함
**원인**: PLY 파일에 노이즈가 많거나 기울기가 심함  
**해결**:
- 먼저 `preview_scan()`으로 데이터 확인
- `z_threshold` 값 조정
- ROI 지정으로 분석 영역 제한

### 문제 4: 메모리 부족 에러
**원인**: 매우 큰 포인트 클라우드 (수백만 개 포인트)  
**해결**:
- ROI를 지정하여 분석 영역 축소
- 시스템 메모리 확인
- 64비트 Python 사용 확인

---

## 📝 파일 형식 상세

### PLY 헤더 예시
```
ply
format binary_little_endian 1.0
element vertex 125436
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
end_header
[바이너리 데이터...]
```

### JSON 리포트 구조
```json
{
  "metadata": {
    "filename": "scan.ply",
    "vertex_count": 125436,
    "scan_time": "2024-01-15 14:30:25",
    "x_range": [0.0, 5.0],
    "y_range": [0.0, 4.0],
    "z_range": [-2.5, 0.5]
  },
  "kcs_results": {
    "grade_name": "B",
    "tolerance_mm": 2.0,
    "exceed_count": 8000,
    "exceed_ratio": 0.084
  },
  "stats": {
    "FLAT": {"count": 87234, "percentage": 91.7},
    "DEPRESSION": {"count": 3456, "percentage": 3.6},
    "PROTRUSION": {"count": 4544, "percentage": 4.8}
  },
  "crack_results": {
    "crack_count": 2145,
    "crack_percentage": 2.25,
    "crack_scores": [0.45, 0.62, ...]
  }
}
```

---

## 📚 참고 자료

- **KCS 41 46 01**: 콘크리트 바닥 평활도 기준
- **RANSAC**: 노이즈에 강인한 평면 검출 알고리즘
- **최소자승법(LSQ)**: 최적 평면 피팅 방법

---

## 📞 기술 지원

프로그램 사용 중 문제가 발생하면:
1. 위 "트러블슈팅" 섹션 확인
2. 콘솔 출력 메시지 및 에러 확인
3. `preview_scan()`으로 입력 데이터 검증

---

## 📄 라이센스

본 프로그램은 건축 현장의 바닥 품질 검사 자동화를 위해 개발되었습니다.

**작성일**: 2024년 1월  
**버전**: 3.2  
**주요 개선사항**:
- v3.0: KCS 기준 등급 판정, 자동 좌표계 보정
- v3.1: 다중신호 크랙 탐지 추가
- v3.2: 강한 필터링 개선, 안정성 향상

---

**마지막 업데이트**: 2024-01-15
