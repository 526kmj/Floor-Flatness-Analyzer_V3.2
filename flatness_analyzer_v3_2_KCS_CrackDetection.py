#!/usr/bin/env python3
"""
건물 바닥 슬래브 평활도 측정 프로그램 v3.2 (KCS 기준 + 크랙 탐지)
Floor Slab Flatness Measurement Algorithm v3.2

핵심 기능:
- 자동 좌표계 기울기 감지
- 자동 회전 행렬 계산 및 좌표 보정
- 강한 필터링 (스캔 최저 Z 기준 상대값)
- 다중신호 크랙 탐지 (곡률, 법선 불연속성, 이방성)
- KCS 41 46 01 기준 자동 등급 판정 (A/B/C/D)
- 편차 기반 직선자 검증 (참고용)
- 다양한 기울기의 PLY 파일 자동 처리
- 보간 없는 실측 기반 히트맵
"""

import numpy as np
import json
import os
from datetime import datetime
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# ============================================================================
# 판정 기준 상수 (KCS 41 46 01)
# ============================================================================

# KCS 41 46 01 기준: 초과율(초과 포인트 비율)에 따른 등급
GRADE_BOUNDARIES = {
    'A': 0.05,      # A등급 (우수): 초과율 ≤ 5%
    'B': 0.15,      # B등급 (양호): 초과율 ≤ 15%
    'C': 0.30,      # C등급 (보통): 초과율 ≤ 30%
    'D': 1.0        # D등급 (불량): 초과율 > 30%
}

# 크랙 탐지 임계값
THRESHOLD_CRACK_SCORE = 0.5    # 크랙 점수 임계값 (0~1)
CURVATURE_K = 15               # k-nearest neighbor 개수
CRACK_THRESHOLD_MM = 3.0       # 크랙으로 판정하는 편차 (mm)

# ============================================================================
# 1. PLY 파일 읽기
# ============================================================================

def read_ply_file(filepath):
    """PLY 파일을 읽고 포인트 클라우드 데이터를 반환"""
    print(f"📁 PLY 파일 읽는 중: {filepath}")
    
    with open(filepath, 'rb') as f:
        header_lines = []
        while True:
            line = f.readline().decode('latin1')
            header_lines.append(line)
            if line.strip() == 'end_header':
                break

        format_type = None
        vertex_count = 0
        properties = []
        current_element = None

        for line in header_lines:
            stripped = line.strip()
            if stripped.startswith('format'):
                format_type = stripped.split()[1]
            elif stripped.startswith('element'):
                parts = stripped.split()
                current_element = parts[1]
                if current_element == 'vertex':
                    vertex_count = int(parts[2])
            elif stripped.startswith('property') and current_element == 'vertex':
                parts = stripped.split()
                if len(parts) >= 3:
                    properties.append((parts[1], parts[2]))

        print(f"   포인트 수: {vertex_count:,}")
        print(f"   포맷: {format_type}")

        if format_type == 'binary_little_endian':
            byte_order = '<'
        elif format_type == 'binary_big_endian':
            byte_order = '>'
        else:
            raise ValueError(
                f"지원하지 않는 PLY 포맷입니다: {format_type} "
                f"(binary_little_endian 또는 binary_big_endian만 지원)"
            )

        type_map = {'float': 'f4', 'uchar': 'u1'}
        dtype_fields = []
        for i, (ptype, pname) in enumerate(properties):
            if ptype not in type_map:
                raise ValueError(f"지원하지 않는 property 타입: {ptype} ({pname})")
            dtype_fields.append((f'f{i}', f'{byte_order}{type_map[ptype]}'))

        raw = np.fromfile(f, dtype=np.dtype(dtype_fields), count=vertex_count)
        data = np.column_stack([raw[f'f{i}'] for i in range(len(properties))])

    data = np.array(data, dtype=np.float32)
    
    metadata = {
        'filename': Path(filepath).name,
        'vertex_count': vertex_count,
        'scan_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'x_range': (float(data[:, 0].min()), float(data[:, 0].max())),
        'y_range': (float(data[:, 1].min()), float(data[:, 1].max())),
        'z_range': (float(data[:, 2].min()), float(data[:, 2].max())),
    }
    
    print(f"✅ 읽기 완료")
    print(f"   X: {metadata['x_range'][0]:.3f} ~ {metadata['x_range'][1]:.3f}m")
    print(f"   Y: {metadata['y_range'][0]:.3f} ~ {metadata['y_range'][1]:.3f}m")
    print(f"   Z: {metadata['z_range'][0]:.3f} ~ {metadata['z_range'][1]:.3f}m")
    
    return data, metadata

# ============================================================================
# 2. 강한 필터링 - 바닥만 추출
# ============================================================================

def strong_floor_filtering(points, z_threshold=0.5):
    """강한 필터링으로 명확한 바닥 영역만 추출"""
    print("\n【 Step 1.2: 강한 필터링 - 명확한 바닥 추출 】")
    
    xyz = points[:, :3]
    original_count = len(points)

    if original_count == 0:
        raise ValueError("입력 포인트 클라우드가 비어 있습니다.")

    # 바닥은 스캔에서 가장 낮은 평면이라고 가정하고, 최저 Z 기준 상대값으로 필터링
    z_min = xyz[:, 2].min()
    z_cutoff = z_min + z_threshold
    mask = xyz[:, 2] < z_cutoff
    floor_candidates = points[mask]

    print(f"   필터링 기준: Z < (최저 Z {z_min:.3f}m + {z_threshold}m) = {z_cutoff:.3f}m")
    print(f"   원본: {original_count:,} → 필터 후: {len(floor_candidates):,}")

    if len(floor_candidates) == 0:
        z_max = xyz[:, 2].max()
        raise ValueError(
            f"Z < {z_cutoff:.3f}m 조건을 만족하는 바닥 포인트가 없습니다 "
            f"(이 스캔의 Z 범위: {z_min:.3f} ~ {z_max:.3f}m). "
            f"Z 임계값을 조정해서 다시 시도해주세요."
        )

    print(f"   유지율: {len(floor_candidates) / original_count * 100:.1f}%")
    
    stats = {
        'original': original_count,
        'after_strong_filter': len(floor_candidates),
        'retention_rate': len(floor_candidates) / original_count * 100,
        'z_threshold': z_threshold,
    }
    
    return floor_candidates, stats

# ============================================================================
# 2-1. 스캔 미리보기 (ROI 좌표 확인용)
# ============================================================================

def preview_scan(ply_filepath, output_dir, z_threshold=0.5, grid_interval=0.1):
    """분석 전 스캔 데이터를 위에서 내려다본 평면도로 저장 (테스트보드 등 ROI 좌표 파악용)"""
    print("\n" + "=" * 80)
    print("📐 스캔 미리보기 (ROI 좌표 확인용)")
    print("=" * 80)

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    points, _ = read_ply_file(ply_filepath)
    floor_candidates, _ = strong_floor_filtering(points, z_threshold=z_threshold)
    xyz = floor_candidates[:, :3]

    x_min, x_max = float(xyz[:, 0].min()), float(xyz[:, 0].max())
    y_min, y_max = float(xyz[:, 1].min()), float(xyz[:, 1].max())

    print(f"   전체 바닥 X 범위: {x_min:.3f} ~ {x_max:.3f}m")
    print(f"   전체 바닥 Y 범위: {y_min:.3f} ~ {y_max:.3f}m")

    fig, ax = plt.subplots(figsize=(10, 10))
    sc = ax.scatter(xyz[:, 0], xyz[:, 1], c=xyz[:, 2], cmap='viridis', s=1)
    plt.colorbar(sc, ax=ax, label='Z (m)', shrink=0.8)

    ax.set_xticks(np.arange(np.floor(x_min / grid_interval) * grid_interval, x_max + grid_interval, grid_interval))
    ax.set_yticks(np.arange(np.floor(y_min / grid_interval) * grid_interval, y_max + grid_interval, grid_interval))
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.5)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title('Top-down Preview (테스트보드 ROI 좌표 확인용)')
    ax.set_aspect('equal')

    plt.tight_layout()
    preview_path = f'{output_dir}/preview_topdown.png'
    plt.savefig(preview_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"   ✅ 저장: {preview_path}")
    print(f"   → 이 이미지에서 테스트보드가 차지하는 X/Y 범위를 확인한 뒤,")
    print(f"     main(..., roi_bounds=(x_min, x_max, y_min, y_max)) 형태로 지정해서 분석하세요.")

    return {
        'x_range': (x_min, x_max),
        'y_range': (y_min, y_max),
    }

# ============================================================================
# 2-2. ROI(관심영역) 크롭 - 테스트보드 등 특정 구역만 분석
# ============================================================================

def crop_to_roi(points, roi_bounds):
    """지정한 X/Y 범위 밖의 점을 제거 (예: 테스트보드 영역만 남기고 주변 바닥 제외)"""
    print("\n【 Step 1.25: ROI 크롭 (지정 구역만 분석) 】")

    x_min, x_max, y_min, y_max = roi_bounds
    xyz = points[:, :3]

    mask = (xyz[:, 0] >= x_min) & (xyz[:, 0] <= x_max) & \
           (xyz[:, 1] >= y_min) & (xyz[:, 1] <= y_max)
    cropped = points[mask]

    print(f"   ROI 범위: X [{x_min:.3f}, {x_max:.3f}]m, Y [{y_min:.3f}, {y_max:.3f}]m")
    print(f"   크롭 전: {len(points):,} → 크롭 후: {len(cropped):,}")

    if len(cropped) == 0:
        raise ValueError(
            f"지정한 ROI 범위(X: {x_min}~{x_max}, Y: {y_min}~{y_max})에 포인트가 없습니다. "
            f"preview_scan()으로 좌표를 다시 확인해주세요."
        )

    return cropped

# ============================================================================
# 3. 자동 좌표계 기울기 감지 및 보정
# ============================================================================

def detect_tilt_and_get_normal(points):
    """
    바닥 평면의 기울기를 자동으로 감지하고 법선벡터 반환
    
    Returns:
        tuple: (법선벡터, 기울기각도)
    """
    print("\n【 Step 1.3: 자동 좌표계 기울기 감지 】")
    
    xyz = points[:, :3]
    
    # 최소제곱 평면 피팅
    X = xyz[:, :2]
    Z = xyz[:, 2]
    ones = np.ones((len(xyz), 1))
    A = np.hstack([X, ones])
    
    coeffs, _, _, _ = np.linalg.lstsq(A, Z, rcond=None)
    a, b, d = -coeffs[0], -coeffs[1], -coeffs[2]
    c = 1.0

    # 정규화 (단위 법선벡터, Z축 위쪽을 향하도록)
    norm = np.sqrt(a**2 + b**2 + c**2)
    normal = np.array([a / norm, b / norm, c / norm])
    
    # 기울기 각도 계산
    ideal_z = np.array([0, 0, 1])
    dot_product = np.dot(normal, ideal_z)
    angle_rad = np.arccos(np.clip(dot_product, -1, 1))
    angle_deg = float(np.degrees(angle_rad))
    
    print(f"   현재 법선벡터: ({normal[0]:.6f}, {normal[1]:.6f}, {normal[2]:.6f})")
    print(f"   기울기 각도: {angle_deg:.2f}°")
    
    if angle_deg < 0.5:
        print(f"   ✅ 판정: 좌표계 정상 (보정 불필요)")
        needs_correction = False
    else:
        print(f"   ⚠️  판정: 기울어짐 감지 (자동 보정 진행)")
        needs_correction = True
    
    return normal, angle_deg, needs_correction

# ============================================================================
# 4. 회전 행렬 계산 및 좌표 보정
# ============================================================================

def calculate_rotation_matrix(normal_vector):
    """
    현재 기울어진 법선벡터를 Z축(0,0,1)으로 회전시키는 회전 행렬 계산
    
    Args:
        normal_vector: 현재 바닥의 법선벡터 (3D)
    
    Returns:
        np.ndarray: 3x3 회전 행렬
    """
    print("\n【 Step 1.4: 회전 행렬 계산 】")
    
    n_current = normal_vector
    n_target = np.array([0, 0, 1])
    
    # 회전축: 외적
    rotation_axis = np.cross(n_current, n_target)
    
    # 회전축이 0인 경우 처리 (이미 수평)
    if np.allclose(rotation_axis, 0):
        print(f"   이미 수평 상태 → 단위 행렬 반환")
        return np.eye(3)
    
    rotation_axis_normalized = rotation_axis / np.linalg.norm(rotation_axis)
    
    # 회전각
    cos_angle = np.dot(n_current, n_target) / (np.linalg.norm(n_current) * np.linalg.norm(n_target))
    angle_rad = np.arccos(np.clip(cos_angle, -1, 1))
    
    # Rodrigues 회전 공식을 사용한 회전 행렬 계산
    rotation = Rotation.from_rotvec(rotation_axis_normalized * angle_rad)
    rotation_matrix = rotation.as_matrix()
    
    print(f"   회전축: ({rotation_axis_normalized[0]:.6f}, {rotation_axis_normalized[1]:.6f}, {rotation_axis_normalized[2]:.6f})")
    print(f"   회전각: {np.degrees(angle_rad):.2f}°")
    print(f"   ✅ 회전 행렬 계산 완료")
    
    return rotation_matrix

def apply_coordinate_correction(points, rotation_matrix):
    """
    모든 포인트에 회전 행렬을 적용하여 좌표 보정
    
    Args:
        points: (N, 6) 또는 (N, 3) 포인트 배열
        rotation_matrix: 3x3 회전 행렬
    
    Returns:
        np.ndarray: 보정된 포인트 배열
    """
    print("\n【 Step 1.5: 좌표 보정 적용 】")
    
    xyz = points[:, :3]
    
    # 회전 적용 (xyz @ R^T)
    corrected_xyz = xyz @ rotation_matrix.T
    
    # RGB 정보가 있으면 포함
    if points.shape[1] > 3:
        corrected_points = np.hstack([corrected_xyz, points[:, 3:]])
    else:
        corrected_points = corrected_xyz
    
    print(f"   보정된 포인트: {len(corrected_points):,}개")
    print(f"   보정 전 Z 범위: {xyz[:, 2].min():.4f} ~ {xyz[:, 2].max():.4f}m")
    print(f"   보정 후 Z 범위: {corrected_xyz[:, 2].min():.4f} ~ {corrected_xyz[:, 2].max():.4f}m")
    print(f"   ✅ 좌표 보정 완료")
    
    return corrected_points

# ============================================================================
# 5. RANSAC 평면 검출
# ============================================================================

def ransac_plane_detection(points, threshold=0.05, max_iterations=500):
    """RANSAC을 이용한 바닥 평면 검출"""
    print("\n【 Step 1.6: RANSAC 평면 검출 】")
    
    xyz = points[:, :3]
    best_inliers = []
    best_plane = None
    
    for iteration in range(max_iterations):
        sample_indices = np.random.choice(len(xyz), 3, replace=False)
        p1, p2, p3 = xyz[sample_indices]
        
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)
        
        if np.allclose(normal, 0):
            continue
        
        normal = normal / np.linalg.norm(normal)
        d = -np.dot(normal, p1)
        plane = np.append(normal, d)
        
        distances = np.abs(np.dot(xyz, normal) + d)
        inliers = np.where(distances < threshold)[0]
        
        if len(inliers) > len(best_inliers):
            best_inliers = inliers
            best_plane = plane
        
        if (iteration + 1) % 100 == 0:
            print(f"   진행: {iteration + 1}/{max_iterations}")
    
    print(f"   최고 inlier: {len(best_inliers):,} ({len(best_inliers) / len(xyz) * 100:.1f}%)")
    
    return best_plane, best_inliers

# ============================================================================
# 6. 최소제곱 평면 피팅
# ============================================================================

def least_squares_plane_fitting(points, plane_initial, inlier_indices):
    """최소제곱 방식으로 기준 평면 재계산"""
    print("\n【 Step 2.1: 최소제곱 평면 피팅 】")
    
    xyz = points[:, :3]
    inlier_points = xyz[inlier_indices]
    
    X = inlier_points[:, :2]
    Z = inlier_points[:, 2]
    ones = np.ones((len(inlier_points), 1))
    A = np.hstack([X, ones])
    
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A, Z, rcond=None)
        a, b, d = -coeffs[0], -coeffs[1], -coeffs[2]
        c = 1.0

        norm = np.sqrt(a**2 + b**2 + c**2)
        plane_final = np.array([a / norm, b / norm, c / norm, d / norm])
        
    except:
        plane_final = plane_initial
    
    print(f"   기준평면: {plane_final[0]:.6f}x + {plane_final[1]:.6f}y + {plane_final[2]:.6f}z + {plane_final[3]:.6f} = 0")
    
    return plane_final

# ============================================================================
# 7. 편차 계산 및 3색상 분류
# ============================================================================

def calculate_deviations_and_classify(points, plane):
    """기준평면으로부터의 편차 계산 및 3색상 분류"""
    print("\n【 Step 2.2-2.4: 편차 계산 및 3색상 분류 】")
    
    xyz = points[:, :3]
    
    signed_distances = np.dot(xyz, plane[:3]) + plane[3]
    deviations_mm = signed_distances * 1000
    
    FLAT_RANGE = 2.0
    
    flat_mask = (deviations_mm >= -FLAT_RANGE) & (deviations_mm <= FLAT_RANGE)
    depression_mask = deviations_mm < -FLAT_RANGE
    protrusion_mask = deviations_mm > FLAT_RANGE
    
    stats = {
        'FLAT': {
            'count': int(np.sum(flat_mask)),
            'percentage': float(np.sum(flat_mask) / len(xyz) * 100),
            'mean_dev': float(deviations_mm[flat_mask].mean() if np.sum(flat_mask) > 0 else 0),
            'std_dev': float(deviations_mm[flat_mask].std() if np.sum(flat_mask) > 0 else 0),
            'min_dev': float(deviations_mm[flat_mask].min() if np.sum(flat_mask) > 0 else 0),
            'max_dev': float(deviations_mm[flat_mask].max() if np.sum(flat_mask) > 0 else 0),
        },
        'DEPRESSION': {
            'count': int(np.sum(depression_mask)),
            'percentage': float(np.sum(depression_mask) / len(xyz) * 100),
            'mean_dev': float(deviations_mm[depression_mask].mean() if np.sum(depression_mask) > 0 else 0),
            'std_dev': float(deviations_mm[depression_mask].std() if np.sum(depression_mask) > 0 else 0),
            'min_dev': float(deviations_mm[depression_mask].min() if np.sum(depression_mask) > 0 else 0),
            'max_dev': float(deviations_mm[depression_mask].max() if np.sum(depression_mask) > 0 else 0),
        },
        'PROTRUSION': {
            'count': int(np.sum(protrusion_mask)),
            'percentage': float(np.sum(protrusion_mask) / len(xyz) * 100),
            'mean_dev': float(deviations_mm[protrusion_mask].mean() if np.sum(protrusion_mask) > 0 else 0),
            'std_dev': float(deviations_mm[protrusion_mask].std() if np.sum(protrusion_mask) > 0 else 0),
            'min_dev': float(deviations_mm[protrusion_mask].min() if np.sum(protrusion_mask) > 0 else 0),
            'max_dev': float(deviations_mm[protrusion_mask].max() if np.sum(protrusion_mask) > 0 else 0),
        }
    }
    
    print(f"   편차 통계 (전체 {len(xyz):,}개):")
    print(f"      최소: {deviations_mm.min():.2f}mm")
    print(f"      평균: {deviations_mm.mean():.2f}mm")
    print(f"      표준편차: {deviations_mm.std():.2f}mm")
    print(f"      최대: {deviations_mm.max():.2f}mm")
    print(f"      RMS: {np.sqrt(np.mean(deviations_mm**2)):.2f}mm")
    print(f"      PV: {deviations_mm.max() - deviations_mm.min():.2f}mm")
    
    print(f"\n   【 초록색 (평탄) 】 -2~+2mm: {stats['FLAT']['count']:,} ({stats['FLAT']['percentage']:.1f}%)")
    print(f"   【 노란색 (함몰) 】 < -2mm: {stats['DEPRESSION']['count']:,} ({stats['DEPRESSION']['percentage']:.1f}%)")
    print(f"   【 빨강색 (요철) 】 > +2mm: {stats['PROTRUSION']['count']:,} ({stats['PROTRUSION']['percentage']:.1f}%)")
    
    return deviations_mm, stats

# ============================================================================
# 7-1. 다중신호 크랙 탐지 (곡률, 법선 불연속성, 이방성)
# ============================================================================

def compute_curvature(points, k=CURVATURE_K):
    """k-nearest neighbor를 이용한 곡률 계산"""
    xyz = points[:, :3]
    k = min(k, len(xyz) - 1)
    curvatures = np.zeros(len(xyz))

    if k < 1:
        return curvatures

    # k-d tree 구성
    tree = cKDTree(xyz)

    # 각 점의 k-nearest neighbor 찾기
    _, indices = tree.query(xyz, k=k+1)  # 자신 포함

    for i in range(len(xyz)):
        neighbor_indices = indices[i, 1:]  # 자신 제외
        neighbors = xyz[neighbor_indices]

        # 이웃 점들의 공분산 행렬
        cov = np.cov(neighbors.T)

        # 고유값 계산 (가장 작은 고유값이 곡률 지표)
        try:
            eigenvalues = np.linalg.eigvalsh(cov)
            eigenvalue_sum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2]
            # 이웃점들이 중복/동일 위치이면 분모가 0이 되어 NaN이 발생하므로 방어
            if eigenvalue_sum > 1e-6:
                # 가장 작은 고유값 (작을수록 곡률 높음)
                curvatures[i] = eigenvalues[0] / eigenvalue_sum
            else:
                curvatures[i] = 0
        except:
            curvatures[i] = 0

    return curvatures

def compute_normal_variance(points, k=CURVATURE_K):
    """법선벡터의 불연속성 계산 (각 점의 로컬 법선을 이웃점들의 로컬 법선과 비교)"""
    xyz = points[:, :3]
    n = len(xyz)
    k = min(k, n - 1)

    if k < 1:
        return np.zeros(n)

    # k-d tree 구성
    tree = cKDTree(xyz)
    _, indices = tree.query(xyz, k=k+1)

    # 각 점의 로컬 법선벡터를 한 번씩만 계산
    normals = np.zeros((n, 3))
    valid = np.zeros(n, dtype=bool)

    for i in range(n):
        neighbors = xyz[indices[i, 1:]]
        centered = neighbors - neighbors.mean(axis=0)
        try:
            _, _, Vt = np.linalg.svd(centered.T @ centered)
            normals[i] = Vt[-1, :]  # 가장 작은 특이값에 대응하는 벡터
            valid[i] = True
        except np.linalg.LinAlgError:
            pass

    # 각 점의 법선을 이웃점들의 법선과 비교 (부호 모호성에 영향받지 않도록 |cos| 사용)
    normal_variance = np.zeros(n)
    for i in range(n):
        if not valid[i]:
            continue
        neighbor_idx = indices[i, 1:]
        neighbor_idx = neighbor_idx[valid[neighbor_idx]]
        if len(neighbor_idx) == 0:
            continue
        cos_sim = np.abs(normals[neighbor_idx] @ normals[i])
        normal_variance[i] = np.mean(1.0 - cos_sim)

    return normal_variance

def compute_anisotropy(points, k=CURVATURE_K):
    """이방성 (점들의 분포 비이방성) 계산"""
    xyz = points[:, :3]
    k = min(k, len(xyz) - 1)
    anisotropy = np.zeros(len(xyz))

    if k < 1:
        return anisotropy

    tree = cKDTree(xyz)
    _, indices = tree.query(xyz, k=k+1)
    
    for i in range(len(xyz)):
        neighbor_indices = indices[i, 1:]
        neighbors = xyz[neighbor_indices]
        
        # 공분산 고유값 기반 이방성
        try:
            cov = np.cov(neighbors.T)
            eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
            
            # 이방성: (λ₁ - λ₂) / λ₁ (평면도) + (λ₂ - λ₃) / λ₁ (선형도)
            if eigenvalues[0] > 1e-6:
                planar = (eigenvalues[0] - eigenvalues[1]) / eigenvalues[0]
                linear = (eigenvalues[1] - eigenvalues[2]) / eigenvalues[0]
                anisotropy[i] = planar + linear
            else:
                anisotropy[i] = 0
        except:
            anisotropy[i] = 0
    
    return anisotropy

def detect_cracks(points, deviations, threshold_score=THRESHOLD_CRACK_SCORE):
    """다중신호 기반 크랙 탐지"""
    print("\n【 크랙 탐지 (다중신호) 】")
    
    # 세 가지 신호 계산
    print("   ① 곡률 계산 중...")
    curvature = compute_curvature(points, k=CURVATURE_K)
    
    print("   ② 법선 불연속성 계산 중...")
    normal_var = compute_normal_variance(points, k=CURVATURE_K)
    
    print("   ③ 이방성 계산 중...")
    anisotropy = compute_anisotropy(points, k=CURVATURE_K)
    
    # 정규화
    curvature_norm = (curvature - curvature.min()) / (curvature.max() - curvature.min() + 1e-6)
    normal_var_norm = (normal_var - normal_var.min()) / (normal_var.max() - normal_var.min() + 1e-6)
    anisotropy_norm = (anisotropy - anisotropy.min()) / (anisotropy.max() - anisotropy.min() + 1e-6)
    
    # 종합 크랙 점수 (0~1)
    crack_scores = (curvature_norm + normal_var_norm + anisotropy_norm) / 3
    
    # 편차도 고려 (절댓값 큰 것 = 하자 가능성)
    deviation_norm = np.abs(deviations) / (np.abs(deviations).max() + 1e-6)
    
    # 종합 점수에 편차 가중치 추가
    final_crack_scores = crack_scores * 0.5 + deviation_norm * 0.5
    
    # 크랙으로 판정
    crack_mask = final_crack_scores >= threshold_score
    
    # 크랙 통계
    crack_count = np.sum(crack_mask)
    crack_percentage = crack_count / len(points) * 100
    
    print(f"   크랙 탐지 결과:")
    print(f"      총 포인트: {len(points):,}개")
    print(f"      크랙 포인트: {crack_count:,}개 ({crack_percentage:.2f}%)")
    print(f"      크랙 점수 평균: {final_crack_scores.mean():.3f}")
    print(f"      크랙 점수 최댓값: {final_crack_scores.max():.3f}")
    
    return {
        'crack_scores': final_crack_scores,
        'crack_mask': crack_mask,
        'crack_count': crack_count,
        'crack_percentage': crack_percentage,
        'curvature': curvature,
        'normal_variance': normal_var,
        'anisotropy': anisotropy
    }

# ============================================================================
# 7-2. KCS 등급 판정
# ============================================================================

def grade_by_kcs(deviations, tolerance_mm=2.0):
    """KCS 41 46 01 기준으로 등급 판정"""
    print("\n【 KCS 41 46 01 기준 등급 판정 】")
    
    # 초과 포인트: 기준값을 벗어난 포인트
    exceed_mask = np.abs(deviations) > tolerance_mm
    exceed_count = np.sum(exceed_mask)
    exceed_ratio = exceed_count / len(deviations)
    
    # 등급 결정
    grade = 'D'  # 기본값
    for g, ratio_threshold in sorted(GRADE_BOUNDARIES.items()):
        if exceed_ratio <= ratio_threshold:
            grade = g
            break
    
    grade_names = {
        'A': '우수 (A)',
        'B': '양호 (B)',
        'C': '보통 (C)',
        'D': '불량 (D)'
    }
    
    print(f"   기준값: ±{tolerance_mm}mm")
    print(f"   초과 포인트: {exceed_count:,} / {len(deviations):,} ({exceed_ratio*100:.2f}%)")
    print(f"   등급: {grade_names[grade]}")
    
    print(f"\n   KCS 등급 기준:")
    print(f"      A (우수): 초과율 ≤ 5%")
    print(f"      B (양호): 초과율 ≤ 15%")
    print(f"      C (보통): 초과율 ≤ 30%")
    print(f"      D (불량): 초과율 > 30%")
    
    return {
        'grade': grade,
        'grade_name': grade_names[grade],
        'exceed_ratio': exceed_ratio,
        'exceed_count': exceed_count,
        'tolerance_mm': tolerance_mm
    }

# ============================================================================
# 8. 직선자 시뮬레이션 (편차 기반)
# ============================================================================

def straightedge_simulation_deviation_based(floor_points, deviations, x_range, y_range):
    """편차 기반 직선자 시뮬레이션"""
    print("\n【 Step 3: 직선자 시뮬레이션 (편차 기반) 】")
    
    x_min, x_max = x_range
    y_min, y_max = y_range
    
    x_length = (x_max - x_min) * 0.8
    y_length = (y_max - y_min) * 0.8
    x_start = x_min + ((x_max - x_min) - x_length) / 2
    y_start = y_min + ((y_max - y_min) - y_length) / 2

    print(f"   직선자 길이: X={x_length:.3f}m, Y={y_length:.3f}m")
    
    straightedge_results = []
    gap_threshold_mm = 3.0
    
    # X 방향
    for i in range(5):
        y_pos = y_min + (y_max - y_min) * (i + 1) / 6
        line_mask = (np.abs(floor_points[:, 1] - y_pos) < 0.01) & \
                    (floor_points[:, 0] >= x_start) & \
                    (floor_points[:, 0] <= x_start + x_length)
        
        if np.sum(line_mask) > 0:
            line_deviations = deviations[line_mask]
            gap = line_deviations.max() - line_deviations.min()
            passed = gap <= gap_threshold_mm
            
            straightedge_results.append({
                'direction': f'X-{i+1}',
                'gap_mm': float(gap),
                'passed': bool(passed),
                'point_count': int(np.sum(line_mask))
            })
            
            status = "✓" if passed else "✗"
            print(f"      X-{i+1}: {gap:.2f}mm {status}")
    
    # Y 방향
    for i in range(3):
        x_pos = x_min + (x_max - x_min) * (i + 1) / 4
        line_mask = (np.abs(floor_points[:, 0] - x_pos) < 0.01) & \
                    (floor_points[:, 1] >= y_start) & \
                    (floor_points[:, 1] <= y_start + y_length)
        
        if np.sum(line_mask) > 0:
            line_deviations = deviations[line_mask]
            gap = line_deviations.max() - line_deviations.min()
            passed = gap <= gap_threshold_mm
            
            straightedge_results.append({
                'direction': f'Y-{i+1}',
                'gap_mm': float(gap),
                'passed': bool(passed),
                'point_count': int(np.sum(line_mask))
            })
            
            status = "✓" if passed else "✗"
            print(f"      Y-{i+1}: {gap:.2f}mm {status}")
    
    if straightedge_results:
        all_passed = all([r['passed'] for r in straightedge_results])
        max_gap = max([r['gap_mm'] for r in straightedge_results])
    else:
        print("   ⚠️  직선자 테스트 라인에 포인트가 없어 검증을 수행하지 못했습니다.")
        all_passed = False
        max_gap = 0.0

    print(f"\n   3m 직선자 기준 (3.0mm 이하):")
    print(f"      테스트 직선: {len(straightedge_results)}개")
    print(f"      최대 틈새: {max_gap:.2f}mm")
    print(f"      판정: {'✅ 합격' if all_passed else '❌ 불합격'}")
    
    return {
        'results': straightedge_results,
        'all_passed': all_passed,
        'max_gap_mm': max_gap,
        'line_count': len(straightedge_results)
    }

# ============================================================================
# 9. 히트맵 생성
# ============================================================================

def generate_heatmap(floor_points, deviations, x_range, y_range, output_path):
    """고정 그리드 히트맵 생성 (보간 없이 실측 데이터만 표시)"""
    print("\n【 Step 4.1: 히트맵 생성 】")

    x_min, x_max = x_range
    y_min, y_max = y_range

    cell_size = 0.01
    nx = int(np.ceil((x_max - x_min) / cell_size))
    ny = int(np.ceil((y_max - y_min) / cell_size))

    print(f"   그리드: {nx} × {ny} = {nx * ny:,} 셀")

    ix = np.floor((floor_points[:, 0] - x_min) / cell_size).astype(np.int64)
    iy = np.floor((floor_points[:, 1] - y_min) / cell_size).astype(np.int64)
    valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)

    # 셀당 합계/개수를 누적한 뒤 마지막에 한 번만 나눠 진짜 평균을 계산
    sum_grid = np.zeros((ny, nx))
    count_grid = np.zeros((ny, nx))
    np.add.at(sum_grid, (iy[valid], ix[valid]), deviations[valid])
    np.add.at(count_grid, (iy[valid], ix[valid]), 1)

    with np.errstate(invalid='ignore', divide='ignore'):
        heatmap_data = sum_grid / count_grid
    heatmap_data[count_grid == 0] = np.nan

    # 보간 없이 실측 셀만 표시 (데이터 없는 셀은 투명 처리)
    masked_data = np.ma.masked_invalid(heatmap_data)

    # 렌더링
    fig, ax = plt.subplots(figsize=(12, 3))

    cmap_colors = ['yellow', 'green', 'red']
    cmap = mcolors.LinearSegmentedColormap.from_list('flatness', cmap_colors, N=256)
    cmap.set_bad(color='white', alpha=0)

    im = ax.imshow(masked_data, cmap=cmap, origin='lower', vmin=-10, vmax=10,
                    extent=[x_min, x_max, y_min, y_max], aspect='auto', interpolation='nearest')
    
    plt.colorbar(im, ax=ax, label='Deviation (mm)', shrink=0.8)
    
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title('Floor Flatness Heatmap (3-Color Classification)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_path}/heatmap.png', dpi=150, bbox_inches='tight')
    print(f"   ✅ 저장: {output_path}/heatmap.png")
    plt.close()
    
    return heatmap_data

# ============================================================================
# 10. 3D 시각화
# ============================================================================

def generate_3d_visualization(floor_points, deviations, output_path):
    """3D 포인트 클라우드 시각화 (ASCII PCD 파일로 직접 저장, open3d 불필요)"""
    print("\n【 Step 4.2: 3D 시각화 】")

    flat_mask = (np.abs(deviations) <= 2.0)
    depression_mask = (deviations < -2.0)
    protrusion_mask = (deviations > 2.0)

    colors = np.zeros((len(floor_points), 3), dtype=np.int32)
    colors[flat_mask] = [0, 200, 100]
    colors[depression_mask] = [255, 200, 0]
    colors[protrusion_mask] = [255, 100, 100]

    print(f"   포인트: {len(floor_points):,}")
    print(f"      초록색(평탄): {np.sum(flat_mask):,}")
    print(f"      노란색(함몰): {np.sum(depression_mask):,}")
    print(f"      빨강색(요철): {np.sum(protrusion_mask):,}")

    n = len(floor_points)
    pcd_path = f'{output_path}/pointcloud_3d.pcd'
    header = (
        "# .PCD v0.7 - Point Cloud Data file format\n"
        "VERSION 0.7\n"
        "FIELDS x y z r g b\n"
        "SIZE 4 4 4 1 1 1\n"
        "TYPE F F F U U U\n"
        "COUNT 1 1 1 1 1 1\n"
        f"WIDTH {n}\n"
        "HEIGHT 1\n"
        "VIEWPOINT 0 0 0 1 0 0 0\n"
        f"POINTS {n}\n"
        "DATA ascii"
    )
    combined = np.hstack([floor_points.astype(np.float64), colors])
    np.savetxt(pcd_path, combined, fmt='%.6f %.6f %.6f %d %d %d', header=header, comments='')
    print(f"   ✅ 저장: {pcd_path}")

# ============================================================================
# 11. 그래프 생성
# ============================================================================

def _plot_histogram(ax, dev_filtered):
    ax.hist(dev_filtered, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2)
    ax.axvline(x=-2, color='orange', linestyle='--', linewidth=1)
    ax.axvline(x=2, color='orange', linestyle='--', linewidth=1)
    ax.set_xlabel('Deviation (mm)')
    ax.set_ylabel('Frequency')
    ax.set_title('Deviation Distribution')
    ax.grid(True, alpha=0.3)

def _plot_pie(ax, stats):
    labels = ['FLAT', 'DEPRESSION', 'PROTRUSION']
    sizes = [stats['FLAT']['percentage'], stats['DEPRESSION']['percentage'], stats['PROTRUSION']['percentage']]
    colors_pie = ['green', 'yellow', 'red']
    ax.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%', startangle=90)
    ax.set_title('Category Distribution')

def _plot_cumulative(ax, dev_filtered):
    sorted_dev = np.sort(dev_filtered)
    cumsum = np.arange(1, len(sorted_dev) + 1) / len(sorted_dev)
    ax.plot(sorted_dev, cumsum, 'b-', linewidth=2)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1)
    ax.axvline(x=-2, color='orange', linestyle='--', linewidth=1)
    ax.axvline(x=2, color='orange', linestyle='--', linewidth=1)
    ax.set_xlabel('Deviation (mm)')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('Cumulative Distribution')
    ax.grid(True, alpha=0.3)

def _plot_boxplot(ax, dev_filtered):
    ax.boxplot([dev_filtered], tick_labels=['Floor'])
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1)
    ax.axhline(y=-2, color='orange', linestyle='--', linewidth=1)
    ax.axhline(y=2, color='orange', linestyle='--', linewidth=1)
    ax.set_ylabel('Deviation (mm)')
    ax.set_title('Box Plot')
    ax.grid(True, alpha=0.3)

def generate_graphs(deviations, stats, output_path):
    """통계 그래프 생성 (통합 이미지 + 대시보드용 개별 이미지)"""
    print("\n【 Step 4.3: 그래프 생성 】")

    dev_filtered = deviations

    # 통합 이미지 (리포트/다운로드용)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    _plot_histogram(axes[0, 0], dev_filtered)
    _plot_pie(axes[0, 1], stats)
    _plot_cumulative(axes[1, 0], dev_filtered)
    _plot_boxplot(axes[1, 1], dev_filtered)
    plt.tight_layout()
    plt.savefig(f'{output_path}/additional_graphs.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✅ 저장: {output_path}/additional_graphs.png")

    # 개별 이미지 (대시보드에서 그래프별 설명과 함께 표시하기 위함)
    individual_plots = [
        ('graph_histogram.png', _plot_histogram, (dev_filtered,)),
        ('graph_pie.png', _plot_pie, (stats,)),
        ('graph_cumulative.png', _plot_cumulative, (dev_filtered,)),
        ('graph_boxplot.png', _plot_boxplot, (dev_filtered,)),
    ]
    for filename, plot_fn, args in individual_plots:
        fig_i, ax_i = plt.subplots(figsize=(7, 5))
        plot_fn(ax_i, *args)
        plt.tight_layout()
        plt.savefig(f'{output_path}/{filename}', dpi=150, bbox_inches='tight')
        plt.close(fig_i)
    print(f"   ✅ 개별 그래프 4개 저장 완료")

# ============================================================================
# 12. 리포트 생성
# ============================================================================

def generate_reports(metadata, filter_stats, tilt_angle, deviations, stats, straightedge_results, output_path, crack_results=None, kcs_results=None):
    """JSON과 텍스트 리포트 생성 (크랙 탐지 및 KCS 등급 포함)"""
    print("\n【 Step 5: 리포트 생성 】")

    dev_filtered = deviations
    
    def to_native(obj):
        if isinstance(obj, (np.floating, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        return obj
    
    # JSON
    report_json = {
        'metadata': {
            'filename': metadata['filename'],
            'scan_time': metadata['scan_time'],
            'x_range_m': list(metadata['x_range']),
            'y_range_m': list(metadata['y_range']),
            'z_range_m': list(metadata['z_range']),
        },
        'coordinate_correction': {
            'tilt_angle_degrees': to_native(tilt_angle),
            'correction_applied': tilt_angle >= 0.5,
        },
        'filtering': {
            'original_points': filter_stats['original'],
            'floor_points': filter_stats['after_strong_filter'],
            'retention_rate_percent': to_native(filter_stats['retention_rate']),
        },
        'statistics': {
            'floor_points': int(len(dev_filtered)),
            'max_deviation_mm': to_native(dev_filtered.max()),
            'min_deviation_mm': to_native(dev_filtered.min()),
            'mean_deviation_mm': to_native(dev_filtered.mean()),
            'std_dev_mm': to_native(dev_filtered.std()),
            'rms_mm': to_native(np.sqrt(np.mean(dev_filtered**2))),
            'pv_mm': to_native(dev_filtered.max() - dev_filtered.min()),
        },
        'categories': {
            'FLAT': {
                'points': stats['FLAT']['count'],
                'percentage': to_native(stats['FLAT']['percentage']),
                'mean_deviation_mm': to_native(stats['FLAT']['mean_dev']),
            },
            'DEPRESSION': {
                'points': stats['DEPRESSION']['count'],
                'percentage': to_native(stats['DEPRESSION']['percentage']),
                'mean_deviation_mm': to_native(stats['DEPRESSION']['mean_dev']),
            },
            'PROTRUSION': {
                'points': stats['PROTRUSION']['count'],
                'percentage': to_native(stats['PROTRUSION']['percentage']),
                'mean_deviation_mm': to_native(stats['PROTRUSION']['mean_dev']),
            }
        },
        'straightedge_validation': {
            'passed': to_native(straightedge_results['all_passed']),
            'max_gap_mm': to_native(straightedge_results['max_gap_mm']),
            'line_count': straightedge_results['line_count'],
        }
    }
    
    # 크랙 탐지 결과 추가
    if crack_results is not None:
        report_json['crack_detection'] = {
            'crack_count': to_native(crack_results['crack_count']),
            'crack_percentage': to_native(crack_results['crack_percentage']),
            'crack_score_mean': to_native(np.mean(crack_results['crack_scores'])),
            'crack_score_max': to_native(np.max(crack_results['crack_scores'])),
        }
    
    # KCS 등급 결과 추가
    if kcs_results is not None:
        report_json['kcs_grading'] = {
            'grade': kcs_results['grade'],
            'grade_name': kcs_results['grade_name'],
            'exceed_ratio_percent': to_native(kcs_results['exceed_ratio'] * 100),
            'exceed_count': to_native(kcs_results['exceed_count']),
            'tolerance_mm': to_native(kcs_results['tolerance_mm']),
        }
    
    with open(f'{output_path}/report.json', 'w', encoding='utf-8') as f:
        json.dump(report_json, f, indent=2, ensure_ascii=False)
    print(f"   ✅ JSON 저장: {output_path}/report.json")
    
    # TEXT
    report_text = f"""
================================================================================
🏢 건물 바닥 슬래브 평활도 측정 리포트 (v3.2 - KCS 기준 + 크랙 탐지)
Floor Slab Flatness Measurement Report
================================================================================

【 1. 스캔 정보 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
파일명:           {metadata['filename']}
스캔 일시:        {metadata['scan_time']}
측정 도구:        SiteScape (스마트폰 LiDAR)

스캔 범위:
  X: {metadata['x_range'][0]:.3f}m ~ {metadata['x_range'][1]:.3f}m
  Y: {metadata['y_range'][0]:.3f}m ~ {metadata['y_range'][1]:.3f}m
  Z: {metadata['z_range'][0]:.3f}m ~ {metadata['z_range'][1]:.3f}m

【 2. 좌표계 보정 정보 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
감지된 기울기 각도: {tilt_angle:.2f}°
보정 적용:         {'✅ 보정됨' if tilt_angle >= 0.5 else '✅ 이미 수평 (보정 불필요)'}

해석:
  - 0° ~ 0.5°: 이미 수평 상태
  - 0.5° ~ 90°: 자동 회전 보정 적용
  - 90° 이상: 좌표계 크게 기울어짐 (보정 적용)

【 3. 필터링 정보 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
원본 포인트:      {filter_stats['original']:,}개
필터 후 바닥:     {filter_stats['after_strong_filter']:,}개
유지율:           {filter_stats['retention_rate']:.1f}%
필터링 방법:      Z < {filter_stats['z_threshold']}m (명확한 바닥만 추출)

【 4. 편차 통계 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
바닥 포인트:      {len(dev_filtered):,}개
최소 편차:        {dev_filtered.min():.2f}mm
최대 편차:        {dev_filtered.max():.2f}mm
평균 편차:        {dev_filtered.mean():.2f}mm
표준편차:         {dev_filtered.std():.2f}mm
RMS:              {np.sqrt(np.mean(dev_filtered**2)):.2f}mm
PV (P-V):         {dev_filtered.max() - dev_filtered.min():.2f}mm

【 5. 카테고리별 분석 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【 초록색 (평탄) 】 -2mm ~ +2mm
  포인트: {stats['FLAT']['count']:,} ({stats['FLAT']['percentage']:.1f}%)
  평균 편차: {stats['FLAT']['mean_dev']:.2f}mm
  판정: ✅ 합격

【 노란색 (함몰) 】 < -2mm
  포인트: {stats['DEPRESSION']['count']:,} ({stats['DEPRESSION']['percentage']:.1f}%)
  평균 편차: {stats['DEPRESSION']['mean_dev']:.2f}mm
  판정: ⚠️  주의

【 빨강색 (요철) 】 > +2mm
  포인트: {stats['PROTRUSION']['count']:,} ({stats['PROTRUSION']['percentage']:.1f}%)
  평균 편차: {stats['PROTRUSION']['mean_dev']:.2f}mm
  판정: ❌ 부적합

【 6. 크랙 탐지 (다중신호 기반) 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
{f"크랙 포인트:      {crack_results['crack_count']:,} / {len(dev_filtered):,} ({crack_results['crack_percentage']:.2f}%)" if crack_results else "크랙 탐지: 미수행"}
{f"크랙 점수 평균:   {np.mean(crack_results['crack_scores']):.3f}" if crack_results else ""}
{f"크랙 점수 최댓값: {np.max(crack_results['crack_scores']):.3f}" if crack_results else ""}
탐지 신호: 곡률 + 법선불연속성 + 이방성

【 7. KCS 41 46 01 등급 판정 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
{f"등급:             {kcs_results['grade_name']}" if kcs_results else "등급 판정: 미수행"}
{f"기준값:           ±{kcs_results['tolerance_mm']}mm" if kcs_results else ""}
{f"초과 포인트:      {kcs_results['exceed_count']:,} / {len(dev_filtered):,} ({kcs_results['exceed_ratio']*100:.2f}%)" if kcs_results else ""}

기준:
  - A 우수: 초과율 ≤ 5%
  - B 양호: 초과율 ≤ 15%
  - C 보통: 초과율 ≤ 30%
  - D 불량: 초과율 > 30%

【 8. 3m 직선자 기준 검증 (참고) 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
기준: 3.0mm 이하
테스트 직선: {straightedge_results['line_count']}개
최대 틈새: {straightedge_results['max_gap_mm']:.2f}mm
판정: {'✅ 합격' if straightedge_results['all_passed'] else '❌ 불합격'}

【 9. 최종 판정 】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
{f"KCS 등급:         {kcs_results['grade_name']}" if kcs_results else "최종 등급: 미판정"}
3m 직선자 기준:   {'✅ 만족' if straightedge_results['all_passed'] else '❌ 미충족'} (참고용)

================================================================================
생성 일시: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Version: 3.2 (자동 좌표계 보정 + KCS 등급 판정 + 다중신호 크랙 탐지)
================================================================================
"""
    
    with open(f'{output_path}/report.txt', 'w', encoding='utf-8') as f:
        f.write(report_text)
    print(f"   ✅ 텍스트 저장: {output_path}/report.txt")

# ============================================================================
# 13. 메인 함수
# ============================================================================

def main(ply_filepath, output_dir=str(Path.home() / 'Desktop' / 'flatness_output'), z_threshold=0.5, roi_bounds=None):
    """메인 실행 함수

    roi_bounds: (x_min, x_max, y_min, y_max) 지정 시 해당 구역(예: 테스트보드)만 분석.
                preview_scan()으로 좌표를 먼저 확인한 뒤 사용.
    """
    print("\n" + "=" * 80)
    print("🏢 건물 바닥 슬래브 평활도 측정 프로그램 v3.2 (KCS 기준 + 크랙 탐지) 시작")
    print("=" * 80)
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    try:
        # ========== Phase 1: 바닥 추출, 좌표계 보정, 평면 검출 ==========
        print("\n" + "=" * 80)
        print("【 Phase 1: 바닥 추출, 좌표계 보정, 평면 검출 】")
        print("=" * 80)
        
        # Step 1.1: 데이터 로드
        points, metadata = read_ply_file(ply_filepath)
        
        # Step 1.2: 강한 필터링
        floor_candidates, filter_stats = strong_floor_filtering(points, z_threshold=z_threshold)

        # Step 1.25: ROI 크롭 (지정 시에만 - 예: 테스트보드 영역만 분석)
        if roi_bounds is not None:
            floor_candidates = crop_to_roi(floor_candidates, roi_bounds)

        # Step 1.3-1.5: 자동 좌표계 보정
        normal_vector, tilt_angle, needs_correction = detect_tilt_and_get_normal(floor_candidates)
        
        if needs_correction:
            rotation_matrix = calculate_rotation_matrix(normal_vector)
            floor_candidates_corrected = apply_coordinate_correction(floor_candidates, rotation_matrix)
        else:
            floor_candidates_corrected = floor_candidates
            rotation_matrix = np.eye(3)
        
        # Step 1.6: RANSAC 평면 검출
        points_xyz = floor_candidates_corrected[:, :3]
        plane_ransac, inlier_indices = ransac_plane_detection(points_xyz, threshold=0.05)
        
        # ========== Phase 2: 평활도 계산 ==========
        print("\n" + "=" * 80)
        print("【 Phase 2: 평활도 계산 】")
        print("=" * 80)
        
        # Step 2.1: LSQ 평면 피팅
        plane_final = least_squares_plane_fitting(points_xyz, plane_ransac, inlier_indices)
        
        # Step 2.2-2.4: 편차 계산 및 분류
        deviations, stats = calculate_deviations_and_classify(points_xyz, plane_final)
        
        # ========== Phase 3: 직선자 검증 ==========
        print("\n" + "=" * 80)
        print("【 Phase 3: 직선자 검증 】")
        print("=" * 80)
        
        straightedge_results = straightedge_simulation_deviation_based(
            points_xyz, deviations,
            (points_xyz[:, 0].min(), points_xyz[:, 0].max()),
            (points_xyz[:, 1].min(), points_xyz[:, 1].max())
        )
        
        # ========== Phase 3.5: 크랙 탐지 ==========
        print("\n" + "=" * 80)
        print("【 Phase 3.5: 다중신호 크랙 탐지 】")
        print("=" * 80)
        
        crack_results = detect_cracks(points_xyz, deviations, threshold_score=THRESHOLD_CRACK_SCORE)
        
        # ========== Phase 3.7: KCS 등급 판정 ==========
        print("\n" + "=" * 80)
        print("【 Phase 3.7: KCS 41 46 01 등급 판정 】")
        print("=" * 80)
        
        kcs_results = grade_by_kcs(deviations, tolerance_mm=2.0)
        
        # ========== Phase 4: 시각화 ==========
        print("\n" + "=" * 80)
        print("【 Phase 4: 시각화 】")
        print("=" * 80)
        
        generate_heatmap(points_xyz, deviations, 
                        (points_xyz[:, 0].min(), points_xyz[:, 0].max()),
                        (points_xyz[:, 1].min(), points_xyz[:, 1].max()), output_dir)
        generate_3d_visualization(points_xyz, deviations, output_dir)
        generate_graphs(deviations, stats, output_dir)
        
        # ========== Phase 5: 리포트 ==========
        print("\n" + "=" * 80)
        print("【 Phase 5: 리포트 생성 】")
        print("=" * 80)
        
        generate_reports(metadata, filter_stats, tilt_angle, deviations, stats, straightedge_results, output_dir, crack_results, kcs_results)
        
        # ========== 완료 ==========
        print("\n" + "=" * 80)
        print("✅ 분석 완료!")
        print("=" * 80)
        print(f"\n📁 출력: {output_dir}/")
        print(f"   ✓ report.txt")
        print(f"   ✓ report.json")
        print(f"   ✓ heatmap.png")
        print(f"   ✓ additional_graphs.png")
        print(f"   ✓ pointcloud_3d.pcd")
        
        print(f"\n【 최종 판정 】")
        print(f"KCS 등급: {kcs_results['grade_name']}")
        print(f"좌표계 기울기: {tilt_angle:.2f}° {'(자동 보정됨)' if tilt_angle >= 0.5 else '(정상)'}")
        print(f"3m 직선자 기준: {'✅ 합격' if straightedge_results['all_passed'] else '❌ 불합격'} (참고용)")
        print(f"최대 틈새: {straightedge_results['max_gap_mm']:.2f}mm")
        print(f"평탄 비율: {stats['FLAT']['percentage']:.1f}%")
        print(f"함몰 비율: {stats['DEPRESSION']['percentage']:.1f}%")
        print(f"요철 비율: {stats['PROTRUSION']['percentage']:.1f}%")
        print(f"\n【 크랙 탐지 결과 】")
        print(f"크랙 포인트: {crack_results['crack_count']:,} / {len(points_xyz):,} ({crack_results['crack_percentage']:.2f}%)")
        print(f"크랙 점수 평균: {np.mean(crack_results['crack_scores']):.3f}")

        return {
            'metadata': metadata,
            'filter_stats': filter_stats,
            'tilt_angle': tilt_angle,
            'stats': stats,
            'straightedge_results': straightedge_results,
            'crack_results': crack_results,
            'kcs_results': kcs_results,
            'output_dir': output_dir,
            'points_xyz': points_xyz,
            'deviations': deviations,
        }

    except Exception as e:
        print(f"\n❌ 오류: {e}")
        import traceback
        traceback.print_exc()
        return {'error': str(e)}

if __name__ == '__main__':
    ply_file = str(Path.home() / 'Desktop' / 'ply' / 'S8동 복도_바닥_수평.ply')
    output_dir = str(Path.home() / 'Desktop' / 'flatness_output')

    if os.path.exists(ply_file):
        main(ply_file, output_dir)
    else:
        print(f"❌ 파일 없음: {ply_file}")
