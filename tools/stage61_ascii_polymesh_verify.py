#!/usr/bin/env python3
"""独立读取并验证 Stage 6.1 ASCII OpenFOAM polyMesh，不链接 cartmesh。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


REQUIRED = ("points", "faces", "owner", "neighbour", "boundary")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def list_body(text: str) -> tuple[int, str]:
    match = re.search(r"\n\s*(\d+)\s*\n\s*\(\s*\n", text)
    if match is None:
        raise ValueError("找不到 OpenFOAM 列表头")
    end = text.rfind("\n)")
    if end < match.end():
        raise ValueError("OpenFOAM 列表没有闭合")
    return int(match.group(1)), text[match.end():end]


def parse_points(path: Path) -> np.ndarray:
    expected, body = list_body(path.read_text(encoding="utf-8"))
    rows = re.findall(
        r"\(\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s*\)",
        body,
    )
    points = np.asarray(rows, dtype=np.float64)
    if len(points) != expected or points.shape != (expected, 3):
        raise ValueError("points 数量与列表头不一致")
    if not np.all(np.isfinite(points)):
        raise ValueError("points 包含非有限坐标")
    return points


def parse_faces(path: Path) -> list[np.ndarray]:
    expected, body = list_body(path.read_text(encoding="utf-8"))
    faces: list[np.ndarray] = []
    for line in body.splitlines():
        line = line.strip()
        if not line:
            continue
        match = re.fullmatch(r"(\d+)\s*\(([^)]*)\)", line)
        if match is None:
            raise ValueError(f"无法解析 face：{line[:80]}")
        count = int(match.group(1))
        ids = np.fromstring(match.group(2), sep=" ", dtype=np.int64)
        if len(ids) != count or count < 3:
            raise ValueError("face 顶点数无效")
        faces.append(ids)
    if len(faces) != expected:
        raise ValueError("faces 数量与列表头不一致")
    return faces


def parse_labels(path: Path) -> np.ndarray:
    expected, body = list_body(path.read_text(encoding="utf-8"))
    values = np.fromstring(body, sep=" ", dtype=np.int64)
    if len(values) != expected:
        raise ValueError(f"{path.name} 数量与列表头不一致")
    return values


def parse_boundary(path: Path) -> list[dict[str, int | str]]:
    expected, body = list_body(path.read_text(encoding="utf-8"))
    patches = []
    pattern = re.compile(r"(\w+)\s*\{([^}]*)\}", re.DOTALL)
    for name, block in pattern.findall(body):
        type_match = re.search(r"\btype\s+(\w+)\s*;", block)
        count_match = re.search(r"\bnFaces\s+(\d+)\s*;", block)
        start_match = re.search(r"\bstartFace\s+(\d+)\s*;", block)
        if not (type_match and count_match and start_match):
            raise ValueError(f"boundary patch {name} 字段不完整")
        patches.append({
            "name": name,
            "type": type_match.group(1),
            "nFaces": int(count_match.group(1)),
            "startFace": int(start_match.group(1)),
        })
    if len(patches) != expected:
        raise ValueError("boundary patch 数量与列表头不一致")
    return patches


def polygon_geometry(vertices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    origin = vertices[0]
    area_vector = np.zeros(3)
    weighted_centroid = np.zeros(3)
    area_sum = 0.0
    for index in range(1, len(vertices) - 1):
        cross = np.cross(vertices[index] - origin, vertices[index + 1] - origin)
        triangle_area = 0.5 * float(np.linalg.norm(cross))
        area_vector += 0.5 * cross
        weighted_centroid += triangle_area * (
            origin + vertices[index] + vertices[index + 1]) / 3.0
        area_sum += triangle_area
    if not area_sum > 0.0:
        raise ValueError("发现零面积 face")
    return area_vector, weighted_centroid / area_sum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    poly_mesh = args.case.resolve() / "constant/polyMesh"
    missing = [name for name in REQUIRED if not (poly_mesh / name).is_file()]
    if missing:
        raise SystemExit(f"缺少 polyMesh 文件：{missing}")

    failures: list[str] = []
    points = parse_points(poly_mesh / "points")
    faces = parse_faces(poly_mesh / "faces")
    owner = parse_labels(poly_mesh / "owner")
    neighbour = parse_labels(poly_mesh / "neighbour")
    patches = parse_boundary(poly_mesh / "boundary")

    if len(owner) != len(faces):
        failures.append("owner_face_count_mismatch")
    if len(neighbour) > len(owner):
        failures.append("neighbour_count_exceeds_owner")
    if any(np.any(face < 0) or np.any(face >= len(points)) for face in faces):
        failures.append("point_index_out_of_range")
    if np.any(owner < 0) or np.any(neighbour < 0):
        failures.append("negative_cell_index")
    if len(neighbour) and np.any(owner[:len(neighbour)] >= neighbour):
        failures.append("owner_not_less_than_neighbour")
    if len(neighbour) and np.any(owner[:len(neighbour)][1:] < owner[:len(neighbour)][:-1]):
        failures.append("internal_owner_not_sorted")

    cell_count = int(max(owner.max(initial=-1), neighbour.max(initial=-1)) + 1)
    expected_start = len(neighbour)
    for patch in patches:
        if patch["startFace"] != expected_start:
            failures.append("boundary_patch_not_contiguous")
            break
        expected_start += int(patch["nFaces"])
    if expected_start != len(faces):
        failures.append("boundary_patch_does_not_cover_faces")

    cell_faces: list[list[tuple[int, int]]] = [[] for _ in range(cell_count)]
    for face_id, cell in enumerate(owner):
        cell_faces[int(cell)].append((face_id, 1))
    for face_id, cell in enumerate(neighbour):
        cell_faces[int(cell)].append((face_id, -1))

    nonmanifold_cell_edges = 0
    nonclosed_cells = 0
    nonpositive_volume_cells = 0
    minimum_volume = math.inf
    maximum_closure = 0.0
    for entries in cell_faces:
        edge_counts: Counter[tuple[int, int]] = Counter()
        closure = np.zeros(3)
        volume = 0.0
        for face_id, sign in entries:
            face = faces[face_id]
            vertices = points[face]
            area_vector, centroid = polygon_geometry(vertices)
            closure += sign * area_vector
            volume += sign * float(np.dot(area_vector, centroid)) / 3.0
            for index, first in enumerate(face):
                second = face[(index + 1) % len(face)]
                edge_counts[(min(int(first), int(second)),
                             max(int(first), int(second)))] += 1
        nonmanifold_cell_edges += sum(count != 2 for count in edge_counts.values())
        closure_norm = float(np.linalg.norm(closure))
        maximum_closure = max(maximum_closure, closure_norm)
        scale = sum(float(np.linalg.norm(polygon_geometry(points[faces[f]])[0]))
                    for f, _ in entries)
        if closure_norm > max(1.0e-12, 1.0e-10 * scale):
            nonclosed_cells += 1
        minimum_volume = min(minimum_volume, volume)
        if not volume > 0.0:
            nonpositive_volume_cells += 1

    if nonmanifold_cell_edges:
        failures.append("nonmanifold_cell_edge_incidence")
    if nonclosed_cells:
        failures.append("nonclosed_cells")
    if nonpositive_volume_cells:
        failures.append("nonpositive_cell_volume")

    result = {
        "schema": "cartmesh-stage61-ascii-polymesh-reader-v1",
        "status": "pass" if not failures else "fail",
        "case": str(args.case.resolve()),
        "pointCount": len(points),
        "faceCount": len(faces),
        "internalFaceCount": len(neighbour),
        "boundaryFaceCount": len(faces) - len(neighbour),
        "cellCount": cell_count,
        "boundaryPatches": patches,
        "nonmanifoldCellEdgeCount": nonmanifold_cell_edges,
        "nonclosedCellCount": nonclosed_cells,
        "nonpositiveVolumeCellCount": nonpositive_volume_cells,
        "minimumCellVolume": minimum_volume,
        "maximumAreaClosureResidual": maximum_closure,
        "failures": failures,
        "polyMeshSha256": {name: sha256(poly_mesh / name) for name in REQUIRED},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
