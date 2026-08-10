#!/usr/bin/env python3
"""独立读取并验证 OpenFOAM binary polyMesh，不链接 cartmesh。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import mmap
import os
import platform
import re
import resource
import struct
import sys
import time
from pathlib import Path


LIST_START = re.compile(rb"\s*(\d+)\s*\(")


class MappedFile:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.file = path.open("rb")
        self.data = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)

    def close(self) -> None:
        self.data.close()
        self.file.close()

    def first_list(self) -> tuple[int, int]:
        header_end = self.data.find(b"}\n")
        if header_end < 0:
            raise ValueError(f"{self.path}: missing FoamFile header")
        match = LIST_START.match(self.data, header_end + 2)
        if match is None:
            raise ValueError(f"{self.path}: missing first list")
        return int(match.group(1)), match.end()


def view_i32(data: mmap.mmap, offset: int, count: int) -> memoryview:
    return memoryview(data)[offset : offset + count * 4].cast("i")


def view_f64(data: mmap.mmap, offset: int, count: int) -> memoryview:
    return memoryview(data)[offset : offset + count * 8].cast("d")


def parse_label_list(path: Path) -> tuple[MappedFile, memoryview]:
    mapped = MappedFile(path)
    count, offset = mapped.first_list()
    end = offset + count * 4
    if end > len(mapped.data):
        mapped.close()
        raise ValueError(f"{path}: truncated label list")
    return mapped, view_i32(mapped.data, offset, count)


def parse_points(path: Path) -> tuple[MappedFile, int, memoryview]:
    mapped = MappedFile(path)
    count, offset = mapped.first_list()
    values = view_f64(mapped.data, offset, count * 3)
    return mapped, count, values


def parse_faces(
    path: Path,
) -> tuple[MappedFile, memoryview, memoryview]:
    mapped = MappedFile(path)
    offset_count, offset_start = mapped.first_list()
    offsets = view_i32(mapped.data, offset_start, offset_count)
    offset_end = offset_start + offset_count * 4
    match = LIST_START.match(mapped.data, offset_end + 1)
    if match is None:
        mapped.close()
        raise ValueError(f"{path}: missing face vertex list")
    vertex_count = int(match.group(1))
    vertices = view_i32(mapped.data, match.end(), vertex_count)
    return mapped, offsets, vertices


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_boundary(path: Path) -> list[dict[str, object]]:
    text = path.read_text(encoding="ascii")
    body = text[text.index("}") + 1 :]
    count_match = re.search(r"\s*(\d+)\s*\(", body)
    if count_match is None:
        raise ValueError("boundary: missing patch count")
    patch_count = int(count_match.group(1))
    pattern = re.compile(
        r"([A-Za-z_][A-Za-z0-9_]*)\s*\{\s*"
        r"type\s+([^;]+);\s*nFaces\s+(\d+);\s*"
        r"startFace\s+(\d+);\s*\}",
        re.MULTILINE,
    )
    patches = [
        {
            "name": match.group(1),
            "type": match.group(2).strip(),
            "nFaces": int(match.group(3)),
            "startFace": int(match.group(4)),
        }
        for match in pattern.finditer(body)
    ]
    if len(patches) != patch_count:
        raise ValueError(
            f"boundary: declared {patch_count} patches, parsed {len(patches)}"
        )
    return patches


def diagnose_cell_edges(
    points: memoryview,
    offsets: memoryview,
    vertices: memoryview,
    owner: memoryview,
    neighbour: memoryview,
    cell_count: int,
) -> dict[str, object]:
    def face_area(face_id: int) -> float:
        begin = offsets[face_id]
        end = offsets[face_id + 1]
        origin_id = vertices[begin]
        origin = (
            points[3 * origin_id],
            points[3 * origin_id + 1],
            points[3 * origin_id + 2],
        )
        area_vector = [0.0, 0.0, 0.0]
        for position in range(begin + 1, end - 1):
            first_id = vertices[position]
            second_id = vertices[position + 1]
            first = (
                points[3 * first_id] - origin[0],
                points[3 * first_id + 1] - origin[1],
                points[3 * first_id + 2] - origin[2],
            )
            second = (
                points[3 * second_id] - origin[0],
                points[3 * second_id + 1] - origin[1],
                points[3 * second_id + 2] - origin[2],
            )
            area_vector[0] += first[1] * second[2] - first[2] * second[1]
            area_vector[1] += first[2] * second[0] - first[0] * second[2]
            area_vector[2] += first[0] * second[1] - first[1] * second[0]
        return 0.5 * math.sqrt(sum(value * value for value in area_vector))

    cell_faces: list[list[int]] = [[] for _ in range(cell_count)]
    for face_id, cell in enumerate(owner):
        cell_faces[cell].append(face_id)
    for face_id, cell in enumerate(neighbour):
        cell_faces[cell].append(face_id)
    worst_cell = -1
    worst_edge: tuple[int, int] | None = None
    worst_incidence = 0
    failure_count = 0
    first_failure: dict[str, object] | None = None
    for cell, faces in enumerate(cell_faces):
        edges: dict[tuple[int, int], list[int]] = {}
        for face_id in faces:
            begin = offsets[face_id]
            end = offsets[face_id + 1]
            for position in range(begin, end):
                first = vertices[position]
                second = vertices[begin if position + 1 == end else position + 1]
                key = (first, second) if first < second else (second, first)
                edges.setdefault(key, []).append(face_id)
        for edge, incident in edges.items():
            incidence = len(incident)
            if incidence > worst_incidence:
                worst_cell = cell
                worst_edge = edge
                worst_incidence = incidence
            if incidence != 2:
                failure_count += 1
                if first_failure is None:
                    first_failure = {
                        "cell": cell,
                        "edge": list(edge),
                        "incidence": incidence,
                        "faces": incident,
                        "cellFaces": faces,
                    }
    if first_failure is not None:
        details = []
        for face_id in first_failure["cellFaces"]:
            begin = offsets[face_id]
            end = offsets[face_id + 1]
            details.append(
                {
                    "face": face_id,
                    "owner": owner[face_id],
                    "neighbour": (
                        neighbour[face_id] if face_id < len(neighbour) else None
                    ),
                    "area": face_area(face_id),
                    "pointIds": list(vertices[begin:end]),
                }
            )
        first_failure["cellFaceDetails"] = details
    return {
        "checked": True,
        "nonTwoIncidentCellEdgeCount": failure_count,
        "firstFailure": first_failure,
        "maximumIncidence": worst_incidence,
        "maximumIncidenceCell": worst_cell,
        "maximumIncidenceEdge": list(worst_edge) if worst_edge else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnose-cell-edges", action="store_true")
    args = parser.parse_args()
    start = time.perf_counter()
    poly_mesh = args.case / "constant" / "polyMesh"
    paths = {
        name: poly_mesh / name
        for name in ("points", "faces", "owner", "neighbour", "boundary")
    }
    mapped_files: list[MappedFile] = []
    try:
        points_map, point_count, points = parse_points(paths["points"])
        faces_map, offsets, vertices = parse_faces(paths["faces"])
        owner_map, owner = parse_label_list(paths["owner"])
        neighbour_map, neighbour = parse_label_list(paths["neighbour"])
        mapped_files.extend((points_map, faces_map, owner_map, neighbour_map))
        face_count = len(offsets) - 1
        failures: list[str] = []
        if not all(math.isfinite(value) for value in points):
            failures.append("nonfinite_point_coordinate")
        if face_count < 0 or len(owner) != face_count:
            failures.append("owner_face_count_mismatch")
        if not offsets or offsets[0] != 0 or offsets[-1] != len(vertices):
            failures.append("face_offset_extent_mismatch")
        previous = -1
        for index, offset in enumerate(offsets):
            if offset < previous:
                failures.append("nonmonotone_face_offsets")
                break
            if index != 0 and offset - previous < 3:
                failures.append("face_with_fewer_than_three_vertices")
                break
            previous = offset
        if any(vertex < 0 or vertex >= point_count for vertex in vertices):
            failures.append("face_point_out_of_range")
        maximum_cell = -1
        previous_owner = -1
        for face_id, cell in enumerate(owner):
            if cell < 0:
                failures.append("negative_owner")
                break
            if face_id < len(neighbour) and cell < previous_owner:
                failures.append("owner_not_upper_triangular_sorted")
                break
            if face_id < len(neighbour):
                previous_owner = cell
            maximum_cell = max(maximum_cell, cell)
            if face_id < len(neighbour) and cell >= neighbour[face_id]:
                failures.append("owner_not_less_than_neighbour")
                break
        for cell in neighbour:
            if cell < 0:
                failures.append("negative_neighbour")
                break
            maximum_cell = max(maximum_cell, cell)
        cell_count = maximum_cell + 1
        patches = parse_boundary(paths["boundary"])
        expected_start = len(neighbour)
        for patch in patches:
            if patch["startFace"] != expected_start:
                failures.append("noncontiguous_boundary_patch")
                break
            expected_start += int(patch["nFaces"])
        if expected_start != face_count:
            failures.append("boundary_patch_face_extent_mismatch")
        edge_diagnostics: dict[str, object] = {"checked": False}
        if args.diagnose_cell_edges:
            edge_diagnostics = diagnose_cell_edges(
                points, offsets, vertices, owner, neighbour, cell_count
            )
            if edge_diagnostics["nonTwoIncidentCellEdgeCount"] != 0:
                failures.append("nonmanifold_cell_edge_incidence")
        hashes = {name: file_sha256(path) for name, path in paths.items()}
        elapsed = time.perf_counter() - start
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            peak_rss_bytes = int(rss)
        else:
            peak_rss_bytes = int(rss * 1024)
        report = {
            "schema": "cartmesh-stage6-independent-binary-reader-v1",
            "status": "pass" if not failures else "fail",
            "case": str(args.case),
            "pointCount": point_count,
            "faceCount": face_count,
            "faceVertexReferenceCount": len(vertices),
            "ownerCount": len(owner),
            "internalFaceCount": len(neighbour),
            "boundaryFaceCount": face_count - len(neighbour),
            "cellCount": cell_count,
            "patches": patches,
            "cellEdgeDiagnostics": edge_diagnostics,
            "failures": failures,
            "sha256": hashes,
            "wallClockSeconds": elapsed,
            "peakRssBytes": peak_rss_bytes,
            "runtimeThreads": 1,
            "python": platform.python_version(),
            "platform": platform.platform(),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(report, ensure_ascii=False))
        return 0 if not failures else 2
    finally:
        # memoryview 必须先离开作用域；进程退出会关闭 mmap，这里只关闭
        # 尚未导出 view 的空诊断路径，避免 BufferError 掩盖真实报告。
        if not mapped_files:
            for mapped in mapped_files:
                mapped.close()


if __name__ == "__main__":
    raise SystemExit(main())
