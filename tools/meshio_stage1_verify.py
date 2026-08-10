#!/usr/bin/env python3
"""使用 meshio 独立读取阶段 1 STL/VTU，并复核表面拓扑、分类和诊断标记。"""

from __future__ import annotations

import argparse
import collections
import dataclasses
import json
import math
import pathlib
import warnings
import xml.etree.ElementTree as ET

import meshio
import numpy as np


Coordinate = tuple[float, float, float]
EdgeKey = tuple[int, int]
CoordinateEdge = frozenset[Coordinate]

ISSUE_NAMES = {
    1: "degenerateTriangle",
    2: "duplicateTriangle",
    3: "boundaryEdge",
    4: "nonManifoldEdge",
    5: "orientationConflict",
    6: "nonManifoldVertex",
    7: "componentOrientationMismatch",
}


@dataclasses.dataclass
class SurfaceInspection:
    summary: dict[str, object]
    points: np.ndarray
    triangles: np.ndarray
    triangle_issues: dict[int, set[int]]
    edge_issues: dict[int, dict[CoordinateEdge, tuple[int, ...]]]
    non_manifold_vertices: dict[Coordinate, int]
    component_orientation_mismatches: dict[Coordinate, int]


def one_cell_block(mesh: meshio.Mesh, cell_type: str) -> np.ndarray:
    blocks = [block.data for block in mesh.cells if block.type == cell_type]
    if len(blocks) != 1 or len(mesh.cells) != 1:
        raise ValueError(f"期望唯一 {cell_type} 单元块，实际为 {[b.type for b in mesh.cells]}")
    return np.asarray(blocks[0])


def normalized_coordinate(values: np.ndarray | list[float]) -> Coordinate:
    coordinate = tuple(0.0 if float(value) == 0.0 else float(value) for value in values)
    if len(coordinate) != 3 or not np.isfinite(coordinate).all():
        raise ValueError("坐标必须包含三个有限分量")
    return coordinate  # type: ignore[return-value]


def canonicalize_surface_vertices(
    source_points: np.ndarray, source_triangles: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    vertex_ids: dict[Coordinate, int] = {}
    points: list[Coordinate] = []
    triangles = np.empty(source_triangles.shape, dtype=np.int64)
    for triangle_id, triangle in enumerate(source_triangles):
        for local_id, source_vertex_id in enumerate(triangle):
            coordinate = normalized_coordinate(source_points[int(source_vertex_id), :3])
            vertex_id = vertex_ids.get(coordinate)
            if vertex_id is None:
                vertex_id = len(points)
                vertex_ids[coordinate] = vertex_id
                points.append(coordinate)
            triangles[triangle_id, local_id] = vertex_id
    return np.asarray(points, dtype=np.float64), triangles


def inspect_surface_arrays(
    source_points: np.ndarray, source_triangles: np.ndarray, reader: str
) -> SurfaceInspection:
    source_points = np.asarray(source_points, dtype=np.float64)
    source_triangles = np.asarray(source_triangles)
    if source_points.ndim != 2 or source_points.shape[1] < 3:
        raise ValueError("STL 点坐标数组必须是至少三列的二维数组")
    if not np.isfinite(source_points[:, :3]).all():
        raise ValueError("STL 点坐标包含非有限值")
    if source_triangles.ndim != 2 or source_triangles.shape[1] != 3:
        raise ValueError("STL 三角形连接数组宽度不是 3")
    if source_triangles.size and (
        source_triangles.min() < 0 or source_triangles.max() >= len(source_points)
    ):
        raise ValueError("STL 三角形引用了范围之外的点")

    points, triangles = canonicalize_surface_vertices(source_points, source_triangles)
    if not len(triangles):
        raise ValueError("STL 不得为空")
    extent = points.max(axis=0) - points.min(axis=0)
    diagonal = float(np.linalg.norm(extent))
    scale = max(diagonal, float(np.finfo(np.float64).tiny))
    suggested_length_tolerance = scale * 1.0e-12
    degenerate_area_tolerance = 0.5 * scale * suggested_length_tolerance

    edge_uses: dict[EdgeKey, list[tuple[int, bool]]] = collections.defaultdict(list)
    canonical_triangles: dict[tuple[int, int, int], int] = {}
    incident_triangles: list[list[int]] = [[] for _ in range(len(points))]
    degenerate_ids: set[int] = set()
    duplicate_ids: set[int] = set()
    valid_ids: list[int] = []
    for triangle_id, triangle in enumerate(triangles):
        triangle_points = points[triangle]
        area_vector = np.cross(
            triangle_points[1] - triangle_points[0],
            triangle_points[2] - triangle_points[0],
        )
        area = 0.5 * float(np.linalg.norm(area_vector))
        edge_lengths = (
            float(np.linalg.norm(triangle_points[1] - triangle_points[0])),
            float(np.linalg.norm(triangle_points[2] - triangle_points[1])),
            float(np.linalg.norm(triangle_points[0] - triangle_points[2])),
        )
        triangle_area_tolerance = 0.5 * max(edge_lengths) * suggested_length_tolerance
        if area <= triangle_area_tolerance:
            degenerate_ids.add(triangle_id)
            continue
        valid_ids.append(triangle_id)
        for vertex_id in triangle:
            incident_triangles[int(vertex_id)].append(triangle_id)
        canonical = tuple(sorted(int(value) for value in triangle))
        if canonical in canonical_triangles:
            duplicate_ids.add(triangle_id)
        else:
            canonical_triangles[canonical] = triangle_id
        ids = [int(value) for value in triangle]
        for start, end in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            key = (min(start, end), max(start, end))
            edge_uses[key].append((triangle_id, start == key[0]))
    parent = list(range(len(triangles)))
    rank = [0] * len(triangles)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def unite(lhs: int, rhs: int) -> None:
        lhs = find(lhs)
        rhs = find(rhs)
        if lhs == rhs:
            return
        if rank[lhs] < rank[rhs] or (rank[lhs] == rank[rhs] and lhs > rhs):
            lhs, rhs = rhs, lhs
        parent[rhs] = lhs
        if rank[lhs] == rank[rhs]:
            rank[lhs] += 1

    incident_edges: list[list[EdgeKey]] = [[] for _ in range(len(points))]
    boundary_edges: dict[CoordinateEdge, tuple[int, ...]] = {}
    non_manifold_edges: dict[CoordinateEdge, tuple[int, ...]] = {}
    orientation_conflicts: dict[CoordinateEdge, tuple[int, ...]] = {}
    for edge, uses in sorted(edge_uses.items()):
        incident_edges[edge[0]].append(edge)
        incident_edges[edge[1]].append(edge)
        for use in uses[1:]:
            unite(uses[0][0], use[0])
        coordinate_edge = frozenset(
            (normalized_coordinate(points[edge[0]]), normalized_coordinate(points[edge[1]]))
        )
        triangle_ids = tuple(use[0] for use in uses)
        if len(uses) == 1:
            boundary_edges[coordinate_edge] = triangle_ids
        elif len(uses) > 2:
            non_manifold_edges[coordinate_edge] = triangle_ids
        elif uses[0][1] == uses[1][1]:
            orientation_conflicts[coordinate_edge] = triangle_ids

    non_manifold_vertices: dict[Coordinate, int] = {}
    for vertex_id, vertex_triangles in enumerate(incident_triangles):
        if not vertex_triangles:
            continue
        local_index = {
            triangle_id: index for index, triangle_id in enumerate(vertex_triangles)
        }
        link: list[list[int]] = [[] for _ in vertex_triangles]
        boundary_edge_count = 0
        invalid_edge_incidence = False
        for edge in incident_edges[vertex_id]:
            uses = edge_uses[edge]
            if len(uses) == 1:
                boundary_edge_count += 1
            elif len(uses) == 2:
                first = local_index[uses[0][0]]
                second = local_index[uses[1][0]]
                if first != second:
                    link[first].append(second)
                    link[second].append(first)
            else:
                invalid_edge_incidence = True

        visited = {0}
        pending = collections.deque([0])
        while pending:
            current = pending.popleft()
            for neighbour in link[current]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    pending.append(neighbour)
        degree_one = sum(len(neighbours) == 1 for neighbours in link)
        degree_valid = all(len(neighbours) <= 2 for neighbours in link)
        closed_link = (
            boundary_edge_count == 0
            and degree_valid
            and all(len(neighbours) == 2 for neighbours in link)
        )
        single_triangle_boundary_link = (
            len(vertex_triangles) == 1
            and boundary_edge_count == 2
            and not link[0]
        )
        open_link = boundary_edge_count == 2 and degree_valid and degree_one == 2
        manifold_vertex = (
            not invalid_edge_incidence
            and len(visited) == len(link)
            and (closed_link or open_link or single_triangle_boundary_link)
        )
        if not manifold_vertex:
            non_manifold_vertices[normalized_coordinate(points[vertex_id])] = len(
                vertex_triangles
            )

    roots = sorted({find(triangle_id) for triangle_id in valid_ids})
    component_index = {root: index for index, root in enumerate(roots)}
    component_triangles: list[list[int]] = [[] for _ in roots]
    for triangle_id in valid_ids:
        component_triangles[component_index[find(triangle_id)]].append(triangle_id)
    connected_components = len(component_triangles)
    manifold = not non_manifold_edges and not non_manifold_vertices
    closed = not boundary_edges and manifold
    consistently_oriented = not orientation_conflicts

    def component_contains_point(triangle_ids: list[int], point: np.ndarray) -> bool | None:
        component_points = points[triangles[np.asarray(triangle_ids, dtype=np.int64)].reshape(-1)]
        if np.any(point < component_points.min(axis=0)) or np.any(point > component_points.max(axis=0)):
            return False
        solid_angle = 0.0
        for triangle_id in triangle_ids:
            vectors = points[triangles[triangle_id]].astype(np.longdouble) - point.astype(
                np.longdouble
            )
            lengths = [float(np.sqrt(np.dot(vector, vector))) for vector in vectors]
            if any(length == 0.0 for length in lengths):
                return None
            numerator = float(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))
            denominator = (
                lengths[0] * lengths[1] * lengths[2]
                + float(np.dot(vectors[0], vectors[1])) * lengths[2]
                + float(np.dot(vectors[1], vectors[2])) * lengths[0]
                + float(np.dot(vectors[2], vectors[0])) * lengths[1]
            )
            if numerator == 0.0 and denominator == 0.0:
                return None
            solid_angle += 2.0 * math.atan2(numerator, denominator)
        return abs(solid_angle) > 2.0 * math.pi

    components: list[dict[str, object]] = []
    total_signed_volume = np.longdouble(0.0)
    for component_id, triangle_ids in enumerate(component_triangles):
        reference = points[triangles[triangle_ids[0], 0]].astype(np.longdouble)
        signed_six_volume = np.longdouble(0.0)
        for triangle_id in triangle_ids:
            shifted = points[triangles[triangle_id]].astype(np.longdouble) - reference
            signed_six_volume += np.dot(shifted[0], np.cross(shifted[1], shifted[2]))
        signed_volume = signed_six_volume / np.longdouble(6.0)
        total_signed_volume += signed_volume
        sample_position = points[triangles[triangle_ids[0]]].mean(axis=0)
        components.append(
            {
                "componentId": component_id,
                "triangleCount": len(triangle_ids),
                "signedVolume": float(signed_volume),
                "nestingDepth": 0,
                "expectedOrientationSign": 1,
                "orientationChecked": False,
                "orientationMatchesNesting": False,
                "samplePosition": sample_position.tolist(),
            }
        )

    component_orientation_mismatches: dict[Coordinate, int] = {}
    can_check_components = (
        closed and consistently_oriented and not degenerate_ids and not duplicate_ids
    )
    if can_check_components:
        for component_id, component in enumerate(components):
            sample = np.asarray(component["samplePosition"], dtype=np.float64)
            nesting_depth = 0
            orientation_checked = True
            for other_id, triangle_ids in enumerate(component_triangles):
                if other_id == component_id:
                    continue
                contains = component_contains_point(triangle_ids, sample)
                if contains is None:
                    orientation_checked = False
                    break
                nesting_depth += int(contains)
            expected_sign = 1 if nesting_depth % 2 == 0 else -1
            matches = orientation_checked and (
                float(component["signedVolume"]) > 0.0
                if expected_sign > 0
                else float(component["signedVolume"]) < 0.0
            )
            component["nestingDepth"] = nesting_depth
            component["expectedOrientationSign"] = expected_sign
            component["orientationChecked"] = orientation_checked
            component["orientationMatchesNesting"] = matches
            if not matches:
                component_orientation_mismatches[normalized_coordinate(sample)] = component_id

    valid_for_stage1 = (
        closed
        and manifold
        and consistently_oriented
        and not degenerate_ids
        and not duplicate_ids
        and not non_manifold_vertices
    )
    material_volume = sum(
        int(component["expectedOrientationSign"])
        * abs(float(component["signedVolume"]))
        for component in components
    )
    summary: dict[str, object] = {
        "reader": reader,
        "triangleCount": int(len(triangles)),
        "uniqueVertexCount": int(len(points)),
        "uniqueEdgeCount": int(len(edge_uses)),
        "degenerateTriangleCount": len(degenerate_ids),
        "duplicateTriangleCount": len(duplicate_ids),
        "boundaryEdgeCount": len(boundary_edges),
        "nonManifoldEdgeCount": len(non_manifold_edges),
        "nonManifoldVertexCount": len(non_manifold_vertices),
        "orientationConflictEdgeCount": len(orientation_conflicts),
        "connectedComponentCount": connected_components,
        "componentOrientationMismatchCount": len(component_orientation_mismatches),
        "suggestedLengthTolerance": suggested_length_tolerance,
        "degenerateAreaTolerance": degenerate_area_tolerance,
        "signedVolume": float(total_signed_volume),
        "materialVolume": float(material_volume),
        "closed": closed,
        "manifold": manifold,
        "consistentlyOriented": consistently_oriented,
        "validForStage1Classification": valid_for_stage1,
        "components": components,
    }
    return SurfaceInspection(
        summary=summary,
        points=points,
        triangles=triangles,
        triangle_issues={1: degenerate_ids, 2: duplicate_ids},
        edge_issues={
            3: boundary_edges,
            4: non_manifold_edges,
            5: orientation_conflicts,
        },
        non_manifold_vertices=non_manifold_vertices,
        component_orientation_mismatches=component_orientation_mismatches,
    )


def inspect_surface(path: pathlib.Path) -> SurfaceInspection:
    # meshio 5.3.5 probes an ASCII STL's first four facet bytes as a binary
    # triangle count.  NumPy can warn while multiplying that intentionally
    # nonsensical probe value; the subsequent file-size check still selects
    # the ASCII reader.  Suppress only this known format-probe warning.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="overflow encountered in scalar multiply",
            category=RuntimeWarning,
            module=r"meshio\.stl\._stl",
        )
        surface = meshio.read(path)
    triangles = one_cell_block(surface, "triangle")
    return inspect_surface_arrays(
        np.asarray(surface.points, dtype=np.float64),
        triangles,
        f"meshio {meshio.__version__}",
    )


def analytic_cube_classification(
    points: np.ndarray, cells: np.ndarray, cube_minimum: np.ndarray, cube_maximum: np.ndarray
) -> np.ndarray:
    expected = np.empty(len(cells), dtype=np.int64)
    for cell_id, cell in enumerate(cells):
        cell_points = points[cell]
        cell_minimum = cell_points.min(axis=0)
        cell_maximum = cell_points.max(axis=0)
        intersects_surface = False
        for axis in range(3):
            other_axes = [value for value in range(3) if value != axis]
            overlaps_face_rectangle = all(
                cell_minimum[other] <= cube_maximum[other]
                and cell_maximum[other] >= cube_minimum[other]
                for other in other_axes
            )
            if not overlaps_face_rectangle:
                continue
            if (
                cell_minimum[axis] <= cube_minimum[axis] <= cell_maximum[axis]
                or cell_minimum[axis] <= cube_maximum[axis] <= cell_maximum[axis]
            ):
                intersects_surface = True
                break
        if intersects_surface:
            expected[cell_id] = 2
            continue
        center = 0.5 * (cell_minimum + cell_maximum)
        expected[cell_id] = 1 if np.all((center >= cube_minimum) & (center <= cube_maximum)) else 0
    return expected


def inspect_volume_mesh(
    path: pathlib.Path, analytic_cube: tuple[np.ndarray, np.ndarray] | None
) -> dict[str, object]:
    mesh = meshio.read(path)
    cells = one_cell_block(mesh, "hexahedron")
    points = np.asarray(mesh.points, dtype=np.float64)
    if cells.ndim != 2 or cells.shape[1] != 8:
        raise ValueError("VTU 六面体连接数组宽度不是 8")
    if cells.size and (cells.min() < 0 or cells.max() >= len(mesh.points)):
        raise ValueError("VTU 六面体引用了范围之外的点")
    arrays = mesh.cell_data.get("stl_cell_classification")
    if arrays is None or len(arrays) != 1:
        raise ValueError("VTU 缺少 stl_cell_classification 单元字段")
    values = np.asarray(arrays[0])
    if len(values) != len(cells):
        raise ValueError("分类字段长度与六面体数量不一致")
    rounded = np.rint(values).astype(np.int64)
    if not np.array_equal(values, rounded) or not np.isin(rounded, [0, 1, 2, 3]).all():
        raise ValueError("分类字段包含图例之外的值")
    counts = np.bincount(rounded, minlength=4)
    analytic_match = None
    if analytic_cube is not None:
        expected = analytic_cube_classification(points, cells, *analytic_cube)
        mismatches = np.flatnonzero(expected != rounded)
        if len(mismatches):
            first = int(mismatches[0])
            raise ValueError(
                f"解析立方体逐单元分类不一致：错误数={len(mismatches)}，首个单元={first}，"
                f"期望={expected[first]}，实际={rounded[first]}"
            )
        analytic_match = int(len(cells))

    center_arrays = mesh.cell_data.get("inside_stl_center_sample")
    if center_arrays is None or len(center_arrays) != 1:
        raise ValueError("VTU 缺少 inside_stl_center_sample 单元字段")
    center_inside = np.asarray(center_arrays[0])
    if len(center_inside) != len(cells) or not np.isin(center_inside, [0.0, 1.0]).all():
        raise ValueError("中心采样字段必须与单元数一致且只包含 0/1")
    return {
        "pointCount": int(len(mesh.points)),
        "cellCount": int(len(cells)),
        "classificationCounts": {
            "outside": int(counts[0]),
            "inside": int(counts[1]),
            "intersected": int(counts[2]),
            "conflict": int(counts[3]),
        },
        "centerInsideCount": int(np.count_nonzero(center_inside == 1.0)),
        "analyticCubeCellMatches": analytic_match,
        "cellDataFields": sorted(mesh.cell_data),
    }


def coordinate_edge_from_report(example: dict[str, object]) -> CoordinateEdge:
    return frozenset(
        (
            normalized_coordinate(example["start"]),  # type: ignore[arg-type]
            normalized_coordinate(example["end"]),  # type: ignore[arg-type]
        )
    )


def validate_examples(
    diagnostics: dict[str, object], surface: SurfaceInspection
) -> tuple[np.ndarray, np.ndarray]:
    expected_codes: list[int] = []
    expected_positions: list[Coordinate] = []

    triangle_specs = (
        (1, "degenerateTriangleExamples", "degenerateTriangleCount"),
        (2, "duplicateTriangleExamples", "duplicateTriangleCount"),
    )
    for code, example_name, count_name in triangle_specs:
        examples = diagnostics[example_name]
        if not isinstance(examples, list):
            raise ValueError(f"{example_name} 必须是数组")
        if int(diagnostics[count_name]) > 0 and not examples:
            raise ValueError(f"{count_name} 非零但没有位置示例")
        if len(examples) > int(diagnostics[count_name]):
            raise ValueError(f"{example_name} 数量超过问题总数")
        for raw_triangle_id in examples:
            triangle_id = int(raw_triangle_id)
            if triangle_id not in surface.triangle_issues[code]:
                raise ValueError(f"{example_name} 含非独立检查结果中的三角形 {triangle_id}")
            centroid = surface.points[surface.triangles[triangle_id]].mean(axis=0)
            expected_codes.append(code)
            expected_positions.append(normalized_coordinate(centroid))

    edge_specs = (
        (3, "boundaryEdgeExamples", "boundaryEdgeCount"),
        (4, "nonManifoldEdgeExamples", "nonManifoldEdgeCount"),
        (5, "orientationConflictExamples", "orientationConflictEdgeCount"),
    )
    for code, example_name, count_name in edge_specs:
        examples = diagnostics[example_name]
        if not isinstance(examples, list):
            raise ValueError(f"{example_name} 必须是数组")
        if int(diagnostics[count_name]) > 0 and not examples:
            raise ValueError(f"{count_name} 非零但没有位置示例")
        if len(examples) > int(diagnostics[count_name]):
            raise ValueError(f"{example_name} 数量超过问题总数")
        for example in examples:
            if not isinstance(example, dict):
                raise ValueError(f"{example_name} 的元素必须是对象")
            edge = coordinate_edge_from_report(example)
            uses = surface.edge_issues[code].get(edge)
            if uses is None:
                raise ValueError(f"{example_name} 含不属于该问题类型的边")
            expected_first = uses[0]
            expected_second = uses[1] if len(uses) > 1 else uses[0]
            if int(example["firstTriangle"]) != expected_first or int(
                example["secondTriangle"]
            ) != expected_second:
                raise ValueError(f"{example_name} 的关联三角形与独立拓扑不一致")
            start = np.asarray(example["start"], dtype=np.float64)
            end = np.asarray(example["end"], dtype=np.float64)
            expected_codes.append(code)
            expected_positions.append(normalized_coordinate(0.5 * (start + end)))

    vertex_examples = diagnostics["nonManifoldVertexExamples"]
    if not isinstance(vertex_examples, list):
        raise ValueError("nonManifoldVertexExamples 必须是数组")
    if int(diagnostics["nonManifoldVertexCount"]) > 0 and not vertex_examples:
        raise ValueError("nonManifoldVertexCount 非零但没有位置示例")
    if len(vertex_examples) > int(diagnostics["nonManifoldVertexCount"]):
        raise ValueError("nonManifoldVertexExamples 数量超过问题总数")
    for example in vertex_examples:
        if not isinstance(example, dict):
            raise ValueError("nonManifoldVertexExamples 的元素必须是对象")
        position = normalized_coordinate(example["position"])  # type: ignore[arg-type]
        incident_count = surface.non_manifold_vertices.get(position)
        if incident_count is None:
            raise ValueError("报告中的非流形顶点不在独立 link 检查结果中")
        if int(example["incidentTriangleCount"]) != incident_count:
            raise ValueError("非流形顶点 incidentTriangleCount 与独立检查不一致")
        expected_codes.append(6)
        expected_positions.append(position)

    components = diagnostics["components"]
    if not isinstance(components, list):
        raise ValueError("components 必须是数组")
    for component in components:
        if not isinstance(component, dict):
            raise ValueError("components 的元素必须是对象")
        if int(diagnostics["componentOrientationMismatchCount"]) > 0 and not bool(
            component["orientationMatchesNesting"]
        ):
            position = normalized_coordinate(component["samplePosition"])  # type: ignore[arg-type]
            component_id = surface.component_orientation_mismatches.get(position)
            if component_id is None or component_id != int(component["componentId"]):
                raise ValueError("分量方向 marker 与独立嵌套检查不一致")
            expected_codes.append(7)
            expected_positions.append(position)

    positions = np.asarray(expected_positions, dtype=np.float64)
    if not expected_positions:
        positions = np.empty((0, 3), dtype=np.float64)
    return np.asarray(expected_codes, dtype=np.int64), positions


def compare_surface_diagnostics(
    diagnostics: dict[str, object], surface: SurfaceInspection
) -> tuple[np.ndarray, np.ndarray]:
    fields = (
        "triangleCount",
        "uniqueVertexCount",
        "uniqueEdgeCount",
        "degenerateTriangleCount",
        "duplicateTriangleCount",
        "boundaryEdgeCount",
        "nonManifoldEdgeCount",
        "nonManifoldVertexCount",
        "orientationConflictEdgeCount",
        "connectedComponentCount",
        "componentOrientationMismatchCount",
        "closed",
        "manifold",
        "consistentlyOriented",
        "validForStage1Classification",
    )
    for name in fields:
        if diagnostics[name] != surface.summary[name]:
            raise ValueError(
                f"报告与独立 STL 检查的 {name} 不一致："
                f"{diagnostics[name]} != {surface.summary[name]}"
            )
    for name in (
        "suggestedLengthTolerance",
        "degenerateAreaTolerance",
        "signedVolume",
        "materialVolume",
    ):
        if not np.isclose(
            float(diagnostics[name]),
            float(surface.summary[name]),
            rtol=1.0e-12,
            atol=float(np.finfo(np.float64).tiny),
        ):
            raise ValueError(f"报告与独立 STL 检查的 {name} 不一致")
    reported_components = diagnostics.get("components")
    expected_components = surface.summary["components"]
    if not isinstance(reported_components, list) or not isinstance(expected_components, list):
        raise ValueError("分量诊断必须是数组")
    if len(reported_components) != len(expected_components):
        raise ValueError("报告与独立检查的分量数量不一致")
    exact_component_fields = (
        "componentId",
        "triangleCount",
        "nestingDepth",
        "expectedOrientationSign",
        "orientationChecked",
        "orientationMatchesNesting",
    )
    for reported, expected in zip(reported_components, expected_components):
        if not isinstance(reported, dict) or not isinstance(expected, dict):
            raise ValueError("分量诊断元素必须是对象")
        for name in exact_component_fields:
            if reported[name] != expected[name]:
                raise ValueError(f"分量诊断字段 {name} 与独立检查不一致")
        if not np.isclose(
            float(reported["signedVolume"]),
            float(expected["signedVolume"]),
            rtol=1.0e-12,
            atol=float(np.finfo(np.float64).tiny),
        ):
            raise ValueError("分量 signedVolume 与独立检查不一致")
        if not np.allclose(
            np.asarray(reported["samplePosition"], dtype=np.float64),
            np.asarray(expected["samplePosition"], dtype=np.float64),
            rtol=0.0,
            atol=float(diagnostics["suggestedLengthTolerance"]),
        ):
            raise ValueError("分量 samplePosition 与独立检查不一致")
    return validate_examples(diagnostics, surface)


def numbers(element: ET.Element, converter: type[int] | type[float]) -> list[int] | list[float]:
    return [converter(token) for token in (element.text or "").split()]


def named_arrays(parent: ET.Element) -> dict[str, ET.Element]:
    return {
        array.attrib["Name"]: array
        for array in parent.findall("DataArray")
        if "Name" in array.attrib
    }


def inspect_diagnostic_marker(
    path: pathlib.Path, expected_codes: np.ndarray, expected_positions: np.ndarray, tolerance: float
) -> dict[str, object]:
    root = ET.parse(path).getroot()
    if root.tag != "VTKFile" or root.attrib.get("type") != "PolyData":
        raise ValueError("诊断 marker 根节点必须是 type=PolyData 的 VTKFile")
    piece = root.find("./PolyData/Piece")
    if piece is None:
        raise ValueError("诊断 marker 缺少 PolyData/Piece")
    point_count = int(piece.attrib["NumberOfPoints"])
    vertex_count = int(piece.attrib["NumberOfVerts"])
    if point_count != vertex_count:
        raise ValueError("诊断 marker 必须由一一对应的点和 VTK_VERTEX 单元组成")
    if any(int(piece.attrib.get(name, "0")) != 0 for name in ("NumberOfLines", "NumberOfStrips", "NumberOfPolys")):
        raise ValueError("诊断 marker 不得包含线、条带或多边形单元")

    point_array = piece.find("./Points/DataArray")
    if point_array is None or point_array.attrib.get("NumberOfComponents") != "3":
        raise ValueError("诊断 marker 缺少三分量点坐标")
    coordinates = np.asarray(numbers(point_array, float), dtype=np.float64)
    if len(coordinates) != point_count * 3 or not np.isfinite(coordinates).all():
        raise ValueError("诊断 marker 点坐标数量错误或包含非有限值")
    coordinates = coordinates.reshape((-1, 3))

    vertices = piece.find("Verts")
    if vertices is None:
        raise ValueError("诊断 marker 缺少 Verts")
    vertex_arrays = named_arrays(vertices)
    if set(vertex_arrays) != {"connectivity", "offsets"}:
        raise ValueError("诊断 marker Verts 必须只包含 connectivity 和 offsets")
    connectivity = numbers(vertex_arrays["connectivity"], int)
    offsets = numbers(vertex_arrays["offsets"], int)
    if connectivity != list(range(point_count)) or offsets != list(range(1, point_count + 1)):
        raise ValueError("诊断 marker 的 VTK_VERTEX 连接或偏移不连续")

    cell_data = piece.find("CellData")
    if cell_data is None:
        raise ValueError("诊断 marker 缺少 CellData")
    arrays = named_arrays(cell_data)
    if "issue_code" not in arrays:
        raise ValueError("诊断 marker 缺少 issue_code")
    issue_array = arrays["issue_code"]
    if issue_array.attrib.get("type") != "UInt8":
        raise ValueError("diagnostic marker 的 issue_code 必须是 UInt8")
    issue_codes = np.asarray(numbers(issue_array, int), dtype=np.int64)
    if len(issue_codes) != point_count or not np.isin(issue_codes, list(ISSUE_NAMES)).all():
        raise ValueError("diagnostic marker 的 issue_code 数量错误或超出图例")
    if not np.array_equal(issue_codes, expected_codes):
        raise ValueError(
            f"diagnostic marker 的 issue_code 与报告示例不一致："
            f"实际={issue_codes.tolist()}，期望={expected_codes.tolist()}"
        )
    if coordinates.shape != expected_positions.shape or not np.allclose(
        coordinates, expected_positions, rtol=0.0, atol=max(tolerance, 1.0e-14)
    ):
        raise ValueError("diagnostic marker 的位置与报告中的缺陷位置不一致")
    code_counts = collections.Counter(int(value) for value in issue_codes)
    return {
        "reader": "Python xml.etree.ElementTree VTK PolyData checker",
        "path": str(path),
        "pointCount": point_count,
        "vertexCellCount": vertex_count,
        "issueCodeLegend": {str(code): name for code, name in ISSUE_NAMES.items()},
        "issueCodeCounts": {
            ISSUE_NAMES[code]: code_counts.get(code, 0) for code in sorted(ISSUE_NAMES)
        },
        "positionsMatchReport": True,
        "status": "pass",
    }


def compare_report(
    report_path: pathlib.Path, surface: SurfaceInspection, volume_mesh: dict[str, object]
) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    compare_surface_diagnostics(report["surfaceDiagnostics"], surface)
    if report["pointCount"] != volume_mesh["pointCount"]:
        raise ValueError("报告点数与 meshio 读取结果不一致")
    if report["cellCount"] != volume_mesh["cellCount"]:
        raise ValueError("报告单元数与 meshio 读取结果不一致")
    if report["classificationCounts"] != volume_mesh["classificationCounts"]:
        raise ValueError("报告分类计数与 VTU 字段不一致")
    if report["classificationCounts"]["conflict"] != 0:
        raise ValueError("阶段 1 外部验证案例存在分类冲突")
    if report["centerPointCounts"]["inside"] != volume_mesh["centerInsideCount"]:
        raise ValueError("报告中心内部计数与 VTU 中心采样字段不一致")
    return {
        "projectStage": report["projectStage"],
        "status": report["status"],
        "resultHashFnv1a64": report["resultHashFnv1a64"],
    }


def inspect_rejected_report(
    report_path: pathlib.Path, marker_path: pathlib.Path | None, surface: SurfaceInspection
) -> tuple[dict[str, object], dict[str, object]]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("projectStage") != 1 or report.get("status") != "rejected_invalid_surface":
        raise ValueError("无效几何报告必须是阶段 1 rejected_invalid_surface")
    expected_codes, expected_positions = compare_surface_diagnostics(
        report["surfaceDiagnostics"], surface
    )
    if report["surfaceDiagnostics"]["validForStage1Classification"]:
        raise ValueError("被拒绝的表面不得标记为可用于阶段 1 分类")
    reported_marker = report.get("diagnosticMarkerVtp")
    if marker_path is None:
        if not isinstance(reported_marker, str) or not reported_marker:
            raise ValueError("无效几何报告缺少 diagnosticMarkerVtp")
        marker_path = pathlib.Path(reported_marker)
    if not marker_path.is_absolute() and not marker_path.exists():
        candidate = report_path.parent / marker_path.name
        if candidate.exists():
            marker_path = candidate
    marker = inspect_diagnostic_marker(
        marker_path,
        expected_codes,
        expected_positions,
        float(surface.summary["suggestedLengthTolerance"]),
    )
    return (
        {
            "projectStage": report["projectStage"],
            "status": report["status"],
            "diagnosticMarkerVtp": str(marker_path),
        },
        marker,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("surface", type=pathlib.Path)
    parser.add_argument("volume_mesh", type=pathlib.Path, nargs="?")
    parser.add_argument("report", type=pathlib.Path, nargs="?")
    parser.add_argument(
        "--axis-aligned-cube",
        type=float,
        nargs=6,
        metavar=("XMIN", "YMIN", "ZMIN", "XMAX", "YMAX", "ZMAX"),
        help="用解析轴对齐立方体逐单元独立复核分类",
    )
    parser.add_argument(
        "--diagnostic-report",
        type=pathlib.Path,
        help="验证 rejected_invalid_surface JSON；该模式不需要 VTU",
    )
    parser.add_argument(
        "--diagnostic-marker",
        type=pathlib.Path,
        help="显式指定诊断 marker VTP；默认从 rejected report 读取",
    )
    parser.add_argument("--output", type=pathlib.Path, help="另存机器可读验证结果")
    arguments = parser.parse_args()
    try:
        surface = inspect_surface(arguments.surface)
        if arguments.diagnostic_report is not None:
            if arguments.volume_mesh is not None or arguments.report is not None:
                raise ValueError("诊断模式不得再提供 volume_mesh 或普通 report 位置参数")
            rejected_report, marker = inspect_rejected_report(
                arguments.diagnostic_report, arguments.diagnostic_marker, surface
            )
            result: dict[str, object] = {
                "surface": surface.summary,
                "diagnosticReport": rejected_report,
                "diagnosticMarker": marker,
                "status": "pass",
            }
        else:
            if arguments.volume_mesh is None or arguments.report is None:
                raise ValueError("正常分类模式必须提供 volume_mesh 和 report")
            if arguments.diagnostic_marker is not None:
                raise ValueError("--diagnostic-marker 只能与 --diagnostic-report 一起使用")
            analytic_cube = None
            if arguments.axis_aligned_cube:
                values = np.asarray(arguments.axis_aligned_cube, dtype=np.float64)
                analytic_cube = (values[:3], values[3:])
                if not np.all(analytic_cube[1] > analytic_cube[0]):
                    raise ValueError("解析立方体最大坐标必须逐轴大于最小坐标")
            volume_mesh = inspect_volume_mesh(arguments.volume_mesh, analytic_cube)
            report = compare_report(arguments.report, surface, volume_mesh)
            result = {
                "surface": surface.summary,
                "volumeMesh": volume_mesh,
                "report": report,
                "status": "pass",
            }
        payload = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
        if arguments.output:
            arguments.output.parent.mkdir(parents=True, exist_ok=True)
            arguments.output.write_text(payload + "\n", encoding="utf-8")
        print(payload)
        return 0
    except (ET.ParseError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(
            json.dumps(
                {"status": "fail", "error": str(error)}, ensure_ascii=False, sort_keys=True
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
