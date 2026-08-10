#!/usr/bin/env python3
"""独立复核阶段三凸 STL Cut-cell 几何、拓扑和 VTK 输出。"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import warnings
import xml.etree.ElementTree as ET

import meshio
import numpy as np


def polygon_area_vector(vertices: np.ndarray) -> np.ndarray:
    if len(vertices) < 3:
        return np.zeros(3)
    reference = vertices[0]
    return 0.5 * sum(
        (np.cross(vertices[i] - reference, vertices[i + 1] - reference)
         for i in range(1, len(vertices) - 1)),
        start=np.zeros(3),
    )


def boundary_edge_imbalance(cell: dict, tolerance: float) -> int:
    points: list[np.ndarray] = []
    edges: list[tuple[int, int]] = []

    def point_id(value: list[float]) -> int:
        point = np.asarray(value, dtype=np.float64)
        for index, existing in enumerate(points):
            if np.linalg.norm(point - existing) <= tolerance:
                return index
        points.append(point)
        return len(points) - 1

    def append_loop(loop: list[list[float]]) -> None:
        if len(loop) < 3:
            return
        ids = [point_id(value) for value in loop]
        edges.extend(
            (ids[index], ids[(index + 1) % len(ids)])
            for index in range(len(ids))
            if ids[index] != ids[(index + 1) % len(ids)]
        )

    for face in cell["cartesianFaces"]:
        for loop in face["orientedBoundaryLoops"]:
            append_loop(loop)
    for face in cell["embeddedBoundaryFaces"]:
        append_loop(face["vertices"])

    balances: dict[tuple[int, int], int] = {}
    for first_id, second_id in edges:
        first = points[first_id]
        second = points[second_id]
        direction = second - first
        length_squared = float(np.dot(direction, direction))
        split = [(0.0, first_id), (1.0, second_id)]
        for candidate_id, candidate in enumerate(points):
            if candidate_id in (first_id, second_id):
                continue
            parameter = float(np.dot(candidate - first, direction) / length_squared)
            if 0.0 < parameter < 1.0:
                projection = first + parameter * direction
                if np.linalg.norm(candidate - projection) <= tolerance:
                    split.append((parameter, candidate_id))
        split.sort()
        for (_, a), (_, b) in zip(split, split[1:]):
            if a == b:
                continue
            key = (min(a, b), max(a, b))
            balances[key] = balances.get(key, 0) + (1 if a < b else -1)
    return sum(balance != 0 for balance in balances.values())


def read_vtp(path: Path) -> tuple[np.ndarray, list[np.ndarray], dict[str, np.ndarray]]:
    root = ET.parse(path).getroot()
    piece = root.find("./PolyData/Piece")
    if piece is None:
        raise ValueError("VTP 缺少 PolyData/Piece")
    point_array = piece.find("./Points/DataArray")
    if point_array is None:
        raise ValueError("VTP 缺少点坐标")
    points = np.fromstring(point_array.text or "", sep=" ", dtype=np.float64).reshape(-1, 3)
    connectivity_node = piece.find("./Polys/DataArray[@Name='connectivity']")
    offsets_node = piece.find("./Polys/DataArray[@Name='offsets']")
    if connectivity_node is None or offsets_node is None:
        raise ValueError("VTP 缺少 polygon connectivity/offsets")
    connectivity = np.fromstring(connectivity_node.text or "", sep=" ", dtype=np.int64)
    offsets = np.fromstring(offsets_node.text or "", sep=" ", dtype=np.int64)
    polygons: list[np.ndarray] = []
    begin = 0
    for offset in offsets:
        polygons.append(connectivity[begin:int(offset)])
        begin = int(offset)
    if begin != len(connectivity):
        raise ValueError("VTP polygon offsets 未覆盖完整 connectivity")
    cell_data: dict[str, np.ndarray] = {}
    for node in piece.findall("./CellData/DataArray"):
        name = node.attrib.get("Name")
        if not name:
            continue
        components = int(node.attrib.get("NumberOfComponents", "1"))
        dtype = np.uint64 if node.attrib.get("type") == "UInt64" else np.float64
        values = np.fromstring(node.text or "", sep=" ", dtype=dtype)
        if components != 1:
            values = values.reshape(-1, components)
        cell_data[name] = values
    return points, polygons, cell_data


def flatten_cell_data(mesh: meshio.Mesh, name: str) -> np.ndarray:
    if name not in mesh.cell_data:
        raise ValueError(f"VTU 缺少单元字段 {name}")
    return np.concatenate([np.asarray(block) for block in mesh.cell_data[name]])


def expected_cube_cell(
    cell_id: int,
    dimensions: np.ndarray,
    domain_minimum: np.ndarray,
    spacing: np.ndarray,
    solid_minimum: np.ndarray,
    solid_maximum: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    nx, ny, _ = (int(value) for value in dimensions)
    k, remainder = divmod(cell_id, nx * ny)
    j, i = divmod(remainder, nx)
    cell_minimum = domain_minimum + spacing * np.array([i, j, k])
    cell_maximum = cell_minimum + spacing
    overlap_minimum = np.maximum(cell_minimum, solid_minimum)
    overlap_maximum = np.minimum(cell_maximum, solid_maximum)
    overlap_extent = np.maximum(overlap_maximum - overlap_minimum, 0.0)
    solid_volume = float(np.prod(overlap_extent))
    cell_volume = float(np.prod(spacing))
    fluid_volume = cell_volume - solid_volume
    cell_centroid = 0.5 * (cell_minimum + cell_maximum)
    if fluid_volume > 0.0 and solid_volume > 0.0:
        solid_centroid = 0.5 * (overlap_minimum + overlap_maximum)
        fluid_centroid = (cell_volume * cell_centroid - solid_volume * solid_centroid) / fluid_volume
    else:
        fluid_centroid = cell_centroid
    return fluid_volume, fluid_centroid, cell_minimum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--boundary", type=Path, required=True)
    parser.add_argument("--polyhedra", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--shape", choices=("cube", "l_prism", "none"), default="none")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    report = json.loads(arguments.report.read_text())
    geometry = json.loads(arguments.geometry.read_text())
    failures: list[str] = []
    tolerance = 2.0e-12

    if report.get("stage3GeometryTopologyComplete") is not True:
        failures.append("阶段三几何拓扑不变量未通过")
    if report.get("stage3Complete") is not False:
        failures.append("外部 CFD checker 尚未通过时不得标记 stage3Complete=true")
    if (report.get("solverReadyCutCellMesh") is not False or
            report.get("externalCfdCheckerAccepted") is not False):
        failures.append("缺少外部 CFD checker 时不得标记为求解器可用")
    if geometry.get("solverReadyCutCellMesh") is not False:
        failures.append("几何文件不得在完整求解器体网格尚未导出时标记 solverReady")

    dimensions = np.asarray(report["dimensions"], dtype=np.int64)
    spacing = np.asarray(report["spacing"], dtype=np.float64)
    domain_minimum = np.asarray(report["domain"]["minimum"], dtype=np.float64)
    domain_maximum = np.asarray(report["domain"]["maximum"], dtype=np.float64)
    cell_volume = float(np.prod(spacing))
    background_count = int(np.prod(dimensions))

    vtu = meshio.read(arguments.mesh)
    vtk_cell_count = sum(len(block.data) for block in vtu.cells)
    if vtk_cell_count != background_count:
        failures.append("VTU 背景单元数与报告不一致")
    fractions = flatten_cell_data(vtu, "fluid_volume_fraction").astype(np.float64)
    cut_flags = flatten_cell_data(vtu, "cut_cell").astype(np.float64)
    small_flags = flatten_cell_data(vtu, "small_cut_cell").astype(np.float64)
    if np.any(fractions < -tolerance) or np.any(fractions > 1.0 + tolerance):
        failures.append("VTU 流体体积分数超出 [0,1]")
    vtu_fluid_volume = float(np.sum(fractions) * cell_volume)
    if not math.isclose(vtu_fluid_volume, report["totalFluidVolume"], abs_tol=tolerance):
        failures.append("VTU 体积分数积分与报告总流体体积不一致")
    if int(np.count_nonzero(cut_flags == 1.0)) != report["cutCellCount"]:
        failures.append("VTU Cut-cell 数与报告不一致")
    small_threshold = float(report["smallCellThreshold"])
    expected_small = (cut_flags == 1.0) & (fractions < small_threshold)
    if not np.array_equal(small_flags == 1.0, expected_small):
        failures.append("VTU 小 Cut-cell 标记与体积分数阈值不一致")
    if int(np.count_nonzero(expected_small)) != int(report["smallCutCellCount"]):
        failures.append("VTU 小 Cut-cell 数与报告不一致")
    cut_fractions = fractions[cut_flags == 1.0]
    expected_minimum_fraction = (
        float(np.min(cut_fractions)) if len(cut_fractions) else 1.0
    )
    if not math.isclose(
        expected_minimum_fraction,
        float(report["minimumCutCellVolumeFraction"]),
        abs_tol=tolerance,
    ):
        failures.append("最小 Cut-cell 体积分数与 VTU 不一致")
    reported_small_ids = {
        int(cell["backgroundCellId"]) for cell in report["smallCutCells"]
    }
    expected_small_ids = set(np.flatnonzero(expected_small).tolist())
    if reported_small_ids != expected_small_ids:
        failures.append("小 Cut-cell 位置列表与 VTU 标记不一致")

    cells = geometry["fluidCells"]
    cell_by_id = {int(cell["backgroundCellId"]): cell for cell in cells}
    if len(cell_by_id) != len(cells) or len(cells) != report["fluidCellCount"]:
        failures.append("几何 JSON 流体单元 ID 不唯一或计数不一致")
    for small_cell in report["smallCutCells"]:
        background_id = int(small_cell["backgroundCellId"])
        geometry_cell = cell_by_id.get(background_id)
        if geometry_cell is None:
            failures.append(f"小 Cut-cell {background_id} 不在几何 JSON 中")
            continue
        expected_boundary_ids = sorted({
            int(face["boundaryId"])
            for face in geometry_cell["embeddedBoundaryFaces"]
        })
        if expected_boundary_ids != list(small_cell["boundaryIds"]):
            failures.append(f"小 Cut-cell {background_id} 的边界 ID 不一致")
        if not np.allclose(
            geometry_cell["centroid"], small_cell["centroid"], atol=tolerance
        ):
            failures.append(f"小 Cut-cell {background_id} 的位置不一致")
    geometry_volume = sum(float(cell["volume"]) for cell in cells)
    if not math.isclose(geometry_volume, report["totalFluidVolume"], abs_tol=tolerance):
        failures.append("几何 JSON 单元体积和与报告不一致")

    maximum_recomputed_closure = 0.0
    polygon_area_mismatches = 0
    boundary_edge_failure_cells = 0
    maximum_boundary_edge_imbalance = 0
    for cell in cells:
        area_vector_sum = np.zeros(3)
        for face in cell["cartesianFaces"]:
            normal = np.asarray(face["outwardNormal"], dtype=np.float64)
            area = float(face["area"])
            area_vector_sum += normal * area
            signed_loop_area = 0.0
            for loop in face["orientedBoundaryLoops"]:
                signed_loop_area += float(np.dot(polygon_area_vector(np.asarray(loop)), normal))
            if not math.isclose(signed_loop_area, area, abs_tol=tolerance):
                polygon_area_mismatches += 1
        for face in cell["embeddedBoundaryFaces"]:
            normal = np.asarray(face["outwardNormal"], dtype=np.float64)
            area = float(face["area"])
            vertices = np.asarray(face["vertices"], dtype=np.float64)
            vector = polygon_area_vector(vertices)
            if not math.isclose(float(np.linalg.norm(vector)), area, abs_tol=tolerance):
                polygon_area_mismatches += 1
            if float(np.dot(vector, normal)) < -tolerance:
                polygon_area_mismatches += 1
            area_vector_sum += normal * area
        closure = float(np.linalg.norm(area_vector_sum))
        maximum_recomputed_closure = max(maximum_recomputed_closure, closure)
        if closure > tolerance:
            failures.append(f"单元 {cell['backgroundCellId']} 面积向量不闭合")
            break
        edge_imbalance = boundary_edge_imbalance(cell, tolerance)
        maximum_boundary_edge_imbalance = max(maximum_boundary_edge_imbalance, edge_imbalance)
        if edge_imbalance != int(cell.get("boundaryEdgeImbalanceCount", edge_imbalance)):
            failures.append(f"单元 {cell['backgroundCellId']} 边链计数与报告不一致")
            break
        if edge_imbalance:
            boundary_edge_failure_cells += 1
    if polygon_area_mismatches:
        failures.append(f"显式多边形面积或定向不匹配数={polygon_area_mismatches}")
    if boundary_edge_failure_cells:
        failures.append(f"边界有向边链不闭合单元数={boundary_edge_failure_cells}")

    polyhedra = meshio.read(arguments.polyhedra)
    polyhedron_piece_count = sum(len(block.data) for block in polyhedra.cells)
    polyhedron_volume_mismatches = 0
    polyhedron_edge_mismatches = 0
    polyhedron_volume_sum = 0.0
    piece_volumes = flatten_cell_data(polyhedra, "piece_volume").astype(np.float64)
    piece_background_ids = flatten_cell_data(polyhedra, "background_cell_id").astype(np.uint64)
    piece_cursor = 0
    per_background_piece_volume: dict[int, float] = {}
    for block in polyhedra.cells:
        if not block.type.startswith("polyhedron"):
            failures.append(f"显式流体分解包含非 polyhedron 单元类型 {block.type}")
            continue
        for faces in block.data:
            unique_ids = sorted({int(point) for face in faces for point in face})
            reference = np.mean(polyhedra.points[unique_ids], axis=0)
            signed_volume = 0.0
            edge_balances: dict[tuple[int, int], int] = {}
            for face in faces:
                face = np.asarray(face, dtype=np.int64)
                first = polyhedra.points[face[0]]
                for index in range(1, len(face) - 1):
                    second = polyhedra.points[face[index]]
                    third = polyhedra.points[face[index + 1]]
                    signed_volume += float(
                        np.dot(first - reference, np.cross(second - reference, third - reference))
                    ) / 6.0
                for index, a in enumerate(face):
                    b = int(face[(index + 1) % len(face)])
                    a = int(a)
                    key = (min(a, b), max(a, b))
                    edge_balances[key] = edge_balances.get(key, 0) + (1 if a < b else -1)
            expected_piece_volume = float(piece_volumes[piece_cursor])
            if signed_volume <= 0.0 or not math.isclose(
                signed_volume, expected_piece_volume, abs_tol=tolerance
            ):
                polyhedron_volume_mismatches += 1
            if any(balance != 0 for balance in edge_balances.values()):
                polyhedron_edge_mismatches += 1
            polyhedron_volume_sum += signed_volume
            background_id = int(piece_background_ids[piece_cursor])
            per_background_piece_volume[background_id] = (
                per_background_piece_volume.get(background_id, 0.0) + signed_volume
            )
            piece_cursor += 1
    for background_id, piece_volume in per_background_piece_volume.items():
        if not math.isclose(piece_volume, float(cell_by_id[background_id]["volume"]), abs_tol=tolerance):
            polyhedron_volume_mismatches += 1
    expected_piece_count = sum(len(cell.get("fluidPolyhedronPieces", [])) for cell in cells)
    if polyhedron_piece_count != expected_piece_count:
        failures.append("polyhedron VTU 单元数与几何 JSON 分解片数不一致")
    if polyhedron_volume_mismatches:
        failures.append(f"显式 polyhedron 体积不匹配数={polyhedron_volume_mismatches}")
    if polyhedron_edge_mismatches:
        failures.append(f"显式 polyhedron 边闭合不匹配数={polyhedron_edge_mismatches}")

    topology_mismatches = 0
    seen_connections: set[tuple[int, int]] = set()
    for connection in geometry["internalFaces"]:
        first_id = int(connection["firstBackgroundCellId"])
        second_id = int(connection["secondBackgroundCellId"])
        first = cells[int(connection["firstFluidCellIndex"])]
        second = cells[int(connection["secondFluidCellIndex"])]
        first_face = first["cartesianFaces"][int(connection["firstLocalFace"])]
        second_face = second["cartesianFaces"][int(connection["secondLocalFace"])]
        if first["backgroundCellId"] != first_id or second["backgroundCellId"] != second_id:
            topology_mismatches += 1
        if not math.isclose(first_face["area"], second_face["area"], abs_tol=tolerance):
            topology_mismatches += 1
        if not np.allclose(first_face["centroid"], second_face["centroid"], atol=tolerance):
            topology_mismatches += 1
        if not np.allclose(first_face["outwardNormal"], -np.asarray(second_face["outwardNormal"]), atol=tolerance):
            topology_mismatches += 1
        seen_connections.add((first_id, second_id))
    nx, ny, nz = (int(value) for value in dimensions)
    strides = (1, nx, nx * ny)
    positive_faces = (1, 3, 5)
    negative_faces = (0, 2, 4)
    for cell_id, cell in cell_by_id.items():
        i = cell_id % nx
        j = (cell_id // nx) % ny
        k = cell_id // (nx * ny)
        coordinates = (i, j, k)
        limits = (nx, ny, nz)
        for axis in range(3):
            if coordinates[axis] + 1 >= limits[axis]:
                continue
            neighbor_id = cell_id + strides[axis]
            area = float(cell["cartesianFaces"][positive_faces[axis]]["area"])
            neighbor = cell_by_id.get(neighbor_id)
            neighbor_area = 0.0 if neighbor is None else float(
                neighbor["cartesianFaces"][negative_faces[axis]]["area"]
            )
            if not math.isclose(area, neighbor_area, abs_tol=tolerance):
                topology_mismatches += 1
            if area > tolerance and neighbor is not None and (cell_id, neighbor_id) not in seen_connections:
                topology_mismatches += 1
    if topology_mismatches:
        failures.append(f"cell-face-neighbor 拓扑不匹配数={topology_mismatches}")

    points, polygons, boundary_data = read_vtp(arguments.boundary)
    required_boundary_fields = {"background_cell_id", "boundary_id", "area", "fluid_outward_normal"}
    if not required_boundary_fields.issubset(boundary_data):
        failures.append("VTP 缺少嵌入边界字段")
    vtp_area = 0.0
    vtp_orientation_mismatches = 0
    for index, polygon in enumerate(polygons):
        vector = polygon_area_vector(points[polygon])
        area = float(np.linalg.norm(vector))
        vtp_area += area
        if "area" in boundary_data and not math.isclose(area, float(boundary_data["area"][index]), abs_tol=tolerance):
            vtp_orientation_mismatches += 1
        if "fluid_outward_normal" in boundary_data and float(
            np.dot(vector, boundary_data["fluid_outward_normal"][index])
        ) < -tolerance:
            vtp_orientation_mismatches += 1
    if not math.isclose(vtp_area, report["totalEmbeddedBoundaryArea"], abs_tol=tolerance):
        failures.append("VTP 多边形面积和与报告不一致")
    if vtp_orientation_mismatches:
        failures.append(f"VTP 面积/法向不匹配数={vtp_orientation_mismatches}")

    analytic: dict[str, float | int] = {}
    if arguments.shape == "cube":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            surface = meshio.read(arguments.surface)
        solid_minimum = np.min(surface.points, axis=0)
        solid_maximum = np.max(surface.points, axis=0)
        exact_solid_volume = float(np.prod(solid_maximum - solid_minimum))
        exact_fluid_volume = float(np.prod(domain_maximum - domain_minimum)) - exact_solid_volume
        ext = solid_maximum - solid_minimum
        exact_area = float(2.0 * (ext[0] * ext[1] + ext[0] * ext[2] + ext[1] * ext[2]))
        analytic_volume_mismatches = 0
        analytic_centroid_mismatches = 0
        expected_cut = 0
        expected_solid = 0
        expected_fluid = 0
        for cell_id in range(background_count):
            fluid_volume, fluid_centroid, _ = expected_cube_cell(
                cell_id, dimensions, domain_minimum, spacing, solid_minimum, solid_maximum
            )
            if fluid_volume <= tolerance:
                expected_solid += 1
                if cell_id in cell_by_id:
                    analytic_volume_mismatches += 1
                continue
            if math.isclose(fluid_volume, cell_volume, abs_tol=tolerance):
                expected_fluid += 1
            else:
                expected_cut += 1
            actual = cell_by_id.get(cell_id)
            if actual is None or not math.isclose(actual["volume"], fluid_volume, abs_tol=tolerance):
                analytic_volume_mismatches += 1
            elif not np.allclose(actual["centroid"], fluid_centroid, atol=tolerance):
                analytic_centroid_mismatches += 1
        if not math.isclose(report["totalFluidVolume"], exact_fluid_volume, abs_tol=tolerance):
            failures.append("立方体解析总流体体积不匹配")
        if not math.isclose(report["totalEmbeddedBoundaryArea"], exact_area, abs_tol=tolerance):
            failures.append("立方体解析嵌入边界面积不匹配")
        if (expected_cut != report["cutCellCount"] or
                expected_solid != report["fullSolidCellCount"] or
                expected_fluid != report["fullFluidCellCount"]):
            failures.append("立方体解析单元分类计数不匹配")
        if analytic_volume_mismatches or analytic_centroid_mismatches:
            failures.append(
                f"立方体逐单元真值不匹配：volume={analytic_volume_mismatches}, "
                f"centroid={analytic_centroid_mismatches}"
            )
        analytic = {
            "exactSolidVolume": exact_solid_volume,
            "exactFluidVolume": exact_fluid_volume,
            "exactEmbeddedBoundaryArea": exact_area,
            "expectedCutCellCount": expected_cut,
            "expectedFullSolidCellCount": expected_solid,
            "expectedFullFluidCellCount": expected_fluid,
            "cellVolumeMismatchCount": analytic_volume_mismatches,
            "cellCentroidMismatchCount": analytic_centroid_mismatches,
        }
    elif arguments.shape == "l_prism":
        solid_boxes = [
            (np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])),
            (np.array([0.0, 1.0, 0.0]), np.array([1.0, 2.0, 1.0])),
            (np.array([1.0, 0.0, 0.0]), np.array([2.0, 1.0, 1.0])),
        ]
        expected_cut = expected_solid = expected_fluid = 0
        volume_mismatches = centroid_mismatches = 0
        nx, ny, _ = (int(value) for value in dimensions)
        for cell_id in range(background_count):
            k, remainder = divmod(cell_id, nx * ny)
            j, i = divmod(remainder, nx)
            cell_minimum = domain_minimum + spacing * np.array([i, j, k])
            cell_maximum = cell_minimum + spacing
            solid_volume = 0.0
            solid_first_moment = np.zeros(3)
            for box_minimum, box_maximum in solid_boxes:
                overlap_minimum = np.maximum(cell_minimum, box_minimum)
                overlap_maximum = np.minimum(cell_maximum, box_maximum)
                overlap_extent = np.maximum(overlap_maximum - overlap_minimum, 0.0)
                volume = float(np.prod(overlap_extent))
                solid_volume += volume
                if volume > 0.0:
                    solid_first_moment += volume * 0.5 * (overlap_minimum + overlap_maximum)
            fluid_volume = cell_volume - solid_volume
            cell_centroid = 0.5 * (cell_minimum + cell_maximum)
            fluid_centroid = cell_centroid
            if fluid_volume > tolerance and solid_volume > tolerance:
                fluid_centroid = (
                    cell_volume * cell_centroid - solid_first_moment
                ) / fluid_volume
            actual = cell_by_id.get(cell_id)
            if fluid_volume <= tolerance:
                expected_solid += 1
                if actual is not None:
                    volume_mismatches += 1
                continue
            if math.isclose(fluid_volume, cell_volume, abs_tol=tolerance):
                expected_fluid += 1
            else:
                expected_cut += 1
            if actual is None or not math.isclose(actual["volume"], fluid_volume, abs_tol=tolerance):
                volume_mismatches += 1
            elif not np.allclose(actual["centroid"], fluid_centroid, atol=tolerance):
                centroid_mismatches += 1
        exact_solid_volume = 3.0
        exact_fluid_volume = float(np.prod(domain_maximum - domain_minimum)) - 3.0
        exact_area = 14.0
        if not math.isclose(report["totalSolidVolume"], exact_solid_volume, abs_tol=tolerance):
            failures.append("L 棱柱解析总固体体积不匹配")
        if not math.isclose(report["totalFluidVolume"], exact_fluid_volume, abs_tol=tolerance):
            failures.append("L 棱柱解析总流体体积不匹配")
        if not math.isclose(report["totalEmbeddedBoundaryArea"], exact_area, abs_tol=tolerance):
            failures.append("L 棱柱解析表面积不匹配")
        if (expected_cut != report["cutCellCount"] or
                expected_solid != report["fullSolidCellCount"] or
                expected_fluid != report["fullFluidCellCount"]):
            failures.append("L 棱柱解析单元计数不匹配")
        if volume_mismatches or centroid_mismatches:
            failures.append(
                f"L 棱柱逐单元真值不匹配：volume={volume_mismatches}, "
                f"centroid={centroid_mismatches}"
            )
        analytic = {
            "exactSolidVolume": exact_solid_volume,
            "exactFluidVolume": exact_fluid_volume,
            "exactEmbeddedBoundaryArea": exact_area,
            "expectedCutCellCount": expected_cut,
            "expectedFullSolidCellCount": expected_solid,
            "expectedFullFluidCellCount": expected_fluid,
            "cellVolumeMismatchCount": volume_mismatches,
            "cellCentroidMismatchCount": centroid_mismatches,
        }

    result = {
        "schema": "cartmesh-stage3-independent-verification-v1",
        "status": "pass" if not failures else "fail",
        "geometryTopologyValidated": not failures,
        "solverReadyCutCellMesh": False,
        "externalCfdCheckerAccepted": False,
        "meshioVersion": meshio.__version__,
        "backgroundCellCount": vtk_cell_count,
        "geometryFluidCellCount": len(cells),
        "embeddedPolygonCount": len(polygons),
        "fluidPolyhedronPieceCount": polyhedron_piece_count,
        "fluidPolyhedronVolumeSum": polyhedron_volume_sum,
        "polyhedronVolumeMismatchCount": polyhedron_volume_mismatches,
        "polyhedronEdgeMismatchCount": polyhedron_edge_mismatches,
        "vtuIntegratedFluidVolume": vtu_fluid_volume,
        "smallCellThreshold": small_threshold,
        "minimumCutCellVolumeFraction": expected_minimum_fraction,
        "smallCutCellCount": int(np.count_nonzero(expected_small)),
        "vtpIntegratedArea": vtp_area,
        "maximumRecomputedAreaClosureResidual": maximum_recomputed_closure,
        "boundaryEdgeFailureCellCount": boundary_edge_failure_cells,
        "maximumBoundaryEdgeImbalanceCount": maximum_boundary_edge_imbalance,
        "polygonAreaMismatchCount": polygon_area_mismatches,
        "topologyMismatchCount": topology_mismatches,
        "vtpOrientationMismatchCount": vtp_orientation_mismatches,
        "analytic": analytic,
        "failures": failures,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
