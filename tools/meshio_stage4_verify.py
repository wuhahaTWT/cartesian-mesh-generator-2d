#!/usr/bin/env python3
"""使用 meshio/NumPy 独立检查阶段四多区域 Cut-cell 输出。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import meshio
import numpy as np


def fail(message: str) -> None:
    raise RuntimeError(message)


def values(mesh: meshio.Mesh, name: str) -> np.ndarray:
    blocks = mesh.cell_data.get(name)
    if blocks is None:
        fail(f"缺少 cell_data: {name}")
    return np.concatenate([np.asarray(block).reshape(-1) for block in blocks])


def vtp_values(path: Path, name: str) -> np.ndarray:
    root = ET.parse(path).getroot()
    for array in root.findall(".//CellData/DataArray"):
        if array.attrib.get("Name") == name:
            return np.fromstring(array.text or "", sep=" ")
    fail(f"VTP 缺少 CellData: {name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--boundary", type=Path, required=True)
    parser.add_argument("--polyhedra", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shape", choices=("thin-shell", "generic"),
                        default="generic")
    parser.add_argument("--expected-region-count", type=int)
    args = parser.parse_args()

    report = json.loads(args.report.read_text(encoding="utf-8"))
    geometry = json.loads(args.geometry.read_text(encoding="utf-8"))
    background = meshio.read(args.background)
    polyhedra = meshio.read(args.polyhedra)

    if report.get("schema") != "cartmesh-stage4-cutcell-v1":
        fail("项目报告不是阶段四机器契约")
    if report.get("stage4GeometryRobustnessComplete") is not True:
        fail("阶段四几何鲁棒性不变量未通过")
    if report.get("stage4Complete") is not False:
        fail("外部 CFD checker 尚未通过时不得标记 stage4Complete=true")
    if (report.get("status") != "geometry_pass_external_cfd_pending" or
            report.get("solverReadyCutCellMesh") is not False or
            report.get("externalCfdCheckerAccepted") is not False):
        fail("阶段四报告没有如实标记外部 CFD checker 门禁")
    if report.get("geometryRepairApplied") is not False:
        fail("默认路径不得静默修改输入几何")
    expected_actions = (["explicit_geometric_length_tolerance"]
                        if report.get("geometricTolerance", 0.0) > 0.0 else [])
    if report.get("toleranceActions") != expected_actions:
        fail("容差处理未在 toleranceActions 中逐项记录")

    region_count = int(report["globalFluidRegionCount"])
    region_volumes = np.array(
        [entry["volume"] for entry in report["globalFluidRegions"]], dtype=float)
    if region_count != len(region_volumes) or region_count == 0:
        fail("全局流体 region 计数与列表不一致")
    if (args.expected_region_count is not None and
            region_count != args.expected_region_count):
        fail("全局流体 region 计数与案例真值不一致")
    if not np.isclose(region_volumes.sum(), report["totalFluidVolume"],
                      rtol=2e-12, atol=2e-12):
        fail("全局 region 体积之和不等于总流体体积")

    geometry_regions = geometry.get("globalFluidRegions", [])
    if [entry["regionId"] for entry in geometry_regions] != list(range(region_count)):
        fail("完整几何 JSON 的 region ID 不连续")
    if not np.allclose([entry["volume"] for entry in geometry_regions],
                       region_volumes, rtol=2e-12, atol=2e-12):
        fail("报告与完整几何 JSON 的 region 体积不一致")

    fluid_cells = geometry["fluidCells"]
    for cell in fluid_cells:
        local_regions = cell["fluidComponentRegionIds"]
        if len(local_regions) != cell["fluidComponentCount"]:
            fail("单元局部分量数与 region ID 数不一致")
        if any(region < 0 or region >= region_count for region in local_regions):
            fail("单元 region ID 越界")
        for piece in cell["fluidPolyhedronPieces"]:
            if piece["globalRegionId"] != local_regions[piece["componentId"]]:
                fail("显式多面体片的 region ID 与局部分量不一致")
    for face in geometry["componentInternalFaces"]:
        first = fluid_cells[face["firstFluidCellIndex"]]
        second = fluid_cells[face["secondFluidCellIndex"]]
        region = face["globalRegionId"]
        if (first["fluidComponentRegionIds"][face["firstComponentId"]] != region or
                second["fluidComponentRegionIds"][face["secondComponentId"]] != region):
            fail("分量邻接面跨越了不同全局 region")

    background_regions = values(background, "fluid_region_id")
    if np.any((background_regions < -1) | (background_regions >= region_count)):
        fail("背景 VTU 的 fluid_region_id 越界")
    multiple = values(background, "multiple_fluid_components")
    if not np.all(np.isin(multiple, (0.0, 1.0))):
        fail("multiple_fluid_components 必须是 0/1")
    small_field = values(background, "small_cut_cell")
    small_entries = report["smallCutCells"]
    if int(np.count_nonzero(small_field == 1.0)) != report["smallCutCellCount"]:
        fail("背景 VTU 小 Cut-cell 标记数与报告不一致")
    if len(small_entries) != report["smallCutCellCount"]:
        fail("小 Cut-cell 位置列表与计数不一致")
    for entry in small_entries:
        if (entry["volumeFraction"] >= report["smallCellThreshold"] or
                len(entry["centroid"]) != 3 or
                entry["backgroundCellId"] < 0):
            fail("小 Cut-cell 条目缺少阈值、位置或背景 ID 证据")

    boundary_ids = vtp_values(args.boundary, "boundary_id").astype(np.uint64)
    boundary_areas = vtp_values(args.boundary, "area")
    if not np.isclose(boundary_areas.sum(), report["totalEmbeddedBoundaryArea"],
                      rtol=2e-12, atol=2e-12):
        fail("边界 VTP 面积与报告不一致")
    named_ids = {int(item["boundaryId"]) for item in report["boundaries"]}
    if named_ids and set(map(int, np.unique(boundary_ids))) != named_ids:
        fail("边界 VTP 未完整保留已命名 boundary ID")

    poly_regions = values(polyhedra, "global_region_id")
    if np.any((poly_regions < 0) | (poly_regions >= region_count)):
        fail("显式多面体 VTU 的 global_region_id 越界")
    if len(values(polyhedra, "component_id")) != len(poly_regions):
        fail("多面体 component/region 数据长度不一致")

    analytic = {}
    if args.shape == "thin-shell":
        expected_regions = np.array([0.9 ** 3, 1.2 ** 3 - 1.0])
        expected_regions.sort()
        if region_count != 2 or not np.allclose(
                np.sort(region_volumes), expected_regions,
                rtol=2e-11, atol=2e-11):
            fail("薄壁内腔/外部 region 体积不符解析真值")
        expected_solid = 1.0 - 0.9 ** 3
        expected_area = 6.0 + 6.0 * 0.9 ** 2
        if {entry.get("name") for entry in report["globalFluidRegions"]} != {
                "exterior", "cavity"}:
            fail("薄壁案例未保留 exterior/cavity 流体区命名")
        if {entry.get("name") for entry in report["boundaries"]} != {
                "outer_wall", "cavity_wall"}:
            fail("薄壁案例未保留内外壁 boundary 命名")
        if not np.isclose(report["totalSolidVolume"], expected_solid,
                          rtol=2e-11, atol=2e-11):
            fail("薄壁材料体积不符解析真值")
        if not np.isclose(report["totalEmbeddedBoundaryArea"], expected_area,
                          rtol=2e-11, atol=2e-11):
            fail("薄壁内外表面积不符解析真值")
        analytic = {"expectedSolidVolume": expected_solid,
                    "expectedBoundaryArea": expected_area,
                    "expectedRegionVolumes": expected_regions.tolist()}

    result = {
        "status": "pass",
        "geometryRobustnessValidated": True,
        "solverReadyCutCellMesh": False,
        "externalCfdCheckerAccepted": False,
        "meshioVersion": meshio.__version__,
        "shape": args.shape,
        "globalFluidRegionCount": region_count,
        "globalFluidRegionVolumes": region_volumes.tolist(),
        "globalFluidRegionNames": [entry.get("name", "")
                                   for entry in report["globalFluidRegions"]],
        "componentInternalFaceCount": len(geometry["componentInternalFaces"]),
        "backgroundCellCount": int(sum(len(block.data) for block in background.cells)),
        "embeddedBoundaryPolygonCount": int(len(boundary_ids)),
        "explicitFluidPolyhedronCount": int(len(poly_regions)),
        "smallCutCellCount": report["smallCutCellCount"],
        "boundaryIds": sorted(map(int, np.unique(boundary_ids))),
        "discardedNumericalPieceCount": report["discardedNumericalPieceCount"],
        "discardedNumericalPieceVolume": report["discardedNumericalPieceVolume"],
        "toleranceActions": report["toleranceActions"],
        "analytic": analytic,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
