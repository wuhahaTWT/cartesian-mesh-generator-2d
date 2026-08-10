#!/usr/bin/env python3
"""独立验证均匀 Cartesian ASCII VTU 的文件结构和六面体拓扑。"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import xml.etree.ElementTree as ET


def numbers(element: ET.Element, converter: type[int] | type[float]) -> list[int] | list[float]:
    return [converter(token) for token in (element.text or "").split()]


def named_arrays(parent: ET.Element) -> dict[str, ET.Element]:
    return {array.attrib["Name"]: array for array in parent.findall("DataArray")}


def validate(path: pathlib.Path) -> dict[str, object]:
    root = ET.parse(path).getroot()
    if root.tag != "VTKFile" or root.attrib.get("type") != "UnstructuredGrid":
        raise ValueError("根节点必须是 type=UnstructuredGrid 的 VTKFile")
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("缺少 UnstructuredGrid/Piece 节点")

    point_count = int(piece.attrib["NumberOfPoints"])
    cell_count = int(piece.attrib["NumberOfCells"])
    point_array = piece.find("./Points/DataArray")
    if point_array is None:
        raise ValueError("缺少点坐标数组")
    coordinates = numbers(point_array, float)
    if len(coordinates) != point_count * 3:
        raise ValueError(f"点元组数量不一致：期望 {point_count}，实际 {len(coordinates) // 3}")
    if not all(math.isfinite(value) for value in coordinates):
        raise ValueError("点坐标数组包含非有限数值")

    cells = piece.find("Cells")
    if cells is None:
        raise ValueError("缺少 Cells 节点")
    arrays = named_arrays(cells)
    required = {"connectivity", "offsets", "types"}
    if set(arrays) != required:
        raise ValueError(f"单元数组为 {sorted(arrays)}，期望 {sorted(required)}")
    connectivity = numbers(arrays["connectivity"], int)
    offsets = numbers(arrays["offsets"], int)
    cell_types = numbers(arrays["types"], int)
    if len(connectivity) != cell_count * 8:
        raise ValueError("每个六面体必须正好具有 8 个连接索引")
    if offsets != list(range(8, cell_count * 8 + 1, 8)):
        raise ValueError("偏移数组不符合连续六面体的预期格式")
    if cell_types != [12] * cell_count:
        raise ValueError("所有单元都必须是类型 12：VTK_HEXAHEDRON")
    if connectivity and (min(connectivity) < 0 or max(connectivity) >= point_count):
        raise ValueError("连接数组引用了声明范围之外的点")

    sample_ids = sorted({0, cell_count // 2, cell_count - 1}) if cell_count else []
    for cell_id in sample_ids:
        node_ids = connectivity[cell_id * 8 : (cell_id + 1) * 8]
        if len(set(node_ids)) != 8:
            raise ValueError(f"单元 {cell_id} 包含重复顶点")
        points = [coordinates[node * 3 : node * 3 + 3] for node in node_ids]
        edge_x = points[1][0] - points[0][0]
        edge_y = points[3][1] - points[0][1]
        edge_z = points[4][2] - points[0][2]
        if not (edge_x > 0.0 and edge_y > 0.0 and edge_z > 0.0):
            raise ValueError(f"单元 {cell_id} 的 Cartesian 边长不是正数")

    fields: dict[str, int] = {}
    cell_data = piece.find("CellData")
    if cell_data is not None:
        for array in cell_data.findall("DataArray"):
            name = array.attrib.get("Name")
            if not name:
                raise ValueError("单元数据数组缺少名称")
            values = numbers(array, float)
            if len(values) != cell_count:
                raise ValueError(f"单元数据字段 {name!r} 的长度与单元数不一致")
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"单元数据字段 {name!r} 包含非有限数值")
            fields[name] = len(values)

    return {
        "path": str(path),
        "pointCount": point_count,
        "cellCount": cell_count,
        "cellType": "VTK_HEXAHEDRON",
        "cellDataFields": fields,
        "status": "pass",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="独立检查均匀 Cartesian VTU 文件")
    parser.add_argument("path", type=pathlib.Path)
    arguments = parser.parse_args()
    try:
        print(json.dumps(validate(arguments.path), ensure_ascii=False, sort_keys=True))
        return 0
    except (ET.ParseError, OSError, ValueError) as error:
        print(json.dumps({"path": str(arguments.path), "status": "fail", "error": str(error)}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
