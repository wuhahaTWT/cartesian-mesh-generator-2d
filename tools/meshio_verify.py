#!/usr/bin/env python3
"""使用外部 meshio 读取器检查阶段 0 VTU，并可通过 Matplotlib 渲染。"""

from __future__ import annotations

import argparse
import json
import pathlib

import meshio
import numpy as np


HEX_FACES = (
    (0, 1, 2, 3),
    (4, 7, 6, 5),
    (0, 4, 5, 1),
    (1, 5, 6, 2),
    (2, 6, 7, 3),
    (3, 7, 4, 0),
)


def render(
    points: np.ndarray,
    cells: np.ndarray,
    screenshot: pathlib.Path,
    title: str,
) -> int:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    boundary: dict[tuple[int, ...], tuple[int, ...]] = {}
    for cell in cells:
        for local_face in HEX_FACES:
            face = tuple(int(cell[index]) for index in local_face)
            key = tuple(sorted(face))
            if key in boundary:
                del boundary[key]
            else:
                boundary[key] = face
    faces = [points[list(face)] for face in boundary.values()]
    figure = plt.figure(figsize=(9, 8), dpi=140)
    axes = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(
        faces, facecolor="#3ba7d8", edgecolor="#17364a", linewidth=0.08, alpha=0.96
    )
    axes.add_collection3d(collection)
    used_points = points[np.unique(cells)]
    minimum = used_points.min(axis=0)
    maximum = used_points.max(axis=0)
    axes.set_xlim(minimum[0], maximum[0])
    axes.set_ylim(minimum[1], maximum[1])
    axes.set_zlim(minimum[2], maximum[2])
    axes.set_box_aspect(maximum - minimum)
    axes.view_init(elev=24, azim=38)
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_zlabel("z")
    axes.set_title(title)
    figure.tight_layout()
    screenshot.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(screenshot)
    plt.close(figure)
    return len(faces)


def main() -> int:
    parser = argparse.ArgumentParser(description="使用 meshio 独立读取并检查阶段 0 VTU")
    parser.add_argument("mesh", type=pathlib.Path)
    parser.add_argument("--field", help="保留指定单元字段值不小于 0.5 的单元")
    parser.add_argument("--screenshot", type=pathlib.Path)
    parser.add_argument("--title", default="阶段 0 单元中心采样 Cartesian 网格")
    arguments = parser.parse_args()

    mesh = meshio.read(arguments.mesh)
    hexahedra = [block for block in mesh.cells if block.type == "hexahedron"]
    if len(hexahedra) != 1 or len(mesh.cells) != 1:
        raise RuntimeError("必须正好存在一个六面体单元块")
    cells = hexahedra[0].data
    if cells.shape[1] != 8:
        raise RuntimeError("六面体连接数组的宽度不是 8")
    if cells.size and (cells.min() < 0 or cells.max() >= len(mesh.points)):
        raise RuntimeError("meshio 发现超出范围的连接索引")

    selected = cells
    if arguments.field:
        arrays = mesh.cell_data.get(arguments.field)
        if arrays is None or len(arrays) != 1:
            raise RuntimeError(f"缺少单块单元数据字段：{arguments.field}")
        if len(arrays[0]) != len(cells):
            raise RuntimeError("单元数据长度与 meshio 读取的单元数不同")
        selected = cells[np.asarray(arrays[0]) >= 0.5]
    boundary_faces = None
    if arguments.screenshot:
        boundary_faces = render(
            mesh.points, selected, arguments.screenshot, arguments.title
        )

    print(
        json.dumps(
            {
                "reader": f"meshio {meshio.__version__}",
                "pointCount": len(mesh.points),
                "cellCount": len(cells),
                "selectedCellCount": len(selected),
                "cellType": "hexahedron",
                "cellDataFields": sorted(mesh.cell_data),
                "renderedBoundaryFaces": boundary_faces,
                "screenshot": str(arguments.screenshot) if arguments.screenshot else None,
                "status": "pass",
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
