#!/usr/bin/env python3
"""用 meshio/Matplotlib 独立渲染相交单元外表面和分类切片。"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import pathlib
import tempfile
import warnings
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR", str(pathlib.Path(tempfile.gettempdir()) / "cartmesh-matplotlib-cache")
)
pathlib.Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Patch, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import meshio
import numpy as np


HEX_FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)
COLORS = ("#e1e7ef", "#2fb56f", "#1494df", "#df3038")
LABELS = ("outside", "inside", "intersected", "conflict")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def one_block(mesh: meshio.Mesh, cell_type: str) -> np.ndarray:
    blocks = [block.data for block in mesh.cells if block.type == cell_type]
    if len(blocks) != 1 or len(mesh.cells) != 1:
        raise ValueError(f"期望唯一 {cell_type} 块，实际为 {[block.type for block in mesh.cells]}")
    return np.asarray(blocks[0], dtype=np.int64)


def read_mesh(path: pathlib.Path) -> meshio.Mesh:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="overflow encountered in scalar multiply",
            category=RuntimeWarning,
            module=r"meshio\.stl\._stl",
        )
        return meshio.read(path)


def classification(mesh: meshio.Mesh, cell_count: int) -> np.ndarray:
    arrays = mesh.cell_data.get("stl_cell_classification")
    if arrays is None or len(arrays) != 1:
        raise ValueError("VTU 缺少 stl_cell_classification")
    values = np.asarray(arrays[0])
    rounded = np.rint(values).astype(np.int64)
    if (
        len(values) != cell_count
        or not np.array_equal(values, rounded)
        or not np.isin(rounded, (0, 1, 2, 3)).all()
    ):
        raise ValueError("stl_cell_classification 长度或图例值无效")
    return rounded


def selected_boundary_faces(cells: np.ndarray, selected: np.ndarray) -> list[np.ndarray]:
    uses: dict[tuple[int, int, int, int], tuple[int, np.ndarray]] = {}
    for cell_id in selected:
        cell = cells[int(cell_id)]
        for local_face in HEX_FACES:
            face = cell[np.asarray(local_face, dtype=np.int64)]
            key = tuple(sorted(int(value) for value in face))
            previous = uses.get(key)
            if previous is None:
                uses[key] = (1, face)
            else:
                uses[key] = (previous[0] + 1, previous[1])
    return [face for count, face in uses.values() if count == 1]


def set_equal_3d_axes(axis: Any, points: np.ndarray) -> None:
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = 0.5 * (minimum + maximum)
    radius = 0.52 * float(np.max(maximum - minimum))
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1.0, 1.0, 1.0))


def render_overview(
    points: np.ndarray,
    cells: np.ndarray,
    values: np.ndarray,
    surface_points: np.ndarray,
    surface_triangles: np.ndarray,
    output: pathlib.Path,
) -> dict[str, Any]:
    selected = np.flatnonzero(values == 2)
    faces = selected_boundary_faces(cells, selected)
    figure = plt.figure(figsize=(10.5, 8.4), dpi=150)
    axis = figure.add_subplot(111, projection="3d")
    cell_collection = Poly3DCollection(
        [points[face] for face in faces],
        facecolor=COLORS[2],
        edgecolor="#0a3959",
        linewidth=0.12,
        alpha=0.58,
    )
    axis.add_collection3d(cell_collection)
    surface_collection = Poly3DCollection(
        surface_points[surface_triangles],
        facecolor="#d2d7de",
        edgecolor="#2b3139",
        linewidth=0.20,
        alpha=0.32,
    )
    axis.add_collection3d(surface_collection)
    set_equal_3d_axes(axis, points)
    axis.view_init(elev=24.0, azim=38.0)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.set_title("Stage 1 exact intersected Cartesian cells with STL overlay")
    axis.grid(False)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    return {
        "path": str(output),
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "selectedIntersectedCellCount": int(len(selected)),
        "selectedUnionBoundaryFaceCount": int(len(faces)),
    }


def triangle_plane_segments(
    points: np.ndarray, triangles: np.ndarray, z_value: float, tolerance: float
) -> list[np.ndarray]:
    segments: list[np.ndarray] = []
    for triangle in points[triangles]:
        signed = triangle[:, 2] - z_value
        candidates: list[np.ndarray] = []
        for first, second in ((0, 1), (1, 2), (2, 0)):
            a = triangle[first]
            b = triangle[second]
            da = signed[first]
            db = signed[second]
            if abs(da) <= tolerance:
                candidates.append(a[:2])
            if da * db < 0.0:
                fraction = da / (da - db)
                candidates.append((a + fraction * (b - a))[:2])
        unique: list[np.ndarray] = []
        for candidate in candidates:
            if not any(np.linalg.norm(candidate - value) <= tolerance for value in unique):
                unique.append(candidate)
        if len(unique) == 2:
            segments.append(np.asarray(unique))
        elif len(unique) >= 3 and np.all(np.abs(signed) <= tolerance):
            for index in range(3):
                segments.append(np.asarray((triangle[index, :2], triangle[(index + 1) % 3, :2])))
    return segments


def render_slice(
    points: np.ndarray,
    cells: np.ndarray,
    values: np.ndarray,
    surface_points: np.ndarray,
    surface_triangles: np.ndarray,
    z_value: float,
    output: pathlib.Path,
) -> dict[str, Any]:
    cell_points = points[cells]
    minimum = cell_points.min(axis=1)
    maximum = cell_points.max(axis=1)
    global_maximum_z = float(points[:, 2].max())
    in_slice = (minimum[:, 2] <= z_value) & (
        (z_value < maximum[:, 2])
        | (np.isclose(z_value, global_maximum_z) & np.isclose(maximum[:, 2], z_value))
    )
    selected = np.flatnonzero(in_slice)
    if not len(selected):
        raise ValueError("指定 z 平面没有穿过任何单元")
    patches = [
        Rectangle(
            (minimum[cell_id, 0], minimum[cell_id, 1]),
            maximum[cell_id, 0] - minimum[cell_id, 0],
            maximum[cell_id, 1] - minimum[cell_id, 1],
        )
        for cell_id in selected
    ]
    patch_collection = PatchCollection(
        patches,
        facecolor=[COLORS[int(values[cell_id])] for cell_id in selected],
        edgecolor="#55606d",
        linewidth=0.18,
    )
    scale = max(1.0, float(np.max(surface_points.max(axis=0) - surface_points.min(axis=0))))
    segments = triangle_plane_segments(
        surface_points, surface_triangles, z_value, 64.0 * np.finfo(float).eps * scale
    )
    figure, axis = plt.subplots(figsize=(9.2, 8.4), dpi=160)
    axis.add_collection(patch_collection)
    if segments:
        axis.add_collection(LineCollection(segments, colors="#111820", linewidths=1.35))
    axis.set_xlim(points[:, 0].min(), points[:, 0].max())
    axis.set_ylim(points[:, 1].min(), points[:, 1].max())
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(f"Stage 1 cell classification slice at z={z_value:.9g}")
    axis.legend(
        handles=[Patch(facecolor=COLORS[index], label=f"{index}: {LABELS[index]}") for index in range(4)],
        loc="upper right",
        framealpha=0.95,
    )
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    counts = collections.Counter(int(values[cell_id]) for cell_id in selected)
    return {
        "path": str(output),
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "z": z_value,
        "sliceCellCount": int(len(selected)),
        "sliceClassificationCounts": {
            LABELS[index]: int(counts[index]) for index in range(4)
        },
        "surfaceIntersectionSegmentCount": int(len(segments)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", required=True, type=pathlib.Path)
    parser.add_argument("--mesh", required=True, type=pathlib.Path)
    parser.add_argument("--overview", required=True, type=pathlib.Path)
    parser.add_argument("--slice", required=True, type=pathlib.Path)
    parser.add_argument("--report", required=True, type=pathlib.Path)
    parser.add_argument("--slice-z", required=True, type=float)
    arguments = parser.parse_args()
    try:
        volume = read_mesh(arguments.mesh)
        cells = one_block(volume, "hexahedron")
        points = np.asarray(volume.points, dtype=np.float64)
        values = classification(volume, len(cells))
        surface = read_mesh(arguments.surface)
        surface_triangles = one_block(surface, "triangle")
        surface_points = np.asarray(surface.points, dtype=np.float64)
        if not np.isfinite(points).all() or not np.isfinite(surface_points).all():
            raise ValueError("输入包含非有限坐标")
        counts = np.bincount(values, minlength=4)
        overview = render_overview(
            points,
            cells,
            values,
            surface_points,
            surface_triangles,
            arguments.overview,
        )
        slice_result = render_slice(
            points,
            cells,
            values,
            surface_points,
            surface_triangles,
            arguments.slice_z,
            arguments.slice,
        )
        result = {
            "schemaVersion": 1,
            "projectStage": 1,
            "reader": f"meshio {meshio.__version__}",
            "renderer": f"Matplotlib {matplotlib.__version__} Agg",
            "surfaceTriangleCount": int(len(surface_triangles)),
            "pointCount": int(len(points)),
            "cellCount": int(len(cells)),
            "classificationField": "stl_cell_classification",
            "classificationCounts": {
                LABELS[index]: int(counts[index]) for index in range(4)
            },
            "overview": overview,
            "slice": slice_result,
            "solverReadyCutCellMesh": False,
            "status": "pass",
        }
        serialized = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
        arguments.report.parent.mkdir(parents=True, exist_ok=True)
        arguments.report.write_text(serialized + "\n", encoding="utf-8")
        print(serialized)
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(json.dumps({"status": "fail", "error": str(error)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
