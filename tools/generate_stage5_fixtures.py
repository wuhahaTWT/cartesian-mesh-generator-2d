#!/usr/bin/env python3
"""生成阶段 5 的确定性体素构造解析 STL 对。

这里的体素只用于构造封闭的解析测试表面；被测输出仍是程序生成的
自适应 Cut-cell 网格，不能把输入构造方式称为 Cut-cell。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Iterable


Point = tuple[float, float, float]
Triangle = tuple[Point, Point, Point]


def _face_vertices(i: int, j: int, k: int, direction: tuple[int, int, int],
                   nx: int, ny: int, nz: int) -> tuple[Point, Point, Point, Point]:
    x0, x1 = i / nx, (i + 1) / nx
    y0, y1 = j / ny, (j + 1) / ny
    z0, z1 = k / nz, (k + 1) / nz
    if direction == (-1, 0, 0):
        return ((x0, y0, z0), (x0, y0, z1), (x0, y1, z1), (x0, y1, z0))
    if direction == (1, 0, 0):
        return ((x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1))
    if direction == (0, -1, 0):
        return ((x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1))
    if direction == (0, 1, 0):
        return ((x0, y1, z0), (x0, y1, z1), (x1, y1, z1), (x1, y1, z0))
    if direction == (0, 0, -1):
        return ((x0, y0, z0), (x0, y1, z0), (x1, y1, z0), (x1, y0, z0))
    if direction == (0, 0, 1):
        return ((x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1))
    raise ValueError(f"unsupported direction {direction}")


def voxel_surface(nx: int, ny: int, nz: int,
                  solid: Callable[[int, int, int], bool]) -> list[Triangle]:
    directions = ((-1, 0, 0), (1, 0, 0), (0, -1, 0),
                  (0, 1, 0), (0, 0, -1), (0, 0, 1))
    triangles: list[Triangle] = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if not solid(i, j, k):
                    continue
                for di, dj, dk in directions:
                    ni, nj, nk = i + di, j + dj, k + dk
                    neighbor_solid = (
                        0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz
                        and solid(ni, nj, nk)
                    )
                    if neighbor_solid:
                        continue
                    vertices = _face_vertices(i, j, k, (di, dj, dk), nx, ny, nz)
                    triangles.append((vertices[0], vertices[1], vertices[2]))
                    triangles.append((vertices[0], vertices[2], vertices[3]))
    return triangles


def normal(triangle: Triangle) -> Point:
    a, b, c = triangle
    u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    cross = (u[1] * v[2] - u[2] * v[1],
             u[2] * v[0] - u[0] * v[2],
             u[0] * v[1] - u[1] * v[0])
    length = sum(value * value for value in cross) ** 0.5
    return tuple(value / length for value in cross)  # type: ignore[return-value]


def write_ascii_stl(path: Path, name: str, triangles: Iterable[Triangle]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as output:
        output.write(f"solid {name}\n")
        for triangle in triangles:
            n = normal(triangle)
            output.write(f"  facet normal {n[0]:.17g} {n[1]:.17g} {n[2]:.17g}\n")
            output.write("    outer loop\n")
            for point in triangle:
                output.write(
                    f"      vertex {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}\n"
                )
            output.write("    endloop\n  endfacet\n")
        output.write(f"endsolid {name}\n")


def make_fixture(output: Path, name: str, removed_xy: set[tuple[int, int]]) -> None:
    nx, ny, nz = 10, 10, 6
    triangles = voxel_surface(
        nx, ny, nz,
        lambda i, j, _k: (i, j) not in removed_xy,
    )
    write_ascii_stl(output / f"{name}.stl", name, triangles)


def make_minimal_coplanar_regression(path: Path) -> None:
    nx, ny, nz = 4, 4, 2
    removed_xy = {(1, 1), (2, 1), (1, 2)}
    triangles = voxel_surface(
        nx, ny, nz,
        lambda i, j, _k: (i, j) not in removed_xy,
    )
    write_ascii_stl(path, "stage5_coplanar_tunnel_corner", triangles)


def generate(output: Path) -> None:
    small_center = {(4, 4), (4, 5), (5, 4), (5, 5)}
    large_center = {(i, j) for i in range(3, 7) for j in range(3, 7)}
    shifted = {(5, 4), (5, 5), (6, 4), (6, 5)}
    local_lobe = small_center | {(6, 5)}

    make_fixture(output, "hole_diameter_old", small_center)
    make_fixture(output, "hole_diameter_new", large_center)
    make_fixture(output, "hole_position_old", small_center)
    make_fixture(output, "hole_position_new", shifted)
    make_fixture(output, "local_contour_old", small_center)
    make_fixture(output, "local_contour_new", local_lobe)
    make_minimal_coplanar_regression(
        Path("tests/data/stage5_coplanar_tunnel_corner_ascii.stl")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/analytic/stage5"),
    )
    args = parser.parse_args()
    generate(args.output)


if __name__ == "__main__":
    main()
