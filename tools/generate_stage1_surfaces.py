#!/usr/bin/env python3
"""生成阶段 1 的确定性封闭 ASCII STL 球体和圆柱解析基准。"""

from __future__ import annotations

import argparse
import math
import pathlib
from collections.abc import Iterable

Point = tuple[float, float, float]
Triangle = tuple[Point, Point, Point]


def subtract(lhs: Point, rhs: Point) -> Point:
    return (lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2])


def cross(lhs: Point, rhs: Point) -> Point:
    return (
        lhs[1] * rhs[2] - lhs[2] * rhs[1],
        lhs[2] * rhs[0] - lhs[0] * rhs[2],
        lhs[0] * rhs[1] - lhs[1] * rhs[0],
    )


def normal(triangle: Triangle) -> Point:
    first = subtract(triangle[1], triangle[0])
    second = subtract(triangle[2], triangle[0])
    value = cross(first, second)
    length = math.sqrt(sum(component * component for component in value))
    return tuple(component / length for component in value)  # type: ignore[return-value]


def sphere_triangles(latitudes: int, longitudes: int) -> list[Triangle]:
    top: Point = (0.0, 0.0, 1.0)
    bottom: Point = (0.0, 0.0, -1.0)
    rings: list[list[Point]] = []
    for latitude in range(1, latitudes):
        phi = math.pi * latitude / latitudes
        ring = []
        for longitude in range(longitudes):
            theta = 2.0 * math.pi * longitude / longitudes
            ring.append(
                (math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi))
            )
        rings.append(ring)
    triangles: list[Triangle] = []
    for longitude in range(longitudes):
        following = (longitude + 1) % longitudes
        triangles.append((top, rings[0][longitude], rings[0][following]))
        for ring in range(len(rings) - 1):
            triangles.append(
                (rings[ring][longitude], rings[ring + 1][longitude], rings[ring + 1][following])
            )
            triangles.append(
                (rings[ring][longitude], rings[ring + 1][following], rings[ring][following])
            )
        triangles.append((bottom, rings[-1][following], rings[-1][longitude]))
    return triangles


def cylinder_triangles(segments: int) -> list[Triangle]:
    bottom_center: Point = (0.0, 0.0, -1.0)
    top_center: Point = (0.0, 0.0, 1.0)
    bottom: list[Point] = []
    top: list[Point] = []
    for segment in range(segments):
        theta = 2.0 * math.pi * segment / segments
        bottom.append((math.cos(theta), math.sin(theta), -1.0))
        top.append((math.cos(theta), math.sin(theta), 1.0))
    triangles: list[Triangle] = []
    for segment in range(segments):
        following = (segment + 1) % segments
        triangles.extend(
            (
                (bottom_center, bottom[following], bottom[segment]),
                (top_center, top[segment], top[following]),
                (bottom[segment], bottom[following], top[following]),
                (bottom[segment], top[following], top[segment]),
            )
        )
    return triangles


def cube_triangles(shift: Point, scale: float, reverse: bool = False) -> list[Triangle]:
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    ]
    transformed = [
        tuple(shift[axis] + scale * point[axis] for axis in range(3))
        for point in points
    ]
    result: list[Triangle] = []
    for first, second, third in faces:
        result.append(
            (transformed[first], transformed[third], transformed[second])
            if reverse
            else (transformed[first], transformed[second], transformed[third])
        )
    return result


def sliver_tetrahedron_triangles() -> list[Triangle]:
    points: list[Point] = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0e-20, 0.0),
        (0.0, 0.0, 1.0),
    ]
    center = tuple(sum(point[axis] for point in points) / 4.0 for axis in range(3))
    result: list[Triangle] = []
    for first, second, third in ((0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)):
        triangle = (points[first], points[second], points[third])
        direction = subtract(triangle[0], center)  # type: ignore[arg-type]
        if sum(a * b for a, b in zip(cross(subtract(triangle[1], triangle[0]), subtract(triangle[2], triangle[0])), direction)) < 0.0:
            triangle = (points[first], points[third], points[second])
        result.append(triangle)
    return result


def tube_triangles(segments: int, inner_radius: float) -> list[Triangle]:
    outer_bottom: list[Point] = []
    outer_top: list[Point] = []
    inner_bottom: list[Point] = []
    inner_top: list[Point] = []
    for segment in range(segments):
        theta = 2.0 * math.pi * segment / segments
        cosine = math.cos(theta)
        sine = math.sin(theta)
        outer_bottom.append((cosine, sine, -1.0))
        outer_top.append((cosine, sine, 1.0))
        inner_bottom.append((inner_radius * cosine, inner_radius * sine, -1.0))
        inner_top.append((inner_radius * cosine, inner_radius * sine, 1.0))
    triangles: list[Triangle] = []
    for segment in range(segments):
        following = (segment + 1) % segments
        triangles.extend(
            (
                (outer_bottom[segment], outer_bottom[following], outer_top[following]),
                (outer_bottom[segment], outer_top[following], outer_top[segment]),
                (inner_bottom[segment], inner_top[following], inner_bottom[following]),
                (inner_bottom[segment], inner_top[segment], inner_top[following]),
                (outer_top[segment], outer_top[following], inner_top[following]),
                (outer_top[segment], inner_top[following], inner_top[segment]),
                (outer_bottom[segment], inner_bottom[following], outer_bottom[following]),
                (outer_bottom[segment], inner_bottom[segment], inner_bottom[following]),
            )
        )
    return triangles


def write_ascii_stl(path: pathlib.Path, name: str, triangles: Iterable[Triangle]) -> int:
    triangle_list = list(triangles)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as output:
        output.write(f"solid {name}\n")
        for triangle in triangle_list:
            surface_normal = normal(triangle)
            output.write(
                "facet normal " + " ".join(format(value, ".17g") for value in surface_normal) + "\n"
            )
            output.write("  outer loop\n")
            for vertex in triangle:
                output.write("    vertex " + " ".join(format(value, ".17g") for value in vertex) + "\n")
            output.write("  endloop\nendfacet\n")
        output.write(f"endsolid {name}\n")
    return len(triangle_list)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=pathlib.Path)
    arguments = parser.parse_args()
    sphere_path = arguments.output_directory / "sphere_lat24_lon48_ascii.stl"
    cylinder_path = arguments.output_directory / "cylinder_96_ascii.stl"
    tube_path = arguments.output_directory / "tube_outer96_inner05_ascii.stl"
    small_hole_tube_path = arguments.output_directory / "tube_outer96_inner008_ascii.stl"
    tiny_cube_path = arguments.output_directory / "cube_scale1e-7_ascii.stl"
    translated_cube_path = arguments.output_directory / "cube_shift1e9_ascii.stl"
    mixed_orientation_path = arguments.output_directory / "two_cubes_one_reversed_ascii.stl"
    sliver_path = arguments.output_directory / "sliver_tetrahedron_ascii.stl"
    sphere_count = write_ascii_stl(sphere_path, "sphere_lat24_lon48", sphere_triangles(24, 48))
    cylinder_count = write_ascii_stl(cylinder_path, "cylinder_96", cylinder_triangles(96))
    tube_count = write_ascii_stl(tube_path, "tube_outer96_inner05", tube_triangles(96, 0.5))
    small_hole_tube_count = write_ascii_stl(
        small_hole_tube_path, "tube_outer96_inner008", tube_triangles(96, 0.08)
    )
    tiny_cube_count = write_ascii_stl(
        tiny_cube_path, "cube_scale1e-7", cube_triangles((0.0, 0.0, 0.0), 1.0e-7)
    )
    translated_cube_count = write_ascii_stl(
        translated_cube_path, "cube_shift1e9", cube_triangles((1.0e9, 1.0e9, 1.0e9), 1.0)
    )
    mixed_orientation_count = write_ascii_stl(
        mixed_orientation_path,
        "two_cubes_one_reversed",
        cube_triangles((0.0, 0.0, 0.0), 1.0)
        + cube_triangles((2.0, 0.0, 0.0), 1.0, reverse=True),
    )
    sliver_count = write_ascii_stl(
        sliver_path, "sliver_tetrahedron", sliver_tetrahedron_triangles()
    )
    print(f"{sphere_path} triangles={sphere_count}")
    print(f"{cylinder_path} triangles={cylinder_count}")
    print(f"{tube_path} triangles={tube_count}")
    print(f"{small_hole_tube_path} triangles={small_hole_tube_count}")
    print(f"{tiny_cube_path} triangles={tiny_cube_count}")
    print(f"{translated_cube_path} triangles={translated_cube_count}")
    print(f"{mixed_orientation_path} triangles={mixed_orientation_count}")
    print(f"{sliver_path} triangles={sliver_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
