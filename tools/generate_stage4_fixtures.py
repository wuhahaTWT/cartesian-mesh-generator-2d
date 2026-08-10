#!/usr/bin/env python3
"""生成阶段四的确定性 ASCII STL 验收几何。"""

from __future__ import annotations

import argparse
from pathlib import Path


BOX_FACES = (
    (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
    (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
)


def box(minimum: tuple[float, float, float],
        maximum: tuple[float, float, float], reverse: bool = False):
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    vertices = (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )
    for face in BOX_FACES:
        indices = (face[0], face[2], face[1]) if reverse else face
        yield tuple(vertices[index] for index in indices)


def write_ascii_stl(path: Path, name: str, triangles) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"solid {name}"]
    for triangle in triangles:
        lines.extend((
            "  facet normal 0 0 0",
            "    outer loop",
            *(f"      vertex {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}"
              for point in triangle),
            "    endloop",
            "  endfacet",
        ))
    lines.append(f"endsolid {name}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-directory", type=Path,
                        default=Path("benchmarks/analytic/stage4"))
    args = parser.parse_args()
    output = args.output_directory

    thin_shell = list(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
    thin_shell.extend(box((0.05, 0.05, 0.05), (0.95, 0.95, 0.95), True))
    write_ascii_stl(output / "thin_shell_wall005_ascii.stl",
                    "thin_shell_wall005", thin_shell)

    two_cubes = list(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
    two_cubes.extend(box((2.0, 0.0, 0.0), (3.0, 1.0, 1.0)))
    write_ascii_stl(output / "two_disjoint_cubes_ascii.stl",
                    "two_disjoint_cubes", two_cubes)

    overlap = [
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ((0.2, 0.2, 0.0), (0.8, 0.2, 0.0), (0.2, 0.8, 0.0)),
    ]
    write_ascii_stl(output / "partial_overlap_ascii.stl",
                    "partial_overlap", overlap)

    crossing = [overlap[0],
                ((0.25, 0.1, -1.0), (0.25, 0.1, 1.0),
                 (0.25, 0.8, 0.0))]
    write_ascii_stl(output / "self_intersection_ascii.stl",
                    "self_intersection", crossing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
