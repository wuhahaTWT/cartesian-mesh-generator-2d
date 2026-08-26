#!/usr/bin/env python3
"""Render an H4-1 debug VTK artifact as deterministic standalone SVG."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path


def read_points_cells(path: Path) -> tuple[list[tuple[float, float]], list[tuple[int, ...]]]:
    tokens = path.read_text(encoding="utf-8").split()
    point_at = tokens.index("POINTS")
    point_count = int(tokens[point_at + 1])
    start = point_at + 3
    points = [
        (float(tokens[start + 3 * i]), float(tokens[start + 3 * i + 1]))
        for i in range(point_count)
    ]
    cell_at = tokens.index("CELLS", start + 3 * point_count)
    cell_count = int(tokens[cell_at + 1])
    cursor = cell_at + 3
    cells: list[tuple[int, ...]] = []
    for _ in range(cell_count):
        count = int(tokens[cursor])
        cursor += 1
        cells.append(tuple(int(value) for value in tokens[cursor : cursor + count]))
        cursor += count
    return points, cells


def fmt(value: float) -> str:
    return f"{value:.8g}"


def render(vtk: Path, report_path: Path, output: Path) -> None:
    points, cells = read_points_cells(vtk)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    strips = report.get("strips", [])
    if report.get("layer_status") != "success" or len(strips) != 1:
        raise ValueError("renderer currently requires one successful H4-1 strip")
    strip = strips[0]
    wall_count = int(strip["wall_vertex_count"])
    layers = int(strip["n_layers"])
    if len(points) != wall_count * (layers + 1):
        raise ValueError("ring-major point layout disagrees with JSON report")

    width, height = 1200, 900
    panel = 300
    left, right, top, bottom = 54.0, 40.0 + panel, 64.0, 56.0
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    span_x = max(max_x - min_x, 1.0e-12)
    span_y = max(max_y - min_y, 1.0e-12)
    scale = min((width - left - right) / span_x, (height - top - bottom) / span_y)
    used_x, used_y = span_x * scale, span_y * scale
    origin_x = left + ((width - left - right) - used_x) * 0.5
    origin_y = top + ((height - top - bottom) - used_y) * 0.5

    def project(point: tuple[float, float]) -> tuple[float, float]:
        return (
            origin_x + (point[0] - min_x) * scale,
            origin_y + (max_y - point[1]) * scale,
        )

    def polyline(ids: list[int], close: bool = True) -> str:
        projected = [project(points[index]) for index in ids]
        if close:
            projected.append(projected[0])
        return " ".join(f"{fmt(x)},{fmt(y)}" for x, y in projected)

    out: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>"
        ".bg{fill:#f8fafc}.cell{stroke:#64748b;stroke-width:.65}.hair{stroke:#64748b;stroke-width:.7;stroke-dasharray:3 3;opacity:.7}"
        ".ring{fill:none;stroke:#2563eb;stroke-width:1.2}.wall{fill:none;stroke:#b91c1c;stroke-width:3}.outer{fill:none;stroke:#0f172a;stroke-width:2.5}"
        ".title{font:600 20px system-ui,sans-serif;fill:#111827}.subtitle{font:13px system-ui,sans-serif;fill:#475569}"
        ".panel-title{font:600 14px system-ui,sans-serif;fill:#111827}.panel{font:12px ui-monospace,monospace;fill:#334155}"
        ".legend{font:12px system-ui,sans-serif;fill:#334155}"
        "</style>",
        f'<rect class="bg" x="0" y="0" width="{width}" height="{height}"/>',
        '<text class="title" x="54" y="31">CartMesh2D H4-1 — body-fitted quad wrapper strip</text>',
        '<text class="subtitle" x="54" y="51">isolated candidate; no Cartesian/Cut-cell remainder integration</text>',
    ]
    colors = ["#dbeafe", "#bfdbfe", "#93c5fd", "#60a5fa", "#3b82f6", "#2563eb"]
    for index, cell in enumerate(cells):
        layer = index // wall_count
        fill = colors[min(layer, len(colors) - 1)]
        out.append(
            f'<polygon class="cell" fill="{fill}" points="{polyline(list(cell), close=False)}"/>'
        )
    for vertex in range(wall_count):
        ids = [ring * wall_count + vertex for ring in range(layers + 1)]
        out.append(f'<polyline class="hair" points="{polyline(ids, close=False)}"/>')
    for ring in range(layers + 1):
        ids = [ring * wall_count + vertex for vertex in range(wall_count)]
        css = "wall" if ring == 0 else "outer" if ring == layers else "ring"
        out.append(f'<polyline class="{css}" points="{polyline(ids)}"/>')

    panel_x = width - panel - 20
    out.append(f'<text class="panel-title" x="{panel_x}" y="92">Verified H4-1 artifact</text>')
    rows = [
        ("patch", str(strip["patch"])),
        ("wall vertices", str(wall_count)),
        ("layers", str(layers)),
        ("quad cells", str(report["layer_cell_count"])),
        ("layer vertices", str(report["layer_vertex_count"])),
        ("first thickness", fmt(float(strip["first_layer_thickness"]))),
        ("total thickness", fmt(float(strip["total_thickness"]))),
        ("growth ratio", fmt(float(strip["growth_ratio"]))),
        ("min cell area", fmt(float(report["min_cell_area"]))),
        ("max cell area", fmt(float(report["max_cell_area"]))),
        ("min hair spacing", fmt(float(report["min_layer_thickness"]))),
        ("max hair spacing", fmt(float(report["max_layer_thickness"]))),
    ]
    for row, (key, value) in enumerate(rows):
        y = 120 + row * 24
        out.append(
            f'<text class="panel" x="{panel_x}" y="{y}">{html.escape(key)}: {html.escape(value)}</text>'
        )
    legend_y = 445
    out.extend(
        [
            f'<line class="wall" x1="{panel_x}" y1="{legend_y}" x2="{panel_x + 38}" y2="{legend_y}"/>',
            f'<text class="legend" x="{panel_x + 48}" y="{legend_y + 4}">wall</text>',
            f'<line class="ring" x1="{panel_x}" y1="{legend_y + 28}" x2="{panel_x + 38}" y2="{legend_y + 28}"/>',
            f'<text class="legend" x="{panel_x + 48}" y="{legend_y + 32}">layer edge</text>',
            f'<line class="hair" x1="{panel_x}" y1="{legend_y + 56}" x2="{panel_x + 38}" y2="{legend_y + 56}"/>',
            f'<text class="legend" x="{panel_x + 48}" y="{legend_y + 60}">hair edge</text>',
            f'<line class="outer" x1="{panel_x}" y1="{legend_y + 84}" x2="{panel_x + 38}" y2="{legend_y + 84}"/>',
            f'<text class="legend" x="{panel_x + 48}" y="{legend_y + 88}">outer envelope</text>',
        ]
    )
    out.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(out) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtk", type=Path)
    parser.add_argument("json", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.vtk, args.json, args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
