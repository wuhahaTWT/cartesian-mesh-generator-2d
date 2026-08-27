#!/usr/bin/env python3
"""Render a deterministic H4-2 hybrid VTK artifact to standalone SVG."""

from __future__ import annotations

import argparse
import html
import json
from collections import defaultdict
from pathlib import Path


def read_vtk(path: Path):
    tokens = path.read_text(encoding="utf-8").split()
    at = tokens.index("POINTS"); count = int(tokens[at + 1]); start = at + 3
    points = [(float(tokens[start + 3*i]), float(tokens[start + 3*i + 1])) for i in range(count)]
    at = tokens.index("CELLS", start + 3*count); cell_count = int(tokens[at + 1]); cursor = at + 3
    cells = []
    for _ in range(cell_count):
        size = int(tokens[cursor]); cursor += 1
        cells.append(tuple(int(v) for v in tokens[cursor:cursor + size])); cursor += size
    at = tokens.index("hybrid_kind", tokens.index("CELL_DATA", cursor))
    at = tokens.index("LOOKUP_TABLE", at) + 2
    kinds = [int(v) for v in tokens[at:at + cell_count]]
    at = tokens.index("layer_index", at + cell_count)
    at = tokens.index("LOOKUP_TABLE", at) + 2
    layers = [int(v) for v in tokens[at:at + cell_count]]
    return points, cells, kinds, layers


def render(vtk: Path, report_path: Path, output: Path) -> None:
    points, cells, kinds, layers = read_vtk(vtk)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    width, height, panel = 1400, 980, 330
    left, right, top, bottom = 52.0, panel + 35.0, 68.0, 48.0
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    span_x, span_y = max(xs)-min(xs), max(ys)-min(ys)
    scale = min((width-left-right)/span_x, (height-top-bottom)/span_y)
    ox = left + ((width-left-right)-span_x*scale)/2
    oy = top + ((height-top-bottom)-span_y*scale)/2
    def project(p): return ox+(p[0]-min(xs))*scale, oy+(max(ys)-p[1])*scale
    def polygon(ids): return " ".join(f"{project(points[v])[0]:.8g},{project(points[v])[1]:.8g}" for v in ids)
    owners = defaultdict(list)
    for cid, cell in enumerate(cells):
        for i, first in enumerate(cell): owners[tuple(sorted((first, cell[(i+1)%len(cell)])))].append(cid)
    interface = [edge for edge, own in owners.items()
                 if len(own) == 2 and ((kinds[own[0]] == 0) != (kinds[own[1]] == 0))]
    wall = [edge for edge, own in owners.items() if len(own) == 1 and kinds[own[0]] == 0]
    colors = ["#dbeafe", "#bfdbfe", "#93c5fd", "#60a5fa", "#3b82f6", "#2563eb"]
    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="52" y="30" font-family="system-ui" font-size="21" font-weight="650" fill="#111827">CartMesh2D H4-2 — conformal boundary-layer / Cut-cell mesh</text>',
        '<text x="52" y="52" font-family="system-ui" font-size="13" fill="#475569">one global topology; highlighted outer envelope is shared by both sides</text>',
    ]
    order = sorted(range(len(cells)), key=lambda i: kinds[i] == 0)
    for cid in order:
        if kinds[cid] == 0:
            fill = colors[min(max(layers[cid], 0), len(colors)-1)]
        elif kinds[cid] == 1:
            fill = "#fed7aa"
        elif kinds[cid] == 3:
            fill = "#fde68a"
        else:
            fill = "#ecfccb"
        out.append(f'<polygon points="{polygon(cells[cid])}" fill="{fill}" stroke="#64748b" stroke-width="0.55"/>')
    for first, second in interface:
        a, b = project(points[first]), project(points[second])
        out.append(f'<line x1="{a[0]:.8g}" y1="{a[1]:.8g}" x2="{b[0]:.8g}" y2="{b[1]:.8g}" stroke="#111827" stroke-width="3"/>')
    for first, second in wall:
        a, b = project(points[first]), project(points[second])
        out.append(f'<line x1="{a[0]:.8g}" y1="{a[1]:.8g}" x2="{b[0]:.8g}" y2="{b[1]:.8g}" stroke="#b91c1c" stroke-width="3.4"/>')
    px = width-panel
    out.append(f'<text x="{px}" y="92" font-family="system-ui" font-size="15" font-weight="650">H4-2 topology audit</text>')
    rows = [
        ("layer cells", report["boundary_layer_cell_count"]),
        ("transition cells", report["transition_polygon_count"]),
        ("remainder cut cells", report["remainder_cut_cell_count"]),
        ("Cartesian cells", report["remainder_cartesian_cell_count"]),
        ("shared interface edges", report["interface_edge_count"]),
        ("single-owner interface", report["single_owner_interface_edges"]),
        ("wrong interface pairs", report["wrong_cell_pair_interface_edges"]),
        ("non-two-valent vertices", report["non_two_valent_interface_vertices"]),
        ("area error", f'{float(report["area_error"]):.4e}'),
        ("topology", "PASS" if report["topology_valid"] else "FAIL"),
    ]
    for i, (key, value) in enumerate(rows):
        out.append(f'<text x="{px}" y="{124+i*25}" font-family="ui-monospace,monospace" font-size="12" fill="#334155">{html.escape(str(key))}: {html.escape(str(value))}</text>')
    legend = [("#b91c1c", "wall"), ("#60a5fa", "boundary layers"), ("#111827", "shared outer envelope"), ("#fde68a", "graded transition fan"), ("#fed7aa", "remainder cut cells"), ("#ecfccb", "Cartesian region")]
    y = 430
    for color, label in legend:
        out.append(f'<rect x="{px}" y="{y-11}" width="34" height="14" fill="{color}"/><text x="{px+46}" y="{y}" font-family="system-ui" font-size="12" fill="#334155">{label}</text>')
        y += 30
    out.append('</svg>')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(out)+"\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtk", type=Path); parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path); args = parser.parse_args()
    render(args.vtk, args.report, args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
