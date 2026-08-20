#!/usr/bin/env python3
"""Render deterministic CM2D v1 topology as standalone SVG.

The renderer is a thin post-processing layer: it never rebuilds or reclassifies
mesh geometry. Optional .viz.json metadata is emitted by cartmesh2d_cli from the
already-computed source Cut-cell/small-cell state.
"""

from __future__ import annotations

import argparse
import html
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Vertex:
    id: int
    x: float
    y: float


@dataclass(frozen=True)
class Edge:
    id: int
    v0: int
    v1: int
    owner: int
    neighbour: int | None
    patch: int


@dataclass(frozen=True)
class Cell:
    id: int
    source_id: int
    source_key: int
    area: float
    vertices: tuple[int, ...]
    edges: tuple[int, ...]

    @property
    def level(self) -> int:
        return self.source_key & 0x3F


@dataclass(frozen=True)
class Mesh:
    vertices: tuple[Vertex, ...]
    edges: tuple[Edge, ...]
    cells: tuple[Cell, ...]
    audit: tuple[int, int, int, int, int, int, int]

    @property
    def valid(self) -> bool:
        return not any(self.audit)


class Cm2dError(ValueError):
    pass


def _require(tokens: list[str], index: int, expected: str) -> int:
    if index >= len(tokens) or tokens[index] != expected:
        got = "<eof>" if index >= len(tokens) else tokens[index]
        raise Cm2dError(f"expected {expected!r}, got {got!r}")
    return index + 1


def read_cm2d(path: Path) -> Mesh:
    tokens = path.read_text(encoding="utf-8-sig").split()
    i = 0
    i = _require(tokens, i, "CM2D")
    if i >= len(tokens) or tokens[i] != "1":
        raise Cm2dError("only CM2D v1 is supported")
    i += 1

    i = _require(tokens, i, "VERTICES")
    vertex_count = int(tokens[i]); i += 1
    vertices: list[Vertex] = []
    for expected_id in range(vertex_count):
        vid = int(tokens[i]); x = float(tokens[i + 1]); y = float(tokens[i + 2]); i += 3
        if vid != expected_id:
            raise Cm2dError("vertex ids must be contiguous and deterministic")
        vertices.append(Vertex(vid, x, y))

    i = _require(tokens, i, "EDGES")
    edge_count = int(tokens[i]); i += 1
    edges: list[Edge] = []
    for expected_id in range(edge_count):
        eid = int(tokens[i]); v0 = int(tokens[i + 1]); v1 = int(tokens[i + 2])
        owner = int(tokens[i + 3]); neighbour_raw = int(tokens[i + 4]); patch = int(tokens[i + 5]); i += 6
        if eid != expected_id:
            raise Cm2dError("edge ids must be contiguous and deterministic")
        if not (0 <= v0 < vertex_count and 0 <= v1 < vertex_count):
            raise Cm2dError("edge references invalid vertex")
        if patch not in (0, 1, 2, 3):
            raise Cm2dError("edge contains unsupported patch value")
        edges.append(Edge(eid, v0, v1, owner, None if neighbour_raw < 0 else neighbour_raw, patch))

    i = _require(tokens, i, "CELLS")
    cell_count = int(tokens[i]); i += 1
    cells: list[Cell] = []
    for expected_id in range(cell_count):
        cid = int(tokens[i]); source_id = int(tokens[i + 1]); source_key = int(tokens[i + 2])
        area = float(tokens[i + 3]); nverts = int(tokens[i + 4]); i += 5
        if cid != expected_id or nverts < 3:
            raise Cm2dError("cell id/vertex-loop is invalid")
        vertex_ids = tuple(int(v) for v in tokens[i:i + nverts]); i += nverts
        if any(v < 0 or v >= vertex_count for v in vertex_ids):
            raise Cm2dError("cell references invalid vertex")
        nedges = int(tokens[i]); i += 1
        edge_ids = tuple(int(e) for e in tokens[i:i + nedges]); i += nedges
        if nedges != nverts or any(e < 0 or e >= edge_count for e in edge_ids):
            raise Cm2dError("cell edge-loop is invalid")
        cells.append(Cell(cid, source_id, source_key, area, vertex_ids, edge_ids))

    i = _require(tokens, i, "AUDIT")
    if i + 7 > len(tokens):
        raise Cm2dError("truncated AUDIT record")
    audit_values = tuple(int(v) for v in tokens[i:i + 7]); i += 7
    i = _require(tokens, i, "END")
    if i != len(tokens):
        raise Cm2dError("unexpected tokens after END")

    for edge in edges:
        if not (0 <= edge.owner < cell_count):
            raise Cm2dError("edge owner is invalid")
        if edge.neighbour is not None and not (0 <= edge.neighbour < cell_count):
            raise Cm2dError("edge neighbour is invalid")

    return Mesh(tuple(vertices), tuple(edges), tuple(cells), audit_values)  # type: ignore[arg-type]


def _load_json(path: Path | None) -> dict | None:
    if path is None:
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def _fmt(value: float) -> str:
    return f"{value:.6g}"


def _hsl_for_level(level: int, min_level: int, max_level: int) -> str:
    if max_level <= min_level:
        hue = 205.0
    else:
        hue = 215.0 - 165.0 * (level - min_level) / (max_level - min_level)
    return f"hsl({_fmt(hue)} 72% 82%)"


def _hsl_for_area(area: float, min_area: float, max_area: float) -> str:
    if max_area <= min_area:
        t = 0.5
    else:
        lo = math.log(max(min_area, 1e-300)); hi = math.log(max(max_area, 1e-300))
        t = (math.log(max(area, 1e-300)) - lo) / max(hi - lo, 1e-300)
    return f"hsl({_fmt(205.0 - 145.0 * t)} 70% 83%)"


def _text(x: float, y: float, value: str, css_class: str = "") -> str:
    cls = f' class="{css_class}"' if css_class else ""
    return f'<text x="{_fmt(x)}" y="{_fmt(y)}"{cls}>{html.escape(value)}</text>'


def render_svg(
    mesh: Mesh,
    output: Path,
    *,
    quality: dict | None = None,
    viz: dict | None = None,
    width: int = 1200,
    height: int = 900,
    labels: bool = False,
    color_by: str = "level",
) -> None:
    if width < 400 or height < 300:
        raise ValueError("SVG dimensions are too small")
    if color_by not in {"level", "area", "none"}:
        raise ValueError("color_by must be level, area, or none")
    if not mesh.vertices or not mesh.cells:
        raise ValueError("mesh is empty")

    xs = [v.x for v in mesh.vertices]; ys = [v.y for v in mesh.vertices]
    xmin, xmax = min(xs), max(xs); ymin, ymax = min(ys), max(ys)
    span_x = max(xmax - xmin, 1e-12); span_y = max(ymax - ymin, 1e-12)

    panel_width = 300 if quality is not None else 0
    left = 48.0; right = 38.0 + panel_width; top = 54.0; bottom = 58.0
    plot_w = width - left - right; plot_h = height - top - bottom
    scale = min(plot_w / span_x, plot_h / span_y)
    used_w = span_x * scale; used_h = span_y * scale
    ox = left + (plot_w - used_w) * 0.5; oy = top + (plot_h - used_h) * 0.5

    def project_xy(x: float, y: float) -> tuple[float, float]:
        return ox + (x - xmin) * scale, oy + (ymax - y) * scale

    def project(vid: int) -> tuple[float, float]:
        v = mesh.vertices[vid]
        return project_xy(v.x, v.y)

    levels = [c.level for c in mesh.cells]
    areas = [c.area for c in mesh.cells if c.area > 0.0]
    min_level, max_level = min(levels), max(levels)
    min_area = min(areas) if areas else 0.0; max_area = max(areas) if areas else 0.0
    cut_like_cells = {edge.owner for edge in mesh.edges if edge.patch == 1}

    source_cells = []
    if viz is not None:
        raw = viz.get("source_cells", [])
        if isinstance(raw, list):
            source_cells = [item for item in raw if isinstance(item, dict)]

    invalid = (not mesh.valid) or (quality is not None and quality.get("valid") is False)

    out: list[str] = []
    out.append('<?xml version="1.0" encoding="UTF-8"?>')
    out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    out.append("<style>"
               ".mesh-bg{fill:#fff}.cell{stroke:none}.cell-cut{stroke:#b91c1c;stroke-width:.7}"
               ".edge-internal{stroke:#64748b;stroke-width:.8;opacity:.62}.edge-embedded{stroke:#c62828;stroke-width:2.4}"
               ".edge-domain{stroke:#111827;stroke-width:2.2}.edge-unclassified{stroke:#f59e0b;stroke-width:2.2;stroke-dasharray:5 4}"
               ".source-bg{fill:none;stroke:#94a3b8;stroke-width:.55;stroke-dasharray:2 3;opacity:.5}"
               ".source-small{fill:none;stroke:#7c3aed;stroke-width:2.2;stroke-dasharray:5 3}.small-dot{fill:#7c3aed;stroke:#fff;stroke-width:1}"
               ".label{font:10px ui-monospace,monospace;fill:#111827;text-anchor:middle;dominant-baseline:middle}"
               ".title{font:600 18px system-ui,sans-serif;fill:#111827}.small{font:12px system-ui,sans-serif;fill:#334155}"
               ".panel-title{font:600 14px system-ui,sans-serif;fill:#111827}.panel{font:12px ui-monospace,monospace;fill:#334155}"
               ".legend{font:11px system-ui,sans-serif;fill:#334155}.invalid{font:700 16px system-ui,sans-serif;fill:#b91c1c}"
               "</style>")
    out.append(f'<rect class="mesh-bg" x="0" y="0" width="{width}" height="{height}"/>')
    out.append(_text(left, 29, "cartmesh2d — stabilized solver topology", "title"))
    if invalid:
        out.append(_text(left, 48, "INVALID TOPOLOGY / QUALITY — inspect audit", "invalid"))

    if source_cells:
        out.append('<g id="source-background">')
        for item in source_cells:
            bounds = item.get("background_bounds")
            if not isinstance(bounds, list) or len(bounds) != 4:
                continue
            x0, y0, x1, y1 = (float(v) for v in bounds)
            p0 = project_xy(x0, y1); p1 = project_xy(x1, y0)
            cls = "source-small" if item.get("small") is True else "source-bg"
            out.append(f'<rect class="{cls}" data-source="{item.get("source_id", "")}" data-kind="{html.escape(str(item.get("kind", "")))}" x="{_fmt(p0[0])}" y="{_fmt(p0[1])}" width="{_fmt(p1[0]-p0[0])}" height="{_fmt(p1[1]-p0[1])}"/>')
            centroid = item.get("centroid")
            if item.get("small") is True and isinstance(centroid, list) and len(centroid) == 2:
                cx, cy = project_xy(float(centroid[0]), float(centroid[1]))
                out.append(f'<circle class="small-dot" cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="4"/>')
        out.append("</g>")

    out.append('<g id="cells">')
    for cell in mesh.cells:
        points = " ".join(f"{_fmt(x)},{_fmt(y)}" for x, y in (project(v) for v in cell.vertices))
        if color_by == "level":
            fill = _hsl_for_level(cell.level, min_level, max_level)
        elif color_by == "area":
            fill = _hsl_for_area(cell.area, min_area, max_area)
        else:
            fill = "#f8fafc"
        classes = "cell cell-cut" if cell.id in cut_like_cells else "cell"
        out.append(f'<polygon class="{classes}" data-cell="{cell.id}" data-source="{cell.source_id}" data-level="{cell.level}" points="{points}" fill="{fill}"/>')
    out.append("</g>")

    groups = {0: [], 1: [], 2: [], 3: []}
    for edge in mesh.edges:
        groups[edge.patch].append(edge)
    out.append('<g id="edges">')
    for patch, css in ((0, "edge-internal"), (3, "edge-unclassified"), (2, "edge-domain"), (1, "edge-embedded")):
        for edge in groups[patch]:
            x1, y1 = project(edge.v0); x2, y2 = project(edge.v1)
            out.append(f'<line class="{css}" data-edge="{edge.id}" x1="{_fmt(x1)}" y1="{_fmt(y1)}" x2="{_fmt(x2)}" y2="{_fmt(y2)}"/>')
    out.append("</g>")

    if labels:
        out.append('<g id="labels">')
        for cell in mesh.cells:
            pts = [project(v) for v in cell.vertices]
            cx = sum(p[0] for p in pts) / len(pts); cy = sum(p[1] for p in pts) / len(pts)
            out.append(_text(cx, cy, str(cell.id), "label"))
        out.append("</g>")

    legend_y = height - 28
    legend_items = [
        ("edge-internal", "internal", 0),
        ("edge-embedded", "embedded boundary", 100),
        ("edge-domain", "domain boundary", 270),
    ]
    for css, label, offset in legend_items:
        out.append(f'<line class="{css}" x1="{left+offset}" y1="{legend_y}" x2="{left+offset+30}" y2="{legend_y}"/>')
        out.append(_text(left + offset + 36, legend_y + 4, label, "legend"))
    if source_cells:
        out.append(f'<rect class="source-small" x="{left+420}" y="{legend_y-7}" width="18" height="14"/>')
        out.append(_text(left + 444, legend_y + 4, "source small-cell (pre-agglomeration)", "legend"))

    if quality is not None:
        px = width - panel_width + 18
        out.append(f'<line x1="{width-panel_width}" y1="20" x2="{width-panel_width}" y2="{height-20}" stroke="#cbd5e1"/>')
        out.append(_text(px, 48, "Quality summary", "panel-title"))
        counts = quality.get("counts", {}) if isinstance(quality.get("counts", {}), dict) else {}
        q = quality.get("quality", {}) if isinstance(quality.get("quality", {}), dict) else {}
        audit = quality.get("topology_audit", {}) if isinstance(quality.get("topology_audit", {}), dict) else {}
        lines = [
            f"valid: {quality.get('valid', '?')}",
            f"vertices: {counts.get('vertices', '?')}",
            f"edges: {counts.get('edges', '?')}",
            f"cells: {counts.get('cells', '?')}",
            f"cut source: {counts.get('source_cut_cells', '?')}",
            f"small source: {counts.get('source_small_cells', '?')}",
            f"min area: {q.get('min_cell_area', '?')}",
            f"min edge: {q.get('min_edge_length', '?')}",
            f"max aspect: {q.get('max_edge_aspect_ratio', '?')}",
            f"max skew: {q.get('max_centroid_skewness', '?')}",
            "", "topology audit:",
        ]
        for key in ("duplicate_vertices", "duplicate_edges", "orphan_internal_edges", "non_manifold_edges", "unclassified_boundary_edges", "open_cell_loops", "area_mismatches"):
            lines.append(f"  {key}: {audit.get(key, '?')}")
        if viz is not None:
            lines.extend(["", f"viz small: {viz.get('source_small_cell_count', '?')}",
                          f"viz merged: {viz.get('merged_small_cell_count', '?')}"])
        for n, line in enumerate(lines):
            out.append(_text(px, 72 + n * 19, line, "panel"))

    out.append(_text(left, height - 7, f"bbox=({_fmt(xmin)},{_fmt(ymin)})..({_fmt(xmax)},{_fmt(ymax)})  levels={min_level}..{max_level}  cells={len(mesh.cells)}  edges={len(mesh.edges)}", "small"))
    out.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render CM2D v1 topology to standalone SVG")
    parser.add_argument("cm2d", type=Path)
    parser.add_argument("svg", type=Path)
    parser.add_argument("--quality", type=Path)
    parser.add_argument("--viz", type=Path, help="optional cartmesh2d-viz-v1 source-cell metadata")
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--labels", action="store_true")
    parser.add_argument("--color-by", choices=("level", "area", "none"), default="level")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    mesh = read_cm2d(args.cm2d)
    quality = _load_json(args.quality)
    viz = _load_json(args.viz)
    if viz is not None and viz.get("format") != "cartmesh2d-viz-v1":
        raise ValueError("unsupported visualization sidecar format")
    render_svg(mesh, args.svg, quality=quality, viz=viz, width=args.width, height=args.height,
               labels=args.labels, color_by=args.color_by)
    print(f"cartmesh2d visualization PASS: {args.svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
