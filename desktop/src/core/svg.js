'use strict';

// SVG -> closed polylines.
//
// SVG is the format a non-CAD user actually has: anything drawn in Illustrator,
// Inkscape or Figma exports one.  Curves are flattened to a chord tolerance the
// same way the DXF converter flattens arcs and splines, so the discretisation
// error entering the mesher is a stated number rather than whatever the drawing
// tool happened to emit.
//
// Deliberately not supported: `transform`, `use`, nested groups with their own
// coordinate systems, and stroke geometry.  Those change coordinates silently,
// and a boundary that is silently in the wrong place is worse than a refusal.

const COMMAND = /([MmLlHhVvCcSsQqTtAaZz])([^MmLlHhVvCcSsQqTtAaZz]*)/g;
const NUMBER = /-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g;

const numbers = text => (text.match(NUMBER) || []).map(Number);

// Recursive de Casteljau until the control points sit within `tolerance` of the
// chord.  Subdividing on measured flatness keeps a gentle curve cheap and a tight
// one dense, which is what a chord-error contract means.
function flattenCubic(p0, p1, p2, p3, tolerance, out, depth = 0) {
  const dx = p3[0] - p0[0];
  const dy = p3[1] - p0[1];
  const chord = Math.hypot(dx, dy);
  const deviation = chord > 0
    ? Math.max(Math.abs((p1[0] - p0[0]) * dy - (p1[1] - p0[1]) * dx) / chord,
               Math.abs((p2[0] - p0[0]) * dy - (p2[1] - p0[1]) * dx) / chord)
    : Math.max(Math.hypot(p1[0] - p0[0], p1[1] - p0[1]),
               Math.hypot(p2[0] - p0[0], p2[1] - p0[1]));
  if (deviation <= tolerance || depth >= 20) {
    out.push(p3);
    return;
  }
  const mid = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
  const a = mid(p0, p1), b = mid(p1, p2), c = mid(p2, p3);
  const d = mid(a, b), e = mid(b, c), f = mid(d, e);
  flattenCubic(p0, a, d, f, tolerance, out, depth + 1);
  flattenCubic(f, e, c, p3, tolerance, out, depth + 1);
}

const quadraticToCubic = (p0, q, p1) => [
  [p0[0] + (2 / 3) * (q[0] - p0[0]), p0[1] + (2 / 3) * (q[1] - p0[1])],
  [p1[0] + (2 / 3) * (q[0] - p1[0]), p1[1] + (2 / 3) * (q[1] - p1[1])]
];

// SVG endpoint-parameterised elliptical arc -> centre parameterisation, then a
// uniform angular sweep whose step is chosen from the sagitta so the chord error
// matches the Bézier path.
function flattenArc(p0, rx, ry, rotationDeg, largeArc, sweep, p1, tolerance, out) {
  if (!(rx > 0) || !(ry > 0)) { out.push(p1); return; }
  const phi = (rotationDeg * Math.PI) / 180;
  const cos = Math.cos(phi), sin = Math.sin(phi);
  const dx2 = (p0[0] - p1[0]) / 2, dy2 = (p0[1] - p1[1]) / 2;
  const x1 = cos * dx2 + sin * dy2;
  const y1 = -sin * dx2 + cos * dy2;
  let ax = rx, ay = ry;
  const oversize = (x1 * x1) / (ax * ax) + (y1 * y1) / (ay * ay);
  if (oversize > 1) { const s = Math.sqrt(oversize); ax *= s; ay *= s; }
  const numerator = ax * ax * ay * ay - ax * ax * y1 * y1 - ay * ay * x1 * x1;
  const denominator = ax * ax * y1 * y1 + ay * ay * x1 * x1;
  const factor = (largeArc === sweep ? -1 : 1) *
    Math.sqrt(Math.max(0, numerator) / (denominator || 1));
  const cx1 = (factor * ax * y1) / ay;
  const cy1 = (-factor * ay * x1) / ax;
  const cx = cos * cx1 - sin * cy1 + (p0[0] + p1[0]) / 2;
  const cy = sin * cx1 + cos * cy1 + (p0[1] + p1[1]) / 2;
  const angleOf = (ux, uy) => Math.atan2(uy, ux);
  const start = angleOf((x1 - cx1) / ax, (y1 - cy1) / ay);
  let sweepAngle = angleOf((-x1 - cx1) / ax, (-y1 - cy1) / ay) - start;
  if (!sweep && sweepAngle > 0) sweepAngle -= 2 * Math.PI;
  if (sweep && sweepAngle < 0) sweepAngle += 2 * Math.PI;
  const radius = Math.max(ax, ay);
  const step = 2 * Math.acos(Math.max(-1, Math.min(1, 1 - tolerance / radius)));
  const count = Math.max(2, Math.ceil(Math.abs(sweepAngle) / Math.max(step, 1e-4)));
  for (let i = 1; i <= count; i++) {
    const angle = start + (sweepAngle * i) / count;
    const ex = ax * Math.cos(angle), ey = ay * Math.sin(angle);
    out.push([cos * ex - sin * ey + cx, sin * ex + cos * ey + cy]);
  }
}

// One `d` attribute may hold several subpaths; each `M` starts a new one and each
// subpath becomes its own loop.  Only closed subpaths are kept: an open stretch is
// not a boundary and the mesher would reject it later with less context.
function walkPath(d, tolerance, loops, warnings) {
  let current = [0, 0];
  let start = [0, 0];
  let points = null;
  let control = null;
  let previous = '';
  const finish = closed => {
    if (points && points.length >= 3) {
      if (closed) loops.push(points);
      else warnings.push(`跳过一条未闭合的 path 子路径（${points.length} 个点）`);
    }
    points = null;
  };
  const move = to => { current = to; points.push(to.slice()); };
  const cubic = (c1, c2, to) => {
    flattenCubic(current, c1, c2, to, tolerance, points);
    control = c2;
    current = to;
  };
  COMMAND.lastIndex = 0;
  let match;
  while ((match = COMMAND.exec(d)) !== null) {
    const raw = match[1];
    const command = raw.toUpperCase();
    const relative = raw !== command;
    const args = numbers(match[2]);
    const point = i => relative
      ? [current[0] + args[i], current[1] + args[i + 1]]
      : [args[i], args[i + 1]];

    if (command === 'Z') {
      finish(true);
      current = start.slice();
      control = null;
      previous = command;
      continue;
    }
    if (command === 'M') {
      finish(false);
      current = point(0);
      start = current.slice();
      points = [current.slice()];
      for (let i = 2; i + 1 < args.length; i += 2) move(point(i));
      control = null;
      previous = command;
      continue;
    }
    if (!points) { points = [current.slice()]; start = current.slice(); }

    if (command === 'L') {
      for (let i = 0; i + 1 < args.length; i += 2) move(point(i));
      control = null;
    } else if (command === 'H') {
      for (const value of args) move([relative ? current[0] + value : value, current[1]]);
      control = null;
    } else if (command === 'V') {
      for (const value of args) move([current[0], relative ? current[1] + value : value]);
      control = null;
    } else if (command === 'C') {
      for (let i = 0; i + 5 < args.length; i += 6) cubic(point(i), point(i + 2), point(i + 4));
    } else if (command === 'S') {
      // The first control point mirrors the previous one only after another cubic;
      // otherwise it coincides with the current point.
      for (let i = 0; i + 3 < args.length; i += 4) {
        const mirrored = (previous === 'C' || previous === 'S') && control
          ? [2 * current[0] - control[0], 2 * current[1] - control[1]]
          : current.slice();
        cubic(mirrored, point(i), point(i + 2));
      }
    } else if (command === 'Q') {
      for (let i = 0; i + 3 < args.length; i += 4) {
        const q = point(i), to = point(i + 2);
        const [c1, c2] = quadraticToCubic(current, q, to);
        flattenCubic(current, c1, c2, to, tolerance, points);
        control = q;
        current = to;
      }
    } else if (command === 'T') {
      for (let i = 0; i + 1 < args.length; i += 2) {
        const q = (previous === 'Q' || previous === 'T') && control
          ? [2 * current[0] - control[0], 2 * current[1] - control[1]]
          : current.slice();
        const to = point(i);
        const [c1, c2] = quadraticToCubic(current, q, to);
        flattenCubic(current, c1, c2, to, tolerance, points);
        control = q;
        current = to;
      }
    } else if (command === 'A') {
      for (let i = 0; i + 6 < args.length; i += 7) {
        const to = relative
          ? [current[0] + args[i + 5], current[1] + args[i + 6]]
          : [args[i + 5], args[i + 6]];
        flattenArc(current, Math.abs(args[i]), Math.abs(args[i + 1]), args[i + 2],
                   args[i + 3] !== 0, args[i + 4] !== 0, to, tolerance, points);
        current = to;
        control = null;
      }
    }
    previous = command;
  }
  finish(false);
}

const attribute = (tag, name) => {
  const match = tag.match(new RegExp(`\\b${name}\\s*=\\s*"([^"]*)"`)) ||
                tag.match(new RegExp(`\\b${name}\\s*=\\s*'([^']*)'`));
  return match ? match[1] : null;
};

const ellipseLoop = (cx, cy, rx, ry, tolerance) => {
  const radius = Math.max(rx, ry);
  const step = 2 * Math.acos(Math.max(-1, Math.min(1, 1 - tolerance / radius)));
  const count = Math.max(8, Math.ceil((2 * Math.PI) / Math.max(step, 1e-4)));
  return Array.from({ length: count }, (_, i) => {
    const angle = (2 * Math.PI * i) / count;
    return [cx + rx * Math.cos(angle), cy + ry * Math.sin(angle)];
  });
};

// A <rect> with rounded corners is a different shape from its corner points, so
// the radii are refused rather than quietly ignored.
function rectLoop(tag, warnings) {
  const value = name => Number(attribute(tag, name) || 0);
  if (Number(attribute(tag, 'rx') || 0) > 0 || Number(attribute(tag, 'ry') || 0) > 0) {
    warnings.push('rect 的圆角半径未被支持，已跳过该元素');
    return null;
  }
  const x = value('x'), y = value('y');
  const w = value('width'), h = value('height');
  if (!(w > 0) || !(h > 0)) return null;
  return [[x, y], [x + w, y], [x + w, y + h], [x, y + h]];
}

// SVG y grows downward.  Flipping it keeps the mesh looking like the drawing; the
// reflection reverses loop orientation, which BoundaryRegion2D::normalizeAlternating
// fixes anyway.
function parseSvgLoops(text, { chordToleranceFraction = 0.002 } = {}) {
  const warnings = [];
  if (/\btransform\s*=/.test(text)) {
    warnings.push('文件里有 transform 属性，本导入器不解析它；坐标可能与视觉位置不符');
  }
  const raw = [];
  const collect = (tag, loop) => { if (loop && loop.length >= 3) raw.push(loop); };

  // Two passes: geometry first at unit tolerance to learn the span, then curves are
  // flattened against a tolerance derived from it.  Without that the tolerance would
  // depend on whatever coordinate range the drawing happens to use.
  const spanProbe = [];
  for (const match of text.matchAll(/<(path|polygon|polyline|circle|ellipse|rect)\b([^>]*)>/gi)) {
    spanProbe.push(...numbers(match[2]));
  }
  const finite = spanProbe.filter(Number.isFinite);
  const extent = finite.length ? Math.max(...finite) - Math.min(...finite) : 1;
  const tolerance = Math.max(1e-9, chordToleranceFraction * (extent || 1));

  for (const match of text.matchAll(/<(path|polygon|polyline|circle|ellipse|rect)\b([^>]*)>/gi)) {
    const [, element, tag] = match;
    const kind = element.toLowerCase();
    if (kind === 'path') {
      const d = attribute(tag, 'd');
      if (d) walkPath(d, tolerance, raw, warnings);
    } else if (kind === 'polygon' || kind === 'polyline') {
      const values = numbers(attribute(tag, 'points') || '');
      const loop = [];
      for (let i = 0; i + 1 < values.length; i += 2) loop.push([values[i], values[i + 1]]);
      if (kind === 'polyline' && loop.length >= 3) {
        const [ax, ay] = loop[0];
        const [bx, by] = loop[loop.length - 1];
        if (Math.hypot(bx - ax, by - ay) > tolerance) {
          warnings.push('polyline 首尾不重合，已按闭合处理');
        }
      }
      collect(tag, loop);
    } else if (kind === 'circle') {
      const r = Number(attribute(tag, 'r') || 0);
      if (r > 0) collect(tag, ellipseLoop(Number(attribute(tag, 'cx') || 0),
                                         Number(attribute(tag, 'cy') || 0), r, r, tolerance));
    } else if (kind === 'ellipse') {
      const rx = Number(attribute(tag, 'rx') || 0);
      const ry = Number(attribute(tag, 'ry') || 0);
      if (rx > 0 && ry > 0) collect(tag, ellipseLoop(Number(attribute(tag, 'cx') || 0),
                                                    Number(attribute(tag, 'cy') || 0), rx, ry, tolerance));
    } else if (kind === 'rect') {
      collect(tag, rectLoop(tag, warnings));
    }
  }
  return { loops: raw.map(loop => loop.map(([x, y]) => [x, -y])), warnings };
}

module.exports = { parseSvgLoops, flattenCubic, flattenArc };
