'use strict';

const path = require('node:path');
const { parseSvgLoops } = require('./svg');

// Geometry input, one place.  Everything that is not DXF is converted to loops here
// and written out as the native .xy; DXF keeps going through the fail-closed C++
// converter because that is where unit handling and entity diagnostics live.

const NUMBER_PAIR = /^\s*(-?[\d.eE+]+)[\s,;]+(-?[\d.eE+]+)\s*$/;

// Shared by .xy, .csv and .txt: a blank line ends a loop, `#` is a comment.  The
// separator is whatever sits between the two numbers, so a comma-delimited export
// and a space-delimited one are the same file to us.
function parseCoordinateText(text) {
  const loops = [];
  let loop = [];
  const warnings = [];
  const lines = text.split(/\r?\n/);
  lines.forEach((line, index) => {
    const trimmed = line.trim();
    if (!trimmed) { if (loop.length) { loops.push(loop); loop = []; } return; }
    if (trimmed.startsWith('#') || trimmed.startsWith('//')) return;
    const match = trimmed.match(NUMBER_PAIR);
    if (!match) {
      // A header row is the normal first line of a CSV export, so skip one silently
      // and complain about anything after that.
      if (index === 0) return;
      warnings.push(`第 ${index + 1} 行不是坐标对，已跳过：${trimmed.slice(0, 40)}`);
      return;
    }
    const x = Number(match[1]);
    const y = Number(match[2]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) {
      warnings.push(`第 ${index + 1} 行坐标不是有限值，已跳过`);
      return;
    }
    loop.push([x, y]);
  });
  if (loop.length) loops.push(loop);
  return { loops, warnings };
}

// Drop a repeated closing vertex: the mesher closes every loop implicitly, and a
// duplicated endpoint would be a zero-length wall segment.
function dropClosingDuplicate(loop) {
  if (loop.length < 2) return loop;
  const [ax, ay] = loop[0];
  const [bx, by] = loop[loop.length - 1];
  const span = Math.max(1e-300, Math.hypot(bx - ax, by - ay));
  const scale = loop.reduce((acc, [x, y]) => Math.max(acc, Math.abs(x), Math.abs(y)), 1);
  return span <= 1e-12 * scale ? loop.slice(0, -1) : loop;
}

function validateLoops(loops) {
  const issues = [];
  const cleaned = loops.map(dropClosingDuplicate).filter(loop => loop.length >= 3);
  if (!cleaned.length) issues.push('文件里没有找到至少 3 个顶点的闭合环。');
  cleaned.forEach((loop, index) => {
    const area = loop.reduce((sum, [x, y], i) => {
      const [nx, ny] = loop[(i + 1) % loop.length];
      return sum + (x * ny - nx * y);
    }, 0) / 2;
    if (!(Math.abs(area) > 0)) issues.push(`第 ${index + 1} 个环的面积为零。`);
  });
  return { loops: cleaned, issues };
}

// The size field asks for everything in body spans, so absolute scale only matters
// for the hybrid path's first-layer thickness.  Normalising an unitless format to a
// 1 m body span therefore keeps every sizing default meaningful.
function normalizeToUnitSpan(loops) {
  const xs = loops.flat().map(p => p[0]);
  const ys = loops.flat().map(p => p[1]);
  const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys));
  if (!(span > 0)) return { loops, scale: 1 };
  const scale = 1 / span;
  return { loops: loops.map(loop => loop.map(([x, y]) => [x * scale, y * scale])), scale };
}

function loopsToXyText(loops, header) {
  const lines = header ? [`# ${header}`] : [];
  loops.forEach((loop, index) => {
    if (index > 0) lines.push('');
    for (const [x, y] of loop) lines.push(`${x.toPrecision(12)} ${y.toPrecision(12)}`);
  });
  return `${lines.join('\n')}\n`;
}

const classify = filePath => {
  const extension = path.extname(filePath).toLowerCase().replace('.', '');
  if (extension === 'dxf') return 'dxf';
  if (extension === 'svg') return 'svg';
  if (extension === 'xy') return 'xy';
  if (extension === 'csv' || extension === 'txt' || extension === 'dat') return 'text';
  return null;
};

// Convert in-process.  Returns null for DXF: that one needs the C++ converter, so
// the caller runs it instead of us guessing at unit handling.
function convertToLoops(filePath, text, options = {}) {
  const kind = classify(filePath);
  if (kind === null) {
    return { kind: null, issues: [`不支持的文件类型：${path.extname(filePath) || '(无扩展名)'}`] };
  }
  if (kind === 'dxf') return { kind, loops: null, issues: [], warnings: [] };

  const parsed = kind === 'svg'
    ? parseSvgLoops(text, options)
    : parseCoordinateText(text);
  const { loops, issues } = validateLoops(parsed.loops);
  if (issues.length) return { kind, issues, warnings: parsed.warnings || [] };

  // .xy is already in metres by contract; the unitless formats get normalised.
  const scaled = kind === 'xy' ? { loops, scale: 1 } : normalizeToUnitSpan(loops);
  return {
    kind,
    loops: scaled.loops,
    scale: scaled.scale,
    normalized: kind !== 'xy',
    issues: [],
    warnings: parsed.warnings || []
  };
}

module.exports = {
  classify,
  convertToLoops,
  parseCoordinateText,
  validateLoops,
  normalizeToUnitSpan,
  loopsToXyText,
  dropClosingDuplicate
};
