'use strict';

const path = require('node:path');
const { parseSvgLoops } = require('./svg');

// Geometry input, one place.  Everything that is not DXF is converted to loops here
// and written out as the native .xy; DXF keeps going through the fail-closed C++
// converter because that is where unit handling and entity diagnostics live.

const NUMBER = '[+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?';
const NUMBER_PAIR = new RegExp(`^\\s*(${NUMBER})[\\s,;]+(${NUMBER})\\s*$`);
const UNIT_TO_METRES = Object.freeze({ m: 1, mm: 0.001, cm: 0.01, in: 0.0254, ft: 0.3048 });

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
      if (index === 0 && /^\s*[a-z_][a-z_\s]*[,;\s]+[a-z_][a-z_\s]*\s*$/i.test(trimmed)) return;
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
  return ax === bx && ay === by ? loop.slice(0, -1) : loop;
}

function validateLoops(loops) {
  const issues = [];
  const cleaned = loops.map(dropClosingDuplicate);
  if (!cleaned.length) issues.push('文件里没有找到至少 3 个顶点的闭合环。');
  cleaned.forEach((loop, index) => {
    if (loop.length < 3) { issues.push(`第 ${index + 1} 个环不足三个顶点。`); return; }
    const [ox, oy] = loop[0];
    const area = loop.reduce((sum, [x, y], i) => {
      const [nx, ny] = loop[(i + 1) % loop.length];
      return sum + ((x-ox) * (ny-oy) - (nx-ox) * (y-oy));
    }, 0) / 2;
    if (!(Math.abs(area) > 0)) issues.push(`第 ${index + 1} 个环的面积为零。`);
  });
  return { loops: cleaned, issues };
}

// Explicit utility for callers that deliberately want a unit-span geometry.
// The physical import path does not call it.
function normalizeToUnitSpan(loops) {
  const span = boundsOfLoops(loops).bodySpan;
  if (!(span > 0)) return { loops, scale: 1 };
  const scale = 1 / span;
  return { loops: loops.map(loop => loop.map(([x, y]) => [x * scale, y * scale])), scale };
}

function loopsToXyText(loops, header) {
  const lines = header ? [`# ${header}`] : [];
  loops.forEach((loop, index) => {
    if (index > 0) lines.push('');
    for (const [x, y] of loop) lines.push(`${x.toPrecision(17)} ${y.toPrecision(17)}`);
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
  if (kind !== 'svg' && parsed.warnings.length) issues.push(...parsed.warnings);
  if (issues.length) return { kind, issues, warnings: parsed.warnings || [] };

  // Unit conversion is physical. Never normalise a drawing silently: the same
  // body exported as XY and CSV must reach the solver at the same size.
  const units = options.sourceUnits && options.sourceUnits !== 'auto' ? options.sourceUnits : 'm';
  const scale = UNIT_TO_METRES[units];
  if (!scale) return { kind, issues: [`未知源单位：${units}`], warnings: [] };
  const scaled = { loops: loops.map(loop => loop.map(([x,y]) => [x*scale,y*scale])), scale };
  return {
    kind,
    loops: scaled.loops,
    scale: scaled.scale,
    normalized: false,
    sourceUnits: units,
    outputUnits: 'm',
    issues: [],
    warnings: [...(parsed.warnings || []), ...((!options.sourceUnits || options.sourceUnits === 'auto') && kind !== 'xy'
      ? ['该格式未提供物理单位，当前按米读取；请在导入设置中确认源单位。'] : [])]
  };
}

function boundsOfLoops(loops) {
  let minX=Infinity, minY=Infinity, maxX=-Infinity, maxY=-Infinity;
  for (const loop of loops) for (const [x,y] of loop) {
    minX=Math.min(minX,x); minY=Math.min(minY,y);
    maxX=Math.max(maxX,x); maxY=Math.max(maxY,y);
  }
  const width=maxX-minX, height=maxY-minY;
  return { minX,minY,maxX,maxY,width,height,centreX:minX+width/2,centreY:minY+height/2,
    bodySpan: Math.max(width,height) };
}

module.exports = {
  classify,
  convertToLoops,
  parseCoordinateText,
  validateLoops,
  normalizeToUnitSpan,
  loopsToXyText,
  dropClosingDuplicate,
  boundsOfLoops
};
