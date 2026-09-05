'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');

const { classify, convertToLoops, parseCoordinateText, loopsToXyText,
        dropClosingDuplicate } = require('../src/core/geometry');
const { parseSvgLoops } = require('../src/core/svg');

test('every advertised extension is recognised, and nothing else is', () => {
  assert.equal(classify('a.xy'), 'xy');
  assert.equal(classify('a.DXF'), 'dxf');
  assert.equal(classify('a.svg'), 'svg');
  assert.equal(classify('a.csv'), 'text');
  assert.equal(classify('a.txt'), 'text');
  assert.equal(classify('a.stl'), null);
  assert.equal(classify('a'), null);
});

test('a blank line separates loops and comments are ignored', () => {
  const { loops } = parseCoordinateText('# header\n0 0\n1 0\n1 1\n\n2 2\n3 2\n3 3\n');
  assert.equal(loops.length, 2);
  assert.deepEqual(loops[0], [[0, 0], [1, 0], [1, 1]]);
});

test('comma and space separated coordinates are the same file', () => {
  const comma = parseCoordinateText('0,0\n1,0\n1,1\n').loops;
  const space = parseCoordinateText('0 0\n1 0\n1 1\n').loops;
  assert.deepEqual(comma, space);
});

test('a CSV header row is skipped but later junk is reported', () => {
  const first = parseCoordinateText('x,y\n0,0\n1,0\n1,1\n');
  assert.equal(first.loops[0].length, 3);
  assert.equal(first.warnings.length, 0);
  const later = parseCoordinateText('0,0\n1,0\nnot a point\n1,1\n');
  assert.equal(later.warnings.length, 1);
});

// The mesher closes loops implicitly, so a repeated first vertex would become a
// zero-length wall segment.
test('a repeated closing vertex is dropped', () => {
  assert.equal(dropClosingDuplicate([[0, 0], [1, 0], [1, 1], [0, 0]]).length, 3);
  assert.equal(dropClosingDuplicate([[0, 0], [1, 0], [1, 1]]).length, 3);
});

test('degenerate input is refused rather than passed on', () => {
  assert.ok(convertToLoops('a.csv', '0,0\n1,0\n').issues.length, 'two vertices is not a loop');
  assert.ok(convertToLoops('a.csv', '0,0\n1,0\n2,0\n').issues.length, 'collinear has zero area');
  assert.ok(convertToLoops('a.stl', '').issues.length);
});

test('a unitless format is normalised to a unit body span, .xy is left alone', () => {
  const csv = convertToLoops('a.csv', '0,0\n200,0\n200,80\n0,80\n');
  assert.equal(csv.normalized, true);
  const xs = csv.loops[0].map(p => p[0]);
  assert.ok(Math.abs(Math.max(...xs) - Math.min(...xs) - 1) < 1e-12);

  const xy = convertToLoops('a.xy', '0 0\n200 0\n200 80\n0 80\n');
  assert.equal(xy.normalized, false);
  assert.equal(xy.loops[0][1][0], 200);
});

test('the .xy writer round-trips through the reader', () => {
  const loops = [[[0, 0], [1, 0], [1, 1]], [[3, 3], [4, 3], [4, 4]]];
  const parsed = parseCoordinateText(loopsToXyText(loops, 'test')).loops;
  assert.equal(parsed.length, 2);
  assert.equal(parsed[1][2][1], 4);
});

test('SVG polygon and rect become loops, y is flipped to point up', () => {
  const { loops } = parseSvgLoops('<svg><polygon points="0,0 10,0 10,4"/><rect x="0" y="0" width="2" height="2"/></svg>');
  assert.equal(loops.length, 2);
  // SVG y grows downward; after the flip the third polygon vertex sits below zero.
  assert.deepEqual(loops[0], [[0, -0], [10, -0], [10, -4]]);
});

test('an SVG circle is flattened to a chord tolerance, not a fixed segment count', () => {
  const coarse = parseSvgLoops('<svg><circle cx="0" cy="0" r="100"/></svg>', { chordToleranceFraction: 0.01 });
  const fine = parseSvgLoops('<svg><circle cx="0" cy="0" r="100"/></svg>', { chordToleranceFraction: 0.0005 });
  assert.ok(fine.loops[0].length > coarse.loops[0].length,
    'a tighter tolerance must produce more vertices');
  for (const [x, y] of fine.loops[0]) {
    assert.ok(Math.abs(Math.hypot(x, y) - 100) < 1e-9, 'vertices must lie on the circle');
  }
});

test('a cubic path is flattened and closed subpaths are kept', () => {
  const { loops } = parseSvgLoops('<svg><path d="M0,0 C 30,0 30,30 0,30 Z"/></svg>');
  assert.equal(loops.length, 1);
  assert.ok(loops[0].length > 4, 'the curve must be subdivided');
});

test('an open subpath is skipped with a warning rather than silently closed', () => {
  const { loops, warnings } = parseSvgLoops('<svg><path d="M0,0 L10,0 L10,10"/></svg>');
  assert.equal(loops.length, 0);
  assert.equal(warnings.length, 1);
});

// transform would move the boundary somewhere other than where it looks, so it is
// called out instead of ignored.
test('transform and rounded rects are reported, not quietly mishandled', () => {
  assert.ok(parseSvgLoops('<svg><g transform="translate(5,5)"><rect x="0" y="0" width="2" height="2"/></g></svg>')
    .warnings.some(w => /transform/.test(w)));
  assert.ok(parseSvgLoops('<svg><rect x="0" y="0" width="2" height="2" rx="0.4"/></svg>')
    .warnings.some(w => /圆角/.test(w)));
});
