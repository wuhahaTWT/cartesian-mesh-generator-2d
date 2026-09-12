'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const { extractContours } = require('../src/core/raster');

function image(width, height, rgba = [255, 255, 255, 255]) {
  const data = new Uint8ClampedArray(width * height * 4);
  for (let i = 0; i < width * height; i++) data.set(rgba, i * 4);
  return { width, height, data };
}

function paint(result, x0, y0, x1, y1, rgba = [0, 0, 0, 255]) {
  for (let y = y0; y < y1; y++) for (let x = x0; x < x1; x++)
    result.data.set(rgba, (y * result.width + x) * 4);
  return result;
}

function area(loop) {
  return loop.reduce((sum, point, index) => {
    const next = loop[(index + 1) % loop.length];
    return sum + point[0] * next[1] - next[0] * point[1];
  }, 0) / 2;
}

function totalArea(loops) { return loops.reduce((sum, loop) => sum + area(loop), 0); }

function pointSegmentDistance(point, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1];
  const denominator = dx * dx + dy * dy;
  const t = denominator ? Math.max(0, Math.min(1,
    ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / denominator)) : 0;
  return Math.hypot(point[0] - a[0] - t * dx, point[1] - a[1] - t * dy);
}

function distanceToLoop(point, loop) {
  return Math.min(...loop.map((a, index) => pointSegmentDistance(point, a, loop[(index + 1) % loop.length])));
}

function orientation(a, b, c) {
  const cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  return Math.sign(cross);
}

function selfIntersects(loop) {
  const on = (a, b, p) => orientation(a, b, p) === 0 && p[0] >= Math.min(a[0], b[0]) &&
    p[0] <= Math.max(a[0], b[0]) && p[1] >= Math.min(a[1], b[1]) && p[1] <= Math.max(a[1], b[1]);
  const intersects = (a, b, c, d) => {
    const o1 = orientation(a, b, c), o2 = orientation(a, b, d);
    const o3 = orientation(c, d, a), o4 = orientation(c, d, b);
    return (o1 !== o2 && o3 !== o4) || (!o1 && on(a, b, c)) || (!o2 && on(a, b, d)) ||
      (!o3 && on(c, d, a)) || (!o4 && on(c, d, b));
  };
  for (let i = 0; i < loop.length; i++) for (let j = i + 1; j < loop.length; j++) {
    if (j === i + 1 || (i === 0 && j === loop.length - 1)) continue;
    if (intersects(loop[i], loop[(i + 1) % loop.length], loop[j], loop[(j + 1) % loop.length])) return true;
  }
  return false;
}

test('a filled rectangle becomes one exact pixel-boundary loop', () => {
  const result = extractContours(paint(image(12, 10), 2, 3, 9, 8), { threshold: null });
  assert.equal(result.mode, 'dark');
  assert.equal(result.threshold, 0);
  assert.equal(result.loops.length, 1);
  assert.equal(result.loops[0].length, 4);
  assert.equal(area(result.loops[0]), 35, 'y-down outer loops have positive area');
  assert.deepEqual(result.loops[0], [[2, 3], [9, 3], [9, 8], [2, 8]]);
});

test('a concave silhouette keeps its notch and exact area', () => {
  const source = image(14, 14);
  paint(source, 2, 2, 11, 11);
  paint(source, 7, 5, 11, 8, [255, 255, 255, 255]);
  const result = extractContours(source);
  assert.equal(result.loops.length, 1);
  assert.equal(totalArea(result.loops), 69);
  assert.ok(result.loops[0].some(([x, y]) => x === 7 && y === 5), 'the inner corner survives');
  assert.equal(selfIntersects(result.loops[0]), false);
});

test('holes are preserved by default and filled only when explicitly requested', () => {
  const source = image(14, 14);
  paint(source, 2, 2, 12, 12);
  paint(source, 5, 5, 9, 9, [255, 255, 255, 255]);
  const preserved = extractContours(source);
  assert.deepEqual(preserved.loops.map(area).sort((a, b) => a - b), [-16, 100]);
  assert.equal(preserved.stats.filled_holes, 0);
  assert.equal(totalArea(preserved.loops), 84);

  const filled = extractContours(source, { fillHoles: true });
  assert.equal(filled.loops.length, 1);
  assert.equal(totalArea(filled.loops), 100);
  assert.equal(filled.stats.filled_holes, 1);
  assert.equal(filled.stats.filled_area_px, 16);
  assert.ok(filled.warnings.some(warning => /显式设置填充 1 个内孔/.test(warning)));
});

test('largest, all, and a seed select connected components deterministically', () => {
  const source = image(24, 14);
  paint(source, 2, 2, 6, 6);
  paint(source, 12, 3, 21, 11);
  const largest = extractContours(source);
  assert.equal(totalArea(largest.loops), 72);
  assert.ok(largest.warnings.some(warning => /最大主体/.test(warning)));
  const all = extractContours(source, { selection: 'all' });
  assert.equal(all.loops.length, 2);
  assert.equal(totalArea(all.loops), 88);
  const seeded = extractContours(source, { seed: { x: 3.8, y: 3.2 } });
  assert.equal(totalArea(seeded.loops), 16);
  assert.throws(() => extractContours(source, { seed: { x: 8, y: 8 } }), /背景/);
});

test('small connected noise is removed with an explicit count and area', () => {
  const source = paint(image(30, 20), 5, 5, 20, 15);
  paint(source, 25, 4, 26, 5);
  paint(source, 24, 12, 25, 13);
  const result = extractContours(source, { selection: 'all', minAreaPx: 4 });
  assert.equal(result.loops.length, 1);
  assert.equal(result.stats.removed_components, 2);
  assert.equal(result.stats.removed_area_px, 2);
  assert.ok(result.warnings.some(warning => /2 个小噪点主体，共 2 像素/.test(warning)));
});

test('transparent coloured foreground is selected by alpha in auto mode', () => {
  const source = image(12, 12, [0, 0, 0, 0]);
  paint(source, 3, 2, 10, 9, [20, 190, 90, 255]);
  const result = extractContours(source);
  assert.equal(result.mode, 'alpha');
  assert.equal(result.threshold, 0);
  assert.equal(totalArea(result.loops), 49);
});

test('light-on-dark and border-colour distance modes select the visible object', () => {
  const lightSource = image(12, 12, [10, 10, 10, 255]);
  paint(lightSource, 3, 4, 9, 10, [245, 245, 245, 255]);
  const light = extractContours(lightSource);
  assert.equal(light.mode, 'light');
  assert.equal(totalArea(light.loops), 36);

  const colourSource = image(12, 12, [40, 140, 220, 255]);
  paint(colourSource, 2, 3, 8, 9, [220, 60, 80, 255]);
  const colour = extractContours(colourSource, { mode: 'background', threshold: 30 });
  assert.equal(totalArea(colour.loops), 36);
  assert.deepEqual(colour.stats.background_rgb, [40, 140, 220]);
  assert.ok(colour.warnings.some(warning => /不是任意场景的语义分割/.test(warning)));
});

test('default simplification stays within one pixel and bounds area drift', () => {
  const source = image(50, 40);
  for (let y = 5; y < 32; y++) paint(source, 5, y, 25 + Math.floor(y / 3), y + 1);
  const exact = extractContours(source, { simplifyTolerance: 0 });
  const simplified = extractContours(source, { simplifyTolerance: 1 });
  assert.ok(simplified.stats.points < exact.stats.points);
  for (const point of exact.loops[0])
    assert.ok(distanceToLoop(point, simplified.loops[0]) <= 1 + 1e-9);
  assert.ok(Math.abs(totalArea(simplified.loops) - totalArea(exact.loops)) <=
    Math.max(1, totalArea(exact.loops) * 0.02));
  assert.equal(selfIntersects(simplified.loops[0]), false);
});

test('diagonal-only contact remains two simple degree-four contours', () => {
  const source = image(9, 9);
  paint(source, 2, 2, 3, 3);
  paint(source, 3, 3, 4, 4);
  const result = extractContours(source, { selection: 'all', minAreaPx: 1, simplifyTolerance: 0 });
  assert.equal(result.loops.length, 2);
  assert.deepEqual(result.loops.map(area), [1, 1]);
  assert.ok(result.loops.every(loop => !selfIntersects(loop)));
  assert.deepEqual(result.loops[0][2], result.loops[1][0], 'the degree-four vertex is shared, not spliced');
});

test('the same bytes and options produce byte-for-byte deterministic output', () => {
  const source = image(22, 18);
  paint(source, 3, 3, 17, 14);
  paint(source, 8, 7, 12, 11, [255, 255, 255, 255]);
  const options = { selection: 'all', threshold: null, simplifyTolerance: 1 };
  assert.equal(JSON.stringify(extractContours(source, options)), JSON.stringify(extractContours(source, options)));
});

test('empty, full, border-cut, and over-complex masks fail explicitly', () => {
  assert.throws(() => extractContours(image(10, 10)), /结果为空/);
  assert.throws(() => extractContours(image(10, 10, [0, 0, 0, 255]), { mode: 'dark' }), /覆盖整张图片/);
  assert.throws(() => extractContours(image(10, 10), { mode: 'dark', threshold: 255 }), /覆盖整张图片/,
    'dark includes pixels equal to the threshold, including the 255 endpoint');
  const border = paint(image(10, 10), 0, 2, 5, 7);
  assert.throws(() => extractContours(border), /接触图片边缘/);

  const noisy = image(50, 50);
  for (let i = 0; i < 65; i++) {
    const x = 2 + (i % 13) * 3, y = 2 + Math.floor(i / 13) * 3;
    paint(noisy, x, y, x + 1, y + 1);
  }
  assert.throws(() => extractContours(noisy,
    { selection: 'all', minAreaPx: 1, simplifyTolerance: 0 }), /超过 64 个闭合环/);
});

test('invalid dimensions, bytes, threshold, and excessive images are bounded', () => {
  assert.throws(() => extractContours({ width: 4, height: 4, data: new Uint8ClampedArray(3) }), /字节长度/);
  assert.throws(() => extractContours(image(4, 4), { threshold: 256 }), /0 到 255/);
  assert.throws(() => extractContours({ width: 2000, height: 2000, data: [] }), /图片过大/);
});

test('the browser-worker build publishes self.CartMeshRaster', () => {
  const source = fs.readFileSync(path.join(__dirname, '../src/core/raster.js'), 'utf8');
  const context = { self: {} };
  vm.runInNewContext(source, context);
  assert.equal(typeof context.self.CartMeshRaster.extractContours, 'function');
});
