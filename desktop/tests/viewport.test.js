'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');

test('sidebar and window resize round trips preserve zoom and world centre', () => {
  const sandbox = { window: { devicePixelRatio: 2 }, ResizeObserver: class { observe() {} } };
  vm.createContext(sandbox);
  vm.runInContext(fs.readFileSync(path.join(__dirname, '../src/renderer/viewport.js'), 'utf8'), sandbox);
  let size = { width: 1000, height: 600 };
  const canvas = { addEventListener() {}, getBoundingClientRect: () => size,
    getContext: () => ({ setTransform() {}, clearRect() {} }) };
  const view = new sandbox.window.MeshView.Viewport(canvas);
  view.fitTo({ minX: -3, maxX: 3, minY: -1, maxY: 1 }, .1);
  const initial = { scale: view.scale, x: view.offset.x, y: view.offset.y };
  for (let i = 0; i < 5; i++) {
    for (const next of [{ width: 1336, height: 600 }, { width: 1000, height: 600 },
                        { width: 500, height: 232 }, { width: 1000, height: 600 }]) {
      size = next;
      view.context();
    }
  }
  assert.ok(Math.abs(view.scale - initial.scale) < 1e-10);
  assert.ok(Math.abs(view.offset.x - initial.x) < 1e-10);
  assert.ok(Math.abs(view.offset.y - initial.y) < 1e-10);
});

function viewportSandbox() {
  let nextFrame = 1;
  const frames = new Map();
  class FakePath2D {
    constructor() { this.commands = []; }
    moveTo(x, y) { this.commands.push(['M', x, y]); }
    lineTo(x, y) { this.commands.push(['L', x, y]); }
    closePath() { this.commands.push(['Z']); }
  }
  const sandbox = {
    Path2D: FakePath2D,
    ResizeObserver: class { observe() {} },
    window: {
      devicePixelRatio: 2,
      requestAnimationFrame(callback) {
        const id = nextFrame++;
        frames.set(id, callback);
        return id;
      },
      cancelAnimationFrame(id) { frames.delete(id); }
    }
  };
  vm.createContext(sandbox);
  vm.runInContext(fs.readFileSync(path.join(__dirname, '../src/renderer/viewport.js'), 'utf8'), sandbox);
  return { sandbox, frames };
}

function fakeCanvas(listeners = {}) {
  const draws = { fill: [], stroke: [], transforms: [] };
  const context = {
    setTransform() {}, clearRect() {}, save() {}, restore() {},
    transform(...values) { draws.transforms.push(values); },
    beginPath() {}, moveTo() {}, lineTo() {}, closePath() {}, fillRect() {},
    fill(path) { draws.fill.push(path); },
    stroke(path) { draws.stroke.push(path); },
    setLineDash() {}, strokeRect() {}, fillText() {}
  };
  return {
    canvas: {
      classList: { add() {}, remove() {} },
      addEventListener(name, callback) { listeners[name] = callback; },
      setPointerCapture() {},
      getBoundingClientRect: () => ({ width: 800, height: 500, left: 0, top: 0 }),
      getContext: () => context
    },
    draws
  };
}

test('mesh paths and boundary groups are built once and reused across view changes', () => {
  const { sandbox } = viewportSandbox();
  const { canvas, draws } = fakeCanvas();
  let vertexReads = 0;
  let cellIterations = 0;
  let edgeIterations = 0;
  const observed = (values, onRead) => new Proxy(values, {
    get(target, key, receiver) {
      if (key === Symbol.iterator) onRead();
      else if (typeof key === 'string' && /^\d+$/.test(key)) vertexReads++;
      return Reflect.get(target, key, receiver);
    }
  });
  const mesh = {
    bounds: { minX: 0, minY: 0, maxX: 2, maxY: 1 }, minLevel: 1, maxLevel: 2,
    vertices: observed([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]], () => {}),
    cells: observed([
      { level: 1, vertices: [0, 1, 2, 3] },
      { level: 2, vertices: [1, 4, 5, 2] }
    ], () => { cellIterations++; }),
    edges: observed([
      { patch: 2, a: 0, b: 1 }, { patch: 3, a: 1, b: 4 }, { patch: 1, a: 2, b: 3 }
    ], () => { edgeIterations++; })
  };
  const view = new sandbox.window.MeshView.Viewport(canvas);
  view.setMesh(mesh);
  const firstVertexReads = vertexReads;
  assert.equal(cellIterations, 1);
  assert.equal(edgeIterations, 1);
  assert.ok(firstVertexReads > 0);
  assert.equal(view.meshCache.levels.length, 2);
  assert.equal(draws.fill.length, 2);
  assert.equal(draws.stroke.length, 5);

  view.scale *= 2;
  view.draw();
  assert.equal(cellIterations, 1);
  assert.equal(edgeIterations, 1);
  assert.equal(vertexReads, firstVertexReads);
  assert.equal(draws.fill.length, 4);
  assert.equal(draws.stroke.length, 10);
  assert.deepEqual(draws.transforms.at(-1), [view.scale, 0, 0, -view.scale,
    -view.offset.x * view.scale, 500 + view.offset.y * view.scale]);

  view.setOutline([[[0, 0], [1, 0], [0, 1]]]);
  assert.equal(view.meshCache, null);
  view.clear();
  assert.equal(view.meshCache, null);
});

test('interactive redraw requests are merged into one animation frame', () => {
  const { sandbox, frames } = viewportSandbox();
  const { canvas } = fakeCanvas();
  const view = new sandbox.window.MeshView.Viewport(canvas);
  let draws = 0;
  const originalDraw = view.draw.bind(view);
  view.draw = () => { draws++; originalDraw(); };

  view.requestDraw();
  view.requestDraw();
  view.requestDraw();
  assert.equal(frames.size, 1);
  const callback = [...frames.values()][0];
  frames.clear();
  callback();
  assert.equal(draws, 1);
  assert.equal(view.pendingDraw, null);

  view.requestDraw();
  assert.equal(frames.size, 1);
  view.draw();
  assert.equal(frames.size, 0);
  assert.equal(draws, 2);
});

test('unchanged view size does not reallocate the canvas backing store', () => {
  const { sandbox } = viewportSandbox();
  let size = { width: 800, height: 500 };
  let width = 0, height = 0, allocations = 0;
  const canvas = {
    addEventListener() {},
    getBoundingClientRect: () => size,
    getContext: () => ({ setTransform() {}, clearRect() {} }),
    get width() { return width; },
    set width(value) { width = value; allocations++; },
    get height() { return height; },
    set height(value) { height = value; allocations++; }
  };
  const view = new sandbox.window.MeshView.Viewport(canvas);
  view.draw();
  view.draw();
  assert.equal(allocations, 2);
  size = { width: 700, height: 500 };
  view.draw();
  assert.equal(allocations, 3);
});

test('zoomed view submits only intersecting cell path chunks', () => {
  const { sandbox } = viewportSandbox();
  const { canvas, draws } = fakeCanvas();
  const cells = Array.from({ length: 64 }, () => ({ level: 0, vertices: [4, 5, 6, 7] }));
  cells.push({ level: 0, vertices: [0, 1, 2, 3] });
  const mesh = {
    bounds: { minX: 0, minY: 0, maxX: 101, maxY: 1 }, minLevel: 0, maxLevel: 0,
    vertices: [[0, 0], [1, 0], [1, 1], [0, 1],
      [100, 0], [101, 0], [101, 1], [100, 1]],
    cells,
    edges: []
  };
  const view = new sandbox.window.MeshView.Viewport(canvas);
  view.setMesh(mesh);
  assert.equal(view.meshCache.levels[0].chunks.length, 2);
  const fillsAtFit = draws.fill.length;
  view.scale = 10;
  view.offset = { x: 0, y: 0 };
  view.draw();
  assert.equal(draws.fill.length - fillsAtFit, 1);
});
