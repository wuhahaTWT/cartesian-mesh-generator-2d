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
