'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const { FLUID_PRESETS, estimateBoundaryLayerYPlus } = require('../src/core/yplus');

test('turbulent Re=1e6 y+=1 estimate matches Schlichting flat-plate physics and FVM dy1=2*yP', () => {
  const res = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 15.0,
    nu: 1.5e-5,
    targetYPlus: 1,
    growthRatio: 1.2,
    regime: 'auto'
  });
  assert.equal(res.effectiveRegime, 'turbulent');
  assert.ok(Math.abs(res.reynolds - 1e6) < 1e-6);
  // Schlichting Cf(1e6) = (2*6 - 0.65)^-2.3 ≈ 0.003745
  assert.ok(Math.abs(res.cf - 0.003745) < 5e-5);
  assert.ok(Math.abs(res.firstLayerHeight - 2 * res.firstCellCenterDistance) < 1e-15);
  assert.ok(res.firstLayerHeight > 4e-5 && res.firstLayerHeight < 5.5e-5);
  assert.ok(res.nLayers >= 24 && res.nLayers <= 28);
  assert.ok(res.coverageAchieved >= 1.0);
  assert.ok(res.recommendedWallRelativeSize > res.firstLayerRelativeSize * 10);
  assert.equal(res.warnings.length, 0);
});

test('buffer-layer y+=12 triggers explicit advisory warning and wall-function y+=30 uses fewer layers', () => {
  const buf = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 20.0,
    nu: FLUID_PRESETS.air20.nu,
    targetYPlus: 12,
    growthRatio: 1.2
  });
  assert.ok(buf.warnings.some(w => w.includes('缓冲层')));

  const wf = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 20.0,
    nu: FLUID_PRESETS.air20.nu,
    targetYPlus: 30,
    growthRatio: 1.2
  });
  const lowRe = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 20.0,
    nu: FLUID_PRESETS.air20.nu,
    targetYPlus: 1,
    growthRatio: 1.2
  });
  assert.ok(wf.firstLayerHeight > lowRe.firstLayerHeight * 25);
  assert.ok(wf.nLayers < lowRe.nLayers);
});

test('laminar auto regime uses Blasius correlation for low Re', () => {
  const res = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 1.0,
    nu: 1e-4,
    targetYPlus: 1,
    growthRatio: 1.2,
    regime: 'auto'
  });
  assert.equal(res.effectiveRegime, 'laminar');
  assert.ok(Math.abs(res.cf - 0.664 / 100) < 1e-9);
  assert.ok(Math.abs(res.delta99 - 4.91 / 100) < 1e-9);
});

test('invalid non-positive inputs fail closed', () => {
  assert.throws(() => estimateBoundaryLayerYPlus({ referenceLength: 0, velocity: 10, nu: 1e-5 }), /参考长度/);
  assert.throws(() => estimateBoundaryLayerYPlus({ referenceLength: 1, velocity: -5, nu: 1e-5 }), /来流速度/);
  assert.throws(() => estimateBoundaryLayerYPlus({ referenceLength: 1, velocity: 10, nu: 0 }), /运动黏度/);
  assert.throws(() => estimateBoundaryLayerYPlus({ referenceLength: 1, velocity: 10, nu: 1e-5, growthRatio: 0.9 }), /增长率/);
});
