'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const {
  FLUID_PRESETS,
  WAKE_REGION_SOURCE,
  buildRecommendedWakeBox,
  upsertYPlusWakeRegion,
  estimateBoundaryLayerYPlus
} = require('../src/core/yplus');
const { validateJob, buildInvocation } = require('../src/core/job');

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
  assert.throws(() => buildRecommendedWakeBox({ referenceLength: 1, wakeLengthSpans: 0 }), /尾迹长度/);
  assert.throws(() => buildRecommendedWakeBox({ referenceLength: 1, wakeAngleDeg: 60 }), /尾迹偏角/);
});

test('exterior high-Re y+ estimate builds aligned wake box, respects deflection/clamping, and skips interior flow', () => {
  const frame = { minX: 0, maxX: 1, minY: -0.06, maxY: 0.06, width: 1.0, height: 0.12, centreX: 0.5, centreY: 0.0, bodySpan: 1.0 };
  const ext = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 90.0,
    nu: 1.5e-5,
    targetYPlus: 1,
    growthRatio: 1.18,
    coverageRatio: 1.5,
    fluidRegion: 'exterior',
    frame,
    wakeLengthSpans: 2.0,
    wakeLevelsBelowWall: 2
  });
  assert.ok(ext.recommendedWakeBox);
  assert.equal(ext.recommendedWakeBox.source, WAKE_REGION_SOURCE);
  assert.equal(ext.recommendedWakeBox.levelsBelowWall, 2);
  // Upstream covers leading edge (xmin < -0.5), downstream extends 2 chords past TE (xmax = 2.5)
  assert.ok(ext.recommendedWakeBox.xmin < -0.55 && ext.recommendedWakeBox.xmin >= -0.75);
  assert.equal(ext.recommendedWakeBox.xmax, 2.5);
  assert.ok(ext.recommendedWakeBox.ymax >= 0.35 && ext.recommendedWakeBox.ymin <= -0.35);
  assert.ok(ext.recommendedFarFieldSpans >= 15);

  // Deflected wake angle widens vertical envelope toward wake direction
  const deflected = buildRecommendedWakeBox({
    referenceLength: 1.0,
    frame,
    wakeLengthSpans: 2.0,
    wakeAngleDeg: -10,
    wakeLevelsBelowWall: 2
  });
  assert.ok(deflected.ymin < ext.recommendedWakeBox.ymin - 0.2);

  // Tight far-field clamps wake box strictly inside domain
  const clamped = buildRecommendedWakeBox({
    referenceLength: 1.0,
    frame,
    farFieldSpans: 1.0,
    wakeLengthSpans: 5.0
  });
  assert.ok(clamped.xmax <= 0.5 + 0.85 * 1.0 + 1e-9);

  // Interior flow must not emit an exterior wake box
  const interior = estimateBoundaryLayerYPlus({
    referenceLength: 1.0,
    velocity: 15.0,
    nu: 1.5e-5,
    fluidRegion: 'interior',
    frame
  });
  assert.equal(interior.recommendedWakeBox, null);
  assert.equal(interior.recommendedFarFieldSpans, null);
});

test('upsertYPlusWakeRegion replaces prior wake box idempotently and wires valid --refine-box CLI args', () => {
  const manualBox = { xmin: -0.2, xmax: 0.3, ymin: -0.1, ymax: 0.1, levelsBelowWall: 1, active: false };
  const wake1 = { xmin: -0.65, xmax: 2.5, ymin: -0.35, ymax: 0.35, levelsBelowWall: 2, source: WAKE_REGION_SOURCE };
  const wake2 = { xmin: -0.65, xmax: 3.5, ymin: -0.4, ymax: 0.4, levelsBelowWall: 3, source: WAKE_REGION_SOURCE };

  const afterFirst = upsertYPlusWakeRegion([manualBox], wake1);
  assert.equal(afterFirst.length, 2);
  assert.equal(afterFirst[0].xmin, -0.2);
  assert.equal(afterFirst[1].source, WAKE_REGION_SOURCE);
  assert.equal(afterFirst[1].active, true);

  const afterSecond = upsertYPlusWakeRegion(afterFirst, wake2);
  assert.equal(afterSecond.length, 2);
  assert.equal(afterSecond[1].xmax, 3.5);
  assert.equal(afterSecond[1].levelsBelowWall, 3);

  // End-to-end job validation + hybrid CLI invocation emits --refine-box accurately
  const frame = { minX: 0, maxX: 1, minY: -0.06, maxY: 0.06, width: 1.0, height: 0.12, centreX: 0.5, centreY: 0.0, bodySpan: 1.0 };
  const { job } = validateJob({
    method: 'hybrid',
    fluidRegion: 'exterior',
    geometryPath: 'airfoil.xy',
    outputDirectory: 'out',
    sizingMode: 'relative',
    referenceMode: 'explicit',
    referenceLength: 1.0,
    wallRelativeSize: 0.012,
    backgroundRelativeSize: 0.6,
    farFieldSpans: 15,
    cellsPerLevel: 3,
    allowUnsafeWallLevel: true,
    nLayers: 30,
    firstLayerRelativeSize: 9e-6,
    growthRatio: 1.2,
    extrusionRelativeSize: 0.05,
    refineBoxes: [wake2]
  });
  const inv = buildInvocation(job, {
    xyPath: 'airfoil.xy',
    prefix: 'out/mesh',
    casePath: 'out/foam',
    frame
  });
  assert.equal(inv.executable, 'cartmesh2d_hybrid_cli');
  const idx = inv.args.indexOf('--refine-box');
  assert.ok(idx > 0);
  // centreX=0.5 + xmin(-0.65) = -0.15; centreX=0.5 + xmax(3.5) = 4.0
  assert.ok(Math.abs(Number(inv.args[idx + 1]) - (-0.15)) < 1e-9);
  assert.ok(Math.abs(Number(inv.args[idx + 2]) - (-0.4)) < 1e-9);
  assert.ok(Math.abs(Number(inv.args[idx + 3]) - 4.0) < 1e-9);
  assert.ok(Math.abs(Number(inv.args[idx + 4]) - 0.4) < 1e-9);
  assert.ok(Number(inv.args[idx + 5]) >= 1);
});
