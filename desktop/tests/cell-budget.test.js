'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');

const {
  BUDGET_PRESETS,
  MIN_TARGET_CELLS,
  MAX_TARGET_CELLS,
  planBudget,
  refineBudget,
  assessBudget
} = require('../src/core/cell-budget');
const { validateJob } = require('../src/core/job');

const request = extra => ({
  targetCells: 15000,
  method: 'cutcell',
  geometryPath: '/input.xy',
  fluidRegion: 'exterior',
  farFieldSpans: 0.5,
  smallAlpha: 0.15,
  ...extra
});

test('the catalog exposes the five bounded approximate-count presets', () => {
  assert.deepEqual(BUDGET_PRESETS.map(option => option.value),
    [5000, 15000, 50000, 150000, 300000]);
  assert.ok(BUDGET_PRESETS.every(option => /约/.test(option.label)));
});

test('the area seed is invariant under physical rescaling with automatic reference length', () => {
  const small = planBudget(request(), { bodySpan: 2, fluidArea: 12 });
  const large = planBudget(request(), { bodySpan: 2000, fluidArea: 12000000 });
  assert.equal(small.wallRelativeSize, large.wallRelativeSize);
  assert.equal(small.backgroundRelativeSize, large.backgroundRelativeSize);
  assert.equal(small.budgetPlan.areaSource, 'fluidArea');
  assert.equal(small.budgetPlan.areaEstimated, false);
});

test('an explicit reference changes only the dimensionless sizes, not the physical seed', () => {
  const automatic = planBudget(request(), { bodySpan: 2, fluidArea: 12 });
  // Keep the physical far-field padding fixed when changing the reference.
  const explicit = planBudget(request({ referenceLength: 4, farFieldSpans: 0.25 }),
    { bodySpan: 2, fluidArea: 12 });
  assert.equal(explicit.backgroundRelativeSize, automatic.backgroundRelativeSize / 2);
  assert.equal(explicit.wallRelativeSize, automatic.wallRelativeSize / 2);
  assert.equal(explicit.budgetPlan.domainSpan, automatic.budgetPlan.domainSpan);
  assert.equal(explicit.referenceLength, 4);
});

test('planning preserves domain, physics, method, and layer controls', () => {
  const original = request({
    targetCells: 50000,
    method: 'hybrid',
    fluidRegion: 'interior',
    farFieldSpans: 0.75,
    referenceLength: 6,
    nLayers: 5,
    growthRatio: 1.17,
    extrusionRelativeSize: 0.02,
    remainderSmallAlpha: 0.1
  });
  const planned = planBudget(original, { bodySpan: 6, fluidArea: 4.5 });
  for (const key of ['method', 'fluidRegion', 'farFieldSpans', 'referenceLength',
    'nLayers', 'growthRatio']) assert.equal(planned[key], original[key]);
  assert.equal(planned.cellsPerLevel, 6);
  assert.equal(planned.firstLayerRelativeSize, planned.wallRelativeSize / 4);
  assert.equal(validateJob(planned).job.fluidRegion, 'interior');
});

test('planning respects the explicit wall ratio and transition-band width', () => {
  const planned = planBudget(request({
    targetCells: 15000,
    wallToBackgroundRatio: 8,
    cellsPerLevel: 12,
    allowUnsafeWallLevel: true
  }), { bodySpan: 1, fluidArea: 1 });
  assert.equal(planned.backgroundRelativeSize / planned.wallRelativeSize, 8);
  assert.equal(planned.cellsPerLevel, 12);
  assert.equal(planned.budgetPlan.wallToBackgroundRatio, 8);
});

test('planning never grants the unsafe wall-depth override', () => {
  const absent = planBudget(request({ targetCells: MAX_TARGET_CELLS, method: 'hybrid' }),
    { bodySpan: 1, fluidArea: 1 });
  assert.equal(Object.hasOwn(absent, 'allowUnsafeWallLevel'), false);
  assert.equal(absent.budgetPlan.limitedBySafeWallLevel, true);

  const refused = planBudget(request({ targetCells: MAX_TARGET_CELLS, method: 'hybrid',
    allowUnsafeWallLevel: false }), { bodySpan: 1, fluidArea: 1 });
  assert.equal(refused.allowUnsafeWallLevel, false);
  assert.ok(refused.wallRelativeSize >= refused.budgetPlan.minimumWallRelativeSize);
  assert.equal(refused.backgroundRelativeSize / refused.wallRelativeSize, 4);

  const explicit = planBudget(request({ targetCells: MAX_TARGET_CELLS, method: 'hybrid',
    allowUnsafeWallLevel: true }), { bodySpan: 1, fluidArea: 1 });
  assert.equal(explicit.allowUnsafeWallLevel, true);
  assert.equal(explicit.budgetPlan.limitedBySafeWallLevel, false);
});

test('fallback area is clearly identified as a conservative bbox estimate', () => {
  const planned = planBudget(request(), {
    bodySpan: 2,
    bounds: { minX: -1, maxX: 1, minY: -0.25, maxY: 0.25 }
  });
  assert.equal(planned.budgetPlan.areaSource, 'boundingBoxEstimate');
  assert.equal(planned.budgetPlan.areaEstimated, true);
  assert.equal(planned.budgetPlan.area, 16, 'exterior native domain is a padded square');
});

test('target and frame inputs fail closed outside their bounds', () => {
  assert.throws(() => planBudget(request({ targetCells: MIN_TARGET_CELLS - 1 }),
    { bodySpan: 1 }), /目标单元数/);
  assert.throws(() => planBudget(request({ targetCells: MAX_TARGET_CELLS + 1 }),
    { bodySpan: 1 }), /目标单元数/);
  assert.throws(() => planBudget(request({ targetCells: 1500.5 }), { bodySpan: 1 }), /目标单元数/);
  assert.throws(() => planBudget(request(), { bodySpan: 0 }), /几何体跨度/);
  assert.throws(() => planBudget(request(), { bodySpan: 1, fluidArea: NaN }), /流体面积/);
  assert.throws(() => planBudget(request({ fluidRegion: 'solid' }), { bodySpan: 1 }), /流体区域/);
  assert.throws(() => planBudget(request({ method: 'hybrid', farFieldSpans: 1000 }),
    { bodySpan: 1 }), /安全壁面层级/);
});

test('assessment reports an honest 70 to 130 percent target interval', () => {
  assert.equal(assessBudget(7000, 10000).withinRange, true);
  assert.equal(assessBudget(13000, 10000).withinRange, true);
  assert.equal(assessBudget(6999, 10000).withinRange, false);
  assert.equal(assessBudget(13001, 10000).withinRange, false);
});

test('observed count correction scales wall, background, and first layer together', () => {
  const planned = planBudget(request({
    targetCells: 50000,
    method: 'hybrid',
    nLayers: 4,
    growthRatio: 1.2,
    extrusionRelativeSize: 0.01,
    cellsPerLevel: 0,
    remainderSmallAlpha: 0.1,
    allowUnsafeWallLevel: true
  }), { bodySpan: 1, fluidArea: 1 });
  const refined = refineBudget(planned, 200000, 50000);
  assert.equal(refined.budgetPlan.correctionFactor, 1.6, 'one correction is bounded');
  assert.equal(refined.wallRelativeSize / planned.wallRelativeSize, 1.6);
  assert.equal(refined.backgroundRelativeSize / planned.backgroundRelativeSize, 1.6);
  assert.equal(refined.firstLayerRelativeSize / planned.firstLayerRelativeSize, 1.6);
  assert.equal(refined.farFieldSpans, planned.farFieldSpans);
  assert.equal(refined.nLayers, planned.nLayers);
  assert.equal(refined.growthRatio, planned.growthRatio);
});

test('a repeated out-of-range count is marked as a plateau and forces a useful step', () => {
  const planned = planBudget(request({ targetCells: 50000, allowUnsafeWallLevel: true,
    cellsPerLevel: 0 }),
    { bodySpan: 1, fluidArea: 1 });
  const once = refineBudget(planned, 20000);
  const twice = refineBudget(once, 20100);
  assert.equal(twice.budgetPlan.plateau, true);
  assert.ok(twice.budgetPlan.correctionFactor <= 0.8);
  const exhausted = refineBudget(twice, 20200);
  assert.equal(exhausted.budgetPlan.attempt, 3);
  assert.equal(exhausted.budgetPlan.exhausted, true);
});

test('the quantized initial uniform floor stays within half the requested budget', () => {
  for (const targetCells of [5000, 15000, 50000, 150000, 300000]) {
    const planned = planBudget(request({ targetCells }),
      { bodySpan: 2, fluidArea: 12.8785 });
    const backgroundSize = planned.backgroundRelativeSize * planned.budgetPlan.referenceLength;
    const uniformFloor = planned.budgetPlan.area / (backgroundSize * backgroundSize);
    assert.ok(uniformFloor <= targetCells / 2,
      `${targetCells}: quantized floor ${uniformFloor} exceeds half-target`);
  }
});

test('the circle 15k plateau narrows the transition band and preserves the physical plan', () => {
  const planned = planBudget(request({ targetCells: 15000, referenceLength: 2,
    cellsPerLevel: 6 }), { bodySpan: 2, fluidArea: 12.8785 });
  // Reproduce the prior quantized seed from the real circle run.
  const oldSeed = {
    ...planned,
    wallRelativeSize: 0.004095,
    backgroundRelativeSize: 0.01638,
    cellsPerLevel: 6,
    budgetPlan: {
      ...planned.budgetPlan,
      area: 12.8785,
      domainSpan: 4,
      referenceLength: 2
    }
  };
  const refined = refineBudget(oldSeed, 21572);
  assert.equal(refined.budgetPlan.correctionKind, 'band-width');
  assert.ok(refined.cellsPerLevel < oldSeed.cellsPerLevel);
  assert.equal(refined.wallRelativeSize, oldSeed.wallRelativeSize);
  assert.equal(refined.backgroundRelativeSize, oldSeed.backgroundRelativeSize);
  assert.equal(refined.budgetPlan.domainSpan, 4);
  assert.equal(refined.budgetPlan.area, 12.8785);
  assert.equal(refined.budgetPlan.targetCells, 15000);
  assert.equal(refined.budgetPlan.minimumCells, oldSeed.budgetPlan.minimumCells);
  assert.equal(refined.budgetPlan.maximumCells, oldSeed.budgetPlan.maximumCells);
});
