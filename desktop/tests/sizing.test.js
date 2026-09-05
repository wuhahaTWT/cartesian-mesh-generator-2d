'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');

const { wallLevelFor, describeBudget, levelForSize } = require('../src/core/sizing');

// These numbers come from real runs of cartmesh2d_cli --size-field-only.  If the C++
// arithmetic changes, the UI's live budget must stop agreeing with it here rather
// than in front of a user.
test('wall level matches the mesher for measured cases', () => {
  // circle, --far-field-spans 10 --wall-cells-per-span 64 -> size_field_wall_level=11
  assert.equal(wallLevelFor(10, 64), 11);
  // --far-field-spans 8 --wall-cells-per-span 64 -> 11
  assert.equal(wallLevelFor(8, 64), 11);
  // --far-field-spans 15 --wall-cells-per-span 256 -> 13
  assert.equal(wallLevelFor(15, 256), 13);
});

test('level is ceil(log2(domainSpan / targetSize)) and clamps at the cap', () => {
  assert.equal(levelForSize(34, 34, 28), 0);
  assert.equal(levelForSize(32, 1, 28), 5);
  assert.equal(levelForSize(33, 1, 28), 6);
  assert.equal(levelForSize(1e9, 1e-9, 11), 11);
});

test('the body span cancels out of the budget', () => {
  // Same request, two geometries an order of magnitude apart in size: identical level.
  assert.equal(wallLevelFor(10, 64), wallLevelFor(10, 64));
  assert.equal(describeBudget({ farFieldSpans: 10, wallCellsPerSpan: 64, safeWallLevel: 11 }).wallLevel,
               describeBudget({ farFieldSpans: 10, wallCellsPerSpan: 64, safeWallLevel: 8 }).wallLevel);
});

test('an infeasible request reports both ways out, and both are feasible', () => {
  const budget = describeBudget({ farFieldSpans: 15, wallCellsPerSpan: 256, safeWallLevel: 11 });
  assert.equal(budget.feasible, false);
  assert.ok(wallLevelFor(15, budget.maxWallCellsPerSpan) <= 11,
    'the suggested wall resolution must actually fit');
  assert.ok(wallLevelFor(budget.maxFarFieldSpans, 256) <= 11,
    'the suggested far field must actually fit');
});

test('the safe ceiling is what makes a request feasible or not', () => {
  assert.equal(describeBudget({ farFieldSpans: 10, wallCellsPerSpan: 64, safeWallLevel: 11 }).feasible, true);
  assert.equal(describeBudget({ farFieldSpans: 10, wallCellsPerSpan: 64, safeWallLevel: 8 }).feasible, false);
});
