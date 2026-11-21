'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const {
  validateAdaptiveRequest,
  buildAdaptiveInvocation,
  parseAdaptiveSuggestion
} = require('../src/core/adaptive-refinement');

test('validateAdaptiveRequest applies defaults and rejects invalid bounds', () => {
  assert.deepEqual(validateAdaptiveRequest(), {
    refineFraction: 0.15,
    maxBoxes: 6,
    padCells: 1.0
  });
  assert.throws(() => validateAdaptiveRequest({ refineFraction: 0 }), /比例/);
  assert.throws(() => validateAdaptiveRequest({ refineFraction: 0.6 }), /比例/);
  assert.throws(() => validateAdaptiveRequest({ maxBoxes: 0 }), /数量/);
  assert.throws(() => validateAdaptiveRequest({ padCells: -1 }), /外扩/);
  assert.throws(() => validateAdaptiveRequest({ wallLevel: 29 }), /层级/);
});

test('buildAdaptiveInvocation constructs deterministic CLI flags', () => {
  const inv = buildAdaptiveInvocation('m.solver.cm2d', 'm.flow.cells.csv', 'm.adaptive.json', {
    refineFraction: 0.2,
    maxBoxes: 4,
    padCells: 1.5,
    wallLevel: 8,
    maxTargetLevel: 8
  });
  assert.equal(inv.executable, 'cartmesh2d_flow_cli');
  assert.deepEqual(inv.args, [
    '--mesh', 'm.solver.cm2d',
    '--suggest-refinement-from', 'm.flow.cells.csv',
    '--suggest-refinement', 'm.adaptive.json',
    '--adapt-fraction', '0.2',
    '--adapt-max-boxes', '4',
    '--adapt-pad-cells', '1.5',
    '--adapt-max-level', '8',
    '--adapt-wall-level', '8'
  ]);
});

test('parseAdaptiveSuggestion maps world boxes into active body-frame regions', () => {
  const payload = {
    format: 'cartmesh2d-adaptive-refinement-v1',
    valid: true,
    cell_count: 1000,
    selected_cell_count: 150,
    selected_area_fraction: 0.08,
    inferred_wall_level: 8,
    max_target_level: 8,
    indicator_max: 0.42,
    boxes: [{
      xmin: 1.5, ymin: -0.5, xmax: 3.5, ymax: 0.5,
      target_level: 7, levels_below_wall: 1,
      body_xmin: 0.5, body_ymin: -0.5, body_xmax: 2.5, body_ymax: 0.5,
      cell_count: 80, mean_indicator: 0.12, max_indicator: 0.42
    }]
  };
  const parsed = parseAdaptiveSuggestion(payload, { cx: 0.5, cy: 0, span: 2.0 });
  assert.equal(parsed.regions.length, 1);
  assert.deepEqual(parsed.regions[0], {
    xmin: 0.5,
    xmax: 1.5,
    ymin: -0.25,
    ymax: 0.25,
    levelsBelowWall: 1,
    active: true,
    targetLevel: 7,
    cellCount: 80,
    meanIndicator: 0.12,
    maxIndicator: 0.42
  });
  assert.throws(() => parseAdaptiveSuggestion({ format: 'bad', valid: true, boxes: [] }), /未知/);
  assert.throws(() => parseAdaptiveSuggestion({
    format: 'cartmesh2d-adaptive-refinement-v1', valid: false, issues: ['NaN'], boxes: []
  }), /NaN/);
});
