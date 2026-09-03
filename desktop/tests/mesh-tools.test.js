'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const {
  getCapabilities,
  validateMeshJob,
  buildMeshInvocation,
  normalizeSummary
} = require('../src/mesh-tools');

const common = {
  dxfPath: '/tmp/body.dxf',
  outputDirectory: '/tmp/out',
  sourceUnits: 'm',
  chordError: 0.001,
  maxLevel: 8,
  minimumLevel: 0,
  cellBudget: 250000
};

test('capability registry separates stable and beta product methods', () => {
  const capabilities = getCapabilities();
  assert.equal(capabilities.schemaVersion, 1);
  assert.deepEqual(capabilities.methods.map(method => [method.id, method.status]), [
    ['cutcell', 'stable'], ['hybrid', 'beta']
  ]);
  assert.equal(capabilities.methods[0].presets.dense.minimumLevel, 0);
});

test('cut-cell tool validates and serializes local refinement controls', () => {
  const { job, estimate } = validateMeshJob({
    ...common,
    method: 'cutcell',
    paddingFraction: 0.25,
    smallAlpha: 0.05,
    distanceBands: [{ distance: 0.1, targetLevel: 7 }],
    refinementBoxes: [{ xmin: 0, ymin: -0.5, xmax: 2, ymax: 0.5, targetLevel: 6 }]
  });
  assert.equal(estimate.globalLeafFloor, 1);
  const invocation = buildMeshInvocation(job, {
    xyPath: '/tmp/body.xy', prefix: '/tmp/body-cutcell', casePath: '/tmp/case'
  });
  assert.equal(invocation.executableName, 'cartmesh2d_cli');
  assert.deepEqual(invocation.args.slice(-9), [
    '--distance-band', '0.1', '7', '--refine-box', '0', '-0.5', '2', '0.5', '6'
  ]);
});

test('hybrid tool has a distinct contract and reports fallback truthfully', () => {
  const { job } = validateMeshJob({
    ...common,
    method: 'hybrid',
    minimumLevel: 3,
    boundaryLevel: 8,
    nLayers: 4,
    firstThickness: 0.02,
    growthRatio: 1.2,
    domainPadding: 1,
    extrusionThickness: 0.01
  });
  const invocation = buildMeshInvocation(job, {
    xyPath: '/tmp/body.xy', prefix: '/tmp/body-hybrid', casePath: '/tmp/case'
  });
  assert.equal(invocation.executableName, 'cartmesh2d_hybrid_cli');
  assert.equal(invocation.args.length, 9);
  assert.match(invocation.cm2dCandidates[0], /hybrid\.solver\.cm2d$/);
  const fallback = normalizeSummary(job, {
    mesh_mode: 'pure_cutcell_fallback', solver_cells: '420', solver_quality: 'pass'
  });
  assert.equal(fallback.actual_method, 'cutcell-fallback');
  assert.equal(fallback.quality_pass, true);
  const rejected = normalizeSummary(job, {
    mesh_mode: 'hybrid', solver_cells: '908', solver_quality: 'pass',
    quality_contract: 'FAIL', openfoam: 'not_requested'
  });
  assert.equal(rejected.quality_status, 'contract fail');
  assert.equal(rejected.openfoam_cells, null);
});

test('global floor budget rejects accidental whole-domain refinement', () => {
  assert.throws(() => validateMeshJob({
    ...common,
    method: 'cutcell',
    maxLevel: 11,
    minimumLevel: 10,
    cellBudget: 250000,
    paddingFraction: 0.25,
    smallAlpha: 0.05
  }), /全域底格下限/);
});

test('job schema rejects unknown units and malformed refinement boxes', () => {
  assert.throws(() => validateMeshJob({
    ...common, method: 'cutcell', sourceUnits: 'yard',
    paddingFraction: 0.25, smallAlpha: 0.05
  }), /未知 DXF 源单位/);
  assert.throws(() => validateMeshJob({
    ...common, method: 'cutcell', paddingFraction: 0.25, smallAlpha: 0.05,
    refinementBoxes: [{ xmin: 2, ymin: 0, xmax: 1, ymax: 1, targetLevel: 5 }]
  }), /xmax > xmin/);
});
