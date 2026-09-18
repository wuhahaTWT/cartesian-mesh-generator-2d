'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const { FLOW_OUTPUT_SUFFIXES, buildFlowInvocation, commitFlowFiles, parseFlowProgress,
        validateFlowOutput, validateFlowRequest } = require('../src/core/flow');
const { exportGuide } = require('../src/core/export-guide');

const summary = {
  format: 'cartmesh2d-flow-summary-v1', case: 'external', nu: 0.01, speed: 1,
  cells: 2, iterations: 27, converged: true, status: 'converged',
  continuity: 1e-10, globalImbalance: -2e-12, globalRelativeImbalance: 2e-12, tolerance: 1e-6, velocityChange: 2e-9,
  pressureChange: 3e-9, momentumResidual: 4e-9
};
const fields = {
  format: 'cartmesh2d-flow-v1',
  cells: [
    { id: 1, u: 0.3, v: 0.4, p: -0.1, speed: 0.5 },
    { id: 0, u: 1, v: 0, p: 0.2, speed: 1 }
  ]
};

test('flow invocation uses the final solver mesh and the small supported parameter set', () => {
  assert.deepEqual(FLOW_OUTPUT_SUFFIXES,
    ['.json', '.fields.json', '.vtk', '.residuals.csv', '.cells.csv', '.faces.csv']);
  const invocation = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { case: 'external', nu: '0.01', speed: '1', maxIterations: '1500' });
  assert.equal(invocation.executable, 'cartmesh2d_flow_cli');
  assert.equal(invocation.request.convection, 'upwind');
  assert.deepEqual(invocation.args, ['--mesh', '/tmp/final.solver.cm2d', '--output', '/tmp/run',
    '--case', 'external', '--nu', '0.01', '--speed', '1', '--max-iterations', '1500',
    '--convection', 'upwind', '--viscous-stress', 'symmetric']);
  const limited = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 10, convection: 'limited-linear' });
  assert.equal(limited.request.convection, 'limited-linear');
  assert.equal(limited.request.viscousStress, 'symmetric');
  assert.deepEqual(limited.args.slice(-4), ['--convection', 'limited-linear', '--viscous-stress', 'symmetric']);
  assert.equal(limited.args.at(-1), 'symmetric');
  assert.throws(() => buildFlowInvocation('/tmp/intermediate.cm2d', '/tmp/run', invocation.request), /solver\.cm2d/);
  assert.throws(() => validateFlowRequest({ case: 'rans', nu: 0.01, speed: 1, maxIterations: 10 }), /未知/);
  assert.throws(() => validateFlowRequest({ case: 'external', nu: 0.01, speed: 1, maxIterations: 10, convection: 'central' }), /对流格式/);
  assert.throws(() => validateFlowRequest({ case: 'cavity', nu: 0, speed: 1, maxIterations: 10 }), /大于 0/);
});

test('progress accepts only complete native flow progress records', () => {
  assert.deepEqual(parseFlowProgress(JSON.stringify({ type: 'flow-progress', iteration: 12,
    continuity: 1e-5, velocityChange: 2e-5, momentumResidual: 3e-5 })), {
    type: 'flow-progress', iteration: 12, continuity: 1e-5,
    velocityChange: 2e-5, momentumResidual: 3e-5
  });
  assert.equal(parseFlowProgress('iteration=12'), null);
  assert.equal(parseFlowProgress('{not json}'), null);
  assert.throws(() => parseFlowProgress(JSON.stringify({ type: 'flow-progress', iteration: 1,
    continuity: null, velocityChange: 0, momentumResidual: 0 })), /finite|有限/);
});

test('flow outputs cover every final cell and are reordered by native cell id', () => {
  const validated = validateFlowOutput(summary, fields, 2);
  assert.deepEqual(validated.fields.cells.map(cell => cell.id), [0, 1]);
  assert.equal(validated.summary.status, 'converged');
  assert.equal(validated.summary.convection, 'upwind');
  assert.equal(validated.summary.convectionInferred, true);
  assert.equal(validated.summary.pressureDiscretization, 'legacy-unspecified');
  assert.equal(validated.summary.pressureDiscretizationInferred, true);
  assert.equal(validated.summary.viscousStress, 'laplacian');
  assert.equal(validated.summary.viscousStressInferred, true);
  assert.throws(() => validateFlowOutput(summary, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' }), /缺少压力离散格式/);
});

test('flow output accepts the limited-linear scheme and rejects unknown or mixed requests', () => {
  const limited = validateFlowOutput({ ...summary, convection: 'limited-linear', pressureDiscretization: 'shared-face-gauss',
    viscousStress: 'symmetric', forceDefinition: 'shared-face-newtonian-traction',
    forceX: 2, forceY: -3, pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 }, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'limited-linear' });
  assert.equal(limited.summary.convection, 'limited-linear');
  assert.equal(limited.summary.pressureDiscretization, 'shared-face-gauss');
  assert.equal(limited.summary.pressureDiscretizationInferred, false);
  assert.throws(() => validateFlowOutput({ ...summary, convection: 'central' }, fields, 2), /对流格式/);
  assert.throws(() => validateFlowOutput({ ...summary, convection: 'limited-linear', pressureDiscretization: 'shared-face-gauss',
    viscousStress: 'symmetric', forceDefinition: 'shared-face-newtonian-traction',
    forceX: 2, forceY: -3, pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 }, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' }), /不一致/);
  assert.throws(() => validateFlowOutput({ ...summary, pressureDiscretization: 'cell-centre' }, fields, 2), /压力离散/);
});

test('new symmetric stress output requires force metadata and consistent totals', () => {
  const current = { ...summary, pressureDiscretization: 'shared-face-gauss', viscousStress: 'symmetric',
    forceDefinition: 'shared-face-newtonian-traction', forceX: 2, forceY: -3,
    pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 };
  const request = { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' };
  const validated = validateFlowOutput(current, fields, 2, request);
  assert.equal(validated.summary.forceX, 2);
  assert.throws(() => validateFlowOutput({ ...current, forceDefinition: 'legacy' }, fields, 2, request), /受力定义/);
  assert.throws(() => validateFlowOutput({ ...current, discreteForceX: 2.1 }, fields, 2, request), /离散力不一致/);
  const withoutWall = { ...current };
  delete withoutWall.wallForceX;
  assert.throws(() => validateFlowOutput(withoutWall, fields, 2, request), /wallForceX/);
});

test('export guide distinguishes recorded shared-face pressure from legacy summaries', () => {
  const result = { counts: { cells: 2 }, gates: {} };
  const legacy = exportGuide({ result, flow: { summary } });
  assert.match(legacy, /压力离散：旧结果未记录/);
  assert.doesNotMatch(legacy, /压力离散：共享面压力/);
  const current = exportGuide({ result, flow: {
    summary: { ...summary, pressureDiscretization: 'shared-face-gauss' }
  } });
  assert.match(current, /压力离散：共享面压力/);
  assert.match(exportGuide({ result, flow: { summary: {
    ...summary, pressureDiscretization: 'shared-face-gauss', viscousStress: 'symmetric'
  } } }), /压力力和黏性力在同一组共享面上积分/);
});

test('iteration-limit output is valid but never reported as converged', () => {
  const valid = validateFlowOutput({ ...summary, status: 'iteration_limit', converged: false }, fields, 2);
  assert.equal(valid.summary.converged, false);
  assert.throws(() => validateFlowOutput({ ...summary, status: 'iteration_limit' }, fields, 2), /矛盾/);
  assert.throws(() => validateFlowOutput({ ...summary, globalRelativeImbalance: 1e-4 }, fields, 2), /未达到/);
});

test('malformed field ids, counts and vector magnitudes fail closed', () => {
  assert.throws(() => validateFlowOutput(summary, { ...fields, cells: [fields.cells[0]] }, 2), /全部单元/);
  assert.throws(() => validateFlowOutput(summary, { ...fields,
    cells: [{ ...fields.cells[0], id: 0 }, { ...fields.cells[1], id: 0 }] }, 2), /不重复/);
  assert.throws(() => validateFlowOutput(summary, { ...fields,
    cells: [{ ...fields.cells[0], id: 0, speed: 9 }, fields.cells[1]] }, 2), /速度分量/);
  assert.throws(() => validateFlowOutput({ ...summary, cells: 3 }, fields, 2), /单元数/);
});

test('a partial flow-file copy removes every top-level committed name', async () => {
  const copied = [];
  const removed = [];
  const fileSystem = {
    async copyFile(source, destination) {
      copied.push([source, destination]);
      if (source.endsWith('.faces.csv')) throw new Error('faces copy failed');
    },
    async rm(file) { removed.push(file); }
  };
  const entries = FLOW_OUTPUT_SUFFIXES.map(suffix => ({
    source: `/pending/flow${suffix}`, destination: `/result/mesh.flow${suffix}`
  }));
  await assert.rejects(commitFlowFiles(fileSystem, entries), /faces copy failed/);
  assert.equal(copied.length, entries.length);
  assert.deepEqual(removed.sort(), entries.map(entry => entry.destination).sort());
});
