'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const { FLOW_OUTPUT_SUFFIXES, buildFlowInvocation, commitFlowFiles, parseFlowProgress,
        validateFlowOutput, validateFlowRequest } = require('../src/core/flow');

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
  assert.deepEqual(invocation.args, ['--mesh', '/tmp/final.solver.cm2d', '--output', '/tmp/run',
    '--case', 'external', '--nu', '0.01', '--speed', '1', '--max-iterations', '1500']);
  assert.throws(() => buildFlowInvocation('/tmp/intermediate.cm2d', '/tmp/run', invocation.request), /solver\.cm2d/);
  assert.throws(() => validateFlowRequest({ case: 'rans', nu: 0.01, speed: 1, maxIterations: 10 }), /未知/);
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
