'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const {
  buildFlowInvocation, parseFlowProgress, validateFlowOutput,
  validateFlowRequest, validateTimeHistory
} = require('../src/core/flow');
const { exportGuide } = require('../src/core/export-guide');

const fields = {
  format: 'cartmesh2d-flow-v1',
  cells: [{ id: 0, u: 0.3, v: 0.4, p: 0, speed: 0.5 }]
};

const transientRequest = {
  case: 'channel', nu: 0.01, speed: 1, maxIterations: 100,
  convection: 'upwind', mode: 'transient', dt: 0.1, steps: 2
};

const transientSummary = {
  format: 'cartmesh2d-flow-summary-v1', case: 'channel', nu: 0.01, speed: 1,
  cells: 1, iterations: 12, converged: true, status: 'converged',
  continuity: 1e-10, globalImbalance: 0, globalRelativeImbalance: 1e-10,
  tolerance: 1e-6, velocityChange: 2e-9, pressureChange: 3e-9, momentumResidual: 4e-9,
  convection: 'upwind', pressurePreconditioner: 'ic0', pressureDiscretization: 'shared-face-gauss', viscousStress: 'symmetric',
  forceDefinition: 'shared-face-newtonian-traction',
  forceX: 2, forceY: -3, pressureForceX: 1, pressureForceY: -1,
  discreteForceX: 2, discreteForceY: -3, wallForceX: 2, wallForceY: -3,
  wallViscousForceX: 1, wallViscousForceY: -2,
  temporalDiscretization: 'backward-euler',
  temporalFaceInterpolation: 'old-and-iteration-flux-defect-skew-corrected-v2',
  time: 0.4, dt: 0.1, acceptedTime: 0.4, requestedSteps: 2,
  completedSteps: 2, maxCourant: 0.5
};

const timeHistoryHeader = 'step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY';
const validHistory = [
  timeHistoryHeader,
  '1,0.3,0.1,1,10,1e-9,1e-10,0.5,0.2,2,-3',
  '2,0.4,0.1,1,12,4e-9,1e-10,0.5,0.2,2,-3'
].join('\n');

test('transient invocation requires dt and steps and puts restart only in backend args', () => {
  assert.throws(() => validateFlowRequest({ ...transientRequest, dt: undefined }), /时间步长/);
  assert.throws(() => validateFlowRequest({ ...transientRequest, steps: undefined }), /时间步数/);
  assert.throws(() => buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { ...transientRequest, resume: true }), /重启状态/);
  const fresh = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run', transientRequest, '/tmp/ignored.checkpoint');
  assert.deepEqual(fresh.args.slice(-4), ['--time-step', '0.1', '--steps', '2',]);
  assert.equal(fresh.args.includes('--restart'), false);
  const resumed = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { ...transientRequest, resume: true }, '/tmp/state.checkpoint');
  assert.deepEqual(resumed.args.slice(-6), ['--time-step', '0.1', '--steps', '2', '--restart', '/tmp/state.checkpoint']);
});

test('physical progress events retain accepted time-step information', () => {
  const accepted = parseFlowProgress(JSON.stringify({ type: 'flow-time-step', time: 0.3, step: 1,
    maxCourant: 0.5, kineticEnergy: 0.2, forceX: 2, forceY: -3 }));
  assert.deepEqual(accepted, { type: 'flow-time-step', time: 0.3, step: 1, maxCourant: 0.5,
    kineticEnergy: 0.2, forceX: 2, forceY: -3 });
  assert.deepEqual(parseFlowProgress(JSON.stringify({ type: 'flow-progress', iteration: 12,
    continuity: 1e-10, velocityChange: 2e-9, momentumResidual: 4e-9, time: 0.4, timeStep: 2 })),
    { type: 'flow-progress', iteration: 12, continuity: 1e-10,
      velocityChange: 2e-9, momentumResidual: 4e-9, time: 0.4, timeStep: 2 });
});

test('physical accepted progress rejects negative or incomplete monitoring values', () => {
  for (const field of ['kineticEnergy', 'forceX', 'forceY']) {
    const value = { type: 'flow-time-step', time: 0.3, step: 1, maxCourant: 0.5,
      kineticEnergy: 0.2, forceX: 2, forceY: -3 };
    delete value[field];
    assert.throws(() => parseFlowProgress(JSON.stringify(value)), /进度|监测|finite|有限/, `missing ${field}`);
  }
  for (const field of ['kineticEnergy', 'forceX', 'forceY']) {
    const value = { type: 'flow-time-step', time: 0.3, step: 1, maxCourant: 0.5,
      kineticEnergy: 0.2, forceX: 2, forceY: -3 };
    value[field] = field === 'kineticEnergy' ? -0.1 : NaN;
    assert.throws(() => parseFlowProgress(JSON.stringify(value)), /进度|监测|finite|有限|负/, `invalid ${field}`);
  }
  assert.throws(() => parseFlowProgress(JSON.stringify({ type: 'flow-time-step', time: 0.3, step: 1,
    maxCourant: 0.5, kineticEnergy: ' 0.2 ', forceX: 2, forceY: -3 })), /进度|finite|有限/);
});

test('a resumed transient summary validates against start time and complete force schema', () => {
  const validated = validateFlowOutput(transientSummary, fields, 1, transientRequest, 0.2);
  assert.equal(validated.summary.acceptedTime, 0.4);
  assert.equal(validated.summary.wallViscousForceY, -2);
  assert.deepEqual(validateTimeHistory(validHistory, transientSummary, 0.2).map(row => row.step), [1, 2]);
});

test('time history rejects failed acceptance, negative values, inconsistent dt, last CFL, and counts', () => {
  const mutations = [
    ['false acceptance', validHistory.replace('1,0.3,0.1,1,10', '1,0.3,0.1,0,10')],
    ['negative CFL', validHistory.replace(',0.5,0.2,2,-3', ',-0.5,0.2,2,-3')],
    ['inconsistent dt', validHistory.replace('2,0.4,0.1', '2,0.4,0.2')],
    ['last CFL', validHistory.replace('2,0.4,0.1,1,12,4e-9,1e-10,0.5', '2,0.4,0.1,1,12,4e-9,1e-10,0.6')]
  ];
  for (const [name, text] of mutations) assert.throws(() => validateTimeHistory(text, transientSummary, 0.2), /时间历史|失败|负/, name);
  assert.throws(() => validateTimeHistory(validHistory, { ...transientSummary, completedSteps: 1 }, 0.2), /最终状态/);
});

test('steady output rejects unexpected time fields', () => {
  const steady = {
    ...transientSummary,
    temporalDiscretization: undefined, temporalFaceInterpolation: undefined,
    time: undefined, dt: undefined, acceptedTime: undefined,
    requestedSteps: undefined, completedSteps: undefined, maxCourant: undefined
  };
  assert.throws(() => validateFlowOutput({ ...steady, time: 0.1 }, fields, 1), /非定常离散格式|时间字段/);
});

test('transient export guide names restart, time history, accepted time, and monitoring', () => {
  const guide = exportGuide({ result: { counts: { cells: 1 }, gates: {} }, flow: { summary: transientSummary } });
  assert.match(guide, /checkpoint/);
  assert.match(guide, /time-history/);
  assert.match(guide, /已接受到 t=0\.4 s/);
  assert.match(guide, /全部物理时间步/);
});

test('progress numeric fields and booleans do not accept coercive encodings', () => {
  assert.throws(() => parseFlowProgress(JSON.stringify({ type: 'flow-time-step', time: '0.3', step: 1,
    maxCourant: 0.5, kineticEnergy: 0.2, forceX: 2, forceY: -3 })), /物理时间|进度|finite|有限/);
  assert.throws(() => parseFlowProgress(JSON.stringify({ type: 'flow-time-step', time: 0.3, step: true,
    maxCourant: 0.5, kineticEnergy: 0.2, forceX: 2, forceY: -3 })), /时间步|进度/);
  assert.throws(() => parseFlowProgress(JSON.stringify({ type: 'flow-time-step', time: 0.3, step: 1,
    maxCourant: ' 0.5 ', kineticEnergy: 0.2, forceX: 2, forceY: -3 })), /CFL|进度|finite|有限/);
});

test('failed candidate summary and history preserve accepted time separately from candidate time', () => {
  const failed = { ...transientSummary, status: 'time_step_not_converged', converged: false,
    time: 0.4, acceptedTime: 0.3, requestedSteps: 2, completedSteps: 1 };
  const candidateHistory = [
    timeHistoryHeader,
    '1,0.3,0.1,1,10,1e-9,1e-10,0.5,0.2,2,-3',
    '2,0.4,0.1,0,12,4e-9,1e-10,0.5,0.2,2,-3'
  ].join('\n');
  const validated = validateFlowOutput(failed, fields, 1, transientRequest, 0.2);
  assert.equal(validated.summary.acceptedTime, 0.3);
  assert.equal(validateTimeHistory(candidateHistory, failed, 0.2).at(-1).accepted, 0);
});

test('three cross-mode contract violations fail closed', () => {
  assert.throws(() => validateFlowRequest({ ...transientRequest, mode: 'steady', resume: true }), /稳态模式/);
  assert.throws(() => validateFlowOutput(transientSummary, fields, 1,
    { ...transientRequest, mode: 'steady' }), /时间模式|非定常时间|时间步长/);
  assert.throws(() => validateFlowOutput({ ...transientSummary, temporalDiscretization: undefined }, fields, 1), /非定常离散格式/);
});
