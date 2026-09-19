'use strict';

const assert = require('node:assert/strict');
const test = require('node:test');
const { parseCm2d } = require('../src/core/cm2d');
const {
  buildThermalInvocation,
  parseThermalProgress,
  thermalBoundaryCsv,
  validateThermalOutput,
  validateThermalRequest
} = require('../src/core/thermal');

const boundaries = {
  wall: { kind: 'value', value: 1, inflowValue: 0 },
  inlet: { kind: 'flux', value: 0, inflowValue: 0 },
  outlet: { kind: 'flux', value: 0, inflowValue: 0 },
  top: { kind: 'flux', value: 0, inflowValue: 0 },
  bottom: { kind: 'flux', value: 0, inflowValue: 0 }
};

function request(overrides = {}) {
  return {
    case: 'channel', nu: 0.1, speed: 1, maxIterations: 100,
    convection: 'limited-linear', pressurePreconditioner: 'ic0',
    mode: 'transient', dt: 0.05, steps: 4,
    diffusivity: 0.1, initial: 0, source: 0,
    scalarConvection: 'upwind', boundaries,
    ...overrides
  };
}

const RECTANGLE = [
  'CM2D 1', 'VERTICES 4',
  '0 0 0', '1 1 0', '2 1 1', '3 0 1',
  'EDGES 4',
  '0 0 1 0 -1 2', '1 1 2 0 -1 2', '2 2 3 0 -1 2', '3 3 0 0 -1 2',
  'CELLS 1', '0 0 0 1 4 0 1 2 3 4 0 1 2 3', 'AUDIT 0'
].join('\n');

const EMBEDDED_RECTANGLE = [
  'CM2D 1', 'VERTICES 6',
  '0 0 0', '1 1 0', '2 1 1', '3 0 1', '4 0.4 0.4', '5 0.6 0.4',
  'EDGES 5',
  '0 0 1 0 -1 2', '1 1 2 0 -1 2', '2 2 3 0 -1 2', '3 3 0 0 -1 2',
  '4 4 5 0 -1 1',
  'CELLS 1', '0 0 0 1 4 0 1 2 3 5 0 1 2 3 4', 'AUDIT 0'
].join('\n');

test('thermal request validation rejects coercion and preserves the physical contract', () => {
  const validated = validateThermalRequest(request());
  assert.equal(validated.diffusivity, 0.1);
  assert.equal(validated.initial, 0);
  assert.equal(validated.boundaries.wall.value, 1);
  for (const [field, value] of [['diffusivity', '0.1'], ['initial', false], ['source', '0']]) {
    assert.throws(() => validateThermalRequest(request({ [field]: value })), /必须是有限数/);
  }
  assert.throws(() => validateThermalRequest(request({ scalarConvection: 'central' })), /未知温度对流格式/);
  assert.throws(() => validateThermalRequest(request({ boundaries: {
    ...boundaries, wall: { kind: 'value', value: '1', inflowValue: 0 }
  } })), /wall 必须是有限数/);
  assert.throws(() => validateThermalRequest(request({ boundaries: {
    ...boundaries, wall: { kind: 'value', value: 1, inflowValue: true }
  } })), /wall 回流温度.*必须是有限数/);
});

test('thermal boundary CSV groups rectangle domain and embedded external wall', () => {
  const rectangle = parseCm2d(RECTANGLE);
  const domainRows = thermalBoundaryCsv(rectangle, request()).trim().split('\n');
  assert.deepEqual(domainRows, [
    'face,type,value,inflowValue', '0,flux,0,0', '1,flux,0,0',
    '2,flux,0,0', '3,flux,0,0'
  ]);

  const external = parseCm2d(EMBEDDED_RECTANGLE);
  const externalRows = thermalBoundaryCsv(external, request({ case: 'external' })).trim().split('\n');
  assert.equal(externalRows.at(-1), '4,value,1,0');
  assert.equal(externalRows.length, 6);
  assert.throws(() => thermalBoundaryCsv(parseCm2d(EMBEDDED_RECTANGLE), request()), /轴对齐矩形外边界/);
});

test('thermal invocation rejects unsafe mesh and restart inputs and emits native arguments', () => {
  assert.throws(() => buildThermalInvocation('/tmp/mesh.cm2d', '/tmp/run', '/tmp/boundary.csv', request()), /solver\.cm2d/);
  assert.throws(() => buildThermalInvocation('/tmp/final.solver.cm2d', '/tmp/run', '/tmp/boundary.csv', request({ resume: true }), null), /联合续算状态/);
  const invocation = buildThermalInvocation('/tmp/final.solver.cm2d', '/tmp/run', '/tmp/boundary.csv', request());
  assert.equal(invocation.executable, 'cartmesh2d_transport_cli');
  assert.deepEqual(invocation.args, [
    '--mesh', '/tmp/final.solver.cm2d', '--output', '/tmp/run', '--boundary', '/tmp/boundary.csv',
    '--evolve-flow', 'channel', '--flow-nu', '0.1', '--flow-speed', '1',
    '--flow-max-iterations', '100', '--flow-convection', 'limited-linear',
    '--pressure-preconditioner', 'ic0', '--diffusivity', '0.1', '--source', '0',
    '--initial', '0', '--convection', 'upwind', '--dt', '0.05', '--steps', '4'
  ]);
});

test('accepted thermal progress requires numeric, converged positive-step records', () => {
  const line = JSON.stringify({ type: 'thermal-time-step', step: 1, time: 0.05, accepted: 1,
    flowIterations: 20, momentumResidual: 1e-10, continuity: 1e-12,
    scalarIterations: 12, scalarResidual: 2e-10, heatContent: 0.1,
    globalBalance: 0, maxCourant: 0.2 });
  assert.deepEqual(parseThermalProgress(line), JSON.parse(line));
  for (const [field, value] of [['accepted', true], ['time', '0.05'], ['step', 1.5], ['maxCourant', -1]]) {
    const parsed = JSON.parse(line);
    parsed[field] = value;
    assert.throws(() => parseThermalProgress(JSON.stringify(parsed)), /监测|有限|步数|CFL/);
  }
  assert.equal(parseThermalProgress('progress=1'), null);
});

function thermalContractFixture() {
  const mesh = parseCm2d(RECTANGLE);
  const input = request({ case: 'external' });
  const summary = {
    format: 'cartmesh2d-scalar-transport-v1', status: 'converged', converged: true,
    evolvingFlow: true, cells: 1, faces: 4, steps: 4,
    time: 0.2, acceptedTime: 0.2, carrierTime: 0.2, timeStep: 0.05,
    diffusivity: 0.1, flowNu: 0.1, flowSpeed: 1, flowCase: 'external',
    convection: 'upwind', flowConvection: 'limited-linear', constantSource: 0,
    initialValue: 0, minValue: 0, maxValue: 0, globalBalance: 0,
    residualNorm: 1e-10, maxDiagonalScaledImbalance: 1e-10
  };
  const cells = [
    'cell,x,y,area,value,previous,sourceIntegral,temporalIntegral,exact',
    '0,0.5,0.5,1,0,0,0,0,0'
  ].join('\n') + '\n';
  const history = [
    'step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant',
    '1,0.05,1,10,1e-10,1e-12,10,1e-10,0,0,0.1',
    '2,0.1,1,10,1e-10,1e-12,10,1e-10,0,0,0.1',
    '3,0.15,1,10,1e-10,1e-12,10,1e-10,0,0,0.1',
    '4,0.2,1,10,1e-10,1e-12,10,1e-10,0,0,0.1'
  ].join('\n') + '\n';
  const joint = [
    'CARTMESH2D_THERMAL_CHECKPOINT 1',
    'COUPLING new-time-flux-Euler-v1',
    'SCALAR 1 0',
    'FLOW',
    'CARTMESH2D_FLOW_CHECKPOINT 1',
    'TIME 0.2',
    'FLUX 4 0 0 0 0',
    ''
  ].join('\n');
  return { summary, mesh, input, cells, history, joint };
}

test('synthetic thermal output contract validates mesh, joint state, values and history', () => {
  const fixture = thermalContractFixture();
  const validated = validateThermalOutput(fixture.summary, fixture.cells, fixture.history,
    fixture.joint, fixture.mesh, fixture.input);
  assert.equal(validated.summary.steps, 4);
  assert.equal(validated.fields.cells.length, fixture.mesh.cells.length);
  assert.equal(validated.history.length, 4);
});

test('thermal output contract rejects tampered joint time/value, summary types, history and request', () => {
  const fixture = thermalContractFixture();
  const cases = [
    ['joint time', { joint: fixture.joint.replace('TIME 0.2', 'TIME 0.25') }],
    ['joint scalar value', { joint: fixture.joint.replace(/^SCALAR 1 0$/m, 'SCALAR 1 9') }],
    ['missing summary field', { summary: (({ flowNu, ...rest }) => rest)(fixture.summary) }],
    ['NaN summary field', { summary: { ...fixture.summary, flowNu: NaN } }],
    ['boolean summary field', { summary: { ...fixture.summary, steps: true } }],
    ['inconsistent history', { history: fixture.history.replace('4,0.2,1,10', '3,0.2,1,10') }],
    ['physical request mismatch', { input: request({ case: 'external', diffusivity: 0.2 }) }]
  ];
  for (const [name, mutation] of cases) {
    assert.throws(() => validateThermalOutput(
      mutation.summary || fixture.summary,
      mutation.cells || fixture.cells,
      mutation.history || fixture.history,
      mutation.joint || fixture.joint,
      fixture.mesh,
      mutation.input || fixture.input
    ), /热输运|温度|联合|历史|请求|有限|不一致/, name);
  }
});
