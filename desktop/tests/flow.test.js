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
  assert.equal(invocation.request.pressurePreconditioner, 'ic0');
  assert.equal(invocation.request.outletBackflow, 'reject');
  assert.deepEqual(invocation.args, ['--mesh', '/tmp/final.solver.cm2d', '--output', '/tmp/run',
    '--case', 'external', '--nu', '0.01', '--speed', '1', '--max-iterations', '1500', '--tolerance', '0.000001',
    '--convection', 'upwind', '--pressure-preconditioner', 'ic0', '--outlet-backflow', 'reject', '--viscous-stress', 'symmetric']);
  const limited = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 10, convection: 'limited-linear' });
  assert.equal(limited.request.convection, 'limited-linear');
  const aggregation = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 10, pressurePreconditioner: 'aggregation' });
  assert.equal(aggregation.request.pressurePreconditioner, 'aggregation');
  assert.deepEqual(aggregation.args.slice(-8), ['--convection', 'upwind', '--pressure-preconditioner', 'aggregation', '--outlet-backflow', 'reject', '--viscous-stress', 'symmetric']);
  assert.equal(limited.request.viscousStress, 'symmetric');
  assert.deepEqual(limited.args.slice(-8), ['--convection', 'limited-linear', '--pressure-preconditioner', 'ic0', '--outlet-backflow', 'reject', '--viscous-stress', 'symmetric']);
  assert.equal(limited.args.at(-1), 'symmetric');
  assert.throws(() => buildFlowInvocation('/tmp/intermediate.cm2d', '/tmp/run', invocation.request), /solver\.cm2d/);
  assert.throws(() => validateFlowRequest({ case: 'rans', nu: 0.01, speed: 1, maxIterations: 10 }), /未知/);
  assert.throws(() => validateFlowRequest({ case: 'external', nu: 0.01, speed: 1, maxIterations: 10, convection: 'central' }), /对流格式/);
  assert.throws(() => validateFlowRequest({ case: 'external', nu: 0.01, speed: 1, maxIterations: 10, pressurePreconditioner: 'amg' }), /压力预条件器/);
  assert.throws(() => validateFlowRequest({ case: 'cavity', nu: 0, speed: 1, maxIterations: 10 }), /大于 0/);
});

test('curved duct can run steady and transient and is retained in output metadata', () => {
  for (const mode of ['steady', 'transient']) {
    const invocation = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
      { case: 'duct', nu: .1, speed: 1, maxIterations: 200, mode, dt: .01, steps: 2 });
    assert.equal(invocation.args[invocation.args.indexOf('--case') + 1], 'duct');
    assert.equal(invocation.args.includes('--time-step'), mode === 'transient');
  }
  assert.equal(validateFlowOutput({ ...summary, case: 'duct' }, fields, 2).summary.case, 'duct');
  assert.match(exportGuide({ result: { counts: { cells: 2 } },
    flow: { summary: { ...summary, case: 'duct' } } }), /几何物面标签并不全是流动壁面/);
});

test('face-frame momentum scheme survives invocation, result validation and export', () => {
  const invocation = buildFlowInvocation('/tmp/final.solver.cm2d', '/tmp/run',
    { case: 'duct', nu: .1, speed: 1, maxIterations: 200, convection: 'face-limited-linear' });
  assert.equal(invocation.args[invocation.args.indexOf('--convection')+1], 'face-limited-linear');
  const result = validateFlowOutput({ ...summary, convection: 'face-limited-linear' }, fields, 2);
  assert.equal(result.summary.convection, invocation.request.convection);
  assert.match(exportGuide({ result: { counts: { cells: 2 } }, flow: result }), /面方向限制/);
});

test('optional outlet backflow diagnostics validate without inventing legacy values', () => {
  const accepted = validateFlowOutput({ ...summary, outletBackflowFaces: 2, outletInflow: 0.25 }, fields, 2);
  assert.equal(accepted.summary.outletBackflowFaces, 2);
  assert.equal(accepted.summary.outletInflow, 0.25);
  const legacy = validateFlowOutput(summary, fields, 2);
  assert.equal(Object.hasOwn(legacy.summary, 'outletBackflowFaces'), false);
  assert.equal(Object.hasOwn(legacy.summary, 'outletInflow'), false);
  assert.throws(() => validateFlowOutput({ ...summary, outletBackflowFaces: -1 }, fields, 2), /回流出口面数/);
  assert.throws(() => validateFlowOutput({ ...summary, outletInflow: -0.1 }, fields, 2), /出口流入量/);
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
  assert.equal(validated.summary.pressurePreconditioner, 'legacy-unspecified');
  assert.equal(validated.summary.pressurePreconditionerInferred, true);
  assert.equal(validated.summary.viscousStress, 'laplacian');
  assert.equal(validated.summary.viscousStressInferred, true);
  assert.throws(() => validateFlowOutput(summary, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' }), /缺少压力离散格式/);
});

test('flow output accepts the limited-linear scheme and rejects unknown or mixed requests', () => {
  const limited = validateFlowOutput({ ...summary, convection: 'limited-linear', pressureDiscretization: 'shared-face-gauss', pressurePreconditioner: 'ic0',
    viscousStress: 'symmetric', forceDefinition: 'shared-face-newtonian-traction',
    forceX: 2, forceY: -3, pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 }, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'limited-linear' });
  assert.equal(limited.summary.convection, 'limited-linear');
  assert.equal(limited.summary.pressureDiscretization, 'shared-face-gauss');
  assert.equal(limited.summary.pressureDiscretizationInferred, false);
  assert.throws(() => validateFlowOutput({ ...summary, convection: 'central' }, fields, 2), /对流格式/);
  assert.throws(() => validateFlowOutput({ ...summary, convection: 'limited-linear', pressureDiscretization: 'shared-face-gauss', pressurePreconditioner: 'ic0',
    viscousStress: 'symmetric', forceDefinition: 'shared-face-newtonian-traction',
    forceX: 2, forceY: -3, pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 }, fields, 2,
    { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' }), /不一致/);
  assert.throws(() => validateFlowOutput({ ...summary, pressureDiscretization: 'cell-centre' }, fields, 2), /压力离散/);
});

test('new symmetric stress output requires force metadata and consistent totals', () => {
  const current = { ...summary, pressureDiscretization: 'shared-face-gauss', pressurePreconditioner: 'ic0', viscousStress: 'symmetric',
    forceDefinition: 'shared-face-newtonian-traction', forceX: 2, forceY: -3,
    pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 };
  const request = { case: 'external', nu: 0.01, speed: 1, maxIterations: 30, convection: 'upwind' };
  const validated = validateFlowOutput(current, fields, 2, request);
  assert.throws(() => validateFlowOutput(current, fields, 2, {...request,tolerance:1e-8}), /不一致/);
  assert.equal(validateFlowOutput({...current,tolerance:1e-8}, fields, 2, {...request,tolerance:1e-8}).summary.tolerance,1e-8);
  assert.equal(validated.summary.forceX, 2);
  assert.throws(() => validateFlowOutput({ ...current, forceDefinition: 'legacy' }, fields, 2, request), /受力定义/);
  assert.throws(() => validateFlowOutput({ ...current, discreteForceX: 2.1 }, fields, 2, request), /离散力不一致/);
  const withoutWall = { ...current };
  delete withoutWall.wallForceX;
  assert.throws(() => validateFlowOutput(withoutWall, fields, 2, request), /wallForceX/);
  assert.throws(() => validateFlowOutput({ ...current, pressurePreconditioner: 'aggregation' }, fields, 2, request), /不一致/);
});

test('a requested live output must record its pressure preconditioner', () => {
  const current = { ...summary, pressureDiscretization: 'shared-face-gauss', viscousStress: 'symmetric',
    forceDefinition: 'shared-face-newtonian-traction', forceX: 2, forceY: -3,
    pressureForceX: 1, pressureForceY: -1, discreteForceX: 2, discreteForceY: -3,
    wallForceX: 2, wallForceY: -3, wallViscousForceX: 1, wallViscousForceY: -2 };
  const request = { case: 'external', nu: 0.01, speed: 1, maxIterations: 30 };
  assert.throws(() => validateFlowOutput(current, fields, 2, request), /缺少压力预条件器/);
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


test('desktop tolerance can only tighten the historical flow default',()=>{
  const request={case:'external',nu:.1,speed:1,maxIterations:500,tolerance:1e-9};
  const invocation=buildFlowInvocation('/tmp/final.solver.cm2d','/tmp/run',request);
  assert.equal(invocation.args[invocation.args.indexOf('--tolerance')+1],'1e-9');
  for(const tolerance of [0,-1,NaN,Infinity,1e-5,1e-13])
    assert.throws(()=>validateFlowRequest({...request,tolerance}));
});

test('steady acceleration is explicit, unavailable in time marching, and bound to recorded diagnostics', () => {
  const request={case:'external',nu:.01,speed:1,maxIterations:30,steadyAcceleration:'anderson'};
  const invocation=buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/flow',request);
  assert.equal(invocation.args[invocation.args.indexOf('--steady-acceleration')+1],'anderson');
  assert.equal(validateFlowRequest({...request,steadyAcceleration:undefined}).steadyAcceleration,'none');
  assert.throws(()=>validateFlowRequest({...request,steadyAcceleration:'unknown'}),/稳态加速/);
  assert.throws(()=>validateFlowRequest({...request,mode:'transient',dt:.02,steps:1}),/稳态加速/);
  const current={...summary,pressureDiscretization:'shared-face-gauss',pressurePreconditioner:'ic0',viscousStress:'symmetric',
    forceDefinition:'shared-face-newtonian-traction',forceX:2,forceY:-3,pressureForceX:1,pressureForceY:-1,
    discreteForceX:2,discreteForceY:-3,wallForceX:2,wallForceY:-3,wallViscousForceX:1,wallViscousForceY:-2,
    steadyAcceleration:'anderson',accelerationCandidates:6,accelerationAccepted:4,accelerationRejected:2};
  assert.equal(validateFlowOutput(current,fields,2,request).summary.steadyAcceleration,'anderson');
  assert.throws(()=>validateFlowOutput({...current,accelerationAccepted:5},fields,2,request),/统计不一致/);
  assert.throws(()=>validateFlowOutput({...current,accelerationCandidates:NaN},fields,2,request),/统计无效/);
  assert.throws(()=>validateFlowOutput({...current,steadyAcceleration:undefined},fields,2,request),/缺少模式/);
  const legacy={...current};for(const key of ['steadyAcceleration','accelerationCandidates','accelerationAccepted','accelerationRejected'])delete legacy[key];
  assert.throws(()=>validateFlowOutput(legacy,fields,2,request),/不一致/);
});

test('linear efficiency controls bind requests and enforce strict final certification',()=>{
  const request={case:'external',nu:.01,speed:1,maxIterations:30,linearPolicy:'adaptive',velocityRelaxation:.8};
  for(const mode of ['steady','transient','adaptive']) {
    const q={...request,mode,dt:.02,steps:1,endTime:.1};
    const invoke=buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/flow',q);
    assert.equal(invoke.args[invoke.args.indexOf('--linear-policy')+1],'adaptive');
    assert.equal(invoke.args[invoke.args.indexOf('--velocity-relaxation')+1],'0.8');
  }
  const current={...summary,pressureDiscretization:'shared-face-gauss',pressurePreconditioner:'ic0',viscousStress:'symmetric',
    forceDefinition:'shared-face-newtonian-traction',forceX:2,forceY:-3,pressureForceX:1,pressureForceY:-1,
    discreteForceX:2,discreteForceY:-3,wallForceX:2,wallForceY:-3,wallViscousForceX:1,wallViscousForceY:-2,
    adaptiveLinear:true,strictLinearFinal:true,velocityRelaxation:.8};
  const validated=validateFlowOutput(current,fields,2,request);
  assert.equal(validated.summary.linearPolicy,'adaptive');assert.equal(validated.summary.velocityRelaxation,.8);
  for(const patch of [{adaptiveLinear:false},{adaptiveLinear:undefined},{adaptiveLinear:'true'},
    {strictLinearFinal:false},{strictLinearFinal:undefined},{velocityRelaxation:.6},{velocityRelaxation:undefined},
    {velocityRelaxation:0},{velocityRelaxation:1.1}])
    assert.throws(()=>validateFlowOutput({...current,...patch},fields,2,request));
  assert.equal(validateFlowOutput({...current,converged:false,status:'iteration_limit',strictLinearFinal:false},fields,2,request).summary.converged,false);
  for(const patch of [{linearPolicy:'unknown'},{velocityRelaxation:0},{velocityRelaxation:1.1},{velocityRelaxation:NaN}])
    assert.throws(()=>validateFlowRequest({...request,...patch}));
});

test('system Cholesky is explicit, platform-bound for new solves, and retained in exported results',()=>{
  const request={case:'external',nu:.01,speed:1,maxIterations:30,pressurePreconditioner:'cholesky'};
  if(process.platform==='darwin') {
    const invocation=buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/flow',request);
    assert.equal(invocation.args[invocation.args.indexOf('--pressure-preconditioner')+1],'cholesky');
  } else assert.throws(()=>validateFlowRequest(request),/macOS/);
  const result=validateFlowOutput({...summary,pressurePreconditioner:'cholesky'},fields,2);
  assert.equal(result.summary.pressurePreconditioner,'cholesky');
  assert.match(exportGuide({result:{counts:{cells:2}},flow:result}),/系统稀疏 Cholesky/);
});
