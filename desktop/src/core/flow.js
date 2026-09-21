'use strict';
const {normalizeInitialVortex,validateInitialVortexOutput}=require('./initial-vortex');
const { validateAdaptiveSummary, validateAttemptHistory } = require('./adaptive-flow');
const { validateWallLoads } = require('./wall-loads');
const { validateBoundaryFluxes } = require('./boundary-flux');
const { normalizeBoundaryDefinition, sameConditions, conditions } = require('./flow-boundaries');

const FLOW_CASES = Object.freeze({
  custom: {
    id: 'custom', label: '命名边界',
    scope: '按最终网格逐面指定速度入口、压力出口、定压开口、对称面或静止／移动壁面。参考速度仅用于归一化。定压开口限水平／竖直方向，允许法向进出；普通压力出口仍拒绝回流。对称面限水平／竖直方向，不穿透且切向无摩擦。'
  },
  external: {
    id: 'external', label: '外流',
    scope: '矩形域：左侧沿 +X 的恒速入口、右侧压力出口、上下滑移，物面无滑移。尾迹加密方向不会旋转入口速度。'
  },
  channel: {
    id: 'channel', label: '通道',
    scope: '矩形内流：左侧抛物线入口、右侧 p=0，其余边界无滑移。'
  },
  duct: {
    id: 'duct', label: '曲壁通道／喷管',
    scope: '内流：最左侧竖直端面为均匀速度入口，最右侧竖直端面为 p=0 出口，其余曲壁及障碍物无滑移。入口和出口必须是明确的竖直端面。'
  },
  cavity: {
    id: 'cavity', label: '顶盖方腔',
    scope: '方形内流：顶盖移动，其余边界无滑移，并固定参考压力。'
  }
});
const FLOW_CONVECTION_SCHEMES = Object.freeze({
  upwind: {
    id: 'upwind', label: '一阶迎风',
    description: '稳健的默认格式，数值扩散较大。'
  },
  'limited-linear': {
    id: 'limited-linear', label: '线性迎风（限制重构）',
    description: '减少数值扩散，但不保证所有工况都更准确。'
  },
  'face-limited-linear': {
    id: 'face-limited-linear', label: '线性迎风（面方向限制）',
    description: '沿网格面的法向、切向限制速度，减少旋转工况的坐标方向影响。'
  }
});
const FLOW_PRESSURE_PRECONDITIONERS = Object.freeze({
  ic0: {
    id: 'ic0', label: '标准（IC0）',
    description: '标准不完全 Cholesky 压力预条件。'
  },
  aggregation: {
    id: 'aggregation', label: '多重网格（试验）',
    description: '试验性聚合多重网格；不保证更快。'
  }
});
const FLOW_OUTLET_BACKFLOW_MODES = Object.freeze({
  reject: { id: 'reject', label: '检测到回流时停止', description: '压力出口保持给定压力；检测到回流时停止。' },
  'normal-inlet': { id: 'normal-inlet', label: '允许法向回流（试验）', description: '压力出口保持给定压力；回流方向垂直出口，适用于当前预设的竖直右出口。' }
});
const LEGACY_PRESSURE_PRECONDITIONER = 'legacy-unspecified';
const PRESSURE_DISCRETIZATION = 'shared-face-gauss';
const LEGACY_PRESSURE_DISCRETIZATION = 'legacy-unspecified';
const VISCOUS_STRESS = 'symmetric';
const LEGACY_VISCOUS_STRESS = 'laplacian';
const FORCE_DEFINITION = 'shared-face-newtonian-traction';
const FLOW_OUTPUT_SUFFIXES = Object.freeze([
  '.json', '.fields.json', '.vtk', '.residuals.csv', '.cells.csv', '.faces.csv'
]);

const finite = (value, name) => {
  if (!['number', 'string'].includes(typeof value) || (typeof value === 'string' && !value.trim())) throw new Error(`${name} 必须是有限数。`);
  const number = Number(value);
  if (!Number.isFinite(number)) throw new Error(`${name} 必须是有限数。`);
  return number;
};
const knownConvection = value => typeof value === 'string'
  && Object.prototype.hasOwnProperty.call(FLOW_CONVECTION_SCHEMES, value);
const knownPressurePreconditioner = value => typeof value === 'string'
  && Object.prototype.hasOwnProperty.call(FLOW_PRESSURE_PRECONDITIONERS, value);
const knownOutletBackflow = value => typeof value === 'string'
  && Object.prototype.hasOwnProperty.call(FLOW_OUTLET_BACKFLOW_MODES, value);

function validateFlowRequest(request = {}) {
  const flowCase = Object.hasOwn(FLOW_CASES, request.case) ? FLOW_CASES[request.case] : null;
  if (!flowCase) throw new Error('未知流动工况。');
  const convection = request.convection === undefined ? 'upwind' : request.convection;
  if (!knownConvection(convection)) throw new Error('未知对流格式。请选择 upwind、limited-linear 或 face-limited-linear。');
  const pressurePreconditioner = request.pressurePreconditioner === undefined
    ? 'ic0' : request.pressurePreconditioner;
  if (!knownPressurePreconditioner(pressurePreconditioner))
    throw new Error('未知压力预条件器。请选择 ic0 或 aggregation。');
  const outletBackflow = request.outletBackflow === undefined ? 'reject' : request.outletBackflow;
  if (!knownOutletBackflow(outletBackflow))
    throw new Error('未知出口回流处理。请选择一种已支持的处理方式。');
  const nu = finite(request.nu, '运动黏度');
  const speed = finite(request.speed, '参考速度');
  const maxIterations = finite(request.maxIterations, '最大迭代数');
  const tolerance = finite(request.tolerance ?? 1e-6, '流动停止容差');
  if (tolerance < 1e-12 || tolerance > 1e-6) throw new Error('流动停止容差须在 1e-12 到 1e-6 之间；只能加严默认停止条件。');
  if (!(nu > 0)) throw new Error('运动黏度必须大于 0。');
  if (!(speed > 0)) throw new Error('参考速度必须大于 0。');
  if (!Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 100000)
    throw new Error('最大迭代数必须是 1 到 100000 的整数。');
  const viscousStress = request.viscousStress === undefined ? VISCOUS_STRESS : request.viscousStress;
  if (viscousStress !== VISCOUS_STRESS) throw new Error('未知黏性应力格式。当前仅支持 symmetric。');
  const mode = request.mode ?? 'steady';
  if (!['steady', 'transient', 'adaptive'].includes(mode)) throw new Error('未知时间模式。');
  const steadyAcceleration=request.steadyAcceleration ?? 'none';
  if (!['none','anderson'].includes(steadyAcceleration) || (steadyAcceleration!=='none' && mode!=='steady'))
    throw new Error('稳态加速仅支持稳态流动的 none 或 anderson。');
  const linearPolicy=request.linearPolicy ?? 'strict';
  if (!['strict','adaptive'].includes(linearPolicy)) throw new Error('线性迭代精度须为 strict 或 adaptive。');
  const velocityRelaxation=finite(request.velocityRelaxation ?? .6,'速度松弛系数');
  if (!(velocityRelaxation>0 && velocityRelaxation<=1)) throw new Error('速度松弛系数须大于0且不大于1。');
  if (request.resume !== undefined && typeof request.resume !== 'boolean') throw new Error('续算选项无效。');
  if (mode === 'steady' && request.resume) throw new Error('稳态模式不能读取非定常重启状态。');
  const normalized = { case: flowCase.id, nu, speed, maxIterations, tolerance, convection,
    pressurePreconditioner, outletBackflow, viscousStress, steadyAcceleration, linearPolicy, velocityRelaxation, mode, resume: Boolean(request.resume) };
  if (flowCase.id === 'custom') {
    if (outletBackflow !== 'reject') throw new Error('命名边界目前只支持检测到出口回流时停止。');
    normalized.boundaryDefinition = normalizeBoundaryDefinition(request.boundaryDefinition);
  }
  if (mode !== 'steady') {
    normalized.dt = finite(request.dt, '时间步长');
    if (!(normalized.dt > 0)) throw new Error('时间步长必须大于 0。');
    if (mode === 'transient') {
      normalized.steps = finite(request.steps, '本次时间步数');
      if (!Number.isInteger(normalized.steps) || normalized.steps < 1 || normalized.steps > 1000000)
        throw new Error('本次时间步数必须是 1 到 1000000 的整数。');
      if (!Number.isFinite(normalized.dt * normalized.steps)) throw new Error('物理时间超出数值范围。');
    } else {
      normalized.endTime=finite(request.endTime,'目标物理时间');
      normalized.minDt=finite(request.minDt ?? normalized.dt/1024,'最小时间步长');
      normalized.maxCourant=finite(request.maxCourant ?? 1,'目标 CFL');
      normalized.maxRetries=finite(request.maxRetries ?? 10,'最大重试次数');
      normalized.maxSteps=finite(request.maxSteps ?? 100000,'最大接受步数');
      if (!(normalized.endTime>0) || !(normalized.minDt>0) || normalized.minDt>normalized.dt || !(normalized.maxCourant>0)
          || !Number.isInteger(normalized.maxRetries) || normalized.maxRetries<0 || normalized.maxRetries>30
          || !Number.isInteger(normalized.maxSteps) || normalized.maxSteps<1 || normalized.maxSteps>1000000)
        throw new Error('自动步长的目标时间、步长范围、CFL 或计算预算无效。');
    }
  }
  if (request.initialVortex!==undefined) {
    if (mode==='steady' || normalized.resume) throw new Error('初始局部涡仅用于新的非定常计算，续算不可重复施加。');
    normalized.initialVortex=normalizeInitialVortex(request.initialVortex);
  }
  return normalized;
}

function buildFlowInvocation(meshPath, outputPrefix, request, restartPath = null, boundaryPath = null) {
  if (typeof meshPath !== 'string' || !meshPath.endsWith('.solver.cm2d'))
    throw new Error('流动求解只能读取本次最终 solver.cm2d。');
  const validated = validateFlowRequest(request);
  if (validated.resume && !restartPath) throw new Error('请先选择重启状态。');
  const temporalArgs = validated.mode === 'transient'
    ? ['--time-step', String(validated.dt), '--steps', String(validated.steps)]
    : validated.mode === 'adaptive' ? ['--time-step',String(validated.dt),'--end-time',String(validated.endTime),
      '--min-time-step',String(validated.minDt),'--max-courant',String(validated.maxCourant),
      '--max-step-retries',String(validated.maxRetries),'--max-time-steps',String(validated.maxSteps)] : [];
  if (validated.resume) temporalArgs.push('--restart', restartPath);
  if (validated.steadyAcceleration!=='none') temporalArgs.push('--steady-acceleration',validated.steadyAcceleration);
  if (validated.linearPolicy!=='strict') temporalArgs.push('--linear-policy',validated.linearPolicy);
  if (validated.velocityRelaxation!==.6) temporalArgs.push('--velocity-relaxation',String(validated.velocityRelaxation));
  if (validated.initialVortex) {
    const v=validated.initialVortex;
    temporalArgs.push('--initial-vortex-x',String(v.centre[0]),'--initial-vortex-y',String(v.centre[1]),
      '--initial-vortex-radius',String(v.radius),'--initial-vortex-speed',String(v.peakSpeed));
  }
  if (validated.case === 'custom') {
    if (typeof boundaryPath !== 'string' || !boundaryPath) throw new Error('缺少命名边界输入路径。');
    temporalArgs.push('--boundary', boundaryPath);
  }
  return {
    executable: 'cartmesh2d_flow_cli',
    request: validated,
    args: ['--mesh', meshPath, '--output', outputPrefix,
      '--case', validated.case, '--nu', String(validated.nu),
      '--speed', String(validated.speed), '--max-iterations', String(validated.maxIterations),
      '--tolerance', String(validated.tolerance),
      '--convection', validated.convection, '--pressure-preconditioner', validated.pressurePreconditioner,
      '--outlet-backflow', validated.outletBackflow, '--viscous-stress', validated.viscousStress, ...temporalArgs]
  };
}

function parseFlowProgress(line) {
  const progressNumber = (value, label) => {
    if (typeof value !== 'number' || !Number.isFinite(value)) throw new Error(`${label} 进度必须为有限数值。`);
    return value;
  };
  if (typeof line !== 'string' || line[0] !== '{') return null;
  let value;
  try { value = JSON.parse(line); } catch { return null; }
  if (value?.type === 'flow-time-retry') {
    for (const key of ['step','attempt','acceptedTime','candidateTime','dt','nextDt','maxCourant']) progressNumber(value[key],key);
    if (!Number.isInteger(value.step) || value.step<1 || !Number.isInteger(value.attempt) || value.attempt<1
        || value.acceptedTime<0 || !(value.dt>value.nextDt) || !(value.nextDt>0) || value.maxCourant<0
        || !near(value.candidateTime,value.acceptedTime+value.dt) || !['courant','nonconverged'].includes(value.reason))
      throw new Error('自动步长重试进度无效。');
    return value;
  }
  if (value?.type === 'flow-time-step') {
    const time = progressNumber(value.time, '物理时间');
    const step = progressNumber(value.step, '时间步');
    const maxCourant = progressNumber(value.maxCourant, 'CFL');
    if (!(time > 0) || !Number.isInteger(step) || step < 1 || maxCourant < 0) throw new Error('时间进度无效。');
    const result = { type: value.type, time, step, maxCourant };
    if (value.kineticEnergy !== undefined || value.forceX !== undefined || value.forceY !== undefined) {
      for (const key of ['kineticEnergy','forceX','forceY']) result[key] = progressNumber(value[key], key);
      if (result.kineticEnergy < 0) throw new Error('动能不能为负。');
    }
    return result;
  }
  if (value?.type !== 'flow-progress') return null;
  const iteration = progressNumber(value.iteration, '进度 iteration');
  if (!Number.isInteger(iteration) || iteration < 0) throw new Error('进度 iteration 无效。');
  const progress = { type: value.type, iteration,
    continuity: progressNumber(value.continuity, '进度 continuity'),
    velocityChange: progressNumber(value.velocityChange, '进度 velocityChange'),
    momentumResidual: progressNumber(value.momentumResidual, '进度 momentumResidual') };
  if (value.time !== undefined || value.timeStep !== undefined) {
    progress.time = progressNumber(value.time, '物理时间');
    progress.timeStep = progressNumber(value.timeStep, '时间步');
    if (!(progress.time > 0) || !Number.isInteger(progress.timeStep) || progress.timeStep < 1)
      throw new Error('物理时间进度无效。');
  }
  if ([progress.continuity, progress.velocityChange, progress.momentumResidual].some(v => v < 0))
    throw new Error('进度残差不能为负。');
  return progress;
}

const near = (a, b) => Math.abs(a - b) <= 1e-12 + 1e-9 * Math.max(Math.abs(a), Math.abs(b));
function flowOutputSuffixes(request) {
  return [...FLOW_OUTPUT_SUFFIXES, ...(['transient','adaptive'].includes(request?.mode) ? ['.checkpoint', '.time-history.csv'] : []),
    ...(request?.mode === 'adaptive' ? ['.attempt-history.csv'] : []),
    ...(request?.initialVortex ? ['.initial.checkpoint'] : []),
    ...(request?.case === 'custom' ? ['.boundaries'] : [])];
}

function validateTimeHistory(text, summary, startTime = 0) {
  const lines = text.trim().split(/\r?\n/);
  const keys = ['step', 'time', 'dt', 'accepted', 'innerIterations', 'momentumResidual', 'continuity', 'maxCourant', 'kineticEnergy', 'forceX', 'forceY'];
  if (lines.shift() !== keys.join(',')) throw new Error('时间历史表头无效。');
  const adaptive=summary.timeStepControl==='adaptive-cfl-retry';
  let lastTime = startTime;
  const rows = lines.map((line, i) => {
    const values = line.split(',');
    if (values.length !== keys.length) throw new Error('时间历史列数不匹配。');
    const row = Object.fromEntries(keys.map((k,j) => [k, finite(values[j], k)]));
    if (row.step !== i + 1 || !(row.time > lastTime) || !near(row.time, lastTime + row.dt) || !(row.dt>0)
        || (!adaptive && !near(row.dt, summary.dt)) || ![0,1].includes(row.accepted)
        || !Number.isInteger(row.innerIterations) || row.innerIterations < 1)
      throw new Error('时间历史顺序或步长不一致。');
    if (row.accepted === 0 && i !== lines.length - 1) throw new Error('失败时间步必须是最后一步。');
    if (['momentumResidual','continuity','maxCourant','kineticEnergy'].some(k => row[k] < 0))
      throw new Error('时间历史包含负的误差或物理指标。');
    if (row.accepted && (row.innerIterations < 10 || row.momentumResidual >= summary.tolerance || row.continuity >= 1e-8))
      throw new Error('时间步未达到停止条件却被接受。');
    lastTime = row.time;
    return row;
  });
  const last = rows.at(-1);
  if (!last || rows.filter(r => r.accepted).length !== summary.completedSteps
      || rows.length !== summary.completedSteps + (summary.converged ? 0 : 1)
      || Boolean(last.accepted) !== summary.converged || !near(last.time, summary.time)
      || last.innerIterations !== summary.iterations || !near(last.dt,summary.dt))
    throw new Error('时间历史与最终状态不一致。');
  for (const name of ['momentumResidual','continuity','maxCourant','forceX','forceY'])
    if (!near(last[name], summary[name])) throw new Error(`时间历史 ${name} 与摘要不一致。`);
  return rows;
}

function validateFlowOutput(summary, fields, expectedCells, expectedRequest = null, startTime = 0) {
  if (!summary || summary.format !== 'cartmesh2d-flow-summary-v1')
    throw new Error('流动摘要格式无效。');
  if (!Object.hasOwn(FLOW_CASES, summary.case)) throw new Error('流动摘要工况无效。');
  if (summary.case === 'custom') {
    const entries = conditions(summary.boundaryConditions);
    const reference = entries.some(b=>b.type==='pressure-opening')
      ? 'explicit pressure opening faces, prescribed static kinematic pressure' : entries.some(b=>b.type==='pressure-outlet')
      ? 'explicit pressure outlet faces, prescribed kinematic pressure' : 'cell 0, kinematic pressure zero';
    if (summary.boundaryFileSuffix !== '.boundaries' || summary.referenceSpeedRole !== 'normalization-only' ||
        summary.outletBackflow !== 'reject' || summary.pressureReference !== reference)
      throw new Error('命名边界结果的压力基准或输入信息不匹配。');
  }
  const transient = summary.temporalDiscretization !== undefined;
  const adaptive=summary.timeStepControl==='adaptive-cfl-retry';
  if (summary.timeStepControl!==undefined && (!adaptive || !transient)) throw new Error('未知时间步控制模式。');
  if (summary.status !== 'converged' && summary.status !== (transient ? 'time_step_not_converged' : 'iteration_limit'))
    throw new Error('流动摘要状态无效。');
  if (typeof summary.converged !== 'boolean' || summary.converged !== (summary.status === 'converged'))
    throw new Error('流动摘要的收敛状态互相矛盾。');
  const convectionInferred = summary.convection === undefined;
  const convection = convectionInferred ? 'upwind' : summary.convection;
  if (!knownConvection(convection)) throw new Error('流动摘要对流格式无效。');
  const pressurePreconditionerInferred = summary.pressurePreconditioner === undefined;
  const pressurePreconditioner = pressurePreconditionerInferred
    ? LEGACY_PRESSURE_PRECONDITIONER : summary.pressurePreconditioner;
  if (!pressurePreconditionerInferred && !knownPressurePreconditioner(pressurePreconditioner))
    throw new Error('流动摘要压力预条件器无效。');
  const outletBackflowInferred = summary.outletBackflow === undefined;
  const outletBackflow = outletBackflowInferred ? 'reject' : summary.outletBackflow;
  if (!knownOutletBackflow(outletBackflow)) throw new Error('流动摘要出口回流处理无效。');
  const pressureDiscretizationInferred = summary.pressureDiscretization === undefined;
  const pressureDiscretization = pressureDiscretizationInferred
    ? LEGACY_PRESSURE_DISCRETIZATION : summary.pressureDiscretization;
  if (!pressureDiscretizationInferred && pressureDiscretization !== PRESSURE_DISCRETIZATION)
    throw new Error('流动摘要压力离散格式无效。');
  const viscousStressInferred = summary.viscousStress === undefined;
  const viscousStress = viscousStressInferred ? LEGACY_VISCOUS_STRESS : summary.viscousStress;
  if (!viscousStressInferred && viscousStress !== VISCOUS_STRESS)
    throw new Error('流动摘要黏性应力格式无效。');
  const iterations = finite(summary.iterations, 'iterations');
  const cells = finite(summary.cells, 'cells');
  if (!Number.isInteger(iterations) || iterations < 1) throw new Error('iterations 无效。');
  if (!Number.isInteger(cells) || cells !== expectedCells) throw new Error('流场单元数与最终网格不一致。');
  const normalizedSummary = {
    ...summary,
    convection,
    convectionInferred,
    pressurePreconditioner,
    pressurePreconditionerInferred,
    outletBackflow,
    outletBackflowInferred,
    pressureDiscretization,
    pressureDiscretizationInferred,
    viscousStress,
    viscousStressInferred,
    nu: finite(summary.nu, 'nu'), speed: finite(summary.speed, 'speed'),
    iterations, cells,
    continuity: finite(summary.continuity, 'continuity'),
    globalImbalance: finite(summary.globalImbalance, 'globalImbalance'),
    globalRelativeImbalance: finite(summary.globalRelativeImbalance, 'globalRelativeImbalance'),
    velocityChange: finite(summary.velocityChange, 'velocityChange'),
    pressureChange: finite(summary.pressureChange, 'pressureChange'),
    momentumResidual: finite(summary.momentumResidual, 'momentumResidual'),
    tolerance: finite(summary.tolerance, 'tolerance')
  };
  if (summary.outletBackflowFaces !== undefined) {
    if (!Number.isSafeInteger(summary.outletBackflowFaces) || summary.outletBackflowFaces < 0)
      throw new Error('流动摘要回流出口面数无效。');
    normalizedSummary.outletBackflowFaces = summary.outletBackflowFaces;
  }
  if (summary.outletInflow !== undefined) {
    normalizedSummary.outletInflow = finite(summary.outletInflow, 'outletInflow');
    if (normalizedSummary.outletInflow < 0) throw new Error('流动摘要出口流入量不能为负。');
  }
  if (!(normalizedSummary.nu > 0) || !(normalizedSummary.speed > 0))
    throw new Error('流动摘要物性或参考速度无效。');
  if (!(normalizedSummary.tolerance > 0)) throw new Error('流动停止条件无效。');
  for (const key of ['continuity', 'globalRelativeImbalance', 'velocityChange', 'pressureChange', 'momentumResidual'])
    if (normalizedSummary[key] < 0) throw new Error('流动误差指标不能为负。');
  if (summary.converged && (iterations < 10 || normalizedSummary.continuity >= 1e-8
      || normalizedSummary.globalRelativeImbalance >= 1e-8
      || ['velocityChange', 'pressureChange', 'momentumResidual'].some(key => normalizedSummary[key] >= normalizedSummary.tolerance)))
    throw new Error('摘要声称收敛，但实际指标未达到停止条件。');

  if (transient) {
    if (summary.temporalDiscretization !== 'backward-euler' || summary.temporalFaceInterpolation !== 'old-and-iteration-flux-defect-skew-corrected-v2')
      throw new Error('未知非定常离散格式。');
    for (const key of ['time','dt','acceptedTime','completedSteps','maxCourant',...(adaptive?[]:['requestedSteps'])])
      normalizedSummary[key] = finite(summary[key], key);
    const q = normalizedSummary;
    if (adaptive) validateAdaptiveSummary(q,startTime);
    if (!(q.dt > 0) || q.maxCourant < 0 || q.acceptedTime < 0 || q.time <= startTime
        || (!adaptive && (!Number.isInteger(q.requestedSteps) || q.requestedSteps < 1
        || !Number.isInteger(q.completedSteps) || q.completedSteps < 0 || q.completedSteps > q.requestedSteps
        || (q.converged && q.completedSteps !== q.requestedSteps)
        || (!q.converged && q.completedSteps >= q.requestedSteps)
        || !near(q.acceptedTime, startTime + q.completedSteps*q.dt)
        || !near(q.time, q.acceptedTime + (q.converged ? 0 : q.dt)))))
      throw new Error('非定常时间、接受状态或步数不一致。');
  } else if (['dt','time','acceptedTime','completedSteps','requestedSteps'].some(k => summary[k] !== undefined)) {
    throw new Error('缺少非定常离散格式，不能解释时间字段。');
  }

  const forceKeys = ['forceX', 'forceY', 'pressureForceX', 'pressureForceY',
    'discreteForceX', 'discreteForceY', 'wallForceX', 'wallForceY',
    'wallViscousForceX', 'wallViscousForceY'];
  if (viscousStress === VISCOUS_STRESS) {
    if (summary.forceDefinition !== FORCE_DEFINITION)
      throw new Error('对称黏性应力结果的受力定义不匹配。');
    for (const key of forceKeys) {
      if (summary[key] !== undefined) normalizedSummary[key] = finite(summary[key], key);
    }
    for (const key of forceKeys)
      if (normalizedSummary[key] === undefined) throw new Error(`对称黏性应力结果缺少 ${key}。`);
    const forceTolerance = 1e-10 * Math.max(1, Math.abs(normalizedSummary.forceX), Math.abs(normalizedSummary.forceY),
      Math.abs(normalizedSummary.discreteForceX), Math.abs(normalizedSummary.discreteForceY));
    if (Math.abs(normalizedSummary.forceX - normalizedSummary.discreteForceX) > forceTolerance
        || Math.abs(normalizedSummary.forceY - normalizedSummary.discreteForceY) > forceTolerance)
      throw new Error('总力与共享面离散力不一致。');
  }

  validateWallLoads(normalizedSummary);
  validateBoundaryFluxes(normalizedSummary);
  const acceleration=summary.steadyAcceleration ?? 'none';
  if (summary.steadyAcceleration===undefined && ['accelerationCandidates','accelerationAccepted','accelerationRejected'].some(key=>Object.hasOwn(summary,key)))
    throw new Error('稳态加速统计缺少模式。');
  if (!['none','anderson'].includes(acceleration) || (acceleration!=='none' && transient))
    throw new Error('稳态加速结果模式无效。');
  normalizedSummary.steadyAcceleration=acceleration;
  if (summary.steadyAcceleration!==undefined) {
    for (const key of ['accelerationCandidates','accelerationAccepted','accelerationRejected'])
      if (!Number.isInteger(summary[key]) || summary[key]<0 || summary[key]>iterations)
        throw new Error('稳态加速统计无效。');
    if (summary.accelerationCandidates!==summary.accelerationAccepted+summary.accelerationRejected ||
        (acceleration==='none' && summary.accelerationCandidates!==0)) throw new Error('稳态加速统计不一致。');
  }
  validateInitialVortexOutput(summary,null,startTime);

  if (summary.adaptiveLinear!==undefined && typeof summary.adaptiveLinear!=='boolean')
    throw new Error('原生线性迭代精度模式无效。');
  if (summary.strictLinearFinal!==undefined && typeof summary.strictLinearFinal!=='boolean')
    throw new Error('原生严格线性复核标记无效。');
  const linearPolicy=summary.adaptiveLinear===true?'adaptive':'strict';
  const velocityRelaxation=summary.velocityRelaxation===undefined?.6:finite(summary.velocityRelaxation,'速度松弛系数');
  if (!(velocityRelaxation>0 && velocityRelaxation<=1)) throw new Error('原生速度松弛系数无效。');
  if (linearPolicy==='adaptive' && normalizedSummary.converged && summary.strictLinearFinal!==true)
    throw new Error('自适应线性结果缺少严格收敛复核。');
  if (linearPolicy==='adaptive' && transient && (!Number.isInteger(summary.strictAcceptedSteps)
      || summary.strictAcceptedSteps!==normalizedSummary.completedSteps))
    throw new Error('自适应线性时间推进缺少逐步严格复核。');
  normalizedSummary.linearPolicy=linearPolicy;
  normalizedSummary.velocityRelaxation=velocityRelaxation;

  if (expectedRequest) {
    const request = validateFlowRequest(expectedRequest);
    validateInitialVortexOutput(summary,request,startTime);
    if (request.case === 'custom' && !sameConditions(request.boundaryDefinition.records, summary.boundaryConditions))
      throw new Error('求解结果的命名边界与请求不一致。');
    if ((request.mode !== 'steady') !== transient || (request.mode==='adaptive')!==adaptive
        || (request.mode==='transient' && (!near(normalizedSummary.dt,request.dt) || normalizedSummary.requestedSteps !== request.steps)))
      throw new Error('原生求解时间模式或步长与请求不一致。');
    if (adaptive) validateAdaptiveSummary(normalizedSummary,startTime,request);
    if (pressureDiscretizationInferred)
      throw new Error('本次新流动结果缺少压力离散格式，不能与请求绑定。');
    if (viscousStressInferred || viscousStress !== request.viscousStress)
      throw new Error('本次新流动结果缺少对称黏性应力格式，不能与请求绑定。');
    if (pressurePreconditionerInferred)
      throw new Error('本次新流动结果缺少压力预条件器，不能与请求绑定。');
    if (normalizedSummary.case !== request.case || normalizedSummary.nu !== request.nu
        || normalizedSummary.speed !== request.speed || normalizedSummary.convection !== request.convection
        || normalizedSummary.pressurePreconditioner !== request.pressurePreconditioner
        || normalizedSummary.outletBackflow !== request.outletBackflow
        || normalizedSummary.tolerance !== request.tolerance
        || normalizedSummary.steadyAcceleration !== request.steadyAcceleration
        || normalizedSummary.linearPolicy !== request.linearPolicy
        || normalizedSummary.velocityRelaxation !== request.velocityRelaxation
        || normalizedSummary.iterations > request.maxIterations)
      throw new Error('原生求解结果与请求工况不一致。');
  }

  if (!fields || fields.format !== 'cartmesh2d-flow-v1' || !Array.isArray(fields.cells))
    throw new Error('流场文件格式无效。');
  if (fields.cells.length !== expectedCells) throw new Error('流场没有覆盖最终网格的全部单元。');
  const normalizedCells = new Array(expectedCells);
  for (const cell of fields.cells) {
    if (!Number.isInteger(cell?.id) || cell.id < 0 || cell.id >= expectedCells || normalizedCells[cell.id])
      throw new Error('流场单元 ID 必须从 0 连续且不重复。');
    const u = finite(cell.u, `cell ${cell.id} u`);
    const v = finite(cell.v, `cell ${cell.id} v`);
    const p = finite(cell.p, `cell ${cell.id} p`);
    const speed = finite(cell.speed, `cell ${cell.id} speed`);
    if (speed < 0 || Math.abs(speed - Math.hypot(u, v)) > 1e-8 * Math.max(1, speed))
      throw new Error(`cell ${cell.id} 的 speed 与速度分量不一致。`);
    normalizedCells[cell.id] = { id: cell.id, u, v, p, speed };
  }
  if (normalizedCells.some(cell => !cell)) throw new Error('流场单元 ID 不连续。');
  return { summary: normalizedSummary, fields: { ...fields, cells: normalizedCells } };
}

async function commitFlowFiles(fileSystem, entries) {
  try {
    for (const entry of entries) await fileSystem.copyFile(entry.source, entry.destination);
  } catch (error) {
    await Promise.allSettled(entries.map(entry => fileSystem.rm(entry.destination, { force: true })));
    throw error;
  }
}

module.exports = {
  FLOW_CASES, FLOW_CONVECTION_SCHEMES, FLOW_PRESSURE_PRECONDITIONERS, FLOW_OUTLET_BACKFLOW_MODES, LEGACY_PRESSURE_PRECONDITIONER, FLOW_OUTPUT_SUFFIXES,
  LEGACY_PRESSURE_DISCRETIZATION, PRESSURE_DISCRETIZATION,
  VISCOUS_STRESS, LEGACY_VISCOUS_STRESS, FORCE_DEFINITION,
  buildFlowInvocation, commitFlowFiles, parseFlowProgress, validateFlowOutput, validateFlowRequest,
  flowOutputSuffixes, validateTimeHistory, validateAttemptHistory
};
