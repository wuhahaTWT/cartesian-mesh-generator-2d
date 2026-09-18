'use strict';

const FLOW_CASES = Object.freeze({
  external: {
    id: 'external', label: '外流',
    scope: '矩形域：左侧恒速入口、右侧压力出口、上下滑移，物面无滑移。'
  },
  channel: {
    id: 'channel', label: '通道',
    scope: '矩形内流：左侧抛物线入口、右侧 p=0，其余边界无滑移。'
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
  }
});
const PRESSURE_DISCRETIZATION = 'shared-face-gauss';
const LEGACY_PRESSURE_DISCRETIZATION = 'legacy-unspecified';
const VISCOUS_STRESS = 'symmetric';
const LEGACY_VISCOUS_STRESS = 'laplacian';
const FORCE_DEFINITION = 'shared-face-newtonian-traction';
const FLOW_OUTPUT_SUFFIXES = Object.freeze([
  '.json', '.fields.json', '.vtk', '.residuals.csv', '.cells.csv', '.faces.csv'
]);

const finite = (value, name) => {
  if (value === null || value === '') throw new Error(`${name} 必须是有限数。`);
  const number = Number(value);
  if (!Number.isFinite(number)) throw new Error(`${name} 必须是有限数。`);
  return number;
};
const knownConvection = value => typeof value === 'string'
  && Object.prototype.hasOwnProperty.call(FLOW_CONVECTION_SCHEMES, value);

function validateFlowRequest(request = {}) {
  const flowCase = FLOW_CASES[request.case];
  if (!flowCase) throw new Error('未知流动工况。');
  const convection = request.convection === undefined ? 'upwind' : request.convection;
  if (!knownConvection(convection)) throw new Error('未知对流格式。请选择 upwind 或 limited-linear。');
  const nu = finite(request.nu, '运动黏度');
  const speed = finite(request.speed, '参考速度');
  const maxIterations = finite(request.maxIterations, '最大迭代数');
  if (!(nu > 0)) throw new Error('运动黏度必须大于 0。');
  if (!(speed > 0)) throw new Error('参考速度必须大于 0。');
  if (!Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 100000)
    throw new Error('最大迭代数必须是 1 到 100000 的整数。');
  const viscousStress = request.viscousStress === undefined ? VISCOUS_STRESS : request.viscousStress;
  if (viscousStress !== VISCOUS_STRESS) throw new Error('未知黏性应力格式。当前仅支持 symmetric。');
  return { case: flowCase.id, nu, speed, maxIterations, convection, viscousStress };
}

function buildFlowInvocation(meshPath, outputPrefix, request) {
  if (typeof meshPath !== 'string' || !meshPath.endsWith('.solver.cm2d'))
    throw new Error('流动求解只能读取本次最终 solver.cm2d。');
  const validated = validateFlowRequest(request);
  return {
    executable: 'cartmesh2d_flow_cli',
    request: validated,
    args: ['--mesh', meshPath, '--output', outputPrefix,
      '--case', validated.case, '--nu', String(validated.nu),
      '--speed', String(validated.speed), '--max-iterations', String(validated.maxIterations),
      '--convection', validated.convection, '--viscous-stress', validated.viscousStress]
  };
}

function parseFlowProgress(line) {
  if (typeof line !== 'string' || line[0] !== '{') return null;
  let value;
  try { value = JSON.parse(line); } catch { return null; }
  if (value?.type !== 'flow-progress') return null;
  const iteration = finite(value.iteration, '进度 iteration');
  if (!Number.isInteger(iteration) || iteration < 0) throw new Error('进度 iteration 无效。');
  return {
    type: 'flow-progress', iteration,
    continuity: finite(value.continuity, '进度 continuity'),
    velocityChange: finite(value.velocityChange, '进度 velocityChange'),
    momentumResidual: finite(value.momentumResidual, '进度 momentumResidual')
  };
}

function validateFlowOutput(summary, fields, expectedCells, expectedRequest = null) {
  if (!summary || summary.format !== 'cartmesh2d-flow-summary-v1')
    throw new Error('流动摘要格式无效。');
  if (!FLOW_CASES[summary.case]) throw new Error('流动摘要工况无效。');
  if (summary.status !== 'converged' && summary.status !== 'iteration_limit')
    throw new Error('流动摘要状态无效。');
  if (typeof summary.converged !== 'boolean' || summary.converged !== (summary.status === 'converged'))
    throw new Error('流动摘要的收敛状态互相矛盾。');
  const convectionInferred = summary.convection === undefined;
  const convection = convectionInferred ? 'upwind' : summary.convection;
  if (!knownConvection(convection)) throw new Error('流动摘要对流格式无效。');
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
  if (!(normalizedSummary.nu > 0) || !(normalizedSummary.speed > 0))
    throw new Error('流动摘要物性或参考速度无效。');
  if (!(normalizedSummary.tolerance > 0)) throw new Error('流动停止条件无效。');
  for (const key of ['continuity', 'globalRelativeImbalance', 'velocityChange', 'pressureChange', 'momentumResidual'])
    if (normalizedSummary[key] < 0) throw new Error('流动误差指标不能为负。');
  if (summary.converged && (iterations < 10 || normalizedSummary.continuity >= 1e-8
      || normalizedSummary.globalRelativeImbalance >= 1e-8
      || ['velocityChange', 'pressureChange', 'momentumResidual'].some(key => normalizedSummary[key] >= normalizedSummary.tolerance)))
    throw new Error('摘要声称收敛，但实际指标未达到停止条件。');

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

  if (expectedRequest) {
    const request = validateFlowRequest(expectedRequest);
    if (pressureDiscretizationInferred)
      throw new Error('本次新流动结果缺少压力离散格式，不能与请求绑定。');
    if (viscousStressInferred || viscousStress !== request.viscousStress)
      throw new Error('本次新流动结果缺少对称黏性应力格式，不能与请求绑定。');
    if (normalizedSummary.case !== request.case || normalizedSummary.nu !== request.nu
        || normalizedSummary.speed !== request.speed || normalizedSummary.convection !== request.convection
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
  FLOW_CASES, FLOW_CONVECTION_SCHEMES, FLOW_OUTPUT_SUFFIXES,
  LEGACY_PRESSURE_DISCRETIZATION, PRESSURE_DISCRETIZATION,
  VISCOUS_STRESS, LEGACY_VISCOUS_STRESS, FORCE_DEFINITION,
  buildFlowInvocation, commitFlowFiles, parseFlowProgress, validateFlowOutput, validateFlowRequest
};
