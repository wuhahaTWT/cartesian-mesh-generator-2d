'use strict';

const CAPABILITIES = Object.freeze({
  schemaVersion: 1,
  methods: [
    {
      id: 'cutcell',
      label: 'Pure Cut-cell',
      status: 'stable',
      statusLabel: '稳定',
      description: '通用外流默认方法；物面局部加密，远场保持稀疏。',
      supports: { distanceBands: true, refinementBoxes: true, boundaryLayers: false, openFoam: true },
      defaults: {
        maxLevel: 6,
        minimumLevel: 0,
        paddingFraction: 0.25,
        smallAlpha: 0.05,
        cellBudget: 250000
      },
      presets: {
        quick: { maxLevel: 5, minimumLevel: 0 },
        standard: { maxLevel: 6, minimumLevel: 0 },
        dense: { maxLevel: 8, minimumLevel: 0 }
      }
    },
    {
      id: 'hybrid',
      label: 'Hybrid + Boundary Layer',
      status: 'beta',
      statusLabel: 'Beta',
      description: '共形壁面层与 Cartesian/Cut-cell 余域；高层级能力仍受 W2 边界约束。',
      supports: { distanceBands: false, refinementBoxes: false, boundaryLayers: true, openFoam: false },
      defaults: {
        maxLevel: 6,
        minimumLevel: 3,
        boundaryLevel: 6,
        nLayers: 4,
        firstThickness: 0.02,
        growthRatio: 1.2,
        domainPadding: 1.0,
        extrusionThickness: 0.01,
        cellBudget: 250000
      },
      presets: {
        quick: { maxLevel: 5, minimumLevel: 2, boundaryLevel: 5 },
        standard: { maxLevel: 6, minimumLevel: 3, boundaryLevel: 6 },
        dense: { maxLevel: 8, minimumLevel: 3, boundaryLevel: 8 }
      }
    }
  ]
});

function methodById(methodId) {
  return CAPABILITIES.methods.find(method => method.id === methodId);
}

function finiteNumber(value, name, { min = -Infinity, max = Infinity, integer = false } = {}) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed < min || parsed > max ||
      (integer && !Number.isInteger(parsed))) {
    throw new Error(`${name} 超出允许范围。`);
  }
  return parsed;
}

function estimateMeshJob(job) {
  const minimumLevel = Number(job.minimumLevel);
  const maxLevel = Number(job.maxLevel);
  const globalLeafFloor = Number.isInteger(minimumLevel) && minimumLevel >= 0
    ? Math.pow(4, minimumLevel)
    : null;
  const boundaryScale = Number.isInteger(maxLevel) && maxLevel >= 0
    ? Math.pow(2, maxLevel)
    : null;
  return {
    kind: 'lower-bound',
    globalLeafFloor,
    boundaryScale,
    note: '这是全域底格下限；物面、距离带、尾迹区和壁面层会继续增加单元。'
  };
}

function validateMeshJob(request) {
  const method = methodById(request.method);
  if (!method) throw new Error('未知网格方法。');
  const sourceUnits = request.sourceUnits || 'auto';
  if (!['auto', 'mm', 'cm', 'm', 'in', 'ft'].includes(sourceUnits)) {
    throw new Error('未知 DXF 源单位。');
  }
  const job = {
    method: method.id,
    dxfPath: request.dxfPath,
    outputDirectory: request.outputDirectory,
    sourceUnits,
    chordError: finiteNumber(request.chordError, '曲线弦高误差', { min: 1e-12 }),
    maxLevel: finiteNumber(request.maxLevel, '物面最高层级', { min: 1, max: 28, integer: true }),
    minimumLevel: finiteNumber(request.minimumLevel, '全域最低层级', { min: 0, max: 28, integer: true }),
    cellBudget: finiteNumber(request.cellBudget, '单元预算', { min: 100, max: 100000000, integer: true })
  };
  if (!job.dxfPath || !job.outputDirectory) throw new Error('请选择 DXF 文件和输出目录。');
  if (job.minimumLevel > job.maxLevel) throw new Error('全域最低层级不能高于物面最高层级。');

  if (method.id === 'cutcell') {
    job.paddingFraction = finiteNumber(request.paddingFraction, '计算域留白比例', { min: 0.01, max: 10 });
    job.smallAlpha = finiteNumber(request.smallAlpha, '小单元阈值', { min: 0.001, max: 0.999 });
    job.distanceBands = (request.distanceBands || []).map((band, index) => ({
      distance: finiteNumber(band.distance, `距离带 ${index + 1} 的距离`, { min: 1e-12 }),
      targetLevel: finiteNumber(band.targetLevel, `距离带 ${index + 1} 的层级`,
        { min: 0, max: job.maxLevel, integer: true })
    }));
    job.refinementBoxes = (request.refinementBoxes || []).map((box, index) => {
      const normalized = {
        xmin: finiteNumber(box.xmin, `区域 ${index + 1} xmin`),
        ymin: finiteNumber(box.ymin, `区域 ${index + 1} ymin`),
        xmax: finiteNumber(box.xmax, `区域 ${index + 1} xmax`),
        ymax: finiteNumber(box.ymax, `区域 ${index + 1} ymax`),
        targetLevel: finiteNumber(box.targetLevel, `区域 ${index + 1} 的层级`,
          { min: 0, max: job.maxLevel, integer: true })
      };
      if (!(normalized.xmax > normalized.xmin) || !(normalized.ymax > normalized.ymin)) {
        throw new Error(`区域 ${index + 1} 必须满足 xmax > xmin 且 ymax > ymin。`);
      }
      return normalized;
    });
  } else {
    job.boundaryLevel = finiteNumber(request.boundaryLevel, '壁面层连接层级',
      { min: job.minimumLevel, max: job.maxLevel, integer: true });
    job.nLayers = finiteNumber(request.nLayers, '壁面层数', { min: 1, max: 64, integer: true });
    job.firstThickness = finiteNumber(request.firstThickness, '首层厚度', { min: 1e-12 });
    job.growthRatio = finiteNumber(request.growthRatio, '增长率', { min: 1, max: 5 });
    job.domainPadding = finiteNumber(request.domainPadding, '计算域留白距离', { min: 1e-12 });
    job.extrusionThickness = finiteNumber(request.extrusionThickness, 'OpenFOAM 挤出厚度', { min: 1e-12 });
  }

  const estimate = estimateMeshJob(job);
  if (estimate.globalLeafFloor > job.cellBudget) {
    throw new Error(`仅全域底格下限就约有 ${estimate.globalLeafFloor.toLocaleString()} 个叶格，已超过 ${job.cellBudget.toLocaleString()} 的单元预算。请降低全域最低层级或提高预算。`);
  }
  return { job, method, estimate };
}

function buildMeshInvocation(job, paths) {
  if (job.method === 'cutcell') {
    const args = [paths.xyPath, paths.prefix, String(job.maxLevel),
      String(job.paddingFraction), String(job.smallAlpha), 'exterior', paths.casePath,
      String(job.minimumLevel), '0'];
    for (const band of job.distanceBands) {
      args.push('--distance-band', String(band.distance), String(band.targetLevel));
    }
    for (const box of job.refinementBoxes) {
      args.push('--refine-box', String(box.xmin), String(box.ymin), String(box.xmax),
        String(box.ymax), String(box.targetLevel));
    }
    return {
      executableName: 'cartmesh2d_cli',
      args,
      progress: '正在生成 Pure Cut-cell 网格…',
      cm2dCandidates: [`${paths.prefix}.cm2d`]
    };
  }
  return {
    executableName: 'cartmesh2d_hybrid_cli',
    args: [paths.xyPath, paths.prefix, String(job.maxLevel), String(job.minimumLevel),
      String(job.boundaryLevel), String(job.nLayers), String(job.firstThickness),
      String(job.growthRatio), String(job.domainPadding)],
    progress: '正在生成 Hybrid + Boundary Layer 网格（Beta）…',
    cm2dCandidates: [
      `${paths.prefix}.hybrid.solver.cm2d`,
      `${paths.prefix}.fallback.solver.cm2d`
    ]
  };
}

function normalizeSummary(job, raw) {
  if (job.method === 'cutcell') {
    return { ...raw, requested_method: 'cutcell', actual_method: 'cutcell' };
  }
  const fallback = raw.mesh_mode === 'pure_cutcell_fallback';
  const contractPass = fallback || raw.quality_contract === 'PASS';
  return {
    ...raw,
    requested_method: 'hybrid',
    actual_method: fallback ? 'cutcell-fallback' : 'hybrid',
    stabilized_cells: raw.solver_cells || raw.cells,
    openfoam_cells: raw.openfoam === 'written' ? (raw.solver_cells || raw.cells) : null,
    quality_status: contractPass ? (raw.solver_quality || 'unknown') : 'contract fail',
    quality_pass: contractPass && String(raw.solver_quality || '').toLowerCase() === 'pass'
  };
}

module.exports = {
  getCapabilities: () => CAPABILITIES,
  estimateMeshJob,
  validateMeshJob,
  buildMeshInvocation,
  normalizeSummary
};
