'use strict';

const { methodById, SAFE_WALL_LEVEL } = require('./capabilities');
const { describeBudget, wallLevelFor } = require('./sizing');

// A job is validated once, here, and then turned into an argv.  The CLI takes up to
// nine positional arguments whose meaning depends on position, so building that
// string anywhere else is how a frontend silently passes small-alpha as a fluid
// region.

const number = (value, name, { min = -Infinity, max = Infinity, integer = false } = {}) => {
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed < min || parsed > max ||
      (integer && !Number.isInteger(parsed))) {
    throw new Error(`${name} 超出允许范围（${min} … ${max}）。`);
  }
  return parsed;
};

function validateSizeField(request, safeWallLevel) {
  const field = {
    farFieldSpans: number(request.farFieldSpans, '远场距离（体长倍数）', { min: 0.25, max: 1000 }),
    wallCellsPerSpan: number(request.wallCellsPerSpan, '壁面单元数（体长 / n）', { min: 1, max: 1 << 20 }),
    cellsPerLevel: number(request.cellsPerLevel, '每级带宽（单元数）', { min: 0, max: 64, integer: true }),
    farLevel: number(request.farLevel ?? 0, '远场最低层级', { min: 0, max: 28, integer: true }),
    curvatureCellsPerRadius: number(request.curvatureCellsPerRadius ?? 0,
      '曲率细化（每半径单元数）', { min: 0, max: 1024 }),
    gapCells: number(request.gapCells ?? 0, '间隙细化（间隙内单元数）', { min: 0, max: 1024 }),
    allowUnsafeWallLevel: Boolean(request.allowUnsafeWallLevel),
    safeWallLevel: number(request.safeWallLevel ?? safeWallLevel, '安全壁面层级上限',
      { min: 1, max: 28, integer: true })
  };
  if (request.wake) {
    field.wake = {
      angleOfAttackDeg: number(request.wake.angleOfAttackDeg, '来流角（度）', { min: -180, max: 180 }),
      downstreamSpans: number(request.wake.downstreamSpans, '尾迹长度（体长倍数）', { min: 0.1, max: 1000 }),
      halfWidthSpans: number(request.wake.halfWidthSpans, '尾迹半宽（体长倍数）', { min: 0.05, max: 100 }),
      levelsBelowWall: number(request.wake.levelsBelowWall, '尾迹低于壁面的级数', { min: 0, max: 20, integer: true })
    };
  }
  // Hand-placed regions.  Stated in body spans about the body's own centre, because
  // that is the frame every other sizing number already uses and it keeps a region
  // meaningful when the geometry is rescaled.  Depth is given relative to the wall so a
  // region does not silently become the finest thing in the mesh.
  field.refineBoxes = (request.refineBoxes || []).map((box, index) => {
    const name = `加密区 ${index + 1}`;
    const normalized = {
      xmin: number(box.xmin, `${name} 起点 x`, { min: -1000, max: 1000 }),
      xmax: number(box.xmax, `${name} 终点 x`, { min: -1000, max: 1000 }),
      ymin: number(box.ymin, `${name} 下边界 y`, { min: -1000, max: 1000 }),
      ymax: number(box.ymax, `${name} 上边界 y`, { min: -1000, max: 1000 }),
      levelsBelowWall: number(box.levelsBelowWall, `${name} 低于壁面级数`,
        { min: 0, max: 20, integer: true })
    };
    if (!(normalized.xmax > normalized.xmin) || !(normalized.ymax > normalized.ymin)) {
      throw new Error(`${name} 需要满足 x 终点 > 起点、y 上边界 > 下边界。`);
    }
    return normalized;
  });

  const budget = describeBudget(field);
  if (!budget.feasible && !field.allowUnsafeWallLevel) {    throw new Error(
      `远场 ${field.farFieldSpans} 倍体长配壁面 1/${field.wallCellsPerSpan} 需要 level ` +
      `${budget.wallLevel}，超过实测安全上限 ${field.safeWallLevel}。` +
      `把远场降到 ${budget.maxFarFieldSpans.toFixed(1)} 倍，或把壁面降到 ` +
      `1/${budget.maxWallCellsPerSpan}，或显式勾选越过上限。`);
  }
  return { field, budget };
}

function validateJob(request) {
  const method = methodById(request.method);
  if (!method) throw new Error('未知网格方法。');
  const fluidRegion = request.fluidRegion === 'interior' ? 'interior' : 'exterior';
  const job = {
    method: method.id,
    fluidRegion,
    geometryPath: request.geometryPath,
    outputDirectory: request.outputDirectory,
    chordError: number(request.chordError ?? 0.001, '曲线弦高误差', { min: 1e-12, max: 1 }),
    sourceUnits: request.sourceUnits || 'auto'
  };
  if (!job.geometryPath) throw new Error('请先选择几何文件或内置样例。');
  if (!job.outputDirectory) throw new Error('请选择输出目录。');

  if (method.id === 'cutcell') {
    job.smallAlpha = number(request.smallAlpha ?? method.defaults.smallAlpha,
      '小单元阈值', { min: 0.001, max: 0.999 });
    const { field, budget } = validateSizeField(request, method.safeWallLevel);
    job.sizeField = field;
    job.budget = budget;
  } else {
    job.maxLevel = number(request.maxLevel, '余域最高层级', { min: 1, max: method.safeWallLevel, integer: true });
    job.minimumLevel = number(request.minimumLevel, '全域最低层级', { min: 0, max: job.maxLevel, integer: true });
    job.boundaryLevel = number(request.boundaryLevel, '壁面连接层级',
      { min: job.minimumLevel, max: job.maxLevel, integer: true });
    job.nLayers = number(request.nLayers, '壁面层数', { min: 1, max: 64, integer: true });
    job.firstThickness = number(request.firstThickness, '首层厚度', { min: 1e-12, max: 1e6 });
    job.growthRatio = number(request.growthRatio, '增长率', { min: 1, max: 5 });
    job.domainPadding = number(request.domainPadding, '计算域留白距离', { min: 1e-9, max: 1e6 });
    job.extrusionThickness = number(request.extrusionThickness ?? method.defaults.extrusionThickness,
      'OpenFOAM 挤出厚度', { min: 1e-12, max: 1e6 });
  }
  return { job, method };
}

// While --size-field is active the CLI ignores positional max-level and
// padding-fraction, but still range-checks them, so they are passed as inert
// in-range placeholders rather than left out.
const SIZE_FIELD_PLACEHOLDER_LEVEL = 8;
const SIZE_FIELD_PLACEHOLDER_PADDING = 0.25;

function sizeFieldArgs(field, { dryRun = false } = {}) {
  const args = [dryRun ? '--size-field-only' : '--size-field',
    '--far-field-spans', String(field.farFieldSpans),
    '--wall-cells-per-span', String(field.wallCellsPerSpan),
    '--cells-per-level', String(field.cellsPerLevel),
    '--far-level', String(field.farLevel),
    '--max-safe-wall-level', String(field.safeWallLevel)];
  if (field.curvatureCellsPerRadius > 0) {
    args.push('--curvature-cells-per-radius', String(field.curvatureCellsPerRadius));
  }
  if (field.gapCells > 0) args.push('--gap-cells', String(field.gapCells));
  if (field.wake) {
    args.push('--wake', String(field.wake.angleOfAttackDeg), String(field.wake.downstreamSpans),
      String(field.wake.halfWidthSpans), String(field.wake.levelsBelowWall));
  }
  if (field.allowUnsafeWallLevel) args.push('--allow-unsafe-wall-level');
  return args;
}

// --refine-box takes absolute coordinates, so the body frame converts.  A box coarser
// than level 1 is not expressible (refine() rejects level 0), and one deeper than the
// wall is refused by the CLI, so both ends are clamped here.
function refineBoxArgs(field, frame) {
  if (!field.refineBoxes || !field.refineBoxes.length || !frame) return [];
  const wallLevel = wallLevelFor(field.farFieldSpans, field.wallCellsPerSpan);
  const args = [];
  for (const box of field.refineBoxes) {
    const level = Math.min(wallLevel, Math.max(1, wallLevel - box.levelsBelowWall));
    args.push('--refine-box',
      String(frame.centreX + box.xmin * frame.bodySpan),
      String(frame.centreY + box.ymin * frame.bodySpan),
      String(frame.centreX + box.xmax * frame.bodySpan),
      String(frame.centreY + box.ymax * frame.bodySpan),
      String(level));
  }
  return args;
}

// An OpenFOAM case directory is always requested.  The CLI only builds the solver
// partition, evaluates solver quality and rates the Q1 contract inside
// `if (openFoamCase)`, so asking for no case means asking for no quality report.
function buildInvocation(job, paths, options = {}) {
  if (job.method === 'cutcell') {
    return {
      executable: 'cartmesh2d_cli',
      args: [paths.xyPath, paths.prefix,
        String(SIZE_FIELD_PLACEHOLDER_LEVEL), String(SIZE_FIELD_PLACEHOLDER_PADDING),
        String(job.smallAlpha), job.fluidRegion,
        options.dryRun ? '-' : paths.casePath,
        '0', '0',
        ...sizeFieldArgs(job.sizeField, options),
        ...(options.dryRun ? [] : refineBoxArgs(job.sizeField, paths.frame))],
      cm2dCandidates: [`${paths.prefix}.cm2d`]
    };
  }
  return {
    executable: 'cartmesh2d_hybrid_cli',
    args: [paths.xyPath, paths.prefix, String(job.maxLevel), String(job.minimumLevel),
      String(job.boundaryLevel), String(job.nLayers), String(job.firstThickness),
      String(job.growthRatio), String(job.domainPadding),
      paths.casePath, String(job.extrusionThickness)],
    cm2dCandidates: [`${paths.prefix}.hybrid.solver.cm2d`,
      `${paths.prefix}.fallback.solver.cm2d`]
  };
}

module.exports = { validateJob, validateSizeField, buildInvocation, sizeFieldArgs,
                   refineBoxArgs, SAFE_WALL_LEVEL, wallLevelFor };
