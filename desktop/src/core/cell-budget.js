'use strict';

const { methodById } = require('./capabilities');

const MIN_TARGET_CELLS = 1000;
const MAX_TARGET_CELLS = 500000;
const TARGET_LOWER_RATIO = 0.7;
const TARGET_UPPER_RATIO = 1.3;
const MAX_ATTEMPTS = 3;
const BAND_RESERVE_FRACTION = 0.5;

const BUDGET_PRESETS = Object.freeze([
  Object.freeze({ value: 5000, label: '快速预览 · 约 5 千' }),
  Object.freeze({ value: 15000, label: '基础检查 · 约 1.5 万' }),
  Object.freeze({ value: 50000, label: '细化检查 · 约 5 万' }),
  Object.freeze({ value: 150000, label: '高密网格 · 约 15 万' }),
  Object.freeze({ value: 300000, label: '规模试验 · 约 30 万' })
]);

function finiteNumber(value, name, { min = -Infinity, max = Infinity, integer = false } = {}) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed < min || parsed > max ||
      (integer && !Number.isInteger(parsed))) {
    throw new Error(`${name} 超出允许范围（${min} … ${max}）。`);
  }
  return parsed;
}

function targetCells(value) {
  return finiteNumber(value, '目标单元数', {
    min: MIN_TARGET_CELLS,
    max: MAX_TARGET_CELLS,
    integer: true
  });
}

function frameExtent(frame, axis) {
  const direct = Number(frame?.[axis]);
  if (Number.isFinite(direct) && direct > 0) return direct;
  const bounds = frame?.bounds;
  const nested = Number(bounds?.[axis]);
  if (Number.isFinite(nested) && nested > 0) return nested;
  const lo = Number(bounds?.[axis === 'width' ? 'minX' : 'minY'] ??
                    frame?.[axis === 'width' ? 'minX' : 'minY']);
  const hi = Number(bounds?.[axis === 'width' ? 'maxX' : 'maxY'] ??
                    frame?.[axis === 'width' ? 'maxX' : 'maxY']);
  return Number.isFinite(lo) && Number.isFinite(hi) && hi > lo ? hi - lo : undefined;
}

function areaFor(request, frame, referenceLength, bodySpan, farFieldSpans) {
  if (frame?.fluidArea !== undefined) {
    return {
      area: finiteNumber(frame.fluidArea, '流体面积', { min: Number.MIN_VALUE }),
      source: 'fluidArea',
      estimated: false
    };
  }

  const width = frameExtent(frame, 'width') ?? bodySpan;
  const height = frameExtent(frame, 'height') ?? bodySpan;
  if (request.fluidRegion === 'interior') {
    return { area: width * height, source: 'boundingBoxEstimate', estimated: true };
  }

  // The native exterior domain is square. Using the longest body span is a
  // conservative fallback when the caller has not measured the true fluid area.
  const domainSpan = bodySpan + 2 * farFieldSpans * referenceLength;
  return { area: domainSpan * domainSpan, source: 'boundingBoxEstimate', estimated: true };
}

function safeMinimumWall(request, method, bodySpan, referenceLength, farFieldSpans) {
  if (request.allowUnsafeWallLevel === true) return 0;
  const safeLevel = finiteNumber(request.safeWallLevel ?? method.safeWallLevel,
    '安全壁面层级上限', { min: 1, max: 28, integer: true });
  const domainSpan = bodySpan + 2 * farFieldSpans * referenceLength;
  return domainSpan / Math.pow(2, safeLevel) / referenceLength;
}

function estimateFromBackground(area, backgroundRelativeSize, referenceLength) {
  const backgroundSize = backgroundRelativeSize * referenceLength;
  const uniformCells = area / (backgroundSize * backgroundSize);
  return Math.round(uniformCells / (1 - BAND_RESERVE_FRACTION));
}

function assessBudget(actual, target) {
  const actualCells = finiteNumber(actual, '实际单元数', { min: 1, integer: true });
  const wanted = targetCells(target);
  const minimumCells = Math.ceil(wanted * TARGET_LOWER_RATIO);
  const maximumCells = Math.floor(wanted * TARGET_UPPER_RATIO);
  const withinRange = actualCells >= minimumCells && actualCells <= maximumCells;
  return {
    actualCells,
    targetCells: wanted,
    minimumCells,
    maximumCells,
    ratio: actualCells / wanted,
    withinRange,
    withinTarget: withinRange,
    isEstimate: true
  };
}

function planBudget(request, frame) {
  if (!request || typeof request !== 'object') throw new Error('缺少网格请求。');
  const wanted = targetCells(request.targetCells);
  const method = methodById(request.method);
  if (!method) throw new Error('未知网格方法。');
  if (request.fluidRegion !== 'interior' && request.fluidRegion !== 'exterior') {
    throw new Error('流体区域必须是 interior 或 exterior。');
  }

  const bodySpan = finiteNumber(frame?.bodySpan, '几何体跨度', { min: Number.MIN_VALUE });
  const referenceLength = request.referenceLength === undefined || request.referenceLength === null ||
      request.referenceLength === ''
    ? bodySpan
    : finiteNumber(request.referenceLength, '参考长度', { min: Number.MIN_VALUE });
  const farFieldSpans = finiteNumber(request.farFieldSpans, '计算域留白 / 参考长度', {
    min: Number.MIN_VALUE,
    max: 1000
  });
  const measuredArea = areaFor(request, frame, referenceLength, bodySpan, farFieldSpans);
  const wallToBackgroundRatio = finiteNumber(request.wallToBackgroundRatio ?? 4,
    '背景 / 壁面尺寸比', { min: 2, max: 8, integer: true });
  if (![2, 4, 8].includes(wallToBackgroundRatio)) {
    throw new Error('背景 / 壁面尺寸比只允许 2、4 或 8。');
  }

  // Reserve half the budget for near-wall refinement. Choose a *coarser*
  // attainable background level: native ceil rounding otherwise turns a small
  // requested-size change into up to four times as many background cells.
  const domainSpan = bodySpan + 2 * farFieldSpans * referenceLength;
  const continuousSize = Math.sqrt(measuredArea.area /
    (wanted * (1 - BAND_RESERVE_FRACTION)));
  const backgroundLevel = Math.max(0, Math.floor(Math.log2(domainSpan / continuousSize)));
  const backgroundSize = domainSpan / Math.pow(2, backgroundLevel) * (1 + 1e-12);
  let backgroundRelativeSize = backgroundSize / referenceLength;
  let wallRelativeSize = backgroundRelativeSize / wallToBackgroundRatio;
  const minimumWallRelativeSize = safeMinimumWall(request, method, bodySpan,
    referenceLength, farFieldSpans);
  if (minimumWallRelativeSize > 1) {
    throw new Error('当前计算域无法在安全壁面层级内规划；请缩小留白，或显式允许更高构造深度。');
  }
  const limitedBySafeWallLevel = wallRelativeSize < minimumWallRelativeSize;
  if (limitedBySafeWallLevel) {
    wallRelativeSize = minimumWallRelativeSize;
    backgroundRelativeSize = wallRelativeSize * wallToBackgroundRatio;
  }

  wallRelativeSize = Math.min(1, Math.max(1e-12, wallRelativeSize));
  backgroundRelativeSize = wallRelativeSize * wallToBackgroundRatio;
  const cellsPerLevel = finiteNumber(request.cellsPerLevel ?? (wanted < 50000 ? 3 : 6),
    '每级过渡带单元数', { min: 0, max: 64, integer: true });
  const minimumCells = Math.ceil(wanted * TARGET_LOWER_RATIO);
  const maximumCells = Math.floor(wanted * TARGET_UPPER_RATIO);
  const seed = {
    ...request,
    sizingMode: 'relative',
    wallRelativeSize,
    backgroundRelativeSize,
    cellsPerLevel
  };
  if (method.id === 'hybrid') seed.firstLayerRelativeSize = wallRelativeSize / 4;

  seed.budgetPlan = {
    targetCells: wanted,
    minimumCells,
    maximumCells,
    estimatedCells: estimateFromBackground(measuredArea.area, backgroundRelativeSize,
      referenceLength),
    isEstimate: true,
    area: measuredArea.area,
    domainSpan: bodySpan + 2 * farFieldSpans * referenceLength,
    referenceLength,
    areaSource: measuredArea.source,
    areaEstimated: measuredArea.estimated,
    bandReserveFraction: BAND_RESERVE_FRACTION,
    wallToBackgroundRatio,
    attempt: 1,
    maxAttempts: MAX_ATTEMPTS,
    limitedBySafeWallLevel,
    minimumWallRelativeSize
  };
  return seed;
}

function correctionBounds(previous) {
  const values = [previous.wallRelativeSize, previous.backgroundRelativeSize];
  if (previous.firstLayerRelativeSize !== undefined) values.push(previous.firstLayerRelativeSize);
  let minimum = 0.6;
  let maximum = 1.6;
  for (const value of values) {
    minimum = Math.max(minimum, 1e-12 / value);
    maximum = Math.min(maximum, (value === previous.backgroundRelativeSize ? 1e6 : 1) / value);
  }
  if (previous.allowUnsafeWallLevel !== true) {
    minimum = Math.max(minimum,
      (previous.budgetPlan?.minimumWallRelativeSize ?? 0) / previous.wallRelativeSize);
  }
  return { minimum, maximum };
}

function refineBudget(previous, actual, target = previous?.budgetPlan?.targetCells ?? previous?.targetCells) {
  if (!previous || typeof previous !== 'object') throw new Error('缺少上一次预算计划。');
  const assessment = assessBudget(actual, target);
  const oldPlan = previous.budgetPlan || {};
  const oldAttempt = finiteNumber(oldPlan.attempt ?? 1, '预算尝试次数', {
    min: 1,
    max: MAX_ATTEMPTS,
    integer: true
  });
  if (assessment.withinRange || oldAttempt >= MAX_ATTEMPTS) {
    return {
      ...previous,
      budgetPlan: {
        ...oldPlan,
        ...assessment,
        correctionFactor: 1,
        exhausted: !assessment.withinRange && oldAttempt >= MAX_ATTEMPTS
      }
    };
  }

  // Background sizes are quantized by the quadtree. While the desired count
  // exceeds that floor, change the width of the fine band instead of repeatedly
  // requesting different floating-point sizes that resolve to the same levels.
  const domainSpan = oldPlan.domainSpan;
  const reference = oldPlan.referenceLength;
  const level = Math.ceil(Math.log2(domainSpan / (previous.backgroundRelativeSize * reference)));
  const actualBackground = domainSpan / Math.pow(2, Math.max(0, level));
  const baseCells = oldPlan.area / (actualBackground * actualBackground);
  const extra = assessment.actualCells - baseCells;
  if (Number.isFinite(baseCells) && baseCells < assessment.targetCells && extra > 0 &&
      previous.cellsPerLevel > 0) {
    const nextBand = Math.max(0, Math.min(64, Math.round(previous.cellsPerLevel *
      (assessment.targetCells-baseCells)/extra)));
    if (nextBand !== previous.cellsPerLevel &&
        (assessment.actualCells > assessment.targetCells ? nextBand < previous.cellsPerLevel : nextBand > previous.cellsPerLevel)) {
      return { ...previous, cellsPerLevel: nextBand,
        budgetPlan: { ...oldPlan, ...assessment, attempt: oldAttempt+1,
          lastActualCells: assessment.actualCells, correctionFactor: 1,
          correctionKind: 'band-width', previousBand: previous.cellsPerLevel,
          estimatedBackgroundCells: Math.round(baseCells), exhausted: false } };
    }
  }
  let factor = Math.sqrt(assessment.actualCells / assessment.targetCells);
  const lastActual = Number(oldPlan.lastActualCells);
  const plateau = Number.isFinite(lastActual) &&
    Math.abs(assessment.actualCells - lastActual) / Math.max(assessment.actualCells, lastActual) <= 0.02;
  if (plateau) factor = assessment.actualCells < assessment.targetCells
    ? Math.min(factor, 0.8)
    : Math.max(factor, 1.25);
  const bounds = correctionBounds(previous);
  factor = Math.min(bounds.maximum, Math.max(bounds.minimum, factor));

  const refined = {
    ...previous,
    wallRelativeSize: previous.wallRelativeSize * factor,
    backgroundRelativeSize: previous.backgroundRelativeSize * factor
  };
  if (previous.firstLayerRelativeSize !== undefined) {
    refined.firstLayerRelativeSize = previous.firstLayerRelativeSize * factor;
  }
  refined.budgetPlan = {
    ...oldPlan,
    ...assessment,
    estimatedCells: Math.round((oldPlan.estimatedCells ?? assessment.actualCells) / (factor * factor)),
    attempt: oldAttempt + 1,
    lastActualCells: assessment.actualCells,
    correctionFactor: factor,
    plateau,
    constrained: factor === bounds.minimum || factor === bounds.maximum,
    exhausted: false
  };
  return refined;
}

module.exports = {
  BUDGET_PRESETS,
  MIN_TARGET_CELLS,
  MAX_TARGET_CELLS,
  planBudget,
  refineBudget,
  assessBudget
};
