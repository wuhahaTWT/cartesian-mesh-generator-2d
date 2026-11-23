'use strict';

// High-Re boundary-layer y+ sizing estimator for cell-centered 2D FVM meshes
// (ANSYS Fluent, OpenFOAM, CartMesh2D). Uses standard flat-plate correlations
// (Schlichting turbulent / Blasius laminar) to compute first-layer total cell
// height dy1 = 2 * y_P so the wall-adjacent cell center sits at the target y+.

const FLUID_PRESETS = Object.freeze({
  air20: {
    id: 'air20',
    label: '标准空气（20 °C）',
    rho: 1.204,
    mu: 1.825e-5,
    nu: 1.825e-5 / 1.204
  },
  water20: {
    id: 'water20',
    label: '液态水（20 °C）',
    rho: 998.2,
    mu: 1.002e-3,
    nu: 1.002e-3 / 998.2
  }
});

const TURBULENT_RE_THRESHOLD = 5e5;
const MIN_RECOMMENDED_LAYERS = 4;
const MAX_SUPPORTED_LAYERS = 64;
const WAKE_REGION_SOURCE = 'yplus-wake';

function roundCoord(value) {
  return Math.round(value * 1000) / 1000;
}

function buildRecommendedWakeBox(input = {}) {
  const fluidRegion = input.fluidRegion === 'interior' ? 'interior' : 'exterior';
  if (fluidRegion === 'interior') return null;

  const referenceLength = positiveFinite(input.referenceLength, '参考长度 Lref');
  const totalThickness = Number.isFinite(Number(input.totalThickness)) ? Math.max(0, Number(input.totalThickness)) : 0;
  const delta99 = Number.isFinite(Number(input.delta99)) ? Math.max(0, Number(input.delta99)) : 0;
  const wakeLengthSpans = positiveFinite(input.wakeLengthSpans ?? 2.0, '尾迹长度');
  if (wakeLengthSpans < 0.2 || wakeLengthSpans > 50) {
    throw new Error('尾迹长度需在 0.2 至 50 Lref 之间。');
  }
  const wakeAngleDeg = Number(input.wakeAngleDeg ?? 0);
  if (!Number.isFinite(wakeAngleDeg) || wakeAngleDeg < -45 || wakeAngleDeg > 45) {
    throw new Error('尾迹偏角需在 -45° 至 45° 之间。');
  }
  const wakeLevelsBelowWall = Math.round(Number(input.wakeLevelsBelowWall ?? 2));
  if (!Number.isInteger(wakeLevelsBelowWall) || wakeLevelsBelowWall < 0 || wakeLevelsBelowWall > 10) {
    throw new Error('尾迹低于壁面级数需为 0 至 10 的整数。');
  }

  const frame = input.frame && typeof input.frame === 'object' ? input.frame : null;
  const bodyHalfX = frame && Number(frame.width) > 0
    ? (0.5 * Number(frame.width)) / referenceLength
    : 0.5;
  const bodyHalfY = frame && Number(frame.height) > 0
    ? (0.5 * Number(frame.height)) / referenceLength
    : 0.15;

  const blThicknessRel = Math.max(totalThickness, delta99) / referenceLength;
  const blMargin = Math.max(0.12, 2.5 * blThicknessRel);
  const upstreamMargin = Math.min(0.25, blMargin);
  const halfY = Math.max(0.35, bodyHalfY + blMargin + 0.12);
  const dyWake = wakeLengthSpans * Math.tan((wakeAngleDeg * Math.PI) / 180);

  let xmin = -(bodyHalfX + upstreamMargin);
  let xmax = bodyHalfX + wakeLengthSpans;
  let ymin = Math.min(-halfY, -halfY + dyWake);
  let ymax = Math.max(halfY, halfY + dyWake);

  // Clamp strictly inside the active far-field domain when a finite farFieldSpans
  // is known, so the wake box never protrudes outside the background tree.
  const farFieldSpans = Number(input.farFieldSpans);
  if (Number.isFinite(farFieldSpans) && farFieldSpans > 0) {
    const maxAbsX = bodyHalfX + 0.85 * farFieldSpans;
    const maxAbsY = bodyHalfY + 0.85 * farFieldSpans;
    xmin = Math.max(-maxAbsX, xmin);
    xmax = Math.min(maxAbsX, xmax);
    ymin = Math.max(-maxAbsY, ymin);
    ymax = Math.min(maxAbsY, ymax);
    if (!(xmax > xmin) || !(ymax > ymin)) return null;
  }

  return {
    xmin: roundCoord(xmin),
    xmax: roundCoord(xmax),
    ymin: roundCoord(ymin),
    ymax: roundCoord(ymax),
    levelsBelowWall: wakeLevelsBelowWall,
    source: WAKE_REGION_SOURCE
  };
}

function upsertYPlusWakeRegion(regions, wakeBox) {
  const retained = (Array.isArray(regions) ? regions : [])
    .filter(box => box && box.source !== WAKE_REGION_SOURCE)
    .map(box => ({ ...box, active: false }));
  if (!wakeBox) return retained;
  retained.push({
    xmin: wakeBox.xmin,
    xmax: wakeBox.xmax,
    ymin: wakeBox.ymin,
    ymax: wakeBox.ymax,
    levelsBelowWall: wakeBox.levelsBelowWall,
    source: WAKE_REGION_SOURCE,
    active: true
  });
  return retained;
}

function positiveFinite(value, label) {
  const n = Number(value);
  if (!Number.isFinite(n) || !(n > 0)) {
    throw new Error(`${label}必须为正有限数值。`);
  }
  return n;
}

function estimateBoundaryLayerYPlus(input = {}) {
  const referenceLength = positiveFinite(input.referenceLength, '参考长度 Lref');
  const velocity = positiveFinite(input.velocity, '来流速度 U∞');
  const nu = positiveFinite(input.nu, '运动黏度 ν');
  const targetYPlus = positiveFinite(input.targetYPlus ?? 1, '目标 y+');
  const growthRatio = positiveFinite(input.growthRatio ?? 1.2, '增长率');
  if (growthRatio < 1 || growthRatio > 5) {
    throw new Error('增长率需在 1 至 5 之间。');
  }
  const coverageRatio = positiveFinite(input.coverageRatio ?? 1.0, '边界层覆盖倍数');
  if (coverageRatio < 0.1 || coverageRatio > 5) {
    throw new Error('边界层覆盖倍数需在 0.1 至 5 之间。');
  }
  const requestedRegime = input.regime || 'auto';
  if (!['auto', 'turbulent', 'laminar'].includes(requestedRegime)) {
    throw new Error(`未知流动状态：${requestedRegime}`);
  }

  const reynolds = (velocity * referenceLength) / nu;
  const effectiveRegime = requestedRegime === 'auto'
    ? (reynolds >= TURBULENT_RE_THRESHOLD ? 'turbulent' : 'laminar')
    : requestedRegime;

  const warnings = [];
  let cf;
  let delta99;
  let correlation;

  if (effectiveRegime === 'turbulent') {
    const denom = 2 * Math.log10(Math.max(reynolds, 100)) - 0.65;
    if (!(denom > 0.5)) {
      throw new Error('雷诺数过低，无法使用 Schlichting 湍流经验公式，请切换为层流或自动模式。');
    }
    cf = Math.pow(denom, -2.3);
    delta99 = 0.37 * referenceLength * Math.pow(reynolds, -0.2);
    correlation = 'Schlichting 湍流平板（Cf = [2 log10(Re) - 0.65]^-2.3）';
    if (reynolds < 1e4) {
      warnings.push(`Re = ${reynolds.toExponential(2)} 低于湍流经验式适用下限（建议切换自动或层流模式）。`);
    }
    if (targetYPlus > 5 && targetYPlus < 30) {
      warnings.push('目标 y+ 位于缓冲层（5 < y+ < 30）：壁面解析模型建议 y+ ≤ 1，标准壁面函数建议 y+ ≥ 30。');
    }
  } else {
    cf = 0.664 / Math.sqrt(reynolds);
    delta99 = 4.91 * referenceLength / Math.sqrt(reynolds);
    correlation = 'Blasius 层流平板（Cf = 0.664 Re^-0.5）';
    if (reynolds > 1e6 && requestedRegime === 'laminar') {
      warnings.push(`Re = ${reynolds.toExponential(2)} 较高，实际工程流动通常已转捩为湍流。`);
    }
  }

  const uTau = velocity * Math.sqrt(cf / 2);
  const firstCellCenterDistance = (targetYPlus * nu) / uTau;
  // Cell-centered FVM places the wall-adjacent cell centroid at half the first
  // layer thickness, so dy1 = 2 * y_P achieves target y+ at the cell center.
  const firstLayerHeight = 2 * firstCellCenterDistance;
  const firstLayerRelativeSize = firstLayerHeight / referenceLength;
  const targetTotalThickness = coverageRatio * delta99;

  let rawLayerCount;
  if (Math.abs(growthRatio - 1) < 1e-9) {
    rawLayerCount = Math.ceil(targetTotalThickness / firstLayerHeight);
  } else {
    const arg = 1 + (targetTotalThickness / firstLayerHeight) * (growthRatio - 1);
    rawLayerCount = Math.ceil(Math.log(Math.max(arg, 1)) / Math.log(growthRatio));
  }
  const layersCapped = rawLayerCount > MAX_SUPPORTED_LAYERS;
  const nLayers = input.explicitLayers != null
    ? Math.min(MAX_SUPPORTED_LAYERS, Math.max(1, Math.round(Number(input.explicitLayers))))
    : Math.min(MAX_SUPPORTED_LAYERS, Math.max(MIN_RECOMMENDED_LAYERS, rawLayerCount));
  if (layersCapped && input.explicitLayers == null) {
    warnings.push(`覆盖目标 δ99 需约 ${rawLayerCount} 层，已按支持上限截断为 ${MAX_SUPPORTED_LAYERS} 层；可适当提高增长率（如 1.25）或采用壁面函数 y+ ≈ 30。`);
  }

  const sumFactor = Math.abs(growthRatio - 1) < 1e-9
    ? nLayers
    : (Math.pow(growthRatio, nLayers) - 1) / (growthRatio - 1);
  const totalThickness = firstLayerHeight * sumFactor;
  const lastLayerHeight = firstLayerHeight * Math.pow(growthRatio, nLayers - 1);
  const coverageAchieved = totalThickness / delta99;

  // Match outer Cartesian wall cell size smoothly (~2x last BL layer height,
  // bounded to [2*dy1, 0.1*Lref]) so the BL-to-remainder volume jump stays mild.
  const recommendedWallSize = Math.min(
    0.1 * referenceLength,
    Math.max(2 * firstLayerHeight, 2.0 * lastLayerHeight)
  );
  const recommendedWallRelativeSize = recommendedWallSize / referenceLength;
  const recommendedBackgroundRelativeSize = Math.max(
    recommendedWallRelativeSize * 8,
    Math.min(1.0, recommendedWallRelativeSize * 32)
  );
  const estimatedFirstLayerAspectRatio = recommendedWallSize / firstLayerHeight;
  if (firstLayerRelativeSize < 1e-7) {
    warnings.push('首层高度相对 Lref 极薄（< 1e-7），请确保几何折线平滑并关注浮点精度与求解器收敛。');
  }
  if (totalThickness > 0.25 * referenceLength) {
    warnings.push('估算边界层总厚度超过 0.25 Lref（低 Re 厚边界层），复杂窄缝/凹角可能触发局部降层。');
  }

  const fluidRegion = input.fluidRegion === 'interior' ? 'interior' : 'exterior';
  const wakeLengthSpans = input.wakeLengthSpans ?? 2.0;
  const recommendedFarFieldSpans = fluidRegion === 'interior'
    ? null
    : (reynolds >= 1e5
        ? Math.max(15, Math.ceil(Number(wakeLengthSpans) * 5))
        : Math.max(8, Math.ceil(Number(wakeLengthSpans) * 3)));
  const recommendedWakeBox = buildRecommendedWakeBox({
    fluidRegion,
    referenceLength,
    totalThickness,
    delta99,
    frame: input.frame,
    farFieldSpans: input.farFieldSpans,
    wakeLengthSpans,
    wakeAngleDeg: input.wakeAngleDeg ?? 0,
    wakeLevelsBelowWall: input.wakeLevelsBelowWall ?? 2
  });

  return {
    referenceLength,
    velocity,
    nu,
    targetYPlus,
    growthRatio,
    coverageRatio,
    fluidRegion,
    reynolds,
    effectiveRegime,
    correlation,
    cf,
    uTau,
    firstCellCenterDistance,
    firstLayerHeight,
    firstLayerRelativeSize,
    delta99,
    nLayers,
    rawLayerCount,
    layersCapped,
    totalThickness,
    totalRelativeThickness: totalThickness / referenceLength,
    lastLayerHeight,
    coverageAchieved,
    recommendedWallSize,
    recommendedWallRelativeSize,
    recommendedBackgroundRelativeSize,
    recommendedFarFieldSpans,
    recommendedWakeBox,
    estimatedFirstLayerAspectRatio,
    warnings
  };
}

module.exports = {
  FLUID_PRESETS,
  TURBULENT_RE_THRESHOLD,
  MAX_SUPPORTED_LAYERS,
  WAKE_REGION_SOURCE,
  buildRecommendedWakeBox,
  upsertYPlusWakeRegion,
  estimateBoundaryLayerYPlus
};
