'use strict';

const finite = (value, name) => {
  const n = Number(value);
  if (!Number.isFinite(n)) throw new Error(`${name} 必须是有限数。`);
  return n;
};

function validateAdaptiveRequest(request = {}) {
  const refineFraction = finite(request.refineFraction ?? 0.15, '自适应加密单元比例');
  if (!(refineFraction > 0 && refineFraction <= 0.5)) {
    throw new Error('自适应加密单元比例须在 (0, 0.5] 范围内。');
  }
  const maxBoxes = finite(request.maxBoxes ?? 6, '最大加密区数量');
  if (!Number.isInteger(maxBoxes) || maxBoxes < 1 || maxBoxes > 32) {
    throw new Error('最大加密区数量须为 1 到 32 的整数。');
  }
  const padCells = finite(request.padCells ?? 1.0, '加密区外扩单元数');
  if (padCells < 0 || padCells > 8) {
    throw new Error('加密区外扩单元数须在 [0, 8] 范围内。');
  }
  const normalized = { refineFraction, maxBoxes, padCells };
  if (request.wallLevel !== undefined && request.wallLevel !== null && request.wallLevel !== '') {
    const wallLevel = finite(request.wallLevel, '参考壁面层级');
    if (!Number.isInteger(wallLevel) || wallLevel < 0 || wallLevel > 28) {
      throw new Error('参考壁面层级须为 0 到 28 的整数。');
    }
    normalized.wallLevel = wallLevel;
  }
  if (request.maxTargetLevel !== undefined && request.maxTargetLevel !== null && request.maxTargetLevel !== '') {
    const maxTargetLevel = finite(request.maxTargetLevel, '自适应最高目标层级');
    if (!Number.isInteger(maxTargetLevel) || maxTargetLevel < 0 || maxTargetLevel > 28) {
      throw new Error('自适应最高目标层级须为 0 到 28 的整数。');
    }
    normalized.maxTargetLevel = maxTargetLevel;
  }
  return normalized;
}

function buildAdaptiveInvocation(cm2dPath, cellsCsvPath, outputJsonPath, request = {}) {
  if (!cm2dPath || !cellsCsvPath || !outputJsonPath) {
    throw new Error('缺少自适应加密所需的网格、流场或输出路径。');
  }
  const normalized = validateAdaptiveRequest(request);
  const args = [
    '--mesh', cm2dPath,
    '--suggest-refinement-from', cellsCsvPath,
    '--suggest-refinement', outputJsonPath,
    '--adapt-fraction', String(normalized.refineFraction),
    '--adapt-max-boxes', String(normalized.maxBoxes),
    '--adapt-pad-cells', String(normalized.padCells)
  ];
  if (normalized.maxTargetLevel !== undefined) {
    args.push('--adapt-max-level', String(normalized.maxTargetLevel));
  }
  if (normalized.wallLevel !== undefined) {
    args.push('--adapt-wall-level', String(normalized.wallLevel));
  }
  if (request.scalarCsvPath) {
    args.push('--adapt-scalar-csv', request.scalarCsvPath);
  }
  if (request.vtkPath) {
    args.push('--adapt-vtk', request.vtkPath);
  }
  if (request.fieldsJsonPath) {
    args.push('--adapt-fields-json', request.fieldsJsonPath);
  }
  return { executable: 'cartmesh2d_flow_cli', args, request: normalized };
}

function validateAdaptiveFields(payload, expectedCells) {
  if (!payload || payload.format !== 'cartmesh2d-adaptive-fields-v1' || !Array.isArray(payload.cells)) {
    throw new Error('自适应梯度场文件格式无效。');
  }
  if (expectedCells !== undefined && payload.cells.length !== expectedCells) {
    throw new Error('自适应梯度场单元数与当前最终网格不匹配。');
  }
  payload.cells.forEach((cell, index) => {
    if (!cell || cell.id !== index || !Number.isFinite(cell.indicator)) {
      throw new Error(`自适应梯度场第 ${index} 号单元无效。`);
    }
  });
  return payload;
}

function roundBoxPair(lo, hi) {
  const rLo = Number(lo.toFixed(4));
  const rHi = Number(hi.toFixed(4));
  return rHi > rLo ? [rLo, rHi] : [lo, hi];
}

function parseAdaptiveSuggestion(payload, frame = null) {
  if (!payload || typeof payload !== 'object') {
    throw new Error('自适应加密结果格式无效。');
  }
  if (payload.format !== 'cartmesh2d-adaptive-refinement-v1') {
    throw new Error(`未知自适应加密格式：${payload.format || '缺失'}。`);
  }
  if (payload.valid !== true) {
    const reason = Array.isArray(payload.issues) && payload.issues.length
      ? payload.issues.join('；') : '原生自适应模块返回无效状态';
    throw new Error(`流场自适应加密失败：${reason}`);
  }
  if (!Array.isArray(payload.boxes)) {
    throw new Error('自适应加密结果缺少 boxes 列表。');
  }
  const useFrame = frame && Number.isFinite(frame.cx) && Number.isFinite(frame.cy) &&
    Number.isFinite(frame.span) && frame.span > 0;
  const regions = payload.boxes.map((box, index) => {
    const xminW = finite(box.xmin, `加密区 ${index + 1} xmin`);
    const yminW = finite(box.ymin, `加密区 ${index + 1} ymin`);
    const xmaxW = finite(box.xmax, `加密区 ${index + 1} xmax`);
    const ymaxW = finite(box.ymax, `加密区 ${index + 1} ymax`);
    if (!(xmaxW > xminW) || !(ymaxW > yminW)) {
      throw new Error(`加密区 ${index + 1} 坐标范围无效。`);
    }
    const levelsBelowWall = finite(box.levels_below_wall ?? 0, `加密区 ${index + 1} 低于壁面级数`);
    if (!Number.isInteger(levelsBelowWall) || levelsBelowWall < 0 || levelsBelowWall > 20) {
      throw new Error(`加密区 ${index + 1} 低于壁面级数超出范围。`);
    }
    const rawXmin = useFrame ? (xminW - frame.cx) / frame.span : finite(box.body_xmin, `加密区 ${index + 1} body_xmin`);
    const rawXmax = useFrame ? (xmaxW - frame.cx) / frame.span : finite(box.body_xmax, `加密区 ${index + 1} body_xmax`);
    const rawYmin = useFrame ? (yminW - frame.cy) / frame.span : finite(box.body_ymin, `加密区 ${index + 1} body_ymin`);
    const rawYmax = useFrame ? (ymaxW - frame.cy) / frame.span : finite(box.body_ymax, `加密区 ${index + 1} body_ymax`);
    const [xmin, xmax] = roundBoxPair(rawXmin, rawXmax);
    const [ymin, ymax] = roundBoxPair(rawYmin, rawYmax);
    return {
      xmin,
      xmax,
      ymin,
      ymax,
      levelsBelowWall,
      active: true,
      targetLevel: Number(box.target_level ?? 0),
      cellCount: Number(box.cell_count ?? 0),
      meanIndicator: Number(box.mean_indicator ?? 0),
      maxIndicator: Number(box.max_indicator ?? 0)
    };
  });
  return {
    cellCount: Number(payload.cell_count ?? 0),
    selectedCellCount: Number(payload.selected_cell_count ?? 0),
    selectedAreaFraction: Number(payload.selected_area_fraction ?? 0),
    inferredWallLevel: Number(payload.inferred_wall_level ?? 0),
    maxTargetLevel: Number(payload.max_target_level ?? 0),
    indicatorMax: Number(payload.indicator_max ?? 0),
    referenceScalarSpan: Number(payload.reference_scalar_span ?? 0),
    cellTargetCount: Array.isArray(payload.cell_targets) ? payload.cell_targets.length : 0,
    regions
  };
}

module.exports = {
  validateAdaptiveRequest,
  buildAdaptiveInvocation,
  validateAdaptiveFields,
  parseAdaptiveSuggestion
};
