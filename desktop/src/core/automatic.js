'use strict';

// A bounded search over generation parameters, never over acceptance thresholds.
// Historical hybrid seeds: tools/verification/generate_q0_baselines.py.
function candidates(request, sample, frame) {
  if (!request.automatic) return [request];
  const dense = request.density === 'dense';
  if (request.method === 'hybrid') {
    if (request.fluidRegion === 'interior') throw new Error('内部网格当前请使用纯 Cut-cell；贴体边界层只支持外流。');
    const sharp = ['sharp_trailing_edge', 'narrow_gap'].includes(sample?.id);
    const superellipse = sample?.id === 'superellipse';
    const scale = frame.bodySpan;
    const seed = { ...request, maxLevel: sharp ? 8 : 6, minimumLevel: 3,
      boundaryLevel: sharp ? 8 : 6, nLayers: superellipse ? 3 : 4,
      firstThickness: sharp ? 0.012 : superellipse ? 0.015 : 0.01 * scale,
      growthRatio: sharp || superellipse ? 1.15 : 1.2,
      domainPadding: sharp || superellipse ? 1 : 0.5 * scale, extrusionThickness: 0.01 * scale };
    const preferred = dense && !sharp ? { ...seed, maxLevel: 7, boundaryLevel: 7 } : seed;
    return [preferred, { ...seed, firstThickness: seed.firstThickness * 0.65 },
      { ...seed, maxLevel: 8, boundaryLevel: 8, firstThickness: seed.firstThickness * 0.5 }, seed]
      .filter((v,i,a) => a.findIndex(x => JSON.stringify(x) === JSON.stringify(v)) === i);
  }
  const field = sample?.sizeField || { farFieldSpans: 6, wallCellsPerSpan: 32, cellsPerLevel: 3 };
  const seed = { ...request, ...field, smallAlpha: sample?.smallAlpha ?? 0.25,
    farFieldSpans: request.fluidRegion === 'interior' ? 1 : field.farFieldSpans,
    wallCellsPerSpan: field.wallCellsPerSpan * (dense ? 2 : 1),
    farLevel: 0, curvatureCellsPerRadius: 0, gapCells: sample?.gapCells || 0,
    wake: request.fluidRegion === 'interior' ? null : sample?.wake || null,
    refineBoxes: [], allowUnsafeWallLevel: false };
  return [seed, { ...seed, smallAlpha: 0.45 },
    { ...seed, cellsPerLevel: 4, smallAlpha: 0.35 },
    { ...seed, wallCellsPerSpan: field.wallCellsPerSpan, smallAlpha: 0.15 }]
    .filter((v,i,a) => a.findIndex(x => JSON.stringify(x) === JSON.stringify(v)) === i);
}

function estimateSeconds(request, historySeconds) {
  if (historySeconds > 0) return { seconds: historySeconds, source: '同类历史耗时' };
  const level = request.method === 'hybrid' ? request.maxLevel :
    Math.ceil(Math.log2((1 + 2 * request.farFieldSpans) * request.wallCellsPerSpan));
  return { seconds: request.method === 'hybrid' ? (level >= 8 ? 60 : 10) : (level >= 11 ? 30 : 5),
    source: '首次粗估，可能提前或超时' };
}
module.exports = { candidates, estimateSeconds };
