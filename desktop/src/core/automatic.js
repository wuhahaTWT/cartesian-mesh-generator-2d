'use strict';
const { wallLevelFor, wallCellsPerSpanCeiling } = require('./sizing');
const { SAFE_WALL_LEVEL } = require('./capabilities');

// A bounded search over generation parameters, never over acceptance thresholds.
// Historical hybrid seeds: tools/verification/generate_q0_baselines.py.
function candidates(request, sample, frame) {
  if (!request.automatic) return [request];
  const dense = request.density === 'dense';
  if (request.method === 'hybrid') {
    if (request.fluidRegion === 'interior') {
      const scale = frame.bodySpan;
      const level = dense ? 7 : 6;
      const seed = { ...request, maxLevel: level, boundaryLevel: level, minimumLevel: level - 1,
        nLayers: 3, remainderSmallAlpha: sample?.id === 'nozzle' ? 0.65 : 0.45, firstThickness: scale * (sample?.id === 'nozzle' ? 0.0025 : 0.01), growthRatio: 1.15,
        domainPadding: scale / 12, extrusionThickness: 0.01 * scale };
      return [seed, { ...seed, firstThickness: seed.firstThickness * 1.5 },
        { ...seed, firstThickness: seed.firstThickness * 0.65 },
        { ...seed, maxLevel: 6, boundaryLevel: 6, minimumLevel: 5, firstThickness: seed.firstThickness }]
        .filter((v,i,a) => a.findIndex(x => JSON.stringify(x) === JSON.stringify(v)) === i);
    }
    const sharp = ['sharp_trailing_edge', 'narrow_gap'].includes(sample?.id);
    const superellipse = sample?.id === 'superellipse';
    const scale = frame.bodySpan;
    const seed = { ...request, maxLevel: sharp ? 8 : 6, minimumLevel: 3,
      boundaryLevel: sharp ? 8 : 6, nLayers: superellipse ? 3 : 4,
      firstThickness: sharp ? 0.012 : superellipse ? 0.015 : 0.01 * scale,
      growthRatio: sharp || superellipse ? 1.15 : 1.2,
      domainPadding: sharp || superellipse ? 1 : 0.5 * scale, extrusionThickness: 0.01 * scale };
    const preferred = dense && !sharp ? { ...seed, maxLevel: 7, boundaryLevel: 7 } : seed;
    return [preferred, { ...preferred, firstThickness: seed.firstThickness * 0.65 },
      { ...seed, maxLevel: 8, boundaryLevel: 8, firstThickness: seed.firstThickness * 0.5 }, seed]
      .filter((v,i,a) => a.findIndex(x => JSON.stringify(x) === JSON.stringify(v)) === i);
  }
  const interiorProfile = request.fluidRegion === 'interior' && sample?.interiorSizeField;
  const field = interiorProfile || sample?.sizeField || { farFieldSpans: 6, wallCellsPerSpan: 32, cellsPerLevel: 3 };
  const seed = { ...request, ...field, smallAlpha: (interiorProfile ? sample.interiorSmallAlpha : sample?.smallAlpha) ?? 0.25,
    farFieldSpans: request.fluidRegion === 'interior' && !interiorProfile ? 1 : field.farFieldSpans,
    wallCellsPerSpan: field.wallCellsPerSpan * (dense ? 2 : 1),
    farLevel: interiorProfile ? field.farLevel + (dense ? 1 : 0) : 0, curvatureCellsPerRadius: 0, gapCells: sample?.gapCells || 0,
    wake: request.fluidRegion === 'interior' ? null : sample?.wake || null,
    refineBoxes: [], allowUnsafeWallLevel: false };
  const requestedWall = seed.wallCellsPerSpan;
  seed.wallCellsPerSpan = Math.min(requestedWall, Math.floor(wallCellsPerSpanCeiling(seed.farFieldSpans, SAFE_WALL_LEVEL)));
  // At the depth ceiling, widen the fine band to add actual cells. The domain
  // and physical boundary stay fixed; repeating the same forbidden depth adds nothing.
  if (dense && wallLevelFor(seed.farFieldSpans, seed.wallCellsPerSpan) === wallLevelFor(seed.farFieldSpans, field.wallCellsPerSpan))
    seed.cellsPerLevel = Math.max(5, field.cellsPerLevel);
  return [seed, { ...seed, smallAlpha: 0.45 },
    { ...seed, cellsPerLevel: 4, smallAlpha: 0.35 },
    { ...seed, wallCellsPerSpan: field.wallCellsPerSpan,
      farLevel: interiorProfile ? field.farLevel : 0,
      smallAlpha: interiorProfile ? sample.interiorSmallAlpha : 0.15 }]
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
