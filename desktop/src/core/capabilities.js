'use strict';

// What the product actually offers, and where each method's evidence stops.  The
// renderer builds its form from this and the main process refuses anything the
// registry does not describe, so there is one place to correct a claim.

const SAFE_WALL_LEVEL = 11;

const METHODS = Object.freeze({
  cutcell: {
    id: 'cutcell',
    label: '纯 Cut-cell',
    executable: 'cartmesh2d_cli',
    status: 'stable',
    statusLabel: '稳定',
    summary: '背景笛卡尔网格 + 真实 cut-cell 多边形。外流默认语义。',
    // Measured in tools/verification/alignment_sensitivity.py: 20/20 pass at wall
    // level 7-11, 11/20 at 12.  Not a taste setting.
    safeWallLevel: SAFE_WALL_LEVEL,
    supports: { sizeField: true, boundaryLayer: false, openFoam: true },
    defaults: {
      farFieldSpans: 10,
      wallCellsPerSpan: 64,
      cellsPerLevel: 3,
      farLevel: 0,
      curvatureCellsPerRadius: 0,
      gapCells: 0,
      wake: null,
      smallAlpha: 0.15
    }
  },
  hybrid: {
    id: 'hybrid',
    label: '贴体边界层 + Cut-cell',
    executable: 'cartmesh2d_hybrid_cli',
    status: 'beta',
    statusLabel: 'Beta',
    summary: '壁面共形四边形层，外围接 Cartesian / Cut-cell 余域。',
    // docs/CURRENT_STATE_CN.md section 2: level 9 fails on faceWeight/volRatio/
    // nonOrtho because wall tangential resolution is pinned to the input polyline.
    safeWallLevel: 8,
    supports: { sizeField: false, boundaryLayer: true, openFoam: true },
    defaults: {
      maxLevel: 6,
      minimumLevel: 3,
      boundaryLevel: 6,
      nLayers: 4,
      firstThickness: 0.02,
      growthRatio: 1.2,
      domainPadding: 1.0,
      extrusionThickness: 0.01
    }
  }
});

// Named starting points for the size field.  Each one is a (far field, wall cells)
// pair that resolves on or under the safe wall level; the labels describe the mesh
// they produce rather than a vague "quality".
const PRESETS = Object.freeze([
  {
    id: 'preview',
    label: '预览',
    hint: '几秒出图，用来确认几何和分区',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 16, cellsPerLevel: 2 }
  },
  {
    id: 'balanced',
    label: '常规',
    hint: '远场 10 倍体长，壁面 1/64 体长，每级 3 格',
    sizeField: { farFieldSpans: 10, wallCellsPerSpan: 64, cellsPerLevel: 3 }
  },
  {
    id: 'wall',
    label: '壁面加密',
    hint: '壁面 1/128 体长；远场收到 6 倍以换取深度',
    sizeField: { farFieldSpans: 6, wallCellsPerSpan: 128, cellsPerLevel: 3 }
  },
  {
    id: 'farfield',
    label: '远场优先',
    hint: '远场 30 倍体长，适合看尾迹和阻塞比',
    sizeField: { farFieldSpans: 30, wallCellsPerSpan: 32, cellsPerLevel: 3 }
  }
]);

const GEOMETRY_FORMATS = Object.freeze([
  { extension: 'xy', label: '原生折线 (.xy)', native: true },
  { extension: 'dxf', label: 'AutoCAD DXF (.dxf)', native: false },
  { extension: 'svg', label: 'SVG 路径 (.svg)', native: false },
  { extension: 'csv', label: '逗号分隔坐标 (.csv)', native: false },
  { extension: 'txt', label: '空格分隔坐标 (.txt)', native: false }
]);

const methodById = id => METHODS[id];

module.exports = { METHODS, PRESETS, GEOMETRY_FORMATS, SAFE_WALL_LEVEL, methodById };
