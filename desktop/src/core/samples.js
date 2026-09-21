'use strict';

// Built-in geometry so the app is usable without hunting for a DXF.  Every entry
// is a file that ships in runtime/samples, a fluid-region semantics, and a sizing
// request that has been run and produces a mesh.

const SAMPLES = Object.freeze([
  {
    id: 'circle',
    label: '圆柱',
    file: 'circle.xy',
    fluidRegion: 'exterior',
    note: '32 段折线圆。最经典的外流验收几何。圆看不出来流方向，所以不带尾迹加密——'
      + '尾迹只在形状本身指明来流方向时才有意义。',
    sizeField: { farFieldSpans: 10, wallCellsPerSpan: 64, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'naca2412',
    label: 'NACA 2412 翼型',
    file: 'naca2412_dense.xy',
    fluidRegion: 'exterior',
    note: '密采样翼型。后缘曲率半径很小，开曲率细化会要求很深的层级。',
    sizeField: { farFieldSpans: 6, wallCellsPerSpan: 48, cellsPerLevel: 3 },
    smallAlpha: 0.45,
    wake: { angleOfAttackDeg: 4, downstreamSpans: 8, halfWidthSpans: 0.7, levelsBelowWall: 4 }
  },
  {
    id: 'thick_airfoil',
    label: '厚弯度翼型',
    file: 'thick_cambered_airfoil.xy',
    fluidRegion: 'exterior',
    note: '钝后缘厚翼型，可观察尾缘与尾迹细化。过粗数量档可能无法通过质量门；自动模式最多细化重试两次。',
    sizeField: { farFieldSpans: 10, wallCellsPerSpan: 64, cellsPerLevel: 3 },
    smallAlpha: 0.35,
    wake: { angleOfAttackDeg: 0, downstreamSpans: 8, halfWidthSpans: 0.7, levelsBelowWall: 4 }
  },
  {
    id: 'two_obstacles',
    label: '双柱串列',
    file: 'two_obstacles.xy',
    fluidRegion: 'exterior',
    note: '两个互不相连的固壁环。多体外流，检验 even-odd 语义。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 64, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'narrow_gap',
    label: '窄缝',
    file: 'narrow_gap.xy',
    fluidRegion: 'exterior',
    note: '两段壁面隔 0.08 相对。用来看间隙细化：不开的话缝里只有一两格。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.15,
    gapCells: 4
  },
  {
    id: 'sharp_trailing_edge',
    label: '尖后缘',
    file: 'sharp_trailing_edge.xy',
    fluidRegion: 'exterior',
    note: '楔形尖尾。曲率细化会在尖点要求很深层级，是观察安全上限的例子。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'gear_star',
    label: '星形齿轮',
    file: 'gear_star.xy',
    fluidRegion: 'exterior',
    note: '交替凹凸的强曲率轮廓，检验 cut-cell 在凹角上的构造。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'superellipse',
    label: '超椭圆',
    file: 'superellipse_24.xy',
    fluidRegion: 'exterior',
    note: '对称超椭圆。对称切点容易擦过格点，是 W1 焦合预算的回归几何。',
    sizeField: { farFieldSpans: 10, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.25
  },
  {
    id: 'serpentine',
    label: '蛇形体',
    file: 'serpentine_body.xy',
    fluidRegion: 'exterior',
    note: '细长弯曲物体，用于检验曲线壁面切分与贴体层连接。支持自动常规和更密档。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 24, cellsPerLevel: 3 },
    smallAlpha: 0.45
  },
  {
    id: 'annulus',
    label: '方形环隙（内流）',
    file: 'annulus.xy',
    fluidRegion: 'interior',
    note: '方形外边界加一个方孔，流体在两者之间。闭域没有左右开口，不能套用通道／喷管预设。圆形旋转验证请选“同心圆旋转环隙”。',
    sizeField: { farFieldSpans: 1, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'rotating_annulus',
    label: '同心圆旋转环隙',
    file: 'rotating_annulus.xy',
    fluidRegion: 'interior',
    note: '半径0.5/1m的64段同心圆。选择命名边界，再生成“同心圆环”模板：内壁逆时针运动，外壁静止；几何网格保持固定。',
    sizeField: { farFieldSpans: .25, wallCellsPerSpan: 16, cellsPerLevel: 3 },
    interiorSizeField: { farFieldSpans: .25, wallCellsPerSpan: 16, cellsPerLevel: 3, farLevel: 4 },
    interiorSmallAlpha: .05,
    smallAlpha: .05
  },
  {
    id: 'rectangle', label: '矩形通道（内流）', file: 'rectangle.xy', fluidRegion: 'interior',
    note: '2×1m矩形。支持速度入口或两端静压；选择“半通道”边界模板时，顶部为对称面、底部无滑移。',
    sizeField: { farFieldSpans: .25, wallCellsPerSpan: 16, cellsPerLevel: 3 },
    interiorSizeField: { farFieldSpans: .25, wallCellsPerSpan: 16, cellsPerLevel: 3, farLevel: 4 },
    interiorSmallAlpha: .1, smallAlpha: .1
  },
  {
    id: 'nozzle',
    label: '喷管型线',
    file: 'nozzle_profile.xy',
    fluidRegion: 'interior',
    note: '收缩扩张喷管，默认生成内部流体网格。不可压计算选择“曲壁通道／喷管”：左均匀速度入口、右零压出口、曲壁无滑移。',
    sizeField: { farFieldSpans: 6, wallCellsPerSpan: 16, cellsPerLevel: 3 },
    interiorSizeField: { farFieldSpans: 0.5, wallCellsPerSpan: 128, cellsPerLevel: 3, farLevel: 7 },
    interiorSmallAlpha: 0.45,
    smallAlpha: 0.35
  }
]);

const sampleById = id => SAMPLES.find(sample => sample.id === id);

module.exports = { SAMPLES, sampleById };
