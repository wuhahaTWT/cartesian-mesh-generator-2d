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
    note: '32 段折线圆。最经典的外流验收几何，也是加密阶梯的基准。',
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
    note: '钝后缘厚翼型，比尖后缘容易收敛，适合先试尾迹加密。',
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
    note: '细长弯曲物体，自身前后段互相靠近。间隙细化默认关闭：开了会在弯段外侧把边界歪斜推到 4.3–5.7，越过硬限 4。',
    sizeField: { farFieldSpans: 8, wallCellsPerSpan: 24, cellsPerLevel: 3 },
    smallAlpha: 0.45
  },
  {
    id: 'annulus',
    label: '圆环（内流）',
    file: 'annulus.xy',
    fluidRegion: 'interior',
    note: '外环加一个孔。这是内流语义，流体在两环之间，不是默认外流。',
    sizeField: { farFieldSpans: 1, wallCellsPerSpan: 32, cellsPerLevel: 3 },
    smallAlpha: 0.15
  },
  {
    id: 'nozzle',
    label: '喷管型线',
    file: 'nozzle_profile.xy',
    fluidRegion: 'exterior',
    note: '收缩扩张型线，长宽比大，能看出方形计算域怎么处理非方几何。',
    sizeField: { farFieldSpans: 6, wallCellsPerSpan: 16, cellsPerLevel: 3 },
    smallAlpha: 0.35
  }
]);

const sampleById = id => SAMPLES.find(sample => sample.id === id);

module.exports = { SAMPLES, sampleById };
