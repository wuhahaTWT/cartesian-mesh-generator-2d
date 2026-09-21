'use strict';

// This adapter is only for drawing complete Cartesian leaves. It deliberately
// produces no CM2D file and is rejected by every flow entry point.
function parseBackgroundGrid(text) {
  const data = JSON.parse(text);
  const require = (ok, message) => { if (!ok) throw new Error(`背景网格：${message}`); };
  require(data.format === 'cartmesh2d-background-v1' && data.solver_ready === false &&
    data.retains_solid_interior === true, '格式或物理语义不匹配。');
  require(JSON.stringify(data.classification_names) === JSON.stringify(['outside','inside','intersected']), '分类定义不匹配。');
  require(['uniform', 'adaptive'].includes(data.mode), '未知网格模式。');
  const d = data.domain;
  require(Array.isArray(d) && d.length === 4 && d.every(Number.isFinite) && d[2] > d[0] && d[3] > d[1], '计算域无效。');
  require(Array.isArray(data.cells) && data.cells.length > 0 && data.cells.length <= 1048576, '单元数超出读取范围。');
  require(Array.isArray(data.boundary_loops) && data.boundary_loops.length > 0 &&
    data.boundary_loops.every(loop => Array.isArray(loop) && loop.length >= 3 &&
      loop.every(p => Array.isArray(p) && p.length === 2 && p.every(Number.isFinite))), '输入轮廓无效。');
  const vertices = [], cells = [], edges = [], counts = [0, 0, 0], nodes = new Set();
  let minLevel = Infinity, maxLevel = 0, areaFraction = 0;
  for (const cell of data.cells) {
    const { level, ix, iy, classification: kind, bounds: b } = cell;
    require(cell.id === cells.length && Number.isInteger(level) && level >= 0 && level <= 12 &&
      Number.isInteger(ix) && Number.isInteger(iy) && ix >= 0 && iy >= 0 && ix < 2**level && iy < 2**level, '单元编号或格坐标无效。');
    const key = `${level}/${ix}/${iy}`;
    require(!nodes.has(key), '重复单元。'); nodes.add(key);
    require(Number.isInteger(kind) && kind >= 0 && kind <= 2, '单元分类无效。');
    const expected = [d[0]+(d[2]-d[0])*ix/2**level, d[1]+(d[3]-d[1])*iy/2**level,
      d[0]+(d[2]-d[0])*(ix+1)/2**level, d[1]+(d[3]-d[1])*(iy+1)/2**level];
    const tolerance = 32*Number.EPSILON*Math.max(...d.map(Math.abs), d[2]-d[0], d[3]-d[1]);
    require(Array.isArray(b) && b.length === 4 && b.every((x,i) => Number.isFinite(x) && Math.abs(x-expected[i]) <= tolerance), '单元不是指定格坐标的完整矩形。');
    const first = vertices.length;
    vertices.push([b[0],b[1]], [b[2],b[1]], [b[2],b[3]], [b[0],b[3]]);
    cells.push({ id: cell.id, level, keyLevel: level, classification: kind,
      area: (b[2]-b[0])*(b[3]-b[1]), vertices: [first,first+1,first+2,first+3] });
    counts[kind]++; minLevel = Math.min(minLevel,level); maxLevel = Math.max(maxLevel,level);
    areaFraction += 4**(-level);
  }
  require(areaFraction === 1, '单元未覆盖完整计算域。');
  for (const cell of data.cells) for (let l=0;l<cell.level;l++)
    require(!nodes.has(`${l}/${cell.ix >> (cell.level-l)}/${cell.iy >> (cell.level-l)}`), '单元重叠。');
  require(Array.isArray(data.counts) && data.counts.length === 3 && counts.every((n,i) => n === data.counts[i]), '分类统计不一致。');
  require(data.mode !== 'uniform' || minLevel === maxLevel, '均匀模式含不同层级。');
  const outline = (loop, patch) => {
    const start=vertices.length; for(const point of loop)vertices.push(point);
    for(let i=0;i<loop.length;i++) edges.push({a:start+i,b:start+(i+1)%loop.length,owner:-1,neighbour:-1,patch});
  };
  // Display-only geometry edges; never numerical wall faces.
  for(const loop of data.boundary_loops) outline(loop,1);
  outline([[d[0],d[1]],[d[2],d[1]],[d[2],d[3]],[d[0],d[3]]],2);
  return { vertices, cells, edges, minLevel, maxLevel, background:true, solverReady:false,
    classificationCounts:counts, backgroundMode:data.mode,
    bounds:{minX:d[0],minY:d[1],maxX:d[2],maxY:d[3]} };
}

function requireFluidMesh(result) {
  if (!result || result.background || result.mesh?.background || !result.cm2dPath)
    throw new Error('请先生成 Cut-cell 或贴体混合流体网格；完整笛卡尔背景网格尚未接入浸入边界求解。');
}
module.exports={parseBackgroundGrid, requireFluidMesh};
