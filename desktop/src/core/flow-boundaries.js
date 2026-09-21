'use strict';

const KINDS = Object.freeze(['velocity-inlet','pressure-outlet','pressure-opening','wall','moving-wall','smooth-moving-wall']);
const fail = message => { throw new Error(`流动边界：${message}`); };
const finite = value => { if (typeof value !== 'number' || !Number.isFinite(value)) fail('数值必须有限。'); return value; };
const numberToken = value => {
  if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/.test(value)) fail('文件数值格式无效。');
  return finite(Number(value));
};
const count = (value, zero = false) => {
  if (!Number.isSafeInteger(value) || value < (zero ? 0 : 1)) fail('数量或面编号无效。');
  return value;
};
function quotedTokens(line) {
  const matches = [...line.matchAll(/"(?:[^"\\]|\\.)*"|[^\s"]+/g)];
  let end = 0;
  const values = matches.map(match => {
    if (line.slice(end, match.index).trim()) fail('文件字符串无效。');
    if (end && match.index === end) fail('文件字段缺少分隔。');
    end = match.index + match[0].length;
    return match[0].startsWith('"') ? match[0].slice(1,-1).replace(/\\(.)/g,'$1') : match[0];
  });
  if (line.slice(end).trim()) fail('未闭合字符串。');
  return values;
}
function condition(entry) {
  if (!entry || typeof entry !== 'object') fail('记录无效。');
  const {face,type,name,u,v,p} = entry;
  count(face, true);
  if (!KINDS.includes(type)) fail('不支持的边界类型。');
  if (typeof name !== 'string' || !name.length || Buffer.byteLength(name,'utf8') > 128 || /[\x00-\x1f\x7f,"]/.test(name)) fail('名称无效。');
  [u,v,p].forEach(finite);
  if (['pressure-outlet','pressure-opening'].includes(type) ? (u !== 0 || v !== 0) : p !== 0) fail('压力与速度条件冲突。');
  if (type === 'wall' && (u !== 0 || v !== 0)) fail('静止壁面速度必须为零。');
  return {face,type,name,u,v,p};
}
function conditions(entries) {
  if (!Array.isArray(entries) || !entries.length) fail('请先生成或导入边界。');
  const seen = new Set(), names = new Map();
  const result = entries.map(entry => {
    const value = condition(entry);
    if (seen.has(value.face)) fail('面编号重复。');
    seen.add(value.face);
    if (names.has(value.name) && names.get(value.name) !== value.type) fail('同名边界不能使用不同类型。');
    names.set(value.name,value.type);
    return value;
  }).sort((a,b)=>a.face-b.face);
  const hasInlet = result.some(b=>b.type==='velocity-inlet'), hasOutlet = result.some(b=>b.type==='pressure-outlet');
  if (!result.some(b=>b.type==='pressure-opening') && hasInlet !== hasOutlet) fail('开放流域需要入口与出口；闭域只允许壁面。');
  return result;
}
function normalizeBoundaryDefinition(definition) {
  if (!definition || typeof definition !== 'object') fail('请先生成或导入边界。');
  const cells = count(definition.cells), faces = count(definition.faces);
  const sorted = conditions(definition.records);
  const byFace = new Map(definition.records.map(b=>[b.face,b]));
  const records = sorted.map(b => {
    const source = byFace.get(b.face);
    if (b.face >= faces) fail('面编号超过网格范围。');
    const owner = count(source.owner,true);
    if (owner >= cells) fail('单元编号超过网格范围。');
    const geometry = Object.fromEntries(['x','y','sx','sy'].map(k=>[k,finite(source[k])]));
    if (!(Math.hypot(geometry.sx,geometry.sy)>0)) fail('边界面长度必须大于零。');
    return {...b,owner,...geometry};
  });
  return {cells,faces,records};
}
function parseBoundaryDefinition(text) {
  if (typeof text !== 'string' || text.length > 32*1024*1024) fail('文件内容无效或超过32MB。');
  const lines = text.split(/\r?\n/).filter(l=>l.trim()).map(quotedTokens);
  if (lines[0]?.join(' ') !== 'CARTMESH2D_FLOW_BOUNDARIES 1' || lines[1]?.[0] !== 'COUNTS' ||
      lines[1].length !== 4 || lines.at(-1)?.join(' ') !== 'END') fail('文件格式无效。');
  const [cells,faces,n] = lines[1].slice(1).map(numberToken);
  count(n);
  if (lines.length !== n+3) fail('边界记录数不符。');
  return normalizeBoundaryDefinition({cells,faces,records:lines.slice(2,-1).map(t=>{
    if (t.length !== 12 || t[0] !== 'BOUNDARY') fail('边界行格式无效。');
    return {face:numberToken(t[1]),owner:numberToken(t[2]),x:numberToken(t[3]),y:numberToken(t[4]),sx:numberToken(t[5]),sy:numberToken(t[6]),
      type:t[7],name:t[8],u:numberToken(t[9]),v:numberToken(t[10]),p:numberToken(t[11])};
  })});
}
const quote = text => '"'+text.replace(/["\\]/g,'\\$&')+'"';
function serializeBoundaryDefinition(value) {
  const d = normalizeBoundaryDefinition(value);
  return ['CARTMESH2D_FLOW_BOUNDARIES 1',`COUNTS ${d.cells} ${d.faces} ${d.records.length}`,
    ...d.records.map(b=>`BOUNDARY ${b.face} ${b.owner} ${b.x} ${b.y} ${b.sx} ${b.sy} ${b.type} ${quote(b.name)} ${b.u} ${b.v} ${b.p}`),
    'END',''].join('\n');
}
function validateBoundaryMesh(value, mesh, speed) {
  const d = normalizeBoundaryDefinition(value);
  finite(speed);
  if (!(speed>0) || d.cells !== mesh.cells.length || d.faces !== mesh.edges.length) fail('配置与最终网格不匹配。');
  const boundaryIds = mesh.edges.flatMap((e,i)=>e.neighbour < 0 ? [i] : []);
  if (boundaryIds.length !== d.records.length) fail('未覆盖全部边界面。');
  const bound = mesh.bounds;
  const positionTolerance = 1e-12+1e-10*Math.max(bound.maxX-bound.minX,bound.maxY-bound.minY);
  for (const b of d.records) {
    const edge=mesh.edges[b.face];
    if (edge.neighbour>=0 || edge.owner!==b.owner) fail('边界所属单元或面类型不匹配。');
    const ids=mesh.cells[edge.owner].vertices;
    let a=edge.a,z=edge.b;
    if (ids[(ids.indexOf(a)+1)%ids.length]!==z) [a,z]=[z,a];
    if (ids.indexOf(a)<0 || ids[(ids.indexOf(a)+1)%ids.length]!==z) fail('无法确定边界外法向。');
    const [ax,ay]=mesh.vertices[a],[zx,zy]=mesh.vertices[z];
    const sx=zy-ay,sy=ax-zx, length=Math.hypot(sx,sy), vectorTolerance=1e-12+1e-10*length;
    if (Math.abs(b.x-(ax+zx)/2)>positionTolerance || Math.abs(b.y-(ay+zy)/2)>positionTolerance ||
        Math.abs(b.sx-sx)>vectorTolerance || Math.abs(b.sy-sy)>vectorTolerance) fail('边界位置或法向与最终网格不匹配；请重新生成边界。');
    if (b.type==='pressure-opening' && Math.min(Math.abs(sx),Math.abs(sy))>(1e-12+1e-10)*length)
      fail(`定压开口 ${b.name} 必须水平或竖直。`);
    const flux=finite(b.u*sx+b.v*sy);
    if (b.type==='velocity-inlet' && !(flux<0)) fail(`入口 ${b.name} 在面 ${b.face} 上没有指向流体内部。`);
    if (['moving-wall','smooth-moving-wall'].includes(b.type) && Math.abs(flux)>(1e-12+1e-10*Math.max(speed,Math.hypot(b.u,b.v)))*length)
      fail(`移动壁面 ${b.name} 的速度必须沿壁面切向。`);
  }
  return d;
}
// Build a pressure-driven counterpart of the native axis-aligned duct preset.
// The initial static pressure difference is 1 m2/s2 and remains editable.
function pressureDrivenBoundaryDefinition(value) {
  const d=normalizeBoundaryDefinition(value);
  if (!d.records.some(b=>b.type==='velocity-inlet') || !d.records.some(b=>b.type==='pressure-outlet'))
    fail('压差模板需要原通道入口与出口。');
  return normalizeBoundaryDefinition({...d,records:d.records.map(b=>
    ['velocity-inlet','pressure-outlet'].includes(b.type)
      ? {...b,type:'pressure-opening',u:0,v:0,p:b.type==='velocity-inlet'?1:0} : b)});
}
function sameConditions(a,b) { return JSON.stringify(conditions(a))===JSON.stringify(conditions(b)); }
module.exports={KINDS,quotedTokens,condition,conditions,normalizeBoundaryDefinition,parseBoundaryDefinition,
  serializeBoundaryDefinition,validateBoundaryMesh,pressureDrivenBoundaryDefinition,sameConditions};
