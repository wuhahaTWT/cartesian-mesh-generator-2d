'use strict';
const { validateFlowRequest } = require('./flow');
const GROUPS = ['wall','inlet','outlet','top','bottom'];
const SUFFIXES = ['.json','.vtk','.cells.csv','.faces.csv','.history.csv','.thermal-history.csv','.thermal.checkpoint','.carrier.checkpoint'];
const requireValue = (ok, message) => { if (!ok) throw new Error(`热输运：${message}`); };
const finite = (v, name) => {
  requireValue(typeof v === 'number' && Number.isFinite(v), `${name} 必须是有限数。`);
  return v;
};
const near = (a,b) => Math.abs(a-b) <= 1e-12 + 1e-9*Math.max(Math.abs(a),Math.abs(b));
function validateThermalRequest(input) {
  requireValue(input && typeof input === 'object','缺少配置。');
  const flow = validateFlowRequest({ ...input, mode:'transient' });
  const r = { ...flow, diffusivity:finite(input.diffusivity,'热扩散率'), initial:finite(input.initial,'初温'),
    source:finite(input.source,'温度源'), scalarConvection:input.scalarConvection, boundaries:{} };
  requireValue(r.diffusivity>0 && r.initial>=0,'热扩散率须为正、初温不得低于 0 K。');
  requireValue(['upwind','limited-linear'].includes(r.scalarConvection),'未知温度对流格式。');
  for (const group of GROUPS) {
    const b=input.boundaries?.[group];
    requireValue(b && ['value','flux'].includes(b.kind),`${group} 边界类型无效。`);
    r.boundaries[group]={kind:b.kind,value:finite(b.value,group),inflowValue:finite(b.inflowValue,`${group} 回流温度`)};
    requireValue(b.inflowValue>=0 && (b.kind!=='value'||b.value>=0),'温度不得低于 0 K。');
  }
  return r;
}
// Match the native flow cases. A nonrectangular boundary outside duct fails;
// assigning the closest screen side would silently invent a physical condition.
function thermalBoundaryCsv(mesh, request) {
  const r=validateThermalRequest(request), b=mesh.bounds;
  const eps=1e-10*Math.max(b.maxX-b.minX,b.maxY-b.minY);
  const equal=(x,y)=>Math.abs(x-y)<=eps;
  const rows=['face,type,value,inflowValue'];
  mesh.edges.forEach((e,id)=>{
    if(e.neighbour>=0)return;
    const a=mesh.vertices[e.a],z=mesh.vertices[e.b];
    let group;
    if(r.case==='external' && e.patch===1)group='wall';
    else if(equal(a[0],b.minX)&&equal(z[0],b.minX))group='inlet';
    else if(equal(a[0],b.maxX)&&equal(z[0],b.maxX))group='outlet';
    else if(r.case==='duct')group='wall';
    else if(equal(a[1],b.maxY)&&equal(z[1],b.maxY))group='top';
    else if(equal(a[1],b.minY)&&equal(z[1],b.minY))group='bottom';
    requireValue(group,'所选流动工况要求轴对齐矩形外边界。');
    const v=r.boundaries[group];rows.push([id,v.kind,v.value,v.inflowValue].join(','));
  });
  requireValue(rows.length>1,'没有边界面。');
  return rows.join('\n')+'\n';
}
function buildThermalInvocation(mesh,prefix,boundary,input,restart=null) {
  const r=validateThermalRequest(input);
  requireValue(typeof mesh==='string'&&mesh.endsWith('.solver.cm2d'),'必须使用最终 solver.cm2d。');
  requireValue(!r.resume||restart,'没有可用的联合续算状态。');
  const args=['--mesh',mesh,'--output',prefix,'--boundary',boundary,'--evolve-flow',r.case,
    '--flow-nu',String(r.nu),'--flow-speed',String(r.speed),'--flow-max-iterations',String(r.maxIterations),
    '--flow-convection',r.convection,'--pressure-preconditioner',r.pressurePreconditioner,'--outlet-backflow',r.outletBackflow,
    '--diffusivity',String(r.diffusivity),'--source',String(r.source),'--initial',String(r.initial),
    '--convection',r.scalarConvection,'--dt',String(r.dt),'--steps',String(r.steps)];
  if(r.resume)args.push('--restart',restart);
  return {executable:'cartmesh2d_transport_cli',request:r,args};
}
const METRICS=['step','time','flowIterations','momentumResidual','continuity','scalarIterations','scalarResidual','heatContent','globalBalance','maxCourant'];
function validateRow(row) {
  for(const k of METRICS)finite(row[k],k);
  requireValue(row.accepted===1 && row.time>0,'监测只接受完成的正时间状态。');
  for(const k of ['step','flowIterations','scalarIterations'])requireValue(Number.isSafeInteger(row[k])&&row[k]>0,'步数/迭代数无效。');
  for(const k of ['momentumResidual','continuity','scalarResidual','maxCourant'])requireValue(row[k]>=0,'残差或 CFL 为负。');
  requireValue(row.momentumResidual<1e-8&&row.continuity<1e-8,'流动未达接受条件。');
  return row;
}
function parseThermalProgress(line) {
  let r;try{r=JSON.parse(line);}catch{return null;}
  return r?.type==='thermal-time-step'?validateRow(r):null;
}
function csvRows(text,header) {
  const lines=text.trim().split(/\r?\n/);
  requireValue(lines.shift()===header,'输出 CSV 表头不匹配。');
  return lines.map(line=>{const f=line.split(',');requireValue(f.length===header.split(',').length,'CSV 列数错误。');return f;});
}
// Metadata only; native restart still validates the complete geometry and state.
function thermalCheckpointTime(text) {
  requireValue(/^CARTMESH2D_THERMAL_CHECKPOINT 1\nCOUPLING new-time-flux-Euler-v1\n/.test(text),'续算文件格式错误。');
  const parts=text.split('\nFLOW\n');
  requireValue(parts.length===2 && /^(?:CARTMESH2D_FLOW_CHECKPOINT 1|CARTMESH2D_FLOW_CHECKPOINT 2)\n/.test(parts[1]),'缺少联合流动状态。');
  const m=parts[1].match(/^TIME (\S+)$/m);
  requireValue(m && /^[+\-\d.eE]+$/.test(m[1]),'缺少物理时间。');
  const time=Number(m[1]);finite(time,'续算时间');requireValue(time>=0,'续算时间不能为负。');
  requireValue(/^FLUX \d+ .+$/m.test(parts[1]),'续算通量缺失。');
  return time;
}
function validateThermalOutput(summary,cellsText,historyText,jointText,mesh,input,startTime=0) {
  const r=validateThermalRequest(input);
  requireValue(summary?.format==='cartmesh2d-scalar-transport-v1'&&summary.status==='converged'&&summary.converged===true&&summary.evolvingFlow===true,'不是完整同步热计算结果。');
  for(const key of ['time','acceptedTime','carrierTime','timeStep','diffusivity','flowNu','flowSpeed','constantSource','initialValue','minValue','maxValue','globalBalance','residualNorm','maxDiagonalScaledImbalance'])finite(summary[key],key);
  requireValue(summary.cells===mesh.cells.length&&summary.faces===mesh.edges.length&&summary.steps===r.steps,'网格数量/步数不一致。');
  for(const [key,value] of Object.entries({timeStep:r.dt,diffusivity:r.diffusivity,flowNu:r.nu,flowSpeed:r.speed,constantSource:r.source,initialValue:r.initial}))requireValue(near(summary[key],value),`${key} 与请求不符。`);
  requireValue(summary.flowCase===r.case&&summary.convection===r.scalarConvection&&summary.flowConvection===r.convection
    && (summary.outletBackflow===undefined ? 'reject' : summary.outletBackflow)===r.outletBackflow,'物理工况/格式不一致。');
  requireValue(summary.maxDiagonalScaledImbalance<=1e-9,'温度单元失衡未达停止条件。');
  const t=startTime+r.dt*r.steps;
  for(const key of ['time','acceptedTime','carrierTime'])requireValue(near(summary[key],t),'流动与温度物理时间不同步。');
  requireValue(near(thermalCheckpointTime(jointText),t),'联合保存时间不同步。');
  const scalarLine=jointText.split('\nFLOW\n')[0].split('\n').find(l=>l.startsWith('SCALAR '));
  requireValue(scalarLine,'缺少联合温度场。');
  const jointValues=scalarLine.trim().split(/\s+/).slice(1).map(Number);
  requireValue(jointValues.shift()===mesh.cells.length&&jointValues.length===mesh.cells.length,'联合温度场数量错误。');
  const rows=csvRows(cellsText,'cell,x,y,area,value,previous,sourceIntegral,temporalIntegral,exact');
  requireValue(rows.length===mesh.cells.length,'温度单元数错误。');
  let min=Infinity,max=-Infinity,heat=0;
  const cells=rows.map((row,i)=>{
    requireValue(row.slice(0,8).every(v=>v.trim()!==''&&Number.isFinite(Number(v))),'温度 CSV 含非法值。');
    const [id,x,y,area,theta]=row.map(Number);
    requireValue(id===i&&near(area,mesh.cells[i].area)&&area>0,'温度单元 ID 或面积错误。');
    requireValue(Number.isFinite(jointValues[i])&&theta===jointValues[i],'温度 CSV 与联合状态不同。');
    min=Math.min(min,theta);max=Math.max(max,theta);heat+=area*theta;return{id,theta};
  });
  requireValue(near(min,summary.minValue)&&near(max,summary.maxValue),'温度范围不一致。');
  const history=csvRows(historyText,'step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant').map((row,i)=>{
    requireValue(row.every(v=>v.trim()!==''&&Number.isFinite(Number(v))),'时间历史含非法值。');
    const [step,time,accepted,flowIterations,momentumResidual,continuity,scalarIterations,scalarResidual,heatContent,globalBalance,maxCourant]=row.map(Number);
    requireValue(step===i+1&&near(time,startTime+(i+1)*r.dt),'时间历史次序错误。');
    return validateRow({step,time,accepted,flowIterations,momentumResidual,continuity,scalarIterations,scalarResidual,heatContent,globalBalance,maxCourant});
  });
  requireValue(history.length===r.steps,'时间历史不完整。');
  requireValue(near(history.at(-1).heatContent,heat)&&near(history.at(-1).globalBalance,summary.globalBalance)&&near(history.at(-1).scalarResidual,summary.residualNorm),'历史、场与摘要不一致。');
  return {summary:{...summary,dt:r.dt},fields:{cells},history,request:r};
}
module.exports={GROUPS,SUFFIXES,validateThermalRequest,thermalBoundaryCsv,buildThermalInvocation,parseThermalProgress,thermalCheckpointTime,validateThermalOutput};
