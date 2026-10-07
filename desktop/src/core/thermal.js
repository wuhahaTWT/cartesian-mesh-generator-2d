'use strict';
const { validateFlowRequest } = require('./flow');
const {validateInitialVortexOutput}=require('./initial-vortex');
const GROUPS = ['wall','inlet','outlet','top','bottom'];
const SUFFIXES = ['.json','.vtk','.cells.csv','.faces.csv','.history.csv','.thermal-history.csv','.thermal.checkpoint','.carrier.checkpoint','.heat-history.csv','.boundary-heat-history.csv'];
const requireValue = (ok, message) => { if (!ok) throw new Error(`热输运：${message}`); };
const finite = (v, name) => {
  requireValue(typeof v === 'number' && Number.isFinite(v), `${name} 必须是有限数。`);
  return v;
};
const near = (a,b) => Math.abs(a-b) <= 1e-12 + 1e-9*Math.max(Math.abs(a),Math.abs(b));
function validateThermalRequest(input) {
  requireValue(input && typeof input === 'object','缺少配置。');
  requireValue(input.case!=='custom','命名边界的温度配置尚未支持。');
  const flow = validateFlowRequest({ ...input, mode:input.mode==='adaptive'?'adaptive':'transient' });
  requireValue(flow.linearPolicy==='strict' && flow.velocityRelaxation===.6 && flow.pressureCorrectionPasses===4,
    '桌面温度联算仍使用固定线性精度、默认速度松弛和4次压力校正；这些加速设置仅用于独立层流。');
  // Preserve the tighter historical coupled-flow default even when standalone
  // flow is configured at 1e-6. A smaller user tolerance also applies here.
  flow.tolerance = Math.min(flow.tolerance, 1e-8);
  const r = { ...flow, diffusivity:finite(input.diffusivity,'热扩散率'), initial:finite(input.initial,'初温'),
    source:finite(input.source,'温度源'), scalarConvection:input.scalarConvection, fluxCorrection:input.fluxCorrection ?? 'unrestricted', boundaries:{} };
  requireValue(r.diffusivity>0 && r.initial>=0,'热扩散率须为正、初温不得低于 0 K。');
  requireValue(['upwind','limited-linear'].includes(r.scalarConvection),'未知温度对流格式。');
  requireValue(['unrestricted','bounded','bounded-spatial'].includes(r.fluxCorrection),'未知温度通量修正。');
  for (const group of GROUPS) {
    const b=input.boundaries?.[group];
    requireValue(b && ['value','flux'].includes(b.kind),`${group} 边界类型无效。`);
    r.boundaries[group]={kind:b.kind,value:finite(b.value,group),inflowValue:finite(b.inflowValue,`${group} 回流温度`)};
    requireValue(b.inflowValue>=0 && (b.kind!=='value'||b.value>=0),'温度不得低于 0 K。');
  }
  requireValue(input.events===undefined || Array.isArray(input.events),'时间事件必须为列表。');
  r.events=(input.events || []).map((event,index)=>{
    requireValue(event && typeof event==='object',`事件 ${index+1} 缺少配置。`);
    const time=finite(event.time,'事件时间'),value=finite(event.value,'事件值');
    requireValue(time>=0,'事件时间不得为负。');
    requireValue(['source',...GROUPS].includes(event.target),'未知事件目标。');
    if(event.target==='source') {
      requireValue(event.kind==='source','体积热源事件须使用热源类型。');
      return {time,target:'source',kind:'source',value};
    }
    requireValue(['value','flux'].includes(event.kind),'事件边界类型无效。');
    const inflowValue=finite(event.inflowValue,'事件流入温度');
    requireValue(inflowValue>=0 && (event.kind!=='value'||value>=0),'事件温度不得低于 0 K。');
    return {time,target:event.target,kind:event.kind,value,inflowValue};
  }).sort((a,b)=>a.time-b.time || ['source',...GROUPS].indexOf(a.target)-['source',...GROUPS].indexOf(b.target));
  for(let i=1;i<r.events.length;++i)requireValue(r.events[i].time!==r.events[i-1].time || r.events[i].target!==r.events[i-1].target,'同一时刻不能重复设置同一事件目标。');
  requireValue(!r.events.length || r.mode==='adaptive','时间事件须使用自动步长模式，以精确停在变化时刻。');
  r.timeError=Boolean(input.timeError);
  r.temperatureScale=finite(input.temperatureScale ?? 1,'温升尺度');
  r.timeRtol=finite(input.timeRtol ?? .01,'时间误差相对容差');
  requireValue(r.temperatureScale>0&&r.timeRtol>0,'温升尺度和时间误差容差须为正。');
  return r;
}
// Match the native flow cases. A nonrectangular boundary outside duct fails;
// assigning the closest screen side would silently invent a physical condition.
function thermalBoundaryGroups(mesh, r) {
  const b=mesh.bounds,eps=1e-10*Math.max(b.maxX-b.minX,b.maxY-b.minY);
  const equal=(x,y)=>Math.abs(x-y)<=eps;
  const groups=Object.fromEntries(GROUPS.map(name=>[name,[]]));
  mesh.edges.forEach((e,id)=>{
    if(e.neighbour>=0)return;
    const a=mesh.vertices[e.a],z=mesh.vertices[e.b];let group;
    if(r.case==='external' && e.patch===1)group='wall';
    else if(equal(a[0],b.minX)&&equal(z[0],b.minX))group='inlet';
    else if(equal(a[0],b.maxX)&&equal(z[0],b.maxX))group='outlet';
    else if(r.case==='duct')group='wall';
    else if(equal(a[1],b.maxY)&&equal(z[1],b.maxY))group='top';
    else if(equal(a[1],b.minY)&&equal(z[1],b.minY))group='bottom';
    requireValue(group,'所选流动工况要求轴对齐矩形外边界。');groups[group].push(id);
  });
  return groups;
}
function thermalBoundaryCsv(mesh, request) {
  const r=validateThermalRequest(request),groups=thermalBoundaryGroups(mesh,r),records=[];
  for(const group of GROUPS)for(const id of groups[group]) {
    const v=r.boundaries[group];records.push([id,v.kind,v.value,v.inflowValue]);
  }
  requireValue(records.length>0,'没有边界面。');
  records.sort((a,b)=>a[0]-b[0]);
  return 'face,type,value,inflowValue\n'+records.map(row=>row.join(',')).join('\n')+'\n';
}
function thermalEventsCsv(mesh, request) {
  const r=validateThermalRequest(request),groups=thermalBoundaryGroups(mesh,r);
  const rows=['time,target,type,value,inflowValue'];
  for(const event of r.events) {
    if(event.target==='source')rows.push([event.time,'source','source',event.value,''].join(','));
    else {
      requireValue(groups[event.target].length>0,`事件目标 ${event.target} 在当前网格/工况中没有边界面。`);
      for(const id of groups[event.target])rows.push([event.time,`face:${id}`,event.kind,event.value,event.inflowValue].join(','));
    }
  }
  return rows.join('\n')+'\n';
}
function buildThermalInvocation(mesh,prefix,boundary,input,restart=null,events=null) {
  const r=validateThermalRequest(input);
  requireValue(typeof mesh==='string'&&mesh.endsWith('.solver.cm2d'),'必须使用最终 solver.cm2d。');
  requireValue(!r.resume||restart,'没有可用的联合续算状态。');
  const args=['--mesh',mesh,'--output',prefix,'--boundary',boundary,'--evolve-flow',r.case,
    '--flow-nu',String(r.nu),'--flow-speed',String(r.speed),'--flow-max-iterations',String(r.maxIterations),
    '--flow-tolerance',String(r.tolerance),
    '--flow-convection',r.convection,'--pressure-preconditioner',r.pressurePreconditioner,'--outlet-backflow',r.outletBackflow,
    '--diffusivity',String(r.diffusivity),'--source',String(r.source),'--initial',String(r.initial),
    '--convection',r.scalarConvection,'--flux-correction',r.fluxCorrection,'--dt',String(r.dt)];
  if(r.mode==='adaptive')args.push('--end-time',String(r.endTime),'--min-dt',String(r.minDt),
    '--max-courant',String(r.maxCourant),'--max-step-retries',String(r.maxRetries),'--max-time-steps',String(r.maxSteps),
    '--time-error',r.timeError?'on':'off','--temperature-scale',String(r.temperatureScale),
    '--velocity-scale',String(r.speed),'--time-rtol',String(r.timeRtol));
  else args.push('--steps',String(r.steps));
  if(r.initialVortex) {
    const v=r.initialVortex;args.push('--initial-vortex-x',String(v.centre[0]),'--initial-vortex-y',String(v.centre[1]),
      '--initial-vortex-radius',String(v.radius),'--initial-vortex-speed',String(v.peakSpeed));
  }
  if(r.events.length) {requireValue(typeof events==='string'&&events.length>0,'缺少完整时间事件文件。');args.push('--thermal-events',events);}
  if(r.resume)args.push('--restart',restart);
  return {executable:'cartmesh2d_transport_cli',request:r,args};
}
const METRICS=['step','time','flowIterations','momentumResidual','continuity','scalarIterations','scalarResidual','heatContent','globalBalance','maxCourant'];
function validateRow(row,tolerance=1e-8) {
  for(const k of METRICS)finite(row[k],k);
  requireValue(row.accepted===1 && row.time>0,'监测只接受完成的正时间状态。');
  for(const k of ['step','flowIterations','scalarIterations'])requireValue(Number.isSafeInteger(row[k])&&row[k]>0,'步数/迭代数无效。');
  for(const k of ['momentumResidual','continuity','scalarResidual','maxCourant'])requireValue(row[k]>=0,'残差或 CFL 为负。');
  requireValue(row.momentumResidual<tolerance&&row.continuity<1e-8,'流动未达接受条件。');
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
function thermalCheckpointTime(text, events=[]) {
  // Native text streams use CRLF on Windows; line endings are not part of the
  // physical checkpoint identity. Keep the same version and field checks.
  text=text.replace(/\r\n/g,'\n');
  requireValue(/^CARTMESH2D_THERMAL_CHECKPOINT [1234]\nCOUPLING new-time-flux-Euler-v1\n/.test(text) &&
    (!text.startsWith('CARTMESH2D_THERMAL_CHECKPOINT 2') || events.length>0),'续算文件格式错误。');
  if(/^CARTMESH2D_THERMAL_CHECKPOINT [234]\n/.test(text)) {
    const count=text.match(/^EVENTS (\d+)$/m),times=[...text.matchAll(/^EVENT (\S+)$/gm)].map(m=>Number(m[1]));
    const expected=[...new Set(events.map(event=>event.time))];
    requireValue(count && Number(count[1])===expected.length && times.length===expected.length && times.every((time,i)=>time===expected[i]),
      '时间事件与保存的完整规律不一致；须保留原 desktop-state.json 或项目包。');
  } else requireValue(events.length===0,'旧状态缺少时间事件。');
  if(text.startsWith('CARTMESH2D_THERMAL_CHECKPOINT 4'))requireValue(/^THERMAL_CONFIG \S+ (?:upwind|limited-linear) (?:bounded|bounded-spatial)$/m.test(text)&&/^CONTROLLER_PRESENT [01]$/m.test(text),'有界温度状态配置缺失。');
  const parts=text.split('\nFLOW\n');
  requireValue(parts.length===2 && /^(?:CARTMESH2D_FLOW_CHECKPOINT 1|CARTMESH2D_FLOW_CHECKPOINT 2)\n/.test(parts[1]),'缺少联合流动状态。');
  const m=parts[1].match(/^TIME (\S+)$/m);
  requireValue(m && /^[+\-\d.eE]+$/.test(m[1]),'缺少物理时间。');
  const time=Number(m[1]);finite(time,'续算时间');requireValue(time>0,'初值不是已接受的物理续算状态。');
  requireValue(/^FLUX \d+ .+$/m.test(parts[1]),'续算通量缺失。');
  return time;
}
// Read the native accepted-step stream. A crashed run may end in a partial
// final line; completed projects still require the entire stream and endpoint.
function readThermalHistory(text,input,startTime,endTime,{partial=false}={}) {
  const r=validateThermalRequest(input);
  const truncated=partial && !text.endsWith('\n');
  if(truncated)text=text.slice(0,text.lastIndexOf('\n')+1);
  if(partial && !text.trim())return {rows:[],truncated:true};
  let previous=startTime;
  const rows=csvRows(text,'step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant').map((row,i)=>{
    requireValue(row.every(v=>v.trim()!==''&&Number.isFinite(Number(v))),'时间历史含非法值。');
    const [step,time,accepted,flowIterations,momentumResidual,continuity,scalarIterations,scalarResidual,heatContent,globalBalance,maxCourant]=row.map(Number);
    requireValue(step===i+1 && time>previous && time<=endTime,'时间历史次序错误。');
    requireValue(r.mode==='adaptive' ? time-previous<=r.dt*(1+1e-9) : near(time,startTime+(i+1)*r.dt),'接受时钟或最大步长错误。');
    previous=time;
    return validateRow({step,time,accepted,flowIterations,momentumResidual,continuity,scalarIterations,scalarResidual,heatContent,globalBalance,maxCourant},r.tolerance);
  });
  if(!partial)requireValue(rows.length>0&&rows.at(-1).time===endTime,'时间历史未到已接受终点。');
  return {rows,truncated};
}
// Failure diagnostics report the native rejected attempt. They neither accept
// its fields nor infer a restart checkpoint from an initial/attempted clock.
function thermalFailureMessage(summary,startTime=0) {
  const time=Number.isFinite(summary?.acceptedTime)?summary.acceptedTime:startTime;
  const lines=[`热计算未完成；已接受到 t=${time} s。`],c=summary?.controller;
  const reasons={flow:'流动方程未收敛',scalar:'温度方程未收敛',
    carrier:'载流面通量未满足温度连续性要求','carrier-half':'半步载流面通量未满足温度连续性要求',
    'flow-half':'时间估计的半步流动方程未收敛','scalar-half':'时间估计的半步温度方程未收敛',
    'flow-linear':'流动线性求解达到迭代上限','flow-linear-half':'时间估计的半步流动线性求解达到迭代上限',
    'nonfinite-error':'时间缺陷含非有限数'};
  const number=value=>Number(value.toPrecision(6)).toString();
  if(c?.reason==='time-error'&&Number.isFinite(c.errorRatio))
    lines.push(`联合时间缺陷为 ${number(c.errorRatio)} 倍预算，接受上限为 1。`);
  else if(c?.reason==='courant'&&Number.isFinite(c.courant)&&Number.isFinite(c.maximumCourant))
    lines.push(`CFL 为 ${number(c.courant)}，超过上限 ${number(c.maximumCourant)}。`);
  else if(reasons[c?.reason || summary?.failedStage])lines.push(reasons[c?.reason || summary.failedStage]+'。');
  if(typeof c?.diagnostic==='string'&&c.diagnostic)lines.push(c.diagnostic);
  if(c&&Number.isFinite(c.timeStep)&&Number.isFinite(c.minimumTimeStep)) {
    lines.push(`最后尝试步长 ${number(c.timeStep)} s；设定最小步长 ${number(c.minimumTimeStep)} s。`);
    if(c.timeStep<=c.minimumTimeStep)
      lines.push('步长已到下限；可减小设定的最小步长后重试，保留当前精度要求。');
  }
  if(c&&Number.isSafeInteger(c.attempts)&&Number.isSafeInteger(c.maximumRetries)&&c.attempts-1>=c.maximumRetries)
    lines.push(`重试预算已用尽（${c.attempts} 次尝试）；可增加重试次数后重试。`);
  return lines.join('\n');
}
function validateThermalOutput(summary,cellsText,historyText,jointText,mesh,input,startTime=0) {
  const r=validateThermalRequest(input);
  requireValue(summary?.format==='cartmesh2d-scalar-transport-v1'&&summary.status==='converged'&&summary.converged===true&&summary.evolvingFlow===true,'不是完整同步热计算结果。');
  for(const key of ['time','acceptedTime','carrierTime','timeStep','diffusivity','flowNu','flowSpeed','constantSource','initialValue','minValue','maxValue','globalBalance','residualNorm','maxDiagonalScaledImbalance'])finite(summary[key],key);
  requireValue(summary.cells===mesh.cells.length&&summary.faces===mesh.edges.length&&(r.mode==='adaptive'||summary.steps===r.steps),'网格数量/步数不一致。');
  for(const [key,value] of Object.entries({...(r.mode==='adaptive'?{maximumTimeStep:r.dt}:{timeStep:r.dt}),diffusivity:r.diffusivity,flowNu:r.nu,flowSpeed:r.speed,constantSource:r.source,initialValue:r.initial}))requireValue(near(summary[key],value),`${key} 与请求不符。`);
  requireValue(summary.flowCase===r.case&&summary.convection===r.scalarConvection&&summary.flowConvection===r.convection
    && (summary.outletBackflow===undefined ? 'reject' : summary.outletBackflow)===r.outletBackflow,'物理工况/格式不一致。');
  requireValue((summary.fluxCorrection ?? 'unrestricted')===r.fluxCorrection,'温度通量修正与请求不符。');
  validateInitialVortexOutput(summary,r,startTime);
  if(r.fluxCorrection!=='unrestricted') {
    for(const key of ['lowerBound','upperBound','maxBoundViolation'])finite(summary[key],key);
    requireValue(summary.lowerBound<=summary.upperBound&&summary.maxBoundViolation>=0&&summary.maxBoundViolation<=1e-9,'有界温度场未达停止条件。');
    const savedMode=jointText.replace(/\r\n/g,'\n').match(/^THERMAL_CONFIG \S+ \S+ (\S+)$/m)?.[1];
    requireValue(jointText.startsWith('CARTMESH2D_THERMAL_CHECKPOINT 4')&&savedMode===r.fluxCorrection,'联合状态未绑定请求的有界温度格式。');
  }
  requireValue(summary.flowTolerance===r.tolerance,'流动停止容差与请求不符。');
  requireValue(summary.maxDiagonalScaledImbalance<=1e-9,'温度单元失衡未达停止条件。');
  const adaptive=r.mode==='adaptive';
  const t=adaptive?r.endTime:startTime+r.dt*r.steps;
  if(adaptive) {requireValue(near(summary.temperatureScale,r.temperatureScale)&&near(summary.timeRelativeTolerance,r.timeRtol),'时间误差尺度与请求不符。');}
  if(adaptive)requireValue(summary.timeStepControl===(r.timeError?'joint-cfl-be-error-retry':'joint-cfl-retry'),'联合时间控制模式不符。');
  for(const key of ['time','acceptedTime','carrierTime'])requireValue(near(summary[key],t),'流动与温度物理时间不同步。');
  if(r.events.length)requireValue(summary.thermalEventCount===new Set(r.events.map(e=>e.time)).size,'原生时间事件数量与请求不符。');
  requireValue(near(thermalCheckpointTime(jointText,r.events),t),'联合保存时间不同步。');
  const scalarLine=jointText.split(/\r?\nFLOW\r?\n/)[0].split(/\r?\n/).find(l=>l.startsWith('SCALAR '));
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
  const history=readThermalHistory(historyText,r,startTime,summary.time).rows;
  requireValue(history.length===(adaptive?summary.completedSteps:r.steps),'时间历史不完整。');
  requireValue(near(history.at(-1).heatContent,heat)&&near(history.at(-1).globalBalance,summary.globalBalance)&&near(history.at(-1).scalarResidual,summary.residualNorm),'历史、场与摘要不一致。');
  return {summary:{...summary,dt:r.dt},fields:{cells},history,request:r};
}
module.exports={GROUPS,SUFFIXES,validateThermalRequest,thermalBoundaryCsv,thermalEventsCsv,buildThermalInvocation,parseThermalProgress,thermalCheckpointTime,readThermalHistory,validateThermalOutput,thermalFailureMessage};
