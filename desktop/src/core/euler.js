'use strict';
const SUFFIXES=['.json','.fields.json','.cells.csv','.faces.csv','.history.csv','.vtk','.checkpoint','.boundaries'];
const PHYSICAL=['case','density','u','v','pressure','gamma','gasConstant','split','thermalConductivity','wallThermal','wallValue'];
const NUMERICAL=['endTime','maximumStep','minimumStep','cfl','maximumSteps','maximumSeconds'];
const requireValue=(ok,message)=>{if(!ok)throw new Error(`可压 Euler：${message}`);};
const finite=(v,name)=>{requireValue(typeof v==='number'&&Number.isFinite(v),`${name} 必须是有限数。`);return v;};
const near=(a,b)=>Math.abs(a-b)<=2e-12+2e-10*Math.max(Math.abs(a),Math.abs(b));
function validateEulerRequest(input) {
  requireValue(input&&typeof input==='object'&&!Array.isArray(input),'缺少配置。');
  requireValue(Object.keys(input).every(k=>[...PHYSICAL,...NUMERICAL,'fluxScheme','order','resume'].includes(k)),'存在未知配置。');
  requireValue(['sod','external','uniform','sealed'].includes(input.case),'未知工况。');
  const r={case:input.case,fluxScheme:input.fluxScheme??'rusanov',order:input.order??1,thermalConductivity:input.thermalConductivity??0,wallThermal:input.wallThermal??'insulated',wallValue:input.wallValue??0};
  requireValue(['rusanov','hllc'].includes(r.fluxScheme)&&[1,2].includes(r.order),'通量格式或精度阶数无效。');
  for(const k of [...PHYSICAL.slice(1).filter(k=>!['thermalConductivity','wallThermal','wallValue'].includes(k)),...NUMERICAL])r[k]=finite(input[k],k);
  finite(r.thermalConductivity,'导热系数');finite(r.wallValue,'热边界值');
  requireValue(r.thermalConductivity>=0&&['insulated','temperature','flux'].includes(r.wallThermal),'导热系数或热壁类型无效。');
  requireValue(r.wallThermal!=='temperature'||r.wallValue>0,'定温壁需要正的绝对温度 K。');
  requireValue(r.wallThermal!=='insulated'||r.wallValue===0,'绝热壁数值必须为零。');
  requireValue(r.wallThermal==='insulated'||(r.thermalConductivity>0&&r.case!=='uniform'),'热壁需要正导热系数和有壁面工况。');
  requireValue(r.density>0&&r.pressure>0&&r.gamma>1&&r.gasConstant>0,'密度、绝对压力、气体常数须为正，gamma须大于1。');
  requireValue(r.split>0&&r.split<1&&(r.case==='sod'||r.split===.5),'隔膜位置只适用于冲击管。');
  requireValue(r.case!=='sod'||(r.u===0&&r.v===0),'冲击管初始速度必须为零。');
  requireValue(r.endTime>0&&r.minimumStep>0&&r.maximumStep>=r.minimumStep&&r.cfl>0&&r.cfl<=.45,'时间步与声学/导热组合 CFL 无效。');
  requireValue(Number.isSafeInteger(r.maximumSteps)&&r.maximumSteps>=1&&r.maximumSteps<=1000000,'接受步数须为1至1000000。');
  requireValue(r.maximumSeconds>0&&r.maximumSeconds<=3600,'单次时间预算须为0至3600秒。');
  requireValue(input.resume===undefined||typeof input.resume==='boolean','续算开关须为布尔值。');r.resume=input.resume===true;
  return r;
}
function buildEulerInvocation(mesh,prefix,input,restart=null) {
  const r=validateEulerRequest(input);requireValue(typeof mesh==='string'&&mesh.endsWith('.solver.cm2d')&&!mesh.endsWith('.failed.solver.cm2d'),'需要最终求解网格。');
  requireValue(!r.resume||restart,'没有可用的已接受状态。');
  const args=['--mesh',mesh,'--output',prefix,'--case',r.case,'--flux',r.fluxScheme,'--order',String(r.order),'--wall-thermal',r.wallThermal];
  for(const [key,flag] of Object.entries({thermalConductivity:'conductivity',wallValue:'wall-value',density:'density',u:'u',v:'v',pressure:'pressure',gamma:'gamma',gasConstant:'gas-r',split:'split',endTime:'end-time',maximumStep:'max-step',minimumStep:'min-step',cfl:'cfl',maximumSteps:'max-steps',maximumSeconds:'max-seconds'}))args.push('--'+flag,String(r[key]));
  if(r.resume)args.push('--restart',restart);
  return {request:r,executable:'cartmesh2d_euler_cli',args};
}
function conservativePrimitive(q,gamma) {
  q.forEach(v=>finite(v,'守恒状态'));requireValue(q.length===4&&q[0]>0,'密度状态错误。');
  const u=q[1]/q[0],v=q[2]/q[0],p=(gamma-1)*(q[3]-.5*(q[1]*u+q[2]*v));
  requireValue([u,v,p].every(Number.isFinite)&&p>0,'非正压力或数值范围错误。');return {rho:q[0],u,v,p};
}
// Metadata/state inspection only; the native restart also compares its full
// serialized Fv geometry and every physical boundary before advancing.
function eulerCheckpoint(text,mesh,input) {
  const r=validateEulerRequest(input),lines=text.trimEnd().split(/\r?\n/);
  const version=r.thermalConductivity>0?2:1;requireValue(lines[0]==='CM2D_EULER_CHECKPOINT '+version,'检查点导热物理/格式不符。');
  if(version===2) {
    const conductivity=lines.filter(line=>line.startsWith('CONDUCTIVITY '));requireValue(conductivity.length===1&&Number(conductivity[0].split(' ')[1])===r.thermalConductivity,'检查点导热系数不匹配。');
    const index=lines.findIndex(line=>line.startsWith('BOUNDARIES ')),count=Number(lines[index]?.split(' ')[1]);requireValue(index>=0&&Number.isSafeInteger(count)&&count>0,'检查点热边界缺失。');
    for(const line of lines.slice(index+1,index+1+count)) {
      const tokens=line.match(/"(?:[^"\\]|\\.)*"|\S+/g);requireValue(tokens?.length===10,'检查点热边界行错误。');
      const wall=Number(tokens[1])===0,kind=wall?['insulated','temperature','flux'].indexOf(r.wallThermal):0,value=wall?r.wallValue:0;
      requireValue(Number(tokens[8])===kind&&Number(tokens[9])===value,'检查点热壁与请求不匹配。');
    }
  }
  const gas=lines[1]?.split(' ').map((s,i)=>i?Number(s):s);
  requireValue(gas?.length===3&&gas[0]==='GAS'&&gas[1]===r.gamma&&gas[2]===r.gasConstant,'检查点气体不匹配。');
  requireValue(lines[2]==='PROBLEM '+JSON.stringify(r.case)&&lines[3]==='CELLS '+mesh.cells.length,'检查点工况/网格数量不匹配。');
  const matches=lines.map((s,i)=>s.startsWith('STATE ')?i:-1).filter(i=>i>=0);requireValue(matches.length===1,'检查点状态段错误。');
  const index=matches[0],header=lines[index].split(' '),time=Number(header[1]),steps=Number(header[2]);
  requireValue(header.length===3&&Number.isFinite(time)&&time>=0&&Number.isSafeInteger(steps)&&steps>=0,'检查点时间/步数错误。');
  requireValue(lines.length===index+mesh.cells.length+2&&lines.at(-1)==='END','检查点状态截断或有额外内容。');
  const cells=lines.slice(index+1,-1).map(line=>line.split(' ').map(Number));cells.forEach(q=>conservativePrimitive(q,r.gamma));
  return {time,steps,cells};
}
function parseEulerProgress(line) {
  let r;try{r=JSON.parse(line);}catch{return null;}
  if(r?.type!=='euler-step')return null;
  for(const k of ['step','time','acousticCourant','minimumDensity','minimumPressure','mass','totalEnergy'])finite(r[k],k);
  requireValue(Number.isSafeInteger(r.step)&&r.step>0&&r.time>0&&r.acousticCourant>0&&r.acousticCourant<=.45*(1+1e-12)&&r.minimumDensity>0&&r.minimumPressure>0,'进度状态无效。');
  if(r.thermalCourant!==undefined||r.combinedCourant!==undefined) {
    finite(r.thermalCourant,'导热 CFL');finite(r.combinedCourant,'组合 CFL');
    requireValue(r.thermalCourant>=0&&r.combinedCourant>=Math.max(r.thermalCourant,r.acousticCourant)-2e-15&&r.combinedCourant<=.45*(1+1e-12),'进度导热 CFL 无效。');
  }
  return r;
}
function rows(text,required) {
  const lines=text.trim().split(/\r?\n/),header=lines.shift().split(',');
  requireValue(required.every(k=>header.includes(k))&&new Set(header).size===header.length,'CSV 表头错误。');
  return lines.map(line=>{const fields=line.split(',');requireValue(fields.length===header.length,'CSV 列数错误。');return Object.fromEntries(header.map((key,i)=>[key,fields[i]]));});
}
function numbers(row,keys) {return keys.map(k=>{requireValue(row[k]!==undefined&&row[k].trim()!=='','CSV 缺少数值。');return finite(Number(row[k]),k);});}
function cellGeometry(mesh,cell) {
  const o=mesh.vertices[cell.vertices[0]];let twiceArea=0,x=0,y=0;
  for(let j=0;j<cell.vertices.length;j++) {
    const a=mesh.vertices[cell.vertices[j]],b=mesh.vertices[cell.vertices[(j+1)%cell.vertices.length]];
    const ax=a[0]-o[0],ay=a[1]-o[1],bx=b[0]-o[0],by=b[1]-o[1],cross=ax*by-ay*bx;
    twiceArea+=cross;x+=(ax+bx)*cross;y+=(ay+by)*cross;
  }
  requireValue(twiceArea>0,'非正面积。');return {area:twiceArea/2,x:o[0]+x/(3*twiceArea),y:o[1]+y/(3*twiceArea)};
}
function validateEulerOutput(summary,fields,cellsText,facesText,historyText,checkpointText,mesh,input,startTime=0) {
  const r=validateEulerRequest(input),s=summary;
  const method=(r.order===1?'first-order ':'limited-linear ')+(r.fluxScheme==='hllc'?'HLLC-HLLE':'Rusanov')+(r.order===1?' / forward Euler':' / SSPRK2');
  requireValue(s?.solver==='native 2D ideal-gas Euler'&&s.method===method&&(s.fluxScheme??'rusanov')===r.fluxScheme&&(s.order??1)===r.order,'求解模型或数值格式不匹配。');
  const modern=s.fluxScheme!==undefined,thermal=s.thermalConductivity!==undefined;
  requireValue((s.thermalConductivity??0)===r.thermalConductivity&&(s.wallThermal??'insulated')===r.wallThermal&&(s.wallValue??0)===r.wallValue,'导热物性或热壁与请求不同。');
  if(thermal)requireValue(s.heatDiscretization==='Fourier / least-squares corrected / full-row-norm explicit bound'&&Number.isSafeInteger(s.heatNonMonotoneRows)&&s.heatNonMonotoneRows>=0&&s.heatNonMonotoneRows<=s.cells,'导热离散诊断错误。');
  if(modern) {
    requireValue(s.shockControl===(r.fluxScheme==='hllc'?'multidimensional pressure-ratio cube HLLC/HLLE blend':'none'),'激波保护格式不匹配。');
    for(const key of ['minimumContactRestoration','lastMinimumContactRestoration'])requireValue(finite(s[key],key)>=0&&s[key]<=1&&(r.fluxScheme==='hllc'||s[key]===1),'激波保护权重无效。');
  }
  const count=(value,name)=>{requireValue(Number.isSafeInteger(value)&&value>=0,`${name} 计数无效。`);return value;};
  for(const key of ['hllcFallbackEvaluations','lastHllcFallbackEvaluations','reconstructionFallbackCells','lastReconstructionFallbackCells'])if(modern)count(s[key],key);
  requireValue(s.status==='target_reached'&&s.targetReached===true&&s.failure==='','未到达目标物理时间。');
  requireValue(s.cells===mesh.cells.length&&s.faces===mesh.edges.length&&s.case===r.case,'网格或工况不匹配。');
  for(const [key,value] of Object.entries({gamma:r.gamma,gasConstant:r.gasConstant,split:r.split,requestedEndTime:r.endTime,time:r.endTime,initialTime:startTime,cflLimit:r.cfl,maximumStep:r.maximumStep,minimumStep:r.minimumStep,maximumSteps:r.maximumSteps,maximumSeconds:r.maximumSeconds}))requireValue(s[key]===value,`${key} 与请求不同。`);
  for(const [key,value] of Object.entries({rho:r.density,u:r.u,v:r.v,p:r.pressure}))requireValue(s.referenceState?.[key]===value,'参考状态与请求不同。');
  requireValue(Number.isSafeInteger(s.acceptedSteps)&&s.acceptedSteps>0&&s.acceptedSteps<=r.maximumSteps&&Number.isSafeInteger(s.steps)&&s.steps>=s.acceptedSteps&&s.lastStep>0,'接受时间步无效。');
  finite(s.lastStep,'最后时间步');
  requireValue(fields?.format==='cartmesh2d-euler-v1'&&fields.time===s.time&&fields.cells?.length===s.cells,'显示场格式或时间不匹配。');
  const checkpoint=eulerCheckpoint(checkpointText,mesh,r);requireValue(checkpoint.time===s.time&&checkpoint.steps===s.steps,'检查点终态不匹配。');
  const raw=rows(cellsText,['cell','x','y','area','rho','u','v','p','rhoU','rhoV','rhoE','temperature','mach','previousRho','previousRhoU','previousRhoV','previousRhoE']);
  requireValue(raw.length===s.cells,'单元 CSV 数量不匹配。');
  const current=[],previous=[],areas=[],heatRates=[];
  for(let i=0;i<raw.length;i++) {
    const row=raw[i],field=fields.cells[i],geometry=cellGeometry(mesh,mesh.cells[i]);
    requireValue(Number(row.cell)===i&&field.id===i,'单元 ID 错位。');
    requireValue(field.rho>0&&field.p>0&&field.temperature>0&&field.mach>=0&&Number(row.p)>0&&Number(row.temperature)>0&&Number(row.mach)>=0,'显示状态非正。');
    for(const k of ['x','y','area'])requireValue(near(Number(row[k]),geometry[k]),'单元几何不匹配。');areas.push(geometry.area);
    const q=numbers(row,['rho','rhoU','rhoV','rhoE']),old=numbers(row,['previousRho','previousRhoU','previousRhoV','previousRhoE']);
    const primitive=conservativePrimitive(q,r.gamma);conservativePrimitive(old,r.gamma);
    const expected={...primitive,rhoE:q[3],temperature:primitive.p/(primitive.rho*r.gasConstant),mach:Math.hypot(primitive.u,primitive.v)/Math.sqrt(r.gamma*primitive.p/primitive.rho)};
    for(const [key,value] of Object.entries(expected))requireValue(near(finite(field[key],key),value)&&near(Number(row[key]),value),'显示场/CSV/理想气体状态不一致。');
    requireValue(q.every((v,k)=>v===checkpoint.cells[i][k]),'显示场与检查点不一致。');
    current.push(q);previous.push(old);
    const heatRate=thermal?numbers(row,['heatRate'])[0]:0;requireValue(heatRate>=0&&(r.thermalConductivity>0||heatRate===0),'导热步长系数错误。');heatRates.push(heatRate);
  }
  const faceRows=rows(facesText,['face','owner','neighbour','x','y','sx','sy','waveSpeed','mass','momentumX','momentumY','energy']);requireValue(faceRows.length===s.faces,'面 CSV 数量不匹配。');
  const residual=mesh.cells.map(()=>[0,0,0,0]),absolute=mesh.cells.map(()=>[0,0,0,0]),spectral=mesh.cells.map(()=>0);
  let lastFallback=0,boundaryHeat=0;
  faceRows.forEach((row,i)=>{
    const mask=count(Number(row.hllcFallbackStages??(modern?NaN:0)),'通量回退');
    requireValue(mask<=(r.order===1?1:3)&&(r.fluxScheme==='hllc'||mask===0),'通量回退阶段错误。');
    if(Number(row.partner)<0||i<Number(row.partner))lastFallback+=(mask&1)+((mask>>1)&1);
    const e=mesh.edges[i];requireValue(Number(row.face)===i&&Number(row.owner)===e.owner&&Number(row.neighbour)===e.neighbour,'面关联不匹配。');
    const vertices=mesh.cells[e.owner].vertices,local=vertices.indexOf(e.a);requireValue(local>=0,'owner 面端点缺失。');
    const forward=vertices[(local+1)%vertices.length]===e.b;
    requireValue(forward||vertices[(local+vertices.length-1)%vertices.length]===e.b,'owner 面不在多边形上。');
    const a=mesh.vertices[forward?e.a:e.b],b=mesh.vertices[forward?e.b:e.a],sx=b[1]-a[1],sy=a[0]-b[0];
    for(const [key,value] of Object.entries({x:(a[0]+b[0])/2,y:(a[1]+b[1])/2,sx,sy}))requireValue(near(Number(row[key]),value),'面几何不匹配。');
    const wave=numbers(row,['waveSpeed'])[0];requireValue(wave>0,'非正声学波速。');const flux=numbers(row,['mass','momentumX','momentumY','energy']);
    if(thermal) {
      const [heat,convective]=numbers(row,['heatFlux','convectiveEnergy']);requireValue(near(heat+convective,flux[3])&&(r.thermalConductivity>0||heat===0),'总能量/导热通量分解错误。');
      if(e.neighbour<0&&Number(row.partner)<0)boundaryHeat+=heat;
      if(row.kind==='slip-wall')requireValue(flux[0]===0&&near(flux[3],heat),'热壁质量或能量通量错误。');
    }
    for(const [id,sign] of [[e.owner,1],...(e.neighbour>=0?[[e.neighbour,-1]]:[])]) {
      spectral[id]+=wave*Math.hypot(sx,sy);
      for(let k=0;k<4;k++){residual[id][k]+=sign*flux[k];absolute[id][k]+=Math.abs(flux[k]);}
    }
  });
  if(modern)requireValue(lastFallback===s.lastHllcFallbackEvaluations,'最后步通量回退计数不匹配。');
  let maximumCellBalanceRelative=0,acousticCourant=0,thermalCourant=0,combinedCourant=0;
  for(let i=0;i<current.length;i++) {
    acousticCourant=Math.max(acousticCourant,s.lastStep*spectral[i]/areas[i]);
    thermalCourant=Math.max(thermalCourant,s.lastStep*heatRates[i]);combinedCourant=Math.max(combinedCourant,s.lastStep*(spectral[i]/areas[i]+heatRates[i]));
    for(let k=0;k<4;k++) {
      const error=areas[i]*(current[i][k]-previous[i][k])+s.lastStep*residual[i][k];
      const scale=areas[i]*(Math.abs(current[i][k])+Math.abs(previous[i][k]))+s.lastStep*absolute[i][k];
      maximumCellBalanceRelative=Math.max(maximumCellBalanceRelative,Math.abs(error)/(scale||1));
    }
  }
  requireValue(maximumCellBalanceRelative<1e-12&&acousticCourant<=r.cfl*(1+1e-12),'最后时间步的守恒或声学 CFL 不通过。');
  requireValue(combinedCourant<=r.cfl*(1+1e-12),'声学/导热组合 CFL 超限。');
  if(thermal)requireValue(near(s.thermalCourant,thermalCourant)&&near(s.combinedCourant,combinedCourant)&&near(s.boundaryHeat,boundaryHeat),'导热步长或边界热量摘要错误。');
  const history=rows(historyText,['step','time','dt','acousticCourant','minimumDensity','minimumPressure','mass','totalEnergy']).map(row=>Object.fromEntries(Object.keys(row).map(k=>[k,finite(Number(row[k]),k)])));
  requireValue(history.length===s.acceptedSteps,'时间历史不完整。');let time=startTime,step=s.steps-s.acceptedSteps,previousIntegral=null,historyMaximumEnergyBalanceRelative=0;
  for(const row of history) {
    requireValue(row.step===++step&&row.time>time&&near(row.time-time,row.dt)&&row.dt>0&&row.acousticCourant<=r.cfl*(1+1e-12)&&row.minimumDensity>0&&row.minimumPressure>0,'时间历史不满足接受条件。');
    if(thermal&&r.case==='sealed')requireValue(row.boundaryMass===0&&near(row.boundaryEnergy,row.boundaryHeat),'封闭历史的总能量通量与壁面热流不匹配。');
    if(previousIntegral)for(const [key,flux] of [['mass','boundaryMass'],['totalEnergy','boundaryEnergy']]) {
      const transfer=row.dt*row[flux],before=previousIntegral[key],after=row[key];
      const relative=Math.abs(after-before+transfer)/(Math.abs(before)+Math.abs(after)+Math.abs(transfer));
      requireValue(relative<1e-12,'相邻历史的质量/总能量收支错误。');
      if(key==='totalEnergy')historyMaximumEnergyBalanceRelative=Math.max(historyMaximumEnergyBalanceRelative,relative);
    }
    previousIntegral=row;
    if(thermal)requireValue(row.thermalCourant>=0&&row.combinedCourant>=row.thermalCourant&&row.combinedCourant<=r.cfl*(1+1e-12),'历史组合 CFL 无效。');time=row.time;
  }
  requireValue(time===s.time,'历史终点不匹配。');
  for(const [k,key] of ['mass','momentumX','momentumY','totalEnergy'].entries())requireValue(near(history.at(-1)[key],current.reduce((sum,q,i)=>sum+areas[i]*q[k],0)),'历史终态积分不匹配。');
  const minDensity=current.reduce((value,q)=>Math.min(value,q[0]),Infinity);
  const minPressure=current.reduce((value,q)=>Math.min(value,conservativePrimitive(q,r.gamma).p),Infinity);
  requireValue(near(history.at(-1).minimumDensity,minDensity)&&near(history.at(-1).minimumPressure,minPressure)&&near(history.at(-1).acousticCourant,acousticCourant),'历史终态范围或 CFL 不匹配。');
  if(thermal)requireValue(near(history.at(-1).thermalCourant,thermalCourant)&&near(history.at(-1).combinedCourant,combinedCourant)&&near(history.at(-1).boundaryHeat,boundaryHeat),'历史终态热通量或组合 CFL 不匹配。');
  if(modern) {
    for(const row of history)requireValue(row.minimumContactRestoration>=0&&row.minimumContactRestoration<=1,'历史激波保护权重无效。');
    requireValue(history.reduce((value,row)=>Math.min(value,row.minimumContactRestoration),1)===s.minimumContactRestoration&&history.at(-1).minimumContactRestoration===s.lastMinimumContactRestoration,'历史激波保护权重不匹配。');
  }
  if(modern)for(const [key,lastKey,limit] of [['hllcFallbackEvaluations','lastHllcFallbackEvaluations',s.faces*r.order],['reconstructionFallbackCells','lastReconstructionFallbackCells',s.cells*r.order]]) {
    const total=history.reduce((sum,row)=>{const value=count(row[key],key);requireValue(value<=limit,'回退计数超出网格范围。');return sum+value;},0);
    requireValue(total===s[key]&&history.at(-1)[key]===s[lastKey],'历史回退计数不匹配。');
  }
  return {summary:s,fields,history,request:r,audit:{maximumCellBalanceRelative,historyMaximumEnergyBalanceRelative,acousticCourant,thermalCourant,combinedCourant,boundaryHeat,scope:'last accepted step finite-volume balance, EOS, mesh/field/checkpoint binding; full selected numerical flux and RK stages audited separately'}};
}
module.exports={SUFFIXES,PHYSICAL,validateEulerRequest,buildEulerInvocation,eulerCheckpoint,parseEulerProgress,validateEulerOutput};
