'use strict';
const fs=require('node:fs/promises');
const path=require('node:path');
const {SUFFIXES,validateThermalRequest,thermalBoundaryCsv,buildThermalInvocation,parseThermalProgress,thermalCheckpointTime,validateThermalOutput}=require('./thermal');

// A run owns a new directory. Publish the in-memory binding only after every
// output is checked; the previous complete result never participates in writes.
async function runThermalJob({currentResult,mesh,request,executable,runProcess,signal,onProgress,log}) {
  const normalized=validateThermalRequest(request);
  const selected=normalized.resume?currentResult.thermalRestart:null;
  if(normalized.resume&&!selected)throw new Error('没有可用的联合续算状态。');
  const boundary=thermalBoundaryCsv(mesh,normalized);
  if(selected) {
    for(const key of ['case','nu','speed','convection','outletBackflow','diffusivity','source','scalarConvection','boundaries'])
      if(JSON.stringify(normalized[key])!==JSON.stringify(selected.metadata.request[key]))
        throw new Error('联合续算必须保持工况、物性、温度源、边界和对流格式。');
  }
  const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'thermal-run-'));
  const prefix=path.join(directory,'thermal'), boundaryPath=path.join(directory,'boundary.csv');
  const previousRestart=currentResult.thermalRestart;
  let startTime=0;
  const readRestart=async()=>({path:`${prefix}.thermal.checkpoint`,metadata:{
    time:thermalCheckpointTime(await fs.readFile(`${prefix}.thermal.checkpoint`,'utf8')),request:normalized}});
  await fs.writeFile(path.join(directory,'desktop-state.json'),JSON.stringify({status:'running',request:normalized}));
  try {
    await fs.writeFile(boundaryPath,boundary);
    let restart=null;
    if(selected) {
      restart=path.join(directory,'input.thermal.checkpoint');
      await fs.copyFile(selected.path,restart);
      startTime=thermalCheckpointTime(await fs.readFile(restart,'utf8'));
      if(startTime!==selected.metadata.time)throw new Error('续算文件已发生变化。');
    }
    const invocation=buildThermalInvocation(currentResult.cm2dPath,prefix,boundaryPath,normalized,restart);
    const processResult=await runProcess(executable(invocation.executable),invocation.args,(line,isError)=>{
      let progress=null;
      if(!isError)try{progress=parseThermalProgress(line);}catch(e){log(`忽略无效热进度：${e.message}`);}
      if(progress)onProgress(progress);else log(line);
    },signal,0,[0,2]);
    signal.throwIfAborted();
    const summary=JSON.parse(await fs.readFile(`${prefix}.json`,'utf8'));
    if(processResult.code!==0||summary.converged!==true)throw new Error(`热计算未完成；已接受到 t=${summary.acceptedTime ?? startTime} s。`);
    await Promise.all(SUFFIXES.map(suffix=>fs.stat(prefix+suffix)));
    const [cells,history,joint]=await Promise.all(['.cells.csv','.thermal-history.csv','.thermal.checkpoint'].map(suffix=>fs.readFile(prefix+suffix,'utf8')));
    const validated=validateThermalOutput(summary,cells,history,joint,mesh,normalized,startTime);
    const restartState=await readRestart();
    signal.throwIfAborted();
    const files=Object.fromEntries(SUFFIXES.map(suffix=>[suffix,path.relative(currentResult.outputDirectory,prefix+suffix)]));
    const payload={...validated,files};
    await fs.writeFile(path.join(directory,'desktop-state.json.tmp'),JSON.stringify({status:'complete',request:normalized,time:summary.time},null,2));
    await fs.rename(path.join(directory,'desktop-state.json.tmp'),path.join(directory,'desktop-state.json'));
    signal.throwIfAborted();
    currentResult.thermal=payload;
    currentResult.thermalRestart=restartState;
    return payload;
  } catch(error) {
    try{currentResult.thermalRestart=await readRestart();}catch{currentResult.thermalRestart=previousRestart;}
    await fs.writeFile(path.join(directory,'desktop-state.json'),JSON.stringify({status:signal.aborted?'cancelled':'failed',
      request:normalized,acceptedTime:currentResult.thermalRestart?.metadata.time ?? null,error:String(error.message)},null,2)).catch(()=>{});
    if(/outlet backflow unsupported/.test(error.message))error.message='当前设置遇到出口回流会停止。请检查计算域，或在流动设置中选择试验性法向回流。\n'+error.message;
    error.message+=`\n诊断及最后联合状态保留在 ${directory}`;
    throw error;
  }
}
module.exports={runThermalJob};
