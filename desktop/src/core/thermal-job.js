'use strict';
const fs=require('node:fs/promises');
const path=require('node:path');
const {checkpointDigest,historyForRestart}=require('./thermal-history');
const {SUFFIXES,validateThermalRequest,thermalBoundaryCsv,thermalEventsCsv,buildThermalInvocation,parseThermalProgress,thermalCheckpointTime,validateThermalOutput,thermalFailureMessage}=require('./thermal');

// A run owns a new directory. Publish the in-memory binding only after every
// output is checked; the previous complete result never participates in writes.
async function runThermalJob({currentResult,mesh,request,executable,runProcess,signal,onProgress,log,onPrepared}) {
  const normalized=validateThermalRequest(request);
  const selected=normalized.resume?currentResult.thermalRestart:null;
  if(normalized.resume&&!selected)throw new Error('没有可用的联合续算状态。');
  const boundary=thermalBoundaryCsv(mesh,normalized);
  const events=normalized.events.length?thermalEventsCsv(mesh,normalized):null;
  if(selected) {
    if(JSON.stringify(normalized.events)!==JSON.stringify(validateThermalRequest(selected.metadata.request).events))
      throw new Error('联合续算必须保持完整热源和边界时间规律。');
    if(normalized.fluxCorrection!==(selected.metadata.request.fluxCorrection ?? 'unrestricted'))
      throw new Error('联合续算必须保持温度通量修正。');
    for(const key of ['case','nu','speed','convection','outletBackflow','diffusivity','source','scalarConvection','boundaries'])
      if(JSON.stringify(normalized[key])!==JSON.stringify(selected.metadata.request[key]))
        throw new Error('联合续算必须保持工况、物性、温度源、边界和对流格式。');
  }
  const suffixes=[...SUFFIXES,...(normalized.initialVortex?['.initial.checkpoint']:[]),...(normalized.mode==='adaptive'?['.attempt-history.csv']:[]),...(normalized.events.length?['.events.csv']:[])];
  const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'thermal-run-'));
  const prefix=path.join(directory,'thermal'), boundaryPath=path.join(directory,'boundary.csv');
  const previousRestart=currentResult.thermalRestart;
  const previousHistory=normalized.resume?historyForRestart(currentResult):[];
  let startTime=0;
  const readRestart=async(complete=false)=>{
    const text=await fs.readFile(`${prefix}.thermal.checkpoint`,'utf8');
    const time=thermalCheckpointTime(text,normalized.events);
    if(time<=startTime)throw new Error('本段尚未接受新的时间步。');
    return {path:`${prefix}.thermal.checkpoint`,metadata:{
    time,request:normalized},sha256:checkpointDigest(text),history:[...previousHistory,{file:`${prefix}.thermal-history.csv`,startTime,endTime:time,request:normalized,complete}]};
  };
  await fs.writeFile(path.join(directory,'desktop-state.json'),JSON.stringify({status:'running',request:normalized}));
  try {
    await fs.writeFile(boundaryPath,boundary);
    if(events)await fs.writeFile(prefix+'.events.csv',events);
    let restart=null;
    if(selected) {
      restart=path.join(directory,'input.thermal.checkpoint');
      await fs.copyFile(selected.path,restart);
      const inputText=await fs.readFile(restart,'utf8');
      if(selected.sha256 && selected.sha256!==checkpointDigest(inputText))throw new Error('续算文件与已保存的历史链不同。');
      startTime=thermalCheckpointTime(inputText,normalized.events);
      if(startTime!==selected.metadata.time)throw new Error('续算文件已发生变化。');
    }
    await onPrepared?.({directory:path.relative(currentResult.outputDirectory,directory).split(path.sep).join('/'),startTime,request:normalized});
    const invocation=buildThermalInvocation(currentResult.cm2dPath,prefix,boundaryPath,normalized,restart,events?prefix+'.events.csv':null);
    const processResult=await runProcess(executable(invocation.executable),invocation.args,(line,isError)=>{
      let progress=null;
      if(!isError)try{progress=parseThermalProgress(line);}catch(e){log(`忽略无效热进度：${e.message}`);}
      if(progress)onProgress(progress);else log(line);
    },signal,0,[0,2],{onSpawn:pid=>{
      const sync=require('node:fs'),file=path.join(directory,'native-process.json');
      sync.writeFileSync(file+'.tmp',JSON.stringify({pid}));sync.renameSync(file+'.tmp',file);
    }});
    signal.throwIfAborted();
    const summary=JSON.parse(await fs.readFile(`${prefix}.json`,'utf8'));
    if(processResult.code!==0||summary.converged!==true)throw new Error(thermalFailureMessage(summary,startTime));
    await Promise.all(suffixes.map(suffix=>fs.stat(prefix+suffix)));
    const [cells,history,joint]=await Promise.all(['.cells.csv','.thermal-history.csv','.thermal.checkpoint'].map(suffix=>fs.readFile(prefix+suffix,'utf8')));
    const validated=validateThermalOutput(summary,cells,history,joint,mesh,normalized,startTime);
    const restartState=await readRestart(true);
    signal.throwIfAborted();
    const files=Object.fromEntries(suffixes.map(suffix=>[suffix,path.relative(currentResult.outputDirectory,prefix+suffix)]));
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
    const retained=currentResult.thermalRestart?.metadata.time;
    error.message+=retained>0?`\n已保留 t=${retained} s 的联合续算状态。`:'\n尚无已接受的联合续算状态。';
    error.message+=`\n诊断文件保留在 ${directory}`;
    throw error;
  }
}
async function checkThermalRestart({currentResult,checkpoint,request,executable,runProcess,signal}) {
  // The checkpoint already contains the carrier; original startup metadata
  // is retained in the project but must not reapply a fresh vortex on validation.
  const r=validateThermalRequest({...request,initialVortex:undefined,resume:true});
  const mesh=currentResult.mesh;
  const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'thermal-check-'));
  const boundary=path.join(directory,'boundary.csv'),events=r.events.length?path.join(directory,'events.csv'):null;
  await fs.writeFile(boundary,thermalBoundaryCsv(mesh,r));
  if(events)await fs.writeFile(events,thermalEventsCsv(mesh,r));
  const invocation=buildThermalInvocation(currentResult.cm2dPath,path.join(directory,'check'),boundary,r,checkpoint,events);
  const result=await runProcess(executable(invocation.executable),[...invocation.args,'--check-restart','on'],()=>{},signal);
  const meta=JSON.parse(result.stdout);
  if(meta.format!=='cartmesh2d-thermal-restart-check-v1'||meta.status!=='valid'||meta.cells!==mesh.cells.length)
    throw new Error('原生联合检查点核对未通过。');
  return meta;
}
module.exports={runThermalJob,checkThermalRestart};
