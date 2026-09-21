'use strict';
const fs=require('node:fs/promises');
const path=require('node:path');
const {createHash}=require('node:crypto');
const {SUFFIXES,PHYSICAL,validateEulerRequest,buildEulerInvocation,eulerCheckpoint,parseEulerProgress,validateEulerOutput}=require('./euler');
const hash=data=>createHash('sha256').update(data).digest('hex');
const manifestFormat='cartmesh2d-euler-run-v1';
async function readCheckpoint(file,mesh,request) {
  const info=await fs.stat(file);if(!info.isFile()||info.size>256*1024**2)throw new Error('Euler 检查点过大或不是文件。');
  const text=await fs.readFile(file,'utf8'),state=eulerCheckpoint(text,mesh,request);
  return {path:file,metadata:{time:state.time,steps:state.steps,request,sha256:hash(text)}};
}
async function importEulerRestart(file,mesh,meshPath) {
  const stat=await fs.stat(file);if(stat.size>65536)throw new Error('Euler 续算清单过大。');
  const record=JSON.parse(await fs.readFile(file,'utf8'));
  if(record?.format!==manifestFormat||!['complete','failed','cancelled'].includes(record.status)||record.checkpoint!=='euler.checkpoint')throw new Error('请选择 Euler 结果目录中的 desktop-state.json。');
  const request=validateEulerRequest(record.request);
  if(record.meshSha256!==hash(await fs.readFile(meshPath)))throw new Error('Euler 续算需要同一最终网格。');
  const state=await readCheckpoint(path.join(path.dirname(file),'euler.checkpoint'),mesh,request);
  if(state.metadata.sha256!==record.checkpointSha256||state.metadata.time!==record.acceptedTime)throw new Error('Euler 检查点与清单不一致。');
  return state;
}
async function runEulerJob({currentResult,mesh,request,executable,runProcess,signal,onProgress=()=>{},log=()=>{}}) {
  const normalized=validateEulerRequest(request),selected=normalized.resume?currentResult.eulerRestart:null;
  if(normalized.resume&&!selected)throw new Error('没有可用的可压续算状态。');
  if(selected) {
    for(const key of PHYSICAL)if(normalized[key]!==selected.metadata.request[key])throw new Error('续算须保持原工况、气体与物理参数；可改目标时间与数值控制。');
    if(normalized.endTime<=selected.metadata.time)throw new Error('目标时间须晚于已接受时间。');
  }
  const meshSha256=hash(await fs.readFile(currentResult.cm2dPath));
  const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'euler-run-')),prefix=path.join(directory,'euler');
  const manifest=path.join(directory,'desktop-state.json'),previousRestart=currentResult.eulerRestart;
  const record={format:manifestFormat,status:'running',request:normalized,meshSha256,checkpoint:'euler.checkpoint'};
  const writeRecord=async(data)=>{await fs.writeFile(manifest+'.tmp',JSON.stringify({...record,...data},null,2));await fs.rename(manifest+'.tmp',manifest);};
  await writeRecord({});let startTime=0;
  try {
    let restart=null;
    if(selected) {
      const source=await readCheckpoint(selected.path,mesh,normalized);
      if(source.metadata.sha256!==selected.metadata.sha256)throw new Error('已选择的 Euler 检查点发生变化。');
      startTime=source.metadata.time;restart=path.join(directory,'input.checkpoint');await fs.copyFile(selected.path,restart);
    }
    const invocation=buildEulerInvocation(currentResult.cm2dPath,prefix,normalized,restart);
    const result=await runProcess(executable(invocation.executable),invocation.args,(line,isError)=>{
      let progress=null;
      if(!isError)try{progress=parseEulerProgress(line);}catch(e){log(e.message);}
      if(progress)onProgress(progress);else log(line);
    },signal,0,[0,2]);
    signal.throwIfAborted();
    if(result.code!==0)throw new Error(result.stderr?.trim()||'Euler 未到达目标时间。');
    await Promise.all(SUFFIXES.map(suffix=>fs.stat(prefix+suffix)));
    const files=await Promise.all(['.json','.fields.json','.cells.csv','.faces.csv','.history.csv','.checkpoint'].map(suffix=>fs.readFile(prefix+suffix,'utf8')));
    const validated=validateEulerOutput(JSON.parse(files[0]),JSON.parse(files[1]),...files.slice(2),mesh,normalized,startTime);
    if(hash(await fs.readFile(currentResult.cm2dPath))!==meshSha256)throw new Error('计算期间网格发生变化。');
    const nextRestart=await readCheckpoint(prefix+'.checkpoint',mesh,normalized);
    signal.throwIfAborted();
    const payload={...validated,files:Object.fromEntries(SUFFIXES.map(suffix=>[suffix,path.relative(currentResult.outputDirectory,prefix+suffix)])),manifest:path.relative(currentResult.outputDirectory,manifest)};
    await writeRecord({status:'complete',acceptedTime:nextRestart.metadata.time,checkpointSha256:nextRestart.metadata.sha256});
    signal.throwIfAborted();currentResult.euler=payload;currentResult.eulerRestart=nextRestart;return payload;
  }catch(error) {
    try{currentResult.eulerRestart=await readCheckpoint(prefix+'.checkpoint',mesh,normalized);}catch{currentResult.eulerRestart=previousRestart;}
    const next=currentResult.eulerRestart?.path===prefix+'.checkpoint'?currentResult.eulerRestart:null;
    await writeRecord({status:signal.aborted?'cancelled':'failed',acceptedTime:next?.metadata.time??null,checkpointSha256:next?.metadata.sha256??null,error:String(error.message)}).catch(()=>{});
    error.message+=`\n诊断与已接受状态保留在 ${directory}`;throw error;
  }
}
module.exports={runEulerJob,readCheckpoint,importEulerRestart};
