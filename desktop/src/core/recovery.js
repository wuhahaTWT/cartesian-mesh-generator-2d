'use strict';
const fs=require('node:fs/promises'),path=require('node:path');
const {randomUUID}=require('node:crypto');
const {writeProjectManifest,readProject}=require('./project');
const {thermalCheckpointTime,validateThermalRequest}=require('./thermal');
const {validateFlowRequest}=require('./flow');
const {readAcceptedCheckpointMetadata}=require('./flow-checkpoint');
const {historyForRestart:flowHistoryForRestart,readFlowTimeline}=require('./flow-history');
const {checkpointDigest,historyForRestart,readThermalTimeline}=require('./thermal-history');
const alive=pid=>{if(!Number.isSafeInteger(pid)||pid<=0)return false;try{process.kill(pid,0);return true;}catch(e){return e.code==='EPERM';}};
async function atomicJson(file,value) {await fs.writeFile(file+'.tmp',JSON.stringify(value,null,2)+'\n');await fs.rename(file+'.tmp',file);}
function inside(root,name) {
  if(typeof name!=='string'||path.isAbsolute(name)||/[\\\0:]/.test(name)||name.split('/').some(p=>p==='..'||!p))throw new Error('恢复记录路径无效。');
  return path.join(root,...name.split('/'));
}
class RecoveryStore {
  constructor(base,session) {this.base=path.resolve(base);this.session=path.resolve(session);this.records=new Map();}
  async checkpoint(current,ui) {
    if(!current || current.incomplete)return;
    await writeProjectManifest(current,ui);
    const root=path.resolve(current.outputDirectory);
    let record=this.records.get(root);
    if(!record)record={format:'cartmesh2d-session-recovery-v1',id:randomUUID(),root,ownerPid:process.pid};
    record={...record,label:path.basename(current.job.geometryPath),updatedAt:new Date().toISOString(),
      thermalTime:current.thermalRestart?.metadata.time ?? null,flowTime:current.flowRestart?.metadata.time ?? null,
      resumeKind:ui?.inputs?.flowResume&&!ui?.inputs?.thermalResume?'flow':'thermal',activeThermal:null,activeFlow:null};
    await atomicJson(path.join(this.session,'recovery-'+record.id+'.json'),record);this.records.set(root,record);
  }
  beginThermal(current,pending) {return this.begin(current,'activeThermal',pending);}
  beginFlow(current,pending) {return this.begin(current,'activeFlow',pending);}
  async begin(current,key,pending) {
    const root=path.resolve(current.outputDirectory),record=this.records.get(root);
    if(!record)throw new Error('计算前的项目恢复点尚未写入。');
    const next={...record,updatedAt:new Date().toISOString(),[key]:pending};
    await atomicJson(path.join(this.session,'recovery-'+record.id+'.json'),next);this.records.set(root,next);
  }
  async readRecord(id) {
    if(typeof id!=='string'||!/^session-[\w-]+\/recovery-[\da-f-]+\.json$/.test(id))throw new Error('恢复会话标识无效。');
    const record=JSON.parse(await fs.readFile(inside(this.base,id),'utf8'));
    if(record?.format!=='cartmesh2d-session-recovery-v1'||typeof record.root!=='string'||!path.isAbsolute(record.root))throw new Error('恢复会话记录无效。');
    return record;
  }
  async active(record) {
    if(record.ownerPid!==process.pid && alive(record.ownerPid))return true;
    for(const pending of [record.activeThermal,record.activeFlow].filter(Boolean)) {
      try {
        const job=inside(record.root,pending.directory);
        const processInfo=JSON.parse(await fs.readFile(path.join(job,'native-process.json'),'utf8'));
        if(alive(processInfo.pid))return true;
      } catch(e) {if(e.code!=='ENOENT')throw e;}
    }
    return false;
  }
  async list() {
    const entries=[];
    for(const session of await fs.readdir(this.base,{withFileTypes:true})) {
      if(!session.isDirectory()||!session.name.startsWith('session-'))continue;
      for(const file of await fs.readdir(path.join(this.base,session.name))) {
        if(!/^recovery-[\da-f-]+\.json$/.test(file))continue;
        const id=session.name+'/'+file;
        try {
          const record=await this.readRecord(id);
          // An idle session in this process is recoverable too (e.g. after a
          // macOS window close/reopen). Main-process exclusivity prevents writes.
          entries.push({id,label:record.label,updatedAt:record.updatedAt,thermalTime:record.thermalTime,flowTime:record.flowTime ?? null,
            kind:record.activeFlow?'flow':record.activeThermal?'thermal':record.resumeKind || 'thermal',
            interrupted:Boolean(record.activeThermal||record.activeFlow),active:await this.active(record)});
        } catch(error) {entries.push({id,label:file,active:true,error:error.message});}
      }
    }
    return entries.sort((a,b)=>String(b.updatedAt).localeCompare(String(a.updatedAt)));
  }
  async restore(id,validateNative,signal) {
    const record=await this.readRecord(id);
    if(await this.active(record))throw new Error('该会话或原生程序仍在运行，不能同时恢复。');
    const root=path.join(await fs.mkdtemp(path.join(this.session,'restored-')),'project');
    await fs.cp(record.root,root,{recursive:true,dereference:false});signal?.throwIfAborted();
    const current=await readProject(root,signal);
    const kind=record.activeFlow?'flow':record.activeThermal?'thermal':record.resumeKind || 'thermal';
    const pending=kind==='flow'?record.activeFlow:record.activeThermal;
    if(pending) {
      const directory=inside(root,pending.directory),prefix=path.join(directory,kind);
      try {
        const request=(kind==='flow'?validateFlowRequest:validateThermalRequest)(pending.request);
        const file=prefix+(kind==='flow'?'.checkpoint':'.thermal.checkpoint'),text=await fs.readFile(file,'utf8');
        const metadata=kind==='flow'?await readAcceptedCheckpointMetadata(file):{time:thermalCheckpointTime(text,request.events),request};
        const time=metadata.time;
        if(time>pending.startTime) {
          const checked=await validateNative(current,file,request,signal,kind);
          if(checked && checked.time!==time)throw new Error('原生检查点时钟与恢复记录不一致。');
          const historyOf=kind==='flow'?flowHistoryForRestart:historyForRestart,key=kind+'Restart';
          const previous=path.resolve(current[key]?.path||'')===file?historyOf(current).slice(0,-1)
            :request.resume?historyOf(current):[];
          const candidate={...current,[key]:{path:file,sha256:checkpointDigest(text),metadata,history:[...previous,
            {file:prefix+(kind==='flow'?'.time-history.csv':'.thermal-history.csv'),startTime:pending.startTime,endTime:time,request,complete:false}]}};
          await (kind==='flow'?readFlowTimeline:readThermalTimeline)(candidate);
          current[key]=candidate[key];
          current.recoveryNotice=`中断运行已找回 t=${time} s 的原生接受状态；场图保留上次完整结果。`;
        } else current.recoveryNotice='本段尚未接受新的时间步，已恢复计算前的项目。';
      } catch(error) {
        signal?.throwIfAborted();
        current.recoveryNotice='最新运行未能完整恢复，保留此前项目和全部诊断：'+error.message;
      }
    } else current.recoveryNotice='已从保留的会话恢复项目，源会话文件未改动。';
    signal?.throwIfAborted();
    if(current[kind+'Restart'])current.projectUi={...current.projectUi,inputs:{...current.projectUi?.inputs,thermalResume:kind==='thermal',flowResume:kind==='flow'}};
    await this.checkpoint(current,current.projectUi);
    return current;
  }
}
module.exports={RecoveryStore,atomicJson};
