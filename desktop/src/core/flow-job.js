'use strict';
const fs=require('node:fs/promises'),path=require('node:path');
const {FLOW_CASES,buildFlowInvocation,validateFlowRequest,parseFlowProgress,validateFlowOutput,validateTimeHistory,validateAttemptHistory}=require('./flow');
const {readCheckpointMetadata,readAcceptedCheckpointMetadata}=require('./flow-checkpoint');
const {validateBoundaryMesh,serializeBoundaryDefinition,parseBoundaryDefinition,sameConditions}=require('./flow-boundaries');
const {checkpointDigest,historyForRestart}=require('./flow-history');
const readJson=async file=>JSON.parse(await fs.readFile(file,'utf8'));

// Each call owns its output directory. Only a completely read and validated
// result replaces the displayed binding; accepted restart state is independent.
async function runFlowJob({currentResult,mesh,request,executable,runProcess,signal,onProgress,log,onPrepared}) {
    request=validateFlowRequest(request);
    const selectedRestart = request?.resume ? currentResult.flowRestart : null;
    // Validate before creating outputs; paths come only from this main process.
    buildFlowInvocation(currentResult.cm2dPath, 'pending', request, selectedRestart?.path, 'pending.boundaries');
    const boundaryDefinition=request.case==='custom' ? validateBoundaryMesh(request.boundaryDefinition,mesh,Number(request.speed)) : null;
    const incompleteDirectory = await fs.mkdtemp(path.join(currentResult.outputDirectory, 'flow-run-'));
    const pendingPrefix = path.join(incompleteDirectory, 'flow');
    let restartPath = null, startTime = 0;
    if (selectedRestart) {
      try {
        restartPath = path.join(incompleteDirectory, 'input.checkpoint');
        await fs.copyFile(selectedRestart.path, restartPath);
        const metadata = await readAcceptedCheckpointMetadata(restartPath);
        if(metadata.time!==selectedRestart.metadata.time ||
          (selectedRestart.sha256 && selectedRestart.sha256!==checkpointDigest(await fs.readFile(restartPath))))
          throw new Error('选入流动状态已改变，不能连接此前历史。');
        startTime = metadata.time;
        for (const key of ['case', 'nu', 'speed', 'convection', 'outletBackflow'])
          if (metadata[key] !== (['nu','speed'].includes(key) ? Number(request[key]) : request[key]))
            throw new Error('续算必须保持原工况、物性和对流格式；可调整时间步与步数。');
        if (boundaryDefinition && !sameConditions(metadata.boundaryDefinition?.records,boundaryDefinition.records))
          throw new Error('续算必须保持原命名边界的名称、类型和数值。');
      } catch (error) {
        // No solver has started and the selected source remains untouched.
        await fs.rm(incompleteDirectory, { recursive: true, force: true });
        throw error;
      }
    }
    const boundaryPath=boundaryDefinition ? path.join(incompleteDirectory,'input.boundaries') : null;
    if (boundaryPath) await fs.writeFile(boundaryPath,serializeBoundaryDefinition(boundaryDefinition));
    const invocation = buildFlowInvocation(currentResult.cm2dPath, pendingPrefix, request, restartPath, boundaryPath);
    const transient = invocation.request.mode !== 'steady';
    const adaptive = invocation.request.mode === 'adaptive';
    if (adaptive && !(invocation.request.endTime>startTime)) {
      await fs.rm(incompleteDirectory,{recursive:true,force:true});
      throw new Error('目标物理时间必须晚于已接受的重启时间。');
    }
    const previousRestart = currentResult.flowRestart;
    const previousHistory=invocation.request.resume?historyForRestart(currentResult):[];
    const readRestart=async(complete=false)=>{
      const file=pendingPrefix+'.checkpoint',metadata=await readAcceptedCheckpointMetadata(file);
      if(metadata.time<=startTime)throw new Error('本段尚未接受新的时间步。');
      return {path:file,metadata,sha256:checkpointDigest(await fs.readFile(file)),history:[...previousHistory,
        {file:pendingPrefix+'.time-history.csv',startTime,endTime:metadata.time,request:invocation.request,complete}]};
    };
    const preserveIncomplete = async error => {
      if (transient) {
        // Ignore .tmp: only the native atomic accepted-state file is resumable.
        try { currentResult.flowRestart = await readRestart(); }
        catch { currentResult.flowRestart = previousRestart; }
      }
      const report = { format: 'cartmesh2d-flow-incomplete-v1',
        status: signal.aborted ? 'cancelled' : 'failed', request: invocation.request,
        mesh: path.basename(currentResult.cm2dPath),
        acceptedTime: currentResult.flowRestart?.metadata.time ?? null,
        exitCode: Number.isInteger(error.code) ? error.code : null,
        message: String(error.message || error).split('\n')[0] };
      await fs.writeFile(path.join(incompleteDirectory, 'desktop-flow-error.json'), JSON.stringify(report, null, 2));
    };
    log(`正在运行原生二维${transient ? '非定常' : '稳态'}层流：${FLOW_CASES[invocation.request.case].label}…`);
    const onLine = (line, isError) => {
      let progress = null;
      if (!isError) {
        try { progress = parseFlowProgress(line); }
        catch (error) { log(`忽略无效进度：${error.message}`); }
      }
      if (progress) onProgress(progress);
      else log(line);
    };
    try {
      await fs.writeFile(path.join(incompleteDirectory,'desktop-state.json'),JSON.stringify({status:'running',request:invocation.request,startTime}));
      await onPrepared?.({directory:path.relative(currentResult.outputDirectory,incompleteDirectory).split(path.sep).join('/'),startTime,request:invocation.request});
      const processResult = await runProcess(executable(invocation.executable), invocation.args,
        onLine, signal, 0, [0, 2],{onSpawn:pid=>{
          const sync=require('node:fs'),file=path.join(incompleteDirectory,'native-process.json');
          sync.writeFileSync(file+'.tmp',JSON.stringify({pid}));sync.renameSync(file+'.tmp',file);
        }});
      signal.throwIfAborted();
      const outputFiles = { summary: `${pendingPrefix}.json`, fields: `${pendingPrefix}.fields.json`,
        vtk: `${pendingPrefix}.vtk`, residuals: `${pendingPrefix}.residuals.csv`,
        cells: `${pendingPrefix}.cells.csv`, faces: `${pendingPrefix}.faces.csv` };
      if (transient) Object.assign(outputFiles, { checkpoint: `${pendingPrefix}.checkpoint`, timeHistory: `${pendingPrefix}.time-history.csv` });
      if (adaptive) outputFiles.attemptHistory=`${pendingPrefix}.attempt-history.csv`;
      if (invocation.request.initialVortex) outputFiles.initialCheckpoint=`${pendingPrefix}.initial.checkpoint`;
      if (boundaryDefinition) outputFiles.boundaries=`${pendingPrefix}.boundaries`;
      const [summary, fields] = await Promise.all([readJson(outputFiles.summary), readJson(outputFiles.fields),
        ...Object.values(outputFiles).map(file => fs.stat(file))]);
      const validated = validateFlowOutput(summary, fields, mesh.cells.length, invocation.request, startTime);
      if (boundaryDefinition) {
        const exported=validateBoundaryMesh(parseBoundaryDefinition(await fs.readFile(outputFiles.boundaries,'utf8')),mesh,invocation.request.speed);
        if (!sameConditions(exported.records,boundaryDefinition.records)) throw new Error('导出边界与输入不一致。');
      }
      if ((processResult.code === 0) !== validated.summary.converged)
        throw new Error('原生求解器退出码与收敛状态不一致。');
      let history = null, attempts = null, checkpointMetadata = null;
      if (transient) {
        history = validateTimeHistory(await fs.readFile(outputFiles.timeHistory, 'utf8'), validated.summary, startTime);
        if (adaptive) attempts=validateAttemptHistory(await fs.readFile(outputFiles.attemptHistory,'utf8'),validated.summary,history);
        checkpointMetadata = await readCheckpointMetadata(outputFiles.checkpoint);
        if (outputFiles.initialCheckpoint) {
          const initial=await readCheckpointMetadata(outputFiles.initialCheckpoint);
          if (initial.time!==0 || ['case','nu','speed','convection','outletBackflow'].some(k=>initial[k]!==checkpointMetadata[k]))
            throw new Error('初始局部涡检查点的时间或物性与结果不一致。');
        }
        if (Math.abs(checkpointMetadata.time-summary.acceptedTime) > 1e-12+1e-9*Math.abs(summary.acceptedTime))
          throw new Error('重启状态时间与摘要不一致。');
        if (!validated.summary.converged)
          throw Object.assign(new Error(summary.acceptedTime>0 ? `时间步未收敛；已接受到 t=${summary.acceptedTime} s，可继续计算。候选场仅留作诊断。` : '本次没有接受物理时间步；初值不能作为续算状态。候选场仅留作诊断。'), { code: 2 });
      }
      signal.throwIfAborted();
      const saved=Object.fromEntries(Object.entries(outputFiles).map(([kind,file])=>[kind,path.relative(currentResult.outputDirectory,file).split(path.sep).join('/')]));
      const payload={...validated,request:invocation.request,files:saved,history,attempts};
      const restartState=transient?await readRestart(true):null;
      await fs.writeFile(path.join(incompleteDirectory,'desktop-state.json.tmp'),JSON.stringify({status:'complete',request:invocation.request,startTime,time:summary.acceptedTime ?? null},null,2));
      await fs.rename(path.join(incompleteDirectory,'desktop-state.json.tmp'),path.join(incompleteDirectory,'desktop-state.json'));
      signal.throwIfAborted();
      currentResult.flow=payload;currentResult.flowRestart=restartState;
      log(transient ? `非定常计算完成，已接受到 t=${summary.acceptedTime} s。`
        : validated.summary.converged ? '层流求解已收敛。' : '层流求解到达迭代上限，保留诊断结果但未收敛。');
      return payload;
    } catch (error) {
      await preserveIncomplete(error).catch(() => {});
      error.message += `\n未完成诊断保留在 ${incompleteDirectory}`;
      throw error;
    }
}
async function checkFlowRestart({currentResult,checkpoint,request,executable,runProcess,signal}) {
  const r=validateFlowRequest({...request,resume:true,initialVortex:undefined});
  const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'flow-check-'));
  let boundary=null;
  if(r.case==='custom') {
    boundary=path.join(directory,'input.boundaries');
    const definition=validateBoundaryMesh(r.boundaryDefinition,currentResult.mesh,r.speed);
    await fs.writeFile(boundary,serializeBoundaryDefinition(definition));
  }
  const invocation=buildFlowInvocation(currentResult.cm2dPath,path.join(directory,'check'),r,checkpoint,boundary);
  const result=await runProcess(executable(invocation.executable),[...invocation.args,'--check-restart','on'],()=>{},signal);
  const meta=JSON.parse(result.stdout);
  if(meta.format!=='cartmesh2d-flow-restart-check-v1'||meta.status!=='valid'||meta.cells!==currentResult.mesh.cells.length)
    throw new Error('原生流动检查点核对未通过。');
  return meta;
}
module.exports={runFlowJob,checkFlowRestart};
