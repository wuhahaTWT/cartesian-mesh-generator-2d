'use strict';
const fs=require('node:fs/promises'),path=require('node:path'),{createHash}=require('node:crypto');
const {readFlowHistory}=require('./flow');
const checkpointDigest=text=>createHash('sha256').update(text).digest('hex');
function historyForRestart(current) {
  const saved=current?.flowRestart;if(!saved)return [];
  if(saved.history)return saved.history;
  const complete=current.flow;
  if(complete?.files.checkpoint && path.resolve(current.outputDirectory,complete.files.checkpoint)===path.resolve(saved.path)) {
    const rows=complete.history;
    if(rows?.length)return [{file:path.resolve(current.outputDirectory,complete.files.timeHistory),
      startTime:rows[0].time-rows[0].dt,endTime:saved.metadata.time,request:complete.request,complete:true}];
  }
  return [];
}
async function readFlowTimeline(current,readText=file=>fs.readFile(file,'utf8')) {
  const restart=current?.flowRestart,segments=historyForRestart(current),rows=[],gaps=[];
  const root=path.resolve(current.outputDirectory);
  let lastEnd=null,previousSample=null,firstStart=null,uncommittedSamples=0;
  for(const segment of segments) {
    const file=path.resolve(segment.file),relative=path.relative(root,file);
    if(relative==='..'||relative.startsWith('..'+path.sep)||path.isAbsolute(relative))throw new Error('流动历史路径越出项目目录。');
    const {startTime,endTime}=segment;
    if(!Number.isFinite(startTime)||!Number.isFinite(endTime)||startTime<0||endTime<=startTime||
      (lastEnd!==null&&startTime<lastEnd))throw new Error('流动历史区间重叠或时钟无效。');
    if(firstStart===null)firstStart=startTime;
    const parsed=readFlowHistory(await readText(file),segment.request,startTime,endTime,{partial:!segment.complete});
    uncommittedSamples+=parsed.uncommittedSamples;
    if(lastEnd!==null&&startTime>lastEnd)gaps.push({from:lastEnd,to:startTime,reason:'未保存此区间的运行历史'});
    if(previousSample!==null&&previousSample<startTime&&parsed.rows.length)parsed.rows[0].breakBefore=true;
    rows.push(...parsed.rows);
    const last=parsed.rows.at(-1)?.time ?? startTime;
    if(last<endTime)gaps.push({from:last,to:endTime,reason:'检查点已接受，历史尚未完整写出'});
    previousSample=parsed.rows.at(-1)?.time ?? previousSample;lastEnd=endTime;
  }
  if(restart&&lastEnd!==null&&lastEnd!==restart.metadata.time)throw new Error('流动历史终点与检查点不同。');
  return {rows,segments:segments.length,availableSteps:rows.length,gaps,uncommittedSamples,
    startTime:firstStart ?? restart?.metadata.time ?? null,endTime:restart?.metadata.time ?? null};
}
module.exports={checkpointDigest,historyForRestart,readFlowTimeline};
