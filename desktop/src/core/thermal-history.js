'use strict';
const fs = require('node:fs/promises');
const path = require('node:path');
const {createHash} = require('node:crypto');
const {readThermalHistory} = require('./thermal');

const checkpointDigest = text => createHash('sha256').update(text).digest('hex');
function historyForRestart(current) {
  const saved=current?.thermalRestart;
  if (!saved) return [];
  if (saved.history) return saved.history;
  // Old projects have no lineage. Retain only the segment whose native
  // checkpoint is exactly the selected input; never infer parents by filename.
  const complete=current.thermal;
  if (complete && path.resolve(current.outputDirectory,complete.files['.thermal.checkpoint'])===path.resolve(saved.path))
    return [{file:path.resolve(current.outputDirectory,complete.files['.thermal-history.csv']),
      startTime:complete.summary.startTime ?? complete.summary.time-complete.summary.steps*complete.summary.timeStep,
      endTime:saved.metadata.time,request:complete.request,complete:true}];
  return [];
}
async function readThermalTimeline(current, readText= file=>fs.readFile(file,'utf8')) {
  const restart=current?.thermalRestart,segments=historyForRestart(current),rows=[],gaps=[];
  const root=path.resolve(current.outputDirectory);
  let lastEnd=null,previousSample=null,firstStart=null;
  for (const segment of segments) {
    const file=path.resolve(segment.file),relative=path.relative(root,file);
    if (relative==='..'||relative.startsWith('..'+path.sep)||path.isAbsolute(relative))throw new Error('温度历史路径越出项目目录。');
    const {startTime,endTime}=segment;
    if (!Number.isFinite(startTime)||!Number.isFinite(endTime)||startTime<0||endTime<=startTime||
      (lastEnd!==null&&startTime<lastEnd))throw new Error('温度历史区间重叠或时钟无效。');
    if(firstStart===null)firstStart=startTime;
    const parsed=readThermalHistory(await readText(file),segment.request,startTime,endTime,{partial:!segment.complete});
    if(lastEnd!==null&&startTime>lastEnd)gaps.push({from:lastEnd,to:startTime,reason:'未保存此区间的运行历史'});
    if(previousSample!==null&&previousSample<startTime) {
      if(parsed.rows.length)parsed.rows[0].breakBefore=true;
    }
    rows.push(...parsed.rows);
    const last=parsed.rows.at(-1)?.time ?? startTime;
    if(last<endTime)gaps.push({from:last,to:endTime,reason:'检查点已接受，历史尚未完整写出'});
    previousSample=parsed.rows.at(-1)?.time ?? previousSample;lastEnd=endTime;
  }
  if(restart && lastEnd!==null && lastEnd!==restart.metadata.time)throw new Error('温度历史终点与联合检查点不同。');
  return {rows,segments:segments.length,availableSteps:rows.length,gaps,
    startTime:firstStart ?? restart?.metadata.time ?? null,endTime:restart?.metadata.time ?? null};
}
module.exports={checkpointDigest,historyForRestart,readThermalTimeline};
