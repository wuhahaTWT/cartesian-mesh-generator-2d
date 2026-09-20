'use strict';

const near = (a,b) => Math.abs(a-b) <= 1e-12+1e-9*Math.max(Math.abs(a),Math.abs(b));
const numeric = (value,name) => {
  if (typeof value!=='number' || !Number.isFinite(value)) throw new Error(`自动步长 ${name} 必须为有限数。`);
  return value;
};

function validateAdaptiveSummary(q, startTime, request=null) {
  for (const name of ['startTime','targetTime','maximumTimeStep','minimumTimeStep','targetCourant',
    'maximumRetries','maximumAcceptedSteps','attemptCount','rejectedSteps']) numeric(q[name],name);
  if (q.timeStepControl!=='adaptive-cfl-retry' || !q.converged || Object.hasOwn(q,'requestedSteps')
      || !near(q.startTime,startTime) || !(q.targetTime>startTime) || !near(q.time,q.targetTime)
      || !near(q.acceptedTime,q.time) || !(q.minimumTimeStep>0) || q.maximumTimeStep<q.minimumTimeStep
      || !(q.targetCourant>0) || q.maxCourant>q.targetCourant || q.dt>q.maximumTimeStep*(1+1e-12)
      || !Number.isInteger(q.maximumRetries) || q.maximumRetries<0 || q.maximumRetries>30
      || !Number.isInteger(q.maximumAcceptedSteps) || q.maximumAcceptedSteps<1 || q.maximumAcceptedSteps>1000000
      || !Number.isInteger(q.completedSteps) || q.completedSteps<1 || q.completedSteps>q.maximumAcceptedSteps
      || !Number.isInteger(q.rejectedSteps) || q.rejectedSteps<0
      || !Number.isInteger(q.attemptCount) || q.attemptCount!==q.completedSteps+q.rejectedSteps)
    throw new Error('自动步长摘要的时间、控制参数或接受记录不一致。');
  if (request) {
    for (const [key,input] of [['maximumTimeStep','dt'],['minimumTimeStep','minDt'],['targetTime','endTime'],
      ['targetCourant','maxCourant'],['maximumRetries','maxRetries'],['maximumAcceptedSteps','maxSteps']])
      if (!near(q[key],request[input])) throw new Error('自动步长结果与请求控制参数不一致。');
  }
}

function validateAttemptHistory(text, q, history) {
  const keys=['attempt','step','startTime','time','dt','accepted','reason','innerConverged','innerIterations',
    'momentumResidual','continuity','velocityChange','pressureChange','maxCourant'];
  const lines=text.trim().split(/\r?\n/);
  if (lines.shift()!==keys.join(',')) throw new Error('自动步长试算历史表头无效。');
  let time=q.startTime, accepted=0, rejected=0, retries=0, previous=null;
  const rows=lines.map((line,i)=>{
    const values=line.split(',');
    if (values.length!==keys.length || values.some(v=>v.trim()==='')) throw new Error('自动步长试算历史列数或数值无效。');
    const r=Object.fromEntries(keys.map((k,j)=>[k,k==='reason'?values[j]:numeric(Number(values[j]),k)]));
    if (r.attempt!==i+1 || r.step!==accepted+1 || ![0,1].includes(r.accepted) || ![0,1].includes(r.innerConverged)
        || !near(r.startTime,time) || !near(r.time,time+r.dt) || !(r.dt>0) || !(r.time>time)
        || r.dt>q.maximumTimeStep*(1+1e-12) || r.time>q.targetTime+1e-12
        || (r.dt<q.minimumTimeStep && !near(r.dt,q.targetTime-time))
        || !Number.isInteger(r.innerIterations) || r.innerIterations<1
        || ['momentumResidual','continuity','velocityChange','pressureChange','maxCourant'].some(k=>r[k]<0))
      throw new Error('自动步长试算历史的顺序、时间或数值无效。');
    if (r.innerConverged && (r.innerIterations<10 || r.continuity>=1e-8 ||
        ['momentumResidual','velocityChange','pressureChange'].some(k=>r[k]>=q.tolerance)))
      throw new Error('自动步长试算未满足收敛条件。');
    if (previous) {
      const factor=previous.maxCourant>0?Math.min(.5,.8*q.targetCourant/previous.maxCourant):.5;
      const expected=Math.max(Math.min(q.minimumTimeStep,q.targetTime-time),previous.dt*factor);
      if (!(r.dt<previous.dt) || !near(r.dt,expected)) throw new Error('自动步长重试未按记录缩减步长。');
    }
    if (r.accepted) {
      if (r.reason!=='accepted' || !r.innerConverged || r.maxCourant>q.targetCourant)
        throw new Error('自动步长错误接受了未收敛或 CFL 超限的试算。');
      const h=history[accepted];
      if (!h || ['step','time','dt','innerIterations','momentumResidual','continuity','maxCourant'].some(k=>!near(r[k],h[k])))
        throw new Error('自动步长试算与已接受时间历史不一致。');
      accepted++;time=r.time;retries=0;previous=null;
    } else {
      if ((r.reason==='courant' && (!r.innerConverged || r.maxCourant<=q.targetCourant))
          || (r.reason==='nonconverged' && r.innerConverged) || !['courant','nonconverged'].includes(r.reason))
        throw new Error('自动步长拒绝原因与试算状态不一致。');
      rejected++;retries++;previous=r;
      if (retries>q.maximumRetries) throw new Error('自动步长超过重试次数上限。');
    }
    return r;
  });
  if (previous || rows.length!==q.attemptCount || accepted!==q.completedSteps || rejected!==q.rejectedSteps || !near(time,q.targetTime))
    throw new Error('自动步长试算记录不完整。');
  return rows;
}

module.exports={validateAdaptiveSummary,validateAttemptHistory};
