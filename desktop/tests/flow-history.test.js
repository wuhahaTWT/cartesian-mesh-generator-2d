'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const {readFlowHistory}=require('../src/core/flow');
const {readFlowTimeline,historyForRestart}=require('../src/core/flow-history');
const request={case:'channel',nu:.01,speed:1,maxIterations:100,mode:'transient',dt:.125,steps:2};
const header='step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY\n';
const row=(step,time,accepted=1)=>`${step},${time},.125,${accepted},12,1e-9,1e-10,.1,.2,1,-1\n`;
test('flow timeline joins accepted segments, breaks missing tails, and excludes failed candidates',async()=>{
 const files={'/p/first':header+row(1,.125)+row(2,.25),'/p/interrupted':header+row(1,.375)+'2,.5',
  '/p/resume':header+row(1,.625)+row(2,.75,0)};
 const current={outputDirectory:'/p',flowRestart:{metadata:{time:.625},history:[
  {file:'/p/first',startTime:0,endTime:.25,complete:true,request},
  {file:'/p/interrupted',startTime:.25,endTime:.5,complete:false,request},
  {file:'/p/resume',startTime:.5,endTime:.625,complete:false,request}]}};
 const data=await readFlowTimeline(current,async p=>files[p]);
 assert.deepEqual(data.rows.map(r=>r.time),[.125,.25,.375,.625]);
 assert.equal(data.rows.at(-1).breakBefore,true);assert.deepEqual(data.gaps.map(g=>[g.from,g.to]),[[.375,.5]]);
 assert.equal(data.availableSteps,4);
 await assert.rejects(readFlowTimeline({...current,flowRestart:{...current.flowRestart,metadata:{time:.75}}},async p=>files[p]),/检查点/);
 await assert.rejects(readFlowTimeline({...current,flowRestart:{history:[{...current.flowRestart.history[0],file:'/escape'}]}},async()=>''),/越出/);
});
test('legacy pre-checkpoint history tails are visible diagnostics and never accepted beyond restart state',()=>{
 const text=header+row(1,.125)+row(2,.25);
 const parsed=readFlowHistory(text,request,0,.125,{partial:true});
 assert.deepEqual(parsed.rows.map(r=>r.time),[.125]);assert.equal(parsed.uncommittedSamples,1);
 assert.throws(()=>readFlowHistory(text,request,0,.125),/终点/);
 assert.throws(()=>readFlowHistory((header+row(1,.125)).replace('1e-9','1e-3'),request,0,.125,{partial:true}),/停止条件/);
 assert.deepEqual(historyForRestart({outputDirectory:'/p',flowRestart:{path:'/p/input',metadata:{time:.25}},flow:{files:{checkpoint:'different'},history:[{time:.125,dt:.125}]}}),[]);
});
