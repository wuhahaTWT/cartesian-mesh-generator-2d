'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const fs=require('node:fs/promises'),path=require('node:path'),os=require('node:os');
const {readThermalTimeline,historyForRestart}=require('../src/core/thermal-history');
const {readThermalHistory}=require('../src/core/thermal');
const REQUEST={case:'channel',nu:.1,speed:.2,maxIterations:100,convection:'upwind',pressurePreconditioner:'ic0',
  mode:'transient',dt:.1,steps:2,diffusivity:.1,initial:300,source:0,scalarConvection:'upwind',
  boundaries:Object.fromEntries(['wall','inlet','outlet','top','bottom'].map(k=>[k,{kind:'value',value:300,inflowValue:300}]))};
const HEADER='step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant\n';
const csv=times=>HEADER+times.map((time,i)=>`${i+1},${time},1,12,1e-10,1e-12,5,1e-10,600,0,.2\n`).join('');

test('timeline joins accepted runs, keeps failed-run progress and marks a crash history gap',async t=>{
  const root=await fs.mkdtemp(path.join(os.tmpdir(),'thermal history '));t.after(()=>fs.rm(root,{recursive:true,force:true}));
  const history=[];
  for(const [name,startTime,endTime,times,dt,complete] of [
    ['first',0,.2,[.1,.2],.1,true],['failed',.2,.3,[.25],.05,false],['resumed',.3,.4,[.35,.4],.05,true]]) {
    const file=path.join(root,name+'.csv');await fs.writeFile(file,csv(times)+(complete?'':'2,0.3,1,12'));
    history.push({file,startTime,endTime,request:{...REQUEST,dt},complete});
  }
  const current={outputDirectory:root,thermalRestart:{metadata:{time:.4},history}};
  const result=await readThermalTimeline(current);
  assert.deepEqual(result.rows.map(r=>r.time),[.1,.2,.25,.35,.4]);
  assert.equal(result.availableSteps,5);assert.equal(result.segments,3);
  assert.deepEqual(result.gaps,[{from:.25,to:.3,reason:'检查点已接受，历史尚未完整写出'}]);
  assert.equal(result.rows[3].breakBefore,true);
  assert.ok(!result.rows[2].breakBefore);
  await fs.writeFile(history[1].file,csv([.25]).replace(',1,12,',',0,12,'));
  await assert.rejects(readThermalTimeline(current),/监测只接受/,'a completed rejected row is never silently plotted');
});

test('completed history is strict, interrupted tail can be missing, overlaps and outside paths fail',async()=>{
  assert.throws(()=>readThermalHistory(csv([.1]),REQUEST,0,.2),/未到/);
  assert.deepEqual(readThermalHistory('',REQUEST,0,.2,{partial:true}).rows,[]);
  const base={file:'/tmp/project/run.csv',startTime:0,endTime:.2,request:REQUEST,complete:true};
  const state=history=>({outputDirectory:'/tmp/project',thermalRestart:{history,metadata:{time:.2}}});
  await assert.rejects(readThermalTimeline(state([{...base,file:'/tmp/outside.csv'}]),async()=>csv([.1,.2])),/越出/);
  await assert.rejects(readThermalTimeline(state([base,{...base,startTime:.1}]),async()=>csv([.1,.2])),/重叠/);
  await assert.rejects(readThermalTimeline({...state([base]),thermalRestart:{history:[base],metadata:{time:.3}}},async()=>csv([.1,.2])),/终点/);
});

test('legacy projects seed only the history tied to the selected checkpoint, not an equal clock',()=>{
  const current={outputDirectory:'/tmp/project',thermal:{summary:{time:.2,startTime:0},request:REQUEST,
    files:{'.thermal.checkpoint':'run/thermal.checkpoint','.thermal-history.csv':'run/history.csv'}},
    thermalRestart:{path:'/tmp/project/run/thermal.checkpoint',metadata:{time:.2}}};
  assert.equal(historyForRestart(current).length,1);
  current.thermalRestart.path='/tmp/project/other/thermal.checkpoint';
  assert.deepEqual(historyForRestart(current),[]);
});
