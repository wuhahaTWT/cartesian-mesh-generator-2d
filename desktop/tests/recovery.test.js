'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const fs=require('node:fs/promises'),path=require('node:path'),os=require('node:os');
const {RecoveryStore,atomicJson}=require('../src/core/recovery');
const {checkpointDigest,readThermalTimeline}=require('../src/core/thermal-history');
const MESH='CM2D 1\nVERTICES 4\n0 0 0\n1 1 0\n2 1 1\n3 0 1\nEDGES 4\n0 0 1 0 -1 2\n1 1 2 0 -1 2\n2 2 3 0 -1 2\n3 3 0 0 -1 2\nCELLS 1\n0 0 0 1 4 0 1 2 3 4 0 1 2 3\nAUDIT 0\n';
const REQUEST={case:'channel',nu:.1,speed:.2,maxIterations:100,convection:'upwind',pressurePreconditioner:'ic0',mode:'transient',dt:.01,steps:2,diffusivity:.1,initial:300,source:0,scalarConvection:'upwind',boundaries:Object.fromEntries(['wall','inlet','outlet','top','bottom'].map(k=>[k,{kind:'value',value:300,inflowValue:300}]))};
const checkpoint=time=>`CARTMESH2D_THERMAL_CHECKPOINT 1\nCOUPLING new-time-flux-Euler-v1\nSCALAR 1 300\nFLOW\nCARTMESH2D_FLOW_CHECKPOINT 2\nTIME ${time}\nFLUX 4 0 0 0 0\n`;
const history=times=>'step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant\n'+times.map((time,i)=>`${i+1},${time},1,2,1e-10,1e-12,2,1e-10,300,0,.1\n`).join('');
async function fixture(t) {
 const base=await fs.mkdtemp(path.join(os.tmpdir(),'cartmesh-recovery-'));t.after(()=>fs.rm(base,{recursive:true,force:true}));
 const session=await fs.mkdtemp(path.join(base,'session-')),root=path.join(session,'mesh');await fs.mkdir(root);
 const prefix=path.join(root,'rectangle');await fs.writeFile(prefix+'.xy','0 0\n1 0\n1 1\n0 1\n');await fs.writeFile(prefix+'.solver.cm2d',MESH);
 const old=path.join(root,'old.checkpoint'),csv=path.join(root,'old.csv');await fs.writeFile(old,checkpoint(.02));await fs.writeFile(csv,history([.01,.02]));
 const current={outputDirectory:root,prefix,cm2dPath:prefix+'.solver.cm2d',job:{geometryPath:'rectangle.xy',outputDirectory:root,method:'cutcell',fluidRegion:'interior'},result:{counts:{cells:1}},
  thermalRestart:{path:old,sha256:checkpointDigest(checkpoint(.02)),metadata:{time:.02,request:REQUEST},history:[{file:csv,startTime:0,endTime:.02,complete:true,request:REQUEST}]}};
 const writer=new RecoveryStore(base,session);await writer.checkpoint(current,{inputs:{thermalResume:true}});
 const dir=path.join(root,'run');await fs.mkdir(dir);await fs.writeFile(path.join(dir,'thermal.thermal.checkpoint'),checkpoint(.04));await fs.writeFile(path.join(dir,'thermal.thermal-history.csv'),history([.03,.04]));
 await writer.beginThermal(current,{directory:'run',startTime:.02,request:{...REQUEST,resume:true}});
 const file=(await fs.readdir(session)).find(f=>f.startsWith('recovery-'));
 const record=JSON.parse(await fs.readFile(path.join(session,file),'utf8'));record.ownerPid=2147483647;await atomicJson(path.join(session,file),record);
 const next=await fs.mkdtemp(path.join(base,'session-'));return {base,root,session,next,file,record,id:path.basename(session)+'/'+file,reader:new RecoveryStore(base,next)};
}
test('recovery opens a new project copy and joins its native-checked interrupted accepted state',async t=>{
 const f=await fixture(t),source=await fs.readFile(path.join(f.session,f.file));
 let checks=0;const restored=await f.reader.restore(f.id,async(current,file,request)=>{checks++;assert.ok(file.startsWith(f.next));assert.equal(request.resume,true);});
 assert.equal(checks,1);assert.equal(restored.thermalRestart.metadata.time,.04);
 assert.deepEqual((await readThermalTimeline(restored)).rows.map(r=>r.time),[.01,.02,.03,.04]);
 assert.equal(restored.projectUi.inputs.thermalResume,true);
 assert.deepEqual(await fs.readFile(path.join(f.session,f.file)),source,'original recovery record remains intact');
 assert.equal(await fs.readFile(path.join(f.root,'run/thermal.thermal.checkpoint'),'utf8'),checkpoint(.04));
});
test('failed native validation retains the earlier accepted state and exposes the failure',async t=>{
 const f=await fixture(t);const restored=await f.reader.restore(f.id,async()=>{throw new Error('native boundary mismatch');});
 assert.equal(restored.thermalRestart.metadata.time,.02);assert.match(restored.recoveryNotice,/native boundary mismatch/);
 assert.deepEqual((await readThermalTimeline(restored)).rows.map(r=>r.time),[.01,.02]);
});
test('recovery refuses a live owner or live native child before copying any project',async t=>{
 const f=await fixture(t);f.record.ownerPid=process.pid;await atomicJson(path.join(f.session,f.file),f.record);
 await assert.rejects(f.reader.restore(f.id,()=>assert.fail('must not validate a live project')),/仍在运行/);
 f.record.ownerPid=2147483647;await atomicJson(path.join(f.session,f.file),f.record);
 await fs.writeFile(path.join(f.root,'run/native-process.json'),JSON.stringify({pid:process.pid}));
 await assert.rejects(f.reader.restore(f.id,()=>assert.fail('must not read live native files')),/仍在运行/);
 assert.deepEqual(await fs.readdir(f.next),[]);
});
test('cancelling native recovery validation does not publish a fallback session',async t=>{
 const f=await fixture(t),stop=new AbortController();
 await assert.rejects(f.reader.restore(f.id,async()=>{stop.abort();throw new Error('cancelled');},stop.signal),/abort/i);
 assert.equal(f.reader.records.size,0);
 assert.equal(JSON.parse(await fs.readFile(path.join(f.session,f.file),'utf8')).activeThermal.startTime,.02);
});
