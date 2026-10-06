'use strict';
const test=require('node:test'),assert=require('node:assert/strict'),fs=require('node:fs/promises'),path=require('node:path'),os=require('node:os');
const {runFlowJob}=require('../src/core/flow-job');
const {readFlowTimeline}=require('../src/core/flow-history');
const request={case:'channel',nu:.01,speed:1,maxIterations:100,mode:'transient',dt:.125,steps:2};
const checkpoint=time=>`CARTMESH2D_FLOW_CHECKPOINT 2\nDISCRETIZATION Euler-RC-v2\nCONFIG "channel" .01 1 upwind symmetric 0 reject\nCELLS 1\nFACES 4\nTIME ${time}\nFLUX 4 0 0 0 0\nEND\n`;
const header='step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY\n';
const row=(step,time,accepted=1)=>`${step},${time},.125,${accepted},12,4e-9,1e-10,.5,.2,2,-3\n`;
const base={format:'cartmesh2d-flow-summary-v1',case:'channel',nu:.01,speed:1,cells:1,iterations:12,converged:true,status:'converged',
 continuity:1e-10,globalImbalance:0,globalRelativeImbalance:1e-10,tolerance:1e-6,velocityChange:2e-9,pressureChange:3e-9,momentumResidual:4e-9,
 convection:'upwind',pressurePreconditioner:'ic0',pressureDiscretization:'shared-face-gauss',viscousStress:'symmetric',forceDefinition:'shared-face-newtonian-traction',
 forceX:2,forceY:-3,pressureForceX:1,pressureForceY:-1,discreteForceX:2,discreteForceY:-3,wallForceX:2,wallForceY:-3,wallViscousForceX:1,wallViscousForceY:-2,
 temporalDiscretization:'backward-euler',temporalFaceInterpolation:'old-and-iteration-flux-defect-skew-corrected-v2',dt:.125,requestedSteps:2,completedSteps:2,maxCourant:.5};
async function fixture(t){const root=await fs.mkdtemp(path.join(os.tmpdir(),'flow-job-'));t.after(()=>fs.rm(root,{recursive:true,force:true}));return {root,currentResult:{outputDirectory:root,cm2dPath:path.join(root,'mesh.solver.cm2d')}};}
const fake=({start=0,fail=false,stop=null}={})=>async(_cmd,args)=>{
 const prefix=args[args.indexOf('--output')+1],time=start+(fail?.125:.25);
 await fs.writeFile(prefix+'.checkpoint',checkpoint(time));
 await fs.writeFile(prefix+'.time-history.csv',header+row(1,start+.125)+row(2,start+.25,fail?0:1));
 if(fail){stop?.abort();throw new Error(stop?'cancelled':'candidate rejected');}
 await fs.writeFile(prefix+'.json',JSON.stringify({...base,time,acceptedTime:time}));
 await fs.writeFile(prefix+'.fields.json',JSON.stringify({format:'cartmesh2d-flow-v1',cells:[{id:0,u:.3,v:.4,p:0,speed:.5}]}));
 for(const suffix of ['.vtk','.residuals.csv','.cells.csv','.faces.csv'])await fs.writeFile(prefix+suffix,'original output '+time);
 return {code:0};
};
async function run(f,request,process,controller=new AbortController()){return runFlowJob({currentResult:f.currentResult,mesh:{cells:[{}]},request,runProcess:process,signal:controller.signal,executable:n=>n,onProgress:()=>{},log:()=>{}});}
test('completed flow runs retain distinct originals, and cancelled partial state joins without replacing displayed fields',async t=>{
 const f=await fixture(t);const first=await run(f,request,fake());
 const oldFile=path.join(f.root,first.files.fields),old=await fs.readFile(oldFile);
 const second=await run(f,{...request,resume:true},fake({start:.25}));
 assert.notEqual(first.files.fields,second.files.fields);assert.deepEqual(await fs.readFile(oldFile),old);
 assert.deepEqual((await readFlowTimeline(f.currentResult)).rows.map(r=>r.time),[.125,.25,.375,.5]);
 const stop=new AbortController();await assert.rejects(run(f,{...request,resume:true},fake({start:.5,fail:true,stop}),stop),/cancelled/);
 assert.equal(f.currentResult.flow,second);assert.equal(f.currentResult.flowRestart.metadata.time,.625);
 assert.deepEqual((await readFlowTimeline(f.currentResult)).rows.map(r=>r.time),[.125,.25,.375,.5,.625]);
 assert.deepEqual(await fs.readFile(oldFile),old);
});
test('same-clock replacement of selected checkpoint cannot silently attach an unrelated trajectory',async t=>{
 const f=await fixture(t);await run(f,request,fake());const file=f.currentResult.flowRestart.path;
 await fs.appendFile(file,'changed physical state\n');
 await assert.rejects(run(f,{...request,resume:true},()=>assert.fail('must not spawn')),/已改变/);
 assert.equal(f.currentResult.flow.summary.acceptedTime,.25);
});
