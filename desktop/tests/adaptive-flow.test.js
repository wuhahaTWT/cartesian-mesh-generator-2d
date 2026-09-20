'use strict';
const test=require('node:test'), assert=require('node:assert/strict');
const {buildFlowInvocation,validateFlowRequest,validateFlowOutput,validateTimeHistory,validateAttemptHistory,
  flowOutputSuffixes,parseFlowProgress}=require('../src/core/flow');
const {exportGuide}=require('../src/core/export-guide');
const request={case:'cavity',nu:.1,speed:1,maxIterations:100,mode:'adaptive',dt:.04,endTime:.08,minDt:.001,maxCourant:.15};
const q={format:'cartmesh2d-flow-summary-v1',case:'cavity',nu:.1,speed:1,cells:1,iterations:12,
  status:'converged',converged:true,continuity:1e-10,globalImbalance:0,globalRelativeImbalance:0,
  tolerance:1e-6,velocityChange:1e-9,pressureChange:1e-9,momentumResidual:1e-9,convection:'upwind',
  pressurePreconditioner:'ic0',pressureDiscretization:'shared-face-gauss',viscousStress:'symmetric',
  forceDefinition:'shared-face-newtonian-traction',forceX:0,forceY:0,pressureForceX:0,pressureForceY:0,
  discreteForceX:0,discreteForceY:0,wallForceX:0,wallForceY:0,wallViscousForceX:0,wallViscousForceY:0,
  temporalDiscretization:'backward-euler',temporalFaceInterpolation:'old-and-iteration-flux-defect-skew-corrected-v2',
  timeStepControl:'adaptive-cfl-retry',time:.08,acceptedTime:.08,startTime:.02,dt:.02,
  targetTime:.08,maximumTimeStep:.04,minimumTimeStep:.001,targetCourant:.15,
  maximumRetries:10,maximumAcceptedSteps:100000,completedSteps:3,attemptCount:4,rejectedSteps:1,maxCourant:.1};
const fields={format:'cartmesh2d-flow-v1',cells:[{id:0,u:0,v:0,p:0,speed:0}]};
const history='step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY\n'+
  '1,.04,.02,1,12,1e-9,1e-10,.1,0,0,0\n2,.06,.02,1,12,1e-9,1e-10,.1,0,0,0\n3,.08,.02,1,12,1e-9,1e-10,.1,0,0,0\n';
const attempts='attempt,step,startTime,time,dt,accepted,reason,innerConverged,innerIterations,momentumResidual,continuity,velocityChange,pressureChange,maxCourant\n'+
  '1,1,.02,.06,.04,0,courant,1,12,1e-9,1e-10,1e-9,1e-9,.2\n'+
  '2,1,.02,.04,.02,1,accepted,1,12,1e-9,1e-10,1e-9,1e-9,.1\n'+
  '3,2,.04,.06,.02,1,accepted,1,12,1e-9,1e-10,1e-9,1e-9,.1\n'+
  '4,3,.06,.08,.02,1,accepted,1,12,1e-9,1e-10,1e-9,1e-9,.1\n';

test('adaptive request binds absolute target and complete controls without fixed step count',()=>{
  const inv=buildFlowInvocation('/tmp/a.solver.cm2d','/tmp/a',request);
  assert.ok(inv.args.includes('--end-time'));assert.ok(!inv.args.includes('--steps'));
  assert.equal(inv.request.maxRetries,10);assert.equal(inv.request.maxSteps,100000);
  assert.equal(inv.request.steps,undefined);assert.ok(flowOutputSuffixes(request).includes('.attempt-history.csv'));
  for(const update of [{endTime:0},{minDt:.05},{maxCourant:0},{maxRetries:31},{maxSteps:0}])
    assert.throws(()=>validateFlowRequest({...request,...update}),/自动步长/);
});
test('adaptive accepted output and rejected trials bind to the same request and timeline',()=>{
  const result=validateFlowOutput(q,fields,1,request,.02);
  const rows=validateTimeHistory(history,result.summary,.02);
  assert.equal(validateAttemptHistory(attempts,result.summary,rows).length,4);
  const guide=exportGuide({result:{counts:{cells:1},gates:{}},flow:{summary:q}});
  assert.match(guide,/attempt-history/);assert.match(guide,/绝对物理时间/);
  for(const update of [{targetTime:.09},{maximumTimeStep:.05},{startTime:0},{requestedSteps:3},
    {completedSteps:1.5},{timeStepControl:'mystery'},{maxCourant:.2},{attemptCount:3}])
    assert.throws(()=>validateFlowOutput({...q,...update},fields,1,request,.02),/自动步长|时间步控制/);
  for(const wrong of [attempts.replace(',courant,',',accepted,'),attempts.replace('2,1,.02','2,1,.03'),
    attempts.replace(',1e-9,.1',',1e-2,.1'),attempts.split('\n').slice(0,-2).join('\n'),
    attempts.replace('2,1,.02,.04,.02','2,1,.02,.05,.03')])
    assert.throws(()=>validateAttemptHistory(wrong,q,rows),/自动步长/);
});
test('adaptive retry progress never masquerades as an accepted step',()=>{
  const r={type:'flow-time-retry',step:1,attempt:1,acceptedTime:0,candidateTime:.04,dt:.04,nextDt:.02,maxCourant:.2,reason:'courant'};
  assert.deepEqual(parseFlowProgress(JSON.stringify(r)),r);
  for(const update of [{nextDt:.05},{acceptedTime:-1},{reason:'accepted'},{candidateTime:.02}])
    assert.throws(()=>parseFlowProgress(JSON.stringify({...r,...update})),/重试/);
});
