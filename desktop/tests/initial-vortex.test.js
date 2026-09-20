'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const {buildFlowInvocation,validateFlowRequest,flowOutputSuffixes}=require('../src/core/flow');
const {validateInitialVortexOutput}=require('../src/core/initial-vortex');
const {exportGuide}=require('../src/core/export-guide');
const initialVortex={centre:[3,0],radius:1,peakSpeed:-.01};
const request={case:'external',nu:.02,speed:1,maxIterations:100,mode:'transient',dt:.01,steps:2,initialVortex};
const metadata={...initialVortex,definition:'compact-cubic-v1',checkpointSuffix:'.initial.checkpoint'};
test('fresh transient vortex binds all parameters and retains the initial state artifact',()=>{
  const call=buildFlowInvocation('/tmp/mesh.solver.cm2d','/tmp/flow',request);
  assert.deepEqual(call.args.slice(-8),['--initial-vortex-x','3','--initial-vortex-y','0','--initial-vortex-radius','1','--initial-vortex-speed','-0.01']);
  assert.ok(flowOutputSuffixes(request).includes('.initial.checkpoint'));
  assert.ok(!flowOutputSuffixes({...request,initialVortex:undefined}).includes('.initial.checkpoint'));
  for(const update of [{mode:'steady'},{resume:true},{initialVortex:null},
    {initialVortex:{...initialVortex,radius:0}},{initialVortex:{...initialVortex,centre:[3,NaN]}},
    {initialVortex:{...initialVortex,peakSpeed:'0.01'}}])
    assert.throws(()=>validateFlowRequest({...request,...update}),/初始局部涡/);
});
test('initial vortex result cannot silently disappear, change amplitude or reappear on resume',()=>{
  const summary={temporalDiscretization:'backward-euler',initialVortex:metadata};
  const q=validateFlowRequest(request);
  validateInitialVortexOutput(summary,q,0);
  for(const [value,expected,time] of [[{...summary,initialVortex:undefined},q,0],
    [{...summary,initialVortex:{...metadata,peakSpeed:.01}},q,0],[summary,q,.01],
    [summary,{...q,initialVortex:undefined},0],[{...summary,temporalDiscretization:undefined},q,0]])
    assert.throws(()=>validateInitialVortexOutput(value,expected,time),/初始局部涡/);
  const guide=exportGuide({result:{counts:{cells:1},gates:{}},flow:{summary}});
  assert.match(guide,/initial.checkpoint/);assert.match(guide,/续算不再施加扰动/);
});
