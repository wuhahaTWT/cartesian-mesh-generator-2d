'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const {createFlowCaseDocument,serializeFlowCase,parseFlowCaseDocument}=require('../src/core/flow-case');
const {buildFlowInvocation}=require('../src/core/flow');
const mesh=['CM2D 1','VERTICES 4','0 0 0','1 1 0','2 1 1','3 0 1','EDGES 4',
  '0 0 1 0 -1 2','1 1 2 0 -1 2','2 2 3 0 -1 2','3 3 0 0 -1 2',
  'CELLS 1','0 0 6 1 4 0 1 2 3 4 0 1 2 3'].join('\n');
const boundaryDefinition={cells:1,faces:4,records:[
  {face:0,owner:0,x:.5,y:0,sx:0,sy:-1,type:'wall',name:'下壁',u:0,v:0,p:0},
  {face:1,owner:0,x:1,y:.5,sx:1,sy:0,type:'pressure-outlet',name:'出口',u:0,v:0,p:7.25},
  {face:2,owner:0,x:.5,y:1,sx:0,sy:1,type:'wall',name:'上壁',u:0,v:0,p:0},
  {face:3,owner:0,x:0,y:.5,sx:-1,sy:0,type:'velocity-inlet',name:'入口',u:1,v:.1,p:0}]};
const common={case:'custom',nu:.1,speed:1,maxIterations:500,tolerance:1e-9,convection:'face-limited-linear',
  pressurePreconditioner:'aggregation',boundaryDefinition};
const temporal={dt:.04,endTime:1,minDt:.00001,maxCourant:.3,maxRetries:8,maxSteps:2000,steps:2,
  initialVortex:{centre:[.5,.5],radius:.1,peakSpeed:-.05}};
test('flow case restores native solver arguments for each time mode and all named values',()=>{
  for(const mode of ['steady','transient','adaptive']) {
    const input={...common,mode,linearPolicy:'adaptive',velocityRelaxation:.8,...(mode==='steady'?{steadyAcceleration:'anderson'}:temporal)};
    const saved=createFlowCaseDocument(input,Buffer.from(mesh));
    const loaded=parseFlowCaseDocument(serializeFlowCase(saved),mesh);
    assert.deepEqual(loaded,saved);
    const invoke=q=>buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/flow',q,null,'/tmp/input.boundaries');
    assert.deepEqual(invoke(loaded.request),invoke(input));
    assert.equal(loaded.request.boundaryDefinition.records[1].p,7.25);
    assert.equal(loaded.request.boundaryDefinition.records[3].v,.1);
    assert.equal(loaded.mesh.cells,1);assert.match(loaded.mesh.sha256,/^[a-f0-9]{64}$/);
    assert.equal(loaded.request.resume,false);
  }
});
test('case rejects a different mesh even with equal counts and checks boundary geometry independently',()=>{
  const saved=createFlowCaseDocument(common,mesh);
  assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(saved),mesh.replace('1 1 0','1 1.1 0')),/网格不匹配/);
  const corrupt=structuredClone(saved);corrupt.request.boundaryDefinition.records[0].sy=1;
  assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(corrupt),mesh),/法向/);
  corrupt.request.boundaryDefinition.records[0].sy=-1;corrupt.request.boundaryDefinition.records.shift();
  assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(corrupt),mesh),/全部边界/);
});
test('case never silently accepts restart, unknown or inactive solver settings, or invalid values',()=>{
  const saved=createFlowCaseDocument({...common,...temporal,mode:'adaptive'},mesh);
  const mutations=[d=>d.format='cartmesh2d-flow-case-v5',d=>d.request.resume=true,
    d=>d.request.steadyAcceleration='anderson',d=>delete d.request.steadyAcceleration,
    d=>delete d.request.tolerance,d=>d.request.tolerance=1e-4,
    d=>delete d.request.linearPolicy,d=>delete d.request.velocityRelaxation,d=>d.request.linearPolicy='unknown',d=>d.request.velocityRelaxation=0,d=>d.request.velocityRelaxation='0.6',
    d=>delete d.request.nu,d=>d.request.nu='0.1',d=>d.request.nu=-1,d=>d.request.dt=0,
    d=>d.request.turbulence='sst',d=>d.request.steps=10,d=>d.request.command='/bin/sh',
    d=>d.request.boundaryDefinition.records[0].path='/tmp/other',
    d=>d.request.initialVortex.decay=5,d=>d.mesh.cells=2,d=>d.request.maxRetries=100];
  for(const mutate of mutations){const d=structuredClone(saved);mutate(d);assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(d),mesh));}
  assert.throws(()=>createFlowCaseDocument({...common,mode:'transient',dt:.1,steps:1,resume:true},mesh),/取消续算/);
  assert.throws(()=>parseFlowCaseDocument('null',mesh));
  assert.throws(()=>parseFlowCaseDocument('{broken',mesh),/JSON/);
});

test('legacy case files retain the exact old default while v2 preserves tighter stopping controls',()=>{
  const saved=createFlowCaseDocument(common,mesh);assert.equal(saved.format,'cartmesh2d-flow-case-v4');
  const legacy=structuredClone(saved);legacy.format='cartmesh2d-flow-case-v1';delete legacy.request.tolerance;delete legacy.request.steadyAcceleration;delete legacy.request.linearPolicy;delete legacy.request.velocityRelaxation;
  const restored=parseFlowCaseDocument(serializeFlowCase(legacy),mesh);
  assert.equal(restored.request.tolerance,1e-6);assert.equal(restored.format,saved.format);
  assert.equal(restored.request.steadyAcceleration,'none');
  const v2=structuredClone(saved);v2.format='cartmesh2d-flow-case-v2';delete v2.request.steadyAcceleration;delete v2.request.linearPolicy;delete v2.request.velocityRelaxation;
  assert.equal(parseFlowCaseDocument(serializeFlowCase(v2),mesh).request.steadyAcceleration,'none');
  v2.request.steadyAcceleration='anderson';assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(v2),mesh),/旧版工况/);
  assert.equal(parseFlowCaseDocument(serializeFlowCase(saved),mesh).request.tolerance,1e-9);
  const v3=structuredClone(saved);v3.format='cartmesh2d-flow-case-v3';delete v3.request.linearPolicy;delete v3.request.velocityRelaxation;v3.request.steadyAcceleration='anderson';
  const old=parseFlowCaseDocument(serializeFlowCase(v3),mesh).request;assert.equal(old.linearPolicy,'strict');assert.equal(old.velocityRelaxation,.6);assert.equal(old.steadyAcceleration,'anderson');
  for(const key of ['linearPolicy','velocityRelaxation']){const bad=structuredClone(v3);bad.request[key]=saved.request[key];assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(bad),mesh),/旧版工况/);}
  legacy.request.tolerance=1e-9;
  assert.throws(()=>parseFlowCaseDocument(serializeFlowCase(legacy),mesh),/旧版工况/);
});
