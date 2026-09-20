'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs/promises');
const os=require('node:os');
const path=require('node:path');
const {parseBoundaryDefinition,serializeBoundaryDefinition,validateBoundaryMesh,conditions,sameConditions}=require('../src/core/flow-boundaries');
const {readCheckpointMetadata}=require('../src/core/flow-checkpoint');
const {buildFlowInvocation,flowOutputSuffixes,validateFlowRequest}=require('../src/core/flow');
const mesh={vertices:[[0,0],[1,0],[1,1],[0,1]],cells:[{vertices:[0,1,2,3]}],
  edges:[{a:0,b:1,owner:0,neighbour:-1},{a:2,b:1,owner:0,neighbour:-1},
         {a:2,b:3,owner:0,neighbour:-1},{a:3,b:0,owner:0,neighbour:-1}],
  bounds:{minX:0,minY:0,maxX:1,maxY:1}};
const definition={cells:1,faces:4,records:[
  {face:0,owner:0,x:.5,y:0,sx:0,sy:-1,type:'wall',name:'曲壁',u:0,v:0,p:0},
  {face:1,owner:0,x:1,y:.5,sx:1,sy:0,type:'pressure-outlet',name:'outlet',u:0,v:0,p:7.25},
  {face:2,owner:0,x:.5,y:1,sx:0,sy:1,type:'wall',name:'曲壁',u:0,v:0,p:0},
  {face:3,owner:0,x:0,y:.5,sx:-1,sy:0,type:'velocity-inlet',name:'入口 \\ A',u:1,v:0,p:0}]};
const copy=()=>structuredClone(definition);

test('named boundary input preserves names, orientations and nonuniform values',()=>{
  const d=copy();
  assert.deepEqual(parseBoundaryDefinition(serializeBoundaryDefinition(d)),d);
  assert.deepEqual(validateBoundaryMesh(d,mesh,1),d);
  assert(sameConditions(d.records,[...d.records].reverse()));
  assert.throws(()=>parseBoundaryDefinition(serializeBoundaryDefinition(d)+'TRAILING'),/边界/);
  const wrong=copy();wrong.records[0].sy=1;
  assert.throws(()=>validateBoundaryMesh(wrong,mesh,1),/法向/);
  wrong.records[0].sy=-1;wrong.records[0].owner=1;
  assert.throws(()=>validateBoundaryMesh(wrong,mesh,1),/单元/);
});

test('named boundary coverage and physical conditions fail explicitly',()=>{
  const mutations=[d=>d.records.pop(),d=>d.records.push({...d.records[0]}),
    d=>d.records[3].u=-1,d=>d.records[0].u=.1,d=>d.records[1].v=1,
    d=>d.records[3].name='bad,name',d=>d.records[3].u=NaN,
    d=>d.records[3].name='outlet',d=>Object.assign(d.records[0],{type:'moving-wall',name:'moving',v:.1}),
    d=>d.records[3].type='symmetry',d=>d.records[0].face=8];
  for(const mutate of mutations){const d=copy();mutate(d);assert.throws(()=>validateBoundaryMesh(d,mesh,1),/边界/);}
  const moving=copy();Object.assign(moving.records[0],{type:'moving-wall',name:'moving',u:.2});
  validateBoundaryMesh(moving,mesh,1);
});

test('custom requests bind explicit input paths, preserve conditions and export snapshot',()=>{
  const request={case:'custom',nu:.1,speed:1,maxIterations:400,boundaryDefinition:copy()};
  assert.throws(()=>buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/run',request),/路径/);
  const invocation=buildFlowInvocation('/tmp/m.solver.cm2d','/tmp/run',request,null,'/tmp/input.boundaries');
  assert.deepEqual(invocation.args.slice(-2),['--boundary','/tmp/input.boundaries']);
  assert.deepEqual(invocation.request.boundaryDefinition,definition);
  assert(flowOutputSuffixes(request).includes('.boundaries'));
  assert(!flowOutputSuffixes({case:'duct'}).includes('.boundaries'));
  assert.throws(()=>validateFlowRequest({...request,outletBackflow:'normal-inlet'}),/回流/);
});

const checkpoint=()=>[
 'CARTMESH2D_FLOW_CHECKPOINT 4','DISCRETIZATION Euler-RC-v2',
 'CONFIG "custom" 0.1 1 upwind symmetric 0 reject','FACE_VISCOSITY 0','BOUNDARIES 4',
 ...conditions(definition.records).map(b=>`BOUNDARY ${b.face} ${b.type} "${b.name.replace(/\\/g,'\\\\')}" ${b.u} ${b.v} ${b.p}`),
 'CELLS 1','CELL 0 0.5 0.5 1 4 0 1 2 3','FACES 4',
 ...definition.records.map(b=>`FACE ${b.face} ${b.owner} - 1 ${b.x} ${b.y} ${b.sx} ${b.sy} 0 0 2 0`),
 'TIME 0.04','U 1 0','V 1 0','P 1 0','FLUX 4 0 0 0 0','END',''].join('\n');
async function metadata(text){const dir=await fs.mkdtemp(path.join(os.tmpdir(),'custom-checkpoint-'));
 try{const file=path.join(dir,'state.checkpoint');await fs.writeFile(file,text);return await readCheckpointMetadata(file);}
 finally{await fs.rm(dir,{recursive:true,force:true});}}
test('custom checkpoint restores every named condition and geometry for the editor',async()=>{
 const result=await metadata(checkpoint());
 assert.equal(result.case,'custom');assert.equal(result.time,.04);
 assert.deepEqual(result.boundaryDefinition,definition);
 validateBoundaryMesh(result.boundaryDefinition,mesh,result.speed);
 for(const text of [checkpoint().replace('FACE_VISCOSITY 0','FACE_VISCOSITY 4 1 1 1 1'),
   checkpoint().replace('BOUNDARIES 4','BOUNDARIES 3'),checkpoint().replace('FACE 3 0 -','FACE 3 0 1'),
   checkpoint().replace('CARTMESH2D_FLOW_CHECKPOINT 4','CARTMESH2D_FLOW_CHECKPOINT 2'),
   checkpoint().replace('FACE 3 0 - 1 0 0.5 -1 0 0 0 2 0\n','')]){
     await assert.rejects(metadata(text));
   }
});
