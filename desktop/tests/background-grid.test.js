'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const {parseBackgroundGrid,requireFluidMesh}=require('../src/core/background-grid');
const {validateJob,buildInvocation}=require('../src/core/job');
const {exportGuide}=require('../src/core/export-guide');
function fixture(){return {format:'cartmesh2d-background-v1',mode:'uniform',solver_ready:false,retains_solid_interior:true,
  classification_names:['outside','inside','intersected'],domain:[0,0,1,1],boundary_loops:[[[.2,.2],[.8,.2],[.5,.8]]],counts:[1,1,2],cells:[
    {id:0,level:1,ix:0,iy:0,classification:0,bounds:[0,0,.5,.5]},
    {id:1,level:1,ix:1,iy:0,classification:1,bounds:[.5,0,1,.5]},
    {id:2,level:1,ix:0,iy:1,classification:2,bounds:[0,.5,.5,1]},
    {id:3,level:1,ix:1,iy:1,classification:2,bounds:[.5,.5,1,1]}]};}
test('background preview preserves every full cell and class, not fluid topology',()=>{
  const mesh=parseBackgroundGrid(JSON.stringify(fixture()));
  assert.equal(mesh.cells.length,4);assert.equal(mesh.background,true);assert.equal(mesh.solverReady,false);
  assert.equal(mesh.cells.reduce((n,c)=>n+c.area,0),1);assert.deepEqual(mesh.classificationCounts,[1,1,2]);
  assert.equal(mesh.edges.filter(e=>e.patch===1).length,3);
  assert.throws(()=>requireFluidMesh({background:true,mesh,cm2dPath:'forged.cm2d'}),/完整笛卡尔/);
  assert.throws(()=>requireFluidMesh({mesh,cm2dPath:'forged.cm2d'}),/完整笛卡尔/);
  assert.doesNotThrow(()=>requireFluidMesh({cm2dPath:'actual.cm2d'}));
});
test('background reader rejects holes, overlap, wrong semantics, counts and nonrectangular records',()=>{
  for(const mutate of [x=>x.cells.pop(),x=>x.cells.push(x.cells[0]),x=>x.solver_ready=true,
    x=>x.counts[0]++,x=>x.cells[0].bounds[0]=.01,x=>x.cells[0].classification=3,
    x=>x.boundary_loops[0][0][0]=null,x=>x.cells[0].level=28]){
    const x=fixture();mutate(x);assert.throws(()=>parseBackgroundGrid(JSON.stringify(x)),/背景网格/);
  }
});
test('background invocation bypasses fluid clipping and OpenFOAM, limits dense requests',()=>{
  const {job}=validateJob({method:'background',geometryPath:'body.xy',fluidRegion:'interior',backgroundMode:'uniform',backgroundLevel:9});
  assert.equal(job.fluidRegion,'exterior');
  assert.equal(validateJob({method:'background',geometryPath:'body.xy',backgroundMode:'uniform',backgroundLevel:2}).job.minimumLevel,2);
  const call=buildInvocation(job,{xyPath:'body.xy',prefix:'mesh',casePath:'never-foam'});
  assert.ok(call.args.includes('--background-grid'));assert.ok(!call.args.includes('never-foam'));
  assert.deepEqual(call.cm2dCandidates,[]);assert.equal(call.backgroundPath,'mesh.background.json');
  assert.throws(()=>validateJob({method:'background',geometryPath:'body.xy',backgroundMode:'uniform',backgroundLevel:10}));
  assert.throws(()=>validateJob({method:'background',geometryPath:'body.xy',backgroundMode:'bad'}));
  assert.throws(()=>buildInvocation(job,{}, {dryRun:true}));
  assert.match(exportGuide({background:true,result:{counts:{cells:4}}}),/不是流体求解拓扑/);
});
