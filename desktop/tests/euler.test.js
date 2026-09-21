'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const fs=require('node:fs/promises'),path=require('node:path'),os=require('node:os');
const fixture=require('./fixtures/euler-uniform.json');
const {parseCm2d}=require('../src/core/cm2d');
const {validateEulerRequest,buildEulerInvocation,validateEulerOutput,parseEulerProgress,eulerCheckpoint}=require('../src/core/euler');
const {runEulerJob,importEulerRestart}=require('../src/core/euler-job');
const mesh=parseCm2d(fixture.mesh),request=fixture.request;
function validate(files=fixture.files,override=request){return validateEulerOutput(JSON.parse(files['.json']),JSON.parse(files['.fields.json']),files['.cells.csv'],files['.faces.csv'],files['.history.csv'],files['.checkpoint'],mesh,override);}
test('Euler accepts real independently audited native fields and rejects stale/tampered output',()=>{
  assert.equal(fixture.independentAudit.valid,true);const result=validate();assert.equal(result.summary.time,.02);assert.ok(result.audit.maximumCellBalanceRelative<1e-12);
  assert.throws(()=>validate(fixture.files,{...request,pressure:2}),/请求不同/);
  for(const [key,value] of [['p',-1],['rho',2],['rhoE',10],['temperature',12],['mach',5]]) {
    const fields=JSON.parse(fixture.files['.fields.json']);fields.cells[0][key]=value;
    assert.throws(()=>validate({...fixture.files,'.fields.json':JSON.stringify(fields)}),/显示/);
  }
  const fields={...fixture.files};let rows=fields['.faces.csv'].trim().split('\n'),header=rows[0].split(','),values=rows[1].split(',');values[header.indexOf('energy')]=String(Number(values[header.indexOf('energy')])+1);rows[1]=values.join(',');fields['.faces.csv']=rows.join('\n');
  assert.throws(()=>validate(fields),/守恒/);
  assert.throws(()=>validate({...fixture.files,'.checkpoint':fixture.files['.checkpoint']+'extra\n'}),/截断|额外/);
  assert.throws(()=>validate({...fixture.files,'.history.csv':fixture.files['.history.csv'].trim().split('\n').slice(0,-1).join('\n')}),/历史/);
});
test('Euler requests and progress keep physical types and units separate',()=>{
  assert.deepEqual(validateEulerRequest(request),request);
  for(const extra of [{pressure:'1'},{density:0},{gamma:1},{cfl:.5},{case:'channel'},{resume:'yes'},{unknown:1},{case:'sod'}])assert.throws(()=>validateEulerRequest({...request,...extra}));
  const command=buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',request);assert.equal(command.executable,'cartmesh2d_euler_cli');assert.ok(command.args.includes('--gas-r'));
  assert.throws(()=>buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',{...request,resume:true}),/状态/);
  assert.equal(parseEulerProgress('normal output'),null);
  const p={type:'euler-step',step:25,time:.1,acousticCourant:.3,minimumDensity:.9,minimumPressure:.8,mass:1,totalEnergy:3};assert.deepEqual(parseEulerProgress(JSON.stringify(p)),p);
  assert.throws(()=>parseEulerProgress(JSON.stringify({...p,minimumPressure:0})),/进度/);
  assert.equal(eulerCheckpoint(fixture.files['.checkpoint'],mesh,request).time,.02);
});
test('Euler job retains prior complete result on failure and binds imported restart manifests',async()=>{
  const directory=await fs.mkdtemp(path.join(os.tmpdir(),'cm2d-euler-job-'));
  try {
    const meshPath=path.join(directory,'mesh.solver.cm2d');await fs.writeFile(meshPath,fixture.mesh);
    const currentResult={cm2dPath:meshPath,outputDirectory:directory};const signal=new AbortController().signal;
    const runner=async(_exe,args)=>{const prefix=args[args.indexOf('--output')+1];await Promise.all(Object.entries(fixture.files).map(([suffix,text])=>fs.writeFile(prefix+suffix,text)));return {code:0,stderr:''};};
    const payload=await runEulerJob({currentResult,mesh,request,executable:x=>x,runProcess:runner,signal});assert.equal(currentResult.euler,payload);
    const imported=await importEulerRestart(path.join(directory,payload.manifest),mesh,meshPath);assert.equal(imported.metadata.time,.02);
    const fail=async(exe,args)=>{await runner(exe,args);return {code:2,stderr:'deliberate budget failure'};};
    await assert.rejects(()=>runEulerJob({currentResult,mesh,request,executable:x=>x,runProcess:fail,signal}),/budget/);
    assert.equal(currentResult.euler,payload);assert.equal(currentResult.eulerRestart.metadata.time,.02);
    await assert.rejects(()=>runEulerJob({currentResult,mesh,request:{...request,resume:true,endTime:.03,gamma:1.3},executable:x=>x,runProcess:runner,signal}),/物理参数/);
    await fs.appendFile(imported.path,'corruption\n');await assert.rejects(()=>importEulerRestart(path.join(directory,payload.manifest),mesh,meshPath),/检查点/);
  }finally{await fs.rm(directory,{recursive:true,force:true});}
});
