'use strict';
// Record small native fields for the App protocol, not a CFD accuracy study.
// node desktop/scripts/record-euler-controls.js outputs/euler-controls
const fs=require('node:fs/promises'),path=require('node:path'),assert=require('node:assert/strict');
const {createHash}=require('node:crypto');
const root=path.resolve(__dirname,'../..');
const {parseCm2d}=require('../src/core/cm2d');
const {runEulerJob,importEulerRestart}=require('../src/core/euler-job');
const {run}=require('../src/core/process');
const {SUFFIXES}=require('../src/core/euler');
const uniform=require('../tests/fixtures/euler-uniform.json'),viscous=require('../tests/fixtures/euler-viscosity.json');
const sha=data=>createHash('sha256').update(data).digest('hex');
async function main() {
  if(process.argv.length!==3)throw new Error('Provide a new outputs directory.');
  const output=path.resolve(process.argv[2]);await fs.mkdir(output,{recursive:false});
  const executable=path.join(root,'build','cartmesh2d_euler_cli'+(process.platform==='win32'?'.exe':''));
  const cases=[
    ['implicit-uniform',uniform,{integrator:'sdirk2',maximumStep:.2,endTime:.2}],
    ['steady-uniform',uniform,{integrator:'sdirk2',maximumStep:.2,endTime:2,mode:'steady',steadyScale:1,steadyTolerance:1e-5}],
    ['total-channel',uniform,{case:'channel',v:0,outletPressure:1,inletModel:'total',inletTotalPressure:1.03,inletTotalTemperature:1.01,integrator:'sdirk2',fluxScheme:'hllc',order:2,maximumStep:.02,endTime:.04}],
    ['implicit-viscous',viscous,{integrator:'sdirk2',maximumStep:.005}],
  ];
  const record={scope:'Native 4/8-cell App protocol fixtures; no physical qualification',platform:process.platform,executableSha256:sha(await fs.readFile(executable)),sourceSha256:{},cases:{},workflow:{}};
  for(const file of ['apps/cartmesh2d_euler_cli.cpp','src/fv/Euler2D.cpp','src/fv/EulerFlux2D.cpp'])record.sourceSha256[file]=sha(await fs.readFile(path.join(root,file)));
  const save=()=>fs.writeFile(path.join(output,'native-run.json'),JSON.stringify(record,null,2)+'\n');
  await save();
  for(const [name,sample,controls] of cases) {
    const directory=path.join(output,name);await fs.mkdir(directory);
    const meshPath=path.join(directory,'mesh.solver.cm2d');await fs.writeFile(meshPath,sample.mesh);
    const mesh=parseCm2d(sample.mesh),currentResult={cm2dPath:meshPath,outputDirectory:directory},request={...sample.request,...controls};
    const commands=[],progress=[];
    record.workflow[name]={commands,request,status:'running'};
    const runner=async(exe,args,...rest)=>{commands.push([exe,...args]);await save();return run(exe,args,...rest);};
    const invoke=(r)=>runEulerJob({currentResult,mesh,request:r,executable:()=>executable,runProcess:runner,signal:new AbortController().signal,onProgress:p=>progress.push(p)});
    const result=await invoke(request),files={};
    for(const suffix of SUFFIXES)files[suffix]=await fs.readFile(path.join(directory,result.files[suffix]),'utf8');
    record.cases[name]={mesh:sample.mesh,request:result.request,files};
    Object.assign(record.workflow[name],{acceptedProgress:progress.length,summary:result.summary,status:'complete'});
    assert.ok(progress.length>0);
    if(name==='implicit-uniform')assert.ok(result.audit.acousticCourant>result.request.cfl);
    if(name==='steady-uniform'){assert.equal(result.summary.steadyConverged,true);assert.ok(result.summary.time<request.endTime);}
    if(name==='total-channel') {
      const original=await importEulerRestart(path.join(directory,result.manifest),mesh,meshPath);assert.equal(original.metadata.request.inletTotalPressure,1.03);
      const limited={...request,endTime:.08,maximumSteps:1};
      await assert.rejects(()=>invoke(limited),/budget/);assert.equal(currentResult.euler,result);
      assert.equal(currentResult.eulerRestart.metadata.steps,1);
      const accepted=currentResult.eulerRestart;
      const resumed=await invoke({...request,endTime:.08,resume:true});
      const resumedCheckpoint=await fs.readFile(currentResult.eulerRestart.path);
      const continuous=await invoke({...request,endTime:.08});
      assert.equal(sha(resumedCheckpoint),sha(await fs.readFile(currentResult.eulerRestart.path)));
      record.workflow[name].restart={retainedAcceptedTime:accepted.metadata.time,resumedTime:resumed.summary.time,continuousTime:continuous.summary.time,checkpointIdentical:true};
    }
    console.log(name,result.summary.status,result.summary.time,result.audit.combinedCourant);
    await save();
  }
  await fs.writeFile(path.join(output,'fixtures.json'),JSON.stringify({scope:record.scope,cases:record.cases},null,2)+'\n');
  await fs.writeFile(path.join(output,'native-run.json'),JSON.stringify(record,null,2)+'\n');
}
main().catch(error=>{console.error(error);process.exitCode=1;});
