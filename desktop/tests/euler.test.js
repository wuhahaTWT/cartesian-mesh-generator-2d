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
  assert.deepEqual(validateEulerRequest(request),{...request,fluxScheme:'rusanov',order:1,wallGradient:'linear',thermalConductivity:0,wallThermal:'insulated',wallValue:0,dynamicViscosity:0,wallModel:'slip'});
  for(const extra of [{pressure:'1'},{density:0},{gamma:1},{cfl:.5},{case:'channel'},{resume:'yes'},{unknown:1},{case:'sod'},{fluxScheme:'roe'},{order:3},{order:'2'}])assert.throws(()=>validateEulerRequest({...request,...extra}));
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

test('Euler numerical method selection is explicit and legacy requests remain reproducible',()=>{
  const selected={...request,fluxScheme:'hllc',order:2};
  const invocation=buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',selected);
  assert.equal(invocation.args[invocation.args.indexOf('--flux')+1],'hllc');
  assert.equal(invocation.args[invocation.args.indexOf('--order')+1],'2');
  assert.throws(()=>validate(fixture.files,selected),/格式不匹配/);
  const summary=JSON.parse(fixture.files['.json']);summary.fluxScheme='hllc';
  assert.throws(()=>validate({...fixture.files,'.json':JSON.stringify(summary)}),/格式不匹配/);
  // The conserved checkpoint binds physics; numerical controls may change explicitly.
  assert.equal(eulerCheckpoint(fixture.files['.checkpoint'],mesh,selected).time,.02);
});

test('Euler validates real second-order HLLC/HLLE output and every fallback/sensor diagnostic',()=>{
  const sample=require('./fixtures/euler-hllc.json'),m=parseCm2d(sample.mesh);
  const check=(files=sample.files)=>validateEulerOutput(JSON.parse(files['.json']),JSON.parse(files['.fields.json']),files['.cells.csv'],files['.faces.csv'],files['.history.csv'],files['.checkpoint'],m,sample.request);
  assert.equal(sample.independentAudit.valid,true);assert.equal(check().summary.order,2);
  for(const [key,value] of [['lastHllcFallbackEvaluations',1],['hllcFallbackEvaluations',1],['minimumContactRestoration',-1],['lastMinimumContactRestoration',1],['shockControl','none']]) {
    const summary=JSON.parse(sample.files['.json']);summary[key]=value;
    assert.throws(()=>check({...sample.files,'.json':JSON.stringify(summary)}),/回退|保护/);
  }
  const lines=sample.files['.faces.csv'].trim().split('\n'),header=lines[0].split(','),row=lines[1].split(',');
  row[header.indexOf('hllcFallbackStages')]='4';lines[1]=row.join(',');
  assert.throws(()=>check({...sample.files,'.faces.csv':lines.join('\n')}),/回退阶段/);
});

test('Fourier total-energy output binds physical walls and audits heat/combined CFL',()=>{
  const sample=require('./fixtures/euler-conduction.json'),m=parseCm2d(sample.mesh),r=sample.request;
  const check=(files=sample.files,request=r)=>validateEulerOutput(JSON.parse(files['.json']),JSON.parse(files['.fields.json']),files['.cells.csv'],files['.faces.csv'],files['.history.csv'],files['.checkpoint'],m,request);
  assert.equal(sample.independentAudit.valid,true);const result=check();assert.ok(result.audit.thermalCourant>0);assert.ok(result.audit.boundaryHeat<0);
  for(const extra of [{thermalConductivity:.6},{wallThermal:'flux'},{wallValue:3}])assert.throws(()=>check(sample.files,{...r,...extra}),/导热|热壁/);
  for(const [suffix,field] of [['.faces.csv','heatFlux'],['.faces.csv','convectiveEnergy'],['.cells.csv','heatRate'],['.history.csv','combinedCourant'],['.history.csv','boundaryHeat']]) {
    const lines=sample.files[suffix].trim().split('\n'),header=lines[0].split(','),values=lines.at(-1).split(',');
    values[header.indexOf(field)]=String(Number(values[header.indexOf(field)])+1);lines[lines.length-1]=values.join(',');
    assert.throws(()=>check({...sample.files,[suffix]:lines.join('\n')}),/导热|能量|热流|CFL|历史/);
  }
  assert.throws(()=>eulerCheckpoint(sample.files['.checkpoint'].replace('CONDUCTIVITY 0.5','CONDUCTIVITY 0.6'),m,r),/导热系数/);
  assert.throws(()=>eulerCheckpoint(sample.files['.checkpoint'],m,{...r,wallValue:2.1}),/热壁/);
  const command=buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',r);
  for(const [flag,value] of [['--conductivity','0.5'],['--wall-thermal','temperature'],['--wall-value','2']])assert.equal(command.args[command.args.indexOf(flag)+1],value);
  for(const extra of [{thermalConductivity:-1},{thermalConductivity:'1'},{thermalConductivity:0},{wallValue:0},{wallThermal:'insulated'},{case:'uniform'},{wallThermal:'radiation'}])assert.throws(()=>validateEulerRequest({...r,...extra}));
  const progress={type:'euler-step',step:1,time:.1,acousticCourant:.1,thermalCourant:.2,combinedCourant:.3,minimumDensity:1,minimumPressure:1,mass:1,totalEnergy:3};
  assert.deepEqual(parseEulerProgress(JSON.stringify(progress)),progress);
  for(const extra of [{thermalCourant:-1},{combinedCourant:.5},{combinedCourant:.05}])assert.throws(()=>parseEulerProgress(JSON.stringify({...progress,...extra})),/CFL/);
});

test('thermal restart manifests prevent transport changes before invoking native solver',async()=>{
  const sample=require('./fixtures/euler-conduction.json'),m=parseCm2d(sample.mesh),r=sample.request;
  const directory=await fs.mkdtemp(path.join(os.tmpdir(),'cm2d-euler-heat-job-'));
  try {
    const meshPath=path.join(directory,'mesh.solver.cm2d');await fs.writeFile(meshPath,sample.mesh);
    const currentResult={cm2dPath:meshPath,outputDirectory:directory};let calls=0;
    const runner=async(_exe,args)=>{calls++;const prefix=args[args.indexOf('--output')+1];await Promise.all(Object.entries(sample.files).map(([suffix,text])=>fs.writeFile(prefix+suffix,text)));return {code:0,stderr:''};};
    const signal=new AbortController().signal;
    const result=await runEulerJob({currentResult,mesh:m,request:r,executable:x=>x,runProcess:runner,signal});
    const imported=await importEulerRestart(path.join(directory,result.manifest),m,meshPath);assert.equal(imported.metadata.request.thermalConductivity,.5);
    for(const extra of [{thermalConductivity:.6},{wallValue:3},{wallThermal:'flux'}])await assert.rejects(()=>runEulerJob({currentResult,mesh:m,request:{...r,resume:true,endTime:.03,...extra},executable:x=>x,runProcess:runner,signal}),/物理参数/);
    assert.equal(calls,1);assert.equal(currentResult.euler,result);
  }finally{await fs.rm(directory,{recursive:true,force:true});}
});

test('viscous output binds dynamic viscosity, no-slip walls and total-energy work',()=>{
  const sample=require('./fixtures/euler-viscosity.json'),m=parseCm2d(sample.mesh),r=sample.request;
  const check=(files=sample.files,request=r)=>validateEulerOutput(JSON.parse(files['.json']),JSON.parse(files['.fields.json']),files['.cells.csv'],files['.faces.csv'],files['.history.csv'],files['.checkpoint'],m,request);
  assert.equal(sample.independentAudit.valid,true);const result=check();assert.ok(result.audit.viscousCourant>0);assert.equal(result.audit.boundaryViscousWork,0);
  for(const extra of [{dynamicViscosity:.3},{wallModel:'slip'}])assert.throws(()=>check(sample.files,{...r,...extra}),/黏度|壁面/);
  for(const extra of [{dynamicViscosity:0},{dynamicViscosity:-1},{dynamicViscosity:'1'},{wallModel:'moving'},{case:'uniform'}])assert.throws(()=>validateEulerRequest({...r,...extra}));
  for(const [suffix,field] of [['.faces.csv','viscousMomentumX'],['.faces.csv','viscousMomentumY'],['.faces.csv','viscousWork'],['.faces.csv','convectiveMomentumX'],['.cells.csv','viscousRate'],['.history.csv','viscousCourant'],['.history.csv','boundaryViscousWork']]) {
    const lines=sample.files[suffix].trim().split('\n'),header=lines[0].split(','),values=lines.at(-1).split(',');values[header.indexOf(field)]=String(field==='viscousRate'?1e9:Number(values[header.indexOf(field)])+1);lines[lines.length-1]=values.join(',');
    assert.throws(()=>check({...sample.files,[suffix]:lines.join('\n')}),/黏性|能量|CFL|历史|无滑移/,field);
  }
  assert.throws(()=>eulerCheckpoint(sample.files['.checkpoint'].replace(/VISCOSITY [^\n]+/,'VISCOSITY 0.3'),m,r),/动力黏度/);
  assert.throws(()=>eulerCheckpoint(sample.files['.checkpoint'],m,{...r,wallModel:'slip'}),/机械壁面/);
  const command=buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',r);assert.equal(command.args[command.args.indexOf('--viscosity')+1],'0.2');assert.equal(command.args[command.args.indexOf('--wall-model')+1],'no-slip');
  const progress={type:'euler-step',step:1,time:.1,acousticCourant:.1,thermalCourant:.1,viscousCourant:.2,combinedCourant:.3,minimumDensity:1,minimumPressure:1,mass:1,totalEnergy:3};
  assert.deepEqual(parseEulerProgress(JSON.stringify(progress)),progress);assert.throws(()=>parseEulerProgress(JSON.stringify({...progress,viscousCourant:.5})),/CFL/);
});

test('viscous restart manifest blocks physical changes before launching native process',async()=>{
  const sample=require('./fixtures/euler-viscosity.json'),m=parseCm2d(sample.mesh),r=sample.request;
  const directory=await fs.mkdtemp(path.join(os.tmpdir(),'cm2d-euler-viscous-job-'));
  try {
    const meshPath=path.join(directory,'mesh.solver.cm2d');await fs.writeFile(meshPath,sample.mesh);
    const currentResult={cm2dPath:meshPath,outputDirectory:directory};let calls=0;
    const runner=async(_exe,args)=>{calls++;const prefix=args[args.indexOf('--output')+1];await Promise.all(Object.entries(sample.files).map(([suffix,text])=>fs.writeFile(prefix+suffix,text)));return {code:0,stderr:''};};
    const signal=new AbortController().signal;const result=await runEulerJob({currentResult,mesh:m,request:r,executable:x=>x,runProcess:runner,signal});
    const imported=await importEulerRestart(path.join(directory,result.manifest),m,meshPath);assert.equal(imported.metadata.request.dynamicViscosity,.2);
    for(const extra of [{dynamicViscosity:.3},{wallModel:'slip'}])await assert.rejects(()=>runEulerJob({currentResult,mesh:m,request:{...r,resume:true,endTime:.01,...extra},executable:x=>x,runProcess:runner,signal}),/物理参数/);
    assert.equal(calls,1);assert.equal(currentResult.euler,result);
  }finally{await fs.rm(directory,{recursive:true,force:true});}
});

test('quadratic wall recovery binds numerical control and actual prescribed walls',()=>{
  const sample=require('./fixtures/euler-wall-accuracy.json'),m=parseCm2d(sample.mesh),r=sample.request;
  const check=(files=sample.files,request=r)=>validateEulerOutput(JSON.parse(files['.json']),JSON.parse(files['.fields.json']),files['.cells.csv'],files['.faces.csv'],files['.history.csv'],files['.checkpoint'],m,request);
  assert.equal(sample.independentAudit.valid,true);const result=check();assert.ok(result.summary.quadraticHeatWalls>0&&result.summary.quadraticViscousWalls>0);
  const invocation=buildEulerInvocation('/tmp/mesh.solver.cm2d','/tmp/run',r);assert.equal(invocation.args[invocation.args.indexOf('--wall-gradient')+1],'quadratic');
  assert.throws(()=>validateEulerRequest({...r,wallGradient:'cubic'}),/梯度格式/);
  assert.throws(()=>check(sample.files,{...r,wallGradient:'linear'}),/梯度格式/);
  for(const [key,value] of [['quadraticHeatWalls',0],['quadraticViscousWalls',0],['quadraticHeatWalls',.5],['quadraticViscousWalls',-1]]){
    const summary=JSON.parse(sample.files['.json']);summary[key]=value;
    assert.throws(()=>check({...sample.files,'.json':JSON.stringify(summary)}),/二次壁面/);
  }
  // Physical states remain restart-compatible when the numerical wall stencil
  // is changed explicitly; the output must still declare the method it used.
  assert.equal(eulerCheckpoint(sample.files['.checkpoint'],m,{...r,wallGradient:'linear'}).time,r.endTime);
});
