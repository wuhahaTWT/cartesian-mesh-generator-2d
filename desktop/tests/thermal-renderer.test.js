'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../src/renderer/app.js'), 'utf8').replace(/\r\n/g, '\n');
const functionSource = name => {
  const start = source.indexOf(`function ${name}(`);
  const next = source.indexOf('\nfunction ', start + 1);
  return source.slice(start, next);
};
function controls() {
  const elements = Object.fromEntries(Object.entries({flowMode:'transient', flowCase:'external', flowNu:'0.01', flowSpeed:'1', flowConvection:'upwind'})
    .map(([id, value]) => [id, {value}]));
  elements.flowResume = {checked:true}; elements.thermalResume = {checked:false};
  elements.flowInitialVortex = {checked:false};
  const counts = { flow:0, thermal:0 };
  const context = {
    $:id => elements[id], state:{flowRestart:{case:'external', nu:.01, speed:1, convection:'upwind'}},
    clearFlowBinding:() => counts.flow++, clearThermalBinding:() => counts.thermal++,
    updateFlowScope() {}, updateFlowMode() {}, setThermalEvents(events) { elements.events=events; }
  };
  vm.createContext(context);
  vm.runInContext(functionSource('applySharedFlowControls') + '\n' + functionSource('applyRestartControls'), context);
  return {context, elements, counts};
}
test('steady flow mode deselects its hidden resume without disrupting thermal resume', () => {
  const {context, elements, counts} = controls();
  elements.flowMode.value = 'steady'; elements.thermalResume.checked = true;
  context.state.flowRestart.nu = .2;
  context.applyRestartControls();
  assert.equal(elements.flowResume.checked, false);
  assert.equal(elements.thermalResume.checked, true);
  assert.equal(elements.flowNu.value, '0.01');
  assert.deepEqual(counts, {flow:0, thermal:0});
});
test('flow resume replaces shared physics and invalidates incompatible displayed fields', () => {
  const {context, elements, counts} = controls();
  elements.thermalResume.checked = true; context.state.flowRestart.nu = .2;
  context.applyRestartControls();
  assert.equal(elements.thermalResume.checked, false);
  assert.equal(elements.flowNu.value, .2);
  assert.deepEqual(counts, {flow:1, thermal:1});
  context.applyRestartControls();
  assert.deepEqual(counts, {flow:1, thermal:1}, 'unchanged physics preserves current results');
});
test('thermal export uses saved backend fields when renderer has been cleared', async () => {
  for (const lineEnding of ['\n', '\r\n']) {
    let call;
    const context = { document:{fonts:{load:async()=>{},ready:Promise.resolve()}}, window:{
      cartmesh:{exportPreviewData:async()=>({mesh:{id:'saved'},thermal:{id:'saved-temperature'}})},
      CartMeshExport:{renderThermal:(mesh,thermal)=>{call={mesh,thermal};return 'png';}}
    }};
    vm.createContext(context);
    // Git may check out CRLF on Windows; delimit by lines, not host bytes.
    const input = source.replace(/\n/g, lineEnding).replace(/\r\n/g, '\n');
    const start = input.indexOf('window.__exportThermalPreview =');
    const end = input.indexOf('\n\nlet previewSequence', start);
    assert.ok(start >= 0 && end > start, 'preview function bounds found');
    vm.runInContext(input.slice(start, end), context);
    assert.equal(await context.window.__exportThermalPreview(), 'png');
    assert.equal(call.mesh.id, 'saved');
    assert.equal(call.thermal.id, 'saved-temperature');
    context.window.cartmesh.exportPreviewData = async()=>({mesh:{},thermal:null});
    assert.equal(await context.window.__exportThermalPreview(), null);
  }
});

test('thermal resume deselects flow resume, restores thermal physics and preserves editable time controls', () => {
  for(const method of ['bounded','bounded-spatial']) {
    const {context, elements, counts} = controls();
    for (const id of ['thermalDiffusivity','thermalInitial','thermalSource','thermalConvection','thermalFluxCorrection','thermalWallKind','thermalWallValue','thermalWallInflow']) elements[id] = {value:''};
    elements.flowDt = {value:'0.005'}; elements.flowSteps = {value:'7'};
    elements.thermalResume.checked = true;
    context.thermalPatches = ['wall']; context.thermalPatchId = ()=>'thermalWall';
    context.state.thermalRestart = {request:{case:'external',nu:.02,speed:1,convection:'upwind',diffusivity:.1,initial:300,source:2,scalarConvection:'limited-linear',fluxCorrection:method,boundaries:{wall:{kind:'value',value:350,inflowValue:300}}}};
    vm.runInContext(functionSource('applyThermalRestartControls'), context);
    context.applyThermalRestartControls();
    assert.equal(elements.flowResume.checked, false);
    assert.equal(elements.thermalResume.checked, true);
    assert.equal(elements.flowNu.value, .02);
    assert.equal(elements.thermalDiffusivity.value, .1);
    assert.equal(elements.thermalWallValue.value, 350);
    assert.equal(elements.thermalFluxCorrection.value, method);
    assert.equal(elements.flowDt.value, '0.005');
    assert.equal(elements.flowSteps.value, '7');
    assert.deepEqual(counts, {flow:1,thermal:1});
  }
});

test('loading a thermal restart clears a locked vortex and preserves supported adaptive controls', () => {
  const {context, elements} = controls();
  for (const id of ['thermalDiffusivity','thermalInitial','thermalSource','thermalConvection','thermalFluxCorrection']) elements[id] = {value:''};
  elements.flowMode.value = 'adaptive';
  elements.flowInitialVortex.checked = true;
  elements.thermalResume.checked = true;
  elements.flowDt = {value:'0.005'}; elements.flowSteps = {value:'7'};
  context.thermalPatches = [];
  context.state.thermalRestart = {request:{case:'external',nu:.01,speed:1,convection:'upwind',
    diffusivity:.1,initial:300,source:2,scalarConvection:'upwind',boundaries:{},events:[{time:.137,target:'source',kind:'source',value:2}]}};
  vm.runInContext(functionSource('applyThermalRestartControls'), context);
  context.applyThermalRestartControls();
  assert.equal(elements.thermalFluxCorrection.value, 'unrestricted', 'legacy restart keeps its original flux operator');
  assert.equal(elements.flowInitialVortex.checked, false, 'restart uses its saved flow, without a fresh vortex');
  assert.equal(elements.flowMode.value, 'adaptive');
  assert.deepEqual(elements.events,context.state.thermalRestart.request.events);
  assert.equal(elements.flowResume.checked, false);
  assert.equal(elements.thermalResume.checked, true);
  assert.equal(elements.flowDt.value, '0.005');
  assert.equal(elements.flowSteps.value, '7');
});


test('thermal failure status exposes controller cause and preserved state while keeping file paths in the log', async () => {
  const detail = '热计算未完成；已接受到 t=0 s。\n联合时间缺陷为 2.88 倍预算，接受上限为 1。\n最后尝试步长为 0.000001 s，已到最小步长；可减小最小步长并保持当前精度预算。';
  for (const retained of [false, true]) {
    const message = detail + (retained ? '\n已保留 t=0.4 s 的联合续算状态。' : '\n尚无已接受的联合续算状态。') + '\n诊断文件保留在 /example/thermal-run。';
    const elements = {runThermal:{}, thermalResume:{checked:false}};
    const saved = retained ? {thermal:{summary:{time:.3}}, restart:{time:.4}} : {};
    const statuses = [], logs = [];
    const context = {
      $:id=>elements[id], state:{result:{},mesh:{},thermalHistory:[]},
      validThermalInputs:()=>true, thermalRequest:()=>({resume:false}),
      clearThermalBinding(){}, renderThermalMonitor(){}, setBusy(){}, applyThermalRestartControls(){},
      status:(title,text,expanded)=>statuses.push({title,text,expanded}), log:text=>logs.push(text),
      refreshThermalState:async restore=>{assert.equal(restore,true);return saved;},
      window:{cartmesh:{runThermal:async()=>{throw new Error("Error invoking remote method 'thermal:run': Error: " + message);}}}
    };
    vm.createContext(context);
    const start = source.indexOf('async function runThermal(');
    const end = source.indexOf('\nfunction boundedFlowHistory(', start);
    assert.ok(start>=0 && end>start);
    vm.runInContext(source.slice(start,end),context);
    await context.runThermal();
    const status = statuses.at(-1);
    assert.equal(status.title,'温度推进未完成');
    assert.equal(status.expanded,true,'failure explanations must wrap instead of using an ellipsis');
    assert.match(status.text,/时间缺陷为 2\.88 倍预算/);
    assert.match(status.text,/最小步长/);
    assert.match(status.text,retained ? /已保留 t=0\.4 s/ : /尚无已接受/);
    assert.equal(elements.thermalResume.checked,retained);
    assert.equal(status.text.includes('/example/thermal-run'),false);
    if (retained) assert.match(status.text,/当前显示上次完整温度结果/);
    assert.equal(logs.at(-1),message,'full diagnostic path remains available in the log');
  }
});


test('fresh thermal startup allows the existing vortex control and joint resume omits it', () => {
  const elements=new Proxy({}, {get:(target,key)=>target[key] ||= {value:'0',checked:false,hidden:false}});
  elements.flowMode.value='adaptive';elements.flowCase.value='external';elements.flowInitialVortex.checked=true;
  elements.flowVortexX.value='3';elements.flowVortexY.value='0';elements.flowVortexRadius.value='1';elements.flowVortexSpeed.value='.05';
  elements.flowDt.value='.05';elements.flowSteps.value='2';elements.thermalTemperatureScale.value='1';elements.thermalTimeRtol.value='.01';
  const context={$:id=>elements[id],state:{result:{},busy:false},document:{querySelectorAll:()=>[]},
    updateThermalEventEditor(){},thermalPatches:[],thermalEventRows:()=>[]};
  vm.createContext(context);vm.runInContext(functionSource('updateThermalMode')+'\n'+functionSource('thermalRequest'),context);
  context.updateThermalMode();assert.equal(elements.runThermal.disabled,false);
  assert.match(elements.thermalTimeHint.textContent,/仅在零时刻施加/);
  assert.equal(context.thermalRequest().initialVortex.peakSpeed,.05);
  context.state.thermalRestart={time:.5};elements.thermalResume.checked=true;
  context.updateThermalMode();assert.equal(context.thermalRequest().initialVortex,undefined);
  assert.doesNotMatch(elements.thermalTimeHint.textContent,/仅在零时刻施加/);
});
