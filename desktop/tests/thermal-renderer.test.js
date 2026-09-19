'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../src/renderer/app.js'), 'utf8');
const functionSource = name => {
  const start = source.indexOf(`function ${name}(`);
  const next = source.indexOf('\nfunction ', start + 1);
  return source.slice(start, next);
};
function controls() {
  const elements = Object.fromEntries(Object.entries({flowMode:'transient', flowCase:'external', flowNu:'0.01', flowSpeed:'1', flowConvection:'upwind'})
    .map(([id, value]) => [id, {value}]));
  elements.flowResume = {checked:true}; elements.thermalResume = {checked:false};
  const counts = { flow:0, thermal:0 };
  const context = {
    $:id => elements[id], state:{flowRestart:{case:'external', nu:.01, speed:1, convection:'upwind'}},
    clearFlowBinding:() => counts.flow++, clearThermalBinding:() => counts.thermal++,
    updateFlowScope() {}, updateFlowMode() {}
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
  let call;
  const context = { document:{fonts:{load:async()=>{},ready:Promise.resolve()}}, window:{
    cartmesh:{exportPreviewData:async()=>({mesh:{id:'saved'},thermal:{id:'saved-temperature'}})},
    CartMeshExport:{renderThermal:(mesh,thermal)=>{call={mesh,thermal};return 'png';}}
  }};
  vm.createContext(context);
  const start = source.indexOf('window.__exportThermalPreview =');
  vm.runInContext(source.slice(start, source.indexOf('\n\nlet previewSequence', start)), context);
  assert.equal(await context.window.__exportThermalPreview(), 'png');
  assert.equal(call.mesh.id, 'saved');
  assert.equal(call.thermal.id, 'saved-temperature');
  context.window.cartmesh.exportPreviewData = async()=>({mesh:{},thermal:null});
  assert.equal(await context.window.__exportThermalPreview(), null);
});

test('thermal resume deselects flow resume, restores thermal physics and preserves editable time controls', () => {
  const {context, elements, counts} = controls();
  for (const id of ['thermalDiffusivity','thermalInitial','thermalSource','thermalConvection','thermalWallKind','thermalWallValue','thermalWallInflow']) elements[id] = {value:''};
  elements.flowDt = {value:'0.005'}; elements.flowSteps = {value:'7'};
  elements.thermalResume.checked = true;
  context.thermalPatches = ['wall']; context.thermalPatchId = ()=>'thermalWall';
  context.state.thermalRestart = {request:{case:'external',nu:.02,speed:1,convection:'upwind',diffusivity:.1,initial:300,source:2,scalarConvection:'limited-linear',boundaries:{wall:{kind:'value',value:350,inflowValue:300}}}};
  vm.runInContext(functionSource('applyThermalRestartControls'), context);
  context.applyThermalRestartControls();
  assert.equal(elements.flowResume.checked, false);
  assert.equal(elements.thermalResume.checked, true);
  assert.equal(elements.flowNu.value, .02);
  assert.equal(elements.thermalDiffusivity.value, .1);
  assert.equal(elements.thermalWallValue.value, 350);
  assert.equal(elements.flowDt.value, '0.005');
  assert.equal(elements.flowSteps.value, '7');
  assert.deepEqual(counts, {flow:1,thermal:1});
});
