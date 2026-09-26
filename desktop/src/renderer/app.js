'use strict';
document.documentElement.dataset.platform = window.cartmesh.platform;

const $ = id => document.getElementById(id);
const state = {
  catalog: null,
  method: 'cutcell',
  geometryPath: '',
  geometryLabel: '',
  outputDirectory: '',
  mesh: null,
  wallBounds: null,
  result: null,
  flow: null,
  flowRestart: null,
  flowBoundaryDefinition: null,
  flowHistory: [],
  euler: null,
  eulerRestart: null,
  eulerHistory: [],
  thermal: null,
  thermalRestart: null,
  thermalHistory: [],
  // Hand-placed refinement regions, in body spans about the body centre.
  regions: [],
  frame: null
};

const disabledControls = new Map();
const SIDEBAR_STORAGE_KEY = 'cartmesh2d-sidebar-collapsed';
function readSidebarCollapsed() {
  try { return window.localStorage.getItem(SIDEBAR_STORAGE_KEY) === 'true'; } catch (_) { return false; }
}
function setSidebarCollapsed(collapsed, persist = true) {
  document.body.classList.toggle('sidebar-collapsed', collapsed);
  const button = $('toggleSidebar');
  if (button) {
    button.setAttribute('aria-expanded', String(!collapsed));
    button.setAttribute('aria-label', collapsed ? '显示设置侧栏' : '隐藏设置侧栏');
    button.title = collapsed ? '显示设置侧栏（⌘B / Ctrl+B）' : '隐藏设置侧栏（⌘B / Ctrl+B）';
  }
  if (persist) {
    try { window.localStorage.setItem(SIDEBAR_STORAGE_KEY, String(collapsed)); } catch (_) { /* private mode */ }
  }
  window.dispatchEvent(new Event('resize'));
}
function initSidebar() {
  setSidebarCollapsed(readSidebarCollapsed(), false);
  $('toggleSidebar')?.addEventListener('click', () => {
    setSidebarCollapsed(!document.body.classList.contains('sidebar-collapsed'));
  });
  window.addEventListener('keydown', event => {
    if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === 'b' &&
        !event.altKey && !event.shiftKey && !['INPUT', 'TEXTAREA', 'SELECT'].includes(document.activeElement?.tagName)) {
      event.preventDefault();
      setSidebarCollapsed(!document.body.classList.contains('sidebar-collapsed'));
    }
  });
}
initSidebar();
function setBusy(busy) {
  state.busy = busy;
  if (busy) {
    for (const control of document.querySelectorAll('.panel input, .panel select, .panel button')) {
      disabledControls.set(control, control.disabled);
      control.disabled = true;
    }
  } else {
    for (const [control, disabled] of disabledControls) control.disabled = disabled;
    disabledControls.clear();
  }
  $('cancel').hidden = !busy;
  $('cancel').disabled = false;
  $('exportResult').disabled = busy;
  $('returnHome').disabled = busy;
  $('actualToManual').disabled = busy;
  $('runFlow').disabled = busy || !state.result;
  updateFlowMode();
  updateEulerControls();
  updateReady();
}
function validInputs() {
  for (const input of document.querySelectorAll('.panel input[type=number]')) {
    if (!input.closest('#flowBlock, #thermalBlock, #eulerBlock') && !input.disabled && input.getClientRects().length &&
        (!input.value.trim() || !input.checkValidity())) {
      input.reportValidity();
      input.focus();
      status('参数需要调整', '请填写有效数值，并检查范围。');
      return false;
    }
  }
  return true;
}
function validFlowInputs() {
  for (const input of document.querySelectorAll('#flowBlock input[type=number]')) {
    if (input.dataset.varying === 'true' && !input.value.trim()) continue;
    if (!input.disabled && input.getClientRects().length && (!input.value.trim() || !input.checkValidity())) {
      input.reportValidity(); input.focus();
      status('工况参数需要调整', '请填写有效数值，并检查范围。');
      return false;
    }
  }
  return true;
}
function clearFlowBinding({ hidePanel = false } = {}) {
  state.flow = null;
  state.flowHistory=[];
  $('flowTimeline').hidden=true;
  if (hidePanel) {state.flowRestart=null;state.flowBoundaryDefinition=null;view.setBoundaryHighlight([]);renderFlowBoundaries();updateFlowMode();
    $('flowCaseFileInfo').textContent='保存从零起算的流动设置，读取前须生成同一最终网格；续算请使用重启文件。温度设置单独配置。';}
  view.setFlowFields(null);
  $('flowResult').hidden = true;
  $('flowResult').replaceChildren();
  $('flowSpeedOption').hidden = true;
  $('flowPressureOption').hidden = true;
  if (['speed', 'pressure'].includes($('displayMode').value)) {
    $('displayMode').value = 'level'; view.mode = 'level';
    $('canvasWrap').classList.remove('light');
    view.draw();
  }
  if (hidePanel) $('flowBlock').hidden = true;
  if (state.mesh) renderLegend(state.mesh, state.levelBasis);
}
function clearResult() {
  clearEulerBinding(true);
  clearFlowBinding({ hidePanel: true });
  clearThermalBinding({ hidePanel: true });
  state.mesh = null; state.result = null; state.wallBounds = null;
  state.job = null;
  state.selectedRequest = null; state.cellBudget = null;
  $('cellBudgetResult').hidden = true; $('actualToManual').hidden = true;
  view.clear();
  $('exportResult').hidden = true;
  for (const id of ['counters', 'gates', 'histogram']) $(id).replaceChildren();
  $('legend').hidden = true;
  $('resolutionResult').hidden = true;
  $('resolutionResult').replaceChildren();
}
async function returnToStart() {
  if (state.busy) return;
  ++previewSequence;
  clearInterval(progressTimer);
  const exportable = !$('exportResult').hidden;
  clearResult();
  state.geometryPath = ''; state.geometryLabel = ''; state.sampleId = null;
  state.geometryLoading = false; state.frame = null; state.regions = [];
  view.regions = []; view.frame = null;
  $('sample').value = '';
  $('sourceUnits').disabled = false;
  for (const id of ['sampleNote', 'geometryFacts', 'probeResult', 'verifiedPreset', 'verifiedPresetNote']) $(id).hidden = true;
  $('log').textContent = '';
  $('empty').hidden = false;
  $('exportResult').hidden = !exportable;
  renderRegions();
  updateReady();
  status('从一个轮廓开始', exportable ? '预览已释放；上一份结果仍可导出，生成新网格后替换。' : '选择样例或导入几何。');
  try { await window.cartmesh.releasePreview(); }
  catch (error) { log(error.message); }
}
$('returnHome').addEventListener('click', returnToStart);
window.__exportMeshPreview = async () => {
  await document.fonts.load('13px "CartMesh UI"', '生成网格');
  await document.fonts.ready;
  const data = await window.cartmesh.exportPreviewData();
  return window.CartMeshExport.render(data.mesh, data.result);
};

window.__exportThermalPreview = async () => {
  await document.fonts.load('13px "CartMesh UI"', '温度输运');
  await document.fonts.ready;
  const payload = await window.cartmesh.exportPreviewData();
  return payload.thermal ? window.CartMeshExport.renderThermal(payload.mesh, payload.thermal) : null;
};

window.__exportEulerPreview = async () => {
  await document.fonts.load('13px "CartMesh UI"','可压结果');await document.fonts.ready;
  const payload=await window.cartmesh.exportPreviewData();
  return payload.euler ? window.CartMeshExport.renderEuler(payload.mesh,payload.euler) : null;
};

let previewSequence = 0;
const view = new window.MeshView.Viewport($('canvas'));
const { levelColour, RAMP, SPEED_RAMP, PRESSURE_RAMP, TEMPERATURE_RAMP } = window.MeshView;

const fmt = value => Number(value || 0).toLocaleString('en-US');
const log = line => { $('log').textContent += `${line}\n`; $('log').scrollTop = 1e9; };
const status = (title, text) => { $('statusTitle').textContent = title; $('statusText').textContent = text || ''; };

function importSettings() {
  return {
    chordError: Number($('chordError').value),
    sourceUnits: $('sourceUnits').value,
    fluidRegion: $('fluidRegion').value
  };
}

// One place builds the request both `probe-sizing` and `generate` take, so a probe
// can never describe a different job from the one that runs.
function buildRequest() {
  if(state.method==='background') return {method:'background',geometryPath:state.geometryPath,
    outputDirectory:state.outputDirectory,automatic:false,fluidRegion:'exterior',
    backgroundMode:$('backgroundMode').value,backgroundLevel:Number($('backgroundLevel').value),
    backgroundMinimumLevel:Number($('backgroundMinimumLevel').value),backgroundPadding:Number($('backgroundPadding').value),...importSettings()};
  const base = {
    automatic: $('controlMode').value === 'auto',
    density: $('density').value,
    targetCells: $('controlMode').value === 'auto' && !['normal','dense'].includes($('density').value)
      ? Number($('density').value === 'custom' ? $('customTargetCells').value : $('density').value) : undefined,
    wallToBackgroundRatio: Number($('autoWallRatio').value),
    method: state.method,
    geometryPath: state.geometryPath,
    outputDirectory: state.outputDirectory,
    sizingMode: 'relative',
    referenceLength: $('referenceMode').value === 'explicit' ? Number($('referenceLength').value) : undefined,
    wallRelativeSize: Number($('wallRelativeSize').value),
    backgroundRelativeSize: Number($('backgroundRelativeSize').value),
    farFieldSpans: Number($('controlMode').value === 'auto' ? $('autoPadding').value : $('relativePadding').value),
    cellsPerLevel: Number($('controlMode').value === 'auto' ? $('autoBandCells').value : $('relativeBandCells').value),
    allowUnsafeWallLevel: $('relativeAllowUnsafe').checked,
    ...importSettings()
  };
  if (state.method === 'hybrid') {
    return {
      ...base,
      maxLevel: Number($('maxLevel').value),
      minimumLevel: Number($('minimumLevel').value),
      boundaryLevel: Number($('boundaryLevel').value),
      nLayers: Number($('controlMode').value === 'auto' ? $('autoLayers').value : $('nLayers').value),
      firstLayerRelativeSize: Number($('firstLayerRelativeSize').value),
      extrusionRelativeSize: Number($('extrusionRelativeSize').value),
      firstThickness: Number($('firstThickness').value),
      growthRatio: Number($('controlMode').value === 'auto' ? $('autoGrowth').value : $('growthRatio').value),
      domainPadding: Number($('domainPadding').value),
      extrusionThickness: Number($('extrusionThickness').value)
    };
  }
  return {
    ...base,
    smallAlpha: Number($('smallAlpha').value),
    wallCellsPerSpan: Number($('wallCellsPerSpan').value),
    farLevel: Number($('farLevel').value),
    curvatureCellsPerRadius: !base.targetCells && $('useCurvature').checked ? Number($('curvatureCellsPerRadius').value) : 0,
    gapCells: !base.targetCells && $('useGap').checked ? Number($('gapCells').value) : 0,
    wake: !base.targetCells && $('useWake').checked ? {
      angleOfAttackDeg: Number($('wakeAngle').value),
      downstreamSpans: Number($('wakeLength').value),
      halfWidthSpans: Number($('wakeHalfWidth').value),
      levelsBelowWall: Number($('wakeLevels').value)
    } : null,
    refineBoxes: (base.targetCells ? [] : state.regions).map(({ xmin, xmax, ymin, ymax, levelsBelowWall }) =>
      ({ xmin, xmax, ymin, ymax, levelsBelowWall }))
  };
}

// Level = ceil(log2((1 + 2*far) * wallCells)) — the body span cancels, so the whole
// far-field-versus-wall-resolution trade is arithmetic and can be shown live.
function updateBudget() {
  if(state.method==='background') {
    const uniform=$('backgroundMode').value==='uniform';
    $('backgroundLevel').max=uniform?'9':'10';
    $('backgroundMinimumLevel').disabled=uniform||state.busy;
    $('generate').disabled=!state.geometryPath||state.busy||state.geometryLoading;
    return;
  }
  const method = state.catalog.methods[state.method];
  if (!method) return;
  const far = Number($('relativePadding').value);
  const wall = Number($('wallRelativeSize').value);
  const background = Number($('backgroundRelativeSize').value);
  const reference = $('referenceMode').value === 'explicit' ? Number($('referenceLength').value) : state.frame?.bodySpan;
  const span = state.frame?.bodySpan;
  const ceiling = method.safeWallLevel;
  const domain = span + 2*far*reference;
  const level = Math.max(0, Math.ceil(Math.log2(domain/(wall*reference))));
  const over = level > ceiling && !$('relativeAllowUnsafe').checked;
  const invalid = !(wall>0) || !(background>=wall) || !(reference>0);
  $('resolutionBudget').className = `budget${over || invalid ? ' over' : ''}`;
  $('resolutionBudget').textContent = !span ? '选择几何后显示物理尺寸和构造深度。'
    : `参考长度 ${reference.toPrecision(5)} m；壁面目标 ${(wall*reference).toPrecision(5)} m（${(100*wall).toPrecision(4)}% Lref）。` +
      `派生树深 ${level}，当前已验证上限 ${ceiling}。` +
      (over ? ' 此请求超过已验证深度，需要明确选择是否尝试。' : '') +
      (background<wall ? ' 背景尺寸不能小于壁面尺寸。' : '');
  $('generate').disabled = (!state.geometryPath || state.busy || state.geometryLoading) || ($('controlMode').value === 'manual' && (over || invalid));
}

function updateReady() {
  updateBudget();
  schedulePlan();
}

let planTimer;
let planSequence = 0;
function schedulePlan() {
  clearTimeout(planTimer);
  const sequence = ++planSequence;
  state.budgetSuggestion = null;
  $('planToManual').disabled = true;
  if (state.method==='background' || !state.frame || $('controlMode').value !== 'auto') return;
  const request = buildRequest();
  if (request.targetCells == null) {
    $('autoBudgetPreview').textContent = '旧版预设保留兼容行为，不按目标数量调整。';
    return;
  }
  planTimer = setTimeout(async () => {
    try {
      const plan = await window.cartmesh.planBudget({ request, frame: state.frame });
      if (sequence !== planSequence) return;
      state.budgetSuggestion = plan;
      const ref = plan.referenceLength ?? state.frame.bodySpan;
      const depth = Math.ceil(Math.log2((state.frame.bodySpan+2*plan.farFieldSpans*ref)/(plan.wallRelativeSize*ref)));
      $('autoBudgetPreview').textContent = `目标 ${fmt(request.targetCells)} 格，接受范围 ${fmt(Math.round(request.targetCells*.7))}–${fmt(Math.round(request.targetCells*1.3))}。` +
        `初始壁面 ${(100*plan.wallRelativeSize).toPrecision(3)}% Lref，背景 ${(100*plan.backgroundRelativeSize).toPrecision(3)}% Lref。` +
        (plan.budgetPlan.areaEstimated ? '多环面积采用包围盒初估，生成后按实际数量校正。' : '') +
        `预计构造深度 ${depth}。` + (plan.budgetPlan.limitedBySafeWallLevel ? '当前受到默认深度限制，可能达不到数量目标；可勾选下方允许更高构造深度。' : depth > state.catalog.methods[state.method].safeWallLevel && !request.allowUnsafeWallLevel ? '需勾选下方“允许更高构造深度”才可尝试。' : '生成后按实测数量调整。');
      $('planToManual').disabled = Boolean(state.busy);
    } catch (error) {
      if (sequence === planSequence) $('autoBudgetPreview').textContent = error.message.replace(/^Error invoking remote method.*?Error: /, '');
    }
  }, 160);
}
function useManualParameters(request) {
  if (!request || state.busy) return;
  const fields = { wallRelativeSize:'wallRelativeSize', backgroundRelativeSize:'backgroundRelativeSize',
    farFieldSpans:'relativePadding', cellsPerLevel:'relativeBandCells', firstLayerRelativeSize:'firstLayerRelativeSize',
    nLayers:'nLayers', growthRatio:'growthRatio', extrusionRelativeSize:'extrusionRelativeSize', smallAlpha:'smallAlpha' };
  for (const [key,id] of Object.entries(fields)) if (request[key] != null) $(id).value = request[key];
  if (request.referenceLength != null) {
    $('referenceMode').value = 'explicit'; $('referenceLength').value = request.referenceLength;
    $('referenceLengthField').hidden = false;
    $('referenceMode').closest('details').open = true;
  }
  $('relativeAllowUnsafe').checked = Boolean(request.allowUnsafeWallLevel);
  $('controlMode').value = 'manual'; updateControlMode(); syncRegions();
  status('已转为手动设置', '参数已填入，可微调后重新生成。');
}
$('planToManual').addEventListener('click', () => useManualParameters(state.budgetSuggestion));
$('actualToManual').addEventListener('click', () => useManualParameters(state.selectedRequest));
for (const id of ['density','customTargetCells','autoPadding','autoWallRatio','autoBandCells','autoLayers','autoGrowth']) {
  $(id).addEventListener('input', () => { $('customTargetField').hidden = $('density').value !== 'custom'; schedulePlan(); });
  $(id).addEventListener('change', schedulePlan);
}
for (const [id, factor] of [['manualCoarser',Math.SQRT2],['manualFiner',1/Math.SQRT2]]) {
  $(id).addEventListener('click', () => {
    for (const field of ['wallRelativeSize','backgroundRelativeSize','firstLayerRelativeSize'])
      $(field).value = Math.min(field === 'backgroundRelativeSize' ? 1e6 : 1,Number($(field).value)*factor).toPrecision(6);
    updateBudget();
  });
}

async function loadVerifiedPreset() {
  if (state.busy || state.geometryLoading) return;
  const reference = { circle:2, naca2412:1, nozzle:6 }[state.sampleId];
  if (!reference) return;
  $('sourceUnits').value = 'm';
  $('fluidRegion').value = state.sampleId === 'nozzle' ? 'interior' : 'exterior';
  const wall = state.sampleId === 'nozzle' ? .00125 : .0025;
  useManualParameters({ referenceLength:reference, wallRelativeSize:wall,
    backgroundRelativeSize:.0125, farFieldSpans:.5, cellsPerLevel:40,
    firstLayerRelativeSize:wall/4, nLayers:4, growthRatio:1.2,
    extrusionRelativeSize:.02/reference, smallAlpha:.1, allowUnsafeWallLevel:true });
  $('useCurvature').checked = false; $('useGap').checked = false; $('useWake').checked = false;
  $('curvatureCellsPerRadius').disabled = true; $('gapCells').disabled = true;
  $('wakeFields').hidden = true; state.regions = []; renderRegions();
  await drawGeometryOutline();
  status('已载入高密样例参数', '参考长度、计算域与高深度选项已按已验证案例替换；检查后点击生成。');
}
$('verifiedPreset').addEventListener('click', loadVerifiedPreset);

function renderMethods() {
  const container = $('methods');
  container.innerHTML = '';
  for (const method of Object.values(state.catalog.methods)) {
    const button = document.createElement('button');
    button.className = `method${method.id === state.method ? ' active' : ''}`;
    button.innerHTML =
      `<b>${method.label}<span class="tag ${method.status}">${method.statusLabel}</span></b>`;
    button.title = method.summary;
    button.addEventListener('click', () => selectMethod(method.id));
    container.appendChild(button);
  }
}

function selectMethod(id) {
  if (id !== state.method) { if (state.flow) clearFlowBinding(); if (state.thermal) clearThermalBinding(); }
  state.method = id;
  const method = state.catalog.methods[id];
  $('sizingBlock').hidden = !method.supports.sizeField;
  $('hybridBlock').hidden = method.supports.sizeField;
  $('methodNote').textContent = id==='background' ? method.summary : method.supports.sizeField
    ? `实测安全壁面层级上限 ${method.safeWallLevel}；越过要显式勾选。支持 OpenFOAM 导出。`
    : `贴体层和余域共用相对尺寸；构造深度由程序推导，当前验证上限 ${method.safeWallLevel}。`;
  $('smallAlphaField').hidden = !method.supports.sizeField;
  $('smallAlphaNote').hidden = !method.supports.sizeField;
  renderMethods();
  updateControlMode();
  updateReady();
}

function updateControlMode() {
  const background=state.method==='background';
  const automatic = !background && $('controlMode').value === 'auto';
  $('backgroundBlock').hidden=!background;
  $('fluidRegion').closest('label').hidden=background;
  $('controlMode').closest('label').hidden=background;
  $('verifiedPreset').hidden=background||!['circle','naca2412','nozzle'].includes(state.sampleId);
  $('verifiedPresetNote').hidden=$('verifiedPreset').hidden;
  $('controlMode').disabled=background||state.busy;
  $('resolutionBlock').hidden = automatic||background;
  $('sizingBlock').hidden = automatic || state.method !== 'cutcell';
  $('hybridBlock').hidden = automatic || state.method !== 'hybrid';
  $('smallAlphaField').hidden = automatic || state.method !== 'cutcell';
  $('smallAlphaNote').hidden = automatic || state.method !== 'cutcell';
  $('automaticControls').hidden = !automatic;
  $('autoLayerFields').hidden = state.method !== 'hybrid';
  $('customTargetField').hidden = $('density').value !== 'custom';
  $('autoNote').hidden = !automatic;
  updateReady();
}
function fitOverview() {
  if (state.wallBounds) view.fitTo(state.wallBounds, 0.25);
  else if (state.mesh) view.fitTo(state.mesh.bounds, 0.1);
}

function applySizeField(field) {
  if (!field) return;
  if (field.farFieldSpans !== undefined) $('farFieldSpans').value = field.farFieldSpans;
  if (field.wallCellsPerSpan !== undefined) $('wallCellsPerSpan').value = field.wallCellsPerSpan;
  if (field.cellsPerLevel !== undefined) $('cellsPerLevel').value = field.cellsPerLevel;
  if (field.wallCellsPerSpan !== undefined) $('wallRelativeSize').value = 1/field.wallCellsPerSpan;
  if (field.farFieldSpans !== undefined) {
    $('relativePadding').value = field.farFieldSpans;
    // Loading a geometry/preset defines its physical domain in both modes.
    // Leaving auto mode at the generic 0.5 default silently truncated the
    // external-flow domain, even when the sample prescribed ten body spans.
    $('autoPadding').value = field.farFieldSpans;
    $('backgroundRelativeSize').value = (1+2*field.farFieldSpans)/Math.pow(2,field.farLevel ?? 0);
  }
  if (field.cellsPerLevel !== undefined) $('relativeBandCells').value = field.cellsPerLevel;
  $('farLevel').value = field.farLevel ?? 0;
  updateBudget();
}

function renderPresets() {
  const container = $('presets');
  container.innerHTML = '';
  for (const preset of state.catalog.presets) {
    const button = document.createElement('button');
    button.textContent = preset.label;
    button.title = preset.hint;
    button.addEventListener('click', () => {
      [...container.children].forEach(child => child.classList.remove('active'));
      button.classList.add('active');
      applySizeField(preset.sizeField);
      status('已套用预设', preset.hint);
    });
    container.appendChild(button);
  }
}

function renderSamples() {
  const select = $('sample');
  for (const sample of state.catalog.samples) {
    const option = document.createElement('option');
    option.value = sample.id;
    option.textContent = sample.label;
    select.appendChild(option);
  }
}

async function chooseGeometry(path, label, sample) {
  $('sourceUnits').disabled = false;
  state.geometryPath = path;
  state.geometryLabel = label;
  state.sampleId = sample?.id;
  $('verifiedPreset').hidden = state.method==='background'||!['circle','naca2412','nozzle'].includes(state.sampleId);
  $('verifiedPresetNote').hidden = $('verifiedPreset').hidden;
  if (sample) {
    $('sourceUnits').value = 'm';
    $('fluidRegion').value = sample.fluidRegion;
    applySizeField(sample.fluidRegion === 'interior' && sample.interiorSizeField || sample.sizeField);
    // Each sample ships the small-cell threshold it was verified at.  The passing
    // combinations are not monotone in alpha: near the level ceiling, whether a wall
    // vertex grazes a grid line is what decides the solver gate.
    if (sample.smallAlpha !== undefined) $('smallAlpha').value = sample.fluidRegion === 'interior'
      ? sample.interiorSmallAlpha ?? sample.smallAlpha : sample.smallAlpha;
    $('useGap').checked = Boolean(sample.gapCells);
    $('gapCells').disabled = !sample.gapCells;
    if (sample.gapCells) $('gapCells').value = sample.gapCells;
    $('useWake').checked = Boolean(sample.wake);
    $('wakeFields').hidden = !sample.wake;
    if (sample.wake) {
      $('wakeAngle').value = sample.wake.angleOfAttackDeg;
      $('wakeLength').value = sample.wake.downstreamSpans;
      $('wakeHalfWidth').value = sample.wake.halfWidthSpans;
      $('wakeLevels').value = sample.wake.levelsBelowWall;
    }
    $('sampleNote').textContent = sample.note;
    $('sampleNote').hidden = false;
  } else {
    $('sampleNote').hidden = true;
  }
  $('probeResult').hidden = true;
  await drawGeometryOutline();
  updateReady();
}

// Draw the imported boundary before meshing.  This is how the user notices that an
// SVG came in mirrored or that a CSV lost half its loops, while it is still cheap.
async function drawGeometryOutline() {
  const sequence = ++previewSequence;
  clearResult();
  state.geometryLoading = true;
  updateReady();
  try {
    status('正在读取几何', state.geometryLabel);
    const preview = await window.cartmesh.previewGeometry({
      geometryPath: state.geometryPath, ...importSettings()
    });
    if (sequence !== previewSequence) return;
    $('sourceUnits').disabled = preview.kind === 'raster';
    state.frame = preview.frame;
    view.setOutline(preview.loops);
    syncRegions();
    $('empty').hidden = true;
    $('legend').hidden = true;
    const spanX = preview.frame.width;
    const spanY = preview.frame.height;
    $('geometryFacts').hidden = false;
    $('geometryFacts').innerHTML =
      `<span>格式 <b>${preview.kind}</b>　环 <b>${preview.loops.length}</b>　` +
      `顶点 <b>${preview.loops.reduce((sum, loop) => sum + loop.length, 0)}</b></span>` +
      `<span>包围盒 <b>${spanX.toPrecision(4)} × ${spanY.toPrecision(4)}</b></span>`;
    preview.warnings.forEach(warning => log(`注意：${warning}`));
    status('几何就绪', state.geometryLabel);
  } catch (error) {
    if (sequence !== previewSequence) return;
    state.geometryPath = '';
    $('geometryFacts').hidden = true;
    status('几何读取失败', error.message.split('\n')[0]);
    log(error.message);
  } finally {
    if (sequence === previewSequence) { state.geometryLoading = false; updateReady(); }
  }
}

async function probeSizing() {
  if (state.busy || state.geometryLoading || !state.geometryPath || !validInputs()) return;
  setBusy(true);
  try {
    const probe = await window.cartmesh.probeSizing({ ...buildRequest(), outputDirectory: '.' });
    const field = probe.field;
    const rows = [
      ['参考长度（m）', field?.reference_length?.toPrecision(5) ?? '—'],
      ['计算域跨度（m）', field?.domain_span?.toPrecision(5) ?? '—'],
      ['请求壁面尺寸（m）', field?.requested_wall_size?.toPrecision(5) ?? '—'],
      ['取整后背景壁面尺寸（m）', field?.wall_cell_size?.toPrecision(5) ?? '—'],
      ['派生壁面 / 背景层级', `${field?.wall_level ?? '—'} / ${field?.minimum_level ?? '—'}`]
    ];
    $('relativeProbeResult').innerHTML = rows
      .map(([key, value]) => `<span><i class="k">${key}</i> ${value}</span>`).join('') +
      (probe.ok ? '<span class="ok">可行</span>'
                : (probe.issues || []).map(issue => `<span class="bad">${issue}</span>`).join(''));
    $('relativeProbeResult').hidden = false;
    status(probe.ok ? '尺寸场可行' : '尺寸场被拒绝',
      probe.ok ? '尺寸解析通过；不代表最终网格或质量验收通过。' : probe.message);
  } catch (error) {
    $('relativeProbeResult').textContent = error.message;
    $('relativeProbeResult').hidden = false;
  } finally {
    setBusy(false);
  }
}

function renderCounters(result) {
  if(result.background){
    $('counters').innerHTML=[['完整单元',result.counts.cells],['外部',result.counts.classification[0]],['内部',result.counts.classification[1]],['相交（未裁切）',result.counts.classification[2]]]
      .map(([label,n])=>`<div><span>${label}</span><b>${fmt(n)}</b></div>`).join('');return;
  }
  // A run that died before its summary block leaves some counters unknown.  Showing a
  // dash beats showing a confident zero.
  const count = value => (value ? fmt(value) : '—');
  const rows = [
    ['solver cells', count(result.counts.cells)],
    ['vertices', count(result.counts.vertices)],
    ['faces', count(result.counts.faces)],
    ['quadtree leaves', count(result.counts.leaves)],
    ['cut cells', count(result.counts.cutCells)],
    ['openfoam cells', result.openFoam.written ? fmt(result.openFoam.cells) : '未导出']
  ];
  if (result.counts.layerCells) rows.splice(4, 0, ['layer cells', fmt(result.counts.layerCells)]);
  $('counters').innerHTML = rows
    .map(([label, value]) => `<div><span>${label}</span><b>${value}</b></div>`).join('');
}

function renderGates(result) {
  if(result.background){$('gates').innerHTML=gateRow('背景网格','已生成','完整域覆盖与记录一致性已核对。物体内部保留；几何轮廓不是数值壁面，未生成流体求解拓扑或OpenFOAM算例。');return;}
  const parts = [];
  parts.push(gateRow('内部拓扑检查', result.gates.topology.pass === null ? '未确认' : result.gates.topology.pass ? 'PASS' : 'FAIL',
    result.gates.topology.pass ? '生成器内部检查通过；外部 checkMesh 需另外执行。' : '本次未完整成功，不能据此确认通过。'));
  const solver = result.gates.solver;
  if (solver) {
    const worst = solver.rows.filter(row => !row.pass);
    parts.push(gateRow('Solver 质量', solver.valid ? 'PASS' : 'FAIL',
      worst.length ? worst.map(row => `${row.label} ${row.value.toPrecision(4)} (限 ${row.limit})`).join('　') :
        solver.rows.map(row => `${row.label} ${row.value.toPrecision(4)}`).join('　')));
  }
  const directional = result.resolution?.directional_connectivity;
  if (directional) {
    const minimum = Number.isFinite(directional.minimum) ? directional.minimum.toPrecision(5) : '未测得';
    parts.push(gateRow('方向连通质量', directional.valid ? 'PASS' : 'FAIL',
      `最小 determinant ${minimum}，下限 ${directional.threshold}；不合格单元 ${directional.failed_cell_ids.length}。` +
      '适用于平面均匀挤出、前后 empty；外部 checkMesh 需另外执行。'));
  }
  if (result.actualMethod === 'cutcell-fallback')
    parts.push(gateRow('方法', 'WARN', '请求了贴体边界层，实际落到纯 Cut-cell fallback'));
  if (!result.openFoam.written)
    parts.push(gateRow('OpenFOAM 导出', 'WARN', '未写出。solver 门未通过，或导出阶段自身失败；CM2D 与 VTK 仍然可用'));
  $('gates').innerHTML = parts.join('');
}

const gateRow = (label, verdict, detail) =>
  `<div class="gate"><b>${label}</b><span class="verdict ${verdict}">${verdict}</span>` +
  `<details class="gate-detail"><summary>详情</summary><span class="detail">${detail}</span></details></div>`;

// Cells per level is the readout that answers the far-field question directly: it
// shows how few cells the coarse levels actually cost.
//
// Built with DOM calls rather than an HTML string because the strict CSP has no
// style-src 'unsafe-inline': a `style="..."` attribute would be dropped, while
// assigning element.style here is allowed.
function renderHistogram(histogram, mesh, basis) {
  const container = $('histogram');
  container.replaceChildren();
  if (!histogram.length) return;
  const title = document.createElement('div');
  title.className = 'title';
  // The hybrid mesh has no Quadtree level to report, so its bands are cell sizes.
  title.textContent = basis === 'size' ? '每尺寸档单元数（粗→细）' : '每层级单元数';
  container.appendChild(title);
  const peak = Math.max(...histogram.map(row => row.count));
  for (const row of histogram) {
    const bar = document.createElement('div');
    bar.className = 'bar';
    const label = document.createElement('em');
    label.textContent = `${basis === 'size' ? '档' : 'L'}${row.level}`;
    const fill = document.createElement('i');
    fill.style.width = `${Math.max(1, Math.round((row.count / peak) * 108))}px`;
    fill.style.background = levelColour(row.level, mesh.minLevel, mesh.maxLevel);
    const count = document.createElement('b');
    count.textContent = fmt(row.count);
    bar.append(label, fill, count);
    container.appendChild(bar);
  }
}

function renderLegend(mesh, basis) {
  const container = $('legend');
  container.replaceChildren();
  if(mesh.background && view.mode==='classification'){
    for(const [i,label] of ['外部','内部','相交'].entries()){
      const entry=document.createElement('span');entry.textContent=`${label} ${fmt(mesh.classificationCounts[i])}　`;
      entry.style.color=['#79b8d4','#9da9bc','#efac56'][i];container.appendChild(entry);
    }
    container.hidden=false;return;
  }
  const fieldMode = ['speed', 'pressure', 'temperature'].includes(view.mode) || view.mode.startsWith('euler-');
  const coloured = view.mode === 'level' || fieldMode;
  const ramp = document.createElement('div');
  ramp.className = 'ramp';
  const palette = fieldMode ? view.fieldPalette() : RAMP;
  for (const colour of palette) {
    const swatch = document.createElement('span');
    swatch.style.background = colour;
    ramp.appendChild(swatch);
  }
  const ends = document.createElement('div');
  ends.className = 'ends';
  const coarse = document.createElement('span');
  const fine = document.createElement('span');
  if (fieldMode) {
    const range = view.fieldRange;
    const unit = ({'euler-rho':'kg/m³','euler-p':'Pa','euler-temperature':'K','euler-mach':'Mach','euler-speed':'m/s'})[view.mode] || (view.mode === 'temperature' ? 'K' : view.mode === 'speed' ? 'm/s' : 'm²/s²');
    coarse.textContent = range ? `${range.min.toPrecision(4)} ${unit}` : '最小';
    fine.textContent = range ? `${range.max.toPrecision(4)} ${unit}` : '最大';
  } else {
    const prefix = basis === 'size' ? '档' : 'L';
    coarse.textContent = `${prefix}${mesh.minLevel} 粗`;
    fine.textContent = `细 ${prefix}${mesh.maxLevel}`;
  }
  ends.append(coarse, fine);

  const keys = document.createElement('div');
  keys.className = 'keys';
  const theme = view.theme();
  for (const [colour, text] of [[theme.wall, '物面（嵌入边界）'],
                                [theme.domain, '计算域边界']]) {
    const key = document.createElement('span');
    const dash = document.createElement('i');
    dash.style.background = colour;
    key.append(dash, document.createTextNode(text));
    keys.appendChild(key);
  }
  if (coloured) container.append(ramp, ends);
  container.appendChild(keys);
  container.hidden = false;
}

// Both of these failures are alignment events: a wall vertex grazing a grid line.
// Measured on a 48-vertex CSV outline, the solver gate fails at a 8x and 9x far field
// and passes at 10x and 11x with the wall resolution untouched, so nudging the far
// field is a real lever and not folklore.
function advice(message) {
  if (/topology audit|UnclassifiedBoundaryEdge/.test(message)) {
    return '提示：这是壁面顶点擦过格点的构造失败（见 CURRENT_STATE_CN.md 第 4 节 A）。' +
      '把远场距离改动 ±1 倍体长会移动整个格点阵列，通常比调壁面分辨率更有效。';
  }
  if (/solver-quality gate failed/.test(message)) {
    return '提示：把小单元阈值 α 提高一档先合并掉薄片，或把远场距离改动 ±1 倍体长以改变格点对齐。' +
      '两者都不改变你要求的壁面分辨率。';
  }
  return '';
}

async function generate() {
  if (state.busy || state.geometryLoading || !state.geometryPath || !validInputs()) return;
  clearResult();
  setBusy(true);
  $('generate').textContent = '正在生成…';
  $('log').textContent = '';
  status('生成中', state.method==='background'?'几何诊断 → 完整笛卡尔网格 → 分类 → 导出':'几何转换 → 尺寸场 → 加密 → cut-cell → 稳定化 → 质量');
  try {
    const payload = await window.cartmesh.generate(buildRequest());
    state.mesh = payload.mesh;
    state.wallBounds = payload.wallBounds;
    state.result = payload.result;
    state.job = payload.job;
    state.levelBasis = payload.levelBasis;
    state.selectedRequest = payload.selectedRequest || null;
    state.cellBudget = payload.cellBudget || null;
    $('actualToManual').hidden = !state.selectedRequest;
    $('cellBudgetResult').hidden = !payload.cellBudget;
    if (payload.cellBudget) {
      const b = payload.cellBudget;
      $('cellBudgetResult').textContent = `目标约 ${fmt(b.targetCells)}；实际 ${fmt(b.actualCells)} 个单元。${b.reached ? '已达到目标范围（±30%）。' : '未达到目标范围，保留本次最接近目标的可导出结果。'} 共尝试 ${b.attempts} 次。` + (b.stoppedReason ? ` 后续调整停止：${b.stoppedReason}` : '');
    }
    $('classificationOption').hidden=!payload.background;
    if(payload.background){$('displayMode').value='classification';view.mode='classification';}
    else if(view.mode==='classification'){$('displayMode').value='level';view.mode='level';}
    view.setMesh(payload.mesh);
    syncRegions();
    $('empty').hidden = true;
    renderLegend(payload.mesh, payload.levelBasis);
    renderCounters(payload.result);
    renderGates(payload.result);
    renderHistogram(payload.levelHistogram, payload.mesh, payload.levelBasis);
    const resolution = payload.result.resolution;
    $('resolutionResult').hidden = !resolution;
    if (resolution) {
      const actual = resolution.actual;
      const percent = resolution.wall_owner_tangential_exceedance_length_fraction;
      $('resolutionResult').textContent = `最终求解网格实测：Lref=${resolution.reference_length.toPrecision(5)} m；` +
        `单元等效尺寸 h/Lref 中位 ${actual.sqrt_area_over_reference?.p50.toPrecision(4) ?? '—'}；` +
        `壁面单元切向跨度 P95 ${actual.wall_owner_tangential_extent_over_reference?.p95.toPrecision(4) ?? '—'}；` +
        `切向跨度超过壁面目标的壁长占比 ${percent == null ? '未设置目标' : (100*percent).toFixed(2)+'%'}。` +
        '切向跨度与法向高度分开测量；完整尺寸报告随结果包保存。';
      const layers = resolution.boundary_layers;
      if (layers && layers.status !== 'not_requested') {
        const percent = value => value == null ? '未测得' : (100*value).toFixed(2)+'%';
        $('resolutionResult').textContent += layers.status === 'no_layers_retained'
          ? ' 本次为纯网格回退，未保留边界层。'
          : ` 边界层最终保留 ${layers.retained_cells}/${layers.requested_cells} 个请求单元；` +
            `首层覆盖壁长 ${percent(layers.first_layer_wall_length_fraction)}，` +
            `完整请求层数覆盖壁长 ${percent(layers.full_requested_layers_wall_length_fraction)}。` +
            `首层法向高度超过目标的壁长占比 ${percent(layers.first_layer_height_exceedance_wall_length_fraction)}。` +
            (layers.status === 'incomplete' ? ' 层身份或连通性核对不完整，详见报告。' : '');
      }
    }
    fitOverview();
    if (payload.automatic) {
      const job = payload.job;
      $('autoNote').textContent = job.relativeSizing
        ? `${payload.cellBudget ? '参数已按最终选中的结果记录。' : payload.densityReduced ? '更密请求未满足，已降密。' : ''}上次生成的壁面目标 h/Lref=${job.sizeField.wallRelativeSize.toPrecision(4)}，背景 h/Lref=${job.sizeField.backgroundRelativeSize.toPrecision(4)}，留白 ${job.sizeField.farFieldSpans.toPrecision(4)} Lref。实际参数与尝试记录随结果包保存。`
        : `${payload.densityReduced ? '更密参数未通过，已降至可生成的密度。' : ''}本次采用：${job.method === 'cutcell'
        ? `壁面目标体长/${job.sizeField.wallCellsPerSpan}，每级带宽 ${job.sizeField.cellsPerLevel} 格，远场 ${job.sizeField.farFieldSpans} 倍，α ${job.smallAlpha}`
        : `余域 / 壁面 level ${job.maxLevel} / ${job.boundaryLevel}，${job.nLayers} 层，首层 ${job.firstThickness}`}。实际参数与尝试记录随结果包保存。`;
    }
    $('exportResult').hidden = Boolean(payload.incomplete);
    $('flowBlock').hidden = Boolean(payload.incomplete || payload.background);
    $('thermalBlock').hidden = Boolean(payload.incomplete || payload.background);
    $('eulerBlock').hidden = Boolean(payload.incomplete || payload.background);
    if (!payload.incomplete && !payload.background) {
      $('flowCase').value = payload.job.fluidRegion === 'interior' ? 'duct' : 'external';
      $('eulerCase').value=payload.job.fluidRegion==='interior'?'sod':'external';
      updateEulerControls();
      updateFlowScope();
    }
    const seconds = payload.result.timings.total_seconds;
    if (payload.incomplete) {
      status('网格已生成，后续步骤失败', payload.incomplete);
    } else {
      status(payload.cellBudget && !payload.cellBudget.reached ? '网格已生成，数量未达目标' : '生成完成', `${payload.cellBudget ? `已检查 ${payload.attempts.length} 组参数` : payload.automatic ? `自动选参成功（第 ${payload.attempts.length} 组）` : '手动生成成功'}；可导出结果包${seconds ? `　${seconds.toFixed(2)} s` : ''}`);
    }
  } catch (error) {
    status('生成失败', error.message.replace(/^Error invoking remote method '[^']+': Error: /, '').split('\n')[0]);
    log(error.message);
    const hint = advice(error.message);
    if (hint) log(hint);
  } finally {
    clearInterval(progressTimer);
    $('generate').textContent = '生成预览';
    setBusy(false);
  }
}

function updateFlowScope() {
  const selected = state.catalog?.flowCases?.[$('flowCase').value];
  const custom=$('flowCase').value==='custom';
  const expectedRegion = $('flowCase').value === 'external' ? 'exterior' : 'interior';
  const mismatch = !custom && state.job && state.job.fluidRegion !== expectedRegion
    ? ` 当前最终网格是${state.job.fluidRegion === 'interior' ? '内流' : '外流'}语义，与此工况不匹配。`
    : '';
  $('flowScope').textContent = (selected?.scope || '') + mismatch;
  $('flowBoundarySettings').hidden=!custom;
  if (custom) $('flowOutletBackflow').value='reject';
  else view.setBoundaryHighlight([]);
  updateFlowMode();
}

function renderFlowBoundaries() {
  const container=$('flowBoundaryPatches');container.replaceChildren();
  const definition=state.flowBoundaryDefinition;
  $('flowBoundaryInfo').textContent=definition
    ? `${definition.records.length} 个边界面，绑定 ${definition.cells.toLocaleString()} 个单元。点击“显示位置”核对边界。`
    : '边界绑定本次最终网格。重新生成网格后需要重新配置。';
  if (!definition) return;
  const groups=new Map();
  for(const b of definition.records){if(!groups.has(b.name))groups.set(b.name,[]);groups.get(b.name).push(b);}
  const edited=()=>{clearFlowBinding();clearThermalBinding();};
  for(const [name,records] of groups){
    const card=document.createElement('details');card.className='flow-boundary-patch';card.open=records[0].type!=='wall';
    const heading=document.createElement('summary');heading.textContent=`${name} · ${records.length} 个面`;card.append(heading);
    const label=document.createElement('label');label.className='field';label.textContent='边界名称';
    const nameInput=document.createElement('input');nameInput.value=name;nameInput.maxLength=128;
    nameInput.addEventListener('change',()=>{records.forEach(b=>b.name=nameInput.value);edited();renderFlowBoundaries();updateFlowMode();});
    label.append(nameInput);card.append(label);
    const type=document.createElement('select');type.setAttribute('aria-label',`${name} 边界类型`);
    for(const [value,text] of [['velocity-inlet','速度入口'],['pressure-outlet','压力出口'],['pressure-opening','定压开口（法向进出）'],['symmetry','对称／自由滑移'],['wall','静止壁面'],['moving-wall','移动壁面（逐面恒速）'],['smooth-moving-wall','移动壁面（平滑变化）']]){
      const option=document.createElement('option');option.value=value;option.textContent=text;type.append(option);
    }
    type.value=records[0].type;
    type.addEventListener('change',()=>{
      for(const b of records){b.type=type.value;if(!['pressure-outlet','pressure-opening'].includes(type.value))b.p=0;if(['wall','symmetry','pressure-outlet','pressure-opening'].includes(type.value))b.u=b.v=0;}
      edited();renderFlowBoundaries();updateFlowMode();
    });card.append(type);
    if(type.value==='pressure-opening'){const hint=document.createElement('p');hint.className='hint';hint.textContent='仅支持水平或竖直开口。压力差决定流量；入流切向速度为零，法向速度自由。输入静压除以密度，不是总压。';card.append(hint);}
    if(type.value==='symmetry'){const hint=document.createElement('p');hint.className='hint';hint.textContent='仅支持水平或竖直边界。不允许穿过边界，切向无摩擦；对称面没有无滑移壁面的阻力。';card.append(hint);}
    const numeric=['pressure-outlet','pressure-opening'].includes(type.value) ? [['p','静态运动学压力 p / ρ（m²/s²）']] : ['wall','symmetry'].includes(type.value) ? [] : [['u','速度 x'],['v','速度 y']];
    for(const [key,caption] of numeric){
      const field=document.createElement('label');field.className='field';field.textContent=caption;
      const input=document.createElement('input');input.type='number';input.step='any';input.dataset.boundaryKey=key;input.dataset.patchName=name;
      const uniform=records.every(b=>b[key]===records[0][key]);input.value=uniform?String(records[0][key]):'';
      input.dataset.varying=String(!uniform);input.placeholder=uniform?'':'保持逐面分布';
      input.addEventListener('input',()=>{if(!input.value.trim()&&input.dataset.varying==='true')return;
        records.forEach(b=>b[key]=Number(input.value));input.dataset.varying='false';edited();});
      field.append(input);card.append(field);
    }
    const show=document.createElement('button');show.textContent='显示位置';show.type='button';
    show.addEventListener('click',()=>view.setBoundaryHighlight(records.map(b=>b.face)));card.append(show);
    container.append(card);
  }
}
async function prepareFlowBoundaries(source) {
  if(state.busy||!state.result)return;
  setBusy(true);
  try{
    const definition=await window.cartmesh.prepareFlowBoundaries({source,speed:Number($('flowSpeed').value)});
    if(definition){if(source!=='import')$('flowBoundaryPreset').value=source;state.flowBoundaryDefinition=definition;clearFlowBinding();clearThermalBinding();renderFlowBoundaries();
      status('命名边界已载入','按名称编辑条件，点击“显示位置”核对，然后启动计算。');}
  }catch(error){status('边界配置未载入',error.message);}
  finally{setBusy(false);}
}

function flowConvectionLabel(summary) {
  const scheme = state.catalog?.flowConvectionSchemes?.[summary.convection];
  const label = scheme?.label || (summary.convection === 'limited-linear'
    ? '线性迎风（限制重构）' : '一阶迎风');
  return `${label}${summary.convectionInferred ? '（旧摘要缺字段，按一阶迎风推断）' : ''}`;
}

function flowPressurePreconditionerLabel(summary) {
  if (summary.pressurePreconditionerInferred || summary.pressurePreconditioner === 'legacy-unspecified')
    return '旧结果未记录';
  const scheme = state.catalog?.flowPressurePreconditioners?.[summary.pressurePreconditioner];
  const label = scheme?.label || (summary.pressurePreconditioner === 'aggregation'
    ? '多重网格（试验）' : '标准（IC0）');
  return label;
}
function flowOutletBackflowLabel(summary) {
  if (summary.outletBackflow === 'normal-inlet') return '允许法向回流（试验）';
  return '检测到回流时停止';
}

function renderFlowResult(summary) {
  const container = $('flowResult');
  container.replaceChildren();
  const convectionText = flowConvectionLabel(summary);
  const pressurePreconditionerText = flowPressurePreconditionerLabel(summary);
  const pressureText = summary.pressureDiscretization === 'shared-face-gauss'
    ? '共享面压力' : '旧结果未记录';
  const stateLine = document.createElement('div');
  stateLine.className = 'flow-state';
  stateLine.textContent = summary.converged
    ? `已收敛 · ${summary.iterations} 次迭代 · 对流：${convectionText}`
    : `到达 ${summary.iterations} 次迭代上限，结果有效但未收敛 · 对流：${convectionText}`;
  if (summary.temporalDiscretization) stateLine.textContent=`已接受 t=${summary.acceptedTime.toPrecision(6)} s · 本次 ${summary.completedSteps} 步 · 最后一步内迭代 ${summary.iterations} 次`;
  container.appendChild(stateLine);
  const rows = [
    ...(summary.temporalDiscretization ? [['最后一步步长（s）',summary.dt],['最后一步最大 CFL',summary.maxCourant]] : []),
    ...(summary.timeStepControl==='adaptive-cfl-retry' ? [['时间步控制','自动 CFL / 重试'],['拒绝的试算次数',String(summary.rejectedSteps)],['CFL 上限',summary.targetCourant]] : []),
    ['对流格式', convectionText],
    ['压力求解', pressurePreconditionerText],
    ['线性迭代精度',summary.linearPolicy==='adaptive'?(summary.converged?'自适应 · 严格复核通过':'自适应 · 尚未收敛'):'固定精度'],
    ['速度松弛系数',summary.velocityRelaxation ?? .6],
    ['压力校正次数',summary.pressureCorrectionPasses ?? 4],
    ['稳态加速', summary.steadyAcceleration==='anderson'
      ? `历史迭代 · 接受 ${summary.accelerationAccepted} / 舍弃 ${summary.accelerationRejected}` : '关闭'],
    ['出口回流', flowOutletBackflowLabel(summary)],
    ['压力离散', pressureText],
    ['局部连续性（无量纲）', summary.continuity],
    ['全局不平衡（m²/s）', summary.globalImbalance],
    ['全局相对不平衡', summary.globalRelativeImbalance],
    ['速度变化', summary.velocityChange],
    ['流动停止容差', summary.tolerance],
    ['压力变化', summary.pressureChange],
    ['动量残差', summary.momentumResidual]
  ];
  if (summary.outletBackflowFaces !== undefined)
    rows.push(['回流出口面数', String(summary.outletBackflowFaces)]);
  if (summary.outletInflow !== undefined)
    rows.push(['出口流入量（m²/s）', summary.outletInflow]);
  if (summary.namedBoundaryFluxes) {
    rows.push(['边界体积流量', '每单位厚度；净流出为正，净流入为负']);
    for(const flux of summary.namedBoundaryFluxes) {
      rows.push([`${flux.name} · 净流量（m²/s）`,flux.net],
        [`${flux.name} · 流入／流出（m²/s）`,`${flux.inflow.toExponential(3)} / ${flux.outflow.toExponential(3)}`],
        [`${flux.name} · 平均法向速度（m/s）`,flux.normalMeanVelocity]);
    }
  }
  if (summary.namedWallLoads) {
    rows.push(['壁面载荷', '流体作用于壁面；力÷ρ÷单位厚度'],
      ['力矩基准', '原点 (0,0)，逆时针为正']);
    for (const load of summary.namedWallLoads) {
      rows.push([`${load.name} · Fx（m³/s²）`,load.forceX],
        [`${load.name} · Fy（m³/s²）`,load.forceY],
        [`${load.name} · 力矩（m⁴/s²）`,load.torque]);
    }
  }
  for (const [label, value] of rows) {
    const item = document.createElement('div');
    const caption = document.createElement('span'); caption.textContent = label;
    const number = document.createElement('b');
    number.textContent = typeof value === 'number' ? value.toExponential(3) : value;
    item.append(caption, number); container.appendChild(item);
  }
  container.hidden = false;
}

function updateFlowMode() {
  const transient = $('flowMode').value !== 'steady';
  const adaptive = $('flowMode').value === 'adaptive';
  $('flowTimeSettings').hidden = !transient && $('flowCase').value === 'custom';
  $('flowAdaptiveSettings').hidden=!adaptive;
  $('flowStepsField').hidden=adaptive;
  $('flowSteps').disabled=state.busy||adaptive;
  $('flowDtLabel').textContent=adaptive?'最大时间步长':transient?'时间步长':'温度联算步长';
  for (const id of ['flowEndTime','flowMaxCourant','flowMinDt','flowMaxRetries','flowMaxSteps']) $(id).disabled=state.busy||!adaptive;
  $('flowRestartSettings').hidden = !transient;
  $('flowIterationLabel').textContent = transient ? '每个时间步的内迭代上限' : '最大 SIMPLE 迭代';
  const restart = state.flowRestart;
  if (!restart) $('flowResume').checked = false;
  $('flowResume').disabled = state.busy || !restart || !transient;
  const resuming = transient && $('flowResume').checked && restart;
  const thermalResuming = $('thermalResume').checked && state.thermalRestart;
  $('saveFlowCase').disabled = Boolean(state.busy || !state.result || resuming || thermalResuming);
  $('saveFlowCase').title = resuming || thermalResuming ? '请先取消续算，再保存从零起算的设置。' : '';
  $('loadFlowCase').disabled = state.busy || !state.result;
  $('flowInitialSettings').hidden=!transient;
  $('flowInitialVortex').disabled=Boolean(state.busy||!transient||resuming||thermalResuming);
  for(const id of ['flowVortexX','flowVortexY','flowVortexRadius','flowVortexSpeed'])
    $(id).disabled=Boolean(state.busy||!transient||resuming||thermalResuming||!$('flowInitialVortex').checked);
  for (const id of ['flowCase','flowNu','flowSpeed','flowConvection','flowOutletBackflow']) {
    const element=$(id); if (element) element.disabled = state.busy || Boolean(resuming || thermalResuming);
  }
  $('flowPressurePreconditioner').disabled = state.busy;
  $('flowSteadyAcceleration').disabled=state.busy||transient;
  $('flowLinearPolicy').disabled=state.busy;
  $('flowVelocityRelaxation').disabled=state.busy;
  $('flowPressureCorrections').disabled=state.busy;
  if ($('flowCase').value==='custom') $('flowOutletBackflow').disabled=true;
  for(const control of document.querySelectorAll('#flowBoundarySettings input, #flowBoundarySettings select, #flowBoundarySettings button'))
    control.disabled=Boolean(state.busy||resuming||thermalResuming||!state.result);
  $('flowRestartInfo').textContent = restart
    ? `可续算：t=${Number(restart.time).toPrecision(6)} s · ${restart.fileName}。启动时原生核对完整网格与状态。`
    : '每个完成的时间步都会保存；取消后可继续。';
  const dt = Number($('flowDt').value), steps = Number($('flowSteps').value);
  const start = resuming ? restart.time : 0;
  $('flowTimeHint').textContent = !transient ? '稳态流动不使用时间步长；这些设置仅在下方温度联合推进时生效。' : adaptive
    ? `从 ${start.toPrecision(5)} s 推进到目标时间；按实际 CFL 和内迭代结果缩步重试。仍是一阶时间格式，需要另做时间步细化验证。`
    : Number.isFinite(dt*steps) && dt > 0 && steps > 0
    ? `本次 ${start.toPrecision(5)} → ${(start+dt*steps).toPrecision(5)} s。一阶时间格式；时间步越小通常越准确，也更慢。`
    : '请填写正的时间步长和整数步数。';
  if (!state.busy) $('runFlow').textContent = transient ? (resuming ? '继续计算' : $('flowInitialVortex').checked?'从初始局部涡开始计算':'从静止开始计算') : '启动层流求解';
  updateThermalMode();
}
function applySharedFlowControls(request) {
  const fields = { case:'flowCase', nu:'flowNu', speed:'flowSpeed', convection:'flowConvection', outletBackflow:'flowOutletBackflow' };
  const changed = Object.entries(fields).some(([key, id]) => {
    const element=$(id); if (!element) return false;
    const value=key==='outletBackflow' ? (request[key] ?? 'reject') : request[key];
    return ['nu', 'speed'].includes(key) ? Number(element.value) !== Number(value) : element.value !== value;
  });
  if (changed) { clearFlowBinding(); clearThermalBinding(); }
  for (const [key, id] of Object.entries(fields)) {
    const element=$(id); if (element) element.value = key==='outletBackflow' ? (request[key] ?? 'reject') : request[key];
  }
  if(request.case==='custom'&&request.boundaryDefinition){
    state.flowBoundaryDefinition=structuredClone(request.boundaryDefinition);renderFlowBoundaries();
  }
}
function applyRestartControls() {
  const q = state.flowRestart;
  if ($('flowMode').value === 'steady') $('flowResume').checked = false;
  if (q && $('flowResume').checked) {
    $('thermalResume').checked = false;
    applySharedFlowControls(q);
    updateFlowScope();
  }
  updateFlowMode();
}
function renderFlowMonitor() {
  const rows = state.flowHistory || [];
  $('flowTimeline').hidden = rows.length === 0;
  const svg = $('flowMonitor'); svg.replaceChildren();
  if (!rows.length) return;
  const metric=$('flowMonitorMetric').value;
  const data=rows.filter(row => Number.isFinite(row[metric]));
  if (!data.length) return;
  const stride=Math.max(1,Math.ceil(data.length/700));
  const sampled=data.filter((_r,i)=>i%stride===0);
  if (sampled.at(-1)!==data.at(-1)) sampled.push(data.at(-1));
  const t0=data[0].time,t1=data.at(-1).time;
  let lo=Infinity,hi=-Infinity;
  for (const r of data) {lo=Math.min(lo,r[metric]);hi=Math.max(hi,r[metric]);}
  const margin=hi===lo?Math.max(1e-12,Math.abs(hi)*.02):.05*(hi-lo);lo-=margin;hi+=margin;
  const x=t=>94+592*(t-t0)/(t1-t0||1),y=v=>100-82*(v-lo)/(hi-lo);
  const add=(tag,attrs,text)=>{const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [key,value]of Object.entries(attrs))el.setAttribute(key,String(value));if(text!==undefined)el.textContent=text;svg.appendChild(el);};
  add('path',{d:'M94 12 V100 H690',fill:'none',stroke:'currentColor',opacity:.3});
  add('path',{d:sampled.map((r,i)=>`${i?'L':'M'}${x(r.time)},${y(r[metric])}`).join(' '),fill:'none',stroke:'currentColor','stroke-width':1.8});
  const last=sampled.at(-1);add('circle',{cx:x(last.time),cy:y(last[metric]),r:2.5,fill:'currentColor'});
  const digits=Math.min(12,Math.max(3,2+Math.ceil(Math.log10(Math.max(Math.abs(lo),Math.abs(hi))/(hi-lo)||1))));
  for(const [tx,ty,label]of [[90,20,hi.toPrecision(digits)],[90,101,lo.toPrecision(digits)],[135,122,`${t0.toPrecision(4)} s`],[675,122,`${t1.toPrecision(4)} s`]])
    add('text',{x:tx,y:ty,fill:'currentColor','text-anchor':'end','font-size':11},label);
  $('flowMonitorCaption').textContent=`已接受时间步 · 最新 ${last[metric].toPrecision(Math.max(5,digits))} · 曲线按需抽样显示，导出 CSV 保留本次全部步数。`;
}
function renderThermalMonitor() {
  const rows = state.thermalHistory || [];
  $('thermalTimeline').hidden = rows.length === 0;
  const svg = $('thermalMonitor'); svg.replaceChildren();
  if (!rows.length) return;
  const metric=$('thermalMonitorMetric').value;
  const data=rows.filter(row => Number.isFinite(row[metric]));
  if (!data.length) return;
  const stride=Math.max(1,Math.ceil(data.length/700));
  const sampled=data.filter((_r,i)=>i%stride===0);
  if (sampled.at(-1)!==data.at(-1)) sampled.push(data.at(-1));
  const t0=data[0].time,t1=data.at(-1).time;
  let lo=Infinity,hi=-Infinity;
  for (const r of data) {lo=Math.min(lo,r[metric]);hi=Math.max(hi,r[metric]);}
  const margin=hi===lo?Math.max(1e-12,Math.abs(hi)*.02):.05*(hi-lo);lo-=margin;hi+=margin;
  const x=t=>94+592*(t-t0)/(t1-t0||1),y=v=>100-82*(v-lo)/(hi-lo);
  const add=(tag,attrs,text)=>{const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [key,value]of Object.entries(attrs))el.setAttribute(key,String(value));if(text!==undefined)el.textContent=text;svg.appendChild(el);};
  add('path',{d:'M94 12 V100 H690',fill:'none',stroke:'currentColor',opacity:.3});
  add('path',{d:sampled.map((r,i)=>`${i?'L':'M'}${x(r.time)},${y(r[metric])}`).join(' '),fill:'none',stroke:'currentColor','stroke-width':1.8});
  const last=sampled.at(-1);add('circle',{cx:x(last.time),cy:y(last[metric]),r:2.5,fill:'currentColor'});
  const digits=Math.min(12,Math.max(3,2+Math.ceil(Math.log10(Math.max(Math.abs(lo),Math.abs(hi))/(hi-lo)||1))));
  for(const [tx,ty,label]of [[90,20,hi.toPrecision(digits)],[90,101,lo.toPrecision(digits)],[135,122,`${t0.toPrecision(4)} s`],[675,122,`${t1.toPrecision(4)} s`]])
    add('text',{x:tx,y:ty,fill:'currentColor','text-anchor':'end','font-size':11},label);
  $('thermalMonitorCaption').textContent=`已接受联合时间步 · 最新 ${last[metric].toPrecision(Math.max(5,digits))} · 曲线按需抽样显示，导出 CSV 保留本次全部步数。`;
}
const thermalPatches = ['wall', 'inlet', 'outlet', 'top', 'bottom'];
const thermalPatchId = patch => `thermal${patch[0].toUpperCase()}${patch.slice(1)}`;
function clearThermalBinding({ hidePanel = false } = {}) {
  state.thermal = null; state.thermalHistory = [];
  if (hidePanel) { state.thermalRestart = null; $('thermalResume').checked = false; $('thermalBlock').hidden = true; }
  view.setThermalFields(null);
  $('thermalResult').hidden = true; $('thermalResult').replaceChildren();
  $('thermalTimeline').hidden = true; $('thermalOption').hidden = true;
  if ($('displayMode').value === 'temperature') {
    $('displayMode').value = 'level'; view.mode = 'level'; view.draw();
  }
  if (state.mesh) renderLegend(state.mesh, state.levelBasis);
  updateThermalMode();
}
function updateThermalMode() {
  const restart = state.thermalRestart;
  if (!restart) $('thermalResume').checked = false;
  const resuming = Boolean(restart && $('thermalResume').checked);
  $('thermalResume').disabled = Boolean(state.busy || !restart);
  $('pickThermalCheckpoint').disabled=state.busy||!state.result;
  for (const input of document.querySelectorAll('#thermalBlock input[type=number], #thermalBlock select'))
    input.disabled = Boolean(state.busy || resuming);
  const vortex=$('flowMode').value!=='steady' && $('flowInitialVortex').checked && !$('flowResume').checked;
  $('runThermal').disabled = Boolean(state.busy || !state.result || $('thermalBlock').hidden || $('flowCase').value==='custom' || $('flowMode').value==='adaptive' || vortex);
  if (!state.busy) $('runThermal').textContent = resuming ? '继续温度与流动推进' : '启动温度与流动推进';
  $('thermalRestartInfo').textContent = restart
    ? `联合续算状态：t=${Number(restart.time).toPrecision(6)} s。物性、源项与边界锁定；可调整时间步和迭代控制。`
    : '每个流动与温度均收敛的时间步保存联合状态；取消后可续算。';
  const start = resuming ? Number(restart.time) : 0;
  const duration = Number($('flowDt').value) * Number($('flowSteps').value);
  $('thermalTimeHint').textContent = vortex ? '温度联合推进尚不支持初始局部涡；请先关闭该初始条件。' : $('flowMode').value==='adaptive'
    ? '温度联合推进目前使用固定步长；请在上方切换为固定步长模式。'
    : Number.isFinite(duration) && duration > 0
    ? `温度始终非定常：本次 ${start.toPrecision(5)} → ${(start + duration).toPrecision(5)} s。温度积分需乘 ρcp 才是单位深度热量。`
    : '请在上方填写时间步长与本次步数。';
}
function applyThermalRestartControls() {
  const request = state.thermalRestart?.request;
  if (request && $('thermalResume').checked) {
    $('flowResume').checked = false;
    applySharedFlowControls(request);
    for (const [field, id] of Object.entries({ diffusivity:'thermalDiffusivity', initial:'thermalInitial', source:'thermalSource', scalarConvection:'thermalConvection' }))
      $(id).value = request[field];
    for (const patch of thermalPatches) {
      const boundary = request.boundaries[patch], id = thermalPatchId(patch);
      $(`${id}Kind`).value = boundary.kind; $(`${id}Value`).value = boundary.value; $(`${id}Inflow`).value = boundary.inflowValue;
    }
    updateFlowScope();
  }
  updateFlowMode();
}
function thermalRequest() {
  const boundaries = {};
  for (const patch of thermalPatches) {
    const id = thermalPatchId(patch);
    boundaries[patch] = { kind:$(`${id}Kind`).value, value:Number($(`${id}Value`).value), inflowValue:Number($(`${id}Inflow`).value) };
  }
  return { case:$('flowCase').value, nu:Number($('flowNu').value), speed:Number($('flowSpeed').value),
    convection:$('flowConvection').value, pressurePreconditioner:$('flowPressurePreconditioner').value,
    outletBackflow:$('flowOutletBackflow').value,
    maxIterations:Number($('flowMaxIterations').value), tolerance:Number($('flowTolerance').value), dt:Number($('flowDt').value), steps:Number($('flowSteps').value),
    resume:$('thermalResume').checked, diffusivity:Number($('thermalDiffusivity').value), initial:Number($('thermalInitial').value),
    source:Number($('thermalSource').value), scalarConvection:$('thermalConvection').value, boundaries };
}
function validThermalInputs() {
  if($('flowLinearPolicy').value!=='strict' || Number($('flowVelocityRelaxation').value)!==.6 || Number($('flowPressureCorrections').value)!==4) {
    status('温度联算的迭代设置需要调整','请恢复固定线性精度、0.6速度松弛和4次压力校正。');return false;
  }
  if($('flowPressurePreconditioner').value==='cholesky') {
    status('温度联算的压力求解设置需要调整','系统稀疏 Cholesky 暂仅用于独立层流；请选择 IC0 或多重网格。');return false;
  }
  for (const input of document.querySelectorAll('#flowBlock input[type=number], #thermalBlock input[type=number]')) {
    if (!input.disabled && (!input.value.trim() || !input.checkValidity())) {
      const details = input.closest('details'); if (details) details.open = true;
      input.reportValidity(); input.focus(); status('温度工况参数需要调整', '请填写有效数值，并检查范围。'); return false;
    }
  }
  return true;
}
function boundedThermalHistory(rows) {
  const accepted = rows.filter(row => row.accepted !== false && Number.isFinite(row.time));
  const stride = Math.max(1, Math.ceil(accepted.length / 1400));
  return accepted.filter((_row, i) => i % stride === 0 || i === accepted.length - 1);
}
function bindThermal(payload) {
  view.setThermalFields(payload.fields.cells);
  state.thermal = payload; state.thermalHistory = boundedThermalHistory(payload.history || []);
  $('thermalOption').hidden = false; $('displayMode').value = 'temperature'; view.mode = 'temperature';
  $('canvasWrap').classList.remove('light'); view.draw(); renderLegend(state.mesh, state.levelBasis);
  const summary = payload.summary, last = state.thermalHistory.at(-1);
  const container = $('thermalResult'); container.replaceChildren();
  const heading = document.createElement('div'); heading.className = 'flow-state';
  heading.textContent = `温度与流动已接受 t=${Number(summary.acceptedTime ?? summary.time).toPrecision(6)} s · 本次 ${summary.steps} 步 · 单向恒物性输运`;
  container.appendChild(heading);
  for (const [label, value] of [['时间步（s）',summary.dt], ['最低温度（K）',summary.minValue], ['最高温度（K）',summary.maxValue],
    ['热扩散率（m²/s）',summary.diffusivity], ['温度积分（K·m²）',last?.heatContent],
    ['全域平衡（K·m²/s）',last?.globalBalance], ['温度残差',last?.scalarResidual], ['最后一步最大 CFL',last?.maxCourant]]) {
    const item = document.createElement('div'), caption = document.createElement('span'), number = document.createElement('b');
    caption.textContent = label; number.textContent = Number.isFinite(value) ? value.toExponential(4) : '未记录';
    item.append(caption, number); container.appendChild(item);
  }
  container.hidden = false; renderThermalMonitor();
}
async function refreshThermalState(restoreResult = false) {
  const saved = await window.cartmesh.thermalState();
  state.thermalRestart = saved.restart;
  if (restoreResult && saved.thermal) bindThermal(saved.thermal);
  updateThermalMode(); return saved;
}
async function runThermal() {
  if (state.busy || !state.result || !state.mesh || !validThermalInputs()) return;
  const request = thermalRequest();
  clearThermalBinding(); setBusy(true); $('runThermal').textContent = '正在推进…';
  status('温度与流动推进中', '只有流动和温度均收敛才接受时间步；可取消并从联合状态继续。');
  try {
    const payload = await window.cartmesh.runThermal(request);
    bindThermal(payload); await refreshThermalState();
    $('thermalResume').checked = Boolean(state.thermalRestart);
    status('温度时间推进完成', `已接受 t=${Number(payload.summary.acceptedTime ?? payload.summary.time).toPrecision(6)} s；温度色图使用最终网格，可导出真实场和联合续算状态。`);
  } catch (error) {
    const message = error.message.replace(/^Error invoking remote method '[^']+': Error: /, '');
    const saved = await refreshThermalState(true).catch(() => null);
    if (saved?.restart) $('thermalResume').checked = true;
    status(/取消/.test(message) ? '温度推进已取消' : '温度推进未完成', message.split('\n')[0] + (saved?.thermal ? ' 当前显示上次完整温度结果。' : '') + (saved?.restart ? ` 可从 t=${saved.restart.time} s 联合续算。` : ''));
    log(message);
  } finally { setBusy(false); applyThermalRestartControls(); }
}

function bindFlow(payload) {
  state.flow=payload; state.flowHistory=payload.history || [];
  view.setFlowFields(payload.fields.cells);
  $('flowSpeedOption').hidden=false; $('flowPressureOption').hidden=false;
  $('displayMode').value='speed';view.mode='speed';view.draw();
  renderLegend(state.mesh,state.levelBasis);renderFlowResult(payload.summary);renderFlowMonitor();
}
async function refreshFlowState(restoreResult=false) {
  const saved=await window.cartmesh.flowState();
  state.flowRestart=saved.restart;
  if (restoreResult && saved.flow) bindFlow(saved.flow);
  updateFlowMode();
  return saved;
}
function flowRequest() {
  const transient=$('flowMode').value!=='steady';
  const request={ case:$('flowCase').value,nu:Number($('flowNu').value),speed:Number($('flowSpeed').value),
    maxIterations:Number($('flowMaxIterations').value),convection:$('flowConvection').value,
    tolerance:Number($('flowTolerance').value),
    pressurePreconditioner:$('flowPressurePreconditioner').value,
    steadyAcceleration:transient?'none':$('flowSteadyAcceleration').value,
    linearPolicy:$('flowLinearPolicy').value,velocityRelaxation:Number($('flowVelocityRelaxation').value),pressureCorrectionPasses:Number($('flowPressureCorrections').value),
    outletBackflow:$('flowOutletBackflow').value,
    mode:$('flowMode').value,dt:Number($('flowDt').value),steps:Number($('flowSteps').value),
    endTime:Number($('flowEndTime').value),minDt:Number($('flowMinDt').value),maxCourant:Number($('flowMaxCourant').value),
    maxRetries:Number($('flowMaxRetries').value),maxSteps:Number($('flowMaxSteps').value),resume:transient && $('flowResume').checked };
  if(request.case==='custom')request.boundaryDefinition=state.flowBoundaryDefinition;
  if(transient && !request.resume && $('flowInitialVortex').checked)
    request.initialVortex={centre:[Number($('flowVortexX').value),Number($('flowVortexY').value)],
      radius:Number($('flowVortexRadius').value),peakSpeed:Number($('flowVortexSpeed').value)};
  return request;
}
function applyFlowCase(request) {
  $('flowResume').checked=false;
  $('thermalResume').checked=false;
  clearFlowBinding(); clearThermalBinding();
  state.flowBoundaryDefinition = request.boundaryDefinition ? structuredClone(request.boundaryDefinition) : null;
  applySharedFlowControls(request);
  const fields = {mode:'flowMode',maxIterations:'flowMaxIterations',tolerance:'flowTolerance',pressurePreconditioner:'flowPressurePreconditioner',steadyAcceleration:'flowSteadyAcceleration',linearPolicy:'flowLinearPolicy',velocityRelaxation:'flowVelocityRelaxation',pressureCorrectionPasses:'flowPressureCorrections',
    dt:'flowDt',steps:'flowSteps',endTime:'flowEndTime',minDt:'flowMinDt',maxCourant:'flowMaxCourant',
    maxRetries:'flowMaxRetries',maxSteps:'flowMaxSteps'};
  for (const [key,id] of Object.entries(fields)) if (request[key] !== undefined) $(id).value=request[key];
  $('flowInitialVortex').checked=Boolean(request.initialVortex);
  if (request.initialVortex) {
    const q=request.initialVortex;
    $('flowVortexX').value=q.centre[0];$('flowVortexY').value=q.centre[1];
    $('flowVortexRadius').value=q.radius;$('flowVortexSpeed').value=q.peakSpeed;
  }
  view.setBoundaryHighlight([]);renderFlowBoundaries();updateFlowScope();updateFlowMode();
}
async function saveFlowCase() {
  if (state.busy || !state.result || !validFlowInputs() || $('saveFlowCase').disabled) return null;
  const request=flowRequest();
  setBusy(true);
  try {
    const saved=await window.cartmesh.saveFlowCase(request);
    if(saved){$('flowCaseFileInfo').textContent=`已保存 ${saved.fileName}：从零起算的流动设置，绑定当前最终网格。`;status('流动工况已保存','包括当前边界、物性、求解和初始设置。');}
    return saved;
  } catch(error){status('工况保存失败',error.message);log(error.message);return null;}
  finally{setBusy(false);}
}
async function loadFlowCase() {
  if(state.busy || !state.result)return null;
  setBusy(true);
  try {
    const saved=await window.cartmesh.loadFlowCase();
    if(saved){applyFlowCase(saved.document.request);$('flowCaseFileInfo').textContent=`已读取 ${saved.fileName}，绑定当前最终网格。后续编辑以当前表单为准。`;
      status('流动工况已读取','设置已恢复，从零开始计算。结果导出仍保存上次已完成的计算；温度设置需单独核对。');}
    return saved;
  }catch(error){status('工况读取失败',error.message);log(error.message);return null;}
  finally{setBusy(false);}
}
const EULER_CONTROLS={thermalConductivity:'eulerConductivity',wallValue:'eulerWallValue',density:'eulerDensity',u:'eulerU',v:'eulerV',pressure:'eulerPressure',gamma:'eulerGamma',gasConstant:'eulerGasConstant',split:'eulerSplit',endTime:'eulerEndTime',maximumStep:'eulerMaximumStep',minimumStep:'eulerMinimumStep',cfl:'eulerCfl',maximumSteps:'eulerMaximumSteps',maximumSeconds:'eulerMaximumSeconds'};
const EULER_PHYSICAL=['thermalConductivity','wallValue','density','u','v','pressure','gamma','gasConstant','split'];
const EULER_OPTIONS=['eulerDensityOption','eulerPressureOption','eulerTemperatureOption','eulerMachOption','eulerSpeedOption'];
function eulerRequest() {
  const request={case:$('eulerCase').value,resume:$('eulerResume').checked,fluxScheme:$('eulerFluxScheme').value,order:Number($('eulerOrder').value),wallThermal:$('eulerWallThermal').value};
  for(const [key,id] of Object.entries(EULER_CONTROLS))request[key]=Number($(id).value);
  if(request.case!=='sod')request.split=.5;
  if(request.wallThermal==='insulated')request.wallValue=0;
  return request;
}
function updateEulerControls() {
  const restart=state.eulerRestart;if(!restart)$('eulerResume').checked=false;
  const resuming=Boolean(restart&&$('eulerResume').checked),sod=$('eulerCase').value==='sod';
  $('eulerResume').disabled=Boolean(state.busy||!restart);$('eulerCase').disabled=Boolean(state.busy||resuming);
  for(const [key,id] of Object.entries(EULER_CONTROLS))$(id).disabled=Boolean(state.busy||(resuming&&EULER_PHYSICAL.includes(key))||(sod&&['u','v'].includes(key)));
  if(sod){$('eulerU').value=0;$('eulerV').value=0;}
  $('eulerSplitField').hidden=!sod;
  const conduction=Number($('eulerConductivity').value)>0,hasWall=$('eulerCase').value!=='uniform';
  if(!conduction||!hasWall)$('eulerWallThermal').value='insulated';
  const thermalKind=$('eulerWallThermal').value;
  $('eulerWallThermal').disabled=Boolean(state.busy||resuming||!conduction||!hasWall);
  $('eulerWallValue').disabled=Boolean(state.busy||resuming||thermalKind==='insulated');
  $('eulerWallValueField').hidden=thermalKind==='insulated';
  $('eulerWallValueLabel').textContent=thermalKind==='temperature'?'壁温（K）':'向外热流（W/m²；负值加热）';
  if(thermalKind==='temperature')$('eulerWallValue').min='0.000000000001';else $('eulerWallValue').removeAttribute('min');
  $('eulerFluxScheme').disabled=Boolean(state.busy);$('eulerOrder').disabled=Boolean(state.busy);
  $('pickEulerCheckpoint').disabled=Boolean(state.busy||!state.result);
  $('runEuler').disabled=Boolean(state.busy||!state.result||$('eulerBlock').hidden);
  if(!state.busy)$('runEuler').textContent=resuming?'继续可压计算':'启动可压计算';
  $('eulerScope').textContent=sod?'仅限轴对齐矩形。左侧使用下方密度和压力，右侧分别为其0.125倍、0.1倍；初始速度为零，上下自由滑移、左右透射。':$('eulerCase').value==='external'?'下方为初始场与远场状态。实际物面自由滑移，外域施加特征远场；当前不会形成黏性边界层。':$('eulerCase').value==='sealed'?'下方为初始场；所有物面与外边界均为固定自由滑移壁。封闭系统只通过设定的壁面导热交换能量。':'无物面的均匀初始场，所有外边界施加相同特征远场。';
  $('eulerRestartInfo').textContent=restart?`可续算到 t=${restart.time.toPrecision(6)} s 的已接受状态（${restart.steps}步）。物理参数锁定，目标时间需更晚；取消续算则从初始场重新计算。`:'每25步保存一次，正常结束或取消时再次保存。重开App后，先生成同一网格，再选择结果目录的 desktop-state.json。';
}
function applyEulerRestart() {
  if($('eulerResume').checked&&state.eulerRestart) {
    const request=state.eulerRestart.request;$('eulerCase').value=request.case;
    for(const key of EULER_PHYSICAL)$(EULER_CONTROLS[key]).value=request[key]??0;
    $('eulerWallThermal').value=request.wallThermal??'insulated';
    $('eulerFluxScheme').value=request.fluxScheme??'rusanov';$('eulerOrder').value=String(request.order??1);
  }
  updateEulerControls();
}
function clearEulerBinding(hidePanel=false) {
  state.euler=null;state.eulerHistory=[];view.setEulerFields(null);
  $('eulerResult').hidden=true;$('eulerResult').replaceChildren();$('eulerTimeline').hidden=true;
  for(const id of EULER_OPTIONS)$(id).hidden=true;
  if($('displayMode').value.startsWith('euler-')){$('displayMode').value='level';view.mode='level';view.draw();}
  if(hidePanel){state.eulerRestart=null;$('eulerResume').checked=false;$('eulerBlock').hidden=true;}
  if(state.mesh)renderLegend(state.mesh,state.levelBasis);
}
function renderEulerMonitor() {
  const rows=state.eulerHistory||[];$('eulerTimeline').hidden=!rows.length;const svg=$('eulerMonitor');svg.replaceChildren();if(!rows.length)return;
  const metric=$('eulerMonitorMetric').value,values=rows.map(r=>r[metric]);let min=Infinity,max=-Infinity;
  for(const value of values){min=Math.min(min,value);max=Math.max(max,value);}
  const first=rows[0].time,last=rows.at(-1).time,span=max-min||Math.max(Math.abs(max)*.02,1e-12);
  const sampled=rows.filter((_r,i)=>i%Math.max(1,Math.ceil(rows.length/1000))===0||i===rows.length-1);
  const line=document.createElementNS('http://www.w3.org/2000/svg','polyline');
  line.setAttribute('points',sampled.map(r=>`${55+620*(r.time-first)/(last-first||1)},${108-85*(r[metric]-min)/span}`).join(' '));
  line.setAttribute('fill','none');line.setAttribute('stroke','#61b7dd');line.setAttribute('stroke-width','2');svg.appendChild(line);
  for(const [x,y,value] of [[6,20,max.toPrecision(4)],[6,110,min.toPrecision(4)],[55,128,first.toPrecision(4)+' s'],[600,128,last.toPrecision(4)+' s']]) {
    const text=document.createElementNS('http://www.w3.org/2000/svg','text');text.setAttribute('x',x);text.setAttribute('y',y);text.setAttribute('fill','currentColor');text.setAttribute('font-size','11');text.textContent=value;svg.appendChild(text);
  }
  $('eulerMonitorCaption').textContent=`已接受状态 · 最新 ${values.at(-1).toPrecision(6)} · 运行时每25步更新，完成后读回全部历史；导出 CSV 保留每一步。`;
}
function bindEuler(payload) {
  state.euler=payload;state.eulerHistory=payload.history.map(row=>({...row}));view.setEulerFields(payload.fields.cells);
  for(const id of EULER_OPTIONS)$(id).hidden=false;
  $('displayMode').value='euler-rho';view.mode='euler-rho';view.draw();renderLegend(state.mesh,state.levelBasis);
  const container=$('eulerResult');container.replaceChildren();container.hidden=false;
  const title=document.createElement('div');title.className='flow-state';title.textContent=`可压 Euler · 已到达 t=${payload.summary.time.toPrecision(6)} s · 本次 ${payload.summary.acceptedSteps} 步 · ${payload.request.fluxScheme==='hllc'?'HLLC/HLLE':'Rusanov'} / ${payload.request.order===2?'二阶':'一阶'}无黏模型`;container.appendChild(title);
  const last=payload.history.at(-1),lines=[`最小密度 ${last.minimumDensity.toPrecision(6)} kg/m³；最小绝对压力 ${last.minimumPressure.toPrecision(6)} Pa`,
    `最后步声学 CFL ${payload.audit.acousticCourant.toPrecision(4)}；逐格守恒相对误差 ${payload.audit.maximumCellBalanceRelative.toExponential(2)}`,
    `导热系数 ${payload.request.thermalConductivity??0} W/m/K；最后步组合 CFL ${(payload.audit.combinedCourant??payload.audit.acousticCourant).toPrecision(4)}；向外净热流 ${(payload.audit.boundaryHeat??0).toPrecision(6)} W/m`,
    `本次 HLLC 通量回退 ${payload.summary.hllcFallbackEvaluations??0} 次；重构退阶 ${payload.summary.reconstructionFallbackCells??0} 个单元阶段`,
    '到达目标时间不代表稳态或任意工况精度合格；导热通过总能量耦合；没有黏性应力或湍流。'];
  for(const line of lines){const p=document.createElement('p');p.className='note';p.textContent=line;container.appendChild(p);}renderEulerMonitor();
}
async function refreshEulerState(show=false) {
  const saved=await window.cartmesh.eulerState();state.eulerRestart=saved.restart||null;
  if(show&&saved.euler)bindEuler(saved.euler);return saved;
}
async function runEuler() {
  if(state.busy||!state.result||!state.mesh)return;
  for(const input of document.querySelectorAll('#eulerBlock input[type=number]'))if(!input.disabled&&(!input.value.trim()||!input.checkValidity())){input.reportValidity();status('可压参数需要调整','检查气体、时间和计算预算。');return;}
  const request=eulerRequest();clearEulerBinding();setBusy(true);$('runEuler').textContent='可压推进中…';status('可压 Euler 推进中','显式声学时间步；取消后保留已接受状态。');
  try {
    const payload=await window.cartmesh.runEuler(request);bindEuler(payload);await refreshEulerState();$('eulerResume').checked=true;
    status('可压时间推进完成',`t=${payload.summary.time.toPrecision(6)} s；可查看密度、绝对压力、温度、Mach数和速度，并导出实际场。`);
  }catch(error) {
    const message=error.message.replace(/^Error invoking remote method '[^']+': Error: /,'');
    const saved=await refreshEulerState(true).catch(()=>null);$('eulerResume').checked=Boolean(saved?.restart);
    status(/取消|cancelled/.test(message)?'可压计算已取消':'可压计算未完成',message.split('\n')[0]+(saved?.euler?' 显示上次完整结果。':''));log(message);
  }finally{setBusy(false);applyEulerRestart();}
}
$('runEuler').addEventListener('click',runEuler);
$('eulerResume').addEventListener('change',applyEulerRestart);
$('eulerCase').addEventListener('change',()=>{if(state.euler)clearEulerBinding();updateEulerControls();});
$('eulerMonitorMetric').addEventListener('change',renderEulerMonitor);
for(const id of ['eulerFluxScheme','eulerOrder','eulerWallThermal'])$(id).addEventListener('change',()=>{if(state.euler)clearEulerBinding();updateEulerControls();});
for(const input of document.querySelectorAll('#eulerBlock input[type=number]'))input.addEventListener('change',()=>{if(state.euler)clearEulerBinding();updateEulerControls();});
$('pickEulerCheckpoint').addEventListener('click',async()=>{
  if(state.busy||!state.result)return;setBusy(true);
  try{const selected=await window.cartmesh.pickEulerCheckpoint();if(selected){state.eulerRestart=selected;$('eulerResume').checked=true;applyEulerRestart();status('可压续算清单已载入','原生计算开始前还会核对完整网格与物理边界。');}}
  catch(error){status('可压续算载入失败',error.message);log(error.message);}finally{setBusy(false);}
});
window.cartmesh.onEulerProgress(progress=>{
  if(!state.busy||state.euler)return;
  state.eulerHistory.push(progress);if(state.eulerHistory.length>1400)state.eulerHistory=state.eulerHistory.filter((_r,i)=>i%2===0||i===state.eulerHistory.length-1);
  renderEulerMonitor();status('可压 Euler 推进中',`已接受 ${progress.step} 步 · t=${progress.time.toPrecision(6)} s · 声学 CFL ${progress.acousticCourant.toPrecision(4)}`);
});

async function runFlow() {
  if (state.busy || !state.result || !state.mesh || !validFlowInputs()) return;
  const request=flowRequest(),transient=request.mode!=='steady';
  clearFlowBinding();setBusy(true);$('runFlow').textContent='正在求解…';
  status('层流求解中',transient?'按物理时间推进；取消后可从最后接受的时间步继续。':'SIMPLE 速度—压力耦合；可随时取消。');
  try {
    const payload=await window.cartmesh.runFlow(request);
    bindFlow(payload);
    if (transient) status('本次时间推进完成',`已接受到 t=${payload.summary.acceptedTime.toPrecision(6)} s；最大 CFL ${payload.summary.maxCourant.toPrecision(4)}。这不等于达到稳态或已验证物理精度。`);
    else status(payload.summary.converged?'层流求解已收敛':'到达迭代上限，未收敛',`${payload.summary.iterations} 次迭代；可切换速度或压力色图并导出。`);
    await refreshFlowState();
    if (transient) $('flowResume').checked=Boolean(state.flowRestart);
  } catch (error) {
    const message=error.message.replace(/^Error invoking remote method '[^']+': Error: /,'');
    const saved=await refreshFlowState(true).catch(()=>null);
    if (transient && saved?.restart) $('flowResume').checked=true;
    status(/取消/.test(message)?'层流求解已取消':'本次计算未完成',message.split('\n')[0]+(saved?.flow?' 当前显示上次完整结果。':''));
    log(message);
  } finally {setBusy(false);applyRestartControls();}
}

async function importGeometryFile(picked) {
  if (!picked || state.busy) return;
  if (/\.(png|jpe?g)$/i.test(picked)) {
    try {
      const imported = await window.CartMeshRasterImport.open(picked);
      if (!imported) return;
      $('sample').value = '';
      $('sourceUnits').value = 'm';
      $('referenceMode').value = 'bbox';
      $('referenceLengthField').hidden = true;
      await chooseGeometry(imported.geometryPath, imported.label, null);
    } catch (error) { status('图片导入失败', error.message); log(error.message); }
  } else {
    $('sample').value = '';
    await chooseGeometry(picked, picked.split(/[\\/]/).pop(), null);
  }
}
$('pickGeometry').addEventListener('click', async () => {
  try { await importGeometryFile(await window.cartmesh.pickGeometry()); }
  catch (error) { status('导入失败', error.message); log(error.message); }
});

$('sample').addEventListener('change', async event => {
  const sample = state.catalog.samples.find(item => item.id === event.target.value);
  if (!sample) return;
  await chooseGeometry(sample.path, sample.label, sample);
});

for (const id of ['farFieldSpans', 'wallCellsPerSpan', 'cellsPerLevel', 'farLevel']) {
  $(id).addEventListener('input', updateBudget);
}
for (const id of ['wallRelativeSize','backgroundRelativeSize','relativePadding','relativeBandCells','referenceLength'])
  $(id).addEventListener('input', () => { updateBudget(); syncRegions(); });
$('relativeAllowUnsafe').addEventListener('change', () => { updateBudget(); schedulePlan(); });
$('referenceMode').addEventListener('change', () => {
  $('referenceLengthField').hidden = $('referenceMode').value !== 'explicit';
  updateBudget();
  schedulePlan();
  syncRegions();
});
$('allowUnsafe').addEventListener('change', updateBudget);
$('useCurvature').addEventListener('change', event => {
  $('curvatureCellsPerRadius').disabled = !event.target.checked;
});
$('useGap').addEventListener('change', event => {
  $('gapCells').disabled = !event.target.checked;
});
$('useWake').addEventListener('change', event => {
  $('wakeFields').hidden = !event.target.checked;
});
$('probe').addEventListener('click', probeSizing);
$('probeRelative').addEventListener('click', probeSizing);
$('generate').addEventListener('click', generate);
$('runFlow').addEventListener('click', runFlow);
$('saveFlowCase').addEventListener('click', saveFlowCase);
$('loadFlowCase').addEventListener('click', loadFlowCase);
$('runThermal').addEventListener('click', runThermal);
$('pickThermalCheckpoint').addEventListener('click',async()=>{
  if(state.busy||!state.result)return;
  setBusy(true);
  try{const picked=await window.cartmesh.pickThermalCheckpoint();if(picked){state.thermalRestart=picked;$('thermalResume').checked=true;applyThermalRestartControls();status('联合状态已选择','开始计算时原生求解器还会核对完整网格与物理配置。');}}
  catch(error){status('联合状态载入失败',error.message);log(error.message);}
  finally{setBusy(false);}
});
$('thermalResume').addEventListener('change', applyThermalRestartControls);
$('thermalMonitorMetric').addEventListener('change', renderThermalMonitor);
for (const control of document.querySelectorAll('#thermalBlock input[type=number], #thermalBlock select'))
  control.addEventListener('change', () => { if (state.thermal) clearThermalBinding(); });
$('flowMode').addEventListener('change',applyRestartControls);
$('flowInitialVortex').addEventListener('change',updateFlowMode);
$('flowResume').addEventListener('change',applyRestartControls);
for(const id of ['flowDt','flowSteps','flowEndTime','flowMaxCourant','flowMinDt','flowMaxRetries','flowMaxSteps']) $(id).addEventListener('input',updateFlowMode);
$('flowMonitorMetric').addEventListener('change',renderFlowMonitor);
$('pickFlowCheckpoint').addEventListener('click',async()=>{
  if(state.busy)return;
  setBusy(true);
  try {const picked=await window.cartmesh.pickFlowCheckpoint();if(picked){state.flowRestart=picked;$('flowResume').checked=true;applyRestartControls();}}
  catch(error){status('无法读取重启状态',error.message);}
  finally {setBusy(false);applyRestartControls();}
});
$('flowCase').addEventListener('change', () => { clearFlowBinding(); clearThermalBinding(); updateFlowScope(); });
$('generateFlowBoundaries').addEventListener('click',()=>prepareFlowBoundaries($('flowBoundaryPreset').value));
$('importFlowBoundaries').addEventListener('click',()=>prepareFlowBoundaries('import'));
for (const id of ['flowNu', 'flowSpeed', 'flowMaxIterations', 'flowTolerance']) {
  $(id).addEventListener('input', () => { if (state.flow) clearFlowBinding(); if (state.thermal) clearThermalBinding(); });
}
$('flowConvection').addEventListener('change', () => { if (state.flow) clearFlowBinding(); if (state.thermal) clearThermalBinding(); });
$('flowPressurePreconditioner').addEventListener('change', () => { if (state.flow) clearFlowBinding(); if (state.thermal) clearThermalBinding(); });
$('flowSteadyAcceleration').addEventListener('change', () => { if (state.flow) clearFlowBinding(); });
$('flowLinearPolicy').addEventListener('change', () => { if (state.flow) clearFlowBinding(); });
$('flowVelocityRelaxation').addEventListener('input', () => { if (state.flow) clearFlowBinding(); });
$('flowPressureCorrections').addEventListener('input', () => { if (state.flow) clearFlowBinding(); });
$('flowOutletBackflow').addEventListener('change', () => { if (state.flow) clearFlowBinding(); if (state.thermal) clearThermalBinding(); });
$('addRegion').addEventListener('click', addRegion);

$('displayMode').addEventListener('change', event => {
  view.mode = event.target.value;
  $('canvasWrap').classList.toggle('light', view.mode === 'light');
  view.draw();
  if (state.mesh) renderLegend(state.mesh, state.levelBasis);
});
$('toggleRegions').addEventListener('click', () => {
  view.showRegions = !view.showRegions;
  $('toggleRegions').classList.toggle('active', view.showRegions);
  view.draw();
});

$('fitOverview').addEventListener('click', fitOverview);
$('controlMode').addEventListener('change', updateControlMode);
$('fluidRegion').addEventListener('change', () => {
  clearResult();
  selectMethod(state.method);
  if (state.geometryPath) drawGeometryOutline();
});

$('fitDomain').addEventListener('click', () => {
  if (state.mesh) view.fitTo(state.mesh.bounds);
});
$('fitWall').addEventListener('click', () => {
  if (state.wallBounds) view.fitTo(state.wallBounds, 0.18);
});
$('toggleGrid').addEventListener('click', () => {
  view.showGrid = !view.showGrid;
  $('toggleGrid').classList.toggle('active', view.showGrid);
  view.draw();
});
$('cancel').addEventListener('click', () => window.cartmesh.cancel());
$('exportResult').addEventListener('click', async () => {
  if (state.busy) return;
  setBusy(true);
  try {
    const destination = await window.cartmesh.exportResult();
    if (destination) status('结果包已保存', destination);
  } catch (error) { status('保存失败', error.message); }
  finally { setBusy(false); }
});

let progressTimer;
window.cartmesh.onProgress(progress => {
  clearInterval(progressTimer);
  const started = Date.now();
  const update = () => {
    const elapsed = Math.floor((Date.now() - started) / 1000);
    const remaining = Math.ceil(progress.estimatedSeconds - elapsed);
    status('生成中', `第 ${progress.attempt}/${progress.maximum} 组 · 本组已用 ${elapsed} 秒 · ` +
      (remaining > 0 ? `预计还需约 ${remaining} 秒` : '已超出粗估，仍在计算，可取消') +
      `（${progress.estimateSource}）`);
  };
  update();
  progressTimer = setInterval(update, 1000);
});
window.cartmesh.onFlowProgress(progress => {
  if (progress.type==='flow-time-retry') {
    const reason=progress.reason==='courant'?'CFL 超限':'内迭代未收敛';
    const message=`${reason}；已接受到 t=${progress.acceptedTime.toPrecision(6)} s。步长 ${progress.dt.toExponential(3)} → ${progress.nextDt.toExponential(3)} s 后重试。`;
    status('当前试算未接受',message);log(message);return;
  }
  if (progress.type==='flow-time-step') {
    state.flowHistory.push(progress);
    if (state.flowHistory.length>1400) state.flowHistory=state.flowHistory.filter((_r,i)=>i%2===0 || i===state.flowHistory.length-1);
    renderFlowMonitor();
    status('非定常计算中',`已接受第 ${progress.step} 步 · t=${progress.time.toPrecision(6)} s · 最大 CFL ${progress.maxCourant.toPrecision(4)}`);
    return;
  }
  status(progress.time===undefined?'层流求解中':`候选时间步 ${progress.timeStep} · t=${progress.time.toPrecision(6)} s`,
    `内迭代 ${fmt(progress.iteration)} · 连续性 ${Number(progress.continuity).toExponential(2)} · 动量残差 ${Number(progress.momentumResidual).toExponential(2)}`);
});
window.cartmesh.onThermalProgress(progress => {
  if (progress.type !== 'thermal-time-step' || progress.accepted === false) return;
  state.thermalHistory.push(progress);
  if (state.thermalHistory.length > 1400) state.thermalHistory = boundedThermalHistory(state.thermalHistory);
  renderThermalMonitor();
  status('温度与流动推进中', `已接受第 ${progress.step} 步 · t=${Number(progress.time).toPrecision(6)} s · 全域平衡 ${Number(progress.globalBalance).toExponential(3)} K·m²/s`);
});
window.cartmesh.onRunLine(log);
window.addEventListener('resize', () => view.requestDraw());

(async () => {
  state.catalog = await window.cartmesh.catalog();
  $('flowPressurePreconditioner').querySelector('option[value="cholesky"]').disabled=window.cartmesh.platform!=='darwin';
  renderMethods();
  renderPresets();
  renderSamples();
  selectMethod('cutcell');
  renderRegions();
  // Smoke tests drive these same handlers; an optional output override retains fixtures.
  window.__smoke = { state, selectMethod, chooseGeometry, generate, runFlow, runThermal, runEuler, eulerRequest, flowRequest, saveFlowCase, loadFlowCase, prepareFlowBoundaries, setOutput, addRegion, renderRegions, view, loadVerifiedPreset, setSidebarCollapsed, returnToStart, importGeometryFile };
})();

function setOutput(directory) {
  state.outputDirectory = directory;
  updateReady();
}

// ---------------------------------------------------------------------------
// Hand-placed refinement regions.
//
// The size field guarantees the refinement that is not negotiable — the band along
// the wall.  Everything else (a wake, a downstream patch, a region you just want to
// look at) is a judgement call, so it is edited by hand here and drawn on the canvas
// in the same body-span frame the numbers are written in.

const REGION_FIELDS = [
  ['xmin', '起点 x'], ['xmax', '终点 x'], ['ymin', '下边界 y'], ['ymax', '上边界 y']
];

function syncRegions() {
  view.regions = state.regions;
  view.frame = state.frame && { ...state.frame, bodySpan: $('referenceMode').value === 'explicit'
    ? Number($('referenceLength').value) : state.frame.bodySpan };
  view.draw();
}

function addRegion() {
  // A downstream patch is the common case, so the default is one: from just behind the
  // body out to six body lengths, four levels coarser than the wall.
  state.regions.push({ xmin: 0.6, xmax: 6, ymin: -0.8, ymax: 0.8, levelsBelowWall: 4, active: true });
  renderRegions();
}

function renderRegions() {
  const list = $('regionList');
  list.replaceChildren();
  state.regions.forEach((box, index) => {
    const card = document.createElement('div');
    card.className = `region${box.active ? ' active' : ''}`;

    const top = document.createElement('div');
    top.className = 'regionTop';
    const name = document.createElement('b');
    name.textContent = `R${index + 1}`;
    const depth = document.createElement('span');
    depth.className = 'depth';
    const depthInput = document.createElement('input');
    depthInput.type = 'number';
    depthInput.min = '0';
    depthInput.max = '20';
    depthInput.step = '1';
    depthInput.value = String(box.levelsBelowWall);
    depthInput.addEventListener('input', () => {
      box.levelsBelowWall = Number(depthInput.value);
      syncRegions();
    });
    depth.append(document.createTextNode('低于壁面'), depthInput, document.createTextNode('级'));
    const grow = document.createElement('span');
    grow.className = 'grow';
    const drop = document.createElement('button');
    drop.className = 'drop';
    drop.textContent = '×';
    drop.title = '删除这个加密区';
    drop.addEventListener('click', () => {
      state.regions.splice(index, 1);
      renderRegions();
    });
    top.append(name, grow, depth, drop);

    const grid = document.createElement('div');
    grid.className = 'grid4';
    for (const [key, label] of REGION_FIELDS) {
      const field = document.createElement('label');
      const caption = document.createElement('span');
      caption.textContent = label;
      const input = document.createElement('input');
      input.type = 'number';
      input.step = 'any';
      input.value = String(box[key]);
      input.addEventListener('input', () => {
        box[key] = Number(input.value);
        syncRegions();
      });
      input.addEventListener('focus', () => {
        state.regions.forEach(other => { other.active = false; });
        box.active = true;
        [...list.children].forEach((item, i) => item.classList.toggle('active', i === index));
        syncRegions();
      });
      field.append(caption, input);
      grid.appendChild(field);
    }
    card.append(top, grid);
    list.appendChild(card);
  });
  syncRegions();
}

for(const id of ['backgroundMode','backgroundLevel','backgroundMinimumLevel','backgroundPadding']) $(id).addEventListener('change',updateReady);
