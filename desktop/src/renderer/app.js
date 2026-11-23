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
  adaptiveFieldPath: null,
  adaptiveSuggestion: null,
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
  $('suggestFlowRefinement').disabled = busy || !state.result || !state.flow || state.method === 'background';
  updateFlowMode();
  updateReady();
}
function validInputs() {
  for (const input of document.querySelectorAll('.panel input[type=number]')) {
    if (!input.closest('#flowBlock, #thermalBlock, #yplusDetails') && !input.disabled && input.getClientRects().length &&
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
  state.adaptiveFieldPath = null;
  state.adaptiveSuggestion = null;
  updateFlowAdaptNote();
  $('flowTimeline').hidden=true;
  if (hidePanel) {state.flowRestart=null;state.flowBoundaryDefinition=null;view.setBoundaryHighlight([]);renderFlowBoundaries();updateFlowMode();
    $('flowCaseFileInfo').textContent='保存从零起算的流动设置，读取前须生成同一最终网格；续算请使用重启文件。温度设置单独配置。';}
  view.setFlowFields(null);
  view.setAdaptiveFields(null);
  $('flowResult').hidden = true;
  $('flowResult').replaceChildren();
  $('suggestFlowRefinement').disabled = true;
  $('flowSpeedOption').hidden = true;
  $('flowPressureOption').hidden = true;
  $('adaptiveIndicatorOption').hidden = true;
  if (['speed', 'pressure', 'adaptive'].includes($('displayMode').value)) {
    $('displayMode').value = 'level'; view.mode = 'level';
    $('canvasWrap').classList.remove('light');
    view.draw();
  }
  if (hidePanel) $('flowBlock').hidden = true;
  if (state.mesh) renderLegend(state.mesh, state.levelBasis);
}
function clearResult() {
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

let previewSequence = 0;
const view = new window.MeshView.Viewport($('canvas'));
const {
  renderCounters, renderGates, renderHistogram, renderLegend: renderLegendView,
  renderResolution, renderFlowResult: renderFlowResultView, renderThermalResult, renderTimeSeriesMonitor
} = window.CartMeshResultsView;
const renderLegend = (mesh, basis) => renderLegendView(mesh, basis, view);
const renderFlowResult = summary => renderFlowResultView(summary, state.catalog);

const fmt = value => Number(value || 0).toLocaleString('en-US');
const log = line => { $('log').textContent += `${line}\n`; $('log').scrollTop = 1e9; };
const status = (title, text) => { $('statusTitle').textContent = title; $('statusText').textContent = text || ''; };

function importSettings() {
  const rawExclude = $('dxfExcludeLayers')?.value || '';
  const dxfExcludeLayers = rawExclude.split(/[,;\n]+/).map(s => s.trim()).filter(Boolean);
  return {
    chordError: Number($('chordError').value),
    sourceUnits: $('sourceUnits').value,
    fluidRegion: $('fluidRegion').value,
    dxfWeldTolerance: Number($('dxfWeldTolerance')?.value || 1e-10),
    dxfIgnoreAnnotations: Boolean($('dxfIgnoreAnnotations')?.checked),
    dxfExcludeLayers
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
  const useAdaptiveCellField = !base.automatic && !base.targetCells &&
    Boolean(state.adaptiveFieldPath) && $('flowAdaptMode').value === 'field';
  const adaptiveRequest = {
    useAdaptiveCellField,
    adaptiveFieldPath: useAdaptiveCellField ? state.adaptiveFieldPath : undefined,
    refineBoxes: (base.targetCells || useAdaptiveCellField ? [] : state.regions).map(({ xmin, xmax, ymin, ymax, levelsBelowWall }) =>
      ({ xmin, xmax, ymin, ymax, levelsBelowWall }))
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
      extrusionThickness: Number($('extrusionThickness').value),
      ...adaptiveRequest
    };
  }
  return {
    ...base,
    smallAlpha: Number($('smallAlpha').value),
    cutcellStrategy: $('cutcellStrategy').value,
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
    ...adaptiveRequest
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
function applyManualParameters(request) {
  if (!request) return;
  const fields = { wallRelativeSize:'wallRelativeSize', backgroundRelativeSize:'backgroundRelativeSize',
    farFieldSpans:'relativePadding', cellsPerLevel:'relativeBandCells', firstLayerRelativeSize:'firstLayerRelativeSize',
    nLayers:'nLayers', growthRatio:'growthRatio', extrusionRelativeSize:'extrusionRelativeSize', smallAlpha:'smallAlpha', cutcellStrategy:'cutcellStrategy' };
  for (const [key,id] of Object.entries(fields)) if (request[key] != null) $(id).value = request[key];
  if (request.referenceLength != null) {
    $('referenceMode').value = 'explicit'; $('referenceLength').value = request.referenceLength;
    $('referenceLengthField').hidden = false;
    $('referenceMode').closest('details').open = true;
  } else {
    $('referenceMode').value = 'bbox';
    $('referenceLengthField').hidden = true;
  }
  $('relativeAllowUnsafe').checked = Boolean(request.allowUnsafeWallLevel);
  $('controlMode').value = 'manual'; updateControlMode(); syncRegions();
}
function useManualParameters(request) {
  if (!request || state.busy) return;
  applyManualParameters(request);
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
  $('regionsBlock').hidden = automatic || background;
  $('cutcellStrategyField').hidden = state.method !== 'cutcell';
  $('smallAlphaField').hidden = automatic || state.method !== 'cutcell' || $('cutcellStrategy').value === 'shape';
  $('smallAlphaNote').hidden = $('smallAlphaField').hidden;
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
    const esc = value => String(value ?? '').replace(/[&<>"']/g, ch => ({ '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;' }[ch]));
    let dxfExtra = '';
    if (preview.dxfReport) {
      const r = preview.dxfReport;
      const parts = [];
      if (r.effective_units_name) parts.push(`单位 <b>${esc(r.effective_units_name)} → m</b>`);
      if (r.insert_entity_count > 0) parts.push(`INSERT <b>${r.insert_entity_count}</b>`);
      if (r.polyline_entity_count > 0) parts.push(`POLYLINE <b>${r.polyline_entity_count}</b>`);
      if (r.ignored_annotation_count > 0) parts.push(`已忽略标注 <b>${r.ignored_annotation_count}</b>`);
      if (parts.length) dxfExtra += `<span>${parts.join('　')}</span>`;
      const summaries = Array.isArray(r.loop_summaries) && r.loop_summaries.length
        ? r.loop_summaries
        : (Array.isArray(r.embedded_patches) ? r.embedded_patches.map((p, i) => ({ loop_id: i, patch_name: p.name, role: p.role, source_layer: p.source_layer })) : []);
      if (summaries.length) {
        const badges = summaries.slice(0, 12).map(s =>
          `${esc(s.patch_name)}(${esc(s.role)})`).join('、') +
          (summaries.length > 12 ? ` 等 ${summaries.length} 个` : '');
        dxfExtra += `<span>边界区 <b>${badges}</b></span>`;
      }
    }
    $('geometryFacts').innerHTML =
      `<span>格式 <b>${esc(preview.kind)}</b>　环 <b>${preview.loops.length}</b>　` +
      `顶点 <b>${preview.loops.reduce((sum, loop) => sum + loop.length, 0)}</b></span>` +
      `<span>包围盒 <b>${spanX.toPrecision(4)} × ${spanY.toPrecision(4)} m</b></span>` + dxfExtra;
    state.dxfReport = preview.dxfReport || null;
    scheduleYPlus();
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

// Both of these failures are alignment events: a wall vertex grazing a grid line.
// Measured on a 48-vertex CSV outline, the solver gate fails at a 8x and 9x far field
// and passes at 10x and 11x with the wall resolution untouched, so nudging the far
// field is a real lever and not folklore.
function advice(message) {
  if (/topology audit|UnclassifiedBoundaryEdge/.test(message)) {
    return '提示：这是壁面顶点擦过格点的构造失败（见 CURRENT_STATE_CN.md 第 4 节 A）。' +
      '把远场距离改动 ±1 倍体长会移动整个格点阵列，通常比调壁面分辨率更有效。';
  }
  return '';
}

async function generate() {
  if (state.busy || state.geometryLoading || !state.geometryPath || !validInputs()) return;
  const request = buildRequest();
  clearResult();
  setBusy(true);
  $('generate').textContent = '正在生成…';
  $('log').textContent = '';
  status('生成中', state.method==='background'?'几何诊断 → 完整笛卡尔网格 → 分类 → 导出':'几何转换 → 尺寸场 → 加密 → cut-cell → 稳定化 → 质量');
  try {
    const payload = await window.cartmesh.generate(request);
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
    renderResolution(payload.result.resolution);
    fitOverview();
    if (payload.automatic) {
      const job = payload.job;
      $('autoNote').textContent = job.relativeSizing
        ? `${payload.cellBudget ? '参数已按最终选中的结果记录。' : payload.densityReduced ? '更密请求未满足，已降密。' : ''}上次生成的壁面目标 h/Lref=${job.sizeField.wallRelativeSize.toPrecision(4)}，背景 h/Lref=${job.sizeField.backgroundRelativeSize.toPrecision(4)}，留白 ${job.sizeField.farFieldSpans.toPrecision(4)} Lref。实际参数与尝试记录随结果包保存。`
        : `${payload.densityReduced ? '更密参数未通过，已降至可生成的密度。' : ''}本次采用：${job.method === 'cutcell'
        ? `壁面目标体长/${job.sizeField.wallCellsPerSpan}，每级带宽 ${job.sizeField.cellsPerLevel} 格，远场 ${job.sizeField.farFieldSpans} 倍，${job.cutcellStrategy === 'shape' ? '形状聚合＋局部重划分' : 'α ' + job.smallAlpha}`
        : `余域 / 壁面 level ${job.maxLevel} / ${job.boundaryLevel}，${job.nLayers} 层，首层 ${job.firstThickness}`}。实际参数与尝试记录随结果包保存。`;
    }
    $('exportResult').hidden = Boolean(payload.incomplete);
    $('flowBlock').hidden = Boolean(payload.incomplete || payload.background);
    $('thermalBlock').hidden = Boolean(payload.incomplete || payload.background);
    if (!payload.incomplete && !payload.background) {
      $('flowCase').value = payload.job.fluidRegion === 'interior' ? 'duct' : 'external';
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
  $('flowIterationLabel').textContent = transient ? '每个时间步的外迭代上限' : '最大外迭代';
  const pimple=$('flowCoupling').value==='pimple';
  $('flowMomentumCorrectorsField').hidden=!pimple;
  $('flowMomentumCorrectors').disabled=state.busy||!pimple;
  $('flowCoupling').disabled=state.busy;
  $('flowPressureRelaxation').disabled=state.busy;
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
  renderTimeSeriesMonitor(state.flowHistory, {
    timelineId: 'flowTimeline', svgId: 'flowMonitor', metricId: 'flowMonitorMetric',
    captionId: 'flowMonitorCaption', prefix: '已接受时间步'
  });
}
function renderThermalMonitor() {
  renderTimeSeriesMonitor(state.thermalHistory, {
    timelineId: 'thermalTimeline', svgId: 'thermalMonitor', metricId: 'thermalMonitorMetric',
    captionId: 'thermalMonitorCaption', prefix: '已接受联合时间步'
  });
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
  if($('flowCoupling').value!=='simple' || Number($('flowPressureRelaxation').value)!==.25 || $('flowLinearPolicy').value!=='strict' || Number($('flowVelocityRelaxation').value)!==.6 || Number($('flowPressureCorrections').value)!==4) {
    status('温度联算的迭代设置需要调整','请恢复 SIMPLE、0.25压力松弛、固定线性精度、0.6速度松弛和4次压力校正。');return false;
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
  renderThermalResult(payload.summary, state.thermalHistory.at(-1));
  renderThermalMonitor();
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
  $('suggestFlowRefinement').disabled = state.busy || !state.result || state.method === 'background';
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
    coupling:$('flowCoupling').value,pressureRelaxation:Number($('flowPressureRelaxation').value),
    momentumCorrectors:$('flowCoupling').value==='pimple'?Number($('flowMomentumCorrectors').value):2,
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
    coupling:'flowCoupling',pressureRelaxation:'flowPressureRelaxation',momentumCorrectors:'flowMomentumCorrectors',
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

async function runFlow() {
  if (state.busy || !state.result || !state.mesh || !validFlowInputs()) return;
  const request=flowRequest(),transient=request.mode!=='steady';
  clearFlowBinding();setBusy(true);$('runFlow').textContent='正在求解…';
  status('层流求解中',transient?'按物理时间推进；取消后可从最后接受的时间步继续。':`${request.coupling.toUpperCase()} 速度—压力耦合；可随时取消。`);
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

function updateFlowAdaptNote(overrideMessage = '') {
  if (overrideMessage) {
    $('flowAdaptNote').textContent = overrideMessage;
    return;
  }
  const s = state.adaptiveSuggestion;
  if (!s) {
    $('flowAdaptNote').textContent = '完成层流预览后可用；自动继承本次参数并展示梯度热图与参考框。';
    return;
  }
  const pct = (s.selectedAreaFraction * 100).toFixed(1);
  const thermalSuffix = s.referenceScalarSpan > 0 ? '（含温度梯度）' : '';
  const fieldActive = Boolean(state.adaptiveFieldPath) && $('flowAdaptMode').value === 'field';
  if (fieldActive) {
    $('flowAdaptNote').textContent =
      `已就绪紧凑单元级流场 AMR${thermalSuffix}（${s.cellTargetCount || s.selectedCellCount} 个高梯度目标单元，选中 ${s.selectedCellCount}/${s.cellCount}，面积占比 ${pct}%；下方 ${state.regions.length} 个参考框仅用于预览，手动编辑将自动切换为矩形加密框模式），点击“生成预览”即可重划。`;
  } else {
    $('flowAdaptNote').textContent =
      `已启用可编辑矩形加密框模式${thermalSuffix}（共 ${state.regions.length} 个加密区；原流场选中 ${s.selectedCellCount}/${s.cellCount} 单元，面积占比 ${pct}%），微调后点击“生成预览”即可局部重划。`;
  }
}

function switchToBoxAdaptOnManualEdit() {
  if ($('flowAdaptMode').value === 'field' && state.adaptiveFieldPath) {
    $('flowAdaptMode').value = 'boxes';
    updateFlowAdaptNote('检测到手动编辑加密区，已自动切换为“可编辑矩形加密框”模式；点击“生成预览”将按下方矩形框局部重划。');
  }
}

async function suggestFlowRefinement() {
  if (state.busy || !state.result || !state.flow || state.method === 'background') return;
  const refineFraction = Number($('flowAdaptFraction').value);
  const maxBoxes = Number($('flowAdaptMaxBoxes').value);
  setBusy(true);
  status('流场自适应分析中', '正在按当前速度与压力梯度提取局部加密区…');
  try {
    const suggestion = await window.cartmesh.suggestFlowRefinement({ refineFraction, maxBoxes });
    if (suggestion.fields?.cells?.length) {
      view.setAdaptiveFields(suggestion.fields.cells);
      $('adaptiveIndicatorOption').hidden = false;
      $('displayMode').value = 'adaptive';
      view.mode = 'adaptive';
      $('canvasWrap').classList.remove('light');
      view.draw();
      if (state.mesh) renderLegend(state.mesh, state.levelBasis);
    }
    if (!suggestion.regions.length && !(suggestion.cellTargetCount > 0)) {
      state.adaptiveFieldPath = null;
      state.adaptiveSuggestion = null;
      $('flowAdaptNote').textContent = '当前流场梯度接近均匀，未生成新的局部加密区。';
      status('无需局部加密', '当前流场未检测到显著高梯度区。');
      return;
    }
    if ($('controlMode').value !== 'manual') {
      if (state.selectedRequest) {
        applyManualParameters(state.selectedRequest);
      } else {
        $('controlMode').value = 'manual';
        updateControlMode();
      }
    }
    state.adaptiveFieldPath = suggestion.cellTargetCount > 0 && suggestion.adaptiveFieldPath
      ? suggestion.adaptiveFieldPath : null;
    state.adaptiveSuggestion = suggestion;
    state.regions = suggestion.regions.map(({ xmin, xmax, ymin, ymax, levelsBelowWall }) => ({
      xmin, xmax, ymin, ymax, levelsBelowWall, active: true
    }));
    renderRegions();
    const thermalSuffix = suggestion.referenceScalarSpan > 0 ? '（含温度梯度）' : '';
    const fieldActive = Boolean(state.adaptiveFieldPath) && $('flowAdaptMode').value === 'field';
    updateFlowAdaptNote();
    status('自适应加密已就绪', fieldActive
      ? `已按梯度${thermalSuffix}生成 ${suggestion.cellTargetCount} 个单元级 AMR 目标并显示梯度热图，点击“生成预览”即可紧凑重划。`
      : `已按梯度${thermalSuffix}填入 ${state.regions.length} 个局部加密框并显示梯度热图，点击“生成预览”即可局部重划。`);
  } catch (error) {
    const message = error.message.replace(/^Error invoking remote method '[^']+': Error: /, '');
    status('自适应分析失败', message.split('\n')[0]);
    log(message);
  } finally {
    setBusy(false);
  }
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
  $(id).addEventListener('input', () => { updateBudget(); syncRegions(); if (id === 'referenceLength' || id === 'relativePadding') scheduleYPlus(); });
$('relativeAllowUnsafe').addEventListener('change', () => { updateBudget(); schedulePlan(); });
$('referenceMode').addEventListener('change', () => {
  $('referenceLengthField').hidden = $('referenceMode').value !== 'explicit';
  updateBudget();
  schedulePlan();
  scheduleYPlus();
  syncRegions();
});
for (const id of ['chordError', 'sourceUnits', 'dxfWeldTolerance', 'dxfIgnoreAnnotations', 'dxfExcludeLayers']) {
  $(id)?.addEventListener('change', () => {
    if (state.geometryPath && /\.dxf$/i.test(state.geometryPath)) drawGeometryOutline();
  });
}

let yplusSequence = 0;
function activeYPlusReferenceLength() {
  if ($('referenceMode').value === 'explicit') {
    const val = Number($('referenceLength').value);
    if (val > 0) return val;
  }
  return state.frame?.bodySpan > 0 ? state.frame.bodySpan : 1.0;
}

function updateYPlusWakeControls() {
  const exterior = $('fluidRegion')?.value !== 'interior';
  const syncWake = Boolean($('yplusSyncWake')?.checked) && exterior;
  if ($('yplusSyncWakeToggle')) $('yplusSyncWakeToggle').hidden = !exterior;
  if ($('yplusWakeFields')) $('yplusWakeFields').hidden = !syncWake;
  if ($('applyYPlus')) {
    $('applyYPlus').textContent = syncWake ? '一键填入边界层与尾迹区' : '填入边界层参数';
  }
}

async function scheduleYPlus() {
  if (!$('yplusPreview') || !window.cartmesh?.estimateYPlus) return;
  updateYPlusWakeControls();
  const sequence = ++yplusSequence;
  const Lref = activeYPlusReferenceLength();
  const fluidRegion = $('fluidRegion')?.value === 'interior' ? 'interior' : 'exterior';
  const currentFar = Number($('relativePadding')?.value);
  try {
    const est = await window.cartmesh.estimateYPlus({
      referenceLength: Lref,
      velocity: Number($('yplusVelocity').value),
      nu: Number($('yplusNu').value),
      targetYPlus: Number($('yplusTarget').value),
      growthRatio: Number($('yplusGrowth').value),
      coverageRatio: Number($('yplusCoverage')?.value ?? 1),
      regime: $('yplusRegime').value,
      fluidRegion,
      frame: state.frame || undefined,
      farFieldSpans: !$('yplusSyncWall')?.checked && currentFar > 0 ? currentFar : undefined,
      wakeLengthSpans: Number($('yplusWakeLength')?.value ?? 2),
      wakeLevelsBelowWall: Number($('yplusWakeLevelsBelow')?.value ?? 2),
      wakeAngleDeg: Number($('yplusWakeAngle')?.value ?? 0)
    });
    if (sequence !== yplusSequence) return;
    state.yplusEstimate = est;
    const regimeLabel = est.effectiveRegime === 'turbulent' ? '湍流' : '层流';
    const lines = [
      `Lref = ${est.referenceLength.toPrecision(4)} m，Re = ${est.reynolds.toExponential(3)}（${regimeLabel}，Cf = ${est.cf.toExponential(3)}）`,
      `首层总高 Δy₁ = ${est.firstLayerHeight.toExponential(3)} m（h₁/Lref = ${est.firstLayerRelativeSize.toExponential(3)}，单元中心 y_P = ${est.firstCellCenterDistance.toExponential(3)} m）`,
      `建议 ${est.nLayers} 层（增长率 ${est.growthRatio}）→ 总厚 ${(est.totalThickness * 1000).toPrecision(4)} mm ≈ ${(est.coverageAchieved * 100).toFixed(0)}% δ₉₉（δ₉₉ = ${(est.delta99 * 1000).toPrecision(4)} mm）`,
      `最外层高 ${(est.lastLayerHeight * 1000).toPrecision(4)} mm；推荐壁面格 h/Lref ≈ ${est.recommendedWallRelativeSize.toPrecision(4)}（首层展弦比 ≈ ${est.estimatedFirstLayerAspectRatio.toFixed(0)}）`
    ];
    if (est.recommendedWakeBox && $('yplusSyncWake')?.checked) {
      const wb = est.recommendedWakeBox;
      const farHint = est.recommendedFarFieldSpans ? `，外流远场 ≥ ${est.recommendedFarFieldSpans} Lref` : '';
      lines.push(`联动尾迹框（--refine-box）：x ∈ [${wb.xmin}, ${wb.xmax}]，y ∈ [${wb.ymin}, ${wb.ymax}] Lref，低于壁面 ${wb.levelsBelowWall} 级${farHint}`);
    } else if (fluidRegion === 'interior') {
      lines.push('内流通道模式不生成外部尾迹加密框。');
    }
    if (est.warnings.length) lines.push(`提示：${est.warnings.join('；')}`);
    $('yplusPreview').className = `budget${est.warnings.length ? ' over' : ''}`;
    $('yplusPreview').textContent = lines.join('\n');
    $('applyYPlus').disabled = Boolean(state.busy);
  } catch (error) {
    if (sequence !== yplusSequence) return;
    state.yplusEstimate = null;
    $('yplusPreview').className = 'budget over';
    $('yplusPreview').textContent = error.message.replace(/^Error invoking remote method.*?Error: /, '');
    $('applyYPlus').disabled = true;
  }
}

$('yplusFluid')?.addEventListener('change', () => {
  const preset = state.catalog?.fluidPresets?.[$('yplusFluid').value];
  if (preset) $('yplusNu').value = preset.nu.toPrecision(6);
  scheduleYPlus();
});
$('yplusNu')?.addEventListener('input', () => {
  if ($('yplusFluid')) $('yplusFluid').value = 'custom';
  scheduleYPlus();
});
for (const id of ['yplusVelocity', 'yplusTarget', 'yplusGrowth', 'yplusCoverage', 'yplusRegime', 'yplusSyncWall', 'yplusSyncWake', 'yplusWakeLength', 'yplusWakeLevelsBelow', 'yplusWakeAngle']) {
  $(id)?.addEventListener('input', scheduleYPlus);
  $(id)?.addEventListener('change', scheduleYPlus);
}
function applyYPlusToForm() {
  const est = state.yplusEstimate;
  if (!est || state.busy) return false;
  $('controlMode').value = 'manual';
  updateControlMode();
  $('nLayers').value = String(est.nLayers);
  $('firstLayerRelativeSize').value = Number(est.firstLayerRelativeSize.toPrecision(6));
  $('growthRatio').value = String(est.growthRatio);
  const exterior = $('fluidRegion')?.value !== 'interior';
  const applyWake = Boolean($('yplusSyncWake')?.checked) && exterior && Boolean(est.recommendedWakeBox);
  if ($('yplusSyncWall')?.checked) {
    $('wallRelativeSize').value = Number(est.recommendedWallRelativeSize.toPrecision(5));
    $('backgroundRelativeSize').value = Number(est.recommendedBackgroundRelativeSize.toPrecision(5));
    if (exterior && est.recommendedFarFieldSpans) {
      const curFar = Number($('relativePadding').value);
      if (!Number.isFinite(curFar) || curFar < est.recommendedFarFieldSpans) {
        $('relativePadding').value = String(est.recommendedFarFieldSpans);
      }
    }
    const span = state.frame?.bodySpan || est.referenceLength;
    const far = Number($('relativePadding').value);
    const domain = span + 2 * far * est.referenceLength;
    const depth = Math.max(0, Math.ceil(Math.log2(domain / (est.recommendedWallRelativeSize * est.referenceLength))));
    const ceiling = state.catalog?.methods?.[state.method]?.safeWallLevel ?? 8;
    if (depth > ceiling) $('relativeAllowUnsafe').checked = true;
  }
  if (applyWake) {
    switchToBoxAdaptOnManualEdit();
    const retained = state.regions
      .filter(box => box && box.source !== 'yplus-wake')
      .map(box => ({ ...box, active: false }));
    retained.push({
      xmin: est.recommendedWakeBox.xmin,
      xmax: est.recommendedWakeBox.xmax,
      ymin: est.recommendedWakeBox.ymin,
      ymax: est.recommendedWakeBox.ymax,
      levelsBelowWall: est.recommendedWakeBox.levelsBelowWall,
      source: 'yplus-wake',
      active: true
    });
    state.regions = retained;
    renderRegions();
  } else if (!exterior && state.regions.some(box => box && box.source === 'yplus-wake')) {
    state.regions = state.regions.filter(box => box && box.source !== 'yplus-wake');
    renderRegions();
  } else {
    syncRegions();
  }
  updateBudget();
  const wakeNote = applyWake
    ? `，已联动尾迹框 x∈[${est.recommendedWakeBox.xmin},${est.recommendedWakeBox.xmax}]`
    : '';
  status(applyWake ? '已按 y⁺ 填入边界层与尾迹区' : '已按 y⁺ 填入边界层参数',
    `Re=${est.reynolds.toExponential(2)}，y⁺=${est.targetYPlus}，首层 h₁/Lref=${est.firstLayerRelativeSize.toExponential(3)}，${est.nLayers} 层${wakeNote}。`);
  return true;
}
$('applyYPlus')?.addEventListener('click', applyYPlusToForm);
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
$('suggestFlowRefinement').addEventListener('click', suggestFlowRefinement);
$('flowAdaptMode').addEventListener('change', () => updateFlowAdaptNote());
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
$('flowCoupling').addEventListener('change', () => { if (state.flow) clearFlowBinding(); updateFlowMode(); });
for(const id of ['flowPressureRelaxation','flowMomentumCorrectors'])
  $(id).addEventListener('input', () => { if (state.flow) clearFlowBinding(); });
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
$('cutcellStrategy').addEventListener('change', () => { clearResult(); updateControlMode(); });
$('fluidRegion').addEventListener('change', () => {
  clearResult();
  selectMethod(state.method);
  scheduleYPlus();
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
  scheduleYPlus();
  // Smoke tests drive these same handlers; an optional output override retains fixtures.
  window.__smoke = { state, selectMethod, chooseGeometry, generate, runFlow, suggestFlowRefinement, runThermal, flowRequest, saveFlowCase, loadFlowCase, prepareFlowBoundaries, setOutput, addRegion, renderRegions, view, loadVerifiedPreset, setSidebarCollapsed, returnToStart, importGeometryFile, scheduleYPlus, applyYPlusToForm };
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
  switchToBoxAdaptOnManualEdit();
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
    name.textContent = box.source === 'yplus-wake' ? `R${index + 1} · 尾迹` : `R${index + 1}`;
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
      switchToBoxAdaptOnManualEdit();
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
      switchToBoxAdaptOnManualEdit();
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
        switchToBoxAdaptOnManualEdit();
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
