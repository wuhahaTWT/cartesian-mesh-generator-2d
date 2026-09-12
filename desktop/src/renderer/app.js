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
  updateReady();
}
function validInputs() {
  for (const input of document.querySelectorAll('.panel input[type=number]')) {
    if (!input.disabled && input.getClientRects().length &&
        (!input.value.trim() || !input.checkValidity())) {
      input.reportValidity();
      input.focus();
      status('参数需要调整', '请填写有效数值，并检查范围。');
      return false;
    }
  }
  return true;
}
function clearResult() {
  state.mesh = null; state.result = null; state.wallBounds = null;
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

let previewSequence = 0;
const view = new window.MeshView.Viewport($('canvas'));
const { levelColour, RAMP } = window.MeshView;

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
  if (!state.frame || $('controlMode').value !== 'auto') return;
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
  state.method = id;
  const method = state.catalog.methods[id];
  $('sizingBlock').hidden = !method.supports.sizeField;
  $('hybridBlock').hidden = method.supports.sizeField;
  $('methodNote').textContent = method.supports.sizeField
    ? `实测安全壁面层级上限 ${method.safeWallLevel}；越过要显式勾选。支持 OpenFOAM 导出。`
    : `贴体层和余域共用相对尺寸；构造深度由程序推导，当前验证上限 ${method.safeWallLevel}。`;
  $('smallAlphaField').hidden = !method.supports.sizeField;
  $('smallAlphaNote').hidden = !method.supports.sizeField;
  renderMethods();
  updateControlMode();
  updateReady();
}

function updateControlMode() {
  const automatic = $('controlMode').value === 'auto';
  $('resolutionBlock').hidden = automatic;
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
  $('verifiedPreset').hidden = !['circle','naca2412','nozzle'].includes(state.sampleId);
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
  const coloured = view.mode === 'level';
  const ramp = document.createElement('div');
  ramp.className = 'ramp';
  for (const colour of RAMP) {
    const swatch = document.createElement('span');
    swatch.style.background = colour;
    ramp.appendChild(swatch);
  }
  const ends = document.createElement('div');
  ends.className = 'ends';
  const prefix = basis === 'size' ? '档' : 'L';
  const coarse = document.createElement('span');
  coarse.textContent = `${prefix}${mesh.minLevel} 粗`;
  const fine = document.createElement('span');
  fine.textContent = `细 ${prefix}${mesh.maxLevel}`;
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
  status('生成中', '几何转换 → 尺寸场 → 加密 → cut-cell → 稳定化 → 质量');
  try {
    const payload = await window.cartmesh.generate(buildRequest());
    state.mesh = payload.mesh;
    state.wallBounds = payload.wallBounds;
    state.result = payload.result;
    state.levelBasis = payload.levelBasis;
    state.selectedRequest = payload.selectedRequest || null;
    state.cellBudget = payload.cellBudget || null;
    $('actualToManual').hidden = !state.selectedRequest;
    $('cellBudgetResult').hidden = !payload.cellBudget;
    if (payload.cellBudget) {
      const b = payload.cellBudget;
      $('cellBudgetResult').textContent = `目标约 ${fmt(b.targetCells)}；实际 ${fmt(b.actualCells)} 个单元。${b.reached ? '已达到目标范围（±30%）。' : '未达到目标范围，保留本次最接近目标的可导出结果。'} 共尝试 ${b.attempts} 次。` + (b.stoppedReason ? ` 后续调整停止：${b.stoppedReason}` : '');
    }
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
$('addRegion').addEventListener('click', addRegion);

$('displayMode').addEventListener('change', event => {
  view.mode = event.target.value;
  $('canvasWrap').classList.toggle('light', view.mode === 'light');
  if (state.mesh) renderLegend(state.mesh, state.levelBasis);
  view.draw();
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
window.cartmesh.onRunLine(log);
window.addEventListener('resize', () => view.requestDraw());

(async () => {
  state.catalog = await window.cartmesh.catalog();
  renderMethods();
  renderPresets();
  renderSamples();
  selectMethod('cutcell');
  renderRegions();
  // Smoke tests drive these same handlers; an optional output override retains fixtures.
  window.__smoke = { state, selectMethod, chooseGeometry, generate, setOutput, addRegion, renderRegions, view, loadVerifiedPreset, setSidebarCollapsed, returnToStart, importGeometryFile };
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
