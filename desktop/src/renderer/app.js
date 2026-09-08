'use strict';

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
  view.clear();
  $('exportResult').hidden = true;
  for (const id of ['counters', 'gates', 'histogram']) $(id).replaceChildren();
  $('legend').hidden = true;
}
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
    method: state.method,
    geometryPath: state.geometryPath,
    outputDirectory: state.outputDirectory,
    ...importSettings()
  };
  if (state.method === 'hybrid') {
    return {
      ...base,
      maxLevel: Number($('maxLevel').value),
      minimumLevel: Number($('minimumLevel').value),
      boundaryLevel: Number($('boundaryLevel').value),
      nLayers: Number($('nLayers').value),
      firstThickness: Number($('firstThickness').value),
      growthRatio: Number($('growthRatio').value),
      domainPadding: Number($('domainPadding').value),
      extrusionThickness: Number($('extrusionThickness').value)
    };
  }
  return {
    ...base,
    smallAlpha: Number($('smallAlpha').value),
    farFieldSpans: Number($('farFieldSpans').value),
    wallCellsPerSpan: Number($('wallCellsPerSpan').value),
    cellsPerLevel: Number($('cellsPerLevel').value),
    farLevel: Number($('farLevel').value),
    curvatureCellsPerRadius: $('useCurvature').checked ? Number($('curvatureCellsPerRadius').value) : 0,
    gapCells: $('useGap').checked ? Number($('gapCells').value) : 0,
    allowUnsafeWallLevel: $('allowUnsafe').checked,
    wake: $('useWake').checked ? {
      angleOfAttackDeg: Number($('wakeAngle').value),
      downstreamSpans: Number($('wakeLength').value),
      halfWidthSpans: Number($('wakeHalfWidth').value),
      levelsBelowWall: Number($('wakeLevels').value)
    } : null,
    refineBoxes: state.regions.map(({ xmin, xmax, ymin, ymax, levelsBelowWall }) =>
      ({ xmin, xmax, ymin, ymax, levelsBelowWall }))
  };
}

// Level = ceil(log2((1 + 2*far) * wallCells)) — the body span cancels, so the whole
// far-field-versus-wall-resolution trade is arithmetic and can be shown live.
function updateBudget() {
  const method = state.catalog.methods[state.method];
  if (!method || !method.supports.sizeField) return;
  const far = Number($('farFieldSpans').value);
  const cells = Number($('wallCellsPerSpan').value);
  const ceiling = method.safeWallLevel;
  const level = Math.max(0, Math.ceil(Math.log2((1 + 2 * far) * cells)));
  const over = level > ceiling && !$('allowUnsafe').checked;
  const maxCells = Math.pow(2, Math.floor(Math.log2(Math.pow(2, ceiling) / (1 + 2 * far))));
  const maxFar = (Math.pow(2, ceiling) / cells - 1) / 2;
  $('budget').className = `budget${over ? ' over' : ''}`;
  $('budget').innerHTML = over
    ? `需要树深 <b>level ${level}</b>，超过实测安全上限 ${ceiling}。<br>` +
      `远场降到 <em>${maxFar.toFixed(1)}</em> 倍，或壁面降到 <em>体长/${maxCells}</em>。`
    : `树深 <b>level ${level}</b> / 上限 ${ceiling}　计算域 <b>${(1 + 2 * far).toFixed(0)}</b> 倍体长<br>` +
      `壁面单元 <b>体长/${cells}</b>　每级带宽 ${$('cellsPerLevel').value} 格`;
  $('generate').disabled = (!state.geometryPath || state.busy || state.geometryLoading) || ($('controlMode').value === 'manual' && over);
}

function updateReady() {
  if (state.method === 'hybrid') {
    $('generate').disabled = (!state.geometryPath || state.busy || state.geometryLoading);
    return;
  }
  updateBudget();
}

function renderMethods() {
  const container = $('methods');
  container.innerHTML = '';
  for (const method of Object.values(state.catalog.methods)) {
    const button = document.createElement('button');
    button.className = `method${method.id === state.method ? ' active' : ''}`;
    button.innerHTML =
      `<b>${method.label}<span class="tag ${method.status}">${method.statusLabel}</span></b>` +
      `<i>${method.summary}</i>`;
    button.addEventListener('click', () => selectMethod(method.id));
    container.appendChild(button);
  }
}

function selectMethod(id) {
  if (id === 'hybrid' && $('fluidRegion').value === 'interior') {
    status('内部网格使用纯 Cut-cell', '贴体边界层当前仅支持外流。');
    id = 'cutcell';
  }
  state.method = id;
  const method = state.catalog.methods[id];
  $('sizingBlock').hidden = !method.supports.sizeField;
  $('hybridBlock').hidden = method.supports.sizeField;
  $('methodNote').textContent = method.supports.sizeField
    ? `实测安全壁面层级上限 ${method.safeWallLevel}；越过要显式勾选。支持 OpenFOAM 导出。`
    : `实测上限 level ${method.safeWallLevel}。这条路径没有尺寸场，层级要直接给。`;
  $('smallAlphaField').hidden = !method.supports.sizeField;
  $('smallAlphaNote').hidden = !method.supports.sizeField;
  renderMethods();
  updateControlMode();
  updateReady();
}

function updateControlMode() {
  const automatic = $('controlMode').value === 'auto';
  $('sizingBlock').hidden = automatic || state.method !== 'cutcell';
  $('hybridBlock').hidden = automatic || state.method !== 'hybrid';
  $('smallAlphaField').hidden = automatic || state.method !== 'cutcell';
  $('smallAlphaNote').hidden = automatic || state.method !== 'cutcell';
  $('densityField').hidden = !automatic;
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
  state.geometryPath = path;
  state.geometryLabel = label;
  if (sample) {
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
  if ($('fluidRegion').value === 'interior' && state.method === 'hybrid') selectMethod('cutcell');
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
    state.frame = preview.frame;
    view.setOutline(preview.loops);
    syncRegions();
    $('empty').hidden = true;
    $('legend').hidden = true;
    const spanX = Math.max(...preview.loops.flat().map(p => p[0])) - Math.min(...preview.loops.flat().map(p => p[0]));
    const spanY = Math.max(...preview.loops.flat().map(p => p[1])) - Math.min(...preview.loops.flat().map(p => p[1]));
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
    const v = probe.values;
    const rows = [
      ['体长', Number(v.size_field_body_span).toPrecision(5)],
      ['计算域跨度', Number(v.size_field_domain_span).toPrecision(5)],
      ['壁面单元尺寸', Number(v.size_field_wall_cell_size).toPrecision(4)],
      ['壁面要求层级', v.size_field_wall_level],
      ['曲率要求层级', v.size_field_curvature_level],
      ['间隙要求层级', v.size_field_proximity_level],
      ['最终树深', v.size_field_max_level],
      ['距离带 / 分段带 / 尾迹箱',
        `${v.sizing_distance_bands} / ${v.sizing_segment_bands} / ${v.sizing_box_regions}`]
    ];
    $('probeResult').innerHTML = rows
      .map(([key, value]) => `<span><i class="k">${key}</i> ${value}</span>`).join('') +
      (probe.ok ? '<span class="ok">可行</span>'
                : (probe.issues || []).map(issue => `<span class="bad">${issue}</span>`).join(''));
    $('probeResult').hidden = false;
    status(probe.ok ? '尺寸场可行' : '尺寸场被拒绝',
      `壁面 ${v.size_field_wall_level} / 曲率 ${v.size_field_curvature_level} / 间隙 ${v.size_field_proximity_level}`);
  } catch (error) {
    $('probeResult').innerHTML = `<span class="bad">${error.message}</span>`;
    $('probeResult').hidden = false;
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

// The gates are independent by design; collapsing them into one verdict is what made
// "checkMesh OK but Q1 FAIL" unreportable before.
function renderGates(result) {
  const parts = [];
  parts.push(gateRow('内部拓扑检查', result.gates.topology.pass === null ? '未确认' : result.gates.topology.pass ? 'PASS' : 'FAIL',
    result.gates.topology.pass ? '生成器内部检查通过；外部 checkMesh 需另外执行。' : '本次未完整成功，不能据此确认通过。'));

  const solver = result.gates.solver;
  if (solver) {
    const worst = solver.rows.filter(row => !row.pass);
    parts.push(gateRow('Solver 质量', solver.valid ? 'PASS' : 'FAIL',
      worst.length
        ? worst.map(row => `${row.label} ${row.value.toPrecision(4)} (限 ${row.limit})`).join('　')
        : solver.rows.map(row => `${row.label} ${row.value.toPrecision(4)}`).join('　')));
  }

  const contract = result.gates.contract;
  if (contract) {
    const detail = contract.byType
      .map(row => `${row.label} ${row.status}${row.hard ? `(${row.hard} hard)` : ''}`).join('　');
    parts.push(gateRow('Q1 合同', contract.status, detail));
    parts.push(`<div class="gate"><b></b><span class="detail">` +
      `Q1 比 solver 门和 checkMesh 都严，是诊断而不是放行条件。分类型计数在 solver ` +
      `凸划分之后统计，所以归到 cartesian 的项可能含划分碎片。</span></div>`);
  }

  if (result.actualMethod === 'cutcell-fallback') {
    parts.push(gateRow('方法', 'WARN', '请求了贴体边界层，实际落到纯 Cut-cell fallback'));
  }
  if (!result.openFoam.written) {
    parts.push(gateRow('OpenFOAM 导出', 'WARN',
      '未写出。solver 门未通过，或导出阶段自身失败；CM2D 与 VTK 仍然可用'));
  }
  $('gates').innerHTML = parts.join('');
}

const gateRow = (label, verdict, detail) =>
  `<div class="gate"><b>${label}</b><span class="verdict ${verdict}">${verdict}</span>` +
  `<span class="detail">${detail}</span></div>`;

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
    label.textContent = `L${row.level}`;
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
    view.setMesh(payload.mesh);
    syncRegions();
    $('empty').hidden = true;
    renderLegend(payload.mesh, payload.levelBasis);
    renderCounters(payload.result);
    renderGates(payload.result);
    renderHistogram(payload.levelHistogram, payload.mesh, payload.levelBasis);
    fitOverview();
    if (payload.automatic) {
      const job = payload.job;
      $('autoNote').textContent = `${payload.densityReduced ? '更密参数未通过，已降至可生成的密度。' : ''}本次采用：${job.method === 'cutcell'
        ? `壁面体长/${job.sizeField.wallCellsPerSpan}，远场 ${job.sizeField.farFieldSpans} 倍，α ${job.smallAlpha}`
        : `余域 / 壁面 level ${job.maxLevel} / ${job.boundaryLevel}，${job.nLayers} 层，首层 ${job.firstThickness}`}。实际参数与尝试记录随结果包保存。`;
    }
    $('exportResult').hidden = Boolean(payload.incomplete);
    const seconds = payload.result.timings.total_seconds;
    if (payload.incomplete) {
      status('网格已生成，后续步骤失败', payload.incomplete);
    } else {
      status('生成完成', `${payload.automatic ? `自动选参成功（第 ${payload.attempts.length} 组）` : '手动生成成功'}；可导出结果包${seconds ? `　${seconds.toFixed(2)} s` : ''}`);
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

$('pickGeometry').addEventListener('click', async () => {
  const picked = await window.cartmesh.pickGeometry();
  if (!picked) return;
  $('sample').value = '';
  await chooseGeometry(picked, picked.split('/').pop(), null);
});

$('sample').addEventListener('change', async event => {
  const sample = state.catalog.samples.find(item => item.id === event.target.value);
  if (!sample) return;
  await chooseGeometry(sample.path, sample.label, sample);
});

for (const id of ['farFieldSpans', 'wallCellsPerSpan', 'cellsPerLevel', 'farLevel']) {
  $(id).addEventListener('input', updateBudget);
}
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
window.addEventListener('resize', () => view.draw());

(async () => {
  state.catalog = await window.cartmesh.catalog();
  renderMethods();
  renderPresets();
  renderSamples();
  selectMethod('cutcell');
  renderRegions();
  // Smoke tests drive these same handlers; an optional output override retains fixtures.
  window.__smoke = { state, selectMethod, chooseGeometry, generate, setOutput, addRegion, renderRegions, view };
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
  view.frame = state.frame;
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
