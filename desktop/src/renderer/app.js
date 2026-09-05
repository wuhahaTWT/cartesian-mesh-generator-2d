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
  result: null
};

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
    } : null
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
  $('generate').disabled = !(state.geometryPath && state.outputDirectory) || over;
}

function updateReady() {
  if (state.method === 'hybrid') {
    $('generate').disabled = !(state.geometryPath && state.outputDirectory);
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
  state.method = id;
  const method = state.catalog.methods[id];
  $('sizingBlock').hidden = !method.supports.sizeField;
  $('hybridBlock').hidden = method.supports.sizeField;
  $('methodNote').textContent = method.supports.sizeField
    ? `实测安全壁面层级上限 ${method.safeWallLevel}；越过要显式勾选。支持 OpenFOAM 导出。`
    : `实测上限 level ${method.safeWallLevel}。这条路径没有尺寸场，层级要直接给。`;
  renderMethods();
  updateReady();
}

function applySizeField(field) {
  if (!field) return;
  if (field.farFieldSpans !== undefined) $('farFieldSpans').value = field.farFieldSpans;
  if (field.wallCellsPerSpan !== undefined) $('wallCellsPerSpan').value = field.wallCellsPerSpan;
  if (field.cellsPerLevel !== undefined) $('cellsPerLevel').value = field.cellsPerLevel;
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
    applySizeField(sample.sizeField);
    // Each sample ships the small-cell threshold it was verified at.  The passing
    // combinations are not monotone in alpha: near the level ceiling, whether a wall
    // vertex grazes a grid line is what decides the solver gate.
    if (sample.smallAlpha !== undefined) $('smallAlpha').value = sample.smallAlpha;
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
  try {
    status('正在读取几何', state.geometryLabel);
    const preview = await window.cartmesh.previewGeometry({
      geometryPath: state.geometryPath, ...importSettings()
    });
    state.mesh = null;
    view.setOutline(preview.loops);
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
    $('geometryFacts').hidden = true;
    status('几何读取失败', error.message.split('\n')[0]);
    log(error.message);
  }
}

async function probeSizing() {
  $('probe').disabled = true;
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
    $('probe').disabled = false;
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
  parts.push(gateRow('拓扑不变量', result.gates.topology.pass ? 'PASS' : 'FAIL',
    '无重复 / 孤立 / 非流形边，面积守恒。生成器 fail-closed，走到这里即已通过'));

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
function renderHistogram(histogram, mesh) {
  const container = $('histogram');
  container.replaceChildren();
  if (!histogram.length) return;
  const title = document.createElement('div');
  title.className = 'title';
  title.textContent = '每层级单元数';
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

function renderLegend(mesh) {
  const container = $('legend');
  container.replaceChildren();
  const ramp = document.createElement('div');
  ramp.className = 'ramp';
  for (const colour of RAMP) {
    const swatch = document.createElement('span');
    swatch.style.background = colour;
    ramp.appendChild(swatch);
  }
  const ends = document.createElement('div');
  ends.className = 'ends';
  const coarse = document.createElement('span');
  coarse.textContent = `L${mesh.minLevel} 粗`;
  const fine = document.createElement('span');
  fine.textContent = `细 L${mesh.maxLevel}`;
  ends.append(coarse, fine);

  const keys = document.createElement('div');
  keys.className = 'keys';
  for (const [colour, text] of [['#ff5a1f', '物面（嵌入边界）'],
                                ['rgba(150,178,196,.55)', '计算域边界']]) {
    const key = document.createElement('span');
    const dash = document.createElement('i');
    dash.style.background = colour;
    key.append(dash, document.createTextNode(text));
    keys.appendChild(key);
  }
  container.append(ramp, ends, keys);
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

async function generate() {  $('generate').disabled = true;
  $('generate').textContent = '正在生成…';
  $('log').textContent = '';
  status('生成中', '几何转换 → 尺寸场 → 加密 → cut-cell → 稳定化 → 质量');
  try {
    const payload = await window.cartmesh.generate(buildRequest());
    state.mesh = payload.mesh;
    state.wallBounds = payload.wallBounds;
    state.result = payload.result;
    view.setMesh(payload.mesh);
    $('empty').hidden = true;
    renderLegend(payload.mesh);
    renderCounters(payload.result);
    renderGates(payload.result);
    renderHistogram(payload.levelHistogram, payload.mesh);
    $('openOutput').hidden = false;
    const seconds = payload.result.timings.total_seconds;
    if (payload.incomplete) {
      status('网格已生成，后续步骤失败', payload.incomplete);
    } else {
      status('生成完成', `${payload.cm2dPath}${seconds ? `　${seconds.toFixed(2)} s` : ''}`);
    }
  } catch (error) {
    status('生成失败', error.message.split('\n')[0]);
    log(error.message);
    const hint = advice(error.message);
    if (hint) log(hint);
  } finally {
    $('generate').textContent = '生成网格';
    updateReady();
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

$('pickOutput').addEventListener('click', async () => {
  const picked = await window.cartmesh.pickOutput();
  if (!picked) return;
  state.outputDirectory = picked;
  $('outputPath').textContent = picked;
  updateReady();
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
$('openOutput').addEventListener('click', () => window.cartmesh.openPath(state.outputDirectory));

window.cartmesh.onRunLine(log);
window.addEventListener('resize', () => view.draw());

(async () => {
  state.catalog = await window.cartmesh.catalog();
  renderMethods();
  renderPresets();
  renderSamples();
  selectMethod('cutcell');
  // Test hook for `electron . --smoke=<sample>`.  The output directory normally comes
  // from a native dialog, which a headless run cannot answer, so the smoke path sets
  // it here and then goes through the same handlers a click would.
  window.__smoke = { state, selectMethod, chooseGeometry, generate, setOutput };
})();

function setOutput(directory) {
  state.outputDirectory = directory;
  $('outputPath').textContent = directory;
  updateReady();
}
