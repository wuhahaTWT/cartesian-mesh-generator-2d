const state = {
  dxfPath: '', outputDirectory: '', mesh: null, output: '', smokeMode: false,
  capabilities: null, methodId: 'cutcell'
};
const $ = id => document.getElementById(id);

function updateReady() { $('generate').disabled = !(state.dxfPath && state.outputDirectory && state.methodId); }
function compactPath(value) { return value || '尚未选择'; }
function logLine(line) { $('log').textContent += `${line}\n`; $('log').scrollTop = $('log').scrollHeight; }
function activeMethod() { return state.capabilities?.methods.find(method => method.id === state.methodId); }
function numberValue(id) { return Number($(id).value); }

function parseCm2d(text) {
  const lines = text.trim().split(/\r?\n/);
  const vertices = [], cells = [], edges = [];
  let i = 1;
  while (i < lines.length) {
    const header = lines[i].trim().split(/\s+/);
    if (header[0] === 'VERTICES') {
      const count = Number(header[1]);
      for (let n = 0; n < count; n++) {
        const p = lines[++i].trim().split(/\s+/).map(Number);
        vertices[p[0]] = [p[1], p[2]];
      }
    } else if (header[0] === 'EDGES') {
      const count = Number(header[1]);
      for (let n = 0; n < count; n++) {
        const e = lines[++i].trim().split(/\s+/).map(Number);
        edges.push({ a: e[1], b: e[2], kind: e[5] });
      }
    } else if (header[0] === 'CELLS') {
      const count = Number(header[1]);
      for (let n = 0; n < count; n++) {
        const c = lines[++i].trim().split(/\s+/).map(Number);
        const vertexCount = c[4];
        cells.push({ level: c[2] & 63, vertices: c.slice(5, 5 + vertexCount) });
      }
    }
    i++;
  }
  if (!vertices.length || !cells.length) throw new Error('CM2D 文件没有可显示的单元。');
  return { vertices, cells, edges };
}

function drawMesh(mesh) {
  const canvas = $('meshCanvas');
  const ratio = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = Math.max(1, Math.floor(rect.width * ratio));
  canvas.height = Math.max(1, Math.floor(rect.height * ratio));
  const ctx = canvas.getContext('2d');
  ctx.scale(ratio, ratio);
  ctx.clearRect(0, 0, rect.width, rect.height);
  const xs = mesh.vertices.map(v => v[0]);
  const ys = mesh.vertices.map(v => v[1]);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const pad = 34;
  const scale = Math.min((rect.width - 2 * pad) / (maxX - minX), (rect.height - 2 * pad) / (maxY - minY));
  const ox = (rect.width - (maxX - minX) * scale) / 2;
  const oy = (rect.height - (maxY - minY) * scale) / 2;
  const point = v => [ox + (v[0] - minX) * scale, rect.height - oy - (v[1] - minY) * scale];

  for (const cell of mesh.cells) {
    ctx.beginPath();
    cell.vertices.forEach((id, index) => {
      const [x, y] = point(mesh.vertices[id]);
      if (index === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    });
    ctx.closePath();
    const shade = Math.max(0, Math.min(5, cell.level - 3));
    ctx.fillStyle = `hsl(${181 - shade * 6} 48% ${92 - shade * 4}%)`;
    ctx.fill();
    ctx.strokeStyle = '#5e909277';
    ctx.lineWidth = .65;
    ctx.stroke();
  }
  for (const edge of mesh.edges.filter(edge => edge.kind !== 0)) {
    const [ax, ay] = point(mesh.vertices[edge.a]);
    const [bx, by] = point(mesh.vertices[edge.b]);
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
    ctx.strokeStyle = edge.kind === 1 ? '#ef5050' : '#203840';
    ctx.lineWidth = edge.kind === 1 ? 2.1 : 1.4;
    ctx.stroke();
  }
}

function renderMethods() {
  $('methodCards').replaceChildren(...state.capabilities.methods.map(method => {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = `method-card${method.id === state.methodId ? ' active' : ''}`;
    button.dataset.method = method.id;
    button.innerHTML = `<span class="method-card-title">${method.label}<i class="status-${method.status}">${method.statusLabel}</i></span><small>${method.description}</small>`;
    button.addEventListener('click', () => selectMethod(method.id));
    return button;
  }));
}

function setValue(id, value) { if (value !== undefined && $(id)) $(id).value = value; }
function applyMethodDefaults(method) {
  for (const [key, value] of Object.entries(method.defaults)) setValue(key, value);
}

function selectMethod(methodId) {
  state.methodId = methodId;
  const method = activeMethod();
  renderMethods();
  applyMethodDefaults(method);
  $('cutcellParameters').hidden = method.id !== 'cutcell';
  $('hybridParameters').hidden = method.id !== 'hybrid';
  $('localRefinementSection').hidden = !method.supports.distanceBands && !method.supports.refinementBoxes;
  $('methodBoundary').textContent = method.status === 'beta'
    ? 'Beta 边界：高层级和狭窄间隙可能失败或明确退化；当前只生成检查网格，不输出 OpenFOAM case。'
    : '稳定默认路径；质量失败会停止输出，不删除坏单元或隐藏告警。';
  document.querySelectorAll('.preset').forEach(button => {
    button.classList.toggle('active', button.dataset.preset === 'standard');
  });
  updateEstimate();
  updateReady();
}

function applyPreset(name) {
  const preset = activeMethod().presets[name];
  for (const [key, value] of Object.entries(preset)) setValue(key, value);
  document.querySelectorAll('.preset').forEach(button => {
    button.classList.toggle('active', button.dataset.preset === name);
  });
  updateEstimate();
}

function readRows(containerId, keys) {
  return [...$(containerId).querySelectorAll('.parameter-row')].map(row => {
    const value = {};
    for (const key of keys) value[key] = Number(row.querySelector(`[data-key="${key}"]`).value);
    return value;
  });
}

function collectRequest() {
  return {
    method: state.methodId,
    dxfPath: state.dxfPath,
    outputDirectory: state.outputDirectory,
    sourceUnits: $('sourceUnits').value,
    chordError: numberValue('chordError'),
    maxLevel: numberValue('maxLevel'),
    minimumLevel: numberValue('minimumLevel'),
    cellBudget: numberValue('cellBudget'),
    paddingFraction: numberValue('paddingFraction'),
    smallAlpha: numberValue('smallAlpha'),
    boundaryLevel: numberValue('boundaryLevel'),
    nLayers: numberValue('nLayers'),
    firstThickness: numberValue('firstThickness'),
    growthRatio: numberValue('growthRatio'),
    domainPadding: numberValue('domainPadding'),
    extrusionThickness: numberValue('extrusionThickness'),
    distanceBands: readRows('distanceBandRows', ['distance', 'targetLevel']),
    refinementBoxes: readRows('refinementBoxRows', ['xmin', 'ymin', 'xmax', 'ymax', 'targetLevel'])
  };
}

async function updateEstimate() {
  const estimate = await window.cartmesh.estimate({
    maxLevel: numberValue('maxLevel'), minimumLevel: numberValue('minimumLevel')
  });
  if (!Number.isFinite(estimate.globalLeafFloor)) {
    $('globalLeafEstimate').textContent = '参数无效';
    $('refinementCostHint').textContent = estimate.note;
    return;
  }
  $('globalLeafEstimate').textContent = `≥ ${estimate.globalLeafFloor.toLocaleString()} leaves`;
  const budget = numberValue('cellBudget');
  const overBudget = Number.isFinite(budget) && estimate.globalLeafFloor > budget;
  $('globalLeafEstimate').classList.toggle('over-budget', overBudget);
  $('refinementCostHint').textContent = `${estimate.note} 全域最低层级每 +1，底格约 ×4。`;
}

function addDistanceBand(values = {}) {
  const row = document.createElement('div');
  row.className = 'parameter-row distance-row';
  row.innerHTML = `<label>距离（米）<input data-key="distance" type="number" min="0.000000001" step="0.01" value="${values.distance ?? 0.1}"></label><label>目标层级<input data-key="targetLevel" type="number" min="0" max="28" value="${values.targetLevel ?? numberValue('maxLevel')}"></label><button type="button" class="remove-row" title="删除">×</button>`;
  row.querySelector('.remove-row').addEventListener('click', () => row.remove());
  $('distanceBandRows').append(row);
}

function addRefinementBox(values = {}) {
  const row = document.createElement('div');
  row.className = 'parameter-row box-row';
  row.innerHTML = `<label>xmin<input data-key="xmin" type="number" step="0.1" value="${values.xmin ?? 0}"></label><label>ymin<input data-key="ymin" type="number" step="0.1" value="${values.ymin ?? -0.5}"></label><label>xmax<input data-key="xmax" type="number" step="0.1" value="${values.xmax ?? 2}"></label><label>ymax<input data-key="ymax" type="number" step="0.1" value="${values.ymax ?? 0.5}"></label><label>层级<input data-key="targetLevel" type="number" min="0" max="28" value="${values.targetLevel ?? numberValue('maxLevel')}"></label><button type="button" class="remove-row" title="删除">×</button>`;
  row.querySelector('.remove-row').addEventListener('click', () => row.remove());
  $('refinementBoxRows').append(row);
}

$('chooseDxf').addEventListener('click', async () => {
  const value = await window.cartmesh.selectDxf();
  if (value) { state.dxfPath = value; $('dxfPath').textContent = compactPath(value); updateReady(); }
});
$('useExample').addEventListener('click', async () => {
  state.dxfPath = await window.cartmesh.examplePath();
  $('dxfPath').textContent = state.dxfPath;
  updateReady();
});
$('chooseOutput').addEventListener('click', async () => {
  const value = await window.cartmesh.selectOutput();
  if (value) { state.outputDirectory = value; $('outputPath').textContent = compactPath(value); updateReady(); }
});
document.querySelectorAll('.preset').forEach(button => {
  button.addEventListener('click', () => applyPreset(button.dataset.preset));
});
['maxLevel', 'minimumLevel', 'cellBudget'].forEach(id => $(id).addEventListener('input', updateEstimate));
$('addDistanceBand').addEventListener('click', () => addDistanceBand());
$('addRefinementBox').addEventListener('click', () => addRefinementBox());

window.cartmesh.onGenerationLine(logLine);
$('generate').addEventListener('click', async () => {
  $('generate').disabled = true;
  $('generate').textContent = '正在生成…';
  $('log').textContent = '';
  $('statusTitle').textContent = '生成中';
  $('statusText').textContent = '质量门通过后才会输出 OpenFOAM case';
  try {
    const result = await window.cartmesh.generate(collectRequest());
    state.mesh = parseCm2d(result.cm2d);
    state.output = result.outputDirectory;
    drawMesh(state.mesh);
    $('emptyState').hidden = true;
    $('legend').hidden = false;
    $('openOutput').hidden = false;
    const isFallback = result.summary.actual_method === 'cutcell-fallback';
    $('actualMethod').textContent = isFallback ? 'Cut-cell fallback' : result.method.label;
    $('actualMethod').classList.toggle('method-warning', isFallback);
    $('cellCount').textContent = Number(result.summary.stabilized_cells || state.mesh.cells.length).toLocaleString();
    $('foamCells').textContent = result.summary.openfoam_cells
      ? Number(result.summary.openfoam_cells).toLocaleString() : '未输出';
    $('vertexCount').textContent = state.mesh.vertices.length.toLocaleString();
    const qualityStatus = result.summary.quality_status || result.summary.solver_quality || 'PASS';
    const qualityPass = result.summary.quality_pass === undefined
      ? String(qualityStatus).toLowerCase() === 'pass' : result.summary.quality_pass;
    $('quality').textContent = String(qualityStatus).toUpperCase();
    $('quality').style.color = qualityPass ? '#087d65' : '#c23c3c';
    $('statusTitle').textContent = isFallback
      ? '生成完成（发生降级）'
      : qualityPass ? '生成完成' : '已生成检查网格（质量未通过）';
    $('statusText').textContent = isFallback
      ? `请求 Hybrid，实际输出 Pure Cut-cell：${result.cm2dPath}`
      : qualityPass ? result.cm2dPath : `未输出 OpenFOAM case：${result.cm2dPath}`;
    if (state.smokeMode) {
      await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      await window.cartmesh.smokeCapture();
    }
  } catch (error) {
    logLine(error.message);
    $('quality').textContent = 'FAIL';
    $('quality').style.color = '#c23c3c';
    $('statusTitle').textContent = '生成失败';
    $('statusText').textContent = error.message.split('\n')[0];
  } finally {
    $('generate').disabled = false;
    $('generate').textContent = '生成网格';
  }
});
$('openOutput').addEventListener('click', () => window.cartmesh.openPath(state.output));
window.addEventListener('resize', () => { if (state.mesh) drawMesh(state.mesh); });

(async () => {
  state.capabilities = await window.cartmesh.capabilities();
  const smoke = await window.cartmesh.smokeConfig();
  selectMethod(smoke.enabled ? smoke.method : 'cutcell');
  if (!smoke.enabled) return;
  state.smokeMode = true;
  state.dxfPath = await window.cartmesh.examplePath();
  state.outputDirectory = smoke.outputDirectory;
  $('dxfPath').textContent = state.dxfPath;
  $('outputPath').textContent = state.outputDirectory;
  updateReady();
  $('generate').click();
})();
