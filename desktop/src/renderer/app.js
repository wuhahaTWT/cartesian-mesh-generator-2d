const state = { dxfPath: '', outputDirectory: '', mesh: null, output: '', smokeMode: false };
const $ = id => document.getElementById(id);
const presets = {
  quick: { maxLevel: 5, minimumLevel: 0 },
  standard: { maxLevel: 6, minimumLevel: 0 },
  dense: { maxLevel: 8, minimumLevel: 0 }
};

function updateReady() { $('generate').disabled = !(state.dxfPath && state.outputDirectory); }
function compactPath(value) { return value || '尚未选择'; }
function logLine(line) { $('log').textContent += `${line}\n`; $('log').scrollTop = $('log').scrollHeight; }
function updateRefinementCostHint() {
  const minimumLevel = Number($('minimumLevel').value);
  const globalLeaves = Number.isInteger(minimumLevel) && minimumLevel >= 0
    ? Math.pow(4, minimumLevel)
    : null;
  $('refinementCostHint').textContent = minimumLevel > 0 && Number.isFinite(globalLeaves)
    ? `当前全域最低层级会先铺约 ${globalLeaves.toLocaleString()} 个叶格，再叠加物面局部加密；每提高一级约再乘 4。`
    : '建议保持全域最低层级为 0；物面精度由最高层级控制。全域最低层级每提高一级，远场叶格约增至 4 倍。';
}

function parseCm2d(text) {
  const lines = text.trim().split(/\r?\n/);
  const vertices = [];
  const cells = [];
  const edges = [];
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
  for (const edge of mesh.edges.filter(e => e.kind !== 0)) {
    const [ax, ay] = point(mesh.vertices[edge.a]);
    const [bx, by] = point(mesh.vertices[edge.b]);
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
    ctx.strokeStyle = edge.kind === 1 ? '#ef5050' : '#203840';
    ctx.lineWidth = edge.kind === 1 ? 2.1 : 1.4;
    ctx.stroke();
  }
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
document.querySelectorAll('.preset').forEach(button => button.addEventListener('click', () => {
  document.querySelectorAll('.preset').forEach(item => item.classList.remove('active'));
  button.classList.add('active');
  const preset = presets[button.dataset.preset];
  $('maxLevel').value = preset.maxLevel;
  $('minimumLevel').value = preset.minimumLevel;
  updateRefinementCostHint();
}));
$('minimumLevel').addEventListener('input', updateRefinementCostHint);

window.cartmesh.onGenerationLine(logLine);
$('generate').addEventListener('click', async () => {
  $('generate').disabled = true;
  $('generate').textContent = '正在生成…';
  $('log').textContent = '';
  $('statusTitle').textContent = '生成中';
  $('statusText').textContent = '质量门通过后才会输出 OpenFOAM case';
  try {
    const result = await window.cartmesh.generate({
      dxfPath: state.dxfPath,
      outputDirectory: state.outputDirectory,
      sourceUnits: $('sourceUnits').value,
      chordError: Number($('chordError').value),
      maxLevel: Number($('maxLevel').value),
      minimumLevel: Number($('minimumLevel').value),
      paddingFraction: Number($('padding').value),
      smallAlpha: Number($('smallAlpha').value)
    });
    state.mesh = parseCm2d(result.cm2d);
    state.output = result.outputDirectory;
    drawMesh(state.mesh);
    $('emptyState').hidden = true;
    $('legend').hidden = false;
    $('openOutput').hidden = false;
    $('cellCount').textContent = Number(result.summary.stabilized_cells).toLocaleString();
    $('foamCells').textContent = Number(result.summary.openfoam_cells).toLocaleString();
    $('vertexCount').textContent = Number(result.summary.vertices).toLocaleString();
    $('quality').textContent = result.summary.solver_quality || 'PASS';
    $('quality').style.color = '#087d65';
    $('statusTitle').textContent = '生成完成';
    $('statusText').textContent = result.cm2dPath;
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
  const smoke = await window.cartmesh.smokeConfig();
  if (!smoke.enabled) return;
  state.smokeMode = true;
  state.dxfPath = await window.cartmesh.examplePath();
  state.outputDirectory = smoke.outputDirectory;
  $('dxfPath').textContent = state.dxfPath;
  $('outputPath').textContent = state.outputDirectory;
  updateReady();
  $('generate').click();
})();

updateRefinementCostHint();
