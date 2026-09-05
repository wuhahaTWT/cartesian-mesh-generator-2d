'use strict';

const { app, BrowserWindow, dialog, ipcMain, shell } = require('electron');
const { spawn } = require('node:child_process');
const fs = require('node:fs/promises');
const path = require('node:path');

const { METHODS, PRESETS, GEOMETRY_FORMATS } = require('./core/capabilities');
const { SAMPLES, sampleById } = require('./core/samples');
const geometry = require('./core/geometry');
const { validateJob, buildInvocation } = require('./core/job');
const { normalizeResult, parseKeyValues } = require('./core/report');
const { parseCm2d, levelHistogram, embeddedBounds,
        assignKeyLevels, assignSizeBands } = require('./core/cm2d');

let mainWindow;

const resourceRoot = () =>
  app.isPackaged ? process.resourcesPath : path.join(__dirname, '..', 'runtime');
const resourcePath = (...parts) => path.join(resourceRoot(), ...parts);
const executable = name =>
  resourcePath('bin', process.platform === 'win32' ? `${name}.exe` : name);

function run(command, args, onLine) {
  return new Promise((resolve, reject) => {
    const child = spawn(command, args, { windowsHide: true });
    let stdout = '';
    let stderr = '';
    const consume = (chunk, isError) => {
      const text = chunk.toString();
      if (isError) stderr += text; else stdout += text;
      text.split(/\r?\n/).filter(Boolean).forEach(onLine);
    };
    child.stdout.on('data', chunk => consume(chunk, false));
    child.stderr.on('data', chunk => consume(chunk, true));
    child.on('error', reject);
    child.on('close', code => {
      // Both CLIs are fail-closed: a non-zero exit means no mesh was committed, and
      // the reason is on stderr.  Surfacing stdout as the fallback keeps the size
      // field's refusal readable even when it printed its diagnosis first.
      if (code === 0) resolve({ stdout, stderr, code });
      else reject(Object.assign(new Error(stderr.trim() || stdout.trim() || `退出码 ${code}`),
                                { stdout, stderr, code }));
    });
  });
}

const readJson = async file => {
  try { return JSON.parse(await fs.readFile(file, 'utf8')); } catch { return null; }
};

const safeBaseName = filePath =>
  path.basename(filePath, path.extname(filePath)).replace(/[^\w-]+/g, '_') || 'mesh';

async function uniquePrefix(directory, base) {
  const first = path.join(directory, base);
  try {
    await fs.access(`${first}.cm2d`);
  } catch {
    return first;
  }
  return path.join(directory, `${base}-${new Date().toISOString().replace(/[:.]/g, '-')}`);
}

// Every input becomes a native .xy before the mesher sees it.  DXF goes through the
// C++ converter because that is where unit handling and entity diagnostics live;
// everything else is converted in process.
async function prepareGeometry(geometryPath, { chordError, sourceUnits }, xyPath, reportPath, log) {
  const kind = geometry.classify(geometryPath);
  if (kind === null) throw new Error(`不支持的文件类型：${path.extname(geometryPath) || '(无扩展名)'}`);

  if (kind === 'dxf') {
    log('正在读取 DXF、换算单位并离散曲线…');
    const args = [geometryPath, xyPath, String(chordError), reportPath];
    if (sourceUnits && sourceUnits !== 'auto') args.push('1e-10', sourceUnits);
    await run(executable('cartmesh2d_dxf_cli'), args, log);
    return { kind, converter: 'cartmesh2d_dxf_cli', warnings: [] };
  }

  const text = await fs.readFile(geometryPath, 'utf8');
  // The SVG chord tolerance is a fraction of the drawing's own extent because SVG
  // carries no units; .xy and the coordinate formats are taken as given.
  const converted = geometry.convertToLoops(geometryPath, text, { chordToleranceFraction: chordError });
  if (converted.issues.length) throw new Error(converted.issues.join('\n'));
  await fs.writeFile(xyPath, geometry.loopsToXyText(converted.loops,
    `converted from ${path.basename(geometryPath)} by CartMesh2D`));
  const loopSizes = converted.loops.map(loop => loop.length);
  log(`已读入 ${converted.loops.length} 个闭合环（顶点 ${loopSizes.join(' / ')}）`);
  if (converted.normalized) {
    log(`该格式无物理单位，已把最大跨度归一到 1 m（缩放 ${converted.scale.toPrecision(4)}）`);
  }
  converted.warnings.forEach(warning => log(`注意：${warning}`));
  return { kind, converter: 'in-process', warnings: converted.warnings, loops: converted.loops };
}

// Which report files a run produces depends on the method, so collect them by name
// rather than guessing from the prefix.
async function collectReports(method, prefix) {
  if (method === 'hybrid') {
    return {
      hybrid: await readJson(`${prefix}.hybrid.json`),
      contract: await readJson(`${prefix}.hybrid.quality-contract.json`),
      solverQuality: await readJson(`${prefix}.hybrid.solver-quality.json`)
    };
  }
  return {
    contract: await readJson(`${prefix}.quality-contract.json`),
    sizeField: await readJson(`${prefix}.size-field.json`),
    sizing: await readJson(`${prefix}.sizing.json`)
  };
}

// Body bbox centre and span, the frame every sizing number is expressed in.
function bodyFrame(loops) {
  const points = loops.flat();
  const xs = points.map(point => point[0]);
  const ys = points.map(point => point[1]);
  const minX = Math.min(...xs), maxX = Math.max(...xs);
  const minY = Math.min(...ys), maxY = Math.max(...ys);
  return {
    centreX: (minX + maxX) / 2,
    centreY: (minY + maxY) / 2,
    bodySpan: Math.max(maxX - minX, maxY - minY)
  };
}

async function firstReadable(candidates) {
  for (const candidate of candidates) {
    try { return { path: candidate, text: await fs.readFile(candidate, 'utf8') }; } catch { /* next */ }
  }
  return null;
}

async function createWindow() {
  mainWindow = new BrowserWindow({    width: 1440,
    height: 900,
    minWidth: 1120,
    minHeight: 720,
    backgroundColor: '#10161c',
    titleBarStyle: 'hiddenInset',
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false
    }
  });
  // Renderer errors are otherwise invisible from a headless smoke run.
  mainWindow.webContents.on('console-message', (_event, level, message) => {
    if (level >= 2) console.error(`renderer: ${message}`);
  });
  await mainWindow.loadFile(path.join(__dirname, 'renderer', 'index.html'));
}
app.whenReady().then(async () => {
  const log = line => mainWindow?.webContents.send('run-line', line);

  ipcMain.handle('catalog', () => ({
    methods: METHODS,
    presets: PRESETS,
    formats: GEOMETRY_FORMATS,
    samples: SAMPLES.map(sample => ({ ...sample, path: resourcePath('samples', sample.file) }))
  }));

  ipcMain.handle('pick-geometry', async () => {
    const result = await dialog.showOpenDialog(mainWindow, {
      title: '选择二维边界几何',
      filters: [
        { name: '所有支持的格式', extensions: ['xy', 'dxf', 'svg', 'csv', 'txt', 'dat'] },
        { name: '原生折线', extensions: ['xy'] },
        { name: 'AutoCAD DXF', extensions: ['dxf'] },
        { name: 'SVG', extensions: ['svg'] },
        { name: '坐标表', extensions: ['csv', 'txt', 'dat'] }
      ],
      properties: ['openFile']
    });
    return result.canceled ? null : result.filePaths[0];
  });

  ipcMain.handle('pick-output', async () => {
    const result = await dialog.showOpenDialog(mainWindow, {
      title: '选择输出目录',
      properties: ['openDirectory', 'createDirectory']
    });
    return result.canceled ? null : result.filePaths[0];
  });

  ipcMain.handle('open-path', (_event, target) => shell.openPath(target));

  // Read a geometry without meshing it, so the outline can be drawn the moment a
  // file is chosen.  DXF needs the converter, so it writes into a scratch directory.
  ipcMain.handle('preview-geometry', async (_event, { geometryPath, chordError, sourceUnits }) => {
    const scratch = path.join(app.getPath('temp'), 'cartmesh2d-preview');
    await fs.mkdir(scratch, { recursive: true });
    const xyPath = path.join(scratch, `${safeBaseName(geometryPath)}.xy`);
    const info = await prepareGeometry(geometryPath, { chordError, sourceUnits }, xyPath,
      path.join(scratch, 'dxf.json'), () => {});
    const loops = info.loops
      || geometry.convertToLoops(xyPath, await fs.readFile(xyPath, 'utf8')).loops;
    return { loops, kind: info.kind, warnings: info.warnings, frame: bodyFrame(loops) };
  });

  // Resolve the size field and stop.  This is the only way to learn the curvature and
  // proximity depths a geometry asks for, because those depend on the wall polyline
  // and cannot be predicted from the flags.
  ipcMain.handle('probe-sizing', async (_event, request) => {
    const { job } = validateJob(request);
    if (job.method !== 'cutcell') throw new Error('贴体边界层路径没有 size field 预检。');
    const scratch = path.join(app.getPath('temp'), 'cartmesh2d-probe');
    await fs.mkdir(scratch, { recursive: true });
    const prefix = path.join(scratch, safeBaseName(job.geometryPath));
    const xyPath = `${prefix}.xy`;
    await prepareGeometry(job.geometryPath, job, xyPath, `${prefix}.dxf.json`, () => {});
    const invocation = buildInvocation(job, { xyPath, prefix, casePath: '-' }, { dryRun: true });
    try {
      const { stdout } = await run(executable(invocation.executable), invocation.args, () => {});
      return { ok: true, values: parseKeyValues(stdout),
               field: await readJson(`${prefix}.size-field.json`) };
    } catch (error) {
      // A refused request still printed every resolved depth before the diagnosis.
      return {
        ok: false,
        values: parseKeyValues(error.stdout || ''),
        issues: (error.stderr || error.message).split(/\r?\n/)
          .filter(line => line.startsWith('size_field_issue='))
          .map(line => line.replace('size_field_issue=', '')),
        message: error.message
      };
    }
  });

  ipcMain.handle('generate', async (_event, request) => {
    const { job, method } = validateJob(request);
    await fs.mkdir(job.outputDirectory, { recursive: true });
    const prefix = await uniquePrefix(job.outputDirectory, safeBaseName(job.geometryPath));
    const paths = { prefix, xyPath: `${prefix}.xy`, casePath: `${prefix}-openfoam` };

    const prepared = await prepareGeometry(job.geometryPath, job, paths.xyPath,
      `${prefix}.dxf.json`, log);
    // A hand-placed region is stated in body spans about the body centre, so the frame
    // has to come from the same loops the mesher is about to read.
    const loops = prepared.loops
      || geometry.convertToLoops(paths.xyPath, await fs.readFile(paths.xyPath, 'utf8')).loops;
    paths.frame = bodyFrame(loops);
    const invocation = buildInvocation(job, paths);
    log(`正在生成${method.label}网格…`);

    // The generator writes the mesh, the size field and every quality report before it
    // attempts the OpenFOAM export, so a late failure still leaves a usable mesh on
    // disk.  Salvaging it beats reporting only the exception: the user gets the
    // preview, the gates, and a precise statement of what did not get written.
    let stdout = '';
    let failure = null;
    try {
      stdout = (await run(executable(invocation.executable), invocation.args, log)).stdout;
    } catch (error) {
      stdout = error.stdout || '';
      failure = error;
    }

    const reports = await collectReports(job.method, prefix);
    const mesh = await firstReadable(invocation.cm2dCandidates);
    if (!mesh) throw failure || new Error('生成结束但没有找到可预览的 CM2D 网格文件。');
    if (failure) log(`网格已写出，但后续步骤失败：${failure.message.split('\n')[0]}`);
    else log('生成完成。');
    // Parsed here rather than in the renderer: contextIsolation means the renderer
    // cannot require() the reader, and duplicating a format parser is how the two
    // copies drift apart.
    // Only the pure path's sourceKey carries a Quadtree level; the hybrid writes a
    // running index there, so it is banded by cell size instead.
    const parsed = job.method === 'hybrid'
      ? assignSizeBands(parseCm2d(mesh.text))
      : assignKeyLevels(parseCm2d(mesh.text));
    return {
      job,
      prefix,
      outputDirectory: job.outputDirectory,
      cm2dPath: mesh.path,
      mesh: parsed,
      levelBasis: job.method === 'hybrid' ? 'size' : 'level',
      levelHistogram: levelHistogram(parsed),
      wallBounds: embeddedBounds(parsed),
      incomplete: failure ? failure.message.split('\n')[0] : null,
      result: normalizeResult({ method: job.method, stdout, reports, paths, mesh: parsed })
    };
  });

  await createWindow();
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
  if (process.argv.some(item => item.startsWith('--smoke='))) await runSmoke();
});

// `--smoke=<sample-id> [--out=<dir>] [--method=<id>] [--shot=<png>]` drives the real
// renderer through a full run and writes a screenshot.  Clicking the actual controls
// is the only check that covers the renderer, the IPC surface and the CLI together.
async function runSmoke() {
  const argument = name => {
    const found = process.argv.find(item => item.startsWith(`--${name}=`));
    return found ? found.slice(name.length + 3) : null;
  };
  const sampleId = argument('smoke');
  const outputDirectory = argument('out') || path.join(app.getPath('temp'), 'cartmesh2d-smoke');
  const method = argument('method') || 'cutcell';
  const shot = argument('shot');
  await fs.mkdir(outputDirectory, { recursive: true });

  // The renderer's init awaits the catalog over IPC, so the hook appears a moment
  // after the page finishes loading.
  for (let attempt = 0; attempt < 100; attempt++) {
    if (await mainWindow.webContents.executeJavaScript('Boolean(window.__smoke)')) break;
    await new Promise(resolve => setTimeout(resolve, 100));
  }

  await mainWindow.webContents.executeJavaScript(`(async () => {
    const smoke = window.__smoke;
    smoke.setOutput(${JSON.stringify(outputDirectory)});
    smoke.selectMethod(${JSON.stringify(method)});
    if (${JSON.stringify(Boolean(argument('regions')))}) {
      smoke.addRegion();
      smoke.addRegion();
      smoke.state.regions[1].xmin = 6; smoke.state.regions[1].xmax = 14;
      smoke.state.regions[1].ymin = -1.6; smoke.state.regions[1].ymax = 1.6;
      smoke.state.regions[1].levelsBelowWall = 6;
      smoke.renderRegions();
    }
    const sample = smoke.state.catalog.samples.find(item => item.id === ${JSON.stringify(sampleId)});
    if (!sample) throw new Error('unknown sample ' + ${JSON.stringify(sampleId)});
    document.getElementById('sample').value = sample.id;
    await smoke.chooseGeometry(sample.path, sample.label, sample);
    await smoke.generate();
    const mode = ${JSON.stringify(argument('mode') || 'level')};
    if (mode !== 'level') {
      const select = document.getElementById('displayMode');
      select.value = mode;
      select.dispatchEvent(new Event('change'));
    }
    return {
      status: document.getElementById('statusTitle').textContent,
      detail: document.getElementById('statusText').textContent,
      counters: document.getElementById('counters').innerText,
      gates: document.getElementById('gates').innerText,
      histogram: document.getElementById('histogram').innerText,
      log: document.getElementById('log').textContent
    };
  })()`).then(async report => {
    console.log(JSON.stringify(report, null, 2));
    if (shot) {
      await new Promise(resolve => setTimeout(resolve, 400));
      await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
      console.log(`screenshot=${shot}`);
    }
    app.exit(/失败/.test(report.status) ? 1 : 0);
  }).catch(error => {
    console.error(error);
    app.exit(1);
  });
}

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});
