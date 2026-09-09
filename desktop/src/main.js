'use strict';

const { app, BrowserWindow, dialog, ipcMain, shell } = require('electron');
const fs = require('node:fs/promises');
const path = require('node:path');

const { METHODS, PRESETS, GEOMETRY_FORMATS } = require('./core/capabilities');
const { SAMPLES, sampleById } = require('./core/samples');
const geometry = require('./core/geometry');
const { candidates, estimateSeconds } = require('./core/automatic');
const { validateJob, buildInvocation } = require('./core/job');
const { normalizeResult, parseKeyValues } = require('./core/report');
const { parseCm2d, levelHistogram, embeddedBounds,
        assignSizeBands } = require('./core/cm2d');

let mainWindow;
let sessionDirectory;
let currentResult;
let operation = null;
const timingHistory = new Map();
async function exclusive(work) {
  if (operation) throw new Error('已有操作正在进行，请等待或取消。');
  operation = new AbortController();
  try { return await work(); } finally { operation = null; }
}
async function exportPackage(destination) {
  if (!currentResult) throw new Error('请先成功生成网格。');
  const temporary = path.join(sessionDirectory, 'export.zip');
  await fs.rm(temporary, { force: true });
  await run('/usr/bin/ditto', ['-c', '-k', '--norsrc', '--noextattr', '--noqtn', '--keepParent', currentResult.outputDirectory, temporary], () => {});
  await fs.copyFile(temporary, destination);
  await fs.rm(temporary, { force: true });
  return destination;
}

const resourceRoot = () =>
  app.isPackaged ? process.resourcesPath : path.join(__dirname, '..', 'runtime');
const resourcePath = (...parts) => path.join(resourceRoot(), ...parts);
const executable = name =>
  resourcePath('bin', process.platform === 'win32' ? `${name}.exe` : name);

const { run: runProcess } = require('./core/process');
const run = (command, args, onLine) => runProcess(command, args, onLine, operation?.signal, operation?.automatic ? 180000 : 0);

const readJson = async file => {
  try { return JSON.parse(await fs.readFile(file, 'utf8')); } catch { return null; }
};

const safeBaseName = filePath =>
  path.basename(filePath, path.extname(filePath)).replace(/[^\w-]+/g, '_') || 'mesh';

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
  const converted = geometry.convertToLoops(geometryPath, text, { chordToleranceFraction: chordError, sourceUnits });
  if (converted.issues.length) throw new Error(converted.issues.join('\n'));
  await fs.writeFile(xyPath, geometry.loopsToXyText(converted.loops,
    `converted from ${path.basename(geometryPath)} by CartMesh2D`));
  const loopSizes = converted.loops.map(loop => loop.length);
  log(`已读入 ${converted.loops.length} 个闭合环（顶点 ${loopSizes.join(' / ')}）`);
  log(`源单位 ${converted.sourceUnits}，输出单位 m（换算 ${converted.scale}）；保留物理尺寸。`);
  converted.warnings.forEach(warning => log(`注意：${warning}`));
  return { kind, converter: 'in-process', warnings: converted.warnings, loops: converted.loops,
    sourceUnits: converted.sourceUnits, outputUnits: 'm', scale: converted.scale };
}

// Which report files a run produces depends on the method, so collect them by name
// rather than guessing from the prefix.
async function collectReports(method, prefix) {
  if (method === 'hybrid') {
    return {
      resolution: await readJson(`${prefix}.resolution.json`),
      sizeField: await readJson(`${prefix}.size-field.json`),
      hybrid: await readJson(`${prefix}.hybrid.json`),
      contract: await readJson(`${prefix}.hybrid.quality-contract.json`),
      solverQuality: await readJson(`${prefix}.hybrid.solver-quality.json`)
    };
  }
  return {
    resolution: await readJson(`${prefix}.resolution.json`),
    contract: await readJson(`${prefix}.quality-contract.json`),
    sizeField: await readJson(`${prefix}.size-field.json`),
    sizing: await readJson(`${prefix}.sizing.json`)
  };
}

// Body bbox centre and span, the frame every sizing number is expressed in.
function bodyFrame(loops) {
  return geometry.boundsOfLoops(loops);
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
  sessionDirectory = await fs.mkdtemp(path.join(app.getPath('temp'), 'cartmesh2d-session-'));
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

  ipcMain.handle('cancel', () => { operation?.abort(); });
  ipcMain.handle('export-result', () => exclusive(async () => {
    if (!currentResult) throw new Error('请先成功生成网格。');
    const result = await dialog.showSaveDialog(mainWindow, {
      title: '保存网格结果包', defaultPath: `${safeBaseName(currentResult.job.geometryPath)}.zip`,
      filters: [{ name: '网格结果包', extensions: ['zip'] }]
    });
    return result.canceled ? null : exportPackage(result.filePath);
  }));

  ipcMain.handle('open-path', (_event, target) => shell.openPath(target));

  // Read a geometry without meshing it, so the outline can be drawn the moment a
  // file is chosen.  DXF needs the converter, so it writes into a scratch directory.
  ipcMain.handle('preview-geometry', async (_event, { geometryPath, chordError, sourceUnits }) => {
    const scratch = await fs.mkdtemp(path.join(sessionDirectory, 'preview-'));
    try {
    const xyPath = path.join(scratch, `${safeBaseName(geometryPath)}.xy`);
    const info = await prepareGeometry(geometryPath, { chordError, sourceUnits }, xyPath,
      path.join(scratch, 'dxf.json'), () => {});
    const loops = info.loops
      || geometry.convertToLoops(xyPath, await fs.readFile(xyPath, 'utf8')).loops;
    return { loops, kind: info.kind, warnings: info.warnings, frame: bodyFrame(loops) };
    } finally { await fs.rm(scratch, { recursive: true, force: true }); }
  });

  // Resolve the size field and stop.  This is the only way to learn the curvature and
  // proximity depths a geometry asks for, because those depend on the wall polyline
  // and cannot be predicted from the flags.
  ipcMain.handle('probe-sizing', (_event, request) => exclusive(async () => {
    const { job } = validateJob(request);
    if (job.method !== 'cutcell' && !job.relativeSizing) throw new Error('请使用相对尺寸进行预检。');
    const scratch = await fs.mkdtemp(path.join(sessionDirectory, 'probe-'));
    try {
    const prefix = path.join(scratch, safeBaseName(job.geometryPath));
    const xyPath = `${prefix}.xy`;
    const prepared = await prepareGeometry(job.geometryPath, job, xyPath, `${prefix}.dxf.json`, () => {});
    const loops = prepared.loops || geometry.convertToLoops(xyPath, await fs.readFile(xyPath, 'utf8')).loops;
    const invocation = buildInvocation(job, { xyPath, prefix, casePath: '-', frame: bodyFrame(loops) }, { dryRun: true });
    try {
      const { stdout } = await run(executable(invocation.executable), invocation.args, () => {});
      return { ok: true, values: parseKeyValues(stdout),
               field: await readJson(`${prefix}.size-field.json`) };
    } catch (error) {
      // A refused request still printed every resolved depth before the diagnosis.
      return {
        ok: false,
        field: await readJson(`${prefix}.size-field.json`),
        values: parseKeyValues(error.stdout || ''),
        issues: (error.stderr || error.message).split(/\r?\n/)
          .filter(line => line.startsWith('size_field_issue='))
          .map(line => line.replace('size_field_issue=', '')),
        message: error.message
      };
    }
      } finally { await fs.rm(scratch, { recursive: true, force: true }); }
  }));

  async function generateOnce(request) {
    const { job, method } = validateJob(request);
    if (currentResult) await fs.rm(currentResult.outputDirectory, { recursive: true, force: true });
    currentResult = null;
    // Only smoke runs may override the temporary root. Every run owns a new folder.
    const root = process.argv.some(arg => arg.startsWith('--smoke=')) && job.outputDirectory
      ? job.outputDirectory : sessionDirectory;
    await fs.mkdir(root, { recursive: true });
    job.outputDirectory = await fs.mkdtemp(path.join(root, 'mesh-'));
    const prefix = path.join(job.outputDirectory, safeBaseName(job.geometryPath));
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

    if (operation.signal.aborted) throw new Error('操作已取消');
    const reports = await collectReports(job.method, prefix);
    const mesh = await firstReadable(invocation.cm2dCandidates);
    if (!mesh) throw failure || new Error('生成结束但没有找到可预览的 CM2D 网格文件。');
    if (failure) log(`网格已写出，但后续步骤失败：${failure.message.split('\n')[0]}`);
    else log('生成完成。');
    // Parsed here rather than in the renderer: contextIsolation means the renderer
    // cannot require() the reader, and duplicating a format parser is how the two
    // copies drift apart.
    // Final solver partitions use size bands; their keys need not encode tree levels.
    const parsed = assignSizeBands(parseCm2d(mesh.text));
    const payload = {
      job,
      prefix,
      outputDirectory: job.outputDirectory,
      cm2dPath: mesh.path,
      mesh: parsed,
      levelBasis: 'size',
      levelHistogram: levelHistogram(parsed),
      wallBounds: embeddedBounds(parsed),
      incomplete: failure ? failure.message.split('\n')[0] : null,
      result: normalizeResult({ method: job.method, stdout, reports, paths, mesh: parsed, incomplete: Boolean(failure) })
    };
    if (!failure) {
      await fs.writeFile(path.join(job.outputDirectory, 'result.json'), JSON.stringify({
        parameters: { ...job, geometryPath: path.basename(job.geometryPath), outputDirectory: undefined },
        geometry: { sourceUnits: prepared.sourceUnits || 'DXF report', outputUnits: 'm', scale: prepared.scale,
          frame: paths.frame },
        result: { ...payload.result, openFoam: { ...payload.result.openFoam, path: path.basename(paths.casePath) } }
      }, null, 2));
      currentResult = payload;
    }
    return payload;
  }

  ipcMain.handle('generate', (_event, request) => exclusive(async () => {
    operation.automatic = Boolean(request.automatic);
    const sample = SAMPLES.find(item => resourcePath('samples', item.file) === request.geometryPath);
    // Get scale using the same import settings, without trusting renderer geometry.
    const scratch = await fs.mkdtemp(path.join(sessionDirectory, 'auto-'));
    let frame;
    try {
      const xy = path.join(scratch, 'input.xy');
      const imported = await prepareGeometry(request.geometryPath, request, xy, path.join(scratch, 'dxf.json'), () => {});
      frame = bodyFrame(imported.loops || geometry.convertToLoops(xy, await fs.readFile(xy, 'utf8')).loops);
    } finally { await fs.rm(scratch, { recursive: true, force: true }); }
    const choices = candidates(request, sample, frame);
    const attempts = [];
    let lastError;
    for (let index = 0; index < choices.length; index++) {
      if (operation.signal.aborted) throw new Error('操作已取消');
      const choice = choices[index];
      const started = Date.now();
      const historyKey = JSON.stringify({ ...choice, outputDirectory: undefined });
      const estimate = estimateSeconds(choice, timingHistory.get(historyKey));
      mainWindow.webContents.send('run-progress', { attempt: index + 1, maximum: choices.length,
        estimatedSeconds: estimate.seconds, estimateSource: estimate.source });
      log(`${request.automatic ? '自动选参' : '手动生成'}：第 ${index + 1}/${choices.length} 次，` +
        (choice.sizingMode === 'relative'
          ? `壁面 h/Lref=${choice.wallRelativeSize}，背景 h/Lref=${choice.backgroundRelativeSize}` +
            (choice.method === 'hybrid' ? `，首层 h/Lref=${choice.firstLayerRelativeSize}` : '')
          : choice.method === 'hybrid' ? `层级 ${choice.maxLevel}，首层 ${choice.firstThickness}` :
         `壁面体长/${choice.wallCellsPerSpan}，远场 ${choice.farFieldSpans}，α ${choice.smallAlpha}`));
      try {
        const payload = await generateOnce(choice);
        if (payload.incomplete) throw new Error(payload.incomplete);
        if (request.automatic && choice.method === 'hybrid' && payload.result.actualMethod !== 'hybrid')
          throw new Error('此参数只生成了纯 Cut-cell 回退网格，继续寻找贴体边界层参数。');
        timingHistory.set(historyKey, (Date.now() - started) / 1000);
        attempts.push({ parameters: choice, seconds: (Date.now() - started) / 1000, success: true });
        payload.automatic = Boolean(request.automatic);
        payload.densityReduced = Boolean(request.automatic && (choice.method === 'cutcell'
          ? (choice.wallCellsPerSpan < choices[0].wallCellsPerSpan ||
             choice.cellsPerLevel < choices[0].cellsPerLevel || choice.farLevel < choices[0].farLevel)
          : choice.maxLevel < choices[0].maxLevel));
        payload.attempts = attempts;
        await fs.writeFile(path.join(payload.outputDirectory, 'selection.json'), JSON.stringify({
          automatic: payload.automatic, attempts
        }, null, 2));
        return payload;
      } catch (error) {
        lastError = error;
        attempts.push({ parameters: choice, seconds: (Date.now() - started) / 1000,
          success: false, reason: error.message.split('\n')[0] });
        if (operation.signal.aborted) throw new Error('操作已取消');
        log(`本组参数未通过：${error.message.split('\n')[0]}`);
        if (/nested wall loops|invalid original wall region|invalid original wall/.test(error.message)) break;
      }
    }
    currentResult = null;
    throw new Error(`${request.automatic ? `自动尝试 ${attempts.length} 组参数后仍未通过；未降低质量标准。` : ''} ${lastError.message}`);
  }));

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
  const outputDirectory = argument('out') || '';
  const method = argument('method') || 'cutcell';
  const shot = argument('shot');
  // Optional physical sizing for repeatable large-mesh smoke runs. These
  // populate the real form after loading the geometry, through its events.
  const sizingInputs = Object.fromEntries([
    ['wallRelativeSize', 'wall-relative-size'],
    ['backgroundRelativeSize', 'background-relative-size'],
    ['referenceLength', 'reference-length'],
    ['relativePadding', 'padding-relative-size'],
    ['relativeBandCells', 'band-cells'],
    ['firstLayerRelativeSize', 'first-layer-relative-size'],
    ['nLayers', 'layers'],
    ['smallAlpha', 'small-alpha'],
  ].map(([id, flag]) => [id, argument(flag)]).filter(([, value]) => value !== null));
  if (outputDirectory) await fs.mkdir(outputDirectory, { recursive: true });

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
    document.getElementById('controlMode').value = ${JSON.stringify(argument('control') || 'auto')};
    document.getElementById('controlMode').dispatchEvent(new Event('change'));
    document.getElementById('density').value = ${JSON.stringify(argument('density') || 'normal')};
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
    const sizingInputs = ${JSON.stringify(sizingInputs)};
    if ('referenceLength' in sizingInputs) {
      const referenceMode = document.getElementById('referenceMode');
      referenceMode.value = 'explicit';
      referenceMode.dispatchEvent(new Event('change'));
    }
    for (const [id, value] of Object.entries(sizingInputs)) {
      const input = document.getElementById(id);
      input.value = value;
      if (!input.checkValidity() || !Number.isFinite(Number(input.value)))
        throw new Error('Invalid smoke sizing input: ' + id);
      input.dispatchEvent(new Event('input'));
      input.dispatchEvent(new Event('change'));
    }
    if (${JSON.stringify(argument('allow-unsafe') === 'true')}) {
      const input = document.getElementById('relativeAllowUnsafe');
      input.checked = true;
      input.dispatchEvent(new Event('change'));
    }
    const regionInput = document.querySelector('#regionList input');
    if (regionInput) {
      regionInput.focus();
      if (!regionInput.isConnected || document.activeElement !== regionInput) throw new Error('Region input lost focus');
    }
    const pending = smoke.generate();
    if (!document.getElementById('sample').disabled || document.getElementById('cancel').hidden)
      throw new Error('Parameters are not locked during generation');
    await pending;
    if (${JSON.stringify(Boolean(argument('repeat')))}) await smoke.generate();
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
    if (argument('export')) report.exported = await exportPackage(argument('export'));
    mainWindow.setSize(1120, 720);
    await new Promise(resolve => setTimeout(resolve, 200));
    report.layout = await mainWindow.webContents.executeJavaScript(`(() => {
      const panel = document.querySelector('.panel');
      panel.scrollTop = panel.scrollHeight;
      const button = document.getElementById('generate').getBoundingClientRect();
      return { bottomReachable: button.bottom <= innerHeight && button.top >= 0,
        pageHeight: document.documentElement.scrollHeight, windowHeight: innerHeight,
        previewCells: window.__smoke.state.mesh?.cells.length,
        exportedCells: window.__smoke.state.result?.openFoam.cells };
    })()`);
    console.log(JSON.stringify(report, null, 2));
    if (!report.layout.bottomReachable) throw new Error('Sidebar bottom is inaccessible');
    if (shot) {
      await new Promise(resolve => setTimeout(resolve, 400));
      await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
      console.log(`screenshot=${shot}`);
    }
    await fs.rm(sessionDirectory, { recursive: true, force: true });
    app.exit(/失败/.test(report.status) ? 1 : 0);
  }).catch(error => {
    console.error(error);
    app.exit(1);
  });
}

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

app.on('will-quit', () => {
  operation?.abort();
  if (sessionDirectory) require('node:fs').rmSync(sessionDirectory, { recursive: true, force: true });
});
