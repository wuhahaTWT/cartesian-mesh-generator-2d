'use strict';

const { app, BrowserWindow, dialog, ipcMain, shell, Menu } = require('electron');
const fs = require('node:fs/promises');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { calibrateRasterLoops, inspectRaster } = require('./core/raster-geometry');

const { METHODS, PRESETS, GEOMETRY_FORMATS } = require('./core/capabilities');
const { SAMPLES, sampleById } = require('./core/samples');
const geometry = require('./core/geometry');
const { candidates, estimateSeconds } = require('./core/automatic');
const { planBudget, BUDGET_PRESETS } = require('./core/cell-budget');
const { runBudget } = require('./core/budget-runner');
const { validateJob, buildInvocation } = require('./core/job');
const { normalizeResult, parseKeyValues } = require('./core/report');
const { exportGuide } = require('./core/export-guide');
const { zipDirectory } = require('./core/archive');
const { parseCm2d, levelHistogram, embeddedBounds,
        assignSizeBands } = require('./core/cm2d');

let mainWindow;
let sessionDirectory;
let currentResult;
let operation = null;
const rasterSources = new Map();
const rasterImports = new Map();
const timingHistory = new Map();
async function exclusive(work) {
  if (operation) throw new Error('已有操作正在进行，请等待或取消。');
  operation = new AbortController();
  try { return await work(); } finally { operation = null; }
}
async function exportPackage(destination) {
  if (!currentResult) throw new Error('请先成功生成网格。');
  const png = await mainWindow.webContents.executeJavaScript('window.__exportMeshPreview()');
  if (typeof png !== 'string' || !png.startsWith('data:image/png;base64,'))
    throw new Error('网格预览图片生成失败，未写出结果包。');
  await fs.writeFile(path.join(currentResult.outputDirectory, 'mesh-preview.png'), Buffer.from(png.split(',')[1], 'base64'));
  await fs.writeFile(path.join(currentResult.outputDirectory, 'README_CN.md'), exportGuide(currentResult));
  const temporary = path.join(sessionDirectory, 'export.zip');
  await fs.rm(temporary, { force: true });
  await zipDirectory(currentResult.outputDirectory, temporary, operation?.signal);
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
  const raster = rasterImports.get(geometryPath);
  if (raster) {
    await fs.copyFile(geometryPath, xyPath);
    const converted = geometry.convertToLoops(xyPath, await fs.readFile(xyPath, 'utf8'), { sourceUnits: 'm' });
    if (converted.issues.length) throw new Error(converted.issues.join('\n'));
    log('使用已确认的图片轮廓与实际尺寸（m）。');
    return { ...converted, kind: 'raster', converter: 'local-raster-contours', warnings: raster.warnings };
  }
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
      solverQuality: await readJson(`${prefix}.hybrid.solver-quality.json`)
    };
  }
  return {
    resolution: await readJson(`${prefix}.resolution.json`),
    solverQuality: await readJson(`${prefix}-openfoam/solver_quality.json`),
    sizeField: await readJson(`${prefix}.size-field.json`),
    sizing: await readJson(`${prefix}.sizing.json`)
  };
}

// Body bbox centre and span, the frame every sizing number is expressed in.
function bodyFrame(loops) {
  const frame = geometry.boundsOfLoops(loops);
  if (loops.length === 1) {
    const loop = loops[0]; const [ox, oy] = loop[0];
    frame.bodyArea = Math.abs(loop.reduce((sum, p, i) => {
      const q = loop[(i+1)%loop.length];
      return sum + (p[0]-ox)*(q[1]-oy) - (q[0]-ox)*(p[1]-oy);
    }, 0)) / 2;
  }
  return frame;
}

function budgetFrame(request, frame) {
  const reference = request.referenceLength ?? frame.bodySpan;
  const domainSpan = frame.bodySpan + 2 * Number(request.farFieldSpans) * reference;
  return { ...frame, fluidArea: frame.bodyArea > 0
    ? (request.fluidRegion === 'interior' ? frame.bodyArea : domainSpan*domainSpan-frame.bodyArea)
    : undefined };
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
    minWidth: 800,
    minHeight: 560,
    backgroundColor: '#10161c',
    titleBarStyle: process.platform === 'darwin' ? 'hiddenInset' : 'default',
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
  if (process.platform !== 'darwin') Menu.setApplicationMenu(null);
  sessionDirectory = await fs.mkdtemp(path.join(app.getPath('temp'), 'cartmesh2d-session-'));
  const log = line => mainWindow?.webContents.send('run-line', line);

  ipcMain.handle('catalog', () => ({
    methods: METHODS,
    presets: PRESETS,
    budgetPresets: BUDGET_PRESETS,
    formats: GEOMETRY_FORMATS,
    samples: SAMPLES.map(sample => ({ ...sample, path: resourcePath('samples', sample.file) }))
  }));

  ipcMain.handle('plan-budget', (_event, { request, frame }) => planBudget(request, budgetFrame(request, frame)));

  ipcMain.handle('pick-geometry', async () => {
    const result = await dialog.showOpenDialog(mainWindow, {
      title: '选择二维边界几何',
      filters: [
        { name: '所有支持的格式', extensions: ['xy', 'dxf', 'svg', 'csv', 'txt', 'dat', 'png', 'jpg', 'jpeg'] },
        { name: '图片轮廓', extensions: ['png', 'jpg', 'jpeg'] },
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
  ipcMain.handle('release-preview', () => {
    if (operation) throw new Error('请先等待当前操作完成。');
    if (currentResult) { currentResult.mesh = null; currentResult.levelHistogram = null; }
    return Boolean(currentResult);
  });
  ipcMain.handle('export-preview-data', async () => {
    if (!currentResult) throw new Error('没有可导出的网格。');
    return { mesh: currentResult.mesh || parseCm2d(await fs.readFile(currentResult.cm2dPath, 'utf8')),
      result: currentResult.result };
  });
  ipcMain.handle('export-result', () => exclusive(async () => {
    if (!currentResult) throw new Error('请先成功生成网格。');
    const result = await dialog.showSaveDialog(mainWindow, {
      title: '保存网格结果包', defaultPath: `${safeBaseName(currentResult.job.geometryPath)}.zip`,
      filters: [{ name: '网格结果包', extensions: ['zip'] }]
    });
    return result.canceled ? null : exportPackage(result.filePath);
  }));

  ipcMain.handle('read-raster', async (_event, sourcePath) => {
    if (operation) throw new Error('请等待当前操作完成。');
    if (typeof sourcePath !== 'string' || !/\.(png|jpe?g)$/i.test(sourcePath)) throw new Error('请选择 PNG 或 JPG 图片。');
    const stat = await fs.stat(sourcePath);
    if (!stat.isFile() || stat.size > 20 * 1024 * 1024) throw new Error('图片最大支持 20 MB，请先缩小或裁剪。');
    const bytes = await fs.readFile(sourcePath);
    const { format, width, height } = inspectRaster(bytes);
    const directory = await fs.mkdtemp(path.join(sessionDirectory, 'raster-source-'));
    const sourceFile = path.join(directory, `source-image.${format === 'jpeg' ? 'jpg' : 'png'}`);
    await fs.writeFile(sourceFile, bytes);
    rasterSources.set(sourcePath, { directory, sourceFile, format, width, height, name: path.basename(sourcePath),
      sha256: createHash('sha256').update(bytes).digest('hex') });
    return { dataUrl: `data:image/${format};base64,${bytes.toString('base64')}`, name: path.basename(sourcePath), format };
  });
  ipcMain.handle('commit-raster', (_event, input) => exclusive(async () => {
    const source = rasterSources.get(input.sourcePath);
    if (!source) throw new Error('图片读取会话已失效，请重新打开。');
    const calibrated = calibrateRasterLoops(input);
    if (typeof input.overlayDataUrl !== 'string' || !input.overlayDataUrl.startsWith('data:image/png;base64,') || input.overlayDataUrl.length > 16 * 1024 * 1024)
      throw new Error('请先完成轮廓预览。');
    const directory = await fs.mkdtemp(path.join(sessionDirectory, 'raster-geometry-'));
    const geometryPath = path.join(directory, `${safeBaseName(source.name)}-contour.xy`);
    try {
      await fs.writeFile(geometryPath, geometry.loopsToXyText(calibrated.loops, 'Calibrated image outline; units=m'));
      const overlay = Buffer.from(input.overlayDataUrl.split(',')[1], 'base64');
      if (!overlay.subarray(0,8).equals(Buffer.from([137,80,78,71,13,10,26,10]))) throw new Error('轮廓预览图片无效。');
      await fs.writeFile(path.join(directory, 'image-outline.png'), overlay);
      const sourceName = path.basename(source.sourceFile);
      await fs.copyFile(source.sourceFile, path.join(directory, sourceName));
      const warnings = (input.warnings || []).filter(v => typeof v === 'string').slice(0,100);
      const report = { source: source.name, sourceSha256: source.sha256, sourceEncodedPixels: { width: source.width, height: source.height }, processedPixels: { width: input.pixelWidth, height: input.pixelHeight },
        calibration: input.calibration, metresPerProcessedPixel: calibrated.metresPerPixel,
        physicalWidth: calibrated.physicalWidth, physicalHeight: calibrated.physicalHeight,
        outputUnits: 'm', pointCount: calibrated.pointCount, loopCount: calibrated.loops.length,
        settings: input.settings, stats: input.stats, warnings,
        limits: 'Local 2D silhouette extraction, confirmed by user; no perspective correction or 3D reconstruction.' };
      await fs.writeFile(path.join(directory, 'image-import.json'), JSON.stringify(report, null, 2));
      rasterImports.set(geometryPath, { directory, warnings, files: [sourceName, 'image-outline.png', 'image-import.json'] });
      return { geometryPath, label: `${source.name} · 图片轮廓` };
    } catch (error) { await fs.rm(directory, { recursive: true, force: true }); throw error; }
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

  async function generateOnce(request, retainPrevious = false) {
    const { job, method } = validateJob(request);
    if (currentResult && !retainPrevious) await fs.rm(currentResult.outputDirectory, { recursive: true, force: true });
    currentResult = null;
    // Only smoke runs may override the temporary root. Every run owns a new folder.
    const root = process.argv.some(arg => arg.startsWith('--smoke=')) && job.outputDirectory
      ? job.outputDirectory : sessionDirectory;
    await fs.mkdir(root, { recursive: true });
    job.outputDirectory = await fs.mkdtemp(path.join(root, 'mesh-'));
    const prefix = path.join(job.outputDirectory, safeBaseName(job.geometryPath));
    const paths = { prefix, xyPath: `${prefix}.xy`, casePath: `${prefix}-openfoam` };

    const rasterImport = rasterImports.get(job.geometryPath);
    if (rasterImport) job.sourceUnits = 'm';
    const prepared = await prepareGeometry(job.geometryPath, job, paths.xyPath,
      `${prefix}.dxf.json`, log);
    if (rasterImport) for (const file of rasterImport.files)
      await fs.copyFile(path.join(rasterImport.directory, file), path.join(job.outputDirectory, file));
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
      rasterImport: Boolean(rasterImport),
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
    if (request.automatic && request.targetCells != null) {
      let payload;
      try { payload = await runBudget(request, budgetFrame(request, frame), {
        signal: operation.signal,
        generate: choice => generateOnce(choice, true),
        progress: ({ attempt, maximum, parameters }) => {
          mainWindow.webContents.send('run-progress', { attempt, maximum,
            estimatedSeconds: request.targetCells >= 100000 ? 90 : 30, estimateSource: '数量档位粗估' });
          log(`数量目标 ${request.targetCells}：第 ${attempt}/${maximum} 组，壁面 h/Lref=${parameters.wallRelativeSize}，背景=${parameters.backgroundRelativeSize}。`);
        }
      }); } catch (error) { currentResult = null; throw error; }
      currentResult = payload;
      await fs.writeFile(path.join(payload.outputDirectory, 'selection.json'), JSON.stringify({
        automatic: true, cellBudget: payload.cellBudget,
        selectedRequest: payload.selectedRequest, attempts: payload.attempts
      }, null, 2));
      return payload;
    }
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
    ['growthRatio', 'growth-ratio'],
    ['extrusionRelativeSize', 'extrusion-relative-size'],
    ['smallAlpha', 'small-alpha'],
    ['autoPadding', 'auto-padding'],
    ['customTargetCells', 'target-cells'],
  ].map(([id, flag]) => [id, argument(flag)]).filter(([, value]) => value !== null));
  if (outputDirectory) await fs.mkdir(outputDirectory, { recursive: true });
  if (argument('small-window')) mainWindow.setSize(800, 560);

  // The renderer's init awaits the catalog over IPC, so the hook appears a moment
  // after the page finishes loading.
  for (let attempt = 0; attempt < 100; attempt++) {
    if (await mainWindow.webContents.executeJavaScript('Boolean(window.__smoke)')) break;
    await new Promise(resolve => setTimeout(resolve, 100));
  }

  await mainWindow.webContents.executeJavaScript(`(async () => {
    const smoke = window.__smoke;
    smoke.setSidebarCollapsed(false);
    smoke.setOutput(${JSON.stringify(outputDirectory)});
    smoke.selectMethod(${JSON.stringify(method)});
    document.getElementById('controlMode').value = ${JSON.stringify(argument('control') || 'auto')};
    document.getElementById('controlMode').dispatchEvent(new Event('change'));
    document.getElementById('density').value = ${JSON.stringify(argument('target-cells') ? 'custom' : argument('density') || '5000')};
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
    if (${JSON.stringify(Boolean(argument('image')))}) {
      const previousPath = smoke.state.geometryPath;
      const pendingImport = smoke.importGeometryFile(${JSON.stringify(argument('image'))});
      const firstResult = await window.__rasterSmoke.ready;
      if (${JSON.stringify(Boolean(argument('image-reject')))}) {
        if (firstResult || !window.__rasterSmoke.state().confirmDisabled) throw new Error('Invalid image was accepted');
        const reason = window.__rasterSmoke.state().status;
        window.__rasterSmoke.cancel(); await pendingImport;
        if (smoke.state.geometryPath !== previousPath) throw new Error('Cancel discarded previous geometry');
        return { rasterPreview: true, rejected: true, reason, previousGeometryPreserved: true };
      }
      if (!firstResult) throw new Error(window.__rasterSmoke.state().status);
      if (!window.__rasterSmoke.state().confirmDisabled) throw new Error('Image accepted without physical calibration');
      const options = ${JSON.stringify({ mode: argument('image-mode') || 'auto', fillHoles: argument('image-fill') === 'true' })};
      await window.__rasterSmoke.setOptions(options);
      const extracted = window.__rasterSmoke.getResult();
      if (!extracted) throw new Error(window.__rasterSmoke.state().status);
      smoke.state.rasterEvidence = { loopCount: extracted.loops.length,
        vertices: extracted.loops.reduce((n,loop)=>n+loop.length,0), stats: extracted.stats, warnings: extracted.warnings };
      if (${JSON.stringify(Boolean(argument('image-preview-only')))})
        return { rasterPreview: true, ...smoke.state.rasterEvidence, calibrationRequired: window.__rasterSmoke.state().confirmDisabled };
      window.__rasterSmoke.setCalibration(Number(${JSON.stringify(argument('image-width') || '200')}), 'mm');
      const committed = await window.__rasterSmoke.confirm();
      if (!committed) throw new Error(window.__rasterSmoke.state()?.status || 'Raster commit failed');
      await pendingImport;
      const expectedWidth = Number(${JSON.stringify(argument('image-width') || '200')}) / 1000;
      if (Math.abs(smoke.state.frame.width - expectedWidth) > expectedWidth * 1e-10 ||
          !document.getElementById('sourceUnits').disabled) throw new Error('Image calibration was not preserved');
    }
    if (${JSON.stringify(argument('verified-preset') === 'true')}) await smoke.loadVerifiedPreset();
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
    document.getElementById('density').dispatchEvent(new Event('change'));
    const theme = document.getElementById('appTheme');
    theme.value = ${JSON.stringify(argument('theme') || 'modern')};
    theme.dispatchEvent(new Event('change'));
    const regionInput = document.querySelector('#regionList input');
    if (regionInput) {
      regionInput.focus();
      if (!regionInput.isConnected || document.activeElement !== regionInput) throw new Error('Region input lost focus');
    }
    if (${JSON.stringify(Boolean(argument('welcome-shot')))}) return { welcomeOnly: true };
    const pending = smoke.generate();
    if (!document.getElementById('sample').disabled || document.getElementById('cancel').hidden)
      throw new Error('Parameters are not locked during generation');
    await pending;
    if (${JSON.stringify(Boolean(argument('interaction-check')))}) {
      const mesh = smoke.state.mesh;
      const display = document.getElementById('displayMode').value;
      const sidebar = document.getElementById('toggleSidebar');
      const initiallyCollapsed = document.body.classList.contains('sidebar-collapsed');
      sidebar.click();
      if (document.body.classList.contains('sidebar-collapsed') === initiallyCollapsed ||
          sidebar.getAttribute('aria-expanded') !== String(!document.body.classList.contains('sidebar-collapsed')) ||
          localStorage.getItem('cartmesh2d-sidebar-collapsed') !== String(document.body.classList.contains('sidebar-collapsed'))) {
        throw new Error('Sidebar toggle did not update state or persistence');
      }
      sidebar.click();
      if (document.body.classList.contains('sidebar-collapsed') !== initiallyCollapsed)
        throw new Error('Sidebar toggle did not restore state');
      for (const value of ['duet','modern']) {
        window.CartMeshTheme.set(value);
        if (smoke.state.mesh !== mesh || document.getElementById('displayMode').value !== display ||
            localStorage.getItem(window.CartMeshTheme.storageKey) !== value)
          throw new Error('Theme switch changed mesh/display state or failed persistence');
      }
      window.CartMeshTheme.set(${JSON.stringify(argument('theme') || 'modern')});
      sidebar.focus();
      window.dispatchEvent(new KeyboardEvent('keydown', { key: 'b', metaKey: true }));
      if (document.body.classList.contains('sidebar-collapsed') === initiallyCollapsed)
        throw new Error('Sidebar Command+B shortcut did not toggle');
      window.dispatchEvent(new KeyboardEvent('keydown', { key: 'b', ctrlKey: true }));
      if (document.body.classList.contains('sidebar-collapsed') !== initiallyCollapsed)
        throw new Error('Sidebar Control+B shortcut did not restore state');
      const actual = smoke.state.selectedRequest;
      if (actual) {
      document.getElementById('actualToManual').click();
      if (document.getElementById('controlMode').value !== 'manual' ||
          Number(document.getElementById('wallRelativeSize').value) !== actual.wallRelativeSize ||
          Number(document.getElementById('relativeBandCells').value) !== actual.cellsPerLevel ||
          Number(document.getElementById('relativePadding').value) !== actual.farFieldSpans)
        throw new Error('Actual-to-manual transfer lost parameters');
      document.getElementById('controlMode').value = 'auto';
      document.getElementById('controlMode').dispatchEvent(new Event('change'));
      } else if (document.getElementById('controlMode').value !== 'manual') {
        throw new Error('Missing selected automatic parameters');
      }
    }
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
      log: document.getElementById('log').textContent,
      cellBudget: smoke.state.cellBudget,
      raster: smoke.state.rasterEvidence,
      theme: document.documentElement.dataset.theme,
      interactionChecks: ${JSON.stringify(Boolean(argument('interaction-check')))}
    };
  })()`).then(async report => {
    if (report.rasterPreview) {
      await new Promise(resolve => setTimeout(resolve, 200));
      report.dialogLayout = await mainWindow.webContents.executeJavaScript(`(() => {
        const dialog = document.querySelector('.raster-import');
        if (!dialog) return null;
        const confirm = dialog.querySelector('[data-raster="confirm"]');
        const box = dialog.getBoundingClientRect();
        const button = confirm?.getBoundingClientRect();
        return { width: innerWidth, height: innerHeight, visible: box.top >= 0 && box.bottom <= innerHeight,
          confirmReachable: Boolean(button && button.top >= 0 && button.bottom <= innerHeight) };
      })()`);
      if (shot) {
        await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
        await fs.writeFile(shot + '.json', JSON.stringify(report, null, 2));
      }
      console.log(JSON.stringify(report, null, 2));
      await mainWindow.webContents.executeJavaScript('window.__rasterSmoke.cancel()');
      await fs.rm(sessionDirectory, { recursive: true, force: true });
      app.exit(0); return;
    }
    if (report.welcomeOnly) {
      await mainWindow.webContents.executeJavaScript(`window.__smoke.state.mesh = null; document.getElementById('empty').hidden = false; window.__smoke.view.clear();`);
      await mainWindow.webContents.executeJavaScript(`(async () => {
        const artwork = new Image(); artwork.src = 'assets/duet-workshop.png'; await artwork.decode();
      })()`);
      // Wake the compositor before capturing a newly loaded CSS wallpaper.
      await mainWindow.webContents.capturePage(undefined, { stayAwake: true });
      await new Promise(resolve => setTimeout(resolve, 1000));
      await fs.writeFile(argument('welcome-shot'), (await mainWindow.webContents.capturePage(undefined, { stayAwake: true })).toPNG());
      await fs.rm(sessionDirectory, { recursive: true, force: true });
      app.exit(0); return;
    }
    if (argument('export')) report.exported = await exportPackage(argument('export'));
    if (argument('mesh-shot')) await fs.writeFile(argument('mesh-shot'), (await mainWindow.webContents.capturePage()).toPNG());
    if (argument('home-check')) {
      report.previewTiming = await mainWindow.webContents.executeJavaScript(`(() => {
        const view = window.__smoke.view;
        const original = view.scale;
        const elapsed = [];
        for (let i = 0; i < 6; i++) {
          view.scale = original * (1 + i * .01);
          const start = performance.now(); view.draw();
          elapsed.push(performance.now() - start);
        }
        view.scale = original; view.draw();
        return { drawSubmissionMs: elapsed, note: 'Canvas command submission; excludes asynchronous GPU completion' };
      })()`);
      report.home = await mainWindow.webContents.executeJavaScript(`(async () => {
        await window.__smoke.returnToStart();
        const { state, view } = window.__smoke;
        return { emptyVisible: !document.getElementById('empty').hidden,
          released: !state.mesh && !state.result && !view.mesh && !view.meshCache && !view.outline,
          geometryCleared: !state.geometryPath,
          exportAvailable: !document.getElementById('exportResult').hidden };
      })()`);
      report.home.mainPreviewReleased = currentResult?.mesh === null;
      if (!Object.values(report.home).every(Boolean)) throw new Error('Return to start did not release preview or preserve export');
      if (argument('export')) report.home.exported = await exportPackage(argument('export').replace(/\.zip$/, '') + '-home.zip');
    }
    mainWindow.setSize(800, 560);
    await new Promise(resolve => setTimeout(resolve, 200));
    report.layout = await mainWindow.webContents.executeJavaScript(`(() => {
      const panel = document.querySelector('.panel');
      panel.scrollTop = panel.scrollHeight;
      const button = document.getElementById('generate').getBoundingClientRect();
      const bodyRect = document.body.getBoundingClientRect();
      const shellRect = document.querySelector('.shell').getBoundingClientRect();
      const canvasRect = document.getElementById('canvas').getBoundingClientRect();
      const geometry = { bodyHeight: bodyRect.height, innerHeight,
        canvasHeight: canvasRect.height, shellBottom: shellRect.bottom };
      return { bottomReachable: button.bottom <= innerHeight && button.top >= 0,
        pageHeight: document.documentElement.scrollHeight, windowHeight: innerHeight,
        geometry,
        wallpaper: {
          image: getComputedStyle(document.body).backgroundImage,
          canvasBackground: getComputedStyle(document.getElementById('canvasWrap')).backgroundColor,
          welcomeHidden: document.getElementById('empty').hidden
        },
      previewCells: window.__smoke.state.mesh?.cells.length,
      exportedCells: window.__smoke.state.result?.openFoam.cells,
      sidebarCollapsed: document.body.classList.contains('sidebar-collapsed'),
      sidebarExpanded: document.getElementById('toggleSidebar')?.getAttribute('aria-expanded') === 'true' };
    })()`);
    console.log(JSON.stringify(report, null, 2));
    if (!report.layout.bottomReachable) throw new Error('Sidebar bottom is inaccessible');
    if (Math.abs(report.layout.geometry.bodyHeight - report.layout.geometry.innerHeight) > 1 ||
        report.layout.geometry.canvasHeight <= 0 || report.layout.geometry.shellBottom > report.layout.geometry.innerHeight + 1)
      throw new Error('Window geometry is not constrained to the viewport');
    if (shot) {
      await mainWindow.webContents.executeJavaScript("document.querySelector('.panel').scrollTop=0");
      await new Promise(resolve => setTimeout(resolve, 400));
      await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
      console.log(`screenshot=${shot}`);
      await fs.writeFile(shot + '.json', JSON.stringify(report, null, 2));
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
