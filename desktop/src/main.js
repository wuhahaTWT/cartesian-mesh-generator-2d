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
const { parseBackgroundGrid, requireFluidMesh } = require('./core/background-grid');
const { exportGuide } = require('./core/export-guide');
const { zipDirectory } = require('./core/archive');
const { FLOW_CASES, FLOW_CONVECTION_SCHEMES, FLOW_PRESSURE_PRECONDITIONERS, FLOW_OUTLET_BACKFLOW_MODES, buildFlowInvocation, commitFlowFiles,
        parseFlowProgress, validateFlowOutput, validateTimeHistory, validateAttemptHistory, flowOutputSuffixes } = require('./core/flow');
const {validateThermalRequest,thermalCheckpointTime}=require('./core/thermal');
const { runThermalJob } = require('./core/thermal-job');
const { runEulerJob, importEulerRestart } = require('./core/euler-job');
const { readCheckpointMetadata } = require('./core/flow-checkpoint');
const { pressureDrivenBoundaryDefinition, parseBoundaryDefinition, serializeBoundaryDefinition, validateBoundaryMesh, sameConditions } = require('./core/flow-boundaries');
const { MAX_BYTES: FLOW_CASE_MAX_BYTES, createFlowCaseDocument, serializeFlowCase, parseFlowCaseDocument } = require('./core/flow-case');
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
  if (currentResult.thermal) {
    const thermalPng=await mainWindow.webContents.executeJavaScript('window.__exportThermalPreview()');
    if(typeof thermalPng!=='string'||!thermalPng.startsWith('data:image/png;base64,'))throw new Error('温度结果图片生成失败。');
    await fs.writeFile(path.join(currentResult.outputDirectory,'temperature-preview.png'),Buffer.from(thermalPng.split(',')[1],'base64'));
  }
  if (currentResult.euler) {
    const png=await mainWindow.webContents.executeJavaScript('window.__exportEulerPreview()');
    if(typeof png!=='string'||!png.startsWith('data:image/png;base64,'))throw new Error('可压结果图片生成失败。');
    await fs.writeFile(path.join(currentResult.outputDirectory,'euler-preview.png'),Buffer.from(png.split(',')[1],'base64'));
  }
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
    samples: SAMPLES.map(sample => ({ ...sample, path: resourcePath('samples', sample.file) })),
    flowCases: FLOW_CASES,
    flowConvectionSchemes: FLOW_CONVECTION_SCHEMES,
    flowPressurePreconditioners: FLOW_PRESSURE_PRECONDITIONERS,
    flowOutletBackflowModes: FLOW_OUTLET_BACKFLOW_MODES
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
    return { mesh: currentResult.mesh || (currentResult.background ? parseBackgroundGrid(await fs.readFile(currentResult.backgroundPath,'utf8')) : parseCm2d(await fs.readFile(currentResult.cm2dPath, 'utf8'))),
      result: currentResult.result, thermal: currentResult.thermal || null, euler:currentResult.euler || null };
  });
  ipcMain.handle('export-result', () => exclusive(async () => {
    if (!currentResult) throw new Error('请先成功生成网格。');
    const result = await dialog.showSaveDialog(mainWindow, {
      title: '保存网格结果包', defaultPath: `${safeBaseName(currentResult.job.geometryPath)}.zip`,
      filters: [{ name: '网格结果包', extensions: ['zip'] }]
    });
    return result.canceled ? null : exportPackage(result.filePath);
  }));

  const flowState = () => ({ flow: currentResult?.flow || null, restart: currentResult?.flowRestart?.metadata || null });
  ipcMain.handle('flow-state', () => flowState());
  // Only the real smoke harness can supply a dialog replacement. Renderer
  // requests never carry read/write paths for case files.
  const smokeCasePath = () => {
    if (!process.argv.includes('--flow-case-check=true') || !process.argv.some(a => a.startsWith('--smoke='))) return null;
    const out = process.argv.find(a => a.startsWith('--out='))?.slice(6);
    if (!out || !path.isAbsolute(out)) throw new Error('工况验证需要绝对输出目录。');
    return path.join(out, 'saved-flow-case.json');
  };
  ipcMain.handle('save-flow-case', (_event, request) => exclusive(async () => {
    requireFluidMesh(currentResult);
    const document = createFlowCaseDocument(request, await fs.readFile(currentResult.cm2dPath));
    const text = serializeFlowCase(document);
    let file = smokeCasePath();
    if (!file) {
      const picked = await dialog.showSaveDialog(mainWindow, { title: '保存流动工况',
        defaultPath: `${safeBaseName(currentResult.job.geometryPath)}.flow-case.json`,
        filters: [{ name: '流动工况', extensions: ['json'] }] });
      if (picked.canceled) return null;
      file = picked.filePath;
    }
    const temporaryDirectory = await fs.mkdtemp(path.join(path.dirname(file), '.cartmesh2d-case-'));
    try {
      const temporary = path.join(temporaryDirectory, 'case.json');
      await fs.writeFile(temporary, text, { flag: 'wx' }); await fs.rename(temporary, file);
    } finally { await fs.rm(temporaryDirectory, { recursive: true, force: true }); }
    return { fileName: path.basename(file), document };
  }));
  ipcMain.handle('load-flow-case', () => exclusive(async () => {
    requireFluidMesh(currentResult);
    let file = smokeCasePath();
    if (!file) {
      const picked = await dialog.showOpenDialog(mainWindow, { title: '读取流动工况', properties: ['openFile'],
        filters: [{ name: '流动工况', extensions: ['json'] }] });
      if (picked.canceled) return null;
      file = picked.filePaths[0];
    }
    if ((await fs.stat(file)).size > FLOW_CASE_MAX_BYTES) throw new Error('工况文件超过32MB。');
    const document = parseFlowCaseDocument(await fs.readFile(file, 'utf8'), await fs.readFile(currentResult.cm2dPath));
    return { fileName: path.basename(file), document };
  }));
  ipcMain.handle('prepare-flow-boundaries', (_event, request) => exclusive(async () => {
    requireFluidMesh(currentResult);
    const speed=Number(request?.speed);
    if (!(Number.isFinite(speed) && speed>0)) throw new Error('参考速度必须大于零。');
    const mesh=currentResult.mesh || parseCm2d(await fs.readFile(currentResult.cm2dPath,'utf8'));
    let text;
    if (request.source==='import') {
      const picked=await dialog.showOpenDialog(mainWindow,{title:'导入与当前最终网格对应的命名边界',
        properties:['openFile'],filters:[{name:'命名流动边界',extensions:['boundaries']}]});
      if (picked.canceled) return null;
      if ((await fs.stat(picked.filePaths[0])).size>32*1024*1024) throw new Error('边界文件超过32MB。');
      text=await fs.readFile(picked.filePaths[0],'utf8');
    } else {
      if (!['channel','duct','pressure-duct','pressure-half-channel','cavity','annulus'].includes(request.source)) throw new Error('请选择支持的边界预设。');
      const directory=await fs.mkdtemp(path.join(currentResult.outputDirectory,'boundary-input-'));
      const target=path.join(directory,'input.boundaries');
      try {
        await runProcess(executable('cartmesh2d_flow_cli'),['--mesh',currentResult.cm2dPath,'--case',request.source==='pressure-half-channel'?'channel':request.source==='pressure-duct'?'duct':request.source,
          '--speed',String(speed),'--export-boundaries',target],()=>{},operation.signal,30000);
        text=await fs.readFile(target,'utf8');
      } finally { await fs.rm(directory,{recursive:true,force:true}); }
    }
    const definition=parseBoundaryDefinition(text);
    return validateBoundaryMesh(['pressure-duct','pressure-half-channel'].includes(request.source) ? pressureDrivenBoundaryDefinition(definition,request.source==='pressure-half-channel') : definition,mesh,speed);
  }));
  ipcMain.handle('pick-flow-checkpoint', () => exclusive(async () => {
    requireFluidMesh(currentResult);
    const picked = await dialog.showOpenDialog(mainWindow, { title: '选择非定常重启状态',
      properties: ['openFile'], filters: [{ name: '已接受流动状态', extensions: ['checkpoint'] }] });
    if (picked.canceled) return null;
    const file = picked.filePaths[0];
    const metadata = await readCheckpointMetadata(file);
    currentResult.flowRestart = { path: file, metadata };
    return metadata;
  }));
  ipcMain.handle('run-flow', (_event, request) => exclusive(async () => {
    requireFluidMesh(currentResult);
    const mesh = currentResult.mesh
      || assignSizeBands(parseCm2d(await fs.readFile(currentResult.cm2dPath, 'utf8')));
    const selectedRestart = request?.resume ? currentResult.flowRestart : null;
    // Validate before creating outputs; paths come only from this main process.
    buildFlowInvocation(currentResult.cm2dPath, 'pending', request, selectedRestart?.path, 'pending.boundaries');
    const boundaryDefinition=request.case==='custom' ? validateBoundaryMesh(request.boundaryDefinition,mesh,Number(request.speed)) : null;
    const incompleteDirectory = await fs.mkdtemp(path.join(currentResult.outputDirectory, 'flow-incomplete-'));
    const pendingPrefix = path.join(incompleteDirectory, 'flow');
    let restartPath = null, startTime = 0;
    if (selectedRestart) {
      try {
        restartPath = path.join(incompleteDirectory, 'input.checkpoint');
        await fs.copyFile(selectedRestart.path, restartPath);
        const metadata = await readCheckpointMetadata(restartPath);
        startTime = metadata.time;
        for (const key of ['case', 'nu', 'speed', 'convection', 'outletBackflow'])
          if (metadata[key] !== (['nu','speed'].includes(key) ? Number(request[key]) : request[key]))
            throw new Error('续算必须保持原工况、物性和对流格式；可调整时间步与步数。');
        if (boundaryDefinition && !sameConditions(metadata.boundaryDefinition?.records,boundaryDefinition.records))
          throw new Error('续算必须保持原命名边界的名称、类型和数值。');
      } catch (error) {
        // No solver has started and the selected source remains untouched.
        await fs.rm(incompleteDirectory, { recursive: true, force: true });
        throw error;
      }
    }
    const boundaryPath=boundaryDefinition ? path.join(incompleteDirectory,'input.boundaries') : null;
    if (boundaryPath) await fs.writeFile(boundaryPath,serializeBoundaryDefinition(boundaryDefinition));
    const invocation = buildFlowInvocation(currentResult.cm2dPath, pendingPrefix, request, restartPath, boundaryPath);
    const transient = invocation.request.mode !== 'steady';
    const adaptive = invocation.request.mode === 'adaptive';
    if (adaptive && !(invocation.request.endTime>startTime)) {
      await fs.rm(incompleteDirectory,{recursive:true,force:true});
      throw new Error('目标物理时间必须晚于已接受的重启时间。');
    }
    const previousFlow = currentResult.flow;
    const previousRestart = currentResult.flowRestart;
    const preserveIncomplete = async error => {
      if (transient) {
        // Ignore .tmp: only the native atomic accepted-state file is resumable.
        try { currentResult.flowRestart = { path: `${pendingPrefix}.checkpoint`,
          metadata: await readCheckpointMetadata(`${pendingPrefix}.checkpoint`) }; }
        catch { currentResult.flowRestart = previousRestart; }
      }
      const report = { format: 'cartmesh2d-flow-incomplete-v1',
        status: operation.signal.aborted ? 'cancelled' : 'failed', request: invocation.request,
        mesh: path.basename(currentResult.cm2dPath),
        acceptedTime: currentResult.flowRestart?.metadata.time ?? null,
        exitCode: Number.isInteger(error.code) ? error.code : null,
        message: String(error.message || error).split('\n')[0] };
      await fs.writeFile(path.join(incompleteDirectory, 'desktop-flow-error.json'), JSON.stringify(report, null, 2));
    };
    log(`正在运行原生二维${transient ? '非定常' : '稳态'}层流：${FLOW_CASES[invocation.request.case].label}…`);
    const onLine = (line, isError) => {
      let progress = null;
      if (!isError) {
        try { progress = parseFlowProgress(line); }
        catch (error) { log(`忽略无效进度：${error.message}`); }
      }
      if (progress) mainWindow.webContents.send('flow-progress', progress);
      else log(line);
    };
    let backups = [];
    let commitStarted = false;
    try {
      const processResult = await runProcess(executable(invocation.executable), invocation.args,
        onLine, operation.signal, 0, [0, 2]);
      operation.signal.throwIfAborted();
      const outputFiles = { summary: `${pendingPrefix}.json`, fields: `${pendingPrefix}.fields.json`,
        vtk: `${pendingPrefix}.vtk`, residuals: `${pendingPrefix}.residuals.csv`,
        cells: `${pendingPrefix}.cells.csv`, faces: `${pendingPrefix}.faces.csv` };
      if (transient) Object.assign(outputFiles, { checkpoint: `${pendingPrefix}.checkpoint`, timeHistory: `${pendingPrefix}.time-history.csv` });
      if (adaptive) outputFiles.attemptHistory=`${pendingPrefix}.attempt-history.csv`;
      if (invocation.request.initialVortex) outputFiles.initialCheckpoint=`${pendingPrefix}.initial.checkpoint`;
      if (boundaryDefinition) outputFiles.boundaries=`${pendingPrefix}.boundaries`;
      const [summary, fields] = await Promise.all([readJson(outputFiles.summary), readJson(outputFiles.fields),
        ...Object.values(outputFiles).map(file => fs.stat(file))]);
      const validated = validateFlowOutput(summary, fields, mesh.cells.length, invocation.request, startTime);
      if (boundaryDefinition) {
        const exported=validateBoundaryMesh(parseBoundaryDefinition(await fs.readFile(outputFiles.boundaries,'utf8')),mesh,invocation.request.speed);
        if (!sameConditions(exported.records,boundaryDefinition.records)) throw new Error('导出边界与输入不一致。');
      }
      if ((processResult.code === 0) !== validated.summary.converged)
        throw new Error('原生求解器退出码与收敛状态不一致。');
      let history = null, attempts = null, checkpointMetadata = null;
      if (transient) {
        history = validateTimeHistory(await fs.readFile(outputFiles.timeHistory, 'utf8'), validated.summary, startTime);
        if (adaptive) attempts=validateAttemptHistory(await fs.readFile(outputFiles.attemptHistory,'utf8'),validated.summary,history);
        checkpointMetadata = await readCheckpointMetadata(outputFiles.checkpoint);
        if (outputFiles.initialCheckpoint) {
          const initial=await readCheckpointMetadata(outputFiles.initialCheckpoint);
          if (initial.time!==0 || ['case','nu','speed','convection','outletBackflow'].some(k=>initial[k]!==checkpointMetadata[k]))
            throw new Error('初始局部涡检查点的时间或物性与结果不一致。');
        }
        if (Math.abs(checkpointMetadata.time-summary.acceptedTime) > 1e-12+1e-9*Math.abs(summary.acceptedTime))
          throw new Error('重启状态时间与摘要不一致。');
        if (!validated.summary.converged)
          throw Object.assign(new Error(`时间步未收敛；已接受到 t=${summary.acceptedTime} s，可继续计算。候选场仅留作诊断。`), { code: 2 });
      }
      operation.signal.throwIfAborted();
      // Preserve the earlier complete result even if copying the new set fails.
      const allSuffixes = flowOutputSuffixes({ mode: 'adaptive', case:'custom', initialVortex:true });
      for (const suffix of allSuffixes) {
        const destination = `${currentResult.prefix}.flow${suffix}`;
        const backup = path.join(incompleteDirectory, `previous${suffix}`);
        try { await fs.copyFile(destination, backup); backups.push({ backup, destination }); }
        catch (error) { if (error.code !== 'ENOENT') throw error; }
      }
      const entries = Object.entries(outputFiles).map(([kind, source]) => ({ kind, source,
        destination: `${currentResult.prefix}.flow${source.slice(pendingPrefix.length)}` }));
      commitStarted = true;
      await commitFlowFiles(fs, entries);
      operation.signal.throwIfAborted();
      for (const suffix of allSuffixes.filter(suffix => !flowOutputSuffixes(invocation.request).includes(suffix)))
        await fs.rm(`${currentResult.prefix}.flow${suffix}`, { force: true });
      const saved = Object.fromEntries(entries.map(entry => [entry.kind, path.basename(entry.destination)]));
      const payload = { ...validated, request: invocation.request, files: saved, history, attempts };
      currentResult.flow = payload;
      currentResult.flowRestart = transient ? { path: `${currentResult.prefix}.flow.checkpoint`, metadata: { ...checkpointMetadata, fileName: path.basename(`${currentResult.prefix}.flow.checkpoint`) } } : null;
      await fs.rm(incompleteDirectory, { recursive: true, force: true }).catch(error => log(`结果已保存，临时目录清理失败：${error.message}`));
      log(transient ? `非定常计算完成，已接受到 t=${summary.acceptedTime} s。`
        : validated.summary.converged ? '层流求解已收敛。' : '层流求解到达迭代上限，保留诊断结果但未收敛。');
      return payload;
    } catch (error) {
      if (commitStarted) {
        for (const suffix of flowOutputSuffixes({ mode: 'adaptive', case:'custom', initialVortex:true }))
          await fs.rm(`${currentResult.prefix}.flow${suffix}`, { force: true }).catch(() => {});
      }
      const restored = await Promise.allSettled(backups.map(entry => fs.copyFile(entry.backup, entry.destination)));
      if (restored.some(entry => entry.status === 'rejected')) {
        currentResult.flow = null;
        error.message += '\n上次结果恢复失败，完整备份保留在诊断目录。';
      }
      if (restored.every(entry => entry.status === 'fulfilled')) currentResult.flow = previousFlow;
      await preserveIncomplete(error).catch(() => {});
      error.message += `\n未完成诊断保留在 ${incompleteDirectory}`;
      throw error;
    }
  }));

  ipcMain.handle('euler-state', () => ({euler:currentResult?.euler || null,restart:currentResult?.eulerRestart?.metadata || null}));
  ipcMain.handle('run-euler', (_event,request) => exclusive(async()=>{
    requireFluidMesh(currentResult);
    const mesh=currentResult.mesh||parseCm2d(await fs.readFile(currentResult.cm2dPath,'utf8'));
    return runEulerJob({currentResult,mesh,request,executable,runProcess,signal:operation.signal,
      onProgress:progress=>mainWindow.webContents.send('euler-progress',progress),log:line=>mainWindow.webContents.send('run-line',line)});
  }));
  ipcMain.handle('pick-euler-checkpoint',()=>exclusive(async()=>{
    requireFluidMesh(currentResult);
    const selected=await dialog.showOpenDialog(mainWindow,{title:'选择 Euler 结果目录中的 desktop-state.json',filters:[{name:'Euler 续算清单',extensions:['json']}],properties:['openFile']});
    if(selected.canceled)return null;
    const mesh=currentResult.mesh||parseCm2d(await fs.readFile(currentResult.cm2dPath,'utf8'));
    currentResult.eulerRestart=await importEulerRestart(selected.filePaths[0],mesh,currentResult.cm2dPath);
    return currentResult.eulerRestart.metadata;
  }));

  ipcMain.handle('thermal-state', () => ({thermal:currentResult?.thermal || null,restart:currentResult?.thermalRestart?.metadata || null}));
  ipcMain.handle('pick-thermal-checkpoint',()=>exclusive(async()=>{
    requireFluidMesh(currentResult);
    const picked=await dialog.showOpenDialog(mainWindow,{title:'载入导出包中的联合续算状态',properties:['openFile'],filters:[{name:'联合状态',extensions:['checkpoint']}]});
    if(picked.canceled)return null;
    const file=picked.filePaths[0];
    if(!file.endsWith('.thermal.checkpoint'))throw new Error('请选择 thermal.checkpoint，carrier.checkpoint 不能联合续算。');
    const saved=await readJson(path.join(path.dirname(file),'desktop-state.json'));
    if(!saved?.request)throw new Error('请保留同目录的 desktop-state.json，以恢复物性和热边界。');
    const request=validateThermalRequest(saved.request),time=thermalCheckpointTime(await fs.readFile(file,'utf8'));
    currentResult.thermalRestart={path:file,metadata:{time,request,fileName:path.basename(file)}};
    return currentResult.thermalRestart.metadata;
  }));
  ipcMain.handle('run-thermal', (_event,request) => exclusive(async()=>{
    requireFluidMesh(currentResult);
    const mesh=currentResult.mesh || parseCm2d(await fs.readFile(currentResult.cm2dPath,'utf8'));
    log('正在同步推进原生流动与温度；温度不反馈物性或浮力。');
    return runThermalJob({currentResult,mesh,request,executable,runProcess,signal:operation.signal,
      onProgress:progress=>mainWindow.webContents.send('thermal-progress',progress),log});
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
    if(job.method === 'background') {
      if(failure) throw failure;
      const parsed=parseBackgroundGrid(await fs.readFile(invocation.backgroundPath,'utf8'));
      const values=parseKeyValues(stdout);
      const result={requestedMethod:'background',actualMethod:'background',background:true,
        counts:{cells:parsed.cells.length,vertices:4*parsed.cells.length,faces:0,leaves:parsed.cells.length,
          cutCells:0,layerCells:0,classification:parsed.classificationCounts},
        gates:{topology:{pass:null,label:'背景网格'},solver:null},
        openFoam:{written:false,path:null,cells:0},
        timings:{total_seconds:Number(values.total_seconds)},raw:values};
      const payload={job,prefix,outputDirectory:job.outputDirectory,background:true,
        backgroundPath:invocation.backgroundPath,cm2dPath:null,mesh:parsed,
        levelBasis:'level',levelHistogram:levelHistogram(parsed),wallBounds:null,
        incomplete:null,rasterImport:Boolean(rasterImport),result};
      await fs.writeFile(path.join(job.outputDirectory,'result.json'),JSON.stringify({
        parameters:{...job,geometryPath:path.basename(job.geometryPath),outputDirectory:undefined},
        invocation:{executable:invocation.executable,args:invocation.args},result},null,2));
      currentResult=payload;
      log('完整笛卡尔背景网格已生成；内部单元保留，未生成流体求解拓扑。');
      return payload;
    }
    const reports = await collectReports(job.method, prefix);
    const mesh = await firstReadable(invocation.cm2dCandidates);
    if (failure) {
      const failedQuality = await readJson(`${prefix}.failed.solver-quality.json`);
      failure.solverQualityFailure = job.method === 'cutcell' && failure.code === 1 &&
        /solver-quality gate failed/.test(failure.stderr || '') &&
        failedQuality?.quality_class === 'solver_quality' && failedQuality.valid === false;
      failure.outputDirectory = job.outputDirectory;
      await fs.writeFile(path.join(job.outputDirectory, 'generation-failure.json'), JSON.stringify({
        parameters: job, invocation, reason: failure.message,
        solverQualityFailure: failure.solverQualityFailure,
        stdout: failure.stdout || '', stderr: failure.stderr || ''
      }, null, 2));
    }
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
    if(request.method==='background') return generateOnce({...request,automatic:false,targetCells:undefined});
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
      }); } catch (error) {
        currentResult = null;
        if (error.outputDirectory && Array.isArray(error.attempts))
          await fs.writeFile(path.join(error.outputDirectory, 'selection-failure.json'), JSON.stringify({
            automatic: true, targetCells: request.targetCells, attempts: error.attempts
          }, null, 2));
        throw error;
      }
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

// `--smoke=<sample-id> [--out=<dir>] [--method=<id>] [--flow-convection=<scheme>] [--flow-pressure-preconditioner=<id>] [--shot=<png>]` drives the real
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
    ['backgroundLevel','background-level'],
    ['backgroundMinimumLevel','background-minimum-level'],
    ['backgroundPadding','background-padding'],
    ['backgroundMode','background-mode'],
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
    const loadedFonts = await document.fonts.load('13px "CartMesh UI"', '生成网格');
    await document.fonts.ready;
    if (!loadedFonts.length || loadedFonts.some(font => font.status !== 'loaded'))
      throw new Error('Bundled Chinese interface font did not load');
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
    const sampleField = sample.fluidRegion === 'interior' && sample.interiorSizeField || sample.sizeField;
    if (sampleField?.farFieldSpans !== undefined &&
        Number(document.getElementById('autoPadding').value) !== sampleField.farFieldSpans)
      throw new Error('Loading sample did not apply its physical domain to automatic mode');
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
      if (!input.checkValidity() || (input.type==='number' && !Number.isFinite(Number(input.value))))
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
    if (${JSON.stringify(Boolean(argument('outlet-backflow')))}) {
      document.getElementById('flowOutletBackflow').value = ${JSON.stringify(argument('outlet-backflow') || 'reject')};
      document.getElementById('flowOutletBackflow').dispatchEvent(new Event('change'));
    }
    if (${JSON.stringify(Boolean(argument('flow')))}) {
      document.getElementById('flowCase').value = ${JSON.stringify(argument('flow') || 'external')};
      document.getElementById('flowCase').dispatchEvent(new Event('change'));
      document.getElementById('flowMaxIterations').value = ${JSON.stringify(argument('flow-max-iterations') || '20')};
      document.getElementById('flowTolerance').value = ${JSON.stringify(argument('flow-tolerance') || '0.000001')};
      document.getElementById('flowNu').value = ${JSON.stringify(argument('flow-nu') || '0.01')};
      document.getElementById('flowSpeed').value = ${JSON.stringify(argument('flow-speed') || '1')};
      document.getElementById('flowConvection').value = ${JSON.stringify(argument('flow-convection') || 'upwind')};
      document.getElementById('flowConvection').dispatchEvent(new Event('change'));
      document.getElementById('flowPressurePreconditioner').value = ${JSON.stringify(argument('flow-pressure-preconditioner') || 'ic0')};
      document.getElementById('flowPressurePreconditioner').dispatchEvent(new Event('change'));
      document.getElementById('flowSteadyAcceleration').value = ${JSON.stringify(argument('flow-steady-acceleration') || 'none')};
      document.getElementById('flowLinearPolicy').value = ${JSON.stringify(argument('flow-linear-policy') || 'strict')};
      document.getElementById('flowVelocityRelaxation').value = ${JSON.stringify(argument('flow-velocity-relaxation') || '.6')};
      document.getElementById('flowPressureCorrections').value = ${JSON.stringify(argument('flow-pressure-corrections') || '4')};
      document.getElementById('flowSteadyAcceleration').dispatchEvent(new Event('change'));
      if (${JSON.stringify(argument('flow') === 'custom')}) {
        await smoke.prepareFlowBoundaries(${JSON.stringify(argument('flow-boundary-preset') || 'duct')});
        if (!smoke.state.flowBoundaryDefinition) throw new Error('Custom boundary editor did not receive definition');
        const pressure=document.querySelector('[data-boundary-key="p"][data-patch-name="outlet"]') || document.querySelector('[data-boundary-key="p"]');
        if (${JSON.stringify(argument('flow-outlet-pressure') !== null)}) {
          if (!pressure) throw new Error('No explicit pressure patch in editor');
          pressure.value=${JSON.stringify(argument('flow-outlet-pressure') || '0')};
          pressure.dispatchEvent(new Event('input'));
        }
        if (${JSON.stringify(argument('flow-opening-pressure') !== null)}) {
          const inlet=document.querySelector('[data-boundary-key="p"][data-patch-name="inlet"]');
          if(!inlet)throw new Error('No pressure-driven inlet in editor');
          inlet.value=${JSON.stringify(argument('flow-opening-pressure') || '0')};
          inlet.dispatchEvent(new Event('input'));
        }
        document.querySelector('#flowBoundaryPatches button').click();
        if (!smoke.view.boundaryHighlight.length) throw new Error('Patch location was not highlighted');
      }
      if (${JSON.stringify(Boolean(argument('flow-dt')))}) {
        document.getElementById('flowMode').value=${JSON.stringify(argument('flow-end-time') ? 'adaptive' : 'transient')};
        document.getElementById('flowMode').dispatchEvent(new Event('change'));
        document.getElementById('flowDt').value=${JSON.stringify(argument('flow-dt') || '.01')};
        document.getElementById('flowSteps').value=${JSON.stringify(argument('flow-steps') || '2')};
        document.getElementById('flowEndTime').value=${JSON.stringify(argument('flow-end-time') || '1')};
        document.getElementById('flowMaxCourant').value=${JSON.stringify(argument('flow-max-courant') || '1')};
      }
      if (${JSON.stringify(argument('flow-vortex-speed')!==null)}) {
        document.getElementById('flowInitialSettings').open=true;
        document.getElementById('flowInitialVortex').checked=true;
        document.getElementById('flowInitialVortex').dispatchEvent(new Event('change'));
        for(const [id,value] of Object.entries({flowVortexX:${JSON.stringify(argument('flow-vortex-x') || '3')},
          flowVortexY:${JSON.stringify(argument('flow-vortex-y') || '0')},flowVortexRadius:${JSON.stringify(argument('flow-vortex-radius') || '1')},
          flowVortexSpeed:${JSON.stringify(argument('flow-vortex-speed') || '.01')}}))document.getElementById(id).value=value;
        if(document.getElementById('flowVortexRadius').disabled)throw new Error('Initial vortex controls are disabled for fresh transient run');
      }
      if (${JSON.stringify(argument('flow-case-check')==='true')}) {
        const saved=await smoke.saveFlowCase();
        if(!saved)throw new Error('App did not save flow case');
        document.getElementById('flowNu').value='4';
        document.getElementById('flowTolerance').value='0.000001';
        document.getElementById('flowSpeed').value='2';
        document.getElementById('flowMode').value='steady';
        document.getElementById('flowMode').dispatchEvent(new Event('change'));
        document.getElementById('flowConvection').value='upwind';
        document.getElementById('flowSteadyAcceleration').value='none';
        document.getElementById('flowLinearPolicy').value='strict';
        document.getElementById('flowVelocityRelaxation').value='.4';
        document.getElementById('flowPressureCorrections').value='2';
        document.getElementById('flowInitialVortex').checked=false;
        document.getElementById('flowVortexSpeed').value='7';
        smoke.state.flowBoundaryDefinition=null;
        const loaded=await smoke.loadFlowCase();
        const repeated=await smoke.saveFlowCase();
        if(!loaded || !repeated || JSON.stringify(saved.document)!==JSON.stringify(repeated.document))
          throw new Error('Saved flow settings did not survive editor modification and reload');
        if(document.getElementById('flowResume').checked || document.getElementById('thermalResume').checked)
          throw new Error('Loading a zero-time case silently enabled resume');
        smoke.state.flowCaseSmoke={mesh:saved.document.mesh,request:saved.document.request,settingsRoundTrip:true};
      }
      await smoke.runFlow();
      if (${JSON.stringify(argument('flow-case-check')==='true')}) {
        if(!smoke.state.flow)throw new Error('Restored case failed to calculate');
        if(document.getElementById('flowMode').value!=='steady' && !document.getElementById('saveFlowCase').disabled)
          throw new Error('Resume settings could be saved as a fresh case');
        const fields=JSON.stringify(smoke.state.flow.fields),history=JSON.stringify(smoke.state.flow.history);
        if(!await smoke.loadFlowCase())throw new Error('Case reload after calculation failed');
        if(smoke.state.flow || !document.getElementById('flowResult').hidden)
          throw new Error('Old result was incorrectly displayed as loaded case output');
        await smoke.runFlow();
        if(!smoke.state.flow || fields!==JSON.stringify(smoke.state.flow.fields) || history!==JSON.stringify(smoke.state.flow.history))
          throw new Error('Repeated saved case changed its actual field or history');
        smoke.state.flowCaseSmoke.repeatedCalculationIdentical=true;
      }
      if (${JSON.stringify(argument('flow-vortex-speed')!==null)}) {
        if(!smoke.state.flow?.summary.initialVortex || !smoke.state.flow.files.initialCheckpoint)
          throw new Error('Initial vortex summary/checkpoint missing from real App result');
        if(!document.getElementById('flowInitialVortex').disabled || !document.getElementById('flowVortexRadius').disabled)
          throw new Error('Resume would allow duplicate vortex initialization');
        smoke.state.initialVortexSmoke=smoke.state.flow.summary.initialVortex;
      }
      if (${JSON.stringify(Boolean(argument('flow-end-time')))}) {
        if (!smoke.state.flow || smoke.state.flow.summary.timeStepControl!=='adaptive-cfl-retry') throw new Error('Missing adaptive App result');
        const initial=smoke.state.flow.summary;
        const before=initial.acceptedTime;
        const target=before+Number(document.getElementById('flowDt').value);
        document.getElementById('flowEndTime').value=String(target);
        document.getElementById('flowResume').checked=true;
        document.getElementById('flowResume').dispatchEvent(new Event('change'));
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==target) throw new Error('Adaptive absolute-time resume failed');
        const completed=smoke.state.flow;
        const iterations=document.getElementById('flowMaxIterations').value;
        document.getElementById('flowMaxIterations').value='1';
        document.getElementById('flowMaxRetries').value='0';
        document.getElementById('flowEndTime').value=String(target+.04);
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==target || smoke.state.flowRestart?.time!==target)
          throw new Error('Adaptive failure replaced the previous accepted result');
        document.getElementById('flowMaxIterations').value=iterations;
        document.getElementById('flowMaxRetries').value='10';
        document.getElementById('flowEndTime').value=String(target+100);
        const pending=smoke.runFlow(),deadline=Date.now()+90000;
        while(smoke.state.busy && !smoke.state.flowHistory.length && Date.now()<deadline)
          await new Promise(resolve=>setTimeout(resolve,25));
        if(!smoke.state.busy || !smoke.state.flowHistory.length) throw new Error('Adaptive cancellation had no accepted live step');
        document.getElementById('cancel').click();await pending;
        const saved=smoke.state.flowRestart?.time;
        if(!(saved>target) || smoke.state.flow?.summary.acceptedTime!==completed.summary.acceptedTime)
          throw new Error('Adaptive cancellation lost accepted checkpoint or earlier result');
        document.getElementById('flowEndTime').value=String(saved+.02);
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==saved+.02) throw new Error('Adaptive cancelled restart used wrong physical time');
        smoke.state.adaptiveSmoke={initial,failedAcceptedTime:target,cancelledAcceptedTime:saved,
          resumedTime:smoke.state.flow.summary.acceptedTime};
      }
      if (${JSON.stringify(Boolean(argument('flow-resume-steps')))}) {
        const before=smoke.state.flow?.summary.acceptedTime;
        document.getElementById('flowSteps').value=${JSON.stringify(argument('flow-resume-steps') || '2')};
        document.getElementById('flowResume').checked=true;
        document.getElementById('flowResume').dispatchEvent(new Event('change'));
        if (${JSON.stringify(Boolean(argument('flow-resume-pressure-preconditioner')))}) {
          const method=document.getElementById('flowPressurePreconditioner');
          if (method.disabled) throw new Error('Pressure solver selector is locked during idle resume');
          method.value=${JSON.stringify(argument('flow-resume-pressure-preconditioner') || 'ic0')};
          method.dispatchEvent(new Event('change'));
        }
        await smoke.runFlow();
        const expected=before+Number(document.getElementById('flowDt').value)*Number(document.getElementById('flowSteps').value);
        if (!smoke.state.flow || Math.abs(smoke.state.flow.summary.acceptedTime-expected)>1e-10)
          throw new Error('Desktop resume did not advance accepted physical time');
      }
      if (${JSON.stringify(argument('flow-failure-check') === 'true')}) {
        const accepted=smoke.state.flow.summary.acceptedTime;
        const maxIterations=document.getElementById('flowMaxIterations').value;
        document.getElementById('flowMaxIterations').value='1';
        document.getElementById('flowSteps').value='2';
        await smoke.runFlow();
        if (smoke.state.flow?.summary.acceptedTime!==accepted || smoke.state.flowRestart?.time!==accepted)
          throw new Error('Rejected time step replaced the complete result or accepted checkpoint');
        document.getElementById('flowMaxIterations').value=maxIterations;
      }
      if (${JSON.stringify(argument('flow-cancel-check') === 'true')}) {
        const completed=smoke.state.flow;
        document.getElementById('flowSteps').value='10000';
        const pendingFlow=smoke.runFlow();
        const deadline=Date.now()+90000;
        while (smoke.state.busy && !smoke.state.flowHistory.length && Date.now()<deadline)
          await new Promise(resolve=>setTimeout(resolve,50));
        if (!smoke.state.busy || !smoke.state.flowHistory.length) throw new Error('Could not observe live accepted step for cancellation');
        document.getElementById('cancel').click();
        await pendingFlow;
        if (!smoke.state.flowRestart || smoke.state.flowRestart.time<=completed.summary.acceptedTime
            || smoke.state.flow?.summary.acceptedTime!==completed.summary.acceptedTime)
          throw new Error('Cancellation lost checkpoint or previous complete result');
        const resumeTime=smoke.state.flowRestart.time;
        document.getElementById('flowSteps').value='2';
        await smoke.runFlow();
        if (Math.abs(smoke.state.flow.summary.acceptedTime-resumeTime-2*Number(document.getElementById('flowDt').value))>1e-10)
          throw new Error('Resume from cancelled calculation did not use the saved time');
      }
      if (!smoke.state.flow || document.getElementById('flowSpeedOption').hidden ||
          document.getElementById('displayMode').value !== 'speed')
        throw new Error('Native flow result did not reach the renderer');
    }
    if (${JSON.stringify(argument('euler') === 'true')}) {
      const target=Number(${JSON.stringify(argument('euler-end-time')||'.0005')});
      for(const [id,value] of Object.entries({eulerCase:${JSON.stringify(argument('euler-case')||'sod')},eulerDensity:'1.225',eulerPressure:'101325',
        eulerViscosity:${JSON.stringify(argument('euler-viscosity')||'0')},eulerWallModel:${JSON.stringify(argument('euler-wall-model')||'slip')},
        eulerConductivity:${JSON.stringify(argument('euler-conductivity')||'0')},eulerWallThermal:${JSON.stringify(argument('euler-wall-thermal')||'insulated')},eulerWallValue:${JSON.stringify(argument('euler-wall-value')||'0')},
        eulerU:${JSON.stringify(argument('euler-u')||'0')},eulerV:'0',eulerFluxScheme:${JSON.stringify(argument('euler-flux')||'rusanov')},eulerOrder:${JSON.stringify(argument('euler-order')||'1')},eulerGamma:'1.4',eulerGasConstant:'287.05',eulerEndTime:String(target)}))document.getElementById(id).value=value;
      document.getElementById('eulerCase').dispatchEvent(new Event('change'));
      await smoke.runEuler();
      if(!smoke.state.euler||smoke.state.euler.summary.time!==target)throw new Error('Euler result did not reach renderer: '+document.getElementById('statusText').textContent);
      const first=smoke.state.euler;
      if(first.summary.fluxScheme!==document.getElementById('eulerFluxScheme').value||first.summary.order!==Number(document.getElementById('eulerOrder').value))throw new Error('Euler method controls did not reach native solver');
      if(first.summary.thermalConductivity!==Number(document.getElementById('eulerConductivity').value)||first.summary.wallThermal!==document.getElementById('eulerWallThermal').value||first.summary.wallValue!==smoke.eulerRequest().wallValue)throw new Error('Thermal controls did not reach native solver');
      if(first.summary.dynamicViscosity!==Number(document.getElementById('eulerViscosity').value)||first.summary.wallModel!==document.getElementById('eulerWallModel').value)throw new Error('Viscous controls did not reach native solver');
      if(!['eulerViscosity','eulerWallModel','eulerConductivity','eulerWallThermal','eulerWallValue'].every(id=>document.getElementById(id).disabled))throw new Error('Restart thermal parameters are not locked');
      if(first.summary.dynamicViscosity>0) {
        document.getElementById('eulerResume').checked=false;document.getElementById('eulerResume').dispatchEvent(new Event('change'));
        document.getElementById('eulerViscosity').value=String(first.summary.dynamicViscosity*2);
        document.getElementById('eulerViscosity').dispatchEvent(new Event('change'));
        if(smoke.state.euler)throw new Error('Viscosity change left stale displayed fields');
        document.getElementById('eulerViscosity').value=String(first.summary.dynamicViscosity);
      }
      if(first.summary.thermalConductivity>0) {
        document.getElementById('eulerResume').checked=false;document.getElementById('eulerResume').dispatchEvent(new Event('change'));
        document.getElementById('eulerConductivity').value=String(first.summary.thermalConductivity*2);
        document.getElementById('eulerConductivity').dispatchEvent(new Event('change'));
        if(smoke.state.euler)throw new Error('Thermal change left stale displayed fields');
        document.getElementById('eulerConductivity').value=String(first.summary.thermalConductivity);
      }
      document.getElementById('eulerOrder').value=first.summary.order===1?'2':'1';
      document.getElementById('eulerOrder').dispatchEvent(new Event('change'));
      if(smoke.state.euler)throw new Error('Euler method change left a stale displayed result');
      document.getElementById('eulerOrder').value=String(first.summary.order);
      document.getElementById('eulerResume').checked=false;document.getElementById('eulerResume').dispatchEvent(new Event('change'));
      await smoke.runEuler();
      if(JSON.stringify(first.fields)!==JSON.stringify(smoke.state.euler?.fields)||JSON.stringify(first.history)!==JSON.stringify(smoke.state.euler?.history))throw new Error('Repeated Euler calculation differs: fields='+String(JSON.stringify(first.fields)===JSON.stringify(smoke.state.euler?.fields))+', history rows='+first.history.length+'/'+smoke.state.euler?.history.length);
      document.getElementById('eulerEndTime').value=String(2*target);await smoke.runEuler();
      if(smoke.state.euler?.summary.time!==2*target)throw new Error('Euler resume failed');
      const complete=smoke.state.euler;
      document.getElementById('eulerMaximumSteps').value='1';document.getElementById('eulerEndTime').value=String(10*target);await smoke.runEuler();
      if(smoke.state.euler.summary.time!==complete.summary.time||!(smoke.state.eulerRestart.time>complete.summary.time))throw new Error('Failed Euler budget lost prior complete result or accepted state');
      document.getElementById('eulerMaximumSteps').value='1000000';document.getElementById('eulerEndTime').value='100';
      const pending=smoke.runEuler(),deadline=Date.now()+30000;
      while(smoke.state.busy&&!smoke.state.eulerHistory.length&&Date.now()<deadline)await new Promise(r=>setTimeout(r,10));
      if(!smoke.state.busy||!smoke.state.eulerHistory.length)throw new Error('No live Euler progress observed');
      await window.cartmesh.cancel();await pending;
      const cancelledTime=smoke.state.eulerRestart?.time;
      if(!(cancelledTime>complete.summary.time)||smoke.state.euler.summary.time!==complete.summary.time)throw new Error('Euler cancel lost accepted state or previous result');
      document.getElementById('eulerMaximumSteps').value='100000';document.getElementById('eulerEndTime').value=String(cancelledTime+target);await smoke.runEuler();
      if(smoke.state.euler?.summary.time!==cancelledTime+target)throw new Error('Euler cancel/resume time mismatch');
      for(const mode of ['euler-rho','euler-p','euler-temperature','euler-mach','euler-speed']) {
        document.getElementById('displayMode').value=mode;document.getElementById('displayMode').dispatchEvent(new Event('change'));
        if(!smoke.view.fieldRange||!Number.isFinite(smoke.view.fieldRange.min))throw new Error('Euler field map missing: '+mode);
      }
      document.getElementById('displayMode').value=first.summary.thermalConductivity>0?'euler-temperature':'euler-rho';document.getElementById('displayMode').dispatchEvent(new Event('change'));
      document.getElementById('eulerBlock').scrollIntoView({block:'start'});
      smoke.state.eulerSmoke={dynamicViscosity:first.summary.dynamicViscosity,wallModel:first.summary.wallModel,viscousControlsReachedNative:true,viscousRestartLocked:true,viscosityChangeClearedStaleResult:first.summary.dynamicViscosity>0,thermalConductivity:first.summary.thermalConductivity,thermalControlsReachedNative:true,thermalRestartLocked:true,thermalChangeClearedStaleResult:first.summary.thermalConductivity>0,fluxScheme:first.summary.fluxScheme,order:first.summary.order,methodChangeClearedStaleResult:true,repeatedFieldsAndHistoryIdentical:true,resumeChecked:true,failedBudgetPreservedComplete:true,cancelledTime,cancelResumeChecked:true,allFiveFieldMaps:true};
    }
    if (${JSON.stringify(argument('thermal') === 'true')}) {
      for(const [id,value] of Object.entries({flowCase:'external',flowNu:'.1',flowSpeed:'1',flowConvection:${JSON.stringify(argument('flow-convection') || 'limited-linear')},flowPressurePreconditioner:'aggregation',flowMaxIterations:'1500',flowDt:'.05',flowSteps:'2',thermalDiffusivity:'.1'})) document.getElementById(id).value=value;
      await smoke.runThermal();
      if(!smoke.state.thermal || smoke.state.thermal.summary.time!==.1)throw new Error('Thermal result did not reach renderer');
      if(smoke.state.thermal.summary.outletBackflow!==${JSON.stringify(argument('outlet-backflow') || 'reject')})
        throw new Error('Thermal outlet mode did not reach the native solver');
      const before=smoke.state.thermal;
      await smoke.runThermal();
      if(!smoke.state.thermal || Math.abs(smoke.state.thermal.summary.time-.2)>1e-12)throw new Error('Thermal resume did not advance');
      const completed=smoke.state.thermal;
      document.getElementById('flowMaxIterations').value='1';
      await smoke.runThermal();
      if(smoke.state.thermal.summary.time!==completed.summary.time)throw new Error('Failed thermal run replaced complete result');
      document.getElementById('flowMaxIterations').value='1500';
      document.getElementById('flowSteps').value='10000';
      const pending=smoke.runThermal();
      const deadline=Date.now()+30000;
      while(smoke.state.busy&&!smoke.state.thermalHistory.length&&Date.now()<deadline)await new Promise(r=>setTimeout(r,30));
      if(!smoke.state.busy||!smoke.state.thermalHistory.length)throw new Error('No live thermal step observed before cancel');
      await window.cartmesh.cancel();await pending;
      const restartTime=smoke.state.thermalRestart?.time;
      if(!(restartTime>completed.summary.time)||smoke.state.thermal.summary.time!==completed.summary.time)throw new Error('Thermal cancellation lost accepted state');
      document.getElementById('flowSteps').value='1';await smoke.runThermal();
      if(Math.abs(smoke.state.thermal.summary.time-restartTime-.05)>1e-12)throw new Error('Thermal cancel/resume clock differs');
      if(document.getElementById('displayMode').value!=='temperature'||document.getElementById('thermalOption').hidden)throw new Error('Temperature map not displayed');
      document.getElementById('thermalBlock').scrollIntoView({block:'start'});
    }
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
      adaptive: smoke.state.adaptiveSmoke || null,
      initialVortex: smoke.state.initialVortexSmoke || null,
      flowCase: smoke.state.flowCaseSmoke || null,
      bundledChineseFontLoaded: true,
      theme: document.documentElement.dataset.theme,
      interactionChecks: ${JSON.stringify(Boolean(argument('interaction-check')))},
      euler:smoke.state.euler?{summary:smoke.state.euler.summary,request:smoke.state.euler.request,files:smoke.state.euler.files,manifest:smoke.state.euler.manifest,
        fieldCells:smoke.state.euler.fields.cells.length,historyRows:smoke.state.euler.history.length,audit:smoke.state.euler.audit,
        resultText:document.getElementById('eulerResult').innerText,restart:smoke.state.eulerRestart,checks:smoke.state.eulerSmoke,
        displayMode:document.getElementById('displayMode').value,monitorVisible:!document.getElementById('eulerTimeline').hidden}:null,
      thermal: smoke.state.thermal ? {summary:smoke.state.thermal.summary,request:smoke.state.thermal.request,files:smoke.state.thermal.files,
        fieldCells:smoke.state.thermal.fields.cells.length,historyRows:smoke.state.thermal.history.length,
        resultText:document.getElementById('thermalResult').innerText,restart:smoke.state.thermalRestart,
        displayMode:document.getElementById('displayMode').value,monitorVisible:!document.getElementById('thermalTimeline').hidden,
        failureAndCancelResumeChecked:true} : null,
      flow: smoke.state.flow ? {
        summary: smoke.state.flow.summary,
        fieldCells: smoke.state.flow.fields.cells.length,
        displayMode: document.getElementById('displayMode').value,
        resultText: document.getElementById('flowResult').innerText,
        historyRows: smoke.state.flow.history?.length || 0,
        monitorVisible: !document.getElementById('flowTimeline').hidden,
        restart: smoke.state.flowRestart,
        cancellationChecked: ${JSON.stringify(argument('flow-cancel-check') === 'true')},
        failureChecked: ${JSON.stringify(argument('flow-failure-check') === 'true')}
      } : null
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
      if (argument('flow-require-converged')==='true' && !report.flow?.summary?.converged)
      throw new Error('App acceptance requires a converged flow; diagnostic output is retained.');
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
    if(currentResult?.background) {
      report.background = await mainWindow.webContents.executeJavaScript(`(async()=>{
        const rejected=[];
        for(const name of ['runFlow','runThermal','runEuler']) {
          try {await window.cartmesh[name]({}); throw new Error('Background accepted by '+name);}
          catch(error){if(!error.message.includes('完整笛卡尔背景网格'))throw error;rejected.push(name);}
        }
        return {rejected,displayMode:document.getElementById('displayMode').value,
          panelsHidden:['flowBlock','thermalBlock','eulerBlock'].every(id=>document.getElementById(id).hidden),
          legend:document.getElementById('legend').innerText};
      })()`);
      if(!report.background.panelsHidden)throw new Error('Background grid exposed flow controls');
    }
    if (argument('export')) report.exported = await exportPackage(argument('export'));
    if (argument('flow-dt')) {
      await mainWindow.webContents.executeJavaScript("document.getElementById('flowTimeline').scrollIntoView({block:'nearest'}); document.getElementById('flowMode').scrollIntoView({block:'start'});");
      await mainWindow.webContents.capturePage(undefined, { stayAwake: true });
      await new Promise(resolve => setTimeout(resolve, 300));
    }
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
    if(report.euler && shot){
      mainWindow.setSize(1320,900);
      await new Promise(resolve=>setTimeout(resolve,200));
      await mainWindow.webContents.executeJavaScript("document.getElementById('eulerViscosity').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
      await new Promise(resolve=>setTimeout(resolve,200));
      await fs.writeFile(shot.replace(/\.png$/, '-euler-controls.png'),(await mainWindow.webContents.capturePage()).toPNG());
    }
    if(report.thermal && shot){
      mainWindow.setSize(1320,900);
      await new Promise(resolve=>setTimeout(resolve,200));
      await mainWindow.webContents.executeJavaScript("document.getElementById('thermalBlock').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=document.querySelector('.results').scrollHeight;");
      await new Promise(resolve=>setTimeout(resolve,200));
      await fs.writeFile(shot.replace(/\.png$/, '-thermal.png'),(await mainWindow.webContents.capturePage()).toPNG());
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
    if (argument('flow-dt')) {
      report.flow.smallWindowLayout = await mainWindow.webContents.executeJavaScript(`(() => {
        const a=document.getElementById('flowResult').getBoundingClientRect();
        const b=document.getElementById('flowTimeline').getBoundingClientRect();
        document.querySelector('.results').scrollTop=document.querySelector('.results').scrollHeight;
        return { resultBottom:a.bottom, monitorTop:b.top, separated:b.top>=a.bottom-1 };
      })()`);
      if (!report.flow.smallWindowLayout.separated) throw new Error('Transient monitor overlaps the result summary');
    }
    if (argument('flow-require-converged')==='true' && !report.flow?.summary?.converged)
      throw new Error('App acceptance requires a converged flow; diagnostic output is retained.');
    if (argument('flow-linear-policy') || argument('flow-velocity-relaxation') || argument('flow-pressure-corrections')) {
      report.flow.linearControls=await mainWindow.webContents.executeJavaScript(`(() => {
        const summary=window.__smoke.state.flow?.summary,result=document.getElementById('flowResult');
        const policy=document.getElementById('flowLinearPolicy').value;
        const relaxation=Number(document.getElementById('flowVelocityRelaxation').value);
        const pressureCorrections=Number(document.getElementById('flowPressureCorrections').value);
        if(!summary || summary.linearPolicy!==policy || summary.velocityRelaxation!==relaxation || summary.pressureCorrectionPasses!==pressureCorrections)
          throw new Error('Linear efficiency controls do not match actual result');
        if(policy==='adaptive' && (!summary.strictLinearFinal || !result.innerText.includes('严格复核通过')))
          throw new Error('Adaptive linear solve lacks strict final certification in App');
        if(!result.innerText.includes('速度松弛系数'))throw new Error('Relaxation not shown in actual App');
        if(!result.innerText.includes('压力校正次数'))throw new Error('Pressure corrections not shown in App');
        return {policy,relaxation,pressureCorrections,strictLinearFinal:summary.strictLinearFinal,strictAcceptedSteps:summary.strictAcceptedSteps};
      })()`);
    }
    console.log(JSON.stringify(report, null, 2));
    if (!report.layout.bottomReachable) throw new Error('Sidebar bottom is inaccessible');
    if (Math.abs(report.layout.geometry.bodyHeight - report.layout.geometry.innerHeight) > 1 ||
        report.layout.geometry.canvasHeight <= 0 || report.layout.geometry.shellBottom > report.layout.geometry.innerHeight + 1)
      throw new Error('Window geometry is not constrained to the viewport');
    if (shot) {
      await mainWindow.webContents.executeJavaScript("document.querySelector('.panel').scrollTop=0");
      if (argument('flow-adaptive-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowMode').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
        await new Promise(resolve=>setTimeout(resolve,200));
      }
      if (argument('flow-initial-shot')) {
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowInitialSettings').open=true; document.getElementById('flowInitialSettings').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
      }
      if (argument('flow-case-check')) {
        await mainWindow.webContents.executeJavaScript("document.getElementById('saveFlowCase').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
      }
      if (argument('flow-load-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult');
          const loads=window.__smoke.state.flow?.summary.namedWallLoads;
          if (!loads?.length || !result.innerText.includes('力矩')) throw new Error('Wall loads did not reach the actual App');
          for (const load of loads) if (!result.innerText.includes(load.name)) throw new Error('Missing named wall load in App');
          result.scrollIntoView({block:'end'});
        })()`);
      }
      if (argument('flow-flux-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult'), summary=window.__smoke.state.flow?.summary;
          if(!summary?.namedBoundaryFluxes?.length || !result.innerText.includes('净流量'))throw new Error('Boundary flux did not reach App');
          for(const flux of summary.namedBoundaryFluxes)if(!result.innerText.includes(flux.name+' · 净流量'))throw new Error('Missing named boundary flux in App');
          [...result.querySelectorAll('div')].find(row=>row.innerText.includes('边界体积流量'))?.scrollIntoView({block:'start'});
        })()`);
      }
      if (argument('flow-boundary-shot')) {
        mainWindow.setSize(1100,800);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowBoundarySettings').scrollIntoView({block:'start'});");
      }
      if (argument('flow-acceleration-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult'),summary=window.__smoke.state.flow?.summary;
          if(summary?.steadyAcceleration!=='anderson' || !(summary.accelerationAccepted>0) || !result.innerText.includes('历史迭代'))
            throw new Error('Steady acceleration did not reach the actual App result');
          if(document.getElementById('flowSteadyAcceleration').value!=='anderson')throw new Error('Acceleration control differs from result');
          document.getElementById('flowSteadyAcceleration').scrollIntoView({block:'center'});
          result.scrollIntoView({block:'start'});
        })()`);
      }
      if (argument('flow-linear-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowLinearPolicy').scrollIntoView({block:'center'}); document.getElementById('flowResult').scrollIntoView({block:'start'});");
      }
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
