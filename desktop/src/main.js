const { app, BrowserWindow, dialog, ipcMain, shell } = require('electron');
const { spawn } = require('node:child_process');
const fs = require('node:fs/promises');
const path = require('node:path');
const {
  getCapabilities,
  estimateMeshJob,
  validateMeshJob,
  buildMeshInvocation,
  normalizeSummary
} = require('./mesh-tools');

let mainWindow;

function resourcePath(...parts) {
  const root = app.isPackaged ? process.resourcesPath : path.join(__dirname, '..', 'runtime');
  return path.join(root, ...parts);
}

function executable(name) {
  return resourcePath('bin', process.platform === 'win32' ? `${name}.exe` : name);
}

function runProcess(command, args, onLine) {
  return new Promise((resolve, reject) => {
    const child = spawn(command, args, { windowsHide: true });
    let stdout = '';
    let stderr = '';
    const consume = (chunk, target) => {
      const text = chunk.toString();
      if (target === 'out') stdout += text;
      else stderr += text;
      text.split(/\r?\n/).filter(Boolean).forEach(onLine);
    };
    child.stdout.on('data', chunk => consume(chunk, 'out'));
    child.stderr.on('data', chunk => consume(chunk, 'err'));
    child.on('error', reject);
    child.on('close', code => {
      if (code === 0) resolve({ stdout, stderr });
      else reject(new Error(stderr.trim() || stdout.trim() || `Process exited with code ${code}`));
    });
  });
}

function parseSummary(output) {
  const summary = {};
  for (const line of output.split(/\r?\n/)) {
    const match = line.match(/^([a-z0-9_]+)=(.+)$/i);
    if (match) summary[match[1]] = match[2];
  }
  return summary;
}

function safeBaseName(filePath) {
  return path.basename(filePath, path.extname(filePath))
    .replace(/[^a-zA-Z0-9_-]+/g, '_') || 'mesh';
}

async function uniquePrefix(directory, base) {
  const first = path.join(directory, base);
  const occupiedMarkers = [
    `${first}.job.json`, `${first}.cm2d`, `${first}.hybrid.cm2d`, `${first}.fallback.cm2d`
  ];
  for (const marker of occupiedMarkers) {
    try {
      await fs.access(marker);
      const stamp = new Date().toISOString().replace(/[:.]/g, '-');
      return path.join(directory, `${base}-${stamp}`);
    } catch {
      // This marker is free; all markers must be checked before reusing the prefix.
    }
  }
  return first;
}

async function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1320,
    height: 820,
    minWidth: 1040,
    minHeight: 680,
    backgroundColor: '#eef2f3',
    titleBarStyle: 'hiddenInset',
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false
    }
  });
  await mainWindow.loadFile(path.join(__dirname, 'renderer', 'index.html'));
}

app.whenReady().then(() => {
  ipcMain.handle('select-dxf', async () => {
    const result = await dialog.showOpenDialog(mainWindow, {
      title: '选择二维 DXF 边界',
      filters: [{ name: 'DXF 文件', extensions: ['dxf'] }],
      properties: ['openFile']
    });
    return result.canceled ? null : result.filePaths[0];
  });

  ipcMain.handle('select-output', async () => {
    const result = await dialog.showOpenDialog(mainWindow, {
      title: '选择输出目录',
      properties: ['openDirectory', 'createDirectory']
    });
    return result.canceled ? null : result.filePaths[0];
  });

  ipcMain.handle('example-path', () => resourcePath('examples', 'spline_circle_mm.dxf'));

  ipcMain.handle('smoke-config', () => {
    if (!app.commandLine.hasSwitch('smoke-test')) return { enabled: false };
    const requestedMethod = app.commandLine.getSwitchValue('smoke-method');
    const method = requestedMethod === 'hybrid' ? 'hybrid' : 'cutcell';
    return {
      enabled: true,
      method,
      outputDirectory: `/private/tmp/cartmesh2d-product-ui-smoke-${method}`
    };
  });

  ipcMain.handle('smoke-capture', async () => {
    if (!app.commandLine.hasSwitch('smoke-test')) return null;
    const image = await mainWindow.webContents.capturePage();
    const requestedMethod = app.commandLine.getSwitchValue('smoke-method');
    const suffix = requestedMethod === 'hybrid' ? '-hybrid' : '-cutcell';
    const target = `/private/tmp/cartmesh2d-product-app${suffix}.png`;
    await fs.writeFile(target, image.toPNG());
    return target;
  });

  ipcMain.handle('open-path', async (_event, target) => shell.openPath(target));

  ipcMain.handle('mesh-capabilities', () => getCapabilities());
  ipcMain.handle('estimate-mesh-job', (_event, request) => estimateMeshJob(request));

  ipcMain.handle('generate-mesh', async (_event, request) => {
    const { job, method, estimate } = validateMeshJob(request);
    await fs.mkdir(job.outputDirectory, { recursive: true });
    const base = `${safeBaseName(job.dxfPath)}-${method.id}`;
    const prefix = await uniquePrefix(job.outputDirectory, base);
    const xyPath = `${prefix}.xy`;
    const dxfReport = `${prefix}.dxf.json`;
    const casePath = `${prefix}-openfoam`;
    const jobPath = `${prefix}.job.json`;
    const send = line => mainWindow?.webContents.send('generation-line', line);

    await fs.writeFile(jobPath, `${JSON.stringify({
      schemaVersion: getCapabilities().schemaVersion,
      job,
      estimate
    }, null, 2)}\n`);
    const dxfArgs = [job.dxfPath, xyPath, String(job.chordError), dxfReport];
    if (job.sourceUnits && job.sourceUnits !== 'auto') {
      dxfArgs.push('1e-10', job.sourceUnits);
    }
    send('正在读取 DXF、换算单位并离散曲线…');
    const dxfRun = await runProcess(executable('cartmesh2d_dxf_cli'), dxfArgs, send);

    const invocation = buildMeshInvocation(job, { xyPath, prefix, casePath });
    send(invocation.progress);
    const meshRun = await runProcess(
      executable(invocation.executableName), invocation.args, send);
    const summary = normalizeSummary(job, parseSummary(meshRun.stdout));
    let cm2dPath = null;
    for (const candidate of invocation.cm2dCandidates) {
      try {
        await fs.access(candidate);
        cm2dPath = candidate;
        break;
      } catch {
        // Try the next explicit output produced by the selected tool.
      }
    }
    if (!cm2dPath) throw new Error('生成器成功退出，但没有找到可读取的 solver CM2D 输出。');
    const cm2d = await fs.readFile(cm2dPath, 'utf8');
    if (summary.actual_method === 'cutcell-fallback') {
      send('Hybrid 未形成，已明确退化为 Pure Cut-cell fallback。');
    } else if (job.method === 'hybrid' && !summary.quality_pass) {
      send('Hybrid 网格已生成供检查，但严格质量合同未通过，Beta 工具未输出 OpenFOAM case。');
    } else {
      send('生成完成。');
    }
    return {
      prefix,
      outputDirectory: job.outputDirectory,
      cm2dPath,
      dxfReport,
      casePath,
      jobPath,
      method,
      estimate,
      summary,
      cm2d,
      converterOutput: dxfRun.stdout
    };
  });

  createWindow();
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});
