const { app, BrowserWindow, dialog, ipcMain, shell } = require('electron');
const { spawn } = require('node:child_process');
const fs = require('node:fs/promises');
const path = require('node:path');

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
  try {
    await fs.access(`${first}.cm2d`);
  } catch {
    return first;
  }
  const stamp = new Date().toISOString().replace(/[:.]/g, '-');
  return path.join(directory, `${base}-${stamp}`);
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

  ipcMain.handle('smoke-config', () => app.commandLine.hasSwitch('smoke-test') ? {
    enabled: true,
    outputDirectory: '/private/tmp/cartmesh2d-product-ui-smoke'
  } : { enabled: false });

  ipcMain.handle('smoke-capture', async () => {
    if (!app.commandLine.hasSwitch('smoke-test')) return null;
    const image = await mainWindow.webContents.capturePage();
    const target = '/private/tmp/cartmesh2d-product-app.png';
    await fs.writeFile(target, image.toPNG());
    return target;
  });

  ipcMain.handle('open-path', async (_event, target) => shell.openPath(target));

  ipcMain.handle('generate-mesh', async (_event, request) => {
    const { dxfPath, outputDirectory, chordError, sourceUnits, maxLevel,
      minimumLevel, paddingFraction, smallAlpha } = request;
    if (!dxfPath || !outputDirectory) throw new Error('请选择 DXF 文件和输出目录。');
    await fs.mkdir(outputDirectory, { recursive: true });
    const base = safeBaseName(dxfPath);
    const prefix = await uniquePrefix(outputDirectory, base);
    const xyPath = `${prefix}.xy`;
    const dxfReport = `${prefix}.dxf.json`;
    const casePath = `${prefix}-openfoam`;
    const send = line => mainWindow?.webContents.send('generation-line', line);

    const dxfArgs = [dxfPath, xyPath, String(chordError), dxfReport];
    if (sourceUnits && sourceUnits !== 'auto') dxfArgs.push('1e-10', sourceUnits);
    send('正在读取 DXF、换算单位并离散曲线…');
    const dxfRun = await runProcess(executable('cartmesh2d_dxf_cli'), dxfArgs, send);

    send('正在生成 Cartesian / Cut-cell 网格…');
    const meshArgs = [xyPath, prefix, String(maxLevel), String(paddingFraction),
      String(smallAlpha), 'exterior', casePath, String(minimumLevel), '0'];
    const meshRun = await runProcess(executable('cartmesh2d_cli'), meshArgs, send);
    const summary = parseSummary(meshRun.stdout);
    const cm2dPath = `${prefix}.cm2d`;
    const cm2d = await fs.readFile(cm2dPath, 'utf8');
    send('生成完成。');
    return {
      prefix,
      outputDirectory,
      cm2dPath,
      dxfReport,
      casePath,
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
