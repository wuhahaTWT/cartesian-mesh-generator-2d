'use strict';

const fs = require('node:fs/promises');
const { createReadStream } = require('node:fs');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { extractZip } = require('./archive');
const { parseCm2d, assignSizeBands, levelHistogram, embeddedBounds } = require('./cm2d');
const { parseBackgroundGrid } = require('./background-grid');
const { METHODS } = require('./capabilities');
const geometry = require('./geometry');
const { readAcceptedCheckpointMetadata } = require('./flow-checkpoint');
const { validateFlowOutput, validateTimeHistory, validateAttemptHistory } = require('./flow');
const { thermalCheckpointTime, validateThermalRequest, validateThermalOutput } = require('./thermal');
const { validateEulerOutput } = require('./euler');
const {checkpointDigest,readThermalTimeline}=require('./thermal-history');
const { readCheckpoint: readEulerCheckpoint } = require('./euler-job');

const MANIFEST = 'cartmesh2d-project.json';
const FORMAT = 'cartmesh2d-project-v1';
const fail = message => { throw new Error(`项目：${message}`); };
function localPath(root, name) {
  if (typeof name !== 'string' || !name || /[\\\0:]/.test(name) || path.posix.isAbsolute(name) ||
      name.split('/').some(part => part === '..' || part === '')) fail('文件路径越出项目目录。');
  return path.join(root, ...name.split('/'));
}
function relative(root, file) {
  const name = path.relative(root, file).split(path.sep).join('/') || '.';
  localPath(root, name);
  return name;
}
function transform(value, visit) {
  if (Array.isArray(value)) return value.map(item => transform(item, visit));
  if (value && typeof value === 'object')
    return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, transform(item, visit)]));
  return visit(value);
}
async function digest(file, signal) {
  const hash = createHash('sha256');
  let bytes = 0;
  for await (const chunk of createReadStream(file, { signal })) { hash.update(chunk); bytes += chunk.length; }
  return { bytes, sha256: hash.digest('hex') };
}
async function inventory(root, signal) {
  const records = [];
  async function visit(directory) {
    for (const entry of (await fs.readdir(directory, { withFileTypes: true })).sort((a,b) => a.name.localeCompare(b.name))) {
      signal?.throwIfAborted();
      const file = path.join(directory, entry.name), name = relative(root, file);
      if (name === MANIFEST || name === MANIFEST + '.tmp') continue;
      if (entry.isDirectory()) await visit(file);
      else if (entry.isFile()) records.push({ name, ...await digest(file, signal) });
      else fail(`不支持链接或特殊文件：${name}`);
    }
  }
  await visit(root);
  return records;
}

// The existing result ZIP becomes a self-contained project. Keep native files,
// including failed attempts, and store only bindings/UI here; fields are rebuilt
// from those native files on open. No accepted output is rewritten.
async function writeProjectManifest(current, ui, signal) {
  if (!current?.outputDirectory || current.incomplete) fail('需要完整生成的网格。');
  const root = path.resolve(current.outputDirectory);
  const state = structuredClone({ ...current, mesh: undefined, levelHistogram: undefined, projectUi: undefined });
  state.outputDirectory = root;
  state.job.outputDirectory = root;
  state.job.geometryPath = path.resolve(current.prefix + '.xy');
  state.job.sourceUnits = 'm';
  if (state.selectedRequest) {
    state.selectedRequest.geometryPath = state.job.geometryPath;
    state.selectedRequest.sourceUnits = 'm';
    state.selectedRequest.outputDirectory = root;
  }
  await fs.access(state.job.geometryPath);
  for (const kind of ['flow', 'thermal', 'euler']) {
    const restart = state[kind + 'Restart'];
    if (restart) {
      const original = path.resolve(restart.path);
      restart.path = original;
      if (!(original.startsWith(root + path.sep))) {
        const directory = path.join(root, 'project-inputs');
        await fs.mkdir(directory, { recursive: true });
        restart.path = path.join(directory, kind + '.checkpoint');
        await fs.copyFile(original, restart.path);
      }
    }
    if (state[kind]) {
      const { request, files, manifest } = state[kind];
      const portableFiles = Object.fromEntries(Object.entries(files).map(([key, file]) => [key, file.split(path.sep).join('/')]));
      state[kind] = { request, files: portableFiles, ...(manifest ? { manifest: manifest.split(path.sep).join('/') } : {}) };
    }
  }
  // Native output names retain the spelling of the requested output directory,
  // which may be relative. Normalize that known prefix before making bindings
  // portable; otherwise the ZIP silently keeps paths into the old checkout.
  const sourceRoot = path.normalize(current.outputDirectory);
  const portable = transform(state, value => {
    if (typeof value !== 'string') return value;
    let file = value;
    if (!path.isAbsolute(value) && sourceRoot !== '.' &&
        (path.normalize(value) === sourceRoot || path.normalize(value).startsWith(sourceRoot + path.sep)))
      file = path.resolve(value);
    return file === root || file.startsWith(root + path.sep) ? 'project://' + relative(root, file) : value;
  });
  const document = { format: FORMAT, geometryLabel: path.basename(current.job.geometryPath),
    state: portable, ui: ui ?? current.projectUi ?? null, files: await inventory(root, signal) };
  signal?.throwIfAborted();
  const temporary = path.join(root, MANIFEST + '.tmp');
  await fs.writeFile(temporary, JSON.stringify(document, null, 2) + '\n');
  await fs.rename(temporary, path.join(root, MANIFEST));
  return document;
}

async function readProject(root, signal) {
  const document = JSON.parse(await fs.readFile(path.join(root, MANIFEST), 'utf8'));
  if (document?.format !== FORMAT || !document.state || !Array.isArray(document.files) || !document.files.length)
    fail('项目版本或文件清单无效。');
  const names = new Set();
  for (const item of document.files) {
    const file = localPath(root, item.name);
    if (names.has(item.name) || item.name === MANIFEST) fail('文件清单包含重复项。');
    names.add(item.name);
    const actual = await digest(file, signal);
    if (item.bytes !== actual.bytes || item.sha256 !== actual.sha256) fail(`文件损坏或被更改：${item.name}`);
  }
  const state = transform(document.state, value => typeof value === 'string' && value.startsWith('project://')
    ? localPath(root, value.slice('project://'.length)) : value);
  const checked = file => {
    const name = relative(root, file);
    if (!names.has(name)) fail(`引用的文件不在清单中：${name}`);
    return file;
  };
  const text = file => fs.readFile(checked(file), 'utf8');
  const native = (binding, key) => {
    if (!binding?.files?.[key]) fail(`缺少原生结果文件：${key}`);
    return localPath(root, binding.files[key]);
  };
  if (state.outputDirectory !== root || state.job?.outputDirectory !== root || !state.result) fail('项目绑定无效。');
  if (!Object.hasOwn(METHODS, state.job.method) || !['interior','exterior'].includes(state.job.fluidRegion)) fail('网格方法或流体域语义无效。');
  const converted = geometry.convertToLoops(state.job.geometryPath, await text(state.job.geometryPath), {sourceUnits:'m'});
  if (converted.issues.length) fail(converted.issues.join('；'));
  relative(root, state.prefix);
  const mesh = state.background
    ? parseBackgroundGrid(await text(state.backgroundPath)) : assignSizeBands(parseCm2d(await text(state.cm2dPath)));
  if (state.result.counts?.cells !== mesh.cells.length) fail('网格与保存的结果数量不一致。');
  state.mesh = mesh;
  state.levelHistogram = levelHistogram(mesh);
  state.wallBounds = state.background ? null : embeddedBounds(mesh);
  if (state.background && ['flow','thermal','euler','flowRestart','thermalRestart','eulerRestart'].some(key => state[key]))
    fail('完整笛卡尔背景不能载入流体状态。');
  if (state.flow) {
    const saved = state.flow;
    const summary = JSON.parse(await text(native(saved, 'summary')));
    const fields = JSON.parse(await text(native(saved, 'fields')));
    const transient = saved.request.mode !== 'steady';
    const historyText = transient ? await text(native(saved, 'timeHistory')) : null;
    // Each CSV row carries the actual dt, including adaptive first steps.
    const first = historyText?.trim().split(/\r?\n/)[1]?.split(',').map(Number);
    const start = first ? first[1] - first[2] : 0;
    const validated = validateFlowOutput(summary, fields, mesh.cells.length, saved.request, start);
    const history = transient ? validateTimeHistory(historyText, validated.summary, start) : null;
    const attempts = saved.request.mode === 'adaptive'
      ? validateAttemptHistory(await text(native(saved, 'attemptHistory')), validated.summary, history) : null;
    state.flow = { ...validated, request: saved.request, files: saved.files, history, attempts };
  }
  if (state.thermal) {
    const saved = state.thermal;
    const summary = JSON.parse(await text(native(saved, '.json')));
    const validated = validateThermalOutput(summary, await text(native(saved, '.cells.csv')),
      await text(native(saved, '.thermal-history.csv')), await text(native(saved, '.thermal.checkpoint')),
      mesh, saved.request, summary.startTime ?? (summary.time - summary.steps * summary.timeStep));
    state.thermal = { ...validated, files: saved.files };
  }
  if (state.euler) {
    const saved = state.euler;
    const summary = JSON.parse(await text(native(saved, '.json')));
    const fields = JSON.parse(await text(native(saved, '.fields.json')));
    const files = await Promise.all(['.cells.csv','.faces.csv','.history.csv','.checkpoint'].map(key => text(native(saved, key))));
    const validated = validateEulerOutput(summary, fields, ...files, mesh, saved.request, summary.initialTime);
    state.euler = { ...validated, files: saved.files, manifest: saved.manifest };
  }
  if (state.flowRestart) {
    const saved = state.flowRestart, metadata = await readAcceptedCheckpointMetadata(checked(saved.path));
    if (metadata.time !== saved.metadata.time) fail('流动检查点时钟与保存记录不同。');
    state.flowRestart = { path: saved.path, metadata };
  }
  if (state.thermalRestart) {
    const saved = state.thermalRestart, request = validateThermalRequest(saved.metadata.request);
    const checkpoint=await text(saved.path),time=thermalCheckpointTime(checkpoint,request.events);
    if(saved.sha256 && saved.sha256!==checkpointDigest(checkpoint))fail('联合状态与保存历史的身份不同。');
    if (time !== saved.metadata.time) fail('联合检查点时钟与保存记录不同。');
    state.thermalRestart = {...saved,metadata:{time,request}};
    await readThermalTimeline(state,text);
  }
  if (state.eulerRestart) {
    const saved = state.eulerRestart;
    state.eulerRestart = await readEulerCheckpoint(checked(saved.path), mesh, saved.metadata.request);
    if (state.eulerRestart.metadata.time !== saved.metadata.time) fail('可压检查点时钟与保存记录不同。');
  }
  signal?.throwIfAborted();
  state.projectUi = document.ui;
  state.projectLabel = document.geometryLabel;
  return state;
}

// Extract into a new owned directory. The caller publishes only after validation,
// so a failed/cancelled open leaves the existing project and source ZIP untouched.
async function openProject(file, sessionDirectory, signal) {
  const scratch = await fs.mkdtemp(path.join(sessionDirectory, 'project-'));
  try {
    await extractZip(file, scratch, signal);
    const entries = await fs.readdir(scratch, { withFileTypes: true });
    if (entries.length !== 1 || !entries[0].isDirectory()) fail('项目包必须包含一个项目目录。');
    const root = path.join(scratch, entries[0].name);
    try { await fs.access(path.join(root, MANIFEST)); }
    catch { fail('此结果包没有项目清单；旧结果仍可按原流程载入检查点。'); }
    return await readProject(root, signal);
  } catch (error) {
    await fs.rm(scratch, { recursive: true, force: true });
    throw error;
  }
}

module.exports = { MANIFEST, writeProjectManifest, readProject, openProject };
