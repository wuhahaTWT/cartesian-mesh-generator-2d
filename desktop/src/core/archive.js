'use strict';
const fs = require('node:fs');
const fsp = require('node:fs/promises');
const path = require('node:path');
const { pipeline } = require('node:stream/promises');
const { ZipFile } = require('yazl');

// Stream file contents, retain empty OpenFOAM directories, and use UTF-8 ZIP names
// on every OS. Never invoke a shell or depend on a system zip utility.
async function zipDirectory(source, destination, signal) {
  const root = path.resolve(source), target = path.resolve(destination);
  if (target === root || target.startsWith(root + path.sep))
    throw new Error('结果包不能写入正在打包的目录。');
  signal?.throwIfAborted();
  const entries = [];
  async function visit(directory, relative = path.basename(root)) {
    const children = await fsp.readdir(directory, { withFileTypes: true });
    entries.push({ relative, directory: true });
    for (const child of children.sort((a, b) => a.name < b.name ? -1 : a.name > b.name ? 1 : 0)) {
      signal?.throwIfAborted();
      const file = path.join(directory, child.name), name = relative + '/' + child.name;
      if (child.isDirectory()) await visit(file, name);
      else if (child.isFile()) entries.push({ file, relative: name });
      else throw new Error(`结果包包含不支持的链接或特殊文件：${name}`);
    }
  }
  await visit(root);
  const zip = new ZipFile();
  const output = fs.createWriteStream(target, { flags: 'wx' });
  zip.on('error', error => zip.outputStream.destroy(error));
  const writing = pipeline(zip.outputStream, output, { signal });
  try {
    for (const entry of entries) {
      if (entry.directory) zip.addEmptyDirectory(entry.relative);
      else zip.addFile(entry.file, entry.relative);
    }
    zip.end();
    await writing;
  } catch (error) {
    zip.outputStream.destroy(error);
    await writing.catch(() => {});
    // Do not remove an existing destination rejected by flags:'wx'.
    if (error.code !== 'EEXIST') await fsp.rm(target, { force: true });
    throw error;
  }
}

// Replace the user's saved ZIP only after a complete archive exists on the
// same filesystem. Cancellation/errors cannot truncate an earlier project.
async function saveDirectoryArchive(source, destination, signal) {
  const root = path.resolve(source), target = path.resolve(destination);
  if (target === root || target.startsWith(root + path.sep))
    throw new Error('结果包不能写入正在打包的目录。');
  signal?.throwIfAborted();
  const staging = await fsp.mkdtemp(path.join(path.dirname(target), '.cartmesh2d-save-'));
  try {
    const file = path.join(staging, 'project.zip');
    await zipDirectory(root, file, signal);
    signal?.throwIfAborted();
    await fsp.rename(file, target);
  } finally { await fsp.rm(staging, { recursive: true, force: true }); }
}

// Read only into a fresh owned directory. Lazy, streamed extraction avoids
// loading meshes into memory and rejects links, duplicate paths and traversal.
async function extractZip(source, destination, signal) {
  signal?.throwIfAborted();
  const root = path.resolve(destination), names = new Set();
  const zip = await new Promise((resolve, reject) => require('yauzl').open(source,
    { lazyEntries: true, strictFileNames: true, validateEntrySizes: true },
    (error, value) => error ? reject(error) : resolve(value)));
  await new Promise((resolve, reject) => {
    let active = Promise.resolve(), stopped = false;
    const abort = () => fail(signal.reason || new Error('操作已取消'));
    const fail = error => {
      stopped = true; zip.close(); signal?.removeEventListener('abort', abort);
      // Wait for an in-flight pipeline to close its files before caller cleanup.
      active.catch(() => {}).then(() => reject(error));
    };
    signal?.addEventListener('abort', abort, { once: true });
    zip.on('error', fail);
    zip.on('end', () => { signal?.removeEventListener('abort', abort); resolve(); });
    zip.on('entry', entry => {
      active = (async () => {
        signal?.throwIfAborted();
        const directory = entry.fileName.endsWith('/');
        const name = directory ? entry.fileName.slice(0,-1) : entry.fileName;
        const mode = (entry.externalFileAttributes >>> 16) & 0xf000;
        if (!name || /[\\\0:]/.test(name) || path.posix.isAbsolute(name) ||
            name.split('/').some(part => !part || part === '.' || part === '..') ||
            (mode && mode !== (directory ? 0x4000 : 0x8000)) || names.has(name))
          throw new Error('项目包包含无效路径、重复文件或链接。');
        names.add(name);
        const target = path.join(root, ...name.split('/'));
        if (directory) await fsp.mkdir(target, { recursive: true });
        else {
          await fsp.mkdir(path.dirname(target), { recursive: true });
          const input = await new Promise((done, error) => zip.openReadStream(entry,
            (err, stream) => err ? error(err) : done(stream)));
          await pipeline(input, fs.createWriteStream(target, { flags: 'wx' }), { signal });
        }
      })();
      active.then(() => { if (!stopped) zip.readEntry(); }, fail);
    });
    if (signal?.aborted) abort(); else zip.readEntry();
  });
}

module.exports = { zipDirectory, extractZip, saveDirectoryArchive };
