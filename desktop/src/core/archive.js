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
module.exports = { zipDirectory };
