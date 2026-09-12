'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const yauzl = require('yauzl');
const { zipDirectory } = require('../src/core/archive');

function readZip(file) {
  return new Promise((resolve, reject) => {
    yauzl.open(file, { lazyEntries: true }, (error, zip) => {
      if (error) return reject(error);
      const result = {};
      zip.on('error', reject);
      zip.on('end', () => resolve(result));
      zip.on('entry', entry => {
        if (entry.fileName.endsWith('/')) { result[entry.fileName] = null; zip.readEntry(); return; }
        zip.openReadStream(entry, (err, stream) => {
          if (err) return reject(err);
          const chunks = [];
          stream.on('error', reject);
          stream.on('data', chunk => chunks.push(chunk));
          stream.on('end', () => { result[entry.fileName] = Buffer.concat(chunks); zip.readEntry(); });
        });
      });
      zip.readEntry();
    });
  });
}
test('ZIP retains UTF-8 names, binary bytes and empty case directories', async t => {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'mesh archive '));
  t.after(() => fs.rm(root, { recursive: true, force: true }));
  const source = path.join(root, '结果 包');
  await fs.mkdir(path.join(source, 'constant', 'polyMesh'), { recursive: true });
  const bytes = Buffer.from([0, 255, 17, 0, 10]);
  await fs.writeFile(path.join(source, '原图.png'), bytes);
  await fs.writeFile(path.join(source, 'README_CN.md'), '网格说明\n');
  const target = path.join(root, '网格.zip');
  await zipDirectory(source, target);
  const entries = await readZip(target);
  assert.deepEqual(entries['结果 包/原图.png'], bytes);
  assert.equal(entries['结果 包/README_CN.md'].toString(), '网格说明\n');
  assert.equal(entries['结果 包/constant/polyMesh/'], null);
  assert.ok(Object.keys(entries).every(name => !name.includes('\\')));
  await assert.rejects(zipDirectory(source, target), /EEXIST/);
  assert.deepEqual((await readZip(target))['结果 包/原图.png'], bytes);
});
test('ZIP refuses recursive output and cancelled work without damaging input', async t => {
  const source = await fs.mkdtemp(path.join(os.tmpdir(), 'zip-guard-'));
  t.after(() => fs.rm(source, { recursive: true, force: true }));
  await assert.rejects(zipDirectory(source, path.join(source, 'export.zip')), /不能写入/);
  const controller = new AbortController(); controller.abort();
  await assert.rejects(zipDirectory(source, source + '.zip', controller.signal), /abort/i);
});
