'use strict';

const fs = require('node:fs');
const path = require('node:path');
const { version } = require('../package.json');

const args = process.argv.slice(2);
if (args.length && (args.length !== 1 || args[0] !== '--apply')) {
  throw new Error('Usage: node scripts/prune-packages.js [--apply]');
}
const dist = path.join(__dirname, '..', 'dist');
const files = fs.existsSync(dist)
  ? fs.readdirSync(dist, { withFileTypes: true }).filter(entry => entry.isFile()).map(entry => entry.name)
  : [];
const groups = new Map();
for (const name of files) {
  const match = name.match(/^CartMesh2D-(\d+\.\d+\.\d+)(?:-(arm64|x64|ia32))?-(mac|win)\.zip$/)
    || name.match(/^cartmesh2d-desktop-(\d+\.\d+\.\d+)(?:-(arm64|x64|ia32))?\.tar\.gz$/);
  if (!match) continue;
  const key = `${match[3] || 'linux'}-${match[2] || 'x64'}`;
  if (!groups.has(key)) groups.set(key, []);
  groups.get(key).push({ name, version: match[1] });
}
function newestFirst(a, b) {
  const x = a.version.split('.').map(Number), y = b.version.split('.').map(Number);
  return y[0] - x[0] || y[1] - x[1] || y[2] - x[2];
}
const removed = [], retained = [];
for (const entries of groups.values()) {
  entries.sort(newestFirst);
  const newest = [...new Set(entries.map(entry => entry.version))].slice(0, 2);
  const keep = new Set([version, '0.3.0', ...newest]);
  for (const entry of entries) {
    if (keep.has(entry.version)) { retained.push(entry.name); continue; }
    for (const name of [entry.name, `${entry.name}.blockmap`]) {
      if (!files.includes(name)) continue;
      removed.push({ name, bytes: fs.statSync(path.join(dist, name)).size });
    }
  }
}
if (args[0] === '--apply') {
  for (const { name } of removed) fs.unlinkSync(path.join(dist, name));
}
console.log(JSON.stringify({ applied: args[0] === '--apply', retained, removed,
  bytes: removed.reduce((sum, entry) => sum + entry.bytes, 0) }, null, 2));
