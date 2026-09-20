'use strict';

const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const test = require('node:test');
const { readCheckpointMetadata } = require('../src/core/flow-checkpoint');

const SERIALIZED = [
  'CARTMESH2D_FLOW_CHECKPOINT 1',
  'DISCRETIZATION Euler-RC-v2',
  'CONFIG "channel" 0.01 1 limited-linear symmetric 0',
  'CELLS 1',
  'CELL 0 0.5 0.5 1 4 0 1 2 3',
  'FACES 1',
  'FACE 0 0 - 0 0 0 1 0 0 0 1 1',
  'TIME 0.125',
  'U 1 0',
  'V 1 0',
  'P 1 0',
  'FLUX 1 0',
  'END',
  ''
].join('\n');

async function withFile(text, callback) {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'cartmesh2d-checkpoint-'));
  const file = path.join(directory, 'restart.flow');
  try {
    await fs.writeFile(file, text);
    return await callback(file);
  } finally {
    await fs.rm(directory, { recursive: true, force: true });
  }
}

test('reads native checkpoint metadata without loading state vectors', async () => {
  const result = await withFile(SERIALIZED, readCheckpointMetadata);
  assert.deepEqual(result, {
    case: 'channel', nu: 0.01, speed: 1, convection: 'limited-linear',
    viscousStress: 'symmetric', outletBackflow: 'reject', time: 0.125, fileName: 'restart.flow'
  });
});

test('restores the face-frame momentum scheme without changing the legacy scheme', async () => {
  const result = await withFile(SERIALIZED.replace('limited-linear', 'face-limited-linear'), readCheckpointMetadata);
  assert.equal(result.convection, 'face-limited-linear');
});

test('reads v2 outlet backflow mode and preserves legacy v1 default', async () => {
  const v2 = SERIALIZED.replace('CARTMESH2D_FLOW_CHECKPOINT 1', 'CARTMESH2D_FLOW_CHECKPOINT 2')
    .replace('symmetric 0', 'symmetric 0 normal-inlet');
  const result = await withFile(v2, readCheckpointMetadata);
  assert.equal(result.outletBackflow, 'normal-inlet');
  await assert.rejects(withFile(v2.replace(' symmetric 0 normal-inlet', ' symmetric 0'), readCheckpointMetadata), /Invalid flow checkpoint metadata/);
  await assert.rejects(withFile(v2.replace('normal-inlet', 'unsupported'), readCheckpointMetadata), /Invalid flow checkpoint metadata/);
  await assert.rejects(withFile(SERIALIZED.replace('symmetric 0', 'symmetric 0 normal-inlet'), readCheckpointMetadata), /Invalid flow checkpoint metadata/);
});

test('accepts supported external, cavity and curved duct configurations', async () => {
  for (const scenario of ['external', 'cavity', 'duct']) {
    const text = SERIALIZED.replace('"channel"', `"${scenario}"`).replace('TIME 0.125', 'TIME 0');
    const result = await withFile(text, readCheckpointMetadata);
    assert.equal(result.case, scenario);
    assert.equal(result.time, 0);
  }
});

test('rejects unsupported, corrupt, and unsafe metadata', async () => {
  const cases = [
    ['empty', ''],
    ['header', SERIALIZED.replace('CARTMESH2D_FLOW_CHECKPOINT 1', 'TG 1')],
    ['case', SERIALIZED.replace('"channel"', '"taylor-green"')],
    ['manufactured', SERIALIZED.replace('"channel"', '"manufactured"')],
    ['discretization', SERIALIZED.replace('Euler-RC-v2', 'Euler')],
    ['nu', SERIALIZED.replace('0.01 1 ', 'NaN 1 ')],
    ['speed', SERIALIZED.replace('0.01 1 ', '0.01 0 ')],
    ['slope', SERIALIZED.replace('symmetric 0', 'symmetric 1')],
    ['counts', SERIALIZED.replace('CELLS 1', 'CELLS 0')],
    ['time', SERIALIZED.replace('TIME 0.125', 'TIME NaN')],
    ['missing time', SERIALIZED.replace('TIME 0.125', 'STEP 1')]
  ];
  for (const [name, text] of cases) {
    await assert.rejects(withFile(text, readCheckpointMetadata), /Invalid flow checkpoint metadata/, name);
  }
});

test('rejects directories and non-regular paths', async () => {
  const directory = await fs.mkdtemp(path.join(os.tmpdir(), 'cartmesh2d-checkpoint-dir-'));
  try {
    await assert.rejects(readCheckpointMetadata(directory), /regular file/);
  } finally {
    await fs.rm(directory, { recursive: true, force: true });
  }
});
