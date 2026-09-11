'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const { runBudget } = require('../src/core/budget-runner');

const request = extra => ({
  targetCells: 10000,
  method: 'cutcell',
  fluidRegion: 'exterior',
  farFieldSpans: 0.5,
  ...extra
});

const frame = { bodySpan: 1, fluidArea: 1 };

function payload(cells, actualMethod = 'cutcell') {
  return { mesh: { cells: Array.from({ length: cells }, () => ({})) },
    result: { actualMethod } };
}

test('selects the nearest successful payload across out-of-range attempts', async () => {
  const counts = [5000, 16000, 10500];
  let calls = 0;
  const result = await runBudget(request(), frame, {
    generate: async () => payload(counts[calls++])
  });

  assert.equal(calls, 3);
  assert.equal(result.mesh.cells.length, 10500);
  assert.equal(result.cellBudget.actualCells, 10500);
  assert.equal(result.cellBudget.reached, true);
  assert.deepEqual(result.attempts.map(attempt => [attempt.success, attempt.actualCells]),
    [[true, 5000], [true, 16000], [true, 10500]]);
});

test('keeps an earlier successful result when a later attempt fails', async () => {
  let calls = 0;
  const result = await runBudget(request(), frame, {
    generate: async () => {
      calls += 1;
      if (calls === 2) throw new Error('质量检查失败');
      return payload(4000);
    }
  });

  assert.equal(calls, 2);
  assert.equal(result.mesh.cells.length, 4000);
  assert.deepEqual(result.attempts.map(attempt => attempt.success), [true, false]);
  assert.equal(result.attempts[1].reason, '质量检查失败');
});

test('rejects a hybrid result after bounded cutcell fallback retries', async () => {
  let calls = 0;
  await assert.rejects(runBudget(request({ method: 'hybrid' }), frame, {
    generate: async () => {
      calls += 1;
      return payload(10000, 'cutcell');
    }
  }), /未保留所选混合网格方法/);
  assert.equal(calls, 3);
});

test('retries a hybrid layer fallback with a thinner first layer', async () => {
  const choices = [];
  const hybridRequest = request({
    method: 'hybrid',
    nLayers: 5,
    growthRatio: 1.25,
    extrusionRelativeSize: 0.02,
    domainPadding: 0.75
  });
  const result = await runBudget(hybridRequest, frame, {
    generate: async choice => {
      choices.push({ ...choice });
      return choices.length === 1 ? payload(10000, 'cutcell') : payload(10000, 'hybrid');
    }
  });

  assert.equal(choices.length, 2);
  assert.equal(choices[1].firstLayerRelativeSize,
    choices[0].firstLayerRelativeSize * 0.65);
  for (const key of ['method', 'nLayers', 'growthRatio', 'extrusionRelativeSize', 'domainPadding']) {
    assert.equal(choices[1][key], choices[0][key], `${key} must survive a layer retry`);
  }
  assert.deepEqual(result.attempts.map(attempt => attempt.success), [false, true]);
  assert.match(result.attempts[0].reason, /未保留所选混合网格方法/);
  assert.equal(result.attempts[1].actualCells, 10000);
  assert.equal(result.cellBudget.stoppedReason, null);
});

test('rejects after three bounded hybrid layer fallback attempts', async () => {
  const choices = [];
  await assert.rejects(runBudget(request({ method: 'hybrid', nLayers: 4, growthRatio: 1.2 }), frame, {
    generate: async choice => {
      choices.push({ ...choice });
      return payload(10000, 'cutcell');
    }
  }), /未保留所选混合网格方法/);

  assert.equal(choices.length, 3);
  assert.equal(choices[1].firstLayerRelativeSize,
    choices[0].firstLayerRelativeSize * 0.65);
  assert.equal(choices[2].firstLayerRelativeSize,
    choices[0].firstLayerRelativeSize * 1.5);
});

test('does not return a result when cancellation arrives after generation', async () => {
  const controller = new AbortController();
  await assert.rejects(runBudget(request(), frame, {
    signal: controller.signal,
    generate: async () => {
      controller.abort();
      return payload(10000);
    }
  }), /操作已取消/);
});

test('stops after the first in-range result', async () => {
  let calls = 0;
  const progress = [];
  const result = await runBudget(request(), frame, {
    progress: event => progress.push(event),
    generate: async () => {
      calls += 1;
      return payload(10000);
    }
  });

  assert.equal(calls, 1);
  assert.equal(progress.length, 1);
  assert.equal(result.attempts.length, 1);
  assert.equal(result.cellBudget.reached, true);
});
