'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { run } = require('../src/core/process');

test('cancellation stops a running worker and permits a fresh operation', async () => {
  const controller = new AbortController();
  const pending = run(process.execPath, ['-e', 'console.log("ready"); setInterval(()=>{},1000)'],
    () => controller.abort(), controller.signal);
  await assert.rejects(pending, /已取消/);
  assert.match((await run(process.execPath, ['-e', 'console.log("next")'])).stdout, /next/);
});
test('pre-cancelled work does not spawn and process errors retain diagnostics', async () => {
  const controller = new AbortController(); controller.abort();
  await assert.rejects(run('/does-not-exist', [], () => {}, controller.signal), /已取消/);
  await assert.rejects(run(process.execPath, ['-e', 'console.error("bad input");process.exit(2)']), /bad input/);
});
test('an automatic attempt has a time limit and leaves later work usable', async () => {
  await assert.rejects(run(process.execPath, ['-e', 'setInterval(()=>{},1000)'], () => {}, undefined, 80), /超时/);
  assert.match((await run(process.execPath, ['-e','console.log("alive")'])).stdout, /alive/);
});
test('selected nonzero result codes can carry a valid result', async () => {
  const result = await run(process.execPath, ['-e', 'console.log("usable");process.exit(2)'],
    () => {}, undefined, 0, [0, 2]);
  assert.equal(result.code, 2);
  assert.match(result.stdout, /usable/);
});
test('line callbacks wait for complete records across output chunks', async () => {
  const lines = [];
  await run(process.execPath, ['-e', 'process.stdout.write("{\\\"type\\\":");setTimeout(()=>process.stdout.end("\\\"flow-progress\\\"}\\n"),20)'],
    line => lines.push(line));
  assert.deepEqual(lines, ['{"type":"flow-progress"}']);
});
