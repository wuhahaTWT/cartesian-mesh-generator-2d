'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { candidates, estimateSeconds } = require('../src/core/automatic');
const { validateJob } = require('../src/core/job');
const { sampleById } = require('../src/core/samples');
test('automatic nozzle keeps interior semantics and bounded validated choices', () => {
  const sample = sampleById('nozzle');
  assert.equal(sample.fluidRegion, 'interior');
  const options = candidates({ automatic: true, method: 'cutcell', geometryPath: '/nozzle.xy', fluidRegion: 'interior' }, sample, { bodySpan: 6 });
  assert.ok(options.length <= 4);
  for (const option of options) {
    assert.equal(validateJob(option).job.fluidRegion, 'interior');
    assert.equal(option.allowUnsafeWallLevel, false);
    assert.equal(option.wake, null);
  }
});
test('sharp hybrid starts with the documented successful layer recipe', () => {
  const choices = candidates({ automatic: true, method: 'hybrid', fluidRegion: 'exterior', geometryPath: '/sharp.xy' }, sampleById('sharp_trailing_edge'), { bodySpan: 4 });
  const job = validateJob(choices[0]).job;
  assert.deepEqual([job.maxLevel,job.boundaryLevel,job.nLayers,job.firstThickness,job.growthRatio], [8,8,4,0.012,1.15]);
  assert.ok(estimateSeconds(choices[0]).seconds >= 30);
});
test('manual requests stay unchanged and hybrid cannot silently mesh the wrong side', () => {
  const request = { method: 'cutcell', automatic: false };
  assert.deepEqual(candidates(request, null, {}), [request]);
  assert.throws(() => validateJob({method:'hybrid',geometryPath:'/in.xy',fluidRegion:'interior'}), /只支持外流/);
});
