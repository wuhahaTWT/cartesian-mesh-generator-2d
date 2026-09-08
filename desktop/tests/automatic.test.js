'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { candidates, estimateSeconds } = require('../src/core/automatic');
const { validateJob, buildInvocation } = require('../src/core/job');
const { sampleById } = require('../src/core/samples');
test('dense hybrid retries layer thickness before reducing spatial resolution', () => {
  const options = candidates({ automatic: true, method: 'hybrid', density: 'dense',
    geometryPath: '/curve.xy', fluidRegion: 'exterior' }, sampleById('serpentine'), { bodySpan: 6 });
  assert.equal(options[1].maxLevel, options[0].maxLevel);
  assert.equal(options[1].boundaryLevel, options[0].boundaryLevel);
  assert.ok(options[1].firstThickness < options[0].firstThickness);
});
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
  const inner = candidates({automatic:true,method:'hybrid',geometryPath:'/in.xy',fluidRegion:'interior'}, null, {bodySpan:6});
  for (const option of inner) {
    const job = validateJob(option).job;
    assert.ok(buildInvocation(job, {}).args.includes('--fluid-region=interior'));
  }
});

test('nozzle density covers the interior and denser retries retain a consistent floor', () => {
  const request = { automatic: true, method: 'cutcell', geometryPath: '/nozzle.xy', fluidRegion: 'interior' };
  const sample = sampleById('nozzle');
  const normal = candidates(request, sample, { bodySpan: 6 });
  const dense = candidates({ ...request, density: 'dense' }, sample, { bodySpan: 6 });
  assert.equal(normal[0].wallCellsPerSpan, 128);
  assert.equal(dense[0].wallCellsPerSpan, 256);
  for (const choice of [...normal, ...dense]) {
    const { job } = validateJob(choice);
    assert.equal(job.budget.wallLevel - job.sizeField.farLevel, 1);
    assert.ok(job.sizeField.wallCellsPerSpan >= 128);
  }
  const exterior = candidates({ ...request, fluidRegion: 'exterior' }, sample, { bodySpan: 6 })[0];
  assert.equal(exterior.farFieldSpans, 6);
  assert.equal(exterior.farLevel, 0);
});


test('all built-in automatic candidates fit the safe wall budget', () => {
  const { SAMPLES } = require('../src/core/samples');
  for (const sample of SAMPLES) {
    for (const density of ['normal', 'dense']) {
      const options = candidates({ automatic: true, method: 'cutcell', density,
        geometryPath: '/input.xy', fluidRegion: sample.fluidRegion }, sample, { bodySpan: 6 });
      assert.ok(options.length <= 4);
      for (const option of options) assert.ok(validateJob(option).job.budget.feasible);
    }
  }
  const sample = sampleById('circle');
  const option = candidates({ automatic: true, method: 'cutcell', density: 'dense' }, sample, { bodySpan: 2 })[0];
  assert.equal(option.farFieldSpans, sample.sizeField.farFieldSpans);
  assert.ok(option.cellsPerLevel > sample.sizeField.cellsPerLevel);
});
