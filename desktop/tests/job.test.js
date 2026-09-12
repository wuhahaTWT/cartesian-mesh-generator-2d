'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');

const { validateJob, buildInvocation } = require('../src/core/job');

const cutcell = extra => ({
  method: 'cutcell', geometryPath: '/in/body.xy', outputDirectory: '/out',
  farFieldSpans: 10, wallCellsPerSpan: 64, cellsPerLevel: 3, ...extra
});
const paths = { xyPath: '/out/body.xy', prefix: '/out/body', casePath: '/out/body-openfoam' };

test('both methods carry the same explicit reference and relative sizes to native code', () => {
  const request = { sizingMode: 'relative', geometryPath: '/in/body.xy', referenceLength: 0.2,
    wallRelativeSize: 0.01, backgroundRelativeSize: 0.1, farFieldSpans: 0.5, cellsPerLevel: 3,
    nLayers: 3, firstLayerRelativeSize: 0.001, growthRatio: 1.2, extrusionRelativeSize: 0.01 };
  for (const method of ['cutcell', 'hybrid']) {
    const { job } = validateJob({ ...request, method });
    const { args } = buildInvocation(job, paths);
    const value = key => args[args.indexOf(key)+1];
    assert.equal(value('--reference-length'), '0.2');
    assert.equal(value('--wall-relative-size'), '0.01');
    assert.equal(value('--background-relative-size'), '0.1');
    assert.ok(!args.includes('--wall-cells-per-span'));
    if (method === 'hybrid') {
      assert.equal(value('--first-layer-relative-size'), '0.001');
      assert.equal(args[10], '0.002', 'extrusion ratio is converted using the same reference');
    }
    assert.throws(() => validateJob({ ...request, method, referenceLength: 0 }), /参考长度/);
    assert.throws(() => validateJob({ ...request, method, backgroundRelativeSize: 0.001 }), /背景尺寸/);
  }
});

test('the pure path argv is the one that was verified against the CLI', () => {
  const { job } = validateJob(cutcell());
  const { executable, args } = buildInvocation(job, paths);
  assert.equal(executable, 'cartmesh2d_cli');
  assert.deepEqual(args, ['/out/body.xy', '/out/body', '8', '0.25', '0.15', 'exterior',
    '/out/body-openfoam', '0', '0',
    '--size-field', '--far-field-spans', '10', '--wall-cells-per-span', '64',
    '--cells-per-level', '3', '--far-level', '0', '--max-safe-wall-level', '11']);
});

// The pure CLI builds and evaluates its solver partition inside `if (openFoamCase)`.
// Dropping the case directory would also drop the desktop Solver report.
test('a case directory is always requested so the quality report exists', () => {
  const { job } = validateJob(cutcell());
  assert.equal(buildInvocation(job, paths).args.includes('/out/body-openfoam'), true);
});

test('the dry run asks for no case and no mesh', () => {
  const { job } = validateJob(cutcell());
  const { args } = buildInvocation(job, paths, { dryRun: true });
  assert.ok(args.includes('--size-field-only'));
  assert.ok(!args.includes('--size-field'));
  assert.equal(args[6], '-');
});

test('optional sizing sources appear only when asked for', () => {
  const plain = buildInvocation(validateJob(cutcell()).job, paths).args;
  assert.ok(!plain.includes('--gap-cells'));
  assert.ok(!plain.includes('--curvature-cells-per-radius'));
  assert.ok(!plain.includes('--wake'));

  const full = buildInvocation(validateJob(cutcell({
    gapCells: 4, curvatureCellsPerRadius: 8,
    wake: { angleOfAttackDeg: 4, downstreamSpans: 8, halfWidthSpans: 0.7, levelsBelowWall: 4 }
  })).job, paths).args;
  assert.deepEqual(full.slice(full.indexOf('--curvature-cells-per-radius')),
    ['--curvature-cells-per-radius', '8', '--gap-cells', '4', '--wake', '4', '8', '0.7', '4']);
});

test('an infeasible depth is refused before a process starts, with the fix in the message', () => {
  assert.throws(() => validateJob(cutcell({ farFieldSpans: 15, wallCellsPerSpan: 256 })),
    error => /level 13/.test(error.message) && /3\.5/.test(error.message) && /1\/64/.test(error.message));
});

test('the ceiling can be crossed, but only on purpose', () => {
  const { job } = validateJob(cutcell({ farFieldSpans: 15, wallCellsPerSpan: 256, allowUnsafeWallLevel: true }));
  assert.ok(buildInvocation(job, paths).args.includes('--allow-unsafe-wall-level'));
});

test('the hybrid path is capped at its own measured ceiling', () => {
  const hybrid = extra => ({
    method: 'hybrid', geometryPath: '/in/body.xy', outputDirectory: '/out',
    maxLevel: 6, minimumLevel: 3, boundaryLevel: 6, nLayers: 4,
    firstThickness: 0.02, growthRatio: 1.2, domainPadding: 1.0, ...extra
  });
  const { job } = validateJob(hybrid());
  assert.deepEqual(buildInvocation(job, paths).args,
    ['/out/body.xy', '/out/body', '6', '3', '6', '4', '0.02', '1.2', '1',
     '/out/body-openfoam', '0.01', '--fluid-region=exterior', '--small-alpha=0.1']);
  // docs/CURRENT_STATE_CN.md section 2: level 9 fails on this path.
  assert.throws(() => validateJob(hybrid({ maxLevel: 9 })), /余域最高层级/);
  // The layer has to attach at or above the far-field floor and at or below the cap.
  assert.throws(() => validateJob(hybrid({ minimumLevel: 5, boundaryLevel: 4 })), /壁面连接层级/);
  assert.throws(() => validateJob(hybrid({ minimumLevel: 7 })), /全域最低层级/);
});

test('unknown methods and missing paths are refused', () => {
  assert.throws(() => validateJob({ method: 'magic' }), /未知网格方法/);
  assert.throws(() => validateJob(cutcell({ geometryPath: '' })), /几何文件/);
  assert.equal(validateJob(cutcell({ outputDirectory: '' })).job.outputDirectory, '');
  assert.deepEqual(buildInvocation(validateJob(cutcell()).job, paths).cm2dCandidates, ['/out/body.solver.cm2d']);
  assert.throws(() => validateJob(cutcell({ cellsPerLevel: 2.5 })), /每级带宽/);
});

test('hand-placed regions convert from body spans to absolute --refine-box', () => {
  const { job } = validateJob(cutcell({
    refineBoxes: [{ xmin: 0.6, xmax: 6, ymin: -0.8, ymax: 0.8, levelsBelowWall: 4 }]
  }));
  // far 10 / wall 64 resolves to wall level 11, so 4 levels below is 7.
  const frame = { centreX: 10, centreY: -5, bodySpan: 2 };
  const args = buildInvocation(job, { ...paths, frame }).args;
  const at = args.indexOf('--refine-box');
  assert.deepEqual(args.slice(at, at + 6),
    ['--refine-box', '11.2', '-6.6', '22', '-3.4', '7']);
});

test('a region is clamped into the levels the primitive accepts', () => {
  const frame = { centreX: 0, centreY: 0, bodySpan: 1 };
  const deep = validateJob(cutcell({
    refineBoxes: [{ xmin: 0, xmax: 1, ymin: 0, ymax: 1, levelsBelowWall: 0 }]
  })).job;
  const shallow = validateJob(cutcell({
    refineBoxes: [{ xmin: 0, xmax: 1, ymin: 0, ymax: 1, levelsBelowWall: 20 }]
  })).job;
  const level = job => {
    const args = buildInvocation(job, { ...paths, frame }).args;
    return Number(args[args.indexOf('--refine-box') + 5]);
  };
  assert.equal(level(deep), 11, 'never deeper than the wall');
  assert.equal(level(shallow), 1, 'refine() rejects level 0, so it saturates at 1');
});

test('an inverted region is refused with the region named', () => {
  assert.throws(() => validateJob(cutcell({
    refineBoxes: [{ xmin: 6, xmax: 0.6, ymin: -0.8, ymax: 0.8, levelsBelowWall: 4 }]
  })), /加密区 1/);
});

test('the dry run leaves regions out, since it only resolves the field', () => {
  const { job } = validateJob(cutcell({
    refineBoxes: [{ xmin: 0.6, xmax: 6, ymin: -0.8, ymax: 0.8, levelsBelowWall: 4 }]
  }));
  const frame = { centreX: 0, centreY: 0, bodySpan: 1 };
  assert.ok(!buildInvocation(job, { ...paths, frame }, { dryRun: true }).args.includes('--refine-box'));
});
