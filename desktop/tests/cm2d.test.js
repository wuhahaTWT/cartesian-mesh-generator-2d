'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');

const { parseCm2d, embeddedBounds, levelHistogram, assignKeyLevels,
        assignSizeBands, PATCH } = require('../src/core/cm2d');
const { parseKeyValues, summarizeContract, summarizeSolverQuality } = require('../src/core/report');

// Hand-built to the exact layout MeshIO2D.cpp writes:
//   CELLS: id sourceId sourceKey area nv v... ne e...
// sourceKey packs (morton << 6) | level, so 0b..._000110 is level 6 and 519 is level 7.
const MESH = [
  'CM2D 1',
  'VERTICES 4',
  '0 0 0', '1 1 0', '2 1 1', '3 0 1',
  'EDGES 4',
  '0 0 1 0 -1 2',
  '1 1 2 0 1 0',
  '2 2 3 0 -1 1',
  '3 3 0 0 -1 2',
  'CELLS 2',
  '0 0 6 1 4 0 1 2 3 4 0 1 2 3',
  '1 1 519 0.25 3 0 1 2 3 0 1 2',
  'AUDIT 0 0 0 0 0 0 0'
].join('\n');

test('the reader follows the writer field for field', () => {
  const mesh = parseCm2d(MESH);
  assert.equal(mesh.vertices.length, 4);
  assert.equal(mesh.edges.length, 4);
  assert.equal(mesh.cells.length, 2);
  assert.deepEqual(mesh.cells[0].vertices, [0, 1, 2, 3]);
  assert.equal(mesh.cells[1].area, 0.25);
});

test('the level is the low six bits of the quadtree key', () => {
  const mesh = assignKeyLevels(parseCm2d(MESH));
  assert.equal(mesh.cells[0].level, 6);
  assert.equal(mesh.cells[1].level, 7);
  assert.equal(mesh.minLevel, 6);
  assert.equal(mesh.maxLevel, 7);
});

// HybridMesh2D.cpp writes `sourceKey = sourceId`, a running index, so reading a level
// out of it paints id-mod-64 noise over a symmetric mesh.  Size banding is the fallback.
test('size bands come from cell area, one band per factor of two', () => {
  const mesh = assignSizeBands(parseCm2d(MESH));
  // Areas 1 and 0.25: sqrt gives 1 and 0.5, exactly one octave apart.
  assert.equal(mesh.cells[0].level, 0);
  assert.equal(mesh.cells[1].level, 1);
});

test('size banding is blind to a sourceKey that is only an index', () => {
  // Same two cells, keys replaced by running ids as the hybrid path writes them.
  const asIndex = MESH.replace('0 0 6 1 4', '0 0 0 1 4').replace('1 1 519 0.25 3', '1 1 1 0.25 3');
  const byKey = assignKeyLevels(parseCm2d(asIndex));
  const bySize = assignSizeBands(parseCm2d(asIndex));
  assert.deepEqual(byKey.cells.map(cell => cell.level), [0, 1], 'the key now says nothing');
  assert.deepEqual(bySize.cells.map(cell => cell.level), [0, 1]);
  // The key-derived levels only agreed here by coincidence; a level-8 cell proves it.
  const shifted = MESH.replace('1 1 519 0.25 3', '1 1 8 0.25 3');
  assert.equal(assignKeyLevels(parseCm2d(shifted)).cells[1].level, 8);
  assert.equal(assignSizeBands(parseCm2d(shifted)).cells[1].level, 1);
});

test('the wall bounds come from embedded edges only', () => {
  const mesh = parseCm2d(MESH);
  assert.equal(mesh.edges[2].patch, PATCH.EMBEDDED);
  // Edge 2 runs (1,1)-(0,1): the wall box is that segment, not the whole domain.
  assert.deepEqual(embeddedBounds(mesh), { minX: 0, minY: 1, maxX: 1, maxY: 1 });
});

test('a mesh with no wall reports no wall bounds', () => {
  assert.equal(embeddedBounds(parseCm2d(MESH.replace('2 2 3 0 -1 1', '2 2 3 0 -1 2'))), null);
});

test('the histogram is ordered coarse to fine', () => {
  assert.deepEqual(levelHistogram(assignKeyLevels(parseCm2d(MESH))),
    [{ level: 6, count: 1 }, { level: 7, count: 1 }]);
});

test('a truncated or foreign file is refused', () => {
  assert.throws(() => parseCm2d('NOTCM2D 1\n'), /CM2D/);
  assert.throws(() => parseCm2d('CM2D 1\nVERTICES 0\nEDGES 0\nCELLS 0\n'), /没有可显示的单元/);
});

// cartmesh2d_cli prints one pair per line; cartmesh2d_hybrid_cli prints many per
// line.  Both have to reach the same place.
test('key=value output is read in both CLI styles', () => {
  const pure = parseKeyValues('stabilized_cells=3452\nquality_contract=FAIL\nsize_field_wall_level=11\n');
  assert.equal(pure.stabilized_cells, '3452');
  assert.equal(pure.size_field_wall_level, '11');

  const hybrid = parseKeyValues('hybrid_status=success cells=700 solver_cells=728 quality_contract=FAIL');
  assert.equal(hybrid.cells, '700');
  assert.equal(hybrid.solver_cells, '728');
  assert.equal(hybrid.quality_contract, 'FAIL');
});

test('a lone pair keeps a value containing spaces', () => {
  const values = parseKeyValues('cm2d=/Users/a b/mesh.cm2d');
  assert.equal(values.cm2d, '/Users/a b/mesh.cm2d');
});

test('an unrated cell type reports OBSERVED rather than passing silently', () => {
  const summary = summarizeContract({
    status: 'FAIL',
    by_cell_type: {
      cartesian: { status: 'FAIL', rated: true, cell_count: 2224, hard_issue_count: 58, preferred_issue_count: 139 },
      boundary_layer: { status: 'OBSERVED', rated: false, cell_count: 128, hard_issue_count: 0, preferred_issue_count: 0 }
    },
    issues: [
      { level: 'hard', metric: 'volume_ratio' },
      { level: 'hard', metric: 'volume_ratio' },
      { level: 'preferred', metric: 'face_weight' }
    ]
  });
  assert.equal(summary.status, 'FAIL');
  assert.equal(summary.byType.find(row => row.type === 'boundary_layer').status, 'OBSERVED');
  assert.equal(summary.metricCounts['hard:volume_ratio'], 2);
  assert.equal(summary.metricCounts['preferred:face_weight'], 1);
});

test('solver metrics are judged in the direction each limit runs', () => {
  const summary = summarizeSolverQuality({
    max_non_orthogonality_deg: 42.3, min_face_weight: 0.0749, min_volume_ratio: 0.0187
  }, false);
  const byKey = Object.fromEntries(summary.rows.map(row => [row.key, row.pass]));
  assert.equal(byKey.max_non_orthogonality_deg, true);
  assert.equal(byKey.min_face_weight, true, '0.0749 clears the 0.05 solver floor');
  assert.equal(byKey.min_volume_ratio, true, '0.0187 clears the 0.01 solver floor');
  assert.equal(summarizeSolverQuality({ min_face_weight: 0.02 }, false).rows[0].pass, false);
});
