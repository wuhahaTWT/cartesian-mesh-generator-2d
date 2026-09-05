'use strict';

// Reader for the native .cm2d text mesh.  Written by MeshIO2D.cpp as:
//
//   CM2D 1
//   VERTICES n   then  id x y
//   EDGES    n   then  id v0 v1 owner neighbour patch      (neighbour -1 on boundary)
//   CELLS    n   then  id sourceId sourceKey area nv v... ne e...
//   AUDIT ...
//
// sourceKey is the Quadtree key, whose low six bits are the level (Quadtree2D::makeKey
// packs `(morton << 6) | level`).  That is what makes it possible to colour the
// preview by refinement level and therefore to see the size field's gradation.

const PATCH = Object.freeze({ NONE: 0, EMBEDDED: 1, DOMAIN: 2, UNCLASSIFIED: 3 });

function parseCm2d(text) {
  const tokens = text.split(/\s+/).filter(Boolean);
  let cursor = 0;
  const word = () => tokens[cursor++];
  const value = () => Number(tokens[cursor++]);
  const expect = expected => {
    if (word() !== expected) throw new Error(`CM2D 文件缺少 ${expected} 段。`);
  };

  expect('CM2D');
  word();

  expect('VERTICES');
  const vertices = new Array(value());
  for (let i = 0; i < vertices.length; i++) {
    const id = value();
    vertices[id] = [value(), value()];
  }

  expect('EDGES');
  const edges = new Array(value());
  for (let i = 0; i < edges.length; i++) {
    value();
    edges[i] = { a: value(), b: value(), owner: value(), neighbour: value(), patch: value() };
  }

  expect('CELLS');
  const cells = new Array(value());
  let minLevel = Infinity;
  let maxLevel = -Infinity;
  for (let i = 0; i < cells.length; i++) {
    value();
    value();
    const sourceKey = value();
    const area = value();
    const cellVertices = new Array(value());
    for (let v = 0; v < cellVertices.length; v++) cellVertices[v] = value();
    const edgeTotal = value();
    cursor += edgeTotal;
    const level = sourceKey % 64;
    minLevel = Math.min(minLevel, level);
    maxLevel = Math.max(maxLevel, level);
    cells[i] = { level, area, vertices: cellVertices };
  }
  if (!vertices.length || !cells.length) throw new Error('CM2D 文件没有可显示的单元。');
  return { vertices, edges, cells, minLevel, maxLevel, bounds: boundsOf(vertices) };
}

function boundsOf(vertices) {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (const [x, y] of vertices) {
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  return { minX, minY, maxX, maxY };
}

// The wall is the only thing whose extent is known without the input geometry, and
// zooming to it is what the user wants first: with a 10-body-span far field the body
// is one seventeenth of the picture.
function embeddedBounds(mesh) {
  const wall = mesh.edges.filter(edge => edge.patch === PATCH.EMBEDDED);
  if (!wall.length) return null;
  const points = [];
  for (const edge of wall) points.push(mesh.vertices[edge.a], mesh.vertices[edge.b]);
  return boundsOf(points);
}

// Cells per level, so the panel can show where the budget actually went.
function levelHistogram(mesh) {
  const counts = new Map();
  for (const cell of mesh.cells) counts.set(cell.level, (counts.get(cell.level) || 0) + 1);
  return [...counts.entries()].sort((a, b) => a[0] - b[0]).map(([level, count]) => ({ level, count }));
}

module.exports = { parseCm2d, boundsOf, embeddedBounds, levelHistogram, PATCH };
