'use strict';

// Mirror of the arithmetic in src/sizing/SizeField2D.cpp so the UI can show the
// consequence of a sizing request before spending a solve on it.
//
// Only the wall level is predictable from the flags alone: it depends on the far
// field and the wall cell count and nothing else.  The curvature and proximity
// levels depend on the geometry, so those come back from --size-field-only.

// resolveSizeField2D: domainSpan = bodySpan * (1 + 2 * farFieldSpans).
const domainSpansPerBody = farFieldSpans => 1 + 2 * farFieldSpans;

// sizeFieldLevelForSize2D: the smallest level whose cell is no larger than the
// target, i.e. ceil(log2(domainSpan / targetSize)), clamped to the cap.
function levelForSize(domainSpan, targetSize, cap = 28) {
  if (!(domainSpan > 0) || !(targetSize > 0) || !Number.isFinite(targetSize)) return 0;
  if (targetSize >= domainSpan) return 0;
  const exact = Math.log2(domainSpan / targetSize);
  if (!(exact > 0)) return 0;
  return Math.min(cap, Math.ceil(exact));
}

// The body span cancels: level = ceil(log2((1 + 2*far) * wallCellsPerSpan)).  So a
// request is feasible or not before any geometry is loaded, which is what makes a
// live budget readout possible.
function wallLevelFor(farFieldSpans, wallCellsPerSpan, cap = 28) {
  return levelForSize(domainSpansPerBody(farFieldSpans), 1 / wallCellsPerSpan, cap);
}

// Largest wallCellsPerSpan that still lands on or under `ceiling`.
function wallCellsPerSpanCeiling(farFieldSpans, ceiling) {
  return Math.pow(2, ceiling) / domainSpansPerBody(farFieldSpans);
}

// The single fact the far-field question reduces to: tree level is measured from
// the domain, so pushing the boundary out and refining the wall both consume the
// same depth budget.  Report the trade rather than letting the mesher refuse.
function describeBudget({ farFieldSpans, wallCellsPerSpan, safeWallLevel }) {
  const wallLevel = wallLevelFor(farFieldSpans, wallCellsPerSpan);
  const affordableCells = wallCellsPerSpanCeiling(farFieldSpans, safeWallLevel);
  return {
    wallLevel,
    feasible: wallLevel <= safeWallLevel,
    domainSpansPerBody: domainSpansPerBody(farFieldSpans),
    // Both ways out, so the choice stays the user's rather than the tool's.
    maxWallCellsPerSpan: Math.max(1, Math.pow(2, Math.floor(Math.log2(affordableCells)))),
    maxFarFieldSpans: Math.max(0.25, (Math.pow(2, safeWallLevel) / wallCellsPerSpan - 1) / 2)
  };
}

module.exports = {
  domainSpansPerBody,
  levelForSize,
  wallLevelFor,
  wallCellsPerSpanCeiling,
  describeBudget
};
