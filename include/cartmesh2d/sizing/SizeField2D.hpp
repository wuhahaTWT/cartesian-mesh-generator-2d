#pragma once

#include "cartmesh2d/quadtree/Quadtree2D.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

// A priori sizing controls.
//
// The pre-existing refinement policy has one geometric criterion — an
// `Intersected` leaf goes to `boundaryLevel` — so after 2:1 balance the finest
// zone is exactly one cell thick and the cell size then doubles every ring.
// That resolves the geometry but not the flow around it, and it leaves the only
// far-field knob (`minimumLevel`) costing 4x per level for resolution nobody
// needs.  This layer expresses the intent instead and compiles it down onto the
// existing primitives; it never introduces a second refinement mechanism.

// Widens every level into a band of `cellsPerLevel` cells rather than the single
// ring 2:1 balance produces on its own — the 2D analogue of snappyHexMesh's
// nCellsBetweenLevels.  `cellsPerLevel == 0` emits no bands at all and is
// therefore exactly the pre-existing boundary-only behaviour.
struct WallDistanceSizing2D {
    std::size_t cellsPerLevel = 3;
    std::size_t farLevel = 0;
};

// Raises the wall level only where the boundary actually turns, using the
// circumradius R of three consecutive wall vertices and requiring
// h <= R / cellsPerRadius.  Collinear triples request nothing.
struct CurvatureSizing2D {
    double cellsPerRadius = 8.0;
};

// Guarantees `cellsAcrossGap` cells across the gap between two wall stretches
// that face each other.  This is the upstream form of the narrow-gap layer/gap
// collision: sizing the gap beats repairing the cells it produced.
struct ProximitySizing2D {
    double cellsAcrossGap = 4.0;
    // Two stretches face each other across fluid only when the shortest connector
    // between them runs close to normal to both.  On a smooth convex wall the
    // connector between two nearby segments is almost tangential instead, which is
    // how a single closed loop avoids reporting its own chords as gaps — an index
    // window alone cannot distinguish the two cases.
    double maxTangentialCosine = 0.5;
    // Segments sharing a vertex have a zero-length connector with no direction, so
    // they are excluded outright rather than measured.
    std::size_t sharedVertexIndexWindow = 1;
};

// The wake is the one sizing source that genuinely depends on the freestream
// direction, so it is requested rather than inferred.  Everything else here is
// direction-free because a 2D body's perturbation magnitude decays isotropically
// (doublet ~ (a/r)^2, vortex ~ (a/r)).
struct WakeSizing2D {
    double angleOfAttackDeg = 0.0;
    double downstreamSpans = 6.0;
    double halfWidthSpans = 0.6;
    // A wake region is large — six body spans long by more than one wide — so it is
    // sized relative to the body, not the wall.  One level below the wall measures
    // out at ~70k extra leaves on a circle; four is ~1k for the same footprint.
    std::size_t levelsBelowWall = 4;
    // The refinement primitive is axis-aligned, so a wake at a non-zero angle is
    // covered by a staircase of boxes along the rotated centreline instead of one
    // oversized bounding box.
    std::size_t slices = 8;
};

struct SizeFieldPolicy2D {
    // Domain half-extent in body spans.  External-aero practice is 5-10 body
    // lengths upstream and 10-20 downstream; octree gradation makes this close to
    // free, so the default is generous rather than minimal.
    double farFieldSpans = 10.0;
    // Body span divided by the target wall cell size.  The tree level is derived
    // from this, so wall resolution stops changing when the domain does.
    //
    // Unset means "as fine as the construction can currently be trusted", i.e. the
    // wall lands exactly on `maxSafeWallLevel`.  Deriving it that way keeps the
    // defaults consistent for every geometry: fixing both a cell count and a ceiling
    // would leave the pair one rounding step away from refusing itself.
    std::optional<double> wallCellsPerSpan;
    // Deepest wall level this field will resolve to without an explicit opt-in.
    //
    // This is a measured construction limit, not a taste setting.  Holding the
    // sizing request fixed and perturbing the body scale by fractions of a percent
    // (tools/verification/alignment_sensitivity.py) gives a 20/20 pass rate at wall
    // level 8 for both a 1.5x and a 21x domain, but 11/20 at level 12 with the 21x
    // domain and 0/20 at level 12 with the 1.5x one.  The failures are discrete
    // events where a wall vertex grazes a grid line, and they get more likely with
    // depth exactly as docs/CURRENT_STATE_CN.md section 4A predicts.  Level 11 is
    // the deepest the W1 evidence actually covers, so that is the default ceiling.
    std::size_t maxSafeWallLevel = 11;
    // Crossing the ceiling has to be asked for, so a coarse `wallCellsPerSpan`
    // typo cannot quietly produce a mesh that fails half the time.
    bool allowUnsafeWallLevel = false;
    std::optional<WallDistanceSizing2D> wallDistance;
    std::optional<CurvatureSizing2D> curvature;
    std::optional<ProximitySizing2D> proximity;
    std::optional<WakeSizing2D> wake;
};

struct ResolvedSizeField2D {
    // Square by construction: `Quadtree2D` divides width and height by the same
    // power of two, so a non-square domain yields non-square cells at every level
    // (a 1 x 0.12 airfoil at padding 0.3 gives every cell an aspect of ~2.2).
    Domain2D domain;
    double bodySpan = 0.0;
    double domainSpan = 0.0;
    double wallCellSize = 0.0;
    // Derived, not requested: the maximum of the wall, curvature and proximity
    // requirements, capped by `maxLevelCap`.
    std::size_t maxLevel = 0;
    std::size_t wallLevel = 0;
    std::size_t curvatureLevel = 0;
    std::size_t proximityLevel = 0;
    bool levelCapReached = false;
    QuadtreeRefinementPolicy2D refinement;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

[[nodiscard]] ResolvedSizeField2D resolveSizeField2D(
    const SizeFieldPolicy2D& policy, const BoundaryRegion2D& boundary,
    std::size_t maxLevelCap = 28, const TolerancePolicy& tol = {});

// Smallest level whose cell is no larger than `targetSize`, clamped to the cap.
[[nodiscard]] std::size_t sizeFieldLevelForSize2D(
    double domainSpan, double targetSize, std::size_t maxLevelCap) noexcept;

[[nodiscard]] std::string resolvedSizeFieldToJson(
    const ResolvedSizeField2D& resolved, int indentSpaces = 2);

} // namespace cartmesh2d
