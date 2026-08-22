#pragma once

#include "cartmesh2d/quadtree/Quadtree2D.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class CutCellKind {
    Empty,
    Full,
    Cut,
    Unsupported
};

// Physical interpretation of a closed BoundaryLoop.
// Exterior is the product default: the loop is a solid wall/obstacle and
// the CFD fluid occupies Domain2D minus the loop interior, matching cartmesh 3D.
// Interior is available explicitly for internal-flow/duct-style fixtures.
enum class FluidRegion2D {
    Exterior,
    Interior
};

enum class CutCellIssueCode {
    InvalidBackgroundCell,
    InvalidBoundary,
    DegeneratePolygon,
    MissingEmbeddedBoundary,
    MultipleEmbeddedComponents,
    AreaOutOfRange
};

struct CutCellIssue2D {
    CutCellIssueCode code;
    std::string message;
};

struct CutCell2D {
    std::size_t sourceId = 0;
    std::uint64_t sourceKey = 0;
    AABB2D backgroundBounds;
    CutCellKind kind = CutCellKind::Empty;
    Polygon2D fluidPolygon;
    double area = 0.0;
    double areaFraction = 0.0;
    std::optional<Point2D> centroid;
    std::vector<Segment2D> embeddedBoundary;
    std::vector<CutCellIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && kind != CutCellKind::Unsupported;
    }
};

// Default CFD semantics: BoundaryLoop is a solid obstacle and the retained
// fluid is on the EXTERIOR side of the loop.
[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

// Explicit physical-side override. Use Interior only for a deliberately
// internal-flow domain represented by the loop interior.
[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    FluidRegion2D fluidRegion,
    const TolerancePolicy& tol = {});

[[nodiscard]] CutCell2D buildCutCell(
    const QuadtreeLeaf2D& leaf,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] CutCell2D buildCutCell(
    const QuadtreeLeaf2D& leaf,
    const BoundaryLoop& boundary,
    FluidRegion2D fluidRegion,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
