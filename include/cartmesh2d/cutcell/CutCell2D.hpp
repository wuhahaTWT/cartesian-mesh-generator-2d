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

[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] CutCell2D buildCutCell(
    const QuadtreeLeaf2D& leaf,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
