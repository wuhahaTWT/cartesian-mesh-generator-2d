#pragma once

#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"

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
    std::vector<CanonicalizedIntersection2D> canonicalizedIntersections;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && kind != CutCellKind::Unsupported;
    }
};

// Single-component API. This remains convenient for unit tests and simple
// leaves. If a leaf genuinely contains more than one disconnected fluid
// component, this API returns Unsupported rather than silently dropping one.
[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    FluidRegion2D fluidRegion,
    const TolerancePolicy& tol = {});

[[nodiscard]] CutCell2D buildCutCell(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryRegion2D& boundary,
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior,
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

[[nodiscard]] CutCell2D buildCutCell(
    const QuadtreeLeaf2D& leaf,
    const BoundaryRegion2D& boundary,
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior,
    const TolerancePolicy& tol = {});

// Product/solver API. A single Cartesian/Quadtree leaf may contain multiple
// disconnected fluid components for a strongly concave solid. Each component
// must become its own solver cell instead of being discarded or bridged.
// A true local hole (a solid loop fully enclosed by one leaf) is still
// Unsupported in the current simple-polygon topology model and is reported
// explicitly; refinement or future polygon-with-holes support is required.
[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryLoop& boundary,
    FluidRegion2D fluidRegion,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const AABB2D& backgroundBounds,
    CellClass classification,
    const BoundaryRegion2D& boundary,
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const QuadtreeLeaf2D& leaf,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const QuadtreeLeaf2D& leaf,
    const BoundaryLoop& boundary,
    FluidRegion2D fluidRegion,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::vector<CutCell2D> buildCutCells(
    const QuadtreeLeaf2D& leaf,
    const BoundaryRegion2D& boundary,
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
