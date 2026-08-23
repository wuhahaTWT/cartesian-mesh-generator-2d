#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <vector>

namespace cartmesh2d {

struct Domain2D {
    AABB2D bounds;

    [[nodiscard]] bool valid(const TolerancePolicy& tol = {}) const noexcept;
    [[nodiscard]] double width() const noexcept { return bounds.max.x - bounds.min.x; }
    [[nodiscard]] double height() const noexcept { return bounds.max.y - bounds.min.y; }
};

enum class CellClass { Outside, Inside, Intersected };

struct CartesianCell2D {
    std::size_t id = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    AABB2D bounds;
    CellClass classification = CellClass::Outside;

    [[nodiscard]] Point2D center() const noexcept;
};

class UniformCartesianGrid2D {
public:
    UniformCartesianGrid2D(Domain2D domain, std::size_t nx, std::size_t ny);

    [[nodiscard]] static UniformCartesianGrid2D fromTargetSpacing(
        Domain2D domain, double targetDx, double targetDy);

    [[nodiscard]] const Domain2D& domain() const noexcept { return domain_; }
    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }
    [[nodiscard]] double dx() const noexcept { return dx_; }
    [[nodiscard]] double dy() const noexcept { return dy_; }
    [[nodiscard]] const std::vector<CartesianCell2D>& cells() const noexcept { return cells_; }
    [[nodiscard]] std::vector<CartesianCell2D>& cells() noexcept { return cells_; }
    [[nodiscard]] const CartesianCell2D& cell(std::size_t i, std::size_t j) const;

private:
    Domain2D domain_;
    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    double dx_ = 0.0;
    double dy_ = 0.0;
    std::vector<CartesianCell2D> cells_;
};

struct ClassificationSummary {
    std::size_t inside = 0;
    std::size_t outside = 0;
    std::size_t intersected = 0;
};

[[nodiscard]] bool segmentIntersectsClosedAABB(
    const Segment2D& segment, const AABB2D& box,
    const TolerancePolicy& tol = {}) noexcept;

[[nodiscard]] CellClass classifyCartesianCell(
    const CartesianCell2D& cell, const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] ClassificationSummary classifyGrid(
    UniformCartesianGrid2D& grid, const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
