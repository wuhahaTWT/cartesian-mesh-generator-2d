#include "cartmesh2d/grid/CartesianGrid2D.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cartmesh2d {

bool Domain2D::valid(const TolerancePolicy& tol) const noexcept {
    if (!bounds.valid(tol)) return false;
    const double scale = std::max({1.0, std::abs(bounds.min.x), std::abs(bounds.max.x),
                                  std::abs(bounds.min.y), std::abs(bounds.max.y)});
    const double eps = tol.scale(scale);
    return width() > eps && height() > eps;
}

Point2D CartesianCell2D::center() const noexcept {
    return {(bounds.min.x + bounds.max.x) * 0.5,
            (bounds.min.y + bounds.max.y) * 0.5};
}

UniformCartesianGrid2D::UniformCartesianGrid2D(Domain2D domain,
                                               std::size_t nx,
                                               std::size_t ny)
    : domain_(domain), nx_(nx), ny_(ny) {
    if (!domain_.valid()) throw std::invalid_argument("invalid 2D Cartesian domain");
    if (nx_ == 0 || ny_ == 0) throw std::invalid_argument("grid dimensions must be positive");

    dx_ = domain_.width() / static_cast<double>(nx_);
    dy_ = domain_.height() / static_cast<double>(ny_);
    cells_.reserve(nx_ * ny_);
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const double x0 = domain_.bounds.min.x + static_cast<double>(i) * dx_;
            const double y0 = domain_.bounds.min.y + static_cast<double>(j) * dy_;
            const double x1 = (i + 1 == nx_) ? domain_.bounds.max.x : x0 + dx_;
            const double y1 = (j + 1 == ny_) ? domain_.bounds.max.y : y0 + dy_;
            cells_.push_back({j * nx_ + i, i, j, {{x0, y0}, {x1, y1}}, CellClass::Outside});
        }
    }
}

UniformCartesianGrid2D UniformCartesianGrid2D::fromTargetSpacing(
    Domain2D domain, double targetDx, double targetDy) {
    if (!domain.valid()) throw std::invalid_argument("invalid 2D Cartesian domain");
    if (!(targetDx > 0.0) || !(targetDy > 0.0)) {
        throw std::invalid_argument("target spacing must be positive");
    }
    const auto nx = static_cast<std::size_t>(std::ceil(domain.width() / targetDx));
    const auto ny = static_cast<std::size_t>(std::ceil(domain.height() / targetDy));
    return UniformCartesianGrid2D(domain, std::max<std::size_t>(1, nx),
                                  std::max<std::size_t>(1, ny));
}

const CartesianCell2D& UniformCartesianGrid2D::cell(std::size_t i, std::size_t j) const {
    if (i >= nx_ || j >= ny_) throw std::out_of_range("Cartesian cell index out of range");
    return cells_.at(j * nx_ + i);
}

bool segmentIntersectsClosedAABB(const Segment2D& segment, const AABB2D& box,
                                 const TolerancePolicy& tol) noexcept {
    if (!box.valid(tol)) return false;
    if (box.contains(segment.a, tol) || box.contains(segment.b, tol)) return true;

    const Point2D bl{box.min.x, box.min.y};
    const Point2D br{box.max.x, box.min.y};
    const Point2D tr{box.max.x, box.max.y};
    const Point2D tl{box.min.x, box.max.y};
    const Segment2D edges[4] = {{bl, br}, {br, tr}, {tr, tl}, {tl, bl}};
    for (const auto& edge : edges) {
        if (intersectSegments(segment, edge, tol).kind != SegmentIntersectionKind::None) {
            return true;
        }
    }
    return false;
}

CellClass classifyCartesianCell(const CartesianCell2D& cell,
                                const BoundaryLoop& boundary,
                                const TolerancePolicy& tol) noexcept {
    const auto& vertices = boundary.vertices();
    if (vertices.size() < 3) return CellClass::Outside;

    for (std::size_t k = 0; k < vertices.size(); ++k) {
        const Segment2D edge{vertices[k], vertices[(k + 1) % vertices.size()]};
        if (segmentIntersectsClosedAABB(edge, cell.bounds, tol)) {
            return CellClass::Intersected;
        }
    }

    const auto pointClass = classifyPointInPolygon(cell.center(), boundary.polygon(), tol);
    if (pointClass == PointInPolygon::Inside) return CellClass::Inside;
    if (pointClass == PointInPolygon::Boundary) return CellClass::Intersected;
    return CellClass::Outside;
}

ClassificationSummary classifyGrid(UniformCartesianGrid2D& grid,
                                   const BoundaryLoop& boundary,
                                   const TolerancePolicy& tol) {
    if (!boundary.diagnose(tol).valid()) {
        throw std::invalid_argument("cannot classify grid against invalid boundary loop");
    }

    ClassificationSummary summary;
    for (auto& cell : grid.cells()) {
        cell.classification = classifyCartesianCell(cell, boundary, tol);
        switch (cell.classification) {
        case CellClass::Inside: ++summary.inside; break;
        case CellClass::Outside: ++summary.outside; break;
        case CellClass::Intersected: ++summary.intersected; break;
        }
    }
    return summary;
}

} // namespace cartmesh2d
