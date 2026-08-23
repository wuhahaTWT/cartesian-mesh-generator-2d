#include "cartmesh2d/grid/CartesianGrid2D.hpp"
#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include "cartmesh2d/spatial/BoundarySegmentIndex2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace cartmesh2d;

namespace {
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

CellClass bruteClassify(const CartesianCell2D& cell,
                        const BoundaryLoop& boundary,
                        const TolerancePolicy& tol) {
    const auto& vertices = boundary.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (segmentIntersectsClosedAABB(
                {vertices[i], vertices[(i + 1U) % vertices.size()]},
                cell.bounds, tol)) return CellClass::Intersected;
    }
    const auto state = classifyPointInPolygon(cell.center(), boundary.polygon(), tol);
    return state == PointInPolygon::Inside ? CellClass::Inside
         : state == PointInPolygon::Boundary ? CellClass::Intersected
         : CellClass::Outside;
}

} // namespace

int main() {
    constexpr double pi = 3.14159265358979323846;
    std::vector<Point2D> points;
    constexpr std::size_t segmentCount = 256U;
    points.reserve(segmentCount);
    for (std::size_t i = 0; i < segmentCount; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) /
                             static_cast<double>(segmentCount);
        const double radius = 1.0 + 0.08 * std::cos(7.0 * angle);
        points.push_back({radius * std::cos(angle), radius * std::sin(angle)});
    }
    const BoundaryLoop boundary(points);
    const TolerancePolicy tol{};
    const BoundarySegmentIndex2D index(boundary, tol);
    check(index.valid(), "BVH accepts a valid closed boundary");
    check(index.segmentCount() == segmentCount, "BVH preserves every source segment");
    check(index.matches(boundary), "BVH matches its exact source loop");

    auto changedPoints = points;
    changedPoints.front().x += 1.0e-6;
    check(!index.matches(BoundaryLoop(changedPoints)),
          "BVH refuses a different boundary instead of using stale geometry");

    const Domain2D domain{{{-1.5, -1.5}, {1.5, 1.5}}};
    UniformCartesianGrid2D grid(domain, 48U, 48U);
    for (const auto& cell : grid.cells()) {
        const CellClass expected = bruteClassify(cell, boundary, tol);
        CellClass actual = CellClass::Outside;
        if (index.intersects(cell.bounds, tol)) {
            actual = CellClass::Intersected;
        } else {
            const auto state = index.classifyPoint(cell.center(), tol);
            actual = state == PointInPolygon::Inside ? CellClass::Inside
                   : state == PointInPolygon::Boundary ? CellClass::Intersected
                   : CellClass::Outside;
        }
        check(actual == expected, "BVH cell classification matches exhaustive scan");
    }

    for (int j = -30; j <= 30; ++j) {
        for (int i = -30; i <= 30; ++i) {
            const Point2D point{0.05 * static_cast<double>(i),
                                0.05 * static_cast<double>(j)};
            check(index.classifyPoint(point, tol) ==
                      classifyPointInPolygon(point, boundary.polygon(), tol),
                  "BVH winding classification matches exhaustive polygon test");
        }
    }

    for (int j = -8; j <= 8; ++j) {
        for (int i = -8; i <= 8; ++i) {
            const Point2D center{0.18 * static_cast<double>(i),
                                 0.18 * static_cast<double>(j)};
            const AABB2D box{{center.x - 0.025, center.y - 0.04},
                             {center.x + 0.025, center.y + 0.04}};
            const double expected = distanceAABBToBoundary(box, boundary, tol);
            const double actual = index.distanceToAABB(box, tol);
            check(std::abs(actual - expected) <= 1.0e-12,
                  "BVH nearest distance matches exhaustive scan");
        }
    }

    const AABB2D probe{{0.9, -0.2}, {1.15, 0.2}};
    const auto first = index.querySegmentIds(probe, tol);
    const auto second = index.querySegmentIds(probe, tol);
    check(!first.empty(), "BVH overlap query returns local candidates");
    check(first == second, "BVH query ordering is deterministic");

    const BoundaryLoop invalid({{0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {1.0, 0.0}});
    check(!BoundarySegmentIndex2D(invalid, tol).valid(),
          "BVH refuses a self-intersecting boundary");

    if (failures != 0) {
        std::cerr << failures << " spatial index test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d P1 boundary segment BVH tests passed\n";
    return EXIT_SUCCESS;
}
