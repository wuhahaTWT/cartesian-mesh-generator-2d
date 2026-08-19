#include "cartmesh2d/geometry/Geometry2D.hpp"

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

void near(double actual, double expected, double eps, const char* message) {
    check(std::abs(actual - expected) <= eps, message);
}
} // namespace

int main() {
    const TolerancePolicy tol{};

    Polygon2D rectangle{{{0, 0}, {4, 0}, {4, 2}, {0, 2}}};
    near(rectangle.signedArea(), 8.0, 1e-12, "rectangle signed area");
    near(rectangle.area(), 8.0, 1e-12, "rectangle area");
    auto rectangleCentroid = rectangle.centroid(tol);
    check(rectangleCentroid.has_value(), "rectangle centroid exists");
    near(rectangleCentroid->x, 2.0, 1e-12, "rectangle centroid x");
    near(rectangleCentroid->y, 1.0, 1e-12, "rectangle centroid y");

    Polygon2D triangle{{{0, 0}, {4, 0}, {0, 2}}};
    near(triangle.area(), 4.0, 1e-12, "triangle area");

    Polygon2D concaveL{{{0, 0}, {3, 0}, {3, 1}, {1, 1}, {1, 3}, {0, 3}}};
    near(concaveL.area(), 5.0, 1e-12, "concave L area");
    check(classifyPointInPolygon({0.5, 2.5}, concaveL, tol) == PointInPolygon::Inside,
          "L inside");
    check(classifyPointInPolygon({2.0, 2.0}, concaveL, tol) == PointInPolygon::Outside,
          "L notch outside");
    check(classifyPointInPolygon({1.0, 2.0}, concaveL, tol) == PointInPolygon::Boundary,
          "L boundary");

    std::vector<Point2D> circlePoints;
    constexpr int count = 64;
    constexpr double pi = 3.14159265358979323846;
    for (int i = 0; i < count; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / count;
        circlePoints.push_back({std::cos(angle), std::sin(angle)});
    }
    Polygon2D circle{circlePoints};
    check(std::abs(circle.area() - pi) < 0.01, "discretized circle area");
    check(classifyPointInPolygon({0, 0}, circle, tol) == PointInPolygon::Inside,
          "circle center inside");

    auto crossing = intersectSegments({{0, 0}, {2, 2}}, {{0, 2}, {2, 0}}, tol);
    check(crossing.kind == SegmentIntersectionKind::Point, "crossing segment kind");
    near(crossing.point->x, 1.0, 1e-12, "crossing x");
    near(crossing.point->y, 1.0, 1e-12, "crossing y");

    auto endpoint = intersectSegments({{0, 0}, {1, 0}}, {{1, 0}, {1, 2}}, tol);
    check(endpoint.kind == SegmentIntersectionKind::Point, "endpoint touch");
    check(nearlyEqual(*endpoint.point, {1, 0}, tol), "endpoint point");

    auto parallel = intersectSegments({{0, 0}, {1, 0}}, {{0, 1}, {1, 1}}, tol);
    check(parallel.kind == SegmentIntersectionKind::None, "parallel non-intersection");

    auto overlap = intersectSegments({{0, 0}, {3, 0}}, {{1, 0}, {4, 0}}, tol);
    check(overlap.kind == SegmentIntersectionKind::Overlap, "collinear overlap kind");
    check(nearlyEqual(overlap.overlap->a, {1, 0}, tol), "overlap start");
    check(nearlyEqual(overlap.overlap->b, {3, 0}, tol), "overlap end");

    BoundaryLoop validLoop({{0, 0}, {0, 2}, {2, 2}, {2, 0}});
    auto validDiagnostics = validLoop.diagnose(tol);
    check(validDiagnostics.valid(), "clockwise rectangle valid before normalization");
    check(validDiagnostics.orientation == LoopOrientation::Clockwise,
          "clockwise orientation detected");
    check(validLoop.normalizeCounterClockwise(tol), "normalize valid loop");
    check(validLoop.diagnose(tol).orientation == LoopOrientation::CounterClockwise,
          "normalized CCW");

    BoundaryLoop bowTie({{0, 0}, {2, 2}, {0, 2}, {2, 0}});
    auto bowTieDiagnostics = bowTie.diagnose(tol);
    check(!bowTieDiagnostics.valid(), "bow tie invalid");
    bool selfIntersectionFound = false;
    for (const auto& issue : bowTieDiagnostics.issues) {
        selfIntersectionFound |= issue.code == BoundaryIssueCode::SelfIntersection;
    }
    check(selfIntersectionFound, "bow tie self intersection diagnosed");
    check(!bowTie.normalizeCounterClockwise(tol), "bow tie normalization rejected");

    BoundaryLoop duplicate({{0, 0}, {1, 0}, {1, 0}, {0, 1}});
    auto duplicateDiagnostics = duplicate.diagnose(tol);
    check(!duplicateDiagnostics.valid(), "duplicate consecutive invalid");
    bool duplicateFound = false;
    bool zeroLengthFound = false;
    for (const auto& issue : duplicateDiagnostics.issues) {
        duplicateFound |= issue.code == BoundaryIssueCode::DuplicateConsecutiveVertex;
        zeroLengthFound |= issue.code == BoundaryIssueCode::ZeroLengthEdge;
    }
    check(duplicateFound, "duplicate diagnostic");
    check(zeroLengthFound, "zero edge diagnostic");

    if (failures == 0) {
        std::cout << "cartmesh2d Stage 2D-0 geometry tests: PASS\n";
        return EXIT_SUCCESS;
    }

    std::cerr << failures << " test(s) failed\n";
    return EXIT_FAILURE;
}
