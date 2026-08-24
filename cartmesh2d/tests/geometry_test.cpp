#include "cartmesh2d/geometry/Geometry2D.hpp"
#include "cartmesh2d/geometry/BoundarySimplification2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
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

    // The two products in the ordinary binary64 determinant both round to
    // 2^54, although the exact determinant of the input doubles is -1.
    const Point2D cancellationA{0.0, 0.0};
    const Point2D cancellationB{134217729.0, 134217728.0};
    const Point2D cancellationC{134217728.0, 134217727.0};
    check(orientationSign(cancellationA, cancellationB, cancellationC) == -1,
          "adaptive orientation recovers cancellation-hidden negative sign");
    check(orientationSign(cancellationA, cancellationC, cancellationB) == 1,
          "adaptive orientation recovers cancellation-hidden positive sign");
    check(orientationSign({0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}) == 0,
          "adaptive orientation identifies exact collinearity");

    const double subnormal = std::numeric_limits<double>::denorm_min();
    check(orientationSign({0.0, 0.0}, {subnormal, 0.0}, {0.0, subnormal}) == 1,
          "exact fallback survives determinant underflow for subnormal inputs");
    check(orientationSign({0.0, 0.0},
                          {std::numeric_limits<double>::infinity(), 0.0},
                          {0.0, 1.0}) == 0,
          "orientation rejects non-finite input deterministically");

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

    BoundaryLoop nonFinite({{0.0, 0.0},
                            {1.0, 0.0},
                            {std::numeric_limits<double>::quiet_NaN(), 1.0}});
    const auto nonFiniteDiagnostics = nonFinite.diagnose(tol);
    check(!nonFiniteDiagnostics.valid(), "non-finite boundary is rejected");
    check(!nonFiniteDiagnostics.issues.empty() &&
              nonFiniteDiagnostics.issues.front().code == BoundaryIssueCode::NonFiniteCoordinate,
          "non-finite boundary reports its explicit issue code");
    check(!AABB2D{{0.0, 0.0},
                  {std::numeric_limits<double>::infinity(), 1.0}}.valid(tol),
          "non-finite Cartesian box is rejected");

    BoundaryRegion2D nestedRegion({
        BoundaryLoop({{0.0,0.0},{4.0,0.0},{4.0,4.0},{0.0,4.0}}),
        BoundaryLoop({{1.0,1.0},{1.0,3.0},{3.0,3.0},{3.0,1.0}}),
        BoundaryLoop({{1.5,1.5},{2.5,1.5},{2.5,2.5},{1.5,2.5}})
    });
    check(nestedRegion.diagnose(tol).valid(), "three-level nested boundary region is valid");
    check(nestedRegion.classifyPoint({0.5,0.5},tol)==PointInPolygon::Inside,
          "outer shell belongs to even-odd region");
    check(nestedRegion.classifyPoint({1.25,1.25},tol)==PointInPolygon::Outside,
          "nested depth-one loop is a hole");
    check(nestedRegion.classifyPoint({2.0,2.0},tol)==PointInPolygon::Inside,
          "nested depth-two loop is an island");
    near(nestedRegion.area(tol),13.0,1e-12,"nested even-odd region area");
    check(nestedRegion.normalizeAlternating(tol),"nested loops normalize by parity depth");
    const auto normalizedDepths=nestedRegion.nestingDepths(tol);
    check(normalizedDepths==std::vector<std::size_t>({0,1,2}),
          "nesting depth is deterministic and orientation independent");
    check(nestedRegion.loops()[0].diagnose(tol).orientation==LoopOrientation::CounterClockwise &&
          nestedRegion.loops()[1].diagnose(tol).orientation==LoopOrientation::Clockwise &&
          nestedRegion.loops()[2].diagnose(tol).orientation==LoopOrientation::CounterClockwise,
          "alternating normalization orients region material on the left");

    BoundaryRegion2D disconnectedRegion({
        BoundaryLoop({{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}}),
        BoundaryLoop({{2.0,0.0},{3.0,0.0},{3.0,1.0},{2.0,1.0}})
    });
    check(disconnectedRegion.diagnose(tol).valid() &&
          disconnectedRegion.classifyPoint({0.5,0.5},tol)==PointInPolygon::Inside &&
          disconnectedRegion.classifyPoint({2.5,0.5},tol)==PointInPolygon::Inside,
          "two disjoint loops form two region components");
    near(disconnectedRegion.area(tol),2.0,1e-12,"disconnected region area sums components");

    BoundaryRegion2D touchingRegion({
        BoundaryLoop({{0.0,0.0},{2.0,0.0},{2.0,2.0},{0.0,2.0}}),
        BoundaryLoop({{2.0,0.5},{3.0,0.5},{3.0,1.5},{2.0,1.5}})
    });
    check(!touchingRegion.diagnose(tol).valid(),
          "touching loops are rejected instead of assigned ambiguous nesting");
    BoundaryRegion2D emptyRegion(std::vector<BoundaryLoop>{});
    check(!emptyRegion.diagnose(tol).valid() && emptyRegion.nestingDepths(tol).empty(),
          "empty boundary region is rejected without unsafe nesting access");
    BoundaryRegion2D toleranceRegion({
        BoundaryLoop({{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}}),
        BoundaryLoop({{1.0001,0.0},{2.0001,0.0},{2.0001,1.0},{1.0001,1.0}})
    });
    TolerancePolicy looseTol=tol;
    looseTol.absolute=1.0e-3;
    check(toleranceRegion.diagnose(tol).valid() && !toleranceRegion.diagnose(looseTol).valid(),
          "boundary-region diagnostic cache is keyed by tolerance policy");

    BoundaryRegion2D denseCircle{BoundaryLoop(circlePoints)};
    const auto simplifiedCircle=simplifyBoundaryRegion2D(denseCircle,0.01,tol);
    check(simplifiedCircle.valid(),"dense circle simplification preserves valid topology");
    check(simplifiedCircle.report.simplifiedVertexCount<
              simplifiedCircle.report.originalVertexCount,
          "dense circle simplification removes sub-grid vertices");
    check(simplifiedCircle.report.measuredMaxDeviation<=0.01+1.0e-12,
          "simplification reports and respects its deviation bound");
    check(simplifiedCircle.boundary->diagnose(tol).valid(),
          "simplified circle remains independently diagnosable");

    std::rotate(circlePoints.begin(),circlePoints.begin()+17,circlePoints.end());
    const auto rotatedCircle=simplifyBoundaryRegion2D(
        BoundaryRegion2D(BoundaryLoop(circlePoints)),0.01,tol);
    const auto& canonicalCircle=simplifiedCircle.boundary->loops()[0].vertices();
    const auto& canonicalRotated=rotatedCircle.boundary->loops()[0].vertices();
    check(rotatedCircle.valid() && canonicalCircle.size()==canonicalRotated.size() &&
              std::equal(canonicalCircle.begin(),canonicalCircle.end(),
                         canonicalRotated.begin(),[](const Point2D& a,const Point2D& b) {
                             return a.x==b.x && a.y==b.y;
                         }),
          "closed-loop simplification is invariant to cyclic input rotation");

    const auto unchangedRectangle=simplifyBoundaryRegion2D(
        BoundaryRegion2D(BoundaryLoop({{0.0,0.0},{4.0,0.0},{4.0,2.0},{0.0,2.0}})),
        0.1,tol);
    check(unchangedRectangle.valid() &&
              unchangedRectangle.report.originalVertexCount==4 &&
              unchangedRectangle.report.simplifiedVertexCount==4 &&
              unchangedRectangle.report.measuredMaxDeviation==0.0,
          "feature corners are preserved exactly");

    const auto invalidSimplification=simplifyBoundaryRegion2D(emptyRegion,0.1,tol);
    check(!invalidSimplification.valid(),"invalid source simplification fails closed");

    BoundaryRegion2D notchAndIsland({
        BoundaryLoop({{0.0,0.0},{4.0,0.0},{4.0,4.0},{2.5,4.0},
                      {2.5,1.0},{1.5,1.0},{1.5,4.0},{0.0,4.0}}),
        BoundaryLoop({{1.8,2.0},{2.2,2.0},{2.2,3.0},{1.8,3.0}})
    });
    const auto topologyPreserving=simplifyBoundaryRegion2D(notchAndIsland,10.0,tol);
    check(topologyPreserving.valid() &&
              topologyPreserving.boundary->nestingDepths(tol)==
              notchAndIsland.nestingDepths(tol),
          "simplification preserves multi-loop nesting depths");

    if (failures == 0) {
        std::cout << "cartmesh2d Stage 2D-0 geometry tests: PASS\n";
        return EXIT_SUCCESS;
    }

    std::cerr << failures << " test(s) failed\n";
    return EXIT_FAILURE;
}
