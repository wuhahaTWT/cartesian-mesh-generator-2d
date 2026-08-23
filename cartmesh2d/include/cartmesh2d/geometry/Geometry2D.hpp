#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

struct TolerancePolicy {
    double absolute = 1.0e-12;
    double relative = 1.0e-10;

    [[nodiscard]] double scale(double magnitude = 1.0) const noexcept;
    [[nodiscard]] bool nearlyEqual(double a, double b, double magnitude = 1.0) const noexcept;
};

struct Vector2D {
    double x = 0.0;
    double y = 0.0;
};

struct Point2D {
    double x = 0.0;
    double y = 0.0;
};

[[nodiscard]] Vector2D operator-(const Point2D& a, const Point2D& b) noexcept;
[[nodiscard]] Point2D operator+(const Point2D& p, const Vector2D& v) noexcept;
[[nodiscard]] Vector2D operator*(const Vector2D& v, double s) noexcept;
[[nodiscard]] double dot(const Vector2D& a, const Vector2D& b) noexcept;
[[nodiscard]] double cross(const Vector2D& a, const Vector2D& b) noexcept;
[[nodiscard]] double squaredNorm(const Vector2D& v) noexcept;
[[nodiscard]] double orientation(const Point2D& a, const Point2D& b, const Point2D& c) noexcept;
// Adaptive exact sign for the orientation determinant of binary64 inputs.
// The common path uses a floating-point error bound; uncertain cases fall
// back to exact signed-integer arithmetic over the input dyadics.
[[nodiscard]] int orientationSign(const Point2D& a, const Point2D& b,
                                  const Point2D& c) noexcept;
[[nodiscard]] bool nearlyEqual(const Point2D& a, const Point2D& b,
                               const TolerancePolicy& tol = {}) noexcept;

struct Segment2D {
    Point2D a;
    Point2D b;
};

struct AABB2D {
    Point2D min;
    Point2D max;

    [[nodiscard]] bool valid(const TolerancePolicy& tol = {}) const noexcept;
    [[nodiscard]] bool contains(const Point2D& p, const TolerancePolicy& tol = {}) const noexcept;
};

struct Polygon2D {
    std::vector<Point2D> vertices;

    [[nodiscard]] double signedArea() const noexcept;
    [[nodiscard]] double area() const noexcept;
    [[nodiscard]] std::optional<Point2D> centroid(const TolerancePolicy& tol = {}) const noexcept;
    [[nodiscard]] AABB2D bounds() const noexcept;
};

enum class SegmentIntersectionKind { None, Point, Overlap };

struct SegmentIntersection {
    SegmentIntersectionKind kind = SegmentIntersectionKind::None;
    std::optional<Point2D> point;
    std::optional<Segment2D> overlap;
};

[[nodiscard]] bool pointOnSegment(const Point2D& p, const Segment2D& segment,
                                  const TolerancePolicy& tol = {}) noexcept;
[[nodiscard]] SegmentIntersection intersectSegments(const Segment2D& lhs,
                                                    const Segment2D& rhs,
                                                    const TolerancePolicy& tol = {}) noexcept;

enum class PointInPolygon { Outside, Inside, Boundary };
[[nodiscard]] PointInPolygon classifyPointInPolygon(const Point2D& p,
                                                    const Polygon2D& polygon,
                                                    const TolerancePolicy& tol = {}) noexcept;

enum class LoopOrientation { Degenerate, Clockwise, CounterClockwise };

enum class BoundaryIssueCode {
    NonFiniteCoordinate,
    TooFewUniqueVertices,
    DuplicateConsecutiveVertex,
    ZeroLengthEdge,
    SelfIntersection,
    ZeroArea
};

struct BoundaryIssue {
    BoundaryIssueCode code;
    std::size_t firstIndex = 0;
    std::size_t secondIndex = 0;
    std::string message;
};

struct BoundaryDiagnostics {
    std::vector<BoundaryIssue> issues;
    LoopOrientation orientation = LoopOrientation::Degenerate;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && orientation != LoopOrientation::Degenerate;
    }
};

class BoundaryLoop {
public:
    explicit BoundaryLoop(std::vector<Point2D> vertices);

    [[nodiscard]] const std::vector<Point2D>& vertices() const noexcept { return vertices_; }
    [[nodiscard]] Polygon2D polygon() const { return Polygon2D{vertices_}; }
    [[nodiscard]] BoundaryDiagnostics diagnose(const TolerancePolicy& tol = {}) const;
    [[nodiscard]] bool normalizeCounterClockwise(const TolerancePolicy& tol = {});
    [[nodiscard]] bool normalizeClockwise(const TolerancePolicy& tol = {});

private:
    std::vector<Point2D> vertices_;
};

enum class BoundaryRegionIssueCode {
    EmptyRegion,
    InvalidLoop,
    IntersectingLoops
};

struct BoundaryRegionIssue2D {
    BoundaryRegionIssueCode code;
    std::size_t loopA = 0;
    std::size_t loopB = 0;
    std::string message;
};

struct BoundaryRegionDiagnostics2D {
    std::vector<BoundaryRegionIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

// A deterministic even-odd region bounded by one or more disjoint or nested
// simple loops.  Nesting alternates material/hole semantics: depth 0 is region,
// depth 1 is a hole, depth 2 is an island, and so on.  Intersecting or touching
// loops are rejected because their topology is ambiguous without arrangement
// repair.
class BoundaryRegion2D {
public:
    explicit BoundaryRegion2D(std::vector<BoundaryLoop> loops);
    explicit BoundaryRegion2D(const BoundaryLoop& loop);

    [[nodiscard]] const std::vector<BoundaryLoop>& loops() const noexcept { return loops_; }
    [[nodiscard]] BoundaryRegionDiagnostics2D diagnose(
        const TolerancePolicy& tol = {}) const;
    [[nodiscard]] bool normalizeAlternating(const TolerancePolicy& tol = {});
    [[nodiscard]] PointInPolygon classifyPoint(
        const Point2D& point, const TolerancePolicy& tol = {}) const;
    [[nodiscard]] std::vector<std::size_t> nestingDepths(
        const TolerancePolicy& tol = {}) const;
    [[nodiscard]] AABB2D bounds() const noexcept;
    [[nodiscard]] double area(const TolerancePolicy& tol = {}) const;

private:
    std::vector<BoundaryLoop> loops_;
    mutable std::optional<BoundaryRegionDiagnostics2D> diagnosticsCache_;
    mutable std::optional<std::vector<std::size_t>> nestingDepthCache_;
    mutable std::optional<TolerancePolicy> cacheTolerance_;
};

} // namespace cartmesh2d
