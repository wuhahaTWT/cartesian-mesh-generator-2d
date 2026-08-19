#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <utility>

namespace cartmesh2d {

namespace {

double maxMagnitude(std::initializer_list<double> values) noexcept {
    double result = 1.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

bool segmentsAreAdjacent(std::size_t i, std::size_t j, std::size_t edgeCount) noexcept {
    if (i == j) return true;
    if (i + 1 == j || j + 1 == i) return true;
    return (i == 0 && j + 1 == edgeCount) || (j == 0 && i + 1 == edgeCount);
}

} // namespace

double TolerancePolicy::scale(double magnitude) const noexcept {
    return absolute + relative * std::max(1.0, std::abs(magnitude));
}

bool TolerancePolicy::nearlyEqual(double a, double b, double magnitude) const noexcept {
    return std::abs(a - b) <= scale(std::max({magnitude, std::abs(a), std::abs(b)}));
}

Vector2D operator-(const Point2D& a, const Point2D& b) noexcept {
    return {a.x - b.x, a.y - b.y};
}

Point2D operator+(const Point2D& p, const Vector2D& v) noexcept {
    return {p.x + v.x, p.y + v.y};
}

Vector2D operator*(const Vector2D& v, double s) noexcept {
    return {v.x * s, v.y * s};
}

double dot(const Vector2D& a, const Vector2D& b) noexcept {
    return a.x * b.x + a.y * b.y;
}

double cross(const Vector2D& a, const Vector2D& b) noexcept {
    return a.x * b.y - a.y * b.x;
}

double squaredNorm(const Vector2D& v) noexcept {
    return dot(v, v);
}

double orientation(const Point2D& a, const Point2D& b, const Point2D& c) noexcept {
    return cross(b - a, c - a);
}

bool nearlyEqual(const Point2D& a, const Point2D& b, const TolerancePolicy& tol) noexcept {
    const double magnitude = maxMagnitude({a.x, a.y, b.x, b.y});
    return tol.nearlyEqual(a.x, b.x, magnitude) && tol.nearlyEqual(a.y, b.y, magnitude);
}

bool AABB2D::valid(const TolerancePolicy& tol) const noexcept {
    return min.x <= max.x + tol.scale(maxMagnitude({min.x, max.x})) &&
           min.y <= max.y + tol.scale(maxMagnitude({min.y, max.y}));
}

bool AABB2D::contains(const Point2D& p, const TolerancePolicy& tol) const noexcept {
    const double sx = tol.scale(maxMagnitude({min.x, max.x, p.x}));
    const double sy = tol.scale(maxMagnitude({min.y, max.y, p.y}));
    return p.x >= min.x - sx && p.x <= max.x + sx &&
           p.y >= min.y - sy && p.y <= max.y + sy;
}

double Polygon2D::signedArea() const noexcept {
    if (vertices.size() < 3) return 0.0;
    double twiceArea = 0.0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % vertices.size()];
        twiceArea += a.x * b.y - b.x * a.y;
    }
    return 0.5 * twiceArea;
}

double Polygon2D::area() const noexcept {
    return std::abs(signedArea());
}

std::optional<Point2D> Polygon2D::centroid(const TolerancePolicy& tol) const noexcept {
    if (vertices.size() < 3) return std::nullopt;

    double twiceArea = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % vertices.size()];
        const double term = a.x * b.y - b.x * a.y;
        twiceArea += term;
        cx += (a.x + b.x) * term;
        cy += (a.y + b.y) * term;
    }

    if (std::abs(twiceArea) <= tol.scale(std::max(1.0, area()))) return std::nullopt;
    return Point2D{cx / (3.0 * twiceArea), cy / (3.0 * twiceArea)};
}

AABB2D Polygon2D::bounds() const noexcept {
    if (vertices.empty()) return {{0.0, 0.0}, {0.0, 0.0}};
    Point2D lo = vertices.front();
    Point2D hi = vertices.front();
    for (const auto& p : vertices) {
        lo.x = std::min(lo.x, p.x);
        lo.y = std::min(lo.y, p.y);
        hi.x = std::max(hi.x, p.x);
        hi.y = std::max(hi.y, p.y);
    }
    return {lo, hi};
}

bool pointOnSegment(const Point2D& p, const Segment2D& segment,
                    const TolerancePolicy& tol) noexcept {
    const Vector2D ab = segment.b - segment.a;
    const Vector2D ap = p - segment.a;
    const double magnitude = maxMagnitude({segment.a.x, segment.a.y,
                                           segment.b.x, segment.b.y, p.x, p.y});
    const double eps = tol.scale(magnitude);
    if (squaredNorm(ab) <= eps * eps) return nearlyEqual(p, segment.a, tol);

    if (std::abs(cross(ab, ap)) > eps * std::max(1.0, std::sqrt(squaredNorm(ab)))) {
        return false;
    }
    const double projection = dot(ap, ab);
    return projection >= -eps && projection <= squaredNorm(ab) + eps;
}

SegmentIntersection intersectSegments(const Segment2D& lhs, const Segment2D& rhs,
                                      const TolerancePolicy& tol) noexcept {
    const Vector2D r = lhs.b - lhs.a;
    const Vector2D s = rhs.b - rhs.a;
    const double rr = squaredNorm(r);
    const double ss = squaredNorm(s);
    const double magnitude = maxMagnitude({lhs.a.x, lhs.a.y, lhs.b.x, lhs.b.y,
                                           rhs.a.x, rhs.a.y, rhs.b.x, rhs.b.y});
    const double eps = tol.scale(magnitude);

    if (rr <= eps * eps && ss <= eps * eps) {
        if (nearlyEqual(lhs.a, rhs.a, tol)) {
            return {SegmentIntersectionKind::Point, lhs.a, std::nullopt};
        }
        return {};
    }
    if (rr <= eps * eps) {
        if (pointOnSegment(lhs.a, rhs, tol)) {
            return {SegmentIntersectionKind::Point, lhs.a, std::nullopt};
        }
        return {};
    }
    if (ss <= eps * eps) {
        if (pointOnSegment(rhs.a, lhs, tol)) {
            return {SegmentIntersectionKind::Point, rhs.a, std::nullopt};
        }
        return {};
    }

    const Vector2D qp = rhs.a - lhs.a;
    const double rxs = cross(r, s);
    const double qpxr = cross(qp, r);
    const double crossEps = eps * std::max({1.0, std::sqrt(rr), std::sqrt(ss)});

    if (std::abs(rxs) <= crossEps) {
        if (std::abs(qpxr) > crossEps) return {};

        const double t0 = dot(rhs.a - lhs.a, r) / rr;
        const double t1 = dot(rhs.b - lhs.a, r) / rr;
        const double lo = std::max(0.0, std::min(t0, t1));
        const double hi = std::min(1.0, std::max(t0, t1));
        const double tEps = eps / std::max(1.0, std::sqrt(rr));
        if (hi < lo - tEps) return {};
        if (std::abs(hi - lo) <= tEps) {
            return {SegmentIntersectionKind::Point,
                    lhs.a + r * ((lo + hi) * 0.5), std::nullopt};
        }
        return {SegmentIntersectionKind::Overlap, std::nullopt,
                Segment2D{lhs.a + r * lo, lhs.a + r * hi}};
    }

    const double t = cross(qp, s) / rxs;
    const double u = cross(qp, r) / rxs;
    const double paramEps = eps / std::max({1.0, std::sqrt(rr), std::sqrt(ss)});
    if (t < -paramEps || t > 1.0 + paramEps ||
        u < -paramEps || u > 1.0 + paramEps) {
        return {};
    }

    return {SegmentIntersectionKind::Point,
            lhs.a + r * std::clamp(t, 0.0, 1.0), std::nullopt};
}

PointInPolygon classifyPointInPolygon(const Point2D& p, const Polygon2D& polygon,
                                      const TolerancePolicy& tol) noexcept {
    const auto& vertices = polygon.vertices;
    if (vertices.size() < 3) return PointInPolygon::Outside;

    int winding = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Point2D& a = vertices[i];
        const Point2D& b = vertices[(i + 1) % vertices.size()];
        if (pointOnSegment(p, {a, b}, tol)) return PointInPolygon::Boundary;

        const double eps = tol.scale(maxMagnitude({a.x, a.y, b.x, b.y, p.x, p.y}));
        if (a.y <= p.y) {
            if (b.y > p.y && orientation(a, b, p) > eps) ++winding;
        } else if (b.y <= p.y && orientation(a, b, p) < -eps) {
            --winding;
        }
    }
    return winding == 0 ? PointInPolygon::Outside : PointInPolygon::Inside;
}

BoundaryLoop::BoundaryLoop(std::vector<Point2D> vertices) : vertices_(std::move(vertices)) {
    if (vertices_.size() > 1 &&
        vertices_.front().x == vertices_.back().x &&
        vertices_.front().y == vertices_.back().y) {
        vertices_.pop_back();
    }
}

BoundaryDiagnostics BoundaryLoop::diagnose(const TolerancePolicy& tol) const {
    BoundaryDiagnostics result;
    const std::size_t n = vertices_.size();

    std::vector<Point2D> unique;
    for (const auto& p : vertices_) {
        bool found = false;
        for (const auto& q : unique) {
            if (nearlyEqual(p, q, tol)) {
                found = true;
                break;
            }
        }
        if (!found) unique.push_back(p);
    }
    if (unique.size() < 3) {
        result.issues.push_back({BoundaryIssueCode::TooFewUniqueVertices, 0, 0,
                                 "boundary loop has fewer than three unique vertices"});
    }

    if (n >= 2) {
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = (i + 1) % n;
            if (nearlyEqual(vertices_[i], vertices_[j], tol)) {
                result.issues.push_back({BoundaryIssueCode::DuplicateConsecutiveVertex, i, j,
                                         "duplicate consecutive boundary vertex"});
                result.issues.push_back({BoundaryIssueCode::ZeroLengthEdge, i, j,
                                         "zero-length boundary edge"});
            }
        }
    }

    if (n >= 4) {
        for (std::size_t i = 0; i < n; ++i) {
            const Segment2D a{vertices_[i], vertices_[(i + 1) % n]};
            for (std::size_t j = i + 1; j < n; ++j) {
                if (segmentsAreAdjacent(i, j, n)) continue;
                const Segment2D b{vertices_[j], vertices_[(j + 1) % n]};
                if (intersectSegments(a, b, tol).kind != SegmentIntersectionKind::None) {
                    result.issues.push_back({BoundaryIssueCode::SelfIntersection, i, j,
                                             "non-adjacent boundary edges intersect"});
                }
            }
        }
    }

    const double signedArea = polygon().signedArea();
    const auto box = polygon().bounds();
    const double scale = std::max({1.0, box.max.x - box.min.x, box.max.y - box.min.y});
    if (std::abs(signedArea) <= tol.scale(scale * scale)) {
        result.orientation = LoopOrientation::Degenerate;
        result.issues.push_back({BoundaryIssueCode::ZeroArea, 0, 0,
                                 "boundary loop has zero or numerically degenerate signed area"});
    } else {
        result.orientation = signedArea > 0.0 ? LoopOrientation::CounterClockwise
                                               : LoopOrientation::Clockwise;
    }
    return result;
}

bool BoundaryLoop::normalizeCounterClockwise(const TolerancePolicy& tol) {
    const auto diagnostics = diagnose(tol);
    if (!diagnostics.valid()) return false;
    if (diagnostics.orientation == LoopOrientation::Clockwise) {
        std::reverse(vertices_.begin(), vertices_.end());
    }
    return true;
}

} // namespace cartmesh2d
