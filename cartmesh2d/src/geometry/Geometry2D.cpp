#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

namespace cartmesh2d {

namespace {

constexpr std::size_t kCoordinateWords = 36;
constexpr std::size_t kProductWords = 72;

template <std::size_t N>
using BigMagnitude = std::array<std::uint64_t, N>;

template <std::size_t N>
struct SignedBig {
    int sign = 0;
    BigMagnitude<N> magnitude{};
};

template <std::size_t N>
[[nodiscard]] int compareMagnitude(const BigMagnitude<N>& lhs,
                                   const BigMagnitude<N>& rhs) noexcept {
    for (std::size_t i = N; i-- > 0;) {
        if (lhs[i] < rhs[i]) return -1;
        if (lhs[i] > rhs[i]) return 1;
    }
    return 0;
}

template <std::size_t N>
[[nodiscard]] BigMagnitude<N> addMagnitude(const BigMagnitude<N>& lhs,
                                           const BigMagnitude<N>& rhs) noexcept {
    BigMagnitude<N> result{};
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const unsigned __int128 sum = static_cast<unsigned __int128>(lhs[i]) +
                                      rhs[i] + carry;
        result[i] = static_cast<std::uint64_t>(sum);
        carry = static_cast<std::uint64_t>(sum >> 64U);
    }
    return result;
}

template <std::size_t N>
[[nodiscard]] BigMagnitude<N> subtractMagnitude(const BigMagnitude<N>& lhs,
                                                const BigMagnitude<N>& rhs) noexcept {
    BigMagnitude<N> result{};
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const unsigned __int128 subtrahend =
            static_cast<unsigned __int128>(rhs[i]) + borrow;
        result[i] = static_cast<std::uint64_t>(
            static_cast<unsigned __int128>(lhs[i]) - subtrahend);
        borrow = static_cast<unsigned __int128>(lhs[i]) < subtrahend ? 1U : 0U;
    }
    return result;
}

template <std::size_t N>
[[nodiscard]] SignedBig<N> addSigned(const SignedBig<N>& lhs,
                                     const SignedBig<N>& rhs) noexcept {
    if (lhs.sign == 0) return rhs;
    if (rhs.sign == 0) return lhs;
    if (lhs.sign == rhs.sign) {
        return {lhs.sign, addMagnitude(lhs.magnitude, rhs.magnitude)};
    }
    const int comparison = compareMagnitude(lhs.magnitude, rhs.magnitude);
    if (comparison == 0) return {};
    if (comparison > 0) {
        return {lhs.sign, subtractMagnitude(lhs.magnitude, rhs.magnitude)};
    }
    return {rhs.sign, subtractMagnitude(rhs.magnitude, lhs.magnitude)};
}

struct Binary64Dyadic {
    int sign = 0;
    std::uint64_t mantissa = 0;
    int exponent = 0;
};

[[nodiscard]] Binary64Dyadic decomposeBinary64(double value) noexcept {
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    const int sign = (bits >> 63U) == 0U ? 1 : -1;
    const std::uint64_t exponentBits = (bits >> 52U) & 0x7ffU;
    const std::uint64_t fraction = bits & ((std::uint64_t{1} << 52U) - 1U);
    if (exponentBits == 0U) {
        if (fraction == 0U) return {};
        return {sign, fraction, -1074};
    }
    return {sign, (std::uint64_t{1} << 52U) | fraction,
            static_cast<int>(exponentBits) - 1023 - 52};
}

[[nodiscard]] SignedBig<kCoordinateWords> exactCoordinate(
    const Binary64Dyadic& value, int commonExponent) noexcept {
    SignedBig<kCoordinateWords> result;
    if (value.sign == 0) return result;
    result.sign = value.sign;
    const unsigned shift = static_cast<unsigned>(value.exponent - commonExponent);
    const std::size_t word = shift / 64U;
    const unsigned bit = shift % 64U;
    result.magnitude[word] |= value.mantissa << bit;
    if (bit != 0U && word + 1U < result.magnitude.size()) {
        result.magnitude[word + 1U] |= value.mantissa >> (64U - bit);
    }
    return result;
}

[[nodiscard]] SignedBig<kProductWords> multiplyExact(
    const SignedBig<kCoordinateWords>& lhs,
    const SignedBig<kCoordinateWords>& rhs) noexcept {
    SignedBig<kProductWords> result;
    if (lhs.sign == 0 || rhs.sign == 0) return result;
    result.sign = lhs.sign * rhs.sign;
    for (std::size_t i = 0; i < lhs.magnitude.size(); ++i) {
        if (lhs.magnitude[i] == 0U) continue;
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < rhs.magnitude.size(); ++j) {
            const std::size_t out = i + j;
            const unsigned __int128 product =
                static_cast<unsigned __int128>(lhs.magnitude[i]) * rhs.magnitude[j] +
                result.magnitude[out] + carry;
            result.magnitude[out] = static_cast<std::uint64_t>(product);
            carry = static_cast<std::uint64_t>(product >> 64U);
        }
        if (i + rhs.magnitude.size() < result.magnitude.size()) {
            result.magnitude[i + rhs.magnitude.size()] = carry;
        }
    }
    return result;
}

[[nodiscard]] int exactOrientationSign(const Point2D& a, const Point2D& b,
                                       const Point2D& c) noexcept {
    const std::array<Binary64Dyadic, 6> values{
        decomposeBinary64(a.x), decomposeBinary64(a.y),
        decomposeBinary64(b.x), decomposeBinary64(b.y),
        decomposeBinary64(c.x), decomposeBinary64(c.y)};
    int commonExponent = 0;
    bool initialized = false;
    for (const auto& value : values) {
        if (value.sign == 0) continue;
        if (!initialized || value.exponent < commonExponent) {
            commonExponent = value.exponent;
            initialized = true;
        }
    }
    if (!initialized) return 0;

    const auto ax = exactCoordinate(values[0], commonExponent);
    const auto ay = exactCoordinate(values[1], commonExponent);
    const auto bx = exactCoordinate(values[2], commonExponent);
    const auto by = exactCoordinate(values[3], commonExponent);
    const auto cx = exactCoordinate(values[4], commonExponent);
    const auto cy = exactCoordinate(values[5], commonExponent);
    const auto bax = addSigned(bx, SignedBig<kCoordinateWords>{-ax.sign, ax.magnitude});
    const auto bay = addSigned(by, SignedBig<kCoordinateWords>{-ay.sign, ay.magnitude});
    const auto cax = addSigned(cx, SignedBig<kCoordinateWords>{-ax.sign, ax.magnitude});
    const auto cay = addSigned(cy, SignedBig<kCoordinateWords>{-ay.sign, ay.magnitude});
    const auto left = multiplyExact(bax, cay);
    const auto right = multiplyExact(bay, cax);
    const auto determinant = addSigned(
        left, SignedBig<kProductWords>{-right.sign, right.magnitude});
    return determinant.sign;
}

bool segmentsAreAdjacent(std::size_t i, std::size_t j, std::size_t edgeCount) noexcept {
    if (i == j) return true;
    if (i + 1 == j || j + 1 == i) return true;
    return (i == 0 && j + 1 == edgeCount) || (j == 0 && i + 1 == edgeCount);
}

} // namespace

double TolerancePolicy::scale(double magnitude) const noexcept {
    return absolute + relative * std::abs(magnitude);
}

bool TolerancePolicy::nearlyEqual(double a, double b, double magnitude) const noexcept {
    return std::abs(a - b) <= scale(std::max(std::abs(magnitude), std::abs(a - b)));
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

int orientationSign(const Point2D& a, const Point2D& b, const Point2D& c) noexcept {
    if (!std::isfinite(a.x) || !std::isfinite(a.y) ||
        !std::isfinite(b.x) || !std::isfinite(b.y) ||
        !std::isfinite(c.x) || !std::isfinite(c.y)) return 0;
    const double acx = a.x - c.x;
    const double bcx = b.x - c.x;
    const double acy = a.y - c.y;
    const double bcy = b.y - c.y;
    const double left = acx * bcy;
    const double right = acy * bcx;
    const double determinant = left - right;
    const double determinantSum = std::abs(left) + std::abs(right);
    constexpr double errorBound = 3.3306690738754716e-16;
    if (std::isfinite(determinant) &&
        std::abs(determinant) > errorBound * determinantSum) {
        return determinant > 0.0 ? 1 : -1;
    }
    return exactOrientationSign(a, b, c);
}

bool nearlyEqual(const Point2D& a, const Point2D& b, const TolerancePolicy& tol) noexcept {
    const double magnitude = std::max(std::abs(a.x - b.x), std::abs(a.y - b.y));
    return tol.nearlyEqual(a.x, b.x, magnitude) && tol.nearlyEqual(a.y, b.y, magnitude);
}

bool AABB2D::valid(const TolerancePolicy& tol) const noexcept {
    if (!std::isfinite(min.x) || !std::isfinite(min.y) ||
        !std::isfinite(max.x) || !std::isfinite(max.y)) return false;
    return min.x <= max.x + tol.scale(max.x - min.x) &&
           min.y <= max.y + tol.scale(max.y - min.y);
}

bool AABB2D::contains(const Point2D& p, const TolerancePolicy& tol) const noexcept {
    const double sx = tol.scale(std::max({std::abs(max.x - min.x),
                                          std::abs(p.x - min.x), std::abs(p.x - max.x)}));
    const double sy = tol.scale(std::max({std::abs(max.y - min.y),
                                          std::abs(p.y - min.y), std::abs(p.y - max.y)}));
    return p.x >= min.x - sx && p.x <= max.x + sx &&
           p.y >= min.y - sy && p.y <= max.y + sy;
}

double Polygon2D::signedArea() const noexcept {
    if (vertices.size() < 3) return 0.0;
    const Point2D origin = vertices.front();
    long double twiceArea = 0.0L;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % vertices.size()];
        const long double ax = static_cast<long double>(a.x) - origin.x;
        const long double ay = static_cast<long double>(a.y) - origin.y;
        const long double bx = static_cast<long double>(b.x) - origin.x;
        const long double by = static_cast<long double>(b.y) - origin.y;
        twiceArea += ax * by - bx * ay;
    }
    return static_cast<double>(0.5L * twiceArea);
}

double Polygon2D::area() const noexcept {
    return std::abs(signedArea());
}

std::optional<Point2D> Polygon2D::centroid(const TolerancePolicy& tol) const noexcept {
    if (vertices.size() < 3) return std::nullopt;

    const Point2D origin = vertices.front();
    long double twiceArea = 0.0L;
    long double cx = 0.0L;
    long double cy = 0.0L;
    double maxSpan = 0.0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % vertices.size()];
        const long double ax = static_cast<long double>(a.x) - origin.x;
        const long double ay = static_cast<long double>(a.y) - origin.y;
        const long double bx = static_cast<long double>(b.x) - origin.x;
        const long double by = static_cast<long double>(b.y) - origin.y;
        const long double term = ax * by - bx * ay;
        twiceArea += term;
        cx += (ax + bx) * term;
        cy += (ay + by) * term;
        maxSpan = std::max({maxSpan, std::abs(a.x - origin.x),
                            std::abs(a.y - origin.y)});
    }

    const long double roundoff =
        64.0L * std::numeric_limits<double>::epsilon() *
        static_cast<long double>(vertices.size()) * maxSpan * maxSpan;
    const long double areaEps = std::max(
        static_cast<long double>(tol.absolute) * tol.absolute, roundoff);
    if (std::abs(twiceArea) <= 2.0L * areaEps) return std::nullopt;
    return Point2D{origin.x + static_cast<double>(cx / (3.0L * twiceArea)),
                   origin.y + static_cast<double>(cy / (3.0L * twiceArea))};
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
    const double magnitude = std::max(std::sqrt(squaredNorm(ab)),
                                      std::sqrt(squaredNorm(ap)));
    const double eps = tol.scale(magnitude);
    const double length = std::sqrt(squaredNorm(ab));
    if (length <= eps) return nearlyEqual(p, segment.a, tol);

    if (std::abs(cross(ab, ap)) > eps * length) {
        return false;
    }
    const double projection = dot(ap, ab);
    const double projectionEps = eps * length;
    return projection >= -projectionEps &&
           projection <= squaredNorm(ab) + projectionEps;
}

SegmentIntersection intersectSegments(const Segment2D& lhs, const Segment2D& rhs,
                                      const TolerancePolicy& tol) noexcept {
    const Vector2D r = lhs.b - lhs.a;
    const Vector2D s = rhs.b - rhs.a;
    const double rr = squaredNorm(r);
    const double ss = squaredNorm(s);
    const double magnitude = std::max({std::sqrt(rr), std::sqrt(ss),
                                       std::sqrt(squaredNorm(rhs.a - lhs.a))});
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
    const double lenR = std::sqrt(rr);
    const double lenS = std::sqrt(ss);
    const double denominatorEps = eps * (lenR + lenS);

    if (std::abs(rxs) <= denominatorEps) {
        if (std::abs(qpxr) > eps * lenR) return {};

        const double t0 = dot(rhs.a - lhs.a, r) / rr;
        const double t1 = dot(rhs.b - lhs.a, r) / rr;
        const double lo = std::max(0.0, std::min(t0, t1));
        const double hi = std::min(1.0, std::max(t0, t1));
        const double tEps = eps / lenR;
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
    const double tEps = eps / lenR;
    const double uEps = eps / lenS;
    if (t < -tEps || t > 1.0 + tEps ||
        u < -uEps || u > 1.0 + uEps) {
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

        if (a.y <= p.y) {
            if (b.y > p.y && orientationSign(a, b, p) > 0) ++winding;
        } else if (b.y <= p.y && orientationSign(a, b, p) < 0) {
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

    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(vertices_[i].x) || !std::isfinite(vertices_[i].y)) {
            result.issues.push_back({BoundaryIssueCode::NonFiniteCoordinate, i, i,
                                     "boundary loop contains a non-finite coordinate"});
        }
    }
    if (!result.issues.empty()) return result;

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
    const double scale = std::max(box.max.x - box.min.x, box.max.y - box.min.y);
    const double areaEps = std::max(tol.absolute * tol.absolute,
                                    tol.relative * scale * scale);
    if (std::abs(signedArea) <= areaEps) {
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
