#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace cartmesh2d {
namespace {

[[nodiscard]] double backgroundArea(const AABB2D& box) noexcept {
    return (box.max.x - box.min.x) * (box.max.y - box.min.y);
}

[[nodiscard]] double areaTolerance(double referenceArea,
                                   const TolerancePolicy& tol) noexcept {
    return std::max(tol.absolute * tol.absolute,
                    tol.relative * std::abs(referenceArea));
}

[[nodiscard]] std::optional<Point2D> cutPolygonCentroid(
    const Polygon2D& polygon, double areaEps) noexcept {
    if (polygon.vertices.size() < 3) return std::nullopt;

    double twiceArea = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
        const auto& a = polygon.vertices[i];
        const auto& b = polygon.vertices[(i + 1) % polygon.vertices.size()];
        const double term = a.x * b.y - b.x * a.y;
        twiceArea += term;
        cx += (a.x + b.x) * term;
        cy += (a.y + b.y) * term;
    }
    if (std::abs(twiceArea) <= 2.0 * areaEps) return std::nullopt;
    return Point2D{cx / (3.0 * twiceArea), cy / (3.0 * twiceArea)};
}

[[nodiscard]] Polygon2D rectanglePolygon(const AABB2D& box) {
    return {{{box.min.x, box.min.y},
             {box.max.x, box.min.y},
             {box.max.x, box.max.y},
             {box.min.x, box.max.y}}};
}

[[nodiscard]] bool scalarNear(double a, double b, const TolerancePolicy& tol) noexcept {
    return std::abs(a - b) <= tol.scale(std::abs(a - b));
}

[[nodiscard]] bool pointNear(const Point2D& a, const Point2D& b,
                             const TolerancePolicy& tol) noexcept {
    return scalarNear(a.x, b.x, tol) && scalarNear(a.y, b.y, tol);
}

[[nodiscard]] Point2D snapPointToBox(Point2D point, const AABB2D& box,
                                     const TolerancePolicy& tol) noexcept {
    if (scalarNear(point.x, box.min.x, tol)) point.x = box.min.x;
    if (scalarNear(point.x, box.max.x, tol)) point.x = box.max.x;
    if (scalarNear(point.y, box.min.y, tol)) point.y = box.min.y;
    if (scalarNear(point.y, box.max.y, tol)) point.y = box.max.y;
    return point;
}

void removeConsecutiveDuplicates(std::vector<Point2D>& vertices,
                                 const TolerancePolicy& tol) {
    std::vector<Point2D> result;
    result.reserve(vertices.size());
    for (const auto& point : vertices) {
        if (result.empty() || !pointNear(result.back(), point, tol)) result.push_back(point);
    }
    if (result.size() > 1 && pointNear(result.front(), result.back(), tol)) result.pop_back();
    vertices.swap(result);
}

void removeCollinearVertices(std::vector<Point2D>& vertices,
                             const TolerancePolicy& tol) {
    removeConsecutiveDuplicates(vertices, tol);
    if (vertices.size() < 3) return;
    bool changed = true;
    while (changed && vertices.size() >= 3) {
        changed = false;
        std::vector<Point2D> out;
        out.reserve(vertices.size());
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const Point2D& prev = vertices[(i + vertices.size() - 1) % vertices.size()];
            const Point2D& cur = vertices[i];
            const Point2D& next = vertices[(i + 1) % vertices.size()];
            const Vector2D a = cur - prev;
            const Vector2D b = next - cur;
            const double lenA = std::sqrt(squaredNorm(a));
            const double lenB = std::sqrt(squaredNorm(b));
            const double lengthEps = tol.scale(std::max(lenA, lenB));
            if (lenA <= lengthEps || lenB <= lengthEps) {
                out.push_back(cur);
                continue;
            }
            const double sinAngle = std::abs(cross(a, b)) / (lenA * lenB);
            if (sinAngle <= tol.scale(1.0) && dot(a, b) >= 0.0) {
                changed = true;
                continue;
            }
            out.push_back(cur);
        }
        vertices.swap(out);
    }
}

[[nodiscard]] bool edgesAreAdjacent(std::size_t i, std::size_t j,
                                    std::size_t edgeCount) noexcept {
    if (i == j) return true;
    if (i + 1 == j || j + 1 == i) return true;
    return (i == 0 && j + 1 == edgeCount) ||
           (j == 0 && i + 1 == edgeCount);
}

[[nodiscard]] bool localSegmentsIntersect(const Segment2D& lhs,
                                          const Segment2D& rhs,
                                          const TolerancePolicy& tol) noexcept {
    const Vector2D r = lhs.b - lhs.a;
    const Vector2D s = rhs.b - rhs.a;
    const double lenR = std::sqrt(squaredNorm(r));
    const double lenS = std::sqrt(squaredNorm(s));
    const double lengthEps = tol.scale(std::max({lenR, lenS,
        std::sqrt(squaredNorm(rhs.a - lhs.a))}));
    if (lenR <= lengthEps || lenS <= lengthEps) return true;

    const Vector2D q = rhs.a - lhs.a;
    const double denominator = cross(r, s);
    const double angularEps = tol.scale(1.0);
    if (std::abs(denominator) <= angularEps * lenR * lenS) {
        const double distance = std::abs(cross(q, r)) / lenR;
        if (distance > lengthEps) return false;
        const double rr = squaredNorm(r);
        const double t0 = dot(rhs.a - lhs.a, r) / rr;
        const double t1 = dot(rhs.b - lhs.a, r) / rr;
        const double lo = std::max(0.0, std::min(t0, t1));
        const double hi = std::min(1.0, std::max(t0, t1));
        const double paramEps = lengthEps / lenR;
        return hi >= lo - paramEps;
    }

    const double t = cross(q, s) / denominator;
    const double u = cross(q, r) / denominator;
    const double paramEps = std::max(lengthEps / lenR, lengthEps / lenS);
    return t >= -paramEps && t <= 1.0 + paramEps &&
           u >= -paramEps && u <= 1.0 + paramEps;
}

[[nodiscard]] bool simplePolygonLoop(const Polygon2D& polygon,
                                     const TolerancePolicy& tol) noexcept {
    const std::size_t n = polygon.vertices.size();
    if (n < 3) return false;
    for (std::size_t i = 0; i < n; ++i) {
        const Segment2D a{polygon.vertices[i], polygon.vertices[(i + 1) % n]};
        if (pointNear(a.a, a.b, tol)) return false;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (edgesAreAdjacent(i, j, n)) continue;
            const Segment2D b{polygon.vertices[j], polygon.vertices[(j + 1) % n]};
            if (localSegmentsIntersect(a, b, tol)) return false;
        }
    }
    return true;
}

[[nodiscard]] std::optional<Segment2D> clipSegmentToAABB(
    const Segment2D& segment, const AABB2D& box,
    const TolerancePolicy& tol) noexcept {
    double t0 = 0.0;
    double t1 = 1.0;
    const double dx = segment.b.x - segment.a.x;
    const double dy = segment.b.y - segment.a.y;

    auto update = [&](double p, double q) {
        const double eps = tol.scale(std::max(std::abs(p), std::abs(q)));
        if (std::abs(p) <= eps) return q >= -eps;
        const double ratio = q / p;
        if (p < 0.0) t0 = std::max(t0, ratio);
        else t1 = std::min(t1, ratio);
        return t0 <= t1 + eps;
    };

    if (!update(-dx, segment.a.x - box.min.x) ||
        !update(dx, box.max.x - segment.a.x) ||
        !update(-dy, segment.a.y - box.min.y) ||
        !update(dy, box.max.y - segment.a.y)) return std::nullopt;

    t0 = std::clamp(t0, 0.0, 1.0);
    t1 = std::clamp(t1, 0.0, 1.0);
    Segment2D clipped{{segment.a.x + t0 * dx, segment.a.y + t0 * dy},
                      {segment.a.x + t1 * dx, segment.a.y + t1 * dy}};
    clipped.a = snapPointToBox(clipped.a, box, tol);
    clipped.b = snapPointToBox(clipped.b, box, tol);
    if (pointNear(clipped.a, clipped.b, tol)) return std::nullopt;
    return clipped;
}

[[nodiscard]] bool liesOnBoxPerimeter(const Segment2D& segment,
                                      const AABB2D& box,
                                      const TolerancePolicy& tol) noexcept {
    return (scalarNear(segment.a.x, box.min.x, tol) && scalarNear(segment.b.x, box.min.x, tol)) ||
           (scalarNear(segment.a.x, box.max.x, tol) && scalarNear(segment.b.x, box.max.x, tol)) ||
           (scalarNear(segment.a.y, box.min.y, tol) && scalarNear(segment.b.y, box.min.y, tol)) ||
           (scalarNear(segment.a.y, box.max.y, tol) && scalarNear(segment.b.y, box.max.y, tol));
}

[[nodiscard]] bool pointStateIsFluid(PointInPolygon state,
                                     FluidRegion2D fluidRegion) noexcept {
    if (state == PointInPolygon::Boundary) return true;
    if (fluidRegion == FluidRegion2D::Exterior) return state == PointInPolygon::Outside;
    return state == PointInPolygon::Inside;
}

[[nodiscard]] std::vector<Segment2D> collectEmbeddedBoundary(
    const BoundaryRegion2D& boundary, const AABB2D& box,
    FluidRegion2D fluidRegion, const TolerancePolicy& tol) {
    std::vector<Segment2D> fragments;
    for (const auto& loop : boundary.loops()) {
        const auto& vertices = loop.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const Segment2D edge{vertices[i], vertices[(i + 1) % vertices.size()]};
            auto clipped = clipSegmentToAABB(edge, box, tol);
            if (!clipped || liesOnBoxPerimeter(*clipped, box, tol)) continue;
            if (fluidRegion == FluidRegion2D::Exterior) std::swap(clipped->a, clipped->b);
            fragments.push_back(*clipped);
        }
    }
    return fragments;
}

[[nodiscard]] std::size_t countFragmentComponents(
    const std::vector<Segment2D>& fragments,
    const TolerancePolicy& tol) {
    if (fragments.empty()) return 0;
    std::vector<bool> visited(fragments.size(), false);
    std::size_t components = 0;
    for (std::size_t start = 0; start < fragments.size(); ++start) {
        if (visited[start]) continue;
        ++components;
        std::queue<std::size_t> queue;
        queue.push(start);
        visited[start] = true;
        while (!queue.empty()) {
            const std::size_t current = queue.front();
            queue.pop();
            for (std::size_t candidate = 0; candidate < fragments.size(); ++candidate) {
                if (visited[candidate]) continue;
                const auto& a = fragments[current];
                const auto& b = fragments[candidate];
                if (pointNear(a.a, b.a, tol) || pointNear(a.a, b.b, tol) ||
                    pointNear(a.b, b.a, tol) || pointNear(a.b, b.b, tol)) {
                    visited[candidate] = true;
                    queue.push(candidate);
                }
            }
        }
    }
    return components;
}

[[nodiscard]] bool pointOnBoxPerimeter(const Point2D& point, const AABB2D& box,
                                       const TolerancePolicy& tol) noexcept {
    return box.contains(point,tol) &&
           (scalarNear(point.x,box.min.x,tol) || scalarNear(point.x,box.max.x,tol) ||
            scalarNear(point.y,box.min.y,tol) || scalarNear(point.y,box.max.y,tol));
}

[[nodiscard]] bool hasClosedEmbeddedComponent(
    const std::vector<Segment2D>& fragments, const AABB2D& box,
    const TolerancePolicy& tol) {
    std::vector<bool> visited(fragments.size(),false);
    for (std::size_t start=0;start<fragments.size();++start) {
        if (visited[start]) continue;
        bool touchesPerimeter=false;
        std::queue<std::size_t> pending;
        pending.push(start);
        visited[start]=true;
        while (!pending.empty()) {
            const std::size_t current=pending.front();
            pending.pop();
            const auto& fragment=fragments[current];
            touchesPerimeter = touchesPerimeter ||
                pointOnBoxPerimeter(fragment.a,box,tol) ||
                pointOnBoxPerimeter(fragment.b,box,tol);
            for (std::size_t candidate=0;candidate<fragments.size();++candidate) {
                if (visited[candidate]) continue;
                const auto& other=fragments[candidate];
                if (pointNear(fragment.a,other.a,tol) || pointNear(fragment.a,other.b,tol) ||
                    pointNear(fragment.b,other.a,tol) || pointNear(fragment.b,other.b,tol)) {
                    visited[candidate]=true;
                    pending.push(candidate);
                }
            }
        }
        if (!touchesPerimeter) return true;
    }
    return false;
}

struct HalfEdge {
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t twin = 0;
    std::size_t next = 0;
};

[[nodiscard]] std::size_t findOrAddPoint(std::vector<Point2D>& points,
                                         const Point2D& point,
                                         const TolerancePolicy& tol) {
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (pointNear(points[i], point, tol)) return i;
    }
    points.push_back(point);
    return points.size() - 1;
}

void addUndirectedEdge(std::vector<HalfEdge>& edges, std::size_t from, std::size_t to) {
    if (from == to) return;
    for (std::size_t i = 0; i + 1 < edges.size(); i += 2) {
        const auto& edge = edges[i];
        if ((edge.from == from && edge.to == to) ||
            (edge.from == to && edge.to == from)) return;
    }
    const std::size_t first = edges.size();
    edges.push_back({from, to, first + 1, first});
    edges.push_back({to, from, first, first + 1});
}

[[nodiscard]] double sideParameter(const Point2D& p, const Segment2D& side) noexcept {
    const Vector2D d = side.b - side.a;
    const double denom = squaredNorm(d);
    if (denom <= 0.0) return 0.0;
    return dot(p - side.a, d) / denom;
}

struct LocalFluidComponent {
    Polygon2D polygon;
    std::vector<Segment2D> embedded;
};

struct LocalIntersectionResult {
    std::vector<LocalFluidComponent> components;
    std::vector<Segment2D> embedded;
    bool graphInvalid = false;
    bool hasHole = false;
    std::string graphFailure;
};

[[nodiscard]] bool pointOnPolygonBoundary(const Point2D& point,
                                          const Polygon2D& polygon,
                                          const TolerancePolicy& tol) noexcept {
    for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
        const Segment2D edge{polygon.vertices[i],
                             polygon.vertices[(i + 1) % polygon.vertices.size()]};
        if (pointOnSegment(point, edge, tol)) return true;
    }
    return false;
}

[[nodiscard]] std::optional<Point2D> interiorProbe(const Polygon2D& polygon,
                                                   const AABB2D& box,
                                                   const TolerancePolicy& tol) noexcept {
    if (const auto centroid = polygon.centroid(tol)) {
        if (classifyPointInPolygon(*centroid, polygon, tol) == PointInPolygon::Inside) {
            return centroid;
        }
    }

    const double cellScale = std::max(box.max.x - box.min.x, box.max.y - box.min.y);
    const double offset = std::max(64.0 * tol.scale(cellScale),
                                   1.0e-8 * cellScale);
    for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
        const Point2D& a = polygon.vertices[i];
        const Point2D& b = polygon.vertices[(i + 1) % polygon.vertices.size()];
        const Vector2D edge = b - a;
        const double length = std::sqrt(squaredNorm(edge));
        if (length <= offset) continue;
        const Point2D probe{0.5 * (a.x + b.x) - offset * edge.y / length,
                            0.5 * (a.y + b.y) + offset * edge.x / length};
        if (classifyPointInPolygon(probe, polygon, tol) == PointInPolygon::Inside) {
            return probe;
        }
    }
    return std::nullopt;
}

[[nodiscard]] LocalIntersectionResult buildLocalIntersection(
    const BoundaryRegion2D& boundary, const AABB2D& box,
    FluidRegion2D fluidRegion, const TolerancePolicy& tol) {
    LocalIntersectionResult result;
    result.embedded = collectEmbeddedBoundary(boundary, box, fluidRegion, tol);
    if (result.embedded.empty()) return result;
    if (hasClosedEmbeddedComponent(result.embedded,box,tol)) {
        result.hasHole=true;
        return result;
    }

    std::vector<Point2D> points;
    std::vector<HalfEdge> edges;
    for (const auto& fragment : result.embedded) {
        const std::size_t a = findOrAddPoint(points, fragment.a, tol);
        const std::size_t b = findOrAddPoint(points, fragment.b, tol);
        addUndirectedEdge(edges, a, b);
    }

    const std::vector<Segment2D> sides{
        {{box.min.x, box.min.y}, {box.max.x, box.min.y}},
        {{box.max.x, box.min.y}, {box.max.x, box.max.y}},
        {{box.max.x, box.max.y}, {box.min.x, box.max.y}},
        {{box.min.x, box.max.y}, {box.min.x, box.min.y}}
    };
    for (const auto& side : sides) {
        std::vector<Point2D> sidePoints{side.a, side.b};
        for (const auto& fragment : result.embedded) {
            if (pointOnSegment(fragment.a, side, tol)) sidePoints.push_back(fragment.a);
            if (pointOnSegment(fragment.b, side, tol)) sidePoints.push_back(fragment.b);
        }
        std::sort(sidePoints.begin(), sidePoints.end(), [&](const Point2D& lhs, const Point2D& rhs) {
            return sideParameter(lhs, side) < sideParameter(rhs, side);
        });
        std::vector<Point2D> unique;
        for (const auto& p : sidePoints) {
            if (unique.empty() || !pointNear(unique.back(), p, tol)) unique.push_back(p);
        }
        for (std::size_t i = 0; i + 1 < unique.size(); ++i) {
            if (pointNear(unique[i], unique[i + 1], tol)) continue;
            const std::size_t a = findOrAddPoint(points, unique[i], tol);
            const std::size_t b = findOrAddPoint(points, unique[i + 1], tol);
            addUndirectedEdge(edges, a, b);
        }
    }

    if (edges.empty()) return result;

    // Build a deterministic planar half-edge rotation system.  For a directed
    // edge u->v, the predecessor of v->u in CCW order is the next edge around
    // the face retained on the left.  Unlike the previous one-outgoing-edge
    // assumption, this remains well-defined at aligned corners and tangencies.
    std::vector<std::vector<std::size_t>> outgoing(points.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        outgoing[edges[i].from].push_back(i);
    }
    for (auto& fan : outgoing) {
        std::sort(fan.begin(), fan.end(), [&](std::size_t lhs, std::size_t rhs) {
            const Point2D& origin = points[edges[lhs].from];
            const Vector2D a = points[edges[lhs].to] - points[edges[lhs].from];
            const Vector2D b = points[edges[rhs].to] - points[edges[rhs].from];
            const bool upperA = a.y > 0.0 || (a.y == 0.0 && a.x >= 0.0);
            const bool upperB = b.y > 0.0 || (b.y == 0.0 && b.x >= 0.0);
            if (upperA != upperB) return upperA;
            const int turn = orientationSign(origin, points[edges[lhs].to],
                                              points[edges[rhs].to]);
            if (turn != 0) return turn > 0;
            const double lengthA = squaredNorm(a);
            const double lengthB = squaredNorm(b);
            if (lengthA != lengthB) return lengthA < lengthB;
            return lhs < rhs;
        });
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& fan = outgoing[edges[i].to];
        const auto twinIt = std::find(fan.begin(), fan.end(), edges[i].twin);
        if (twinIt == fan.end() || fan.empty()) {
            result.graphInvalid = true;
            result.graphFailure = "half-edge twin is missing from the destination rotation fan";
            return result;
        }
        const std::size_t position = static_cast<std::size_t>(std::distance(fan.begin(), twinIt));
        edges[i].next = fan[(position + fan.size() - 1) % fan.size()];
    }

    std::vector<bool> used(edges.size(), false);
    const double areaEps = areaTolerance(backgroundArea(box), tol);
    for (std::size_t startEdge = 0; startEdge < edges.size(); ++startEdge) {
        if (used[startEdge]) continue;
        std::vector<Point2D> loop;
        std::size_t edgeId = startEdge;
        std::size_t guard = 0;
        while (guard++ <= edges.size()) {
            if (used[edgeId]) {
                if (edgeId == startEdge) break;
                result.graphInvalid = true;
                result.graphFailure = "face traversal reached a half-edge already owned by another face";
                return result;
            }
            used[edgeId] = true;
            const auto edge = edges[edgeId];
            loop.push_back(points[edge.from]);
            edgeId = edge.next;
            if (edgeId == startEdge) break;
        }
        if (guard > edges.size() + 1) {
            result.graphInvalid = true;
            result.graphFailure = "half-edge face traversal did not close";
            return result;
        }
        removeCollinearVertices(loop, tol);
        if (loop.size() < 3) continue;
        Polygon2D polygon{loop};
        if (polygon.area() <= areaEps) continue;
        if (!simplePolygonLoop(polygon, tol)) {
            result.graphInvalid = true;
            result.graphFailure = "half-edge traversal produced a non-simple polygon";
            return result;
        }

        if (polygon.signedArea() < 0.0) continue; // unbounded/right-hand face

        const auto probe = interiorProbe(polygon, box, tol);
        if (!probe) {
            result.graphInvalid = true;
            result.graphFailure = "half-edge face has no certified interior probe (vertices=" +
                                  std::to_string(polygon.vertices.size()) +
                                  ", area=" + std::to_string(polygon.area()) + ")";
            return result;
        }
        const auto state = boundary.classifyPoint(*probe, tol);
        if (!pointStateIsFluid(state, fluidRegion)) continue;

        LocalFluidComponent component;
        component.polygon = std::move(polygon);
        for (const auto& fragment : result.embedded) {
            const Point2D midpoint{0.5 * (fragment.a.x + fragment.b.x),
                                   0.5 * (fragment.a.y + fragment.b.y)};
            if (pointOnPolygonBoundary(midpoint, component.polygon, tol)) {
                component.embedded.push_back(fragment);
            }
        }
        // A true partial Cut-cell component must touch the embedded boundary.
        // This also rejects the outer loop of a disconnected polygon-with-hole
        // arrangement instead of silently filling the hole.
        if (component.embedded.empty()) {
            result.hasHole = true;
            continue;
        }
        result.components.push_back(std::move(component));
    }
    return result;
}

void setFullFluidCell(CutCell2D& result, const AABB2D& box,
                      double fullArea, double areaEps) {
    result.kind = CutCellKind::Full;
    result.fluidPolygon = rectanglePolygon(box);
    result.area = fullArea;
    result.areaFraction = 1.0;
    result.centroid = cutPolygonCentroid(result.fluidPolygon, areaEps);
    result.embeddedBoundary.clear();
}

void setEmptyFluidCell(CutCell2D& result) {
    result.kind = CutCellKind::Empty;
    result.fluidPolygon.vertices.clear();
    result.area = 0.0;
    result.areaFraction = 0.0;
    result.centroid.reset();
    result.embeddedBoundary.clear();
}

[[nodiscard]] CutCell2D unsupportedCell(const AABB2D& box,
                                        CutCellIssueCode code,
                                        std::string message) {
    CutCell2D result;
    result.backgroundBounds = box;
    result.kind = CutCellKind::Unsupported;
    result.issues.push_back({code, std::move(message)});
    return result;
}

[[nodiscard]] std::vector<CutCell2D> buildCutCellsImpl(
    const AABB2D& box, CellClass classification,
    const BoundaryRegion2D& inputBoundary, FluidRegion2D fluidRegion,
    const TolerancePolicy& tol) {
    const double fullArea = backgroundArea(box);
    if (!(fullArea > areaTolerance(fullArea, tol))) {
        return {unsupportedCell(box, CutCellIssueCode::InvalidBackgroundCell,
                                "background cell has non-positive area")};
    }

    const auto inputDiagnostics = inputBoundary.diagnose(tol);
    if (!inputDiagnostics.valid()) {
        return {unsupportedCell(box, CutCellIssueCode::InvalidBoundary,
                                "boundary region failed geometry diagnostics")};
    }

    BoundaryRegion2D boundary = inputBoundary;
    if (!boundary.normalizeAlternating(tol)) {
        return {unsupportedCell(box, CutCellIssueCode::InvalidBoundary,
                                "boundary region could not be normalized by nesting depth")};
    }

    const double areaEps = areaTolerance(fullArea, tol);
    CutCell2D simple;
    simple.backgroundBounds = box;

    if (classification == CellClass::Outside) {
        if (fluidRegion == FluidRegion2D::Exterior) {
            setFullFluidCell(simple, box, fullArea, areaEps);
        } else {
            setEmptyFluidCell(simple);
        }
        return {std::move(simple)};
    }
    if (classification == CellClass::Inside) {
        if (fluidRegion == FluidRegion2D::Exterior) {
            setEmptyFluidCell(simple);
        } else {
            setFullFluidCell(simple, box, fullArea, areaEps);
        }
        return {std::move(simple)};
    }

    const auto local = buildLocalIntersection(boundary, box, fluidRegion, tol);
    if (local.graphInvalid) {
        return {unsupportedCell(box, CutCellIssueCode::DegeneratePolygon,
                                "local Cut-cell boundary graph is invalid: " + local.graphFailure)};
    }
    if (local.hasHole) {
        return {unsupportedCell(box, CutCellIssueCode::MultipleEmbeddedComponents,
                                "local fluid region contains a hole; refine the leaf or add polygon-with-holes support")};
    }

    if (local.components.empty()) {
        const Point2D center{0.5 * (box.min.x + box.max.x),
                             0.5 * (box.min.y + box.max.y)};
        const auto state = boundary.classifyPoint(center, tol);
        if (pointStateIsFluid(state, fluidRegion)) {
            setFullFluidCell(simple, box, fullArea, areaEps);
        } else {
            setEmptyFluidCell(simple);
        }
        return {std::move(simple)};
    }

    const std::size_t fragmentComponentCount = countFragmentComponents(local.embedded, tol);
    if (fragmentComponentCount == 0) {
        return {unsupportedCell(box, CutCellIssueCode::MissingEmbeddedBoundary,
                                "partial fluid polygon has no embedded-boundary fragment")};
    }

    std::vector<CutCell2D> result;
    result.reserve(local.components.size());
    double totalArea = 0.0;
    for (const auto& localComponent : local.components) {
        const auto& polygon = localComponent.polygon;
        CutCell2D component;
        component.backgroundBounds = box;
        component.fluidPolygon = polygon;
        component.embeddedBoundary = localComponent.embedded;
        component.area = component.fluidPolygon.area();
        if (component.area <= areaEps) continue;
        if (component.area > fullArea + areaEps) {
            return {unsupportedCell(box, CutCellIssueCode::AreaOutOfRange,
                                    "local Cut-cell fluid component exceeds background-cell area")};
        }
        component.areaFraction = std::clamp(component.area / fullArea, 0.0, 1.0);
        component.centroid = cutPolygonCentroid(component.fluidPolygon, areaEps);
        if (!component.centroid) {
            return {unsupportedCell(box, CutCellIssueCode::DegeneratePolygon,
                                    "local Cut-cell polygon has no valid centroid")};
        }
        component.kind = CutCellKind::Cut;
        totalArea += component.area;
        result.push_back(std::move(component));
    }

    if (result.empty()) {
        setEmptyFluidCell(simple);
        return {std::move(simple)};
    }
    if (totalArea > fullArea + areaEps) {
        return {unsupportedCell(box, CutCellIssueCode::AreaOutOfRange,
                                "sum of local fluid components exceeds background-cell area")};
    }
    if (result.size() == 1 && fullArea - result.front().area <= areaEps) {
        setFullFluidCell(result.front(), box, fullArea, areaEps);
    }
    return result;
}

[[nodiscard]] CutCell2D buildCutCellImpl(
    const AABB2D& box, CellClass classification,
    const BoundaryRegion2D& inputBoundary, FluidRegion2D fluidRegion,
    const TolerancePolicy& tol) {
    auto components = buildCutCellsImpl(box, classification, inputBoundary, fluidRegion, tol);
    if (components.size() == 1) return std::move(components.front());

    CutCell2D result = unsupportedCell(
        box, CutCellIssueCode::MultipleEmbeddedComponents,
        "multiple disconnected fluid components occur in one leaf; use buildCutCells() for solver topology");
    return result;
}

} // namespace

CutCell2D buildCutCell(const AABB2D& box, CellClass classification,
                       const BoundaryLoop& boundary,
                       const TolerancePolicy& tol) {
    return buildCutCellImpl(box, classification, BoundaryRegion2D(boundary),
                            FluidRegion2D::Exterior, tol);
}

CutCell2D buildCutCell(const AABB2D& box, CellClass classification,
                       const BoundaryLoop& boundary, FluidRegion2D fluidRegion,
                       const TolerancePolicy& tol) {
    return buildCutCellImpl(box, classification, BoundaryRegion2D(boundary), fluidRegion, tol);
}

CutCell2D buildCutCell(const AABB2D& box, CellClass classification,
                       const BoundaryRegion2D& boundary, FluidRegion2D fluidRegion,
                       const TolerancePolicy& tol) {
    return buildCutCellImpl(box, classification, boundary, fluidRegion, tol);
}

CutCell2D buildCutCell(const QuadtreeLeaf2D& leaf,
                       const BoundaryLoop& boundary,
                       const TolerancePolicy& tol) {
    CutCell2D result = buildCutCellImpl(leaf.bounds, leaf.classification, BoundaryRegion2D(boundary),
                                        FluidRegion2D::Exterior, tol);
    result.sourceId = leaf.id;
    result.sourceKey = leaf.key;
    return result;
}

CutCell2D buildCutCell(const QuadtreeLeaf2D& leaf,
                       const BoundaryLoop& boundary, FluidRegion2D fluidRegion,
                       const TolerancePolicy& tol) {
    CutCell2D result = buildCutCellImpl(leaf.bounds, leaf.classification, BoundaryRegion2D(boundary),
                                        fluidRegion, tol);
    result.sourceId = leaf.id;
    result.sourceKey = leaf.key;
    return result;
}

CutCell2D buildCutCell(const QuadtreeLeaf2D& leaf,
                       const BoundaryRegion2D& boundary, FluidRegion2D fluidRegion,
                       const TolerancePolicy& tol) {
    CutCell2D result = buildCutCellImpl(leaf.bounds, leaf.classification, boundary,
                                        fluidRegion, tol);
    result.sourceId = leaf.id;
    result.sourceKey = leaf.key;
    return result;
}

std::vector<CutCell2D> buildCutCells(const AABB2D& box, CellClass classification,
                                     const BoundaryLoop& boundary,
                                     const TolerancePolicy& tol) {
    return buildCutCellsImpl(box, classification, BoundaryRegion2D(boundary), FluidRegion2D::Exterior, tol);
}

std::vector<CutCell2D> buildCutCells(const AABB2D& box, CellClass classification,
                                     const BoundaryLoop& boundary, FluidRegion2D fluidRegion,
                                     const TolerancePolicy& tol) {
    return buildCutCellsImpl(box, classification, BoundaryRegion2D(boundary), fluidRegion, tol);
}

std::vector<CutCell2D> buildCutCells(const AABB2D& box, CellClass classification,
                                     const BoundaryRegion2D& boundary,
                                     FluidRegion2D fluidRegion,
                                     const TolerancePolicy& tol) {
    return buildCutCellsImpl(box, classification, boundary, fluidRegion, tol);
}

std::vector<CutCell2D> buildCutCells(const QuadtreeLeaf2D& leaf,
                                     const BoundaryLoop& boundary,
                                     const TolerancePolicy& tol) {
    auto result = buildCutCellsImpl(leaf.bounds, leaf.classification, BoundaryRegion2D(boundary),
                                    FluidRegion2D::Exterior, tol);
    for (auto& component : result) {
        component.sourceId = leaf.id;
        component.sourceKey = leaf.key;
    }
    return result;
}

std::vector<CutCell2D> buildCutCells(const QuadtreeLeaf2D& leaf,
                                     const BoundaryLoop& boundary, FluidRegion2D fluidRegion,
                                     const TolerancePolicy& tol) {
    auto result = buildCutCellsImpl(leaf.bounds, leaf.classification, BoundaryRegion2D(boundary),
                                    fluidRegion, tol);
    for (auto& component : result) {
        component.sourceId = leaf.id;
        component.sourceKey = leaf.key;
    }
    return result;
}

std::vector<CutCell2D> buildCutCells(const QuadtreeLeaf2D& leaf,
                                     const BoundaryRegion2D& boundary,
                                     FluidRegion2D fluidRegion,
                                     const TolerancePolicy& tol) {
    auto result = buildCutCellsImpl(leaf.bounds, leaf.classification, boundary,
                                    fluidRegion, tol);
    for (auto& component : result) {
        component.sourceId = leaf.id;
        component.sourceKey = leaf.key;
    }
    return result;
}

} // namespace cartmesh2d
