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

[[nodiscard]] Polygon2D rectanglePolygon(const AABB2D& box) {
    return {{{box.min.x, box.min.y},
             {box.max.x, box.min.y},
             {box.max.x, box.max.y},
             {box.min.x, box.max.y}}};
}

[[nodiscard]] bool scalarNear(double a, double b, const TolerancePolicy& tol) noexcept {
    return std::abs(a - b) <= tol.scale(std::max({1.0, std::abs(a), std::abs(b)}));
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
            const double magnitude = std::sqrt(squaredNorm(a) * squaredNorm(b));
            const double eps = tol.scale(std::max(1.0, magnitude));
            if (std::abs(cross(a, b)) <= eps && dot(a, b) >= -eps) {
                changed = true;
                continue;
            }
            out.push_back(cur);
        }
        vertices.swap(out);
    }
}

[[nodiscard]] std::optional<Segment2D> clipSegmentToAABB(
    const Segment2D& segment, const AABB2D& box,
    const TolerancePolicy& tol) noexcept {
    double t0 = 0.0;
    double t1 = 1.0;
    const double dx = segment.b.x - segment.a.x;
    const double dy = segment.b.y - segment.a.y;

    auto update = [&](double p, double q) {
        const double eps = tol.scale(std::max({1.0, std::abs(p), std::abs(q)}));
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

[[nodiscard]] std::vector<Segment2D> collectEmbeddedBoundary(
    const BoundaryLoop& boundary, const AABB2D& box,
    const TolerancePolicy& tol) {
    std::vector<Segment2D> fragments;
    const auto& vertices = boundary.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Segment2D edge{vertices[i], vertices[(i + 1) % vertices.size()]};
        const auto clipped = clipSegmentToAABB(edge, box, tol);
        if (clipped && !liesOnBoxPerimeter(*clipped, box, tol)) fragments.push_back(*clipped);
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

struct DirectedEdge {
    std::size_t from = 0;
    std::size_t to = 0;
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

void addDirectedEdge(std::vector<DirectedEdge>& edges, std::size_t from, std::size_t to) {
    if (from == to) return;
    for (const auto& edge : edges) {
        if (edge.from == from && edge.to == to) return;
    }
    edges.push_back({from, to});
}

[[nodiscard]] double sideParameter(const Point2D& p, const Segment2D& side) noexcept {
    const Vector2D d = side.b - side.a;
    const double denom = squaredNorm(d);
    if (denom <= 0.0) return 0.0;
    return dot(p - side.a, d) / denom;
}

struct LocalIntersectionResult {
    std::vector<Polygon2D> components;
    std::vector<Segment2D> embedded;
    bool graphInvalid = false;
};

[[nodiscard]] LocalIntersectionResult buildLocalIntersection(
    const BoundaryLoop& boundary, const AABB2D& box,
    const TolerancePolicy& tol) {
    LocalIntersectionResult result;
    result.embedded = collectEmbeddedBoundary(boundary, box, tol);
    if (result.embedded.empty()) return result;

    std::vector<Point2D> points;
    std::vector<DirectedEdge> edges;
    for (const auto& fragment : result.embedded) {
        const std::size_t a = findOrAddPoint(points, fragment.a, tol);
        const std::size_t b = findOrAddPoint(points, fragment.b, tol);
        addDirectedEdge(edges, a, b);
    }

    const std::vector<Segment2D> sides{
        {{box.min.x, box.min.y}, {box.max.x, box.min.y}},
        {{box.max.x, box.min.y}, {box.max.x, box.max.y}},
        {{box.max.x, box.max.y}, {box.min.x, box.max.y}},
        {{box.min.x, box.max.y}, {box.min.x, box.min.y}}
    };
    const Polygon2D boundaryPolygon = boundary.polygon();

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
            const Point2D mid{0.5 * (unique[i].x + unique[i + 1].x),
                              0.5 * (unique[i].y + unique[i + 1].y)};
            const auto state = classifyPointInPolygon(mid, boundaryPolygon, tol);
            if (state == PointInPolygon::Inside || state == PointInPolygon::Boundary) {
                const std::size_t a = findOrAddPoint(points, unique[i], tol);
                const std::size_t b = findOrAddPoint(points, unique[i + 1], tol);
                addDirectedEdge(edges, a, b);
            }
        }
    }

    if (edges.empty()) return result;

    std::vector<std::vector<std::size_t>> outgoing(points.size());
    std::vector<std::size_t> indegree(points.size(), 0);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        outgoing[edges[i].from].push_back(i);
        ++indegree[edges[i].to];
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (outgoing[i].size() != indegree[i]) {
            result.graphInvalid = true;
            return result;
        }
        if (!outgoing[i].empty() && outgoing[i].size() != 1) {
            result.graphInvalid = true;
            return result;
        }
    }

    std::vector<bool> used(edges.size(), false);
    const double areaEps = tol.scale(std::max(1.0, backgroundArea(box)));
    for (std::size_t startEdge = 0; startEdge < edges.size(); ++startEdge) {
        if (used[startEdge]) continue;
        std::vector<Point2D> loop;
        std::size_t edgeId = startEdge;
        const std::size_t startNode = edges[startEdge].from;
        std::size_t guard = 0;
        while (guard++ <= edges.size()) {
            if (used[edgeId]) {
                result.graphInvalid = true;
                return result;
            }
            used[edgeId] = true;
            const auto edge = edges[edgeId];
            loop.push_back(points[edge.from]);
            if (edge.to == startNode) break;
            if (outgoing[edge.to].size() != 1) {
                result.graphInvalid = true;
                return result;
            }
            edgeId = outgoing[edge.to].front();
        }
        if (guard > edges.size() + 1) {
            result.graphInvalid = true;
            return result;
        }
        removeCollinearVertices(loop, tol);
        if (loop.size() < 3) continue;
        Polygon2D polygon{loop};
        if (polygon.area() <= areaEps) continue;
        if (polygon.signedArea() < 0.0) std::reverse(polygon.vertices.begin(), polygon.vertices.end());
        BoundaryLoop diagnosticsLoop(polygon.vertices);
        if (!diagnosticsLoop.diagnose(tol).valid()) {
            result.graphInvalid = true;
            return result;
        }
        result.components.push_back(std::move(polygon));
    }
    return result;
}

} // namespace

CutCell2D buildCutCell(const AABB2D& box, CellClass classification,
                       const BoundaryLoop& inputBoundary,
                       const TolerancePolicy& tol) {
    CutCell2D result;
    result.backgroundBounds = box;

    const double fullArea = backgroundArea(box);
    if (!(fullArea > tol.scale(std::max(1.0, std::abs(fullArea))))) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::InvalidBackgroundCell,
                                 "background cell has non-positive area"});
        return result;
    }

    const auto inputDiagnostics = inputBoundary.diagnose(tol);
    if (!inputDiagnostics.valid()) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::InvalidBoundary,
                                 "boundary loop failed geometry diagnostics"});
        return result;
    }

    BoundaryLoop boundary = inputBoundary;
    if (!boundary.normalizeCounterClockwise(tol)) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::InvalidBoundary,
                                 "boundary loop could not be normalized to CCW"});
        return result;
    }

    if (classification == CellClass::Outside) {
        result.kind = CutCellKind::Empty;
        return result;
    }
    if (classification == CellClass::Inside) {
        result.kind = CutCellKind::Full;
        result.fluidPolygon = rectanglePolygon(box);
        result.area = fullArea;
        result.areaFraction = 1.0;
        result.centroid = result.fluidPolygon.centroid(tol);
        return result;
    }

    const auto local = buildLocalIntersection(boundary, box, tol);
    result.embeddedBoundary = local.embedded;
    if (local.graphInvalid) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::DegeneratePolygon,
                                 "local Cut-cell boundary graph is open, branched or degenerate"});
        return result;
    }

    if (local.components.empty()) {
        const Point2D center{0.5 * (box.min.x + box.max.x), 0.5 * (box.min.y + box.max.y)};
        const auto state = classifyPointInPolygon(center, boundary.polygon(), tol);
        if (state == PointInPolygon::Inside || state == PointInPolygon::Boundary) {
            result.kind = CutCellKind::Full;
            result.fluidPolygon = rectanglePolygon(box);
            result.area = fullArea;
            result.areaFraction = 1.0;
            result.centroid = result.fluidPolygon.centroid(tol);
            result.embeddedBoundary.clear();
        } else {
            result.kind = CutCellKind::Empty;
        }
        return result;
    }

    if (local.components.size() > 1) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::MultipleEmbeddedComponents,
                                 "multiple disconnected fluid components occur in one leaf"});
        return result;
    }

    result.fluidPolygon = local.components.front();
    result.area = result.fluidPolygon.area();
    const double areaEps = tol.scale(std::max(1.0, fullArea));
    if (result.area <= areaEps) {
        result.kind = CutCellKind::Empty;
        result.fluidPolygon.vertices.clear();
        result.area = 0.0;
        result.areaFraction = 0.0;
        result.centroid.reset();
        return result;
    }
    if (result.area > fullArea + areaEps) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::AreaOutOfRange,
                                 "local Cut-cell fluid area exceeds background-cell area"});
        return result;
    }

    result.areaFraction = std::clamp(result.area / fullArea, 0.0, 1.0);
    result.centroid = result.fluidPolygon.centroid(tol);
    if (!result.centroid) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::DegeneratePolygon,
                                 "local Cut-cell polygon has no valid centroid"});
        return result;
    }

    if (fullArea - result.area <= areaEps) {
        result.kind = CutCellKind::Full;
        result.fluidPolygon = rectanglePolygon(box);
        result.area = fullArea;
        result.areaFraction = 1.0;
        result.centroid = result.fluidPolygon.centroid(tol);
        result.embeddedBoundary.clear();
        return result;
    }

    const std::size_t componentCount = countFragmentComponents(result.embeddedBoundary, tol);
    if (componentCount == 0) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::MissingEmbeddedBoundary,
                                 "partial fluid polygon has no embedded-boundary fragment"});
        return result;
    }

    result.kind = CutCellKind::Cut;
    return result;
}

CutCell2D buildCutCell(const QuadtreeLeaf2D& leaf,
                       const BoundaryLoop& boundary,
                       const TolerancePolicy& tol) {
    CutCell2D result = buildCutCell(leaf.bounds, leaf.classification, boundary, tol);
    result.sourceId = leaf.id;
    result.sourceKey = leaf.key;
    return result;
}

} // namespace cartmesh2d
