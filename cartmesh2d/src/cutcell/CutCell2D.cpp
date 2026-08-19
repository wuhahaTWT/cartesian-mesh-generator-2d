#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <algorithm>
#include <cmath>
#include <queue>

namespace cartmesh2d {
namespace {

enum class ClipSide { Left, Right, Bottom, Top };

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

[[nodiscard]] bool insideClipSide(const Point2D& point, ClipSide side,
                                  const AABB2D& box,
                                  const TolerancePolicy& tol) noexcept {
    switch (side) {
    case ClipSide::Left:
        return point.x >= box.min.x - tol.scale(std::abs(box.min.x));
    case ClipSide::Right:
        return point.x <= box.max.x + tol.scale(std::abs(box.max.x));
    case ClipSide::Bottom:
        return point.y >= box.min.y - tol.scale(std::abs(box.min.y));
    case ClipSide::Top:
        return point.y <= box.max.y + tol.scale(std::abs(box.max.y));
    }
    return false;
}

[[nodiscard]] Point2D intersectClipSide(const Point2D& a, const Point2D& b,
                                        ClipSide side, const AABB2D& box) noexcept {
    if (side == ClipSide::Left || side == ClipSide::Right) {
        const double x = side == ClipSide::Left ? box.min.x : box.max.x;
        const double dx = b.x - a.x;
        const double t = std::abs(dx) > 0.0 ? (x - a.x) / dx : 0.0;
        return {x, a.y + t * (b.y - a.y)};
    }
    const double y = side == ClipSide::Bottom ? box.min.y : box.max.y;
    const double dy = b.y - a.y;
    const double t = std::abs(dy) > 0.0 ? (y - a.y) / dy : 0.0;
    return {a.x + t * (b.x - a.x), y};
}

void removeConsecutiveDuplicates(std::vector<Point2D>& vertices,
                                 const TolerancePolicy& tol) {
    std::vector<Point2D> result;
    result.reserve(vertices.size());
    for (const auto& point : vertices) {
        if (result.empty() || !pointNear(result.back(), point, tol)) {
            result.push_back(point);
        }
    }
    if (result.size() > 1 && pointNear(result.front(), result.back(), tol)) {
        result.pop_back();
    }
    vertices.swap(result);
}

[[nodiscard]] std::vector<Point2D> clipAgainstSide(
    const std::vector<Point2D>& input, ClipSide side,
    const AABB2D& box, const TolerancePolicy& tol) {
    std::vector<Point2D> output;
    if (input.empty()) return output;

    Point2D previous = input.back();
    bool previousInside = insideClipSide(previous, side, box, tol);
    for (const auto& current : input) {
        const bool currentInside = insideClipSide(current, side, box, tol);
        if (currentInside) {
            if (!previousInside) {
                output.push_back(intersectClipSide(previous, current, side, box));
            }
            output.push_back(current);
        } else if (previousInside) {
            output.push_back(intersectClipSide(previous, current, side, box));
        }
        previous = current;
        previousInside = currentInside;
    }
    removeConsecutiveDuplicates(output, tol);
    return output;
}

[[nodiscard]] Polygon2D clipPolygonToAABB(const Polygon2D& polygon,
                                          const AABB2D& box,
                                          const TolerancePolicy& tol) {
    std::vector<Point2D> vertices = polygon.vertices;
    for (const auto side : {ClipSide::Left, ClipSide::Right,
                            ClipSide::Bottom, ClipSide::Top}) {
        vertices = clipAgainstSide(vertices, side, box, tol);
        if (vertices.empty()) break;
    }
    removeConsecutiveDuplicates(vertices, tol);
    Polygon2D result{vertices};
    if (result.signedArea() < 0.0) {
        std::reverse(result.vertices.begin(), result.vertices.end());
    }
    return result;
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
        if (p < 0.0) {
            t0 = std::max(t0, ratio);
        } else {
            t1 = std::min(t1, ratio);
        }
        return t0 <= t1 + eps;
    };

    if (!update(-dx, segment.a.x - box.min.x) ||
        !update(dx, box.max.x - segment.a.x) ||
        !update(-dy, segment.a.y - box.min.y) ||
        !update(dy, box.max.y - segment.a.y)) {
        return std::nullopt;
    }

    t0 = std::clamp(t0, 0.0, 1.0);
    t1 = std::clamp(t1, 0.0, 1.0);
    Segment2D clipped{{segment.a.x + t0 * dx, segment.a.y + t0 * dy},
                      {segment.a.x + t1 * dx, segment.a.y + t1 * dy}};
    if (pointNear(clipped.a, clipped.b, tol)) return std::nullopt;
    return clipped;
}

[[nodiscard]] bool liesOnBoxPerimeter(const Segment2D& segment,
                                      const AABB2D& box,
                                      const TolerancePolicy& tol) noexcept {
    return (scalarNear(segment.a.x, box.min.x, tol) &&
            scalarNear(segment.b.x, box.min.x, tol)) ||
           (scalarNear(segment.a.x, box.max.x, tol) &&
            scalarNear(segment.b.x, box.max.x, tol)) ||
           (scalarNear(segment.a.y, box.min.y, tol) &&
            scalarNear(segment.b.y, box.min.y, tol)) ||
           (scalarNear(segment.a.y, box.max.y, tol) &&
            scalarNear(segment.b.y, box.max.y, tol));
}

[[nodiscard]] std::vector<Segment2D> collectEmbeddedBoundary(
    const BoundaryLoop& boundary, const AABB2D& box,
    const TolerancePolicy& tol) {
    std::vector<Segment2D> fragments;
    const auto& vertices = boundary.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Segment2D edge{vertices[i], vertices[(i + 1) % vertices.size()]};
        const auto clipped = clipSegmentToAABB(edge, box, tol);
        if (clipped && !liesOnBoxPerimeter(*clipped, box, tol)) {
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

    result.fluidPolygon = clipPolygonToAABB(boundary.polygon(), box, tol);
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
                                 "clipped fluid area exceeds background-cell area"});
        return result;
    }

    result.areaFraction = std::clamp(result.area / fullArea, 0.0, 1.0);
    result.centroid = result.fluidPolygon.centroid(tol);
    if (!result.centroid) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::DegeneratePolygon,
                                 "clipped fluid polygon has no valid centroid"});
        return result;
    }

    if (fullArea - result.area <= areaEps) {
        result.kind = CutCellKind::Full;
        result.fluidPolygon = rectanglePolygon(box);
        result.area = fullArea;
        result.areaFraction = 1.0;
        result.centroid = result.fluidPolygon.centroid(tol);
        return result;
    }

    BoundaryLoop clippedLoop(result.fluidPolygon.vertices);
    const auto clippedDiagnostics = clippedLoop.diagnose(tol);
    if (!clippedDiagnostics.valid()) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::DegeneratePolygon,
                                 "clipped fluid polygon is degenerate or self-intersecting"});
        return result;
    }

    result.embeddedBoundary = collectEmbeddedBoundary(boundary, box, tol);
    const std::size_t componentCount = countFragmentComponents(result.embeddedBoundary, tol);
    if (componentCount == 0) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::MissingEmbeddedBoundary,
                                 "partial fluid polygon has no embedded-boundary fragment"});
        return result;
    }
    if (componentCount > 1) {
        result.kind = CutCellKind::Unsupported;
        result.issues.push_back({CutCellIssueCode::MultipleEmbeddedComponents,
                                 "multiple disconnected embedded-boundary components occur in one leaf"});
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
