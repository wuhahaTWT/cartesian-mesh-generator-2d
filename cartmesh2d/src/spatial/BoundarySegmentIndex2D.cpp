#include "cartmesh2d/spatial/BoundarySegmentIndex2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

[[nodiscard]] AABB2D segmentBounds(const Segment2D& segment) noexcept {
    return {{std::min(segment.a.x, segment.b.x), std::min(segment.a.y, segment.b.y)},
            {std::max(segment.a.x, segment.b.x), std::max(segment.a.y, segment.b.y)}};
}

[[nodiscard]] AABB2D mergedBounds(const AABB2D& a, const AABB2D& b) noexcept {
    return {{std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y)},
            {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y)}};
}

[[nodiscard]] double boxScale(const AABB2D& a, const AABB2D& b) noexcept {
    return std::max({a.max.x - a.min.x, a.max.y - a.min.y,
                     b.max.x - b.min.x, b.max.y - b.min.y});
}

[[nodiscard]] bool boxesOverlap(const AABB2D& a, const AABB2D& b,
                                const TolerancePolicy& tol) noexcept {
    const double eps = tol.scale(boxScale(a, b));
    return a.max.x >= b.min.x - eps && b.max.x >= a.min.x - eps &&
           a.max.y >= b.min.y - eps && b.max.y >= a.min.y - eps;
}

[[nodiscard]] bool segmentIntersectsBox(const Segment2D& segment,
                                        const AABB2D& box,
                                        const TolerancePolicy& tol) noexcept {
    if (box.contains(segment.a, tol) || box.contains(segment.b, tol)) return true;
    const Point2D bl{box.min.x, box.min.y};
    const Point2D br{box.max.x, box.min.y};
    const Point2D tr{box.max.x, box.max.y};
    const Point2D tl{box.min.x, box.max.y};
    const Segment2D sides[4]{{bl, br}, {br, tr}, {tr, tl}, {tl, bl}};
    for (const auto& side : sides) {
        if (intersectSegments(segment, side, tol).kind != SegmentIntersectionKind::None) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] double pointSegmentDistance(const Point2D& point,
                                          const Segment2D& segment) noexcept {
    const Vector2D direction = segment.b - segment.a;
    const double denominator = squaredNorm(direction);
    if (denominator <= 0.0) return std::sqrt(squaredNorm(point - segment.a));
    const double parameter = std::clamp(
        dot(point - segment.a, direction) / denominator, 0.0, 1.0);
    return std::sqrt(squaredNorm(point - (segment.a + direction * parameter)));
}

[[nodiscard]] double pointBoxDistance(const Point2D& point,
                                      const AABB2D& box) noexcept {
    const double dx = point.x < box.min.x ? box.min.x - point.x
                    : point.x > box.max.x ? point.x - box.max.x : 0.0;
    const double dy = point.y < box.min.y ? box.min.y - point.y
                    : point.y > box.max.y ? point.y - box.max.y : 0.0;
    return std::hypot(dx, dy);
}

[[nodiscard]] double boxBoxDistance(const AABB2D& a, const AABB2D& b) noexcept {
    const double dx = a.max.x < b.min.x ? b.min.x - a.max.x
                    : b.max.x < a.min.x ? a.min.x - b.max.x : 0.0;
    const double dy = a.max.y < b.min.y ? b.min.y - a.max.y
                    : b.max.y < a.min.y ? a.min.y - b.max.y : 0.0;
    return std::hypot(dx, dy);
}

[[nodiscard]] double segmentBoxDistance(const Segment2D& segment,
                                        const AABB2D& box,
                                        const TolerancePolicy& tol) noexcept {
    if (segmentIntersectsBox(segment, box, tol)) return 0.0;
    const Point2D corners[4]{{box.min.x, box.min.y}, {box.max.x, box.min.y},
                             {box.max.x, box.max.y}, {box.min.x, box.max.y}};
    double best = std::min(pointBoxDistance(segment.a, box),
                           pointBoxDistance(segment.b, box));
    for (const auto& corner : corners) {
        best = std::min(best, pointSegmentDistance(corner, segment));
    }
    return best;
}

} // namespace

BoundarySegmentIndex2D::BoundarySegmentIndex2D(
    const BoundaryLoop& boundary, const TolerancePolicy& tol)
    : BoundarySegmentIndex2D(BoundaryRegion2D(boundary), tol) {}

BoundarySegmentIndex2D::BoundarySegmentIndex2D(
    const BoundaryRegion2D& boundary, const TolerancePolicy& tol) {
    if (!boundary.diagnose(tol).valid()) return;
    std::size_t segmentCount = 0;
    for (const auto& loop : boundary.loops()) segmentCount += loop.vertices().size();
    segments_.reserve(segmentCount);
    boundaryLoops_.reserve(boundary.loops().size());
    std::size_t sourceId = 0;
    for (const auto& loop : boundary.loops()) {
        boundaryLoops_.push_back(loop.vertices());
        const auto& vertices = boundaryLoops_.back();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const Segment2D segment{vertices[i], vertices[(i + 1U) % vertices.size()]};
            segments_.push_back({segment, segmentBounds(segment), sourceId++});
        }
    }
    order_.resize(segments_.size());
    std::iota(order_.begin(), order_.end(), 0U);
    nodes_.reserve(2U * segments_.size());
    if (!segments_.empty()) {
        const std::size_t root = buildNode(0U, segments_.size());
        (void)root;
        bounds_ = nodes_.front().bounds;
    }
    valid_ = !segments_.empty();
}

std::size_t BoundarySegmentIndex2D::buildNode(std::size_t begin, std::size_t end) {
    const std::size_t nodeId = nodes_.size();
    nodes_.push_back({});
    AABB2D bounds = segments_[order_[begin]].bounds;
    for (std::size_t i = begin + 1U; i < end; ++i) {
        bounds = mergedBounds(bounds, segments_[order_[i]].bounds);
    }
    nodes_[nodeId].bounds = bounds;
    nodes_[nodeId].begin = begin;
    nodes_[nodeId].end = end;
    constexpr std::size_t leafCapacity = 8U;
    if (end - begin <= leafCapacity) {
        nodes_[nodeId].leaf = true;
        return nodeId;
    }

    const bool splitX = (bounds.max.x - bounds.min.x) >=
                        (bounds.max.y - bounds.min.y);
    std::stable_sort(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                     order_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t lhs, std::size_t rhs) {
        const auto& a = segments_[lhs].bounds;
        const auto& b = segments_[rhs].bounds;
        const double centerA = splitX ? a.min.x + a.max.x : a.min.y + a.max.y;
        const double centerB = splitX ? b.min.x + b.max.x : b.min.y + b.max.y;
        if (centerA != centerB) return centerA < centerB;
        return segments_[lhs].sourceId < segments_[rhs].sourceId;
    });
    const std::size_t middle = begin + (end - begin) / 2U;
    const std::size_t left = buildNode(begin, middle);
    const std::size_t right = buildNode(middle, end);
    nodes_[nodeId].left = left;
    nodes_[nodeId].right = right;
    return nodeId;
}

bool BoundarySegmentIndex2D::matches(const BoundaryLoop& boundary) const {
    return matches(BoundaryRegion2D(boundary));
}

bool BoundarySegmentIndex2D::matches(const BoundaryRegion2D& boundary) const noexcept {
    if (boundary.loops().size() != boundaryLoops_.size()) return false;
    for (std::size_t loopId = 0; loopId < boundaryLoops_.size(); ++loopId) {
        const auto& vertices = boundary.loops()[loopId].vertices();
        if (vertices.size() != boundaryLoops_[loopId].size()) return false;
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (vertices[i].x != boundaryLoops_[loopId][i].x ||
                vertices[i].y != boundaryLoops_[loopId][i].y) return false;
        }
    }
    return true;
}

std::vector<std::size_t> BoundarySegmentIndex2D::querySegmentIds(
    const AABB2D& box, const TolerancePolicy& tol) const {
    std::vector<std::size_t> result;
    if (!valid_ || !box.valid(tol)) return result;
    std::vector<std::size_t> stack{0U};
    while (!stack.empty()) {
        const std::size_t nodeId = stack.back();
        stack.pop_back();
        const auto& node = nodes_[nodeId];
        if (!boxesOverlap(node.bounds, box, tol)) continue;
        if (node.leaf) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                const auto& record = segments_[order_[i]];
                if (boxesOverlap(record.bounds, box, tol)) result.push_back(record.sourceId);
            }
        } else {
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool BoundarySegmentIndex2D::intersects(const AABB2D& box,
                                        const TolerancePolicy& tol) const {
    if (!valid_ || !box.valid(tol)) return false;
    std::vector<std::size_t> stack{0U};
    while (!stack.empty()) {
        const std::size_t nodeId = stack.back();
        stack.pop_back();
        const auto& node = nodes_[nodeId];
        if (!boxesOverlap(node.bounds, box, tol)) continue;
        if (node.leaf) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                const auto& record = segments_[order_[i]];
                if (boxesOverlap(record.bounds, box, tol) &&
                    segmentIntersectsBox(record.segment, box, tol)) return true;
            }
        } else {
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }
    return false;
}

PointInPolygon BoundarySegmentIndex2D::classifyPoint(
    const Point2D& point, const TolerancePolicy& tol) const {
    if (!valid_) return PointInPolygon::Outside;
    const double scale = std::max(bounds_.max.x - bounds_.min.x,
                                  bounds_.max.y - bounds_.min.y);
    const double eps = tol.scale(scale);
    const AABB2D horizontalSlab{{bounds_.min.x - eps, point.y - eps},
                                {bounds_.max.x + eps, point.y + eps}};
    const auto candidates = querySegmentIds(horizontalSlab, tol);
    int winding = 0;
    for (const std::size_t id : candidates) {
        const auto& segment = segments_[id].segment;
        if (pointOnSegment(point, segment, tol)) return PointInPolygon::Boundary;
        const Point2D& a = segment.a;
        const Point2D& b = segment.b;
        if (a.y <= point.y) {
            if (b.y > point.y && orientationSign(a, b, point) > 0) ++winding;
        } else if (b.y <= point.y && orientationSign(a, b, point) < 0) {
            --winding;
        }
    }
    return winding % 2 == 0 ? PointInPolygon::Outside : PointInPolygon::Inside;
}

double BoundarySegmentIndex2D::distanceToAABB(
    const AABB2D& box, const TolerancePolicy& tol) const {
    if (!valid_ || !box.valid(tol)) return std::numeric_limits<double>::infinity();
    using QueueItem = std::pair<double, std::size_t>;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>> queue;
    queue.push({boxBoxDistance(box, nodes_.front().bounds), 0U});
    double best = std::numeric_limits<double>::infinity();
    while (!queue.empty()) {
        const auto [lowerBound, nodeId] = queue.top();
        queue.pop();
        if (lowerBound >= best) continue;
        const auto& node = nodes_[nodeId];
        if (node.leaf) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                best = std::min(best,
                    segmentBoxDistance(segments_[order_[i]].segment, box, tol));
            }
        } else {
            const double leftDistance = boxBoxDistance(box, nodes_[node.left].bounds);
            const double rightDistance = boxBoxDistance(box, nodes_[node.right].bounds);
            if (leftDistance < best) queue.push({leftDistance, node.left});
            if (rightDistance < best) queue.push({rightDistance, node.right});
        }
    }
    return best;
}

} // namespace cartmesh2d
