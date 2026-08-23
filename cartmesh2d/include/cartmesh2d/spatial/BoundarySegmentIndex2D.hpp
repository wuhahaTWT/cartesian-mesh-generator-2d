#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <vector>

namespace cartmesh2d {

// Deterministic AABB hierarchy over one closed boundary loop.  The index owns
// a copy of the input segments so repeated Cartesian-cell queries never rescan
// the complete boundary.
class BoundarySegmentIndex2D {
public:
    explicit BoundarySegmentIndex2D(const BoundaryLoop& boundary,
                                    const TolerancePolicy& tol = {});

    [[nodiscard]] bool valid() const noexcept { return valid_; }
    [[nodiscard]] std::size_t segmentCount() const noexcept { return segments_.size(); }
    [[nodiscard]] bool matches(const BoundaryLoop& boundary) const noexcept;

    [[nodiscard]] std::vector<std::size_t> querySegmentIds(
        const AABB2D& box, const TolerancePolicy& tol = {}) const;
    [[nodiscard]] bool intersects(const AABB2D& box,
                                  const TolerancePolicy& tol = {}) const;
    [[nodiscard]] PointInPolygon classifyPoint(
        const Point2D& point, const TolerancePolicy& tol = {}) const;
    [[nodiscard]] double distanceToAABB(
        const AABB2D& box, const TolerancePolicy& tol = {}) const;

private:
    struct SegmentRecord {
        Segment2D segment;
        AABB2D bounds;
        std::size_t sourceId = 0;
    };

    struct Node {
        AABB2D bounds;
        std::size_t left = 0;
        std::size_t right = 0;
        std::size_t begin = 0;
        std::size_t end = 0;
        bool leaf = false;
    };

    [[nodiscard]] std::size_t buildNode(std::size_t begin, std::size_t end);

    std::vector<Point2D> boundaryVertices_;
    std::vector<SegmentRecord> segments_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
    AABB2D bounds_{};
    bool valid_ = false;
};

} // namespace cartmesh2d
