#pragma once
#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"
#include <map>
#include <vector>

namespace cartmesh2d {

// One immutable active-vertex snapshot per topology transaction. A support
// interval is partitioned once and both incident cells consume the same list.
// Unused registry vertices (e.g. from a rejected candidate) are not included.
class SharedEdgePartition2D {
public:
    SharedEdgePartition2D(const IntersectionRegistry2D& registry,
                         std::vector<std::size_t> activeHandles);
    [[nodiscard]] const std::vector<std::size_t>& handles() const { return handles_; }
    [[nodiscard]] const std::vector<Point2D>& points() const { return points_; }
    [[nodiscard]] std::size_t denseId(std::size_t handle) const { return denseIds_.at(handle); }
    [[nodiscard]] const std::vector<std::pair<double,std::size_t>>& partition(
        std::size_t a, std::size_t b, double epsilon,const TolerancePolicy& tol);
    [[nodiscard]] std::size_t partitionCount() const { return partitions_.size(); }
    [[nodiscard]] std::size_t cacheHits() const { return cacheHits_; }
private:
    std::vector<std::size_t> handles_;
    std::vector<Point2D> points_;
    std::map<std::size_t,std::size_t> denseIds_;
    std::map<double,std::map<double,std::size_t>> columns_;
    std::map<std::pair<std::size_t,std::size_t>,std::vector<std::pair<double,std::size_t>>> partitions_;
    std::size_t cacheHits_=0;
};
}
