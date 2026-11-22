#pragma once
#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
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
    [[nodiscard]] std::size_t denseId(std::size_t handle) const {
        const auto id=denseIds_.at(handle);
        if (id==std::numeric_limits<std::size_t>::max())
            throw std::out_of_range("inactive canonical vertex handle");
        return id;
    }
    [[nodiscard]] const std::vector<std::pair<double,std::size_t>>& partition(
        std::size_t a, std::size_t b, double epsilon,const TolerancePolicy& tol);
    [[nodiscard]] std::size_t partitionCount() const { return partitions_.size(); }
    [[nodiscard]] std::size_t cacheHits() const { return cacheHits_; }
private:
    struct PairHash {
        [[nodiscard]] std::size_t operator()(const std::pair<std::size_t,std::size_t>& p) const noexcept {
            std::size_t h=p.first*0x9e3779b97f4a7c15ULL;
            h^=p.second+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);
            return h;
        }
    };
    std::vector<std::size_t> handles_;
    std::vector<Point2D> points_;
    std::vector<std::size_t> denseIds_;
    struct Column {
        double x;
        std::vector<std::pair<double,std::size_t>> byY;
    };
    std::vector<Column> columns_;
    std::unordered_map<std::pair<std::size_t,std::size_t>,
                       std::vector<std::pair<double,std::size_t>>,PairHash> partitions_;
    std::size_t cacheHits_=0;
};
}
