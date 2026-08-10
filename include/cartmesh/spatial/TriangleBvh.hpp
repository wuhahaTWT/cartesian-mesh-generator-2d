#pragma once

#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cartmesh {

struct BvhStatistics {
    std::uint64_t triangle_count{};
    std::uint64_t node_count{};
    std::uint64_t leaf_count{};
    std::uint64_t maximum_depth{};
    std::uint64_t maximum_leaf_triangles{};
};

struct RayIntersectionSummary {
    std::uint64_t raw_triangle_hits{};
    std::uint64_t unique_surface_crossings{};
    bool ambiguous{};
};

class TriangleBvh {
  public:
    explicit TriangleBvh(const SurfaceMesh& surface, std::size_t leaf_size = 8);

    [[nodiscard]] const std::vector<Triangle>& triangles() const noexcept { return triangles_; }
    [[nodiscard]] const AABB& bounds() const noexcept { return nodes_.front().bounds; }
    [[nodiscard]] const BvhStatistics& statistics() const noexcept { return statistics_; }
    [[nodiscard]] std::vector<std::uint64_t> query(const AABB& bounds) const;
    [[nodiscard]] bool intersects_surface(const AABB& bounds) const;
    [[nodiscard]] RayIntersectionSummary trace_ray(const Vec3& origin,
                                                   const Vec3& direction) const;
    [[nodiscard]] bool point_on_surface(const Vec3& point, double tolerance) const;
    [[nodiscard]] double distance_to_surface(const Vec3& point) const;

  private:
    struct Node {
        explicit Node(AABB node_bounds) : bounds(node_bounds) {}

        AABB bounds;
        std::uint64_t first{};
        std::uint32_t count{};
        std::uint32_t left{};
        std::uint32_t right{};
        bool leaf{};
    };

    [[nodiscard]] std::uint32_t build(std::size_t begin, std::size_t end,
                                      std::uint64_t depth, std::size_t leaf_size);

    std::vector<Triangle> triangles_;
    std::vector<std::uint64_t> indices_;
    std::vector<Node> nodes_;
    BvhStatistics statistics_;
};

} // 命名空间 cartmesh
