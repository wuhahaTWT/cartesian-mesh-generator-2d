#pragma once

#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace cartmesh {

enum class PointClassification : std::uint8_t {
    outside = 0,
    inside = 1,
    on_surface = 2,
    conflict = 3,
};

enum class CellClassification : std::uint8_t {
    outside = 0,
    inside = 1,
    intersected = 2,
    conflict = 3,
};

struct PointClassificationResult {
    PointClassification classification{PointClassification::conflict};
    std::uint8_t inside_votes{};
    std::uint8_t outside_votes{};
    std::uint8_t ambiguous_rays{};
};

class SurfaceClassifier {
  public:
    explicit SurfaceClassifier(const TriangleBvh& bvh, double surface_tolerance = 0.0);

    [[nodiscard]] PointClassificationResult classify(const Vec3& point) const;
    [[nodiscard]] CellClassification classify_cell(const AABB& bounds) const;
    [[nodiscard]] CellClassification classify_cell(
        const AABB& bounds, PointClassification center_classification) const;
    [[nodiscard]] double surface_tolerance() const noexcept { return surface_tolerance_; }

  private:
    const TriangleBvh& bvh_;
    double surface_tolerance_;
    std::array<Vec3, 3> directions_;
};

struct UniformClassification {
    std::vector<std::uint8_t> cell_classification;
    std::vector<std::uint8_t> center_point_classification;
    std::uint64_t inside_count{};
    std::uint64_t outside_count{};
    std::uint64_t intersected_count{};
    std::uint64_t conflict_count{};
    std::uint64_t center_inside_count{};
    std::uint64_t center_outside_count{};
    std::uint64_t center_on_surface_count{};
    std::uint64_t center_conflict_count{};
};

struct AdaptiveClassification {
    std::vector<std::uint8_t> cell_classification;
    std::vector<std::uint8_t> center_point_classification;
    std::uint64_t inside_count{};
    std::uint64_t outside_count{};
    std::uint64_t intersected_count{};
    std::uint64_t conflict_count{};
    std::uint64_t center_inside_count{};
    std::uint64_t center_outside_count{};
    std::uint64_t center_on_surface_count{};
    std::uint64_t center_conflict_count{};
    double inside_volume{};
    double inside_plus_intersected_volume{};
};

[[nodiscard]] UniformClassification classify_uniform_cells(const UniformCartesianGrid& grid,
                                                            const SurfaceClassifier& classifier);

[[nodiscard]] AdaptiveClassification classify_octree_leaves(
    const LinearOctree& tree, const SurfaceClassifier& classifier);

[[nodiscard]] constexpr const char* point_classification_name(
    PointClassification classification) noexcept {
    switch (classification) {
    case PointClassification::outside:
        return "outside";
    case PointClassification::inside:
        return "inside";
    case PointClassification::on_surface:
        return "on_surface";
    case PointClassification::conflict:
        return "conflict";
    }
    return "conflict";
}

[[nodiscard]] constexpr const char* cell_classification_name(
    CellClassification classification) noexcept {
    switch (classification) {
    case CellClassification::outside:
        return "outside";
    case CellClassification::inside:
        return "inside";
    case CellClassification::intersected:
        return "intersected";
    case CellClassification::conflict:
        return "conflict";
    }
    return "conflict";
}

} // 命名空间 cartmesh
