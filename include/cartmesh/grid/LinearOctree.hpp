#pragma once

#include "cartmesh/geometry/AABB.hpp"
#include "cartmesh/grid/OctreeNodeCode.hpp"

#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

namespace cartmesh {

enum class FaceDirection : std::uint8_t {
    negative_x = 0,
    positive_x = 1,
    negative_y = 2,
    positive_y = 3,
    negative_z = 4,
    positive_z = 5,
};

struct OctreeLevelStatistics {
    std::vector<std::uint64_t> leaf_count_by_level;
    std::uint8_t minimum_leaf_level{};
    std::uint8_t maximum_leaf_level{};
};

struct OctreeRefinementStatistics {
    std::uint64_t initial_leaf_count{};
    std::uint64_t final_leaf_count{};
    std::uint64_t split_count{};
    std::uint64_t iteration_count{};
};

struct OctreeBalanceReport {
    bool balanced{};
    std::uint64_t violating_face_pair_count{};
    std::uint8_t maximum_level_difference{};
};

struct OctreeBalanceStatistics {
    std::uint64_t initial_leaf_count{};
    std::uint64_t final_leaf_count{};
    std::uint64_t split_count{};
    std::uint64_t iteration_count{};
};

class LinearOctree {
  public:
    using DesiredLevelFunction =
        std::function<std::uint8_t(OctreeNodeCode, const AABB&)>;

    LinearOctree(AABB domain, std::uint8_t base_level, std::uint8_t maximum_level);

    [[nodiscard]] constexpr const AABB& domain() const noexcept { return domain_; }
    [[nodiscard]] constexpr std::uint8_t base_level() const noexcept { return base_level_; }
    [[nodiscard]] constexpr std::uint8_t maximum_level() const noexcept {
        return maximum_level_;
    }
    [[nodiscard]] std::uint64_t leaf_count() const noexcept { return leaves_.size(); }
    [[nodiscard]] std::span<const OctreeNodeCode> leaf_codes() const noexcept {
        return leaves_;
    }
    [[nodiscard]] OctreeNodeCode leaf_code(std::uint64_t leaf_id) const;
    [[nodiscard]] AABB cell_bounds(OctreeNodeCode code) const;
    [[nodiscard]] Vec3 cell_center(OctreeNodeCode code) const;
    [[nodiscard]] double cell_volume(OctreeNodeCode code) const;
    [[nodiscard]] std::optional<std::uint64_t> find_leaf(OctreeNodeCode code) const;
    [[nodiscard]] std::optional<std::uint64_t>
    find_leaf_covering_maximum_level_cell(std::uint32_t x, std::uint32_t y,
                                          std::uint32_t z) const;
    [[nodiscard]] std::vector<OctreeNodeCode> face_neighbors(OctreeNodeCode code,
                                                            FaceDirection direction) const;

    bool refine_leaf(OctreeNodeCode code);
    bool coarsen_parent(OctreeNodeCode parent_code);
    [[nodiscard]] OctreeRefinementStatistics
    refine_to_desired_levels(const DesiredLevelFunction& desired_level);
    [[nodiscard]] OctreeBalanceReport check_face_balance() const;
    [[nodiscard]] OctreeBalanceStatistics balance_faces_2_to_1();

    [[nodiscard]] OctreeLevelStatistics level_statistics() const;
    [[nodiscard]] bool validate_partition() const;
    [[nodiscard]] std::uint64_t compact_storage_bytes() const noexcept {
        return static_cast<std::uint64_t>(leaves_.capacity()) * sizeof(OctreeNodeCode);
    }
    [[nodiscard]] std::uint64_t result_hash_fnv1a64() const noexcept;

  private:
    [[nodiscard]] MortonCode anchor(OctreeNodeCode code) const;
    [[nodiscard]] std::uint64_t span(OctreeNodeCode code) const;
    [[nodiscard]] bool code_is_within_levels(OctreeNodeCode code) const;
    void replace_marked_with_children(const std::vector<OctreeNodeCode>& marked);

    AABB domain_;
    std::uint8_t base_level_;
    std::uint8_t maximum_level_;
    std::vector<OctreeNodeCode> leaves_;
};

} // 命名空间 cartmesh
