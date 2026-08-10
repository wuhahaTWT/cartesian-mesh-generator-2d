#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/grid/LinearOctree.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace cartmesh {

enum class IncrementalTopologyRelation : std::uint8_t {
    preserved = 0,
    rebuilt = 1,
    refined = 2,
    coarsened = 3,
};

struct IncrementalCellMappingEntry {
    OctreeNodeCode old_code{};
    OctreeNodeCode new_code{};
    IncrementalTopologyRelation relation{IncrementalTopologyRelation::rebuilt};
    double background_overlap_volume{};
    double old_fluid_volume{};
    double new_fluid_volume{};
    double shared_fluid_overlap_volume{};
    double old_volume_preserved_fraction{};
    double new_volume_from_old_fraction{};
    bool old_cell_is_fluid{};
    bool new_cell_is_fluid{};
    bool exact_fluid_overlap{};
};

struct IncrementalMeshMapping {
    std::vector<IncrementalCellMappingEntry> entries;
    double total_background_overlap_volume{};
    double total_shared_fluid_overlap_volume{};
    double old_fluid_volume{};
    double new_fluid_volume{};
    double removed_fluid_volume{};
    double created_fluid_volume{};
    std::uint64_t preserved_pair_count{};
    std::uint64_t rebuilt_pair_count{};
    std::uint64_t refined_pair_count{};
    std::uint64_t coarsened_pair_count{};
};

[[nodiscard]] constexpr const char* incremental_topology_relation_name(
    IncrementalTopologyRelation relation) noexcept {
    switch (relation) {
    case IncrementalTopologyRelation::preserved: return "preserved";
    case IncrementalTopologyRelation::rebuilt: return "rebuilt";
    case IncrementalTopologyRelation::refined: return "refined";
    case IncrementalTopologyRelation::coarsened: return "coarsened";
    }
    return "unknown";
}

// 线性八叉树分区在 Morton 区间上做确定性双指针重叠。Cut-cell 流体重叠
// 使用凸片半空间交集计算；新出现的流体体积和被几何吞没的旧流体体积
// 分开报告，不能被伪造成守恒一对一权重。
[[nodiscard]] IncrementalMeshMapping build_incremental_mesh_mapping(
    const LinearOctree& old_tree, const ConvexCutCellMesh& old_mesh,
    const LinearOctree& new_tree, const ConvexCutCellMesh& new_mesh,
    const std::optional<AABB>& affected_bounds,
    double geometric_tolerance = 0.0);

} // namespace cartmesh
