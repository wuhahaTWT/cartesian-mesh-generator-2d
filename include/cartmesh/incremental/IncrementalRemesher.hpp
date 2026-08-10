#pragma once

#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"

#include <cstdint>
#include <optional>

namespace cartmesh {

struct IncrementalOctreeStatistics {
    std::uint64_t old_leaf_count{};
    std::uint64_t new_leaf_count{};
    std::uint64_t coarsened_parent_count{};
    std::uint64_t rule_split_count{};
    std::uint64_t balance_split_count{};
    std::uint64_t preserved_leaf_count{};
    std::uint64_t created_leaf_count{};
    std::uint64_t removed_leaf_count{};
    std::uint64_t local_created_leaf_count{};
    std::uint64_t local_removed_leaf_count{};
    std::uint64_t balance_closure_created_leaf_count{};
    std::uint64_t balance_closure_removed_leaf_count{};
};

struct IncrementalOctreeResult {
    LinearOctree tree;
    OctreeAdaptationStatistics adaptation;
    IncrementalOctreeStatistics statistics;
};

// 从旧树出发，先消除新目标不再需要的细化，再应用新几何规则和 2:1 平衡。
// affected_bounds 只用于审计实际变化属于几何局部还是平衡闭包；算法仍全局扫描
// 紧凑叶数组以保证不会遗漏必须撤销的旧平衡细化。
[[nodiscard]] IncrementalOctreeResult update_octree_incrementally(
    LinearOctree old_tree, const OctreeRefinementEngine& new_engine,
    const std::optional<AABB>& affected_bounds);

} // namespace cartmesh
