#include "cartmesh/incremental/IncrementalRemesher.hpp"

#include <algorithm>
#include <set>
#include <vector>

namespace cartmesh {
namespace {

[[nodiscard]] bool intersects_affected(const AABB& bounds,
                                       const std::optional<AABB>& affected) {
    return affected && affected->intersects(bounds);
}

[[nodiscard]] bool parent_order(OctreeNodeCode first,
                                OctreeNodeCode second) {
    const auto first_node = decode_octree_node(first);
    const auto second_node = decode_octree_node(second);
    if (first_node.level != second_node.level) {
        return first_node.level > second_node.level;
    }
    return first < second;
}

} // namespace

IncrementalOctreeResult update_octree_incrementally(
    LinearOctree old_tree, const OctreeRefinementEngine& new_engine,
    const std::optional<AABB>& affected_bounds) {
    const LinearOctree old_snapshot = old_tree;
    IncrementalOctreeStatistics statistics;
    statistics.old_leaf_count = old_tree.leaf_count();

    while (true) {
        std::vector<OctreeNodeCode> parents;
        parents.reserve(static_cast<std::size_t>(old_tree.leaf_count() / 8U));
        for (const auto leaf : old_tree.leaf_codes()) {
            const auto node = decode_octree_node(leaf);
            if (node.level <= old_tree.base_level()) continue;
            parents.push_back(octree_parent(leaf));
        }
        std::sort(parents.begin(), parents.end(), parent_order);
        parents.erase(std::unique(parents.begin(), parents.end()), parents.end());

        std::uint64_t coarsened_this_iteration = 0;
        for (const auto parent : parents) {
            const auto node = decode_octree_node(parent);
            const AABB bounds = old_tree.cell_bounds(parent);
            const auto desired = new_engine.desired_level(old_tree, parent, bounds);
            if (desired <= node.level && old_tree.coarsen_parent(parent)) {
                ++coarsened_this_iteration;
            }
        }
        statistics.coarsened_parent_count += coarsened_this_iteration;
        if (coarsened_this_iteration == 0) break;
    }

    const auto adaptation = new_engine.apply(old_tree);
    statistics.rule_split_count = adaptation.rule_refinement.split_count;
    statistics.balance_split_count = adaptation.balance.split_count;
    statistics.new_leaf_count = old_tree.leaf_count();

    std::set<OctreeNodeCode> old_codes(old_snapshot.leaf_codes().begin(),
                                       old_snapshot.leaf_codes().end());
    std::set<OctreeNodeCode> new_codes(old_tree.leaf_codes().begin(),
                                       old_tree.leaf_codes().end());
    for (const auto code : old_codes) {
        if (new_codes.contains(code)) {
            ++statistics.preserved_leaf_count;
            continue;
        }
        ++statistics.removed_leaf_count;
        if (intersects_affected(old_snapshot.cell_bounds(code), affected_bounds)) {
            ++statistics.local_removed_leaf_count;
        } else {
            ++statistics.balance_closure_removed_leaf_count;
        }
    }
    for (const auto code : new_codes) {
        if (old_codes.contains(code)) continue;
        ++statistics.created_leaf_count;
        if (intersects_affected(old_tree.cell_bounds(code), affected_bounds)) {
            ++statistics.local_created_leaf_count;
        } else {
            ++statistics.balance_closure_created_leaf_count;
        }
    }
    return {std::move(old_tree), adaptation, statistics};
}

} // namespace cartmesh
