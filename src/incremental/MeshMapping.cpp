#include "cartmesh/incremental/MeshMapping.hpp"

#include "cartmesh/cutcell/ConvexPolyhedron.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

[[nodiscard]] std::vector<const FluidCellGeometry*> fluid_cells_by_leaf(
    const LinearOctree& tree, const ConvexCutCellMesh& mesh) {
    std::vector<const FluidCellGeometry*> result(
        static_cast<std::size_t>(tree.leaf_count()), nullptr);
    for (const auto& cell : mesh.fluid_cells) {
        if (cell.background_cell_id >= tree.leaf_count()) {
            throw std::invalid_argument(
                "增量映射输入的 background ID 超出八叉树叶范围");
        }
        result[static_cast<std::size_t>(cell.background_cell_id)] = &cell;
    }
    return result;
}

[[nodiscard]] AABB overlapping_bounds(const AABB& first,
                                      const AABB& second) {
    const Vec3 minimum{std::max(first.minimum().x, second.minimum().x),
                       std::max(first.minimum().y, second.minimum().y),
                       std::max(first.minimum().z, second.minimum().z)};
    const Vec3 maximum{std::min(first.maximum().x, second.maximum().x),
                       std::min(first.maximum().y, second.maximum().y),
                       std::min(first.maximum().z, second.maximum().z)};
    return AABB(minimum, maximum);
}

[[nodiscard]] std::vector<ConvexPolyhedron> fluid_pieces(
    const FluidCellGeometry* cell, const AABB& bounds) {
    if (!cell) return {};
    if (!cell->cut) return {make_box_polyhedron(bounds)};
    std::vector<ConvexPolyhedron> result;
    result.reserve(cell->fluid_polyhedron_pieces.size());
    for (const auto& piece : cell->fluid_polyhedron_pieces) {
        result.push_back(piece.polyhedron);
    }
    return result;
}

[[nodiscard]] double convex_overlap_volume(
    const ConvexPolyhedron& first, const ConvexPolyhedron& second,
    double length_tolerance) {
    ConvexPolyhedron overlap = first;
    const auto second_geometry = measure_polyhedron(second);
    if (!second_geometry.closed || !second_geometry.positive_volume ||
        second_geometry.faces.size() != second.faces.size()) {
        throw std::runtime_error(
            "增量流体映射遇到无效的新网格凸多面体片");
    }
    for (std::size_t face = 0; face < second.faces.size(); ++face) {
        const auto& record = second.faces[face];
        if (record.vertex_indices.empty()) continue;
        const Vec3 point = second.vertices[record.vertex_indices.front()];
        overlap = clip_convex_polyhedron(
            overlap,
            OrientedHalfSpace(point,
                              second_geometry.faces[face].outward_normal),
            length_tolerance, PolyhedronFaceKind::internal_partition);
        if (overlap.empty()) return 0.0;
    }
    const auto geometry = measure_polyhedron(overlap);
    if (!geometry.closed || !geometry.positive_volume) {
        if (geometry.volume <= length_tolerance * length_tolerance *
                                   length_tolerance) {
            return 0.0;
        }
        throw std::runtime_error(
            "增量流体映射的凸片交集不是闭合正体积");
    }
    return geometry.volume;
}

[[nodiscard]] double exact_fluid_overlap(
    const FluidCellGeometry* old_cell, const AABB& old_bounds,
    const FluidCellGeometry* new_cell, const AABB& new_bounds,
    double length_tolerance) {
    if (!old_cell || !new_cell) return 0.0;
    if (!old_cell->cut && !new_cell->cut) {
        return overlapping_bounds(old_bounds, new_bounds).volume();
    }
    const auto old_pieces = fluid_pieces(old_cell, old_bounds);
    const auto new_pieces = fluid_pieces(new_cell, new_bounds);
    double result = 0.0;
    for (const auto& old_piece : old_pieces) {
        for (const auto& new_piece : new_pieces) {
            result += convex_overlap_volume(old_piece, new_piece,
                                            length_tolerance);
        }
    }
    return result;
}

[[nodiscard]] IncrementalTopologyRelation relation_for_pair(
    OctreeNodeCode old_code, OctreeNodeCode new_code, const AABB& bounds,
    const std::optional<AABB>& affected_bounds) {
    if (old_code == new_code) {
        return affected_bounds && affected_bounds->intersects(bounds)
                   ? IncrementalTopologyRelation::rebuilt
                   : IncrementalTopologyRelation::preserved;
    }
    const auto old_level = decode_octree_node(old_code).level;
    const auto new_level = decode_octree_node(new_code).level;
    return old_level < new_level ? IncrementalTopologyRelation::refined
                                 : IncrementalTopologyRelation::coarsened;
}

} // namespace

IncrementalMeshMapping build_incremental_mesh_mapping(
    const LinearOctree& old_tree, const ConvexCutCellMesh& old_mesh,
    const LinearOctree& new_tree, const ConvexCutCellMesh& new_mesh,
    const std::optional<AABB>& affected_bounds,
    double geometric_tolerance) {
    if (!std::isfinite(geometric_tolerance) || geometric_tolerance < 0.0) {
        throw std::invalid_argument(
            "增量映射几何容差必须是非负有限数");
    }
    if (old_tree.domain().minimum().x != new_tree.domain().minimum().x ||
        old_tree.domain().minimum().y != new_tree.domain().minimum().y ||
        old_tree.domain().minimum().z != new_tree.domain().minimum().z ||
        old_tree.domain().maximum().x != new_tree.domain().maximum().x ||
        old_tree.domain().maximum().y != new_tree.domain().maximum().y ||
        old_tree.domain().maximum().z != new_tree.domain().maximum().z ||
        old_tree.maximum_level() != new_tree.maximum_level()) {
        throw std::invalid_argument(
            "增量映射要求旧、新八叉树使用相同固定域和最大层级");
    }
    const double length_tolerance = std::max(
        geometric_tolerance,
        512.0 * std::numeric_limits<double>::epsilon() *
            norm(old_tree.domain().extent()));
    const auto old_cells = fluid_cells_by_leaf(old_tree, old_mesh);
    const auto new_cells = fluid_cells_by_leaf(new_tree, new_mesh);

    IncrementalMeshMapping result;
    result.old_fluid_volume = old_mesh.total_fluid_volume;
    result.new_fluid_volume = new_mesh.total_fluid_volume;
    std::uint64_t old_id = 0;
    std::uint64_t new_id = 0;
    while (old_id < old_tree.leaf_count() && new_id < new_tree.leaf_count()) {
        const auto old_code = old_tree.leaf_code(old_id);
        const auto new_code = new_tree.leaf_code(new_id);
        const auto old_anchor =
            octree_anchor_morton(old_code, old_tree.maximum_level());
        const auto new_anchor =
            octree_anchor_morton(new_code, new_tree.maximum_level());
        const auto old_end = old_anchor +
            octree_morton_span(old_code, old_tree.maximum_level());
        const auto new_end = new_anchor +
            octree_morton_span(new_code, new_tree.maximum_level());
        if (old_end <= new_anchor) {
            ++old_id;
            continue;
        }
        if (new_end <= old_anchor) {
            ++new_id;
            continue;
        }

        const AABB old_bounds = old_tree.cell_bounds(old_code);
        const AABB new_bounds = new_tree.cell_bounds(new_code);
        const AABB overlap_bounds = overlapping_bounds(old_bounds, new_bounds);
        const auto* old_cell = old_cells[static_cast<std::size_t>(old_id)];
        const auto* new_cell = new_cells[static_cast<std::size_t>(new_id)];
        IncrementalCellMappingEntry entry;
        entry.old_code = old_code;
        entry.new_code = new_code;
        entry.relation = relation_for_pair(
            old_code, new_code, overlap_bounds, affected_bounds);
        entry.background_overlap_volume = overlap_bounds.volume();
        entry.old_fluid_volume = old_cell ? old_cell->volume : 0.0;
        entry.new_fluid_volume = new_cell ? new_cell->volume : 0.0;
        entry.old_cell_is_fluid = old_cell != nullptr;
        entry.new_cell_is_fluid = new_cell != nullptr;
        if (entry.relation == IncrementalTopologyRelation::preserved) {
            if (entry.old_cell_is_fluid != entry.new_cell_is_fluid ||
                (old_cell &&
                 std::abs(old_cell->volume - new_cell->volume) >
                     length_tolerance * length_tolerance * length_tolerance)) {
                throw std::runtime_error(
                    "影响范围外的 preserved 单元几何发生变化");
            }
            entry.shared_fluid_overlap_volume =
                old_cell ? old_cell->volume : 0.0;
        } else {
            entry.shared_fluid_overlap_volume = exact_fluid_overlap(
                old_cell, old_bounds, new_cell, new_bounds,
                length_tolerance);
        }
        entry.exact_fluid_overlap = true;
        if (entry.old_fluid_volume > 0.0) {
            entry.old_volume_preserved_fraction =
                entry.shared_fluid_overlap_volume / entry.old_fluid_volume;
        }
        if (entry.new_fluid_volume > 0.0) {
            entry.new_volume_from_old_fraction =
                entry.shared_fluid_overlap_volume / entry.new_fluid_volume;
        }
        result.total_background_overlap_volume +=
            entry.background_overlap_volume;
        result.total_shared_fluid_overlap_volume +=
            entry.shared_fluid_overlap_volume;
        switch (entry.relation) {
        case IncrementalTopologyRelation::preserved:
            ++result.preserved_pair_count;
            break;
        case IncrementalTopologyRelation::rebuilt:
            ++result.rebuilt_pair_count;
            break;
        case IncrementalTopologyRelation::refined:
            ++result.refined_pair_count;
            break;
        case IncrementalTopologyRelation::coarsened:
            ++result.coarsened_pair_count;
            break;
        }
        result.entries.push_back(entry);

        if (old_end <= new_end) ++old_id;
        if (new_end <= old_end) ++new_id;
    }
    const double volume_tolerance =
        length_tolerance * length_tolerance * length_tolerance;
    const double removed =
        result.old_fluid_volume - result.total_shared_fluid_overlap_volume;
    const double created =
        result.new_fluid_volume - result.total_shared_fluid_overlap_volume;
    result.removed_fluid_volume =
        std::abs(removed) <= volume_tolerance ? 0.0 : std::max(0.0, removed);
    result.created_fluid_volume =
        std::abs(created) <= volume_tolerance ? 0.0 : std::max(0.0, created);
    return result;
}

} // namespace cartmesh
