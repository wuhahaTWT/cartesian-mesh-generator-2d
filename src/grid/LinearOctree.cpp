#include "cartmesh/grid/LinearOctree.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

[[nodiscard]] constexpr std::array<FaceDirection, 6> all_faces() noexcept {
    return {FaceDirection::negative_x, FaceDirection::positive_x,
            FaceDirection::negative_y, FaceDirection::positive_y,
            FaceDirection::negative_z, FaceDirection::positive_z};
}

[[nodiscard]] std::uint8_t absolute_level_difference(std::uint8_t first,
                                                     std::uint8_t second) noexcept {
    return first >= second ? static_cast<std::uint8_t>(first - second)
                           : static_cast<std::uint8_t>(second - first);
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (std::uint32_t shift = 0; shift < 64U; shift += 8U) {
        hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
}

} // 匿名命名空间

LinearOctree::LinearOctree(AABB domain, std::uint8_t base_level,
                           std::uint8_t maximum_level)
    : domain_(domain), base_level_(base_level), maximum_level_(maximum_level) {
    if (!domain_.has_positive_volume()) {
        throw std::invalid_argument("八叉树计算域必须具有正体积");
    }
    if (maximum_level_ > maximum_octree_level || base_level_ > maximum_level_) {
        throw std::invalid_argument("八叉树必须满足 0 <= baseLevel <= maximumLevel <= 21");
    }
    const std::uint32_t count_per_axis = 1U << base_level_;
    const std::uint64_t initial_count =
        std::uint64_t{1} << (3U * static_cast<std::uint32_t>(base_level_));
    if (initial_count > leaves_.max_size()) {
        throw std::length_error("基础八叉树叶数超出容器范围");
    }
    leaves_.reserve(static_cast<std::size_t>(initial_count));
    for (std::uint32_t z = 0; z < count_per_axis; ++z) {
        for (std::uint32_t y = 0; y < count_per_axis; ++y) {
            for (std::uint32_t x = 0; x < count_per_axis; ++x) {
                leaves_.push_back(encode_octree_node(base_level_, x, y, z));
            }
        }
    }
    std::sort(leaves_.begin(), leaves_.end(), [&](OctreeNodeCode lhs, OctreeNodeCode rhs) {
        return anchor(lhs) < anchor(rhs);
    });
}

MortonCode LinearOctree::anchor(OctreeNodeCode code) const {
    return octree_anchor_morton(code, maximum_level_);
}

std::uint64_t LinearOctree::span(OctreeNodeCode code) const {
    return octree_morton_span(code, maximum_level_);
}

bool LinearOctree::code_is_within_levels(OctreeNodeCode code) const {
    try {
        const auto node = decode_octree_node(code);
        return node.level >= base_level_ && node.level <= maximum_level_;
    } catch (const std::exception&) {
        return false;
    }
}

OctreeNodeCode LinearOctree::leaf_code(std::uint64_t leaf_id) const {
    if (leaf_id >= leaves_.size()) {
        throw std::out_of_range("八叉树叶 ID 超出范围");
    }
    return leaves_[static_cast<std::size_t>(leaf_id)];
}

AABB LinearOctree::cell_bounds(OctreeNodeCode code) const {
    if (!code_is_within_levels(code)) {
        throw std::out_of_range("八叉树节点码不属于当前层级范围");
    }
    const auto node = decode_octree_node(code);
    const double count_per_axis = std::ldexp(1.0, node.level);
    const auto extent = domain_.extent();
    const Vec3 spacing{extent.x / count_per_axis, extent.y / count_per_axis,
                       extent.z / count_per_axis};
    const auto root_minimum = domain_.minimum();
    const Vec3 minimum{root_minimum.x + spacing.x * static_cast<double>(node.x),
                       root_minimum.y + spacing.y * static_cast<double>(node.y),
                       root_minimum.z + spacing.z * static_cast<double>(node.z)};
    const Vec3 maximum{
        root_minimum.x + spacing.x * static_cast<double>(node.x + 1U),
        root_minimum.y + spacing.y * static_cast<double>(node.y + 1U),
        root_minimum.z + spacing.z * static_cast<double>(node.z + 1U)};
    return AABB(minimum, maximum);
}

Vec3 LinearOctree::cell_center(OctreeNodeCode code) const {
    return cell_bounds(code).center();
}

double LinearOctree::cell_volume(OctreeNodeCode code) const {
    return cell_bounds(code).volume();
}

std::optional<std::uint64_t> LinearOctree::find_leaf(OctreeNodeCode code) const {
    if (!code_is_within_levels(code)) {
        return std::nullopt;
    }
    const auto target_anchor = anchor(code);
    const auto iterator = std::lower_bound(
        leaves_.begin(), leaves_.end(), target_anchor,
        [&](OctreeNodeCode candidate, MortonCode value) { return anchor(candidate) < value; });
    if (iterator == leaves_.end() || *iterator != code) {
        return std::nullopt;
    }
    return static_cast<std::uint64_t>(std::distance(leaves_.begin(), iterator));
}

std::optional<std::uint64_t> LinearOctree::find_leaf_covering_maximum_level_cell(
    std::uint32_t x, std::uint32_t y, std::uint32_t z) const {
    const std::uint32_t coordinate_limit = 1U << maximum_level_;
    if (x >= coordinate_limit || y >= coordinate_limit || z >= coordinate_limit) {
        return std::nullopt;
    }
    const MortonCode target = encode_morton_3d(x, y, z);
    const auto upper = std::upper_bound(
        leaves_.begin(), leaves_.end(), target,
        [&](MortonCode value, OctreeNodeCode candidate) { return value < anchor(candidate); });
    if (upper == leaves_.begin()) {
        return std::nullopt;
    }
    const auto candidate = std::prev(upper);
    const MortonCode candidate_anchor = anchor(*candidate);
    if (target - candidate_anchor >= span(*candidate)) {
        return std::nullopt;
    }
    return static_cast<std::uint64_t>(std::distance(leaves_.begin(), candidate));
}

std::vector<OctreeNodeCode> LinearOctree::face_neighbors(OctreeNodeCode code,
                                                         FaceDirection direction) const {
    if (!find_leaf(code)) {
        throw std::out_of_range("只能查询当前叶单元的面邻居");
    }
    const auto node = decode_octree_node(code);
    const auto level_shift = static_cast<std::uint8_t>(maximum_level_ - node.level);
    const std::uint64_t leaf_size = std::uint64_t{1} << level_shift;
    const std::array<std::uint64_t, 3> origin = {
        static_cast<std::uint64_t>(node.x) << level_shift,
        static_cast<std::uint64_t>(node.y) << level_shift,
        static_cast<std::uint64_t>(node.z) << level_shift};
    const std::uint64_t coordinate_limit = std::uint64_t{1} << maximum_level_;
    const auto raw_direction = static_cast<std::uint8_t>(direction);
    const std::size_t normal_axis = raw_direction / 2U;
    const bool positive = (raw_direction & 1U) != 0U;
    std::int64_t normal_coordinate = 0;
    if (positive) {
        normal_coordinate = static_cast<std::int64_t>(origin[normal_axis] + leaf_size);
        if (normal_coordinate >= static_cast<std::int64_t>(coordinate_limit)) {
            return {};
        }
    } else {
        if (origin[normal_axis] == 0) {
            return {};
        }
        normal_coordinate = static_cast<std::int64_t>(origin[normal_axis] - 1U);
    }
    const std::array<std::array<std::size_t, 2>, 3> tangential_axes = {
        std::array<std::size_t, 2>{1, 2}, std::array<std::size_t, 2>{0, 2},
        std::array<std::size_t, 2>{0, 1}};
    const auto u_axis = tangential_axes[normal_axis][0];
    const auto v_axis = tangential_axes[normal_axis][1];
    const std::uint64_t u_end = origin[u_axis] + leaf_size;
    const std::uint64_t v_end = origin[v_axis] + leaf_size;
    std::vector<OctreeNodeCode> result;
    std::uint64_t v = origin[v_axis];
    while (v < v_end) {
        std::uint64_t next_v = v_end;
        std::uint64_t u = origin[u_axis];
        while (u < u_end) {
            std::array<std::uint32_t, 3> probe{};
            probe[normal_axis] = static_cast<std::uint32_t>(normal_coordinate);
            probe[u_axis] = static_cast<std::uint32_t>(u);
            probe[v_axis] = static_cast<std::uint32_t>(v);
            const auto neighbor_id = find_leaf_covering_maximum_level_cell(
                probe[0], probe[1], probe[2]);
            if (!neighbor_id) {
                throw std::logic_error("八叉树分区存在面邻居查询空洞");
            }
            const auto neighbor_code = leaf_code(*neighbor_id);
            result.push_back(neighbor_code);
            const auto neighbor = decode_octree_node(neighbor_code);
            const auto neighbor_shift =
                static_cast<std::uint8_t>(maximum_level_ - neighbor.level);
            const std::uint64_t neighbor_size = std::uint64_t{1} << neighbor_shift;
            const std::array<std::uint64_t, 3> neighbor_origin = {
                static_cast<std::uint64_t>(neighbor.x) << neighbor_shift,
                static_cast<std::uint64_t>(neighbor.y) << neighbor_shift,
                static_cast<std::uint64_t>(neighbor.z) << neighbor_shift};
            const std::uint64_t neighbor_u_end = neighbor_origin[u_axis] + neighbor_size;
            const std::uint64_t neighbor_v_end = neighbor_origin[v_axis] + neighbor_size;
            if (neighbor_u_end <= u || neighbor_v_end <= v) {
                throw std::logic_error("八叉树面邻居未覆盖查询子区域");
            }
            u = std::min(u_end, neighbor_u_end);
            next_v = std::min(next_v, std::min(v_end, neighbor_v_end));
        }
        if (next_v <= v) {
            throw std::logic_error("八叉树面邻居平铺没有前进");
        }
        v = next_v;
    }
    std::sort(result.begin(), result.end(), [&](OctreeNodeCode lhs, OctreeNodeCode rhs) {
        return anchor(lhs) < anchor(rhs);
    });
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool LinearOctree::refine_leaf(OctreeNodeCode code) {
    const auto leaf_id = find_leaf(code);
    if (!leaf_id) {
        return false;
    }
    const auto node = decode_octree_node(code);
    if (node.level >= maximum_level_) {
        return false;
    }
    const auto position = leaves_.begin() + static_cast<std::ptrdiff_t>(*leaf_id);
    std::array<OctreeNodeCode, 8> children{};
    for (std::uint8_t child_index = 0; child_index < 8U; ++child_index) {
        children[child_index] = octree_child(code, child_index);
    }
    const auto offset = static_cast<std::size_t>(*leaf_id);
    leaves_.erase(position);
    leaves_.insert(leaves_.begin() + static_cast<std::ptrdiff_t>(offset), children.begin(),
                   children.end());
    return true;
}

bool LinearOctree::coarsen_parent(OctreeNodeCode parent_code) {
    const auto parent = decode_octree_node(parent_code);
    if (parent.level < base_level_ || parent.level >= maximum_level_) {
        return false;
    }
    std::array<std::uint64_t, 8> child_ids{};
    for (std::uint8_t child_index = 0; child_index < 8U; ++child_index) {
        const auto child_id = find_leaf(octree_child(parent_code, child_index));
        if (!child_id) {
            return false;
        }
        child_ids[child_index] = *child_id;
    }
    for (std::uint8_t child_index = 1; child_index < 8U; ++child_index) {
        if (child_ids[child_index] != child_ids[0] + child_index) {
            throw std::logic_error("同一父节点的八个叶子节点未保持 Morton 连续顺序");
        }
    }
    const auto first = leaves_.begin() + static_cast<std::ptrdiff_t>(child_ids[0]);
    leaves_.erase(first, first + 8);
    leaves_.insert(leaves_.begin() + static_cast<std::ptrdiff_t>(child_ids[0]), parent_code);
    return true;
}

void LinearOctree::replace_marked_with_children(
    const std::vector<OctreeNodeCode>& marked) {
    if (marked.empty()) {
        return;
    }
    std::vector<OctreeNodeCode> output;
    output.reserve(leaves_.size() + marked.size() * 7U);
    std::size_t marked_index = 0;
    for (const auto leaf : leaves_) {
        while (marked_index < marked.size() && anchor(marked[marked_index]) < anchor(leaf)) {
            ++marked_index;
        }
        if (marked_index < marked.size() && marked[marked_index] == leaf) {
            for (std::uint8_t child_index = 0; child_index < 8U; ++child_index) {
                output.push_back(octree_child(leaf, child_index));
            }
        } else {
            output.push_back(leaf);
        }
    }
    leaves_.swap(output);
}

OctreeRefinementStatistics LinearOctree::refine_to_desired_levels(
    const DesiredLevelFunction& desired_level) {
    if (!desired_level) {
        throw std::invalid_argument("八叉树目标层级函数不得为空");
    }
    OctreeRefinementStatistics statistics;
    statistics.initial_leaf_count = leaves_.size();
    while (true) {
        std::vector<OctreeNodeCode> marked;
        for (const auto leaf : leaves_) {
            const auto node = decode_octree_node(leaf);
            const auto target = std::min(desired_level(leaf, cell_bounds(leaf)), maximum_level_);
            if (target > node.level) {
                marked.push_back(leaf);
            }
        }
        if (marked.empty()) {
            break;
        }
        replace_marked_with_children(marked);
        statistics.split_count += marked.size();
        ++statistics.iteration_count;
    }
    statistics.final_leaf_count = leaves_.size();
    return statistics;
}

OctreeBalanceReport LinearOctree::check_face_balance() const {
    const auto levels = level_statistics();
    const auto global_difference = absolute_level_difference(
        levels.minimum_leaf_level, levels.maximum_leaf_level);
    if (global_difference <= 1U) {
        return {true, 0, global_difference};
    }
    OctreeBalanceReport report{true, 0, 0};
    for (const auto leaf : leaves_) {
        const auto leaf_level = decode_octree_node(leaf).level;
        for (const auto direction : all_faces()) {
            for (const auto neighbor : face_neighbors(leaf, direction)) {
                if (anchor(neighbor) < anchor(leaf)) {
                    continue;
                }
                const auto difference = absolute_level_difference(
                    leaf_level, decode_octree_node(neighbor).level);
                report.maximum_level_difference =
                    std::max(report.maximum_level_difference, difference);
                if (difference > 1U) {
                    report.balanced = false;
                    ++report.violating_face_pair_count;
                }
            }
        }
    }
    return report;
}

OctreeBalanceStatistics LinearOctree::balance_faces_2_to_1() {
    OctreeBalanceStatistics statistics;
    statistics.initial_leaf_count = leaves_.size();
    const auto levels = level_statistics();
    if (absolute_level_difference(levels.minimum_leaf_level,
                                  levels.maximum_leaf_level) <= 1U) {
        statistics.final_leaf_count = leaves_.size();
        return statistics;
    }
    while (true) {
        std::vector<OctreeNodeCode> marked;
        for (const auto leaf : leaves_) {
            const auto leaf_level = decode_octree_node(leaf).level;
            for (const auto direction : all_faces()) {
                for (const auto neighbor : face_neighbors(leaf, direction)) {
                    const auto neighbor_level = decode_octree_node(neighbor).level;
                    if (leaf_level + 1U < neighbor_level) {
                        marked.push_back(leaf);
                    }
                }
            }
        }
        if (marked.empty()) {
            break;
        }
        std::sort(marked.begin(), marked.end(), [&](OctreeNodeCode lhs, OctreeNodeCode rhs) {
            return anchor(lhs) < anchor(rhs);
        });
        marked.erase(std::unique(marked.begin(), marked.end()), marked.end());
        replace_marked_with_children(marked);
        statistics.split_count += marked.size();
        ++statistics.iteration_count;
    }
    statistics.final_leaf_count = leaves_.size();
    return statistics;
}

OctreeLevelStatistics LinearOctree::level_statistics() const {
    OctreeLevelStatistics statistics;
    statistics.leaf_count_by_level.resize(static_cast<std::size_t>(maximum_level_) + 1U, 0);
    statistics.minimum_leaf_level = maximum_level_;
    statistics.maximum_leaf_level = 0;
    for (const auto leaf : leaves_) {
        const auto level = decode_octree_node(leaf).level;
        ++statistics.leaf_count_by_level[level];
        statistics.minimum_leaf_level = std::min(statistics.minimum_leaf_level, level);
        statistics.maximum_leaf_level = std::max(statistics.maximum_leaf_level, level);
    }
    return statistics;
}

bool LinearOctree::validate_partition() const {
    if (leaves_.empty()) {
        return false;
    }
    MortonCode expected_anchor = 0;
    for (const auto leaf : leaves_) {
        if (!code_is_within_levels(leaf) || anchor(leaf) != expected_anchor) {
            return false;
        }
        expected_anchor += span(leaf);
    }
    const auto expected_total = std::uint64_t{1} <<
                                (3U * static_cast<std::uint32_t>(maximum_level_));
    return expected_anchor == expected_total;
}

std::uint64_t LinearOctree::result_hash_fnv1a64() const noexcept {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_byte(hash, base_level_);
    hash_byte(hash, maximum_level_);
    const std::array<double, 6> bounds = {
        domain_.minimum().x, domain_.minimum().y, domain_.minimum().z,
        domain_.maximum().x, domain_.maximum().y, domain_.maximum().z};
    for (const double value : bounds) {
        hash_u64(hash, std::bit_cast<std::uint64_t>(value));
    }
    for (const auto leaf : leaves_) {
        hash_u64(hash, leaf);
    }
    return hash;
}

} // 命名空间 cartmesh
