#pragma once

#include "cartmesh/grid/MortonCode.hpp"

#include <bit>
#include <cstdint>
#include <stdexcept>

namespace cartmesh {

using OctreeNodeCode = std::uint64_t;
inline constexpr std::uint8_t maximum_octree_level = 21;

struct OctreeCoordinates {
    std::uint8_t level{};
    std::uint32_t x{};
    std::uint32_t y{};
    std::uint32_t z{};
};

[[nodiscard]] inline OctreeNodeCode encode_octree_node(std::uint8_t level, std::uint32_t x,
                                                       std::uint32_t y, std::uint32_t z) {
    if (level > maximum_octree_level) {
        throw std::out_of_range("八叉树层级最多为 21");
    }
    const std::uint32_t coordinate_limit = 1U << level;
    if (x >= coordinate_limit || y >= coordinate_limit || z >= coordinate_limit) {
        throw std::out_of_range("八叉树层级坐标超出当前层级范围");
    }
    const auto sentinel_shift = static_cast<std::uint32_t>(3U * level);
    return (OctreeNodeCode{1} << sentinel_shift) | encode_morton_3d(x, y, z);
}

[[nodiscard]] inline OctreeCoordinates decode_octree_node(OctreeNodeCode code) {
    if (code == 0) {
        throw std::invalid_argument("八叉树节点码不得为零");
    }
    const auto highest_bit = static_cast<std::uint32_t>(std::bit_width(code)) - 1U;
    if ((highest_bit % 3U) != 0U) {
        throw std::invalid_argument("八叉树节点码缺少合法层级哨兵位");
    }
    const auto level = static_cast<std::uint8_t>(highest_bit / 3U);
    if (level > maximum_octree_level) {
        throw std::invalid_argument("八叉树节点码层级超出 21");
    }
    const auto sentinel = OctreeNodeCode{1} << highest_bit;
    const auto coordinates = decode_morton_3d(code ^ sentinel);
    const std::uint32_t coordinate_limit = 1U << level;
    if (coordinates.x >= coordinate_limit || coordinates.y >= coordinate_limit ||
        coordinates.z >= coordinate_limit) {
        throw std::invalid_argument("八叉树节点码包含当前层级之外的坐标位");
    }
    return {level, coordinates.x, coordinates.y, coordinates.z};
}

[[nodiscard]] inline OctreeNodeCode octree_parent(OctreeNodeCode code) {
    const auto node = decode_octree_node(code);
    if (node.level == 0) {
        throw std::out_of_range("根节点没有父节点");
    }
    return encode_octree_node(static_cast<std::uint8_t>(node.level - 1U), node.x >> 1U,
                              node.y >> 1U, node.z >> 1U);
}

[[nodiscard]] inline OctreeNodeCode octree_child(OctreeNodeCode code,
                                                 std::uint8_t child_index) {
    const auto node = decode_octree_node(code);
    if (node.level >= maximum_octree_level) {
        throw std::out_of_range("21 层节点不能继续细分");
    }
    if (child_index >= 8U) {
        throw std::out_of_range("八叉树子节点索引必须位于 0..7");
    }
    return encode_octree_node(
        static_cast<std::uint8_t>(node.level + 1U), node.x * 2U + (child_index & 1U),
        node.y * 2U + ((child_index >> 1U) & 1U),
        node.z * 2U + ((child_index >> 2U) & 1U));
}

[[nodiscard]] inline MortonCode octree_anchor_morton(OctreeNodeCode code,
                                                      std::uint8_t maximum_level) {
    const auto node = decode_octree_node(code);
    if (maximum_level < node.level || maximum_level > maximum_octree_level) {
        throw std::out_of_range("Morton 锚点最大层级无效");
    }
    const auto shift = static_cast<std::uint8_t>(maximum_level - node.level);
    return encode_morton_3d(node.x << shift, node.y << shift, node.z << shift);
}

[[nodiscard]] inline std::uint64_t octree_morton_span(OctreeNodeCode code,
                                                      std::uint8_t maximum_level) {
    const auto node = decode_octree_node(code);
    if (maximum_level < node.level || maximum_level > maximum_octree_level) {
        throw std::out_of_range("Morton 区间最大层级无效");
    }
    return std::uint64_t{1} <<
           (3U * static_cast<std::uint32_t>(maximum_level - node.level));
}

} // 命名空间 cartmesh
