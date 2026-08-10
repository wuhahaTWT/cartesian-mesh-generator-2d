#pragma once

#include <cstdint>
#include <stdexcept>

namespace cartmesh {

using MortonCode = std::uint64_t;
inline constexpr std::uint32_t max_morton_coordinate = (1U << 21U) - 1U;

[[nodiscard]] inline MortonCode encode_morton_3d(std::uint32_t x, std::uint32_t y,
                                                 std::uint32_t z) {
    if (x > max_morton_coordinate || y > max_morton_coordinate || z > max_morton_coordinate) {
        throw std::out_of_range("Morton 坐标的每个方向最多使用 21 位");
    }
    MortonCode code = 0;
    for (std::uint32_t bit = 0; bit < 21U; ++bit) {
        code |= (static_cast<MortonCode>((x >> bit) & 1U) << (3U * bit));
        code |= (static_cast<MortonCode>((y >> bit) & 1U) << (3U * bit + 1U));
        code |= (static_cast<MortonCode>((z >> bit) & 1U) << (3U * bit + 2U));
    }
    return code;
}

struct MortonCoordinates {
    std::uint32_t x{};
    std::uint32_t y{};
    std::uint32_t z{};
};

[[nodiscard]] inline MortonCoordinates decode_morton_3d(MortonCode code) noexcept {
    MortonCoordinates result;
    for (std::uint32_t bit = 0; bit < 21U; ++bit) {
        result.x |= static_cast<std::uint32_t>((code >> (3U * bit)) & 1U) << bit;
        result.y |= static_cast<std::uint32_t>((code >> (3U * bit + 1U)) & 1U) << bit;
        result.z |= static_cast<std::uint32_t>((code >> (3U * bit + 2U)) & 1U) << bit;
    }
    return result;
}

} // 命名空间 cartmesh
