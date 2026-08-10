#pragma once

#include "cartmesh/grid/MortonCode.hpp"

#include <cstdint>

namespace cartmesh {

struct CellKey {
    std::uint8_t level{};
    std::uint32_t i{};
    std::uint32_t j{};
    std::uint32_t k{};

    [[nodiscard]] MortonCode morton_code() const { return encode_morton_3d(i, j, k); }
    [[nodiscard]] constexpr bool operator==(const CellKey&) const noexcept = default;
};

} // 命名空间 cartmesh
