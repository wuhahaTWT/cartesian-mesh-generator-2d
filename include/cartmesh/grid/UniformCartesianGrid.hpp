#pragma once

#include "cartmesh/geometry/AABB.hpp"
#include "cartmesh/grid/CellKey.hpp"

#include <cstdint>

namespace cartmesh {

class UniformCartesianGrid {
  public:
    UniformCartesianGrid(AABB domain, std::uint32_t nx, std::uint32_t ny, std::uint32_t nz);

    [[nodiscard]] constexpr const AABB& domain() const noexcept { return domain_; }
    [[nodiscard]] constexpr std::uint32_t nx() const noexcept { return nx_; }
    [[nodiscard]] constexpr std::uint32_t ny() const noexcept { return ny_; }
    [[nodiscard]] constexpr std::uint32_t nz() const noexcept { return nz_; }
    [[nodiscard]] constexpr Vec3 spacing() const noexcept { return spacing_; }
    [[nodiscard]] std::uint64_t cell_count() const noexcept;
    [[nodiscard]] std::uint64_t point_count() const;
    [[nodiscard]] std::uint64_t linear_id(const CellKey& key) const;
    [[nodiscard]] CellKey cell_key(std::uint64_t linear_id) const;
    [[nodiscard]] AABB cell_bounds(const CellKey& key) const;
    [[nodiscard]] Vec3 cell_center(const CellKey& key) const;
    [[nodiscard]] double cell_volume() const noexcept;

  private:
    void validate_key(const CellKey& key) const;

    AABB domain_;
    std::uint32_t nx_;
    std::uint32_t ny_;
    std::uint32_t nz_;
    Vec3 spacing_;
};

} // 命名空间 cartmesh
