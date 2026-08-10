#include "cartmesh/grid/UniformCartesianGrid.hpp"

#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

[[nodiscard]] std::uint64_t checked_product(std::uint64_t lhs, std::uint64_t rhs) {
    if (rhs != 0 && lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
        throw std::overflow_error("Cartesian 网格规模超出 64 位索引范围");
    }
    return lhs * rhs;
}

} // 匿名命名空间

UniformCartesianGrid::UniformCartesianGrid(AABB domain, std::uint32_t nx, std::uint32_t ny,
                                           std::uint32_t nz)
    : domain_(domain), nx_(nx), ny_(ny), nz_(nz), spacing_{} {
    if (!domain_.has_positive_volume()) {
        throw std::invalid_argument("Cartesian 网格计算域必须具有正体积");
    }
    if (nx_ == 0 || ny_ == 0 || nz_ == 0) {
        throw std::invalid_argument("Cartesian 网格三个方向的单元数都必须为正数");
    }
    static_cast<void>(checked_product(checked_product(nx_, ny_), nz_));
    static_cast<void>(checked_product(checked_product(static_cast<std::uint64_t>(nx_) + 1,
                                                      static_cast<std::uint64_t>(ny_) + 1),
                                      static_cast<std::uint64_t>(nz_) + 1));
    const auto size = domain_.extent();
    spacing_ = {size.x / static_cast<double>(nx_), size.y / static_cast<double>(ny_),
                size.z / static_cast<double>(nz_)};
}

std::uint64_t UniformCartesianGrid::cell_count() const noexcept {
    return static_cast<std::uint64_t>(nx_) * static_cast<std::uint64_t>(ny_) *
           static_cast<std::uint64_t>(nz_);
}

std::uint64_t UniformCartesianGrid::point_count() const {
    return checked_product(
        checked_product(static_cast<std::uint64_t>(nx_) + 1,
                        static_cast<std::uint64_t>(ny_) + 1),
        static_cast<std::uint64_t>(nz_) + 1);
}

void UniformCartesianGrid::validate_key(const CellKey& key) const {
    if (key.level != 0 || key.i >= nx_ || key.j >= ny_ || key.k >= nz_) {
        throw std::out_of_range("CellKey 位于当前 0 层均匀网格之外");
    }
}

std::uint64_t UniformCartesianGrid::linear_id(const CellKey& key) const {
    validate_key(key);
    return static_cast<std::uint64_t>(key.i) +
           static_cast<std::uint64_t>(nx_) *
               (static_cast<std::uint64_t>(key.j) +
                static_cast<std::uint64_t>(ny_) * static_cast<std::uint64_t>(key.k));
}

CellKey UniformCartesianGrid::cell_key(std::uint64_t linear_id_value) const {
    if (linear_id_value >= cell_count()) {
        throw std::out_of_range("线性单元 ID 位于网格范围之外");
    }
    const auto plane = static_cast<std::uint64_t>(nx_) * static_cast<std::uint64_t>(ny_);
    const auto k = linear_id_value / plane;
    const auto remainder = linear_id_value % plane;
    const auto j = remainder / nx_;
    const auto i = remainder % nx_;
    return {0, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j),
            static_cast<std::uint32_t>(k)};
}

AABB UniformCartesianGrid::cell_bounds(const CellKey& key) const {
    validate_key(key);
    const auto base = domain_.minimum();
    const Vec3 minimum{base.x + spacing_.x * static_cast<double>(key.i),
                       base.y + spacing_.y * static_cast<double>(key.j),
                       base.z + spacing_.z * static_cast<double>(key.k)};
    const Vec3 maximum{
        base.x + spacing_.x * static_cast<double>(key.i + 1U),
        base.y + spacing_.y * static_cast<double>(key.j + 1U),
        base.z + spacing_.z * static_cast<double>(key.k + 1U)};
    return AABB(minimum, maximum);
}

Vec3 UniformCartesianGrid::cell_center(const CellKey& key) const {
    return cell_bounds(key).center();
}

double UniformCartesianGrid::cell_volume() const noexcept {
    return spacing_.x * spacing_.y * spacing_.z;
}

} // 命名空间 cartmesh
