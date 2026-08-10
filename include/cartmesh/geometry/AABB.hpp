#pragma once

#include "cartmesh/geometry/Vec3.hpp"

#include <stdexcept>

namespace cartmesh {

class AABB {
  public:
    AABB(Vec3 minimum, Vec3 maximum) : minimum_(minimum), maximum_(maximum) {
        if (!is_finite(minimum_) || !is_finite(maximum_) || maximum_.x < minimum_.x ||
            maximum_.y < minimum_.y || maximum_.z < minimum_.z) {
            throw std::invalid_argument("AABB 必须使用有限坐标，并且各方向范围不得为负");
        }
    }

    [[nodiscard]] constexpr const Vec3& minimum() const noexcept { return minimum_; }
    [[nodiscard]] constexpr const Vec3& maximum() const noexcept { return maximum_; }
    [[nodiscard]] constexpr Vec3 extent() const noexcept { return maximum_ - minimum_; }
    [[nodiscard]] constexpr Vec3 center() const noexcept {
        return minimum_ + extent() * 0.5;
    }
    [[nodiscard]] constexpr double volume() const noexcept {
        const auto size = extent();
        return size.x * size.y * size.z;
    }
    [[nodiscard]] constexpr bool has_positive_volume() const noexcept {
        const auto size = extent();
        return size.x > 0.0 && size.y > 0.0 && size.z > 0.0;
    }
    [[nodiscard]] constexpr bool contains(const Vec3& point) const noexcept {
        return point.x >= minimum_.x && point.x <= maximum_.x && point.y >= minimum_.y &&
               point.y <= maximum_.y && point.z >= minimum_.z && point.z <= maximum_.z;
    }
    [[nodiscard]] constexpr bool intersects(const AABB& other) const noexcept {
        return minimum_.x <= other.maximum_.x && maximum_.x >= other.minimum_.x &&
               minimum_.y <= other.maximum_.y && maximum_.y >= other.minimum_.y &&
               minimum_.z <= other.maximum_.z && maximum_.z >= other.minimum_.z;
    }

  private:
    Vec3 minimum_;
    Vec3 maximum_;
};

} // 命名空间 cartmesh
