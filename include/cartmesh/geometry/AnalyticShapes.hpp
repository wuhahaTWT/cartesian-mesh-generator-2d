#pragma once

#include "cartmesh/geometry/AABB.hpp"

#include <numbers>
#include <stdexcept>

namespace cartmesh {

class AnalyticCube {
  public:
    explicit AnalyticCube(AABB bounds) : bounds_(bounds) {}

    [[nodiscard]] constexpr const AABB& bounds() const noexcept { return bounds_; }
    [[nodiscard]] constexpr bool contains(const Vec3& point) const noexcept {
        return bounds_.contains(point);
    }
    [[nodiscard]] constexpr double volume() const noexcept { return bounds_.volume(); }
    [[nodiscard]] constexpr double surface_area() const noexcept {
        const auto size = bounds_.extent();
        return 2.0 * (size.x * size.y + size.y * size.z + size.z * size.x);
    }

  private:
    AABB bounds_;
};

class AnalyticSphere {
  public:
    AnalyticSphere(Vec3 center, double radius) : center_(center), radius_(radius) {
        if (!is_finite(center_) || !std::isfinite(radius_) || radius_ <= 0.0) {
            throw std::invalid_argument("球体必须具有有限坐标的球心和正半径");
        }
    }

    [[nodiscard]] constexpr const Vec3& center() const noexcept { return center_; }
    [[nodiscard]] constexpr double radius() const noexcept { return radius_; }
    [[nodiscard]] constexpr bool contains(const Vec3& point) const noexcept {
        const auto delta = point - center_;
        return dot(delta, delta) <= radius_ * radius_;
    }
    [[nodiscard]] double signed_distance(const Vec3& point) const noexcept {
        return norm(point - center_) - radius_;
    }
    [[nodiscard]] constexpr double volume() const noexcept {
        return (4.0 / 3.0) * std::numbers::pi_v<double> * radius_ * radius_ * radius_;
    }
    [[nodiscard]] constexpr double surface_area() const noexcept {
        return 4.0 * std::numbers::pi_v<double> * radius_ * radius_;
    }

  private:
    Vec3 center_;
    double radius_;
};

} // 命名空间 cartmesh
