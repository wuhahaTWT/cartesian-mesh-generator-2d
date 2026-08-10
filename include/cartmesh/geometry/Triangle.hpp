#pragma once

#include "cartmesh/geometry/AABB.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace cartmesh {

class Triangle {
  public:
    Triangle(Vec3 a, Vec3 b, Vec3 c) : vertices_{a, b, c} {
        if (!is_finite(a) || !is_finite(b) || !is_finite(c)) {
            throw std::invalid_argument("三角形顶点坐标必须是有限数值");
        }
    }

    [[nodiscard]] constexpr const std::array<Vec3, 3>& vertices() const noexcept {
        return vertices_;
    }
    [[nodiscard]] Vec3 area_vector() const noexcept {
        return cross(vertices_[1] - vertices_[0], vertices_[2] - vertices_[0]) * 0.5;
    }
    [[nodiscard]] double area() const noexcept { return norm(area_vector()); }
    [[nodiscard]] Vec3 centroid() const noexcept {
        return (vertices_[0] + vertices_[1] + vertices_[2]) / 3.0;
    }
    [[nodiscard]] AABB bounds() const {
        const Vec3 minimum{std::min({vertices_[0].x, vertices_[1].x, vertices_[2].x}),
                           std::min({vertices_[0].y, vertices_[1].y, vertices_[2].y}),
                           std::min({vertices_[0].z, vertices_[1].z, vertices_[2].z})};
        const Vec3 maximum{std::max({vertices_[0].x, vertices_[1].x, vertices_[2].x}),
                           std::max({vertices_[0].y, vertices_[1].y, vertices_[2].y}),
                           std::max({vertices_[0].z, vertices_[1].z, vertices_[2].z})};
        return AABB(minimum, maximum);
    }

  private:
    std::array<Vec3, 3> vertices_;
};

} // 命名空间 cartmesh
