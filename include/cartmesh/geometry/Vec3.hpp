#pragma once

#include <cmath>

namespace cartmesh {

struct Vec3 {
    double x{};
    double y{};
    double z{};

    [[nodiscard]] constexpr Vec3 operator+(const Vec3& rhs) const noexcept {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }
    [[nodiscard]] constexpr Vec3 operator-(const Vec3& rhs) const noexcept {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }
    [[nodiscard]] constexpr Vec3 operator*(double scale) const noexcept {
        return {x * scale, y * scale, z * scale};
    }
    [[nodiscard]] constexpr Vec3 operator/(double scale) const noexcept {
        return {x / scale, y / scale, z / scale};
    }
};

[[nodiscard]] constexpr double dot(const Vec3& lhs, const Vec3& rhs) noexcept {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

[[nodiscard]] constexpr Vec3 cross(const Vec3& lhs, const Vec3& rhs) noexcept {
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

[[nodiscard]] inline double norm(const Vec3& value) noexcept { return std::sqrt(dot(value, value)); }

[[nodiscard]] inline bool is_finite(const Vec3& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

} // 命名空间 cartmesh
