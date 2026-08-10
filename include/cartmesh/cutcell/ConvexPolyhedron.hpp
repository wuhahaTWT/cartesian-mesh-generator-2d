#pragma once

#include "cartmesh/geometry/AABB.hpp"

#include <cstdint>
#include <vector>

namespace cartmesh {

enum class PolyhedronFaceKind : std::uint8_t {
    cartesian = 0,
    embedded_boundary = 1,
    internal_partition = 2,
};

struct PolyhedronFace {
    std::vector<std::uint32_t> vertex_indices;
    PolyhedronFaceKind kind{PolyhedronFaceKind::cartesian};
    std::uint64_t source_id{};
};

struct ConvexPolyhedron {
    std::vector<Vec3> vertices;
    std::vector<PolyhedronFace> faces;

    [[nodiscard]] bool empty() const noexcept { return faces.empty(); }
};

// 保留 dot(unit_normal, position - point) <= 0 的闭半空间。
class OrientedHalfSpace {
  public:
    OrientedHalfSpace(Vec3 point, Vec3 outward_normal, std::uint64_t boundary_id = 0);

    [[nodiscard]] constexpr const Vec3& point() const noexcept { return point_; }
    [[nodiscard]] constexpr const Vec3& outward_normal() const noexcept {
        return outward_normal_;
    }
    [[nodiscard]] constexpr std::uint64_t boundary_id() const noexcept {
        return boundary_id_;
    }
    [[nodiscard]] constexpr double signed_distance(const Vec3& position) const noexcept {
        return dot(distance_normal_, position - point_);
    }
    [[nodiscard]] constexpr double distance_scale() const noexcept {
        return distance_scale_;
    }

  private:
    Vec3 point_{};
    Vec3 distance_normal_{};
    Vec3 outward_normal_{};
    double distance_scale_{};
    std::uint64_t boundary_id_{};
};

struct PolygonGeometry {
    double area{};
    Vec3 centroid{};
    Vec3 outward_normal{};
    Vec3 area_vector{};
};

struct PolyhedronGeometry {
    double volume{};
    Vec3 centroid{};
    std::vector<PolygonGeometry> faces;
    Vec3 oriented_area_vector_sum{};
    bool closed{};
    bool positive_volume{};
};

[[nodiscard]] ConvexPolyhedron make_box_polyhedron(const AABB& box);

// 构造具有外向面绕序的非退化四面体；输入顶点顺序可以任意。
[[nodiscard]] ConvexPolyhedron make_tetrahedron_polyhedron(
    const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

// 对凸多面体执行确定性半空间裁剪。新生切面的外法向为
// half_space.outward_normal()，并标记为 embedded_boundary。
[[nodiscard]] ConvexPolyhedron clip_convex_polyhedron(
    const ConvexPolyhedron& polyhedron, const OrientedHalfSpace& half_space,
    double length_tolerance = 0.0,
    PolyhedronFaceKind cut_face_kind = PolyhedronFaceKind::embedded_boundary);

[[nodiscard]] PolyhedronGeometry measure_polyhedron(const ConvexPolyhedron& polyhedron);

} // 命名空间 cartmesh
