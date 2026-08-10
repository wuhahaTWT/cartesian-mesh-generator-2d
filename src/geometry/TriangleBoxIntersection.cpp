#include "cartmesh/geometry/TriangleBoxIntersection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace cartmesh {
namespace {

[[nodiscard]] double coordinate_ulp(double value) noexcept {
    const double upward = std::nextafter(value, std::numeric_limits<double>::infinity());
    const double downward = std::nextafter(value, -std::numeric_limits<double>::infinity());
    double spacing = 0.0;
    if (std::isfinite(upward)) {
        spacing = std::max(spacing, std::abs(upward - value));
    }
    if (std::isfinite(downward)) {
        spacing = std::max(spacing, std::abs(value - downward));
    }
    return spacing;
}

[[nodiscard]] bool projections_overlap(const Vec3& axis, const std::array<Vec3, 3>& vertices,
                                       const Vec3& half_extent,
                                       double length_tolerance) noexcept {
    const double axis_squared = dot(axis, axis);
    if (axis_squared == 0.0) {
        return true;
    }

    const std::array<double, 3> projection = {
        dot(axis, vertices[0]), dot(axis, vertices[1]), dot(axis, vertices[2])};
    const auto [minimum, maximum] = std::minmax_element(projection.begin(), projection.end());
    const double radius = half_extent.x * std::abs(axis.x) +
                          half_extent.y * std::abs(axis.y) +
                          half_extent.z * std::abs(axis.z);
    const double scale = std::max({std::abs(*minimum), std::abs(*maximum), std::abs(radius),
                                   std::numeric_limits<double>::min()});
    const double tolerance =
        std::max(64.0 * std::numeric_limits<double>::epsilon() * scale,
                 length_tolerance * std::sqrt(axis_squared));
    return *minimum <= radius + tolerance && *maximum >= -radius - tolerance;
}

} // 匿名命名空间

bool triangle_intersects_aabb(const Triangle& triangle, const AABB& box) noexcept {
    const Vec3 center = box.center();
    const Vec3 half_extent = box.extent() * 0.5;
    const auto& source = triangle.vertices();
    const std::array<Vec3, 3> vertices = {source[0] - center, source[1] - center,
                                          source[2] - center};
    const std::array<Vec3, 3> edges = {vertices[1] - vertices[0],
                                       vertices[2] - vertices[1],
                                       vertices[0] - vertices[2]};
    const std::array<Vec3, 3> box_axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const std::array<double, 15> input_coordinates = {
        box.minimum().x, box.minimum().y, box.minimum().z, box.maximum().x,
        box.maximum().y, box.maximum().z, source[0].x, source[0].y, source[0].z,
        source[1].x, source[1].y, source[1].z, source[2].x, source[2].y,
        source[2].z};
    double maximum_coordinate_ulp = 0.0;
    for (const double coordinate : input_coordinates) {
        maximum_coordinate_ulp =
            std::max(maximum_coordinate_ulp, coordinate_ulp(coordinate));
    }
    const double local_scale =
        std::max({norm(box.extent()), norm(edges[0]), norm(edges[1]), norm(edges[2]),
                  std::numeric_limits<double>::min()});
    // The vertices are translated by the box center before projection.  Four
    // input-coordinate ULPs cover that subtraction without making tolerance
    // grow by an arbitrary 64 ULPs with distance from the global origin.
    const double length_tolerance =
        std::max(4.0 * maximum_coordinate_ulp,
                 64.0 * std::numeric_limits<double>::epsilon() * local_scale);

    for (const auto& axis : box_axes) {
        if (!projections_overlap(axis, vertices, half_extent, length_tolerance)) {
            return false;
        }
    }

    const Vec3 triangle_normal = cross(edges[0], edges[1]);
    if (!projections_overlap(triangle_normal, vertices, half_extent, length_tolerance)) {
        return false;
    }

    for (const auto& edge : edges) {
        for (const auto& box_axis : box_axes) {
            if (!projections_overlap(cross(edge, box_axis), vertices, half_extent,
                                     length_tolerance)) {
                return false;
            }
        }
    }
    return true;
}

} // 命名空间 cartmesh
