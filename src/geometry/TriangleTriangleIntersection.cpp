#include "cartmesh/geometry/TriangleTriangleIntersection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace cartmesh {
namespace {

struct Vec2 {
    double x{};
    double y{};
};

[[nodiscard]] double cross2(const Vec2& first, const Vec2& second) noexcept {
    return first.x * second.y - first.y * second.x;
}

[[nodiscard]] Vec2 project(const Vec3& point, std::size_t dropped_axis) noexcept {
    if (dropped_axis == 0) return {point.y, point.z};
    if (dropped_axis == 1) return {point.x, point.z};
    return {point.x, point.y};
}

[[nodiscard]] double maximum_edge_length(const Triangle& triangle) noexcept {
    const auto& point = triangle.vertices();
    return std::max({norm(point[1] - point[0]), norm(point[2] - point[1]),
                     norm(point[0] - point[2])});
}

[[nodiscard]] bool point_in_triangle(const Vec3& point,
                                     const Triangle& triangle,
                                     double dimensionless_tolerance) noexcept {
    const auto& vertex = triangle.vertices();
    const Vec3 first = vertex[1] - vertex[0];
    const Vec3 second = vertex[2] - vertex[0];
    const Vec3 local = point - vertex[0];
    const double first_first = dot(first, first);
    const double first_second = dot(first, second);
    const double second_second = dot(second, second);
    const double local_first = dot(local, first);
    const double local_second = dot(local, second);
    const double denominator =
        first_first * second_second - first_second * first_second;
    if (!(denominator > 0.0)) return false;
    const double u =
        (second_second * local_first - first_second * local_second) /
        denominator;
    const double v =
        (first_first * local_second - first_second * local_first) /
        denominator;
    return u >= -dimensionless_tolerance && v >= -dimensionless_tolerance &&
           u + v <= 1.0 + dimensionless_tolerance;
}

void append_unique(std::vector<Vec3>& points, const Vec3& point,
                   double tolerance) noexcept {
    for (const auto& existing : points) {
        if (norm(existing - point) <= tolerance) return;
    }
    points.push_back(point);
}

void collect_plane_intersections(const Triangle& source,
                                 const Triangle& target,
                                 const Vec3& target_normal,
                                 double target_offset, double tolerance,
                                 double dimensionless_tolerance,
                                 std::vector<Vec3>& intersections) noexcept {
    const auto& vertex = source.vertices();
    std::array<double, 3> distance{};
    for (std::size_t index = 0; index < 3; ++index) {
        distance[index] = dot(target_normal, vertex[index]) - target_offset;
        if (std::abs(distance[index]) <= tolerance &&
            point_in_triangle(vertex[index], target, dimensionless_tolerance)) {
            append_unique(intersections, vertex[index], tolerance);
        }
    }
    for (std::size_t edge = 0; edge < 3; ++edge) {
        const std::size_t next = (edge + 1U) % 3U;
        if ((distance[edge] < -tolerance && distance[next] > tolerance) ||
            (distance[edge] > tolerance && distance[next] < -tolerance)) {
            const double fraction =
                distance[edge] / (distance[edge] - distance[next]);
            const Vec3 point =
                vertex[edge] + (vertex[next] - vertex[edge]) * fraction;
            if (point_in_triangle(point, target, dimensionless_tolerance)) {
                append_unique(intersections, point, tolerance);
            }
        }
    }
}

[[nodiscard]] std::vector<Vec2> clip_coplanar_triangle(
    const Triangle& subject_triangle, const Triangle& clip_triangle,
    std::size_t dropped_axis, double area_tolerance) noexcept {
    const Vec3 origin = subject_triangle.vertices()[0];
    std::vector<Vec2> polygon;
    for (const auto& point : subject_triangle.vertices()) {
        polygon.push_back(project(point - origin, dropped_axis));
    }
    std::array<Vec2, 3> clip{};
    for (std::size_t index = 0; index < 3; ++index) {
        clip[index] =
            project(clip_triangle.vertices()[index] - origin, dropped_axis);
    }
    const double orientation =
        cross2({clip[1].x - clip[0].x, clip[1].y - clip[0].y},
               {clip[2].x - clip[0].x, clip[2].y - clip[0].y}) >= 0.0
            ? 1.0
            : -1.0;
    for (std::size_t edge = 0; edge < 3 && !polygon.empty(); ++edge) {
        const Vec2 first = clip[edge];
        const Vec2 second = clip[(edge + 1U) % 3U];
        const Vec2 direction{second.x - first.x, second.y - first.y};
        const auto signed_distance = [&](const Vec2& point) {
            return orientation *
                   cross2(direction, {point.x - first.x, point.y - first.y});
        };
        std::vector<Vec2> clipped;
        for (std::size_t index = 0; index < polygon.size(); ++index) {
            const Vec2 current = polygon[index];
            const Vec2 next = polygon[(index + 1U) % polygon.size()];
            const double current_distance = signed_distance(current);
            const double next_distance = signed_distance(next);
            const bool current_inside = current_distance >= -area_tolerance;
            const bool next_inside = next_distance >= -area_tolerance;
            if (current_inside) clipped.push_back(current);
            if (current_inside != next_inside) {
                const double denominator = current_distance - next_distance;
                if (denominator != 0.0) {
                    const double fraction = current_distance / denominator;
                    clipped.push_back(
                        {current.x + (next.x - current.x) * fraction,
                         current.y + (next.y - current.y) * fraction});
                }
            }
        }
        polygon = std::move(clipped);
    }
    return polygon;
}

[[nodiscard]] double polygon_area(const std::vector<Vec2>& polygon) noexcept {
    double twice_area = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1U) % polygon.size()];
        twice_area += first.x * second.y - first.y * second.x;
    }
    return 0.5 * std::abs(twice_area);
}

} // 匿名命名空间

TriangleTriangleRelation classify_triangle_triangle(
    const Triangle& first, const Triangle& second,
    double length_tolerance) noexcept {
    const double scale =
        std::max({maximum_edge_length(first), maximum_edge_length(second),
                  std::numeric_limits<double>::min()});
    const double tolerance = std::max(
        length_tolerance,
        256.0 * std::numeric_limits<double>::epsilon() * scale);
    const double dimensionless_tolerance = tolerance / scale;
    const Vec3 first_area_vector = first.area_vector();
    const Vec3 second_area_vector = second.area_vector();
    const double first_area = norm(first_area_vector);
    const double second_area = norm(second_area_vector);
    if (first_area <= 0.5 * scale * tolerance ||
        second_area <= 0.5 * scale * tolerance) {
        return TriangleTriangleRelation::indeterminate;
    }
    const Vec3 first_normal = first_area_vector / first_area;
    const Vec3 second_normal = second_area_vector / second_area;
    const double first_offset = dot(first_normal, first.vertices()[0]);
    const double second_offset = dot(second_normal, second.vertices()[0]);
    const auto separated_by_plane = [&](const Triangle& triangle,
                                        const Vec3& normal, double offset) {
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        for (const auto& point : triangle.vertices()) {
            const double distance = dot(normal, point) - offset;
            minimum = std::min(minimum, distance);
            maximum = std::max(maximum, distance);
        }
        return minimum > tolerance || maximum < -tolerance;
    };
    if (separated_by_plane(second, first_normal, first_offset) ||
        separated_by_plane(first, second_normal, second_offset)) {
        return TriangleTriangleRelation::disjoint;
    }

    const double cross_normal_length = norm(cross(first_normal, second_normal));
    if (cross_normal_length <= dimensionless_tolerance) {
        if (std::abs(dot(first_normal, second.vertices()[0]) - first_offset) >
            tolerance) {
            return TriangleTriangleRelation::disjoint;
        }
        const Vec3 absolute_normal{std::abs(first_normal.x),
                                   std::abs(first_normal.y),
                                   std::abs(first_normal.z)};
        std::size_t dropped_axis = 0;
        if (absolute_normal.y > absolute_normal.x) dropped_axis = 1;
        if ((dropped_axis == 0 ? absolute_normal.x : absolute_normal.y) <
            absolute_normal.z) {
            dropped_axis = 2;
        }
        const auto polygon = clip_coplanar_triangle(
            first, second, dropped_axis, tolerance * scale);
        if (polygon.empty()) return TriangleTriangleRelation::disjoint;
        if (polygon_area(polygon) > tolerance * scale) {
            return TriangleTriangleRelation::coplanar_area_overlap;
        }
        return TriangleTriangleRelation::boundary_contact;
    }

    std::vector<Vec3> intersections;
    collect_plane_intersections(first, second, second_normal, second_offset,
                                tolerance, dimensionless_tolerance,
                                intersections);
    collect_plane_intersections(second, first, first_normal, first_offset,
                                tolerance, dimensionless_tolerance,
                                intersections);
    if (intersections.empty()) return TriangleTriangleRelation::disjoint;
    for (std::size_t first_point = 0; first_point < intersections.size();
         ++first_point) {
        for (std::size_t second_point = first_point + 1;
             second_point < intersections.size(); ++second_point) {
            if (norm(intersections[first_point] - intersections[second_point]) >
                tolerance) {
                return TriangleTriangleRelation::proper_intersection;
            }
        }
    }
    return TriangleTriangleRelation::boundary_contact;
}

} // 命名空间 cartmesh
