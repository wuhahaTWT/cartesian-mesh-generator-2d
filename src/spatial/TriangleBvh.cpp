#include "cartmesh/spatial/TriangleBvh.hpp"

#include "cartmesh/geometry/TriangleBoxIntersection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cartmesh {
namespace {

[[nodiscard]] AABB triangle_range_bounds(const std::vector<Triangle>& triangles,
                                         const std::vector<std::uint64_t>& indices,
                                         std::size_t begin, std::size_t end) {
    Vec3 minimum = triangles[static_cast<std::size_t>(indices[begin])].vertices().front();
    Vec3 maximum = minimum;
    for (std::size_t position = begin; position < end; ++position) {
        for (const auto& vertex : triangles[static_cast<std::size_t>(indices[position])].vertices()) {
            minimum.x = std::min(minimum.x, vertex.x);
            minimum.y = std::min(minimum.y, vertex.y);
            minimum.z = std::min(minimum.z, vertex.z);
            maximum.x = std::max(maximum.x, vertex.x);
            maximum.y = std::max(maximum.y, vertex.y);
            maximum.z = std::max(maximum.z, vertex.z);
        }
    }
    return AABB(minimum, maximum);
}

[[nodiscard]] double axis_value(const Vec3& value, std::size_t axis) noexcept {
    if (axis == 0) {
        return value.x;
    }
    if (axis == 1) {
        return value.y;
    }
    return value.z;
}

[[nodiscard]] bool ray_intersects_aabb(const Vec3& origin, const Vec3& direction,
                                       const AABB& bounds) noexcept {
    double near_distance = 0.0;
    double far_distance = std::numeric_limits<double>::infinity();
    const std::array<double, 3> origins{origin.x, origin.y, origin.z};
    const std::array<double, 3> directions{direction.x, direction.y, direction.z};
    const std::array<double, 3> minimum{bounds.minimum().x, bounds.minimum().y,
                                        bounds.minimum().z};
    const std::array<double, 3> maximum{bounds.maximum().x, bounds.maximum().y,
                                        bounds.maximum().z};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (std::abs(directions[axis]) <= std::numeric_limits<double>::min()) {
            if (origins[axis] < minimum[axis] || origins[axis] > maximum[axis]) {
                return false;
            }
            continue;
        }
        const double inverse = 1.0 / directions[axis];
        double first = (minimum[axis] - origins[axis]) * inverse;
        double second = (maximum[axis] - origins[axis]) * inverse;
        if (first > second) {
            std::swap(first, second);
        }
        near_distance = std::max(near_distance, first);
        far_distance = std::min(far_distance, second);
        if (far_distance < near_distance) {
            return false;
        }
    }
    return far_distance >= 0.0;
}

struct TriangleRayHit {
    double distance{};
    bool ambiguous{};
};

[[nodiscard]] double coordinate_ulp(double value) noexcept {
    const double magnitude = std::abs(value);
    return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
}

[[nodiscard]] double coordinate_resolution(const Vec3& value) noexcept {
    return std::max({coordinate_ulp(value.x), coordinate_ulp(value.y),
                     coordinate_ulp(value.z)});
}

[[nodiscard]] bool ray_intersects_triangle(const Vec3& origin, const Vec3& direction,
                                           const Triangle& triangle, TriangleRayHit& hit) noexcept {
    const auto& vertices = triangle.vertices();
    const Vec3 edge1 = vertices[1] - vertices[0];
    const Vec3 edge2 = vertices[2] - vertices[0];
    const Vec3 p = cross(direction, edge2);
    const double determinant = dot(edge1, p);
    const double edge_scale = norm(edge1) * norm(edge2);
    const double determinant_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * edge_scale;
    if (std::abs(determinant) <= determinant_tolerance) {
        return false;
    }
    const double inverse = 1.0 / determinant;
    const Vec3 from_vertex = origin - vertices[0];
    const double u = dot(from_vertex, p) * inverse;
    const Vec3 q = cross(from_vertex, edge1);
    const double v = dot(direction, q) * inverse;
    const double distance = dot(edge2, q) * inverse;
    const double barycentric_tolerance = 256.0 * std::numeric_limits<double>::epsilon();
    const double distance_tolerance = std::max(
        64.0 * std::numeric_limits<double>::epsilon() *
            std::max({norm(edge1), norm(edge2), std::abs(distance)}),
        4.0 * std::max(coordinate_resolution(origin), coordinate_resolution(vertices[0])));
    if (u < -barycentric_tolerance || v < -barycentric_tolerance ||
        u + v > 1.0 + barycentric_tolerance || distance <= distance_tolerance) {
        return false;
    }
    hit.distance = distance;
    hit.ambiguous = u <= barycentric_tolerance || v <= barycentric_tolerance ||
                    1.0 - u - v <= barycentric_tolerance;
    return true;
}

[[nodiscard]] double point_triangle_distance_squared(const Vec3& point,
                                                     const Triangle& triangle) noexcept {
    const auto& vertices = triangle.vertices();
    const Vec3 ab = vertices[1] - vertices[0];
    const Vec3 ac = vertices[2] - vertices[0];
    const Vec3 ap = point - vertices[0];
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        return dot(ap, ap);
    }

    const Vec3 bp = point - vertices[1];
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        return dot(bp, bp);
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double fraction = d1 / (d1 - d3);
        const Vec3 projection = vertices[0] + ab * fraction;
        const Vec3 difference = point - projection;
        return dot(difference, difference);
    }

    const Vec3 cp = point - vertices[2];
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        return dot(cp, cp);
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double fraction = d2 / (d2 - d6);
        const Vec3 projection = vertices[0] + ac * fraction;
        const Vec3 difference = point - projection;
        return dot(difference, difference);
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
        const double fraction = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Vec3 projection = vertices[1] + (vertices[2] - vertices[1]) * fraction;
        const Vec3 difference = point - projection;
        return dot(difference, difference);
    }

    const Vec3 normal = cross(ab, ac);
    const double normal_squared = dot(normal, normal);
    if (normal_squared == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    const double signed_scaled_distance = dot(ap, normal);
    return signed_scaled_distance * signed_scaled_distance / normal_squared;
}

[[nodiscard]] double point_aabb_distance_squared(const Vec3& point,
                                                 const AABB& bounds) noexcept {
    const auto axis_distance = [](double value, double minimum, double maximum) {
        if (value < minimum) {
            return minimum - value;
        }
        if (value > maximum) {
            return value - maximum;
        }
        return 0.0;
    };
    const double dx = axis_distance(point.x, bounds.minimum().x, bounds.maximum().x);
    const double dy = axis_distance(point.y, bounds.minimum().y, bounds.maximum().y);
    const double dz = axis_distance(point.z, bounds.minimum().z, bounds.maximum().z);
    return dx * dx + dy * dy + dz * dz;
}

} // 匿名命名空间

TriangleBvh::TriangleBvh(const SurfaceMesh& surface, std::size_t leaf_size)
    : triangles_(surface.triangles()), indices_(triangles_.size()) {
    if (leaf_size == 0 || leaf_size > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("BVH 叶节点三角形数必须是 32 位正整数");
    }
    std::iota(indices_.begin(), indices_.end(), std::uint64_t{0});
    nodes_.reserve(triangles_.size() * 2);
    static_cast<void>(build(0, indices_.size(), 0, leaf_size));
    statistics_.triangle_count = triangles_.size();
    statistics_.node_count = nodes_.size();
}

std::uint32_t TriangleBvh::build(std::size_t begin, std::size_t end, std::uint64_t depth,
                                 std::size_t leaf_size) {
    const auto bounds_value = triangle_range_bounds(triangles_, indices_, begin, end);
    const auto node_index = static_cast<std::uint32_t>(nodes_.size());
    nodes_.emplace_back(bounds_value);
    statistics_.maximum_depth = std::max(statistics_.maximum_depth, depth);
    const auto count = end - begin;
    if (count <= leaf_size) {
        nodes_[node_index].first = begin;
        nodes_[node_index].count = static_cast<std::uint32_t>(count);
        nodes_[node_index].leaf = true;
        ++statistics_.leaf_count;
        statistics_.maximum_leaf_triangles =
            std::max(statistics_.maximum_leaf_triangles, static_cast<std::uint64_t>(count));
        return node_index;
    }

    Vec3 centroid_minimum = triangles_[static_cast<std::size_t>(indices_[begin])].centroid();
    Vec3 centroid_maximum = centroid_minimum;
    for (std::size_t position = begin + 1; position < end; ++position) {
        const auto centroid = triangles_[static_cast<std::size_t>(indices_[position])].centroid();
        centroid_minimum.x = std::min(centroid_minimum.x, centroid.x);
        centroid_minimum.y = std::min(centroid_minimum.y, centroid.y);
        centroid_minimum.z = std::min(centroid_minimum.z, centroid.z);
        centroid_maximum.x = std::max(centroid_maximum.x, centroid.x);
        centroid_maximum.y = std::max(centroid_maximum.y, centroid.y);
        centroid_maximum.z = std::max(centroid_maximum.z, centroid.z);
    }
    const auto centroid_extent = centroid_maximum - centroid_minimum;
    std::size_t axis = 0;
    if (centroid_extent.y > centroid_extent.x) {
        axis = 1;
    }
    if (axis_value(centroid_extent, 2) > axis_value(centroid_extent, axis)) {
        axis = 2;
    }
    std::stable_sort(indices_.begin() + static_cast<std::ptrdiff_t>(begin),
                     indices_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::uint64_t lhs, std::uint64_t rhs) {
                         const double lhs_value =
                             axis_value(triangles_[static_cast<std::size_t>(lhs)].centroid(), axis);
                         const double rhs_value =
                             axis_value(triangles_[static_cast<std::size_t>(rhs)].centroid(), axis);
                         return lhs_value < rhs_value || (lhs_value == rhs_value && lhs < rhs);
                     });
    const auto middle = begin + count / 2;
    const auto left = build(begin, middle, depth + 1, leaf_size);
    const auto right = build(middle, end, depth + 1, leaf_size);
    nodes_[node_index].left = left;
    nodes_[node_index].right = right;
    return node_index;
}

std::vector<std::uint64_t> TriangleBvh::query(const AABB& bounds_value) const {
    std::vector<std::uint64_t> result;
    std::vector<std::uint32_t> stack{0};
    while (!stack.empty()) {
        const auto node_index = stack.back();
        stack.pop_back();
        const auto& node = nodes_[node_index];
        if (!node.bounds.intersects(bounds_value)) {
            continue;
        }
        if (node.leaf) {
            for (std::uint32_t offset = 0; offset < node.count; ++offset) {
                const auto triangle_index = indices_[static_cast<std::size_t>(node.first + offset)];
                if (triangles_[static_cast<std::size_t>(triangle_index)].bounds().intersects(
                        bounds_value)) {
                    result.push_back(triangle_index);
                }
            }
        } else {
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

bool TriangleBvh::intersects_surface(const AABB& bounds_value) const {
    const double world_scale =
        std::max({std::abs(bounds_value.minimum().x), std::abs(bounds_value.minimum().y),
                  std::abs(bounds_value.minimum().z), std::abs(bounds_value.maximum().x),
                  std::abs(bounds_value.maximum().y), std::abs(bounds_value.maximum().z),
                  std::abs(bounds().minimum().x), std::abs(bounds().minimum().y),
                  std::abs(bounds().minimum().z), std::abs(bounds().maximum().x),
                  std::abs(bounds().maximum().y), std::abs(bounds().maximum().z),
                  norm(bounds_value.extent()), norm(bounds().extent()),
                  std::numeric_limits<double>::min()});
    const double candidate_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * world_scale;
    const Vec3 padding{candidate_tolerance, candidate_tolerance, candidate_tolerance};
    const AABB candidate_bounds(bounds_value.minimum() - padding,
                                bounds_value.maximum() + padding);
    std::vector<std::uint32_t> stack{0};
    while (!stack.empty()) {
        const auto node_index = stack.back();
        stack.pop_back();
        const auto& node = nodes_[node_index];
        if (!node.bounds.intersects(candidate_bounds)) {
            continue;
        }
        if (node.leaf) {
            for (std::uint32_t offset = 0; offset < node.count; ++offset) {
                const auto triangle_index = indices_[static_cast<std::size_t>(node.first + offset)];
                const auto& triangle = triangles_[static_cast<std::size_t>(triangle_index)];
                if (triangle.bounds().intersects(candidate_bounds) &&
                    triangle_intersects_aabb(triangle, bounds_value)) {
                    return true;
                }
            }
        } else {
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }
    return false;
}

RayIntersectionSummary TriangleBvh::trace_ray(const Vec3& origin, const Vec3& direction) const {
    if (!is_finite(origin) || !is_finite(direction) || norm(direction) == 0.0) {
        throw std::invalid_argument("BVH 射线必须具有有限原点和非零有限方向");
    }
    std::vector<TriangleRayHit> hits;
    std::vector<std::uint32_t> stack{0};
    while (!stack.empty()) {
        const auto node_index = stack.back();
        stack.pop_back();
        const auto& node = nodes_[node_index];
        if (!ray_intersects_aabb(origin, direction, node.bounds)) {
            continue;
        }
        if (node.leaf) {
            for (std::uint32_t offset = 0; offset < node.count; ++offset) {
                const auto triangle_index = indices_[static_cast<std::size_t>(node.first + offset)];
                TriangleRayHit hit;
                if (ray_intersects_triangle(origin, direction,
                                            triangles_[static_cast<std::size_t>(triangle_index)],
                                            hit)) {
                    hits.push_back(hit);
                }
            }
        } else {
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }
    std::sort(hits.begin(), hits.end(),
              [](const TriangleRayHit& lhs, const TriangleRayHit& rhs) {
                  return lhs.distance < rhs.distance;
              });
    RayIntersectionSummary result;
    result.raw_triangle_hits = hits.size();
    const double local_scale = norm(bounds().extent());
    const double input_resolution = std::max(
        coordinate_resolution(origin),
        std::max(coordinate_resolution(bounds().minimum()),
                 coordinate_resolution(bounds().maximum())));
    std::size_t begin = 0;
    while (begin < hits.size()) {
        std::size_t end = begin + 1;
        const double tolerance = std::max(
            512.0 * std::numeric_limits<double>::epsilon() *
                std::max(local_scale, std::abs(hits[begin].distance)),
            8.0 * input_resolution);
        bool group_ambiguous = hits[begin].ambiguous;
        while (end < hits.size() && hits[end].distance - hits[begin].distance <= tolerance) {
            group_ambiguous = group_ambiguous || hits[end].ambiguous;
            ++end;
        }
        ++result.unique_surface_crossings;
        result.ambiguous = result.ambiguous || group_ambiguous;
        begin = end;
    }
    return result;
}

bool TriangleBvh::point_on_surface(const Vec3& point, double tolerance) const {
    if (!is_finite(point) || !std::isfinite(tolerance) || tolerance < 0.0) {
        throw std::invalid_argument("表面距离查询需要有限点和非负有限容差");
    }
    const Vec3 padding{tolerance, tolerance, tolerance};
    const auto candidates = query(AABB(point - padding, point + padding));
    const double squared_tolerance = tolerance * tolerance;
    for (const auto triangle_index : candidates) {
        if (point_triangle_distance_squared(
                point, triangles_[static_cast<std::size_t>(triangle_index)]) <=
            squared_tolerance) {
            return true;
        }
    }
    return false;
}

double TriangleBvh::distance_to_surface(const Vec3& point) const {
    if (!is_finite(point)) {
        throw std::invalid_argument("表面最短距离查询需要有限坐标");
    }
    double best_squared = std::numeric_limits<double>::infinity();
    std::vector<std::uint32_t> stack{0};
    while (!stack.empty()) {
        const auto node_index = stack.back();
        stack.pop_back();
        const auto& node = nodes_[node_index];
        if (point_aabb_distance_squared(point, node.bounds) > best_squared) {
            continue;
        }
        if (node.leaf) {
            for (std::uint32_t offset = 0; offset < node.count; ++offset) {
                const auto triangle_index = indices_[static_cast<std::size_t>(node.first + offset)];
                best_squared = std::min(
                    best_squared,
                    point_triangle_distance_squared(
                        point, triangles_[static_cast<std::size_t>(triangle_index)]));
            }
            continue;
        }
        const double left_distance = point_aabb_distance_squared(point, nodes_[node.left].bounds);
        const double right_distance =
            point_aabb_distance_squared(point, nodes_[node.right].bounds);
        if (left_distance <= right_distance) {
            stack.push_back(node.right);
            stack.push_back(node.left);
        } else {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }
    return std::sqrt(best_squared);
}

} // 命名空间 cartmesh
