#include "cartmesh/grid/OctreeRefinement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace cartmesh {
namespace {

[[nodiscard]] double point_aabb_distance_squared(const Vec3& point,
                                                 const AABB& bounds) noexcept {
    const auto axis_distance = [](double value, double minimum, double maximum) {
        return value < minimum ? minimum - value
                               : (value > maximum ? value - maximum : 0.0);
    };
    const double dx = axis_distance(point.x, bounds.minimum().x, bounds.maximum().x);
    const double dy = axis_distance(point.y, bounds.minimum().y, bounds.maximum().y);
    const double dz = axis_distance(point.z, bounds.minimum().z, bounds.maximum().z);
    return dx * dx + dy * dy + dz * dz;
}

[[nodiscard]] double half_diagonal(const AABB& bounds) noexcept {
    return 0.5 * norm(bounds.extent());
}

[[nodiscard]] AABB expanded_bounds(const AABB& bounds, double distance) {
    const Vec3 padding{distance, distance, distance};
    return AABB(bounds.minimum() - padding, bounds.maximum() + padding);
}

[[nodiscard]] bool sphere_intersects_box(const SphereRefinementRegion& sphere,
                                         const AABB& bounds) noexcept {
    return point_aabb_distance_squared(sphere.center, bounds) <= sphere.radius * sphere.radius;
}

[[nodiscard]] bool cylinder_may_intersect_box(const CylinderRefinementRegion& cylinder,
                                              const AABB& bounds) noexcept {
    const Vec3 axis = cylinder.second_axis_point - cylinder.first_axis_point;
    const double length = norm(axis);
    const Vec3 direction = axis / length;
    const Vec3 from_first = bounds.center() - cylinder.first_axis_point;
    const double axial_coordinate = dot(from_first, direction);
    const double padding = half_diagonal(bounds);
    if (axial_coordinate < -padding || axial_coordinate > length + padding) {
        return false;
    }
    const double clamped_axial = std::clamp(axial_coordinate, 0.0, length);
    const Vec3 closest_axis_point = cylinder.first_axis_point + direction * clamped_axial;
    return norm(bounds.center() - closest_axis_point) <= cylinder.radius + padding;
}

[[nodiscard]] Vec3 unit_normal(const Triangle& triangle) noexcept {
    const Vec3 normal = triangle.area_vector();
    return normal / norm(normal);
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
        const Vec3 difference = point - (vertices[0] + ab * (d1 / (d1 - d3)));
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
        const Vec3 difference = point - (vertices[0] + ac * (d2 / (d2 - d6)));
        return dot(difference, difference);
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
        const Vec3 difference =
            point - (vertices[1] + (vertices[2] - vertices[1]) *
                                       ((d4 - d3) / ((d4 - d3) + (d5 - d6))));
        return dot(difference, difference);
    }
    const Vec3 normal = cross(ab, ac);
    const double signed_scaled_distance = dot(ap, normal);
    return signed_scaled_distance * signed_scaled_distance / dot(normal, normal);
}

[[nodiscard]] double segment_segment_distance_squared(const Vec3& first_start,
                                                      const Vec3& first_end,
                                                      const Vec3& second_start,
                                                      const Vec3& second_end) noexcept {
    const Vec3 first_direction = first_end - first_start;
    const Vec3 second_direction = second_end - second_start;
    const Vec3 separation = first_start - second_start;
    const double a = dot(first_direction, first_direction);
    const double e = dot(second_direction, second_direction);
    const double f = dot(second_direction, separation);
    double first_parameter = 0.0;
    double second_parameter = 0.0;
    if (a <= std::numeric_limits<double>::min() &&
        e <= std::numeric_limits<double>::min()) {
        return dot(separation, separation);
    }
    if (a <= std::numeric_limits<double>::min()) {
        second_parameter = std::clamp(f / e, 0.0, 1.0);
    } else {
        const double c = dot(first_direction, separation);
        if (e <= std::numeric_limits<double>::min()) {
            first_parameter = std::clamp(-c / a, 0.0, 1.0);
        } else {
            const double b = dot(first_direction, second_direction);
            const double denominator = a * e - b * b;
            if (denominator != 0.0) {
                first_parameter = std::clamp((b * f - c * e) / denominator, 0.0, 1.0);
            }
            second_parameter = (b * first_parameter + f) / e;
            if (second_parameter < 0.0) {
                second_parameter = 0.0;
                first_parameter = std::clamp(-c / a, 0.0, 1.0);
            } else if (second_parameter > 1.0) {
                second_parameter = 1.0;
                first_parameter = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }
    const Vec3 closest_separation =
        (first_start + first_direction * first_parameter) -
        (second_start + second_direction * second_parameter);
    return dot(closest_separation, closest_separation);
}

[[nodiscard]] double triangle_triangle_distance(const Triangle& first,
                                                const Triangle& second) noexcept {
    double squared = std::numeric_limits<double>::infinity();
    for (const auto& vertex : first.vertices()) {
        squared = std::min(squared, point_triangle_distance_squared(vertex, second));
    }
    for (const auto& vertex : second.vertices()) {
        squared = std::min(squared, point_triangle_distance_squared(vertex, first));
    }
    for (std::size_t first_edge = 0; first_edge < 3U; ++first_edge) {
        for (std::size_t second_edge = 0; second_edge < 3U; ++second_edge) {
            squared = std::min(
                squared,
                segment_segment_distance_squared(
                    first.vertices()[first_edge], first.vertices()[(first_edge + 1U) % 3U],
                    second.vertices()[second_edge], second.vertices()[(second_edge + 1U) % 3U]));
        }
    }
    return std::sqrt(std::max(0.0, squared));
}

[[nodiscard]] bool target_is_valid(std::uint8_t target,
                                   const LinearOctree& tree) noexcept {
    return target >= tree.base_level() && target <= tree.maximum_level();
}

} // 匿名命名空间

OctreeRefinementEngine::OctreeRefinementEngine(
    OctreeRefinementConfiguration configuration, const TriangleBvh* bvh)
    : configuration_(std::move(configuration)), bvh_(bvh) {
    const bool needs_surface = configuration_.surface_target_level.has_value() ||
                               !configuration_.distance_bands.empty() ||
                               configuration_.curvature.has_value() ||
                               configuration_.gap.has_value();
    if (needs_surface && bvh_ == nullptr) {
        throw std::invalid_argument("表面、距离、曲率或狭缝细化必须提供 BVH");
    }
}

void OctreeRefinementEngine::validate_for_tree(const LinearOctree& tree) const {
    const bool needs_surface = configuration_.surface_target_level.has_value() ||
                               !configuration_.distance_bands.empty() ||
                               configuration_.curvature.has_value() ||
                               configuration_.gap.has_value();
    if (needs_surface && bvh_ == nullptr) {
        throw std::invalid_argument("表面、距离、曲率或狭缝细化必须提供 BVH");
    }
    const auto validate_target = [&](std::uint8_t target) {
        if (!target_is_valid(target, tree)) {
            throw std::invalid_argument("细化目标层级必须位于八叉树基础层与最大层之间");
        }
    };
    if (configuration_.surface_target_level) {
        validate_target(*configuration_.surface_target_level);
    }
    for (const auto& band : configuration_.distance_bands) {
        if (!std::isfinite(band.maximum_distance) || band.maximum_distance < 0.0) {
            throw std::invalid_argument("表面距离细化阈值必须是非负有限数");
        }
        validate_target(band.target_level);
    }
    if (configuration_.curvature) {
        const auto& rule = *configuration_.curvature;
        if (!std::isfinite(rule.minimum_normal_angle_degrees) ||
            rule.minimum_normal_angle_degrees < 0.0 ||
            rule.minimum_normal_angle_degrees > 180.0 ||
            !std::isfinite(rule.neighborhood_cell_diagonals) ||
            rule.neighborhood_cell_diagonals <= 0.0) {
            throw std::invalid_argument("曲率细化角度或邻域尺度无效");
        }
        validate_target(rule.target_level);
    }
    if (configuration_.gap) {
        const auto& rule = *configuration_.gap;
        if (!std::isfinite(rule.maximum_search_distance) ||
            rule.maximum_search_distance <= 0.0 ||
            !std::isfinite(rule.maximum_opposing_normal_dot) ||
            rule.maximum_opposing_normal_dot < -1.0 ||
            rule.maximum_opposing_normal_dot > 1.0 ||
            !std::isfinite(rule.minimum_facing_dot) || rule.minimum_facing_dot < 0.0 ||
            rule.minimum_facing_dot > 1.0 || rule.minimum_cells_across_gap == 0) {
            throw std::invalid_argument("狭缝保护参数无效");
        }
    }
    for (const auto& region : configuration_.boxes) {
        validate_target(region.target_level);
    }
    for (const auto& region : configuration_.spheres) {
        if (!is_finite(region.center) || !std::isfinite(region.radius) ||
            region.radius <= 0.0) {
            throw std::invalid_argument("球形细化区必须具有有限圆心和正半径");
        }
        validate_target(region.target_level);
    }
    for (const auto& region : configuration_.cylinders) {
        if (!is_finite(region.first_axis_point) || !is_finite(region.second_axis_point) ||
            !std::isfinite(region.radius) || region.radius <= 0.0 ||
            norm(region.second_axis_point - region.first_axis_point) == 0.0) {
            throw std::invalid_argument("圆柱细化区必须具有非零有限轴线和正半径");
        }
        validate_target(region.target_level);
    }
}

bool OctreeRefinementEngine::curvature_rule_matches(const AABB& bounds) const {
    if (!configuration_.curvature || bvh_ == nullptr) {
        return false;
    }
    const auto& rule = *configuration_.curvature;
    const double radius = rule.neighborhood_cell_diagonals * norm(bounds.extent());
    const auto candidates = bvh_->query(expanded_bounds(bounds, radius));
    const double cosine_threshold =
        std::cos(rule.minimum_normal_angle_degrees * std::numbers::pi / 180.0);
    for (std::size_t first = 0; first < candidates.size(); ++first) {
        const auto& first_triangle =
            bvh_->triangles()[static_cast<std::size_t>(candidates[first])];
        if (first_triangle.area() == 0.0) {
            continue;
        }
        const Vec3 first_normal = unit_normal(first_triangle);
        for (std::size_t second = first + 1; second < candidates.size(); ++second) {
            const auto& second_triangle =
                bvh_->triangles()[static_cast<std::size_t>(candidates[second])];
            if (second_triangle.area() == 0.0) {
                continue;
            }
            const double normal_dot =
                std::clamp(dot(first_normal, unit_normal(second_triangle)), -1.0, 1.0);
            if (normal_dot <= cosine_threshold) {
                return true;
            }
        }
    }
    return false;
}

std::optional<OctreeRefinementEngine::DetectedGap>
OctreeRefinementEngine::detect_gap(const AABB& bounds) const {
    if (!configuration_.gap || bvh_ == nullptr) {
        return std::nullopt;
    }
    const auto& rule = *configuration_.gap;
    const auto candidates = bvh_->query(expanded_bounds(bounds, rule.maximum_search_distance));
    double best = std::numeric_limits<double>::infinity();
    Vec3 best_normal{};
    for (std::size_t first = 0; first < candidates.size(); ++first) {
        const auto& first_triangle =
            bvh_->triangles()[static_cast<std::size_t>(candidates[first])];
        if (first_triangle.area() == 0.0) {
            continue;
        }
        const Vec3 first_normal = unit_normal(first_triangle);
        for (std::size_t second = first + 1; second < candidates.size(); ++second) {
            const auto& second_triangle =
                bvh_->triangles()[static_cast<std::size_t>(candidates[second])];
            if (second_triangle.area() == 0.0) {
                continue;
            }
            const Vec3 second_normal = unit_normal(second_triangle);
            if (dot(first_normal, second_normal) > rule.maximum_opposing_normal_dot) {
                continue;
            }
            const Vec3 separation = second_triangle.centroid() - first_triangle.centroid();
            const double width = triangle_triangle_distance(first_triangle, second_triangle);
            if (width == 0.0 || width > rule.maximum_search_distance) {
                continue;
            }
            const double first_facing = std::abs(dot(separation, first_normal)) / width;
            const double second_facing = std::abs(dot(separation, second_normal)) / width;
            if (first_facing >= rule.minimum_facing_dot &&
                second_facing >= rule.minimum_facing_dot && width < best) {
                const Vec3 opposing_average = first_normal - second_normal;
                const double normal_length = norm(opposing_average);
                if (normal_length == 0.0) {
                    continue;
                }
                best = width;
                best_normal = opposing_average / normal_length;
            }
        }
    }
    if (!std::isfinite(best)) {
        return std::nullopt;
    }
    return DetectedGap{best, best_normal};
}

std::uint32_t OctreeRefinementEngine::required_gap_level(
    const LinearOctree& tree, const DetectedGap& gap) const {
    const auto extent = tree.domain().extent();
    const double root_projection = std::abs(gap.normal.x) * extent.x +
                                   std::abs(gap.normal.y) * extent.y +
                                   std::abs(gap.normal.z) * extent.z;
    const double required_cells =
        root_projection * static_cast<double>(configuration_.gap->minimum_cells_across_gap) /
        gap.width;
    return static_cast<std::uint32_t>(
        std::max(0.0, std::ceil(std::log2(required_cells))));
}

std::uint8_t OctreeRefinementEngine::desired_level(const LinearOctree& tree,
                                                   OctreeNodeCode,
                                                   const AABB& bounds) const {
    std::uint8_t target = tree.base_level();
    if (configuration_.surface_target_level && bvh_->intersects_surface(bounds)) {
        target = std::max(target, *configuration_.surface_target_level);
    }
    if (!configuration_.distance_bands.empty()) {
        const double conservative_distance =
            std::max(0.0, bvh_->distance_to_surface(bounds.center()) - half_diagonal(bounds));
        for (const auto& band : configuration_.distance_bands) {
            if (conservative_distance <= band.maximum_distance) {
                target = std::max(target, band.target_level);
            }
        }
    }
    if (configuration_.curvature && curvature_rule_matches(bounds)) {
        target = std::max(target, configuration_.curvature->target_level);
    }
    if (const auto gap = detect_gap(bounds)) {
        const auto required_level = static_cast<std::uint8_t>(std::clamp(
            required_gap_level(tree, *gap), static_cast<std::uint32_t>(tree.base_level()),
            static_cast<std::uint32_t>(tree.maximum_level())));
        target = std::max(target, required_level);
    }
    for (const auto& region : configuration_.boxes) {
        if (region.bounds.intersects(bounds)) {
            target = std::max(target, region.target_level);
        }
    }
    for (const auto& region : configuration_.spheres) {
        if (sphere_intersects_box(region, bounds)) {
            target = std::max(target, region.target_level);
        }
    }
    for (const auto& region : configuration_.cylinders) {
        if (cylinder_may_intersect_box(region, bounds)) {
            target = std::max(target, region.target_level);
        }
    }
    return target;
}

OctreeAdaptationStatistics OctreeRefinementEngine::apply(LinearOctree& tree) const {
    validate_for_tree(tree);
    OctreeAdaptationStatistics statistics;
    statistics.rule_refinement = tree.refine_to_desired_levels(
        [&](OctreeNodeCode code, const AABB& bounds) {
            return desired_level(tree, code, bounds);
        });
    if (configuration_.enforce_face_2_to_1_balance) {
        statistics.balance = tree.balance_faces_2_to_1();
    } else {
        statistics.balance.initial_leaf_count = tree.leaf_count();
        statistics.balance.final_leaf_count = tree.leaf_count();
    }
    for (const auto code : tree.leaf_codes()) {
        const auto bounds = tree.cell_bounds(code);
        if (configuration_.surface_target_level && bvh_->intersects_surface(bounds)) {
            ++statistics.surface_rule_hits;
        }
        if (!configuration_.distance_bands.empty()) {
            const double conservative_distance =
                std::max(0.0, bvh_->distance_to_surface(bounds.center()) -
                                  half_diagonal(bounds));
            const bool distance_hit = std::any_of(
                configuration_.distance_bands.begin(), configuration_.distance_bands.end(),
                [&](const DistanceRefinementBand& band) {
                    return conservative_distance <= band.maximum_distance;
                });
            statistics.distance_rule_hits += distance_hit ? 1U : 0U;
        }
        if (curvature_rule_matches(bounds)) {
            ++statistics.curvature_rule_hits;
        }
        if (const auto detected_gap = detect_gap(bounds)) {
            ++statistics.gap_rule_hits;
            const auto required_level = required_gap_level(tree, *detected_gap);
            statistics.maximum_required_gap_level =
                std::max(statistics.maximum_required_gap_level, required_level);
            statistics.gap_resolution_failure_count +=
                required_level > tree.maximum_level() ? 1U : 0U;
        }
        bool user_hit = false;
        for (const auto& region : configuration_.boxes) {
            user_hit = user_hit || region.bounds.intersects(bounds);
        }
        for (const auto& region : configuration_.spheres) {
            user_hit = user_hit || sphere_intersects_box(region, bounds);
        }
        for (const auto& region : configuration_.cylinders) {
            user_hit = user_hit || cylinder_may_intersect_box(region, bounds);
        }
        statistics.user_region_rule_hits += user_hit ? 1U : 0U;
    }
    return statistics;
}

} // 命名空间 cartmesh
