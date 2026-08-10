#include "cartmesh/classify/SurfaceClassifier.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

[[nodiscard]] double coordinate_ulp(double value) noexcept {
    const double magnitude = std::abs(value);
    return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
}

[[nodiscard]] double bounds_coordinate_resolution(const AABB& bounds) noexcept {
    const auto& minimum = bounds.minimum();
    const auto& maximum = bounds.maximum();
    return std::max({coordinate_ulp(minimum.x), coordinate_ulp(minimum.y),
                     coordinate_ulp(minimum.z), coordinate_ulp(maximum.x),
                     coordinate_ulp(maximum.y), coordinate_ulp(maximum.z)});
}

[[nodiscard]] Vec3 normalized(Vec3 value) {
    const double length = norm(value);
    if (length == 0.0 || !std::isfinite(length)) {
        throw std::invalid_argument("分类射线方向必须有限且非零");
    }
    return value / length;
}

} // 匿名命名空间

SurfaceClassifier::SurfaceClassifier(const TriangleBvh& bvh, double surface_tolerance)
    : bvh_(bvh),
      surface_tolerance_(surface_tolerance > 0.0
                             ? surface_tolerance
                             : std::max(norm(bvh.bounds().extent()) * 1.0e-12,
                                        4.0 * bounds_coordinate_resolution(bvh.bounds()))),
      directions_{normalized({1.0, 0.372013, 0.139021}),
                  normalized({-0.231017, 1.0, 0.417019}),
                  normalized({0.311023, -0.271013, 1.0})} {
    if (!std::isfinite(surface_tolerance) || surface_tolerance < 0.0) {
        throw std::invalid_argument("分类表面容差必须是非负有限数");
    }
}

PointClassificationResult SurfaceClassifier::classify(const Vec3& point) const {
    if (bvh_.point_on_surface(point, surface_tolerance_)) {
        return {PointClassification::on_surface, 0, 0, 0};
    }
    PointClassificationResult result;
    for (const auto& direction : directions_) {
        const auto intersections = bvh_.trace_ray(point, direction);
        if (intersections.ambiguous) {
            ++result.ambiguous_rays;
            continue;
        }
        if ((intersections.unique_surface_crossings & 1U) != 0U) {
            ++result.inside_votes;
        } else {
            ++result.outside_votes;
        }
    }
    if (result.inside_votes > 0 && result.outside_votes == 0) {
        result.classification = PointClassification::inside;
    } else if (result.outside_votes > 0 && result.inside_votes == 0) {
        result.classification = PointClassification::outside;
    } else {
        result.classification = PointClassification::conflict;
    }
    return result;
}

CellClassification SurfaceClassifier::classify_cell(const AABB& bounds) const {
    return classify_cell(bounds, classify(bounds.center()).classification);
}

CellClassification SurfaceClassifier::classify_cell(
    const AABB& bounds, PointClassification center_classification) const {
    if (bvh_.intersects_surface(bounds)) {
        return CellClassification::intersected;
    }
    if (center_classification == PointClassification::inside) {
        return CellClassification::inside;
    }
    if (center_classification == PointClassification::outside) {
        return CellClassification::outside;
    }
    return CellClassification::conflict;
}

UniformClassification classify_uniform_cells(const UniformCartesianGrid& grid,
                                               const SurfaceClassifier& classifier) {
    UniformClassification result;
    result.cell_classification.resize(static_cast<std::size_t>(grid.cell_count()));
    result.center_point_classification.resize(static_cast<std::size_t>(grid.cell_count()));
    for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
        const auto bounds = grid.cell_bounds(grid.cell_key(id));
        const auto center_classification = classifier.classify(bounds.center()).classification;
        const auto cell_classification = classifier.classify_cell(bounds, center_classification);
        result.cell_classification[static_cast<std::size_t>(id)] =
            static_cast<std::uint8_t>(cell_classification);
        result.center_point_classification[static_cast<std::size_t>(id)] =
            static_cast<std::uint8_t>(center_classification);
        switch (cell_classification) {
        case CellClassification::inside:
            ++result.inside_count;
            break;
        case CellClassification::outside:
            ++result.outside_count;
            break;
        case CellClassification::intersected:
            ++result.intersected_count;
            break;
        case CellClassification::conflict:
            ++result.conflict_count;
            break;
        }
        switch (center_classification) {
        case PointClassification::inside:
            ++result.center_inside_count;
            break;
        case PointClassification::outside:
            ++result.center_outside_count;
            break;
        case PointClassification::on_surface:
            ++result.center_on_surface_count;
            break;
        case PointClassification::conflict:
            ++result.center_conflict_count;
            break;
        }
    }
    return result;
}

AdaptiveClassification classify_octree_leaves(const LinearOctree& tree,
                                               const SurfaceClassifier& classifier) {
    AdaptiveClassification result;
    result.cell_classification.resize(static_cast<std::size_t>(tree.leaf_count()));
    result.center_point_classification.resize(static_cast<std::size_t>(tree.leaf_count()));
    for (std::uint64_t leaf_id = 0; leaf_id < tree.leaf_count(); ++leaf_id) {
        const auto code = tree.leaf_code(leaf_id);
        const auto bounds = tree.cell_bounds(code);
        const auto center_classification = classifier.classify(bounds.center()).classification;
        const auto cell_classification = classifier.classify_cell(bounds, center_classification);
        result.cell_classification[static_cast<std::size_t>(leaf_id)] =
            static_cast<std::uint8_t>(cell_classification);
        result.center_point_classification[static_cast<std::size_t>(leaf_id)] =
            static_cast<std::uint8_t>(center_classification);
        switch (cell_classification) {
        case CellClassification::inside:
            ++result.inside_count;
            result.inside_volume += bounds.volume();
            result.inside_plus_intersected_volume += bounds.volume();
            break;
        case CellClassification::outside:
            ++result.outside_count;
            break;
        case CellClassification::intersected:
            ++result.intersected_count;
            result.inside_plus_intersected_volume += bounds.volume();
            break;
        case CellClassification::conflict:
            ++result.conflict_count;
            break;
        }
        switch (center_classification) {
        case PointClassification::inside:
            ++result.center_inside_count;
            break;
        case PointClassification::outside:
            ++result.center_outside_count;
            break;
        case PointClassification::on_surface:
            ++result.center_on_surface_count;
            break;
        case PointClassification::conflict:
            ++result.center_conflict_count;
            break;
        }
    }
    return result;
}

} // 命名空间 cartmesh
