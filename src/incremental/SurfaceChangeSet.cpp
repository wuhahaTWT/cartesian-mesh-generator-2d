#include "cartmesh/incremental/SurfaceChangeSet.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace cartmesh {
namespace {

using TriangleKey = std::array<std::uint64_t, 9>;

[[nodiscard]] std::uint64_t coordinate_bits(double value) noexcept {
    // +0 和 -0 表示相同坐标，统一成 +0；SurfaceMesh 已拒绝非有限坐标。
    if (value == 0.0) value = 0.0;
    return std::bit_cast<std::uint64_t>(value);
}

[[nodiscard]] bool vertex_less(const Vec3& first, const Vec3& second) noexcept {
    if (first.x != second.x) return first.x < second.x;
    if (first.y != second.y) return first.y < second.y;
    return first.z < second.z;
}

[[nodiscard]] TriangleKey canonical_key(const Triangle& triangle) {
    auto vertices = triangle.vertices();
    std::sort(vertices.begin(), vertices.end(), vertex_less);
    TriangleKey key{};
    for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex) {
        key[3U * vertex] = coordinate_bits(vertices[vertex].x);
        key[3U * vertex + 1U] = coordinate_bits(vertices[vertex].y);
        key[3U * vertex + 2U] = coordinate_bits(vertices[vertex].z);
    }
    return key;
}

struct TriangleRecord {
    std::uint64_t id{};
    AABB bounds;
};

using TriangleMap = std::map<TriangleKey, std::vector<TriangleRecord>>;

[[nodiscard]] TriangleMap build_map(const SurfaceMesh& surface) {
    TriangleMap result;
    for (std::uint64_t id = 0; id < surface.triangles().size(); ++id) {
        const auto& triangle = surface.triangles()[static_cast<std::size_t>(id)];
        result[canonical_key(triangle)].push_back({id, triangle.bounds()});
    }
    return result;
}

void include_bounds(std::optional<AABB>& accumulated, const AABB& addition) {
    if (!accumulated) {
        accumulated = addition;
        return;
    }
    const Vec3 minimum{
        std::min(accumulated->minimum().x, addition.minimum().x),
        std::min(accumulated->minimum().y, addition.minimum().y),
        std::min(accumulated->minimum().z, addition.minimum().z)};
    const Vec3 maximum{
        std::max(accumulated->maximum().x, addition.maximum().x),
        std::max(accumulated->maximum().y, addition.maximum().y),
        std::max(accumulated->maximum().z, addition.maximum().z)};
    accumulated = AABB(minimum, maximum);
}

} // namespace

SurfaceChangeSet detect_surface_changes(const SurfaceMesh& old_surface,
                                        const SurfaceMesh& new_surface) {
    const auto old_map = build_map(old_surface);
    const auto new_map = build_map(new_surface);
    SurfaceChangeSet result;

    auto old_iterator = old_map.begin();
    auto new_iterator = new_map.begin();
    while (old_iterator != old_map.end() || new_iterator != new_map.end()) {
        if (new_iterator == new_map.end() ||
            (old_iterator != old_map.end() &&
             old_iterator->first < new_iterator->first)) {
            for (const auto& record : old_iterator->second) {
                result.old_triangle_ids.push_back(record.id);
                include_bounds(result.bounds, record.bounds);
            }
            ++old_iterator;
            continue;
        }
        if (old_iterator == old_map.end() ||
            new_iterator->first < old_iterator->first) {
            for (const auto& record : new_iterator->second) {
                result.new_triangle_ids.push_back(record.id);
                include_bounds(result.bounds, record.bounds);
            }
            ++new_iterator;
            continue;
        }

        const auto shared_count =
            std::min(old_iterator->second.size(), new_iterator->second.size());
        for (std::size_t index = shared_count;
             index < old_iterator->second.size(); ++index) {
            const auto& record = old_iterator->second[index];
            result.old_triangle_ids.push_back(record.id);
            include_bounds(result.bounds, record.bounds);
        }
        for (std::size_t index = shared_count;
             index < new_iterator->second.size(); ++index) {
            const auto& record = new_iterator->second[index];
            result.new_triangle_ids.push_back(record.id);
            include_bounds(result.bounds, record.bounds);
        }
        ++old_iterator;
        ++new_iterator;
    }
    std::sort(result.old_triangle_ids.begin(), result.old_triangle_ids.end());
    std::sort(result.new_triangle_ids.begin(), result.new_triangle_ids.end());
    return result;
}

std::optional<AABB> conservative_affected_bounds(
    const SurfaceChangeSet& changes, const AABB& domain, double margin) {
    if (!std::isfinite(margin) || margin < 0.0) {
        throw std::invalid_argument("增量影响范围扩展必须是非负有限数");
    }
    if (!changes.bounds) return std::nullopt;
    const Vec3 padding{margin, margin, margin};
    const Vec3 expanded_minimum = changes.bounds->minimum() - padding;
    const Vec3 expanded_maximum = changes.bounds->maximum() + padding;
    const Vec3 minimum{
        std::max(domain.minimum().x, expanded_minimum.x),
        std::max(domain.minimum().y, expanded_minimum.y),
        std::max(domain.minimum().z, expanded_minimum.z)};
    const Vec3 maximum{
        std::min(domain.maximum().x, expanded_maximum.x),
        std::min(domain.maximum().y, expanded_maximum.y),
        std::min(domain.maximum().z, expanded_maximum.z)};
    if (maximum.x < minimum.x || maximum.y < minimum.y ||
        maximum.z < minimum.z) {
        return std::nullopt;
    }
    return AABB(minimum, maximum);
}

} // namespace cartmesh
