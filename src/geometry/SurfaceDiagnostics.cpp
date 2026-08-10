#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/geometry/TriangleTriangleIntersection.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace cartmesh {
namespace {

[[nodiscard]] std::uint64_t coordinate_bits(double value) noexcept {
    if (value == 0.0) {
        value = 0.0;
    }
    return std::bit_cast<std::uint64_t>(value);
}

struct VertexKey {
    std::uint64_t x{};
    std::uint64_t y{};
    std::uint64_t z{};

    [[nodiscard]] bool operator<(const VertexKey& other) const noexcept {
        return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);
    }
};

[[nodiscard]] VertexKey vertex_key(const Vec3& vertex) noexcept {
    return {coordinate_bits(vertex.x), coordinate_bits(vertex.y), coordinate_bits(vertex.z)};
}

struct EdgeKey {
    std::uint64_t first{};
    std::uint64_t second{};

    [[nodiscard]] bool operator<(const EdgeKey& other) const noexcept {
        return std::tie(first, second) < std::tie(other.first, other.second);
    }
};

struct TriangleKey {
    std::array<std::uint64_t, 3> vertices{};

    [[nodiscard]] bool operator<(const TriangleKey& other) const noexcept {
        return vertices < other.vertices;
    }
};

struct EdgeUse {
    std::uint64_t triangle{};
    bool low_to_high{};
};

class DisjointSet {
  public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }

    [[nodiscard]] std::size_t find(std::size_t value) {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }

    void unite(std::size_t lhs, std::size_t rhs) {
        lhs = find(lhs);
        rhs = find(rhs);
        if (lhs == rhs) {
            return;
        }
        if (rank_[lhs] < rank_[rhs] || (rank_[lhs] == rank_[rhs] && lhs > rhs)) {
            std::swap(lhs, rhs);
        }
        parent_[rhs] = lhs;
        if (rank_[lhs] == rank_[rhs]) {
            ++rank_[lhs];
        }
    }

  private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
};

[[nodiscard]] SurfaceDiagnosticLocation location_for(
    const EdgeKey& edge, const std::vector<EdgeUse>& uses, const std::vector<Vec3>& vertices) {
    return {uses.front().triangle,
            uses.size() > 1 ? uses[1].triangle : uses.front().triangle,
            vertices[static_cast<std::size_t>(edge.first)],
            vertices[static_cast<std::size_t>(edge.second)]};
}

[[nodiscard]] long double signed_six_volume_relative(const Triangle& triangle,
                                                     const Vec3& reference) noexcept {
    const auto& points = triangle.vertices();
    const long double ax = static_cast<long double>(points[0].x) - reference.x;
    const long double ay = static_cast<long double>(points[0].y) - reference.y;
    const long double az = static_cast<long double>(points[0].z) - reference.z;
    const long double bx = static_cast<long double>(points[1].x) - reference.x;
    const long double by = static_cast<long double>(points[1].y) - reference.y;
    const long double bz = static_cast<long double>(points[1].z) - reference.z;
    const long double cx = static_cast<long double>(points[2].x) - reference.x;
    const long double cy = static_cast<long double>(points[2].y) - reference.y;
    const long double cz = static_cast<long double>(points[2].z) - reference.z;
    return ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
           az * (bx * cy - by * cx);
}

struct ComponentWork {
    std::vector<std::size_t> triangle_ids;
    Vec3 reference{};
    Vec3 minimum{};
    Vec3 maximum{};
    double surface_area{};
};

[[nodiscard]] bool point_in_bounds(const Vec3& point, const ComponentWork& component) noexcept {
    return point.x >= component.minimum.x && point.x <= component.maximum.x &&
           point.y >= component.minimum.y && point.y <= component.maximum.y &&
           point.z >= component.minimum.z && point.z <= component.maximum.z;
}

[[nodiscard]] std::optional<bool> component_contains_point(
    const SurfaceMesh& mesh, const ComponentWork& component, const Vec3& point) noexcept {
    if (!point_in_bounds(point, component)) {
        return false;
    }
    long double solid_angle = 0.0L;
    for (const auto triangle_id : component.triangle_ids) {
        const auto& vertices = mesh.triangles()[triangle_id].vertices();
        const long double ax = static_cast<long double>(vertices[0].x) - point.x;
        const long double ay = static_cast<long double>(vertices[0].y) - point.y;
        const long double az = static_cast<long double>(vertices[0].z) - point.z;
        const long double bx = static_cast<long double>(vertices[1].x) - point.x;
        const long double by = static_cast<long double>(vertices[1].y) - point.y;
        const long double bz = static_cast<long double>(vertices[1].z) - point.z;
        const long double cx = static_cast<long double>(vertices[2].x) - point.x;
        const long double cy = static_cast<long double>(vertices[2].y) - point.y;
        const long double cz = static_cast<long double>(vertices[2].z) - point.z;
        const long double a_length = std::sqrt(ax * ax + ay * ay + az * az);
        const long double b_length = std::sqrt(bx * bx + by * by + bz * bz);
        const long double c_length = std::sqrt(cx * cx + cy * cy + cz * cz);
        if (a_length == 0.0L || b_length == 0.0L || c_length == 0.0L) {
            return std::nullopt;
        }
        const long double numerator =
            ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
            az * (bx * cy - by * cx);
        const long double denominator =
            a_length * b_length * c_length +
            (ax * bx + ay * by + az * bz) * c_length +
            (bx * cx + by * cy + bz * cz) * a_length +
            (cx * ax + cy * ay + cz * az) * b_length;
        if (numerator == 0.0L && denominator == 0.0L) {
            return std::nullopt;
        }
        solid_angle += 2.0L * std::atan2(numerator, denominator);
    }
    return std::abs(solid_angle) > 2.0L * std::numbers::pi_v<long double>;
}

[[nodiscard]] AABB expanded_bounds(const AABB& bounds, double amount) {
    const Vec3 delta{amount, amount, amount};
    return AABB(bounds.minimum() - delta, bounds.maximum() + delta);
}

[[nodiscard]] SurfaceDiagnosticTrianglePair triangle_pair_location(
    const SurfaceMesh& mesh, std::uint64_t first,
    std::uint64_t second) noexcept {
    return {first, second,
            mesh.triangles()[static_cast<std::size_t>(first)].centroid(),
            mesh.triangles()[static_cast<std::size_t>(second)].centroid()};
}

} // 匿名命名空间

SurfaceDiagnostics diagnose_surface(const SurfaceMesh& mesh, std::size_t max_examples) {
    SurfaceDiagnostics result;
    result.triangle_count = mesh.triangles().size();
    const double diagonal = norm(mesh.bounds().extent());
    const double scale = std::max(diagonal, std::numeric_limits<double>::min());
    result.suggested_length_tolerance = scale * 1.0e-12;
    result.degenerate_area_tolerance =
        0.5 * scale * result.suggested_length_tolerance;

    std::map<VertexKey, std::uint64_t> vertex_ids;
    std::vector<Vec3> vertices;
    std::vector<std::array<std::uint64_t, 3>> triangle_vertices;
    triangle_vertices.reserve(mesh.triangles().size());
    for (const auto& triangle : mesh.triangles()) {
        std::array<std::uint64_t, 3> ids{};
        for (std::size_t local = 0; local < 3; ++local) {
            const auto& vertex = triangle.vertices()[local];
            const auto [iterator, inserted] = vertex_ids.try_emplace(vertex_key(vertex), 0);
            if (inserted) {
                iterator->second = vertices.size();
                vertices.push_back(vertex);
            }
            ids[local] = iterator->second;
        }
        triangle_vertices.push_back(ids);
    }
    result.unique_vertex_count = vertices.size();

    std::map<TriangleKey, std::uint64_t> first_triangle_by_key;
    std::map<EdgeKey, std::vector<EdgeUse>> edge_uses;
    std::vector<std::vector<std::uint64_t>> incident_triangles(vertices.size());
    DisjointSet components(mesh.triangles().size());
    std::vector<bool> valid_triangle(mesh.triangles().size(), true);
    for (std::size_t triangle_index = 0; triangle_index < mesh.triangles().size();
         ++triangle_index) {
        const auto& triangle = mesh.triangles()[triangle_index];
        const auto& points = triangle.vertices();
        const double longest_edge =
            std::max({norm(points[1] - points[0]), norm(points[2] - points[1]),
                      norm(points[0] - points[2])});
        const double triangle_area_tolerance =
            0.5 * longest_edge * result.suggested_length_tolerance;
        if (triangle.area() <= triangle_area_tolerance) {
            valid_triangle[triangle_index] = false;
            ++result.degenerate_triangle_count;
            if (result.degenerate_triangle_examples.size() < max_examples) {
                result.degenerate_triangle_examples.push_back(triangle_index);
            }
            continue;
        }
        const auto& original_ids = triangle_vertices[triangle_index];
        for (const auto vertex_id : original_ids) {
            incident_triangles[static_cast<std::size_t>(vertex_id)].push_back(triangle_index);
        }
        auto sorted_ids = original_ids;
        std::sort(sorted_ids.begin(), sorted_ids.end());
        const TriangleKey triangle_key_value{sorted_ids};
        const auto [duplicate_iterator, inserted] =
            first_triangle_by_key.try_emplace(triangle_key_value, triangle_index);
        if (!inserted) {
            ++result.duplicate_triangle_count;
            if (result.duplicate_triangle_examples.size() < max_examples) {
                result.duplicate_triangle_examples.push_back(triangle_index);
            }
        }
        static_cast<void>(duplicate_iterator);

        for (std::size_t edge = 0; edge < 3; ++edge) {
            const auto from = original_ids[edge];
            const auto to = original_ids[(edge + 1) % 3];
            const EdgeKey key{std::min(from, to), std::max(from, to)};
            edge_uses[key].push_back(
                {static_cast<std::uint64_t>(triangle_index), from == key.first});
        }
    }
    result.unique_edge_count = edge_uses.size();
    std::vector<std::vector<EdgeKey>> incident_edges(vertices.size());

    for (const auto& [edge, uses] : edge_uses) {
        incident_edges[static_cast<std::size_t>(edge.first)].push_back(edge);
        incident_edges[static_cast<std::size_t>(edge.second)].push_back(edge);
        for (std::size_t use = 1; use < uses.size(); ++use) {
            components.unite(static_cast<std::size_t>(uses.front().triangle),
                             static_cast<std::size_t>(uses[use].triangle));
        }
        if (uses.size() == 1) {
            ++result.boundary_edge_count;
            if (result.boundary_edge_examples.size() < max_examples) {
                result.boundary_edge_examples.push_back(location_for(edge, uses, vertices));
            }
        } else if (uses.size() > 2) {
            ++result.non_manifold_edge_count;
            if (result.non_manifold_edge_examples.size() < max_examples) {
                result.non_manifold_edge_examples.push_back(location_for(edge, uses, vertices));
            }
        } else if (uses[0].low_to_high == uses[1].low_to_high) {
            ++result.orientation_conflict_edge_count;
            if (result.orientation_conflict_examples.size() < max_examples) {
                result.orientation_conflict_examples.push_back(location_for(edge, uses, vertices));
            }
        }
    }

    for (std::size_t vertex_id = 0; vertex_id < vertices.size(); ++vertex_id) {
        const auto& triangles = incident_triangles[vertex_id];
        if (triangles.empty()) {
            continue;
        }
        std::map<std::uint64_t, std::size_t> local_index;
        for (std::size_t index = 0; index < triangles.size(); ++index) {
            local_index.emplace(triangles[index], index);
        }
        std::vector<std::vector<std::size_t>> link(triangles.size());
        std::size_t boundary_edges = 0;
        bool invalid_edge_incidence = false;
        for (const auto& edge : incident_edges[vertex_id]) {
            const auto& uses = edge_uses.at(edge);
            if (uses.size() == 1) {
                ++boundary_edges;
            } else if (uses.size() == 2) {
                const auto first = local_index.at(uses[0].triangle);
                const auto second = local_index.at(uses[1].triangle);
                if (first != second) {
                    link[first].push_back(second);
                    link[second].push_back(first);
                }
            } else {
                invalid_edge_incidence = true;
            }
        }

        std::vector<bool> visited(link.size(), false);
        std::queue<std::size_t> pending;
        pending.push(0);
        visited[0] = true;
        std::size_t visited_count = 0;
        while (!pending.empty()) {
            const auto current = pending.front();
            pending.pop();
            ++visited_count;
            for (const auto neighbour : link[current]) {
                if (!visited[neighbour]) {
                    visited[neighbour] = true;
                    pending.push(neighbour);
                }
            }
        }
        const auto degree_one = static_cast<std::size_t>(
            std::count_if(link.begin(), link.end(),
                          [](const auto& neighbours) { return neighbours.size() == 1; }));
        const bool degree_valid = std::all_of(
            link.begin(), link.end(),
            [](const auto& neighbours) { return neighbours.size() <= 2; });
        const bool closed_link =
            boundary_edges == 0 && degree_valid &&
            std::all_of(link.begin(), link.end(),
                        [](const auto& neighbours) { return neighbours.size() == 2; });
        const bool single_triangle_boundary_link =
            triangles.size() == 1 && boundary_edges == 2 && link.front().empty();
        const bool open_link =
            boundary_edges == 2 && degree_valid && degree_one == 2;
        const bool manifold_vertex =
            !invalid_edge_incidence && visited_count == link.size() &&
            (closed_link || open_link || single_triangle_boundary_link);
        if (!manifold_vertex) {
            ++result.non_manifold_vertex_count;
            if (result.non_manifold_vertex_examples.size() < max_examples) {
                result.non_manifold_vertex_examples.push_back(
                    {vertices[vertex_id], static_cast<std::uint64_t>(triangles.size())});
            }
        }
    }

    std::vector<std::size_t> roots;
    roots.reserve(mesh.triangles().size());
    for (std::size_t triangle_index = 0; triangle_index < valid_triangle.size();
         ++triangle_index) {
        if (valid_triangle[triangle_index]) {
            roots.push_back(components.find(triangle_index));
        }
    }
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    result.connected_component_count = roots.size();
    result.manifold =
        result.non_manifold_edge_count == 0 && result.non_manifold_vertex_count == 0;
    result.closed = result.boundary_edge_count == 0 && result.manifold;
    result.consistently_oriented = result.orientation_conflict_edge_count == 0;

    const TriangleBvh intersection_bvh(mesh);
    for (std::size_t first = 0; first < mesh.triangles().size(); ++first) {
        if (!valid_triangle[first]) continue;
        const auto candidates = intersection_bvh.query(expanded_bounds(
            mesh.triangles()[first].bounds(), result.suggested_length_tolerance));
        auto first_vertices = triangle_vertices[first];
        std::sort(first_vertices.begin(), first_vertices.end());
        for (const auto second_id : candidates) {
            if (second_id <= first || !valid_triangle[static_cast<std::size_t>(second_id)]) {
                continue;
            }
            auto second_vertices =
                triangle_vertices[static_cast<std::size_t>(second_id)];
            std::sort(second_vertices.begin(), second_vertices.end());
            if (first_vertices == second_vertices) continue;
            std::size_t shared_vertex_count = 0;
            for (const auto first_vertex : first_vertices) {
                shared_vertex_count += static_cast<std::size_t>(
                    std::find(second_vertices.begin(), second_vertices.end(),
                              first_vertex) != second_vertices.end());
            }
            const auto relation = classify_triangle_triangle(
                mesh.triangles()[first],
                mesh.triangles()[static_cast<std::size_t>(second_id)],
                result.suggested_length_tolerance);
            if (relation == TriangleTriangleRelation::coplanar_area_overlap) {
                ++result.overlapping_triangle_pair_count;
                if (result.overlapping_triangle_examples.size() < max_examples) {
                    result.overlapping_triangle_examples.push_back(
                        triangle_pair_location(mesh, first, second_id));
                }
            } else if (relation ==
                           TriangleTriangleRelation::proper_intersection &&
                       shared_vertex_count < 2) {
                ++result.self_intersection_pair_count;
                if (result.self_intersection_examples.size() < max_examples) {
                    result.self_intersection_examples.push_back(
                        triangle_pair_location(mesh, first, second_id));
                }
            } else if (relation == TriangleTriangleRelation::boundary_contact &&
                       shared_vertex_count == 0) {
                ++result.non_adjacent_contact_pair_count;
                if (result.non_adjacent_contact_examples.size() < max_examples) {
                    result.non_adjacent_contact_examples.push_back(
                        triangle_pair_location(mesh, first, second_id));
                }
            }
        }
    }

    std::map<std::size_t, std::size_t> component_index_by_root;
    std::vector<ComponentWork> component_work;
    component_work.reserve(roots.size());
    for (std::size_t component_index = 0; component_index < roots.size(); ++component_index) {
        component_index_by_root.emplace(roots[component_index], component_index);
        component_work.emplace_back();
    }
    for (std::size_t triangle_index = 0; triangle_index < valid_triangle.size();
         ++triangle_index) {
        if (!valid_triangle[triangle_index]) {
            continue;
        }
        const auto root = components.find(triangle_index);
        const auto component_index = component_index_by_root.at(root);
        auto& work = component_work[component_index];
        const auto& triangle = mesh.triangles()[triangle_index];
        if (work.triangle_ids.empty()) {
            work.reference = triangle.vertices()[0];
            work.minimum = triangle.vertices()[0];
            work.maximum = triangle.vertices()[0];
        }
        work.triangle_ids.push_back(triangle_index);
        work.surface_area += triangle.area();
        for (const auto& point : triangle.vertices()) {
            work.minimum.x = std::min(work.minimum.x, point.x);
            work.minimum.y = std::min(work.minimum.y, point.y);
            work.minimum.z = std::min(work.minimum.z, point.z);
            work.maximum.x = std::max(work.maximum.x, point.x);
            work.maximum.y = std::max(work.maximum.y, point.y);
            work.maximum.z = std::max(work.maximum.z, point.z);
        }
    }

    result.components.reserve(component_work.size());
    long double total_signed_volume = 0.0L;
    for (std::size_t component_index = 0; component_index < component_work.size();
         ++component_index) {
        const auto& work = component_work[component_index];
        long double signed_six_volume = 0.0L;
        for (const auto triangle_id : work.triangle_ids) {
            signed_six_volume +=
                signed_six_volume_relative(mesh.triangles()[triangle_id], work.reference);
        }
        const long double component_volume = signed_six_volume / 6.0L;
        total_signed_volume += component_volume;
        result.components.push_back(
            {static_cast<std::uint64_t>(component_index),
             static_cast<std::uint64_t>(work.triangle_ids.size()),
             work.surface_area, AABB(work.minimum, work.maximum),
             norm(work.maximum - work.minimum),
             static_cast<double>(component_volume), 0, 1, false, false,
             mesh.triangles()[work.triangle_ids.front()].centroid()});
    }
    result.signed_volume = static_cast<double>(total_signed_volume);

    const bool can_check_component_orientation =
        result.closed && result.consistently_oriented && result.degenerate_triangle_count == 0 &&
        result.duplicate_triangle_count == 0;
    if (can_check_component_orientation) {
        for (std::size_t component_index = 0; component_index < component_work.size();
             ++component_index) {
            auto& component = result.components[component_index];
            bool containment_ambiguous = false;
            std::uint64_t nesting_depth = 0;
            for (std::size_t other = 0; other < component_work.size(); ++other) {
                if (other == component_index) {
                    continue;
                }
                const auto contains = component_contains_point(
                    mesh, component_work[other], component.sample_position);
                if (!contains) {
                    containment_ambiguous = true;
                    break;
                }
                nesting_depth += *contains ? 1U : 0U;
            }
            component.nesting_depth = nesting_depth;
            component.expected_orientation_sign = (nesting_depth & 1U) == 0U ? 1 : -1;
            component.orientation_checked = !containment_ambiguous;
            component.orientation_matches_nesting =
                component.orientation_checked &&
                (component.expected_orientation_sign > 0 ? component.signed_volume > 0.0
                                                         : component.signed_volume < 0.0);
            if (!component.orientation_matches_nesting) {
                ++result.component_orientation_mismatch_count;
            }
        }
    }
    long double material_volume = 0.0L;
    for (const auto& component : result.components) {
        material_volume += static_cast<long double>(component.expected_orientation_sign) *
                           std::abs(static_cast<long double>(component.signed_volume));
    }
    result.material_volume = static_cast<double>(material_volume);
    result.small_component_diagonal_threshold = scale * 1.0e-6;
    result.minimum_component_bounding_diagonal =
        std::numeric_limits<double>::infinity();
    result.minimum_component_absolute_volume =
        std::numeric_limits<double>::infinity();
    for (const auto& component : result.components) {
        result.minimum_component_bounding_diagonal =
            std::min(result.minimum_component_bounding_diagonal,
                     component.bounding_diagonal);
        result.minimum_component_absolute_volume =
            std::min(result.minimum_component_absolute_volume,
                     std::abs(component.signed_volume));
        if (component.bounding_diagonal <=
            result.small_component_diagonal_threshold) {
            ++result.small_component_count;
        }
    }
    if (result.components.empty()) {
        result.minimum_component_bounding_diagonal = 0.0;
        result.minimum_component_absolute_volume = 0.0;
    }
    return result;
}

} // 命名空间 cartmesh
