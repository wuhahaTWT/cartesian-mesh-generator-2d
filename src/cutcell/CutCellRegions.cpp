#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace cartmesh {
namespace {

struct Vec2 {
    double x{};
    double y{};
};

struct FaceRegion {
    std::size_t component_id{};
    std::vector<Vec3> polygon;
};

struct PolygonIntersection {
    double area{};
    Vec3 centroid{};
};

class DisjointSets {
  public:
    explicit DisjointSets(std::size_t size) : parent_(size) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }

    [[nodiscard]] std::size_t find(std::size_t value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    void unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first == second) return;
        if (second < first) std::swap(first, second);
        parent_[second] = first;
    }

  private:
    std::vector<std::size_t> parent_;
};

[[nodiscard]] double cross2(const Vec2& first, const Vec2& second) noexcept {
    return first.x * second.y - first.y * second.x;
}

[[nodiscard]] Vec2 project(const Vec3& point, std::size_t axis,
                           const Vec3& origin) noexcept {
    const Vec3 local = point - origin;
    if (axis == 0) return {local.y, local.z};
    if (axis == 1) return {local.x, local.z};
    return {local.x, local.y};
}

[[nodiscard]] Vec3 lift(const Vec2& point, std::size_t axis,
                        const Vec3& origin) noexcept {
    if (axis == 0) return {origin.x, origin.y + point.x, origin.z + point.y};
    if (axis == 1) return {origin.x + point.x, origin.y, origin.z + point.y};
    return {origin.x + point.x, origin.y + point.y, origin.z};
}

[[nodiscard]] std::vector<Vec2> clip_polygon(
    const std::vector<Vec3>& subject, const std::vector<Vec3>& clip,
    std::size_t axis, double signed_area_tolerance) {
    if (subject.size() < 3 || clip.size() < 3) return {};
    const Vec3 origin = subject.front();
    std::vector<Vec2> polygon;
    polygon.reserve(subject.size());
    for (const auto& point : subject) polygon.push_back(project(point, axis, origin));
    std::vector<Vec2> clip_polygon_2d;
    clip_polygon_2d.reserve(clip.size());
    for (const auto& point : clip) clip_polygon_2d.push_back(project(point, axis, origin));
    double clip_twice_area = 0.0;
    for (std::size_t index = 0; index < clip_polygon_2d.size(); ++index) {
        const auto& first = clip_polygon_2d[index];
        const auto& second = clip_polygon_2d[(index + 1U) % clip_polygon_2d.size()];
        clip_twice_area += cross2(first, second);
    }
    const double orientation = clip_twice_area >= 0.0 ? 1.0 : -1.0;
    for (std::size_t edge = 0; edge < clip_polygon_2d.size() && !polygon.empty();
         ++edge) {
        const Vec2 edge_first = clip_polygon_2d[edge];
        const Vec2 edge_second =
            clip_polygon_2d[(edge + 1U) % clip_polygon_2d.size()];
        const Vec2 direction{edge_second.x - edge_first.x,
                             edge_second.y - edge_first.y};
        const auto distance = [&](const Vec2& point) {
            return orientation * cross2(
                direction, {point.x - edge_first.x, point.y - edge_first.y});
        };
        std::vector<Vec2> clipped;
        for (std::size_t index = 0; index < polygon.size(); ++index) {
            const Vec2 first = polygon[index];
            const Vec2 second = polygon[(index + 1U) % polygon.size()];
            const double first_distance = distance(first);
            const double second_distance = distance(second);
            const bool first_inside = first_distance >= -signed_area_tolerance;
            const bool second_inside = second_distance >= -signed_area_tolerance;
            if (first_inside) clipped.push_back(first);
            if (first_inside != second_inside) {
                const double denominator = first_distance - second_distance;
                if (denominator != 0.0) {
                    const double fraction = first_distance / denominator;
                    clipped.push_back(
                        {first.x + (second.x - first.x) * fraction,
                         first.y + (second.y - first.y) * fraction});
                }
            }
        }
        polygon = std::move(clipped);
    }
    return polygon;
}

[[nodiscard]] PolygonIntersection intersect_polygons(
    const std::vector<Vec3>& first, const std::vector<Vec3>& second,
    std::size_t axis, double area_tolerance) {
    const auto polygon = clip_polygon(first, second, axis, area_tolerance);
    if (polygon.size() < 3) return {};
    double twice_area = 0.0;
    Vec2 weighted{};
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& current = polygon[index];
        const auto& next = polygon[(index + 1U) % polygon.size()];
        const double cross_value = cross2(current, next);
        twice_area += cross_value;
        weighted.x += (current.x + next.x) * cross_value;
        weighted.y += (current.y + next.y) * cross_value;
    }
    const double area = 0.5 * std::abs(twice_area);
    if (area <= area_tolerance || twice_area == 0.0) return {};
    const Vec2 centroid{weighted.x / (3.0 * twice_area),
                        weighted.y / (3.0 * twice_area)};
    return {area, lift(centroid, axis, first.front())};
}

[[nodiscard]] std::vector<FaceRegion> face_regions(
    const FluidCellGeometry& cell, std::uint8_t local_face) {
    std::vector<FaceRegion> result;
    if (cell.fluid_component_count == 0) return result;
    if (!cell.fluid_polyhedron_pieces.empty()) {
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            for (const auto& face : piece.polyhedron.faces) {
                if (face.kind != PolyhedronFaceKind::cartesian ||
                    face.source_id != local_face) {
                    continue;
                }
                FaceRegion region;
                region.component_id = piece.component_id;
                for (const auto vertex : face.vertex_indices) {
                    region.polygon.push_back(piece.polyhedron.vertices[vertex]);
                }
                result.push_back(std::move(region));
            }
        }
    } else if (!cell.cartesian_faces[local_face].oriented_boundary_loops.empty()) {
        result.push_back({0, cell.cartesian_faces[local_face].oriented_boundary_loops.front()});
    }
    return result;
}

} // 匿名命名空间

void assign_global_fluid_regions(ConvexCutCellMesh& mesh,
                                 double area_tolerance) {
    double face_area_scale = 0.0;
    std::vector<std::size_t> node_offset(mesh.fluid_cells.size() + 1U, 0);
    for (std::size_t cell = 0; cell < mesh.fluid_cells.size(); ++cell) {
        node_offset[cell + 1U] =
            node_offset[cell] + mesh.fluid_cells[cell].fluid_component_count;
        for (const auto& face : mesh.fluid_cells[cell].cartesian_faces) {
            face_area_scale = std::max(face_area_scale, face.area);
        }
    }
    const double tolerance = std::max(
        area_tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() *
            std::max(face_area_scale, std::numeric_limits<double>::min()));
    DisjointSets components(node_offset.back());
    mesh.component_internal_faces.clear();

    for (const auto& connection : mesh.internal_faces) {
        const auto& first_cell = mesh.fluid_cells[connection.first_fluid_cell_index];
        const auto& second_cell = mesh.fluid_cells[connection.second_fluid_cell_index];
        const auto first_regions = face_regions(first_cell, connection.first_local_face);
        const auto second_regions = face_regions(second_cell, connection.second_local_face);
        using Pair = std::pair<std::size_t, std::size_t>;
        std::map<Pair, std::pair<double, Vec3>> overlap_by_components;
        const std::size_t axis = connection.first_local_face / 2U;
        for (const auto& first_region : first_regions) {
            for (const auto& second_region : second_regions) {
                const auto overlap = intersect_polygons(
                    first_region.polygon, second_region.polygon, axis, tolerance);
                if (overlap.area <= tolerance) continue;
                const Pair key{first_region.component_id,
                               second_region.component_id};
                auto& aggregate = overlap_by_components[key];
                aggregate.second = aggregate.second + overlap.centroid * overlap.area;
                aggregate.first += overlap.area;
            }
        }
        for (const auto& [component_pair, aggregate] : overlap_by_components) {
            if (aggregate.first <= tolerance) continue;
            const auto first_node =
                node_offset[connection.first_fluid_cell_index] + component_pair.first;
            const auto second_node =
                node_offset[connection.second_fluid_cell_index] + component_pair.second;
            components.unite(first_node, second_node);
            mesh.component_internal_faces.push_back(
                {connection.first_background_cell_id,
                 connection.second_background_cell_id,
                 connection.first_fluid_cell_index,
                 connection.second_fluid_cell_index,
                 component_pair.first, component_pair.second, 0,
                 aggregate.first, aggregate.second / aggregate.first,
                 connection.normal});
        }
    }

    std::vector<std::size_t> roots;
    roots.reserve(node_offset.back());
    for (std::size_t node = 0; node < node_offset.back(); ++node) {
        roots.push_back(components.find(node));
    }
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    mesh.global_fluid_region_count = roots.size();
    mesh.global_fluid_region_volumes.assign(roots.size(), 0.0);
    const auto region_id = [&](std::size_t node) {
        const auto root = components.find(node);
        return static_cast<std::uint64_t>(
            std::lower_bound(roots.begin(), roots.end(), root) - roots.begin());
    };

    for (std::size_t cell_index = 0; cell_index < mesh.fluid_cells.size();
         ++cell_index) {
        auto& cell = mesh.fluid_cells[cell_index];
        cell.fluid_component_region_ids.resize(cell.fluid_component_count);
        std::vector<double> local_volume(cell.fluid_component_count, 0.0);
        if (cell.fluid_polyhedron_pieces.empty()) {
            if (!local_volume.empty()) local_volume[0] = cell.volume;
        } else {
            for (auto& piece : cell.fluid_polyhedron_pieces) {
                if (piece.component_id >= local_volume.size()) {
                    throw std::runtime_error(
                        "Cut-cell 局部分量 ID 超出单元分量计数");
                }
                const auto id = region_id(node_offset[cell_index] + piece.component_id);
                piece.global_region_id = id;
                local_volume[piece.component_id] += piece.geometry.volume;
            }
        }
        for (std::size_t component = 0; component < cell.fluid_component_count;
             ++component) {
            const auto id = region_id(node_offset[cell_index] + component);
            cell.fluid_component_region_ids[component] = id;
            mesh.global_fluid_region_volumes[static_cast<std::size_t>(id)] +=
                local_volume[component];
        }
    }
    for (auto& connection : mesh.component_internal_faces) {
        connection.global_region_id =
            mesh.fluid_cells[connection.first_fluid_cell_index]
                .fluid_component_region_ids[connection.first_component_id];
    }
}

} // 命名空间 cartmesh
