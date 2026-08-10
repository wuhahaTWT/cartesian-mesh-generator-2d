#include "cartmesh/cutcell/ConvexSurfaceCutter.hpp"

#include "cartmesh/geometry/SurfaceDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cartmesh {
namespace {

struct CandidatePlane {
    Vec3 point{};
    Vec3 outward_normal{};
};

class DisjointSets {
  public:
    explicit DisjointSets(std::size_t size) : parent_(size) {
        for (std::size_t index = 0; index < size; ++index) {
            parent_[index] = index;
        }
    }

    [[nodiscard]] std::size_t find(std::size_t value) {
        while (parent_[value] != value) {
            parent_[value] = parent_[parent_[value]];
            value = parent_[value];
        }
        return value;
    }

    void unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first != second) {
            parent_[second] = first;
        }
    }

  private:
    std::vector<std::size_t> parent_;
};

[[nodiscard]] bool plane_less(const CandidatePlane& first,
                              const CandidatePlane& second) noexcept {
    if (first.outward_normal.x != second.outward_normal.x) {
        return first.outward_normal.x < second.outward_normal.x;
    }
    if (first.outward_normal.y != second.outward_normal.y) {
        return first.outward_normal.y < second.outward_normal.y;
    }
    if (first.outward_normal.z != second.outward_normal.z) {
        return first.outward_normal.z < second.outward_normal.z;
    }
    const double first_offset = dot(first.outward_normal, first.point);
    const double second_offset = dot(second.outward_normal, second.point);
    if (first_offset != second_offset) {
        return first_offset < second_offset;
    }
    if (first.point.x != second.point.x) {
        return first.point.x < second.point.x;
    }
    if (first.point.y != second.point.y) {
        return first.point.y < second.point.y;
    }
    return first.point.z < second.point.z;
}

[[nodiscard]] bool same_plane(const CandidatePlane& first,
                              const CandidatePlane& second,
                              double length_tolerance) noexcept {
    const double angular_tolerance =
        128.0 * std::numeric_limits<double>::epsilon();
    return norm(first.outward_normal - second.outward_normal) <= angular_tolerance &&
           std::abs(dot(first.outward_normal, second.point - first.point)) <=
               length_tolerance;
}

} // 匿名命名空间

ConvexSurfaceCutter::ConvexSurfaceCutter(const SurfaceMesh& surface,
                                         std::uint64_t boundary_id,
                                         double length_tolerance) {
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("凸表面 Cut-cell 长度容差必须为非负有限数");
    }
    const auto diagnostics = diagnose_surface(surface);
    if (!diagnostics.valid_for_stage1_classification()) {
        throw std::invalid_argument("凸表面 Cut-cell 需要封闭、定向一致且流形的 STL");
    }
    if (diagnostics.connected_component_count != 1) {
        throw std::invalid_argument("当前凸表面 Cut-cell 内核只处理单连通分量");
    }
    length_tolerance_ = std::max(length_tolerance, diagnostics.suggested_length_tolerance);
    boundary_id_ = boundary_id;
    input_orientation_reversed_ = diagnostics.signed_volume < 0.0;
    const double orientation_sign = input_orientation_reversed_ ? -1.0 : 1.0;

    std::vector<CandidatePlane> candidates;
    candidates.reserve(surface.triangles().size());
    for (const auto& triangle : surface.triangles()) {
        const Vec3 area_vector = triangle.area_vector() * orientation_sign;
        const double area = norm(area_vector);
        if (!(area > 0.0)) {
            throw std::invalid_argument("凸表面 Cut-cell 不接受零面积三角形");
        }
        candidates.push_back({triangle.vertices().front(), area_vector / area});
    }

    for (const auto& plane : candidates) {
        for (const auto& triangle : surface.triangles()) {
            for (const auto& vertex : triangle.vertices()) {
                if (dot(plane.outward_normal, vertex - plane.point) >
                    length_tolerance_) {
                    throw std::invalid_argument(
                        "STL 不是凸集边界，不得用凸半空间 Cut-cell 内核重建");
                }
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), plane_less);
    std::vector<CandidatePlane> unique;
    for (const auto& candidate : candidates) {
        if (unique.empty() || !same_plane(unique.back(), candidate, length_tolerance_)) {
            unique.push_back(candidate);
        }
    }
    planes_.reserve(unique.size());
    for (const auto& plane : unique) {
        planes_.emplace_back(plane.point, plane.outward_normal, boundary_id);
    }
}

ConvexSurfaceCutResult ConvexSurfaceCutter::cut_box(const AABB& box) const {
    ConvexSurfaceCutResult result;
    result.solid_polyhedron = make_box_polyhedron(box);
    for (const auto& plane : planes_) {
        result.solid_polyhedron = clip_convex_polyhedron(
            result.solid_polyhedron, plane, length_tolerance_);
        if (result.solid_polyhedron.empty()) {
            break;
        }
    }
    result.solid_geometry = measure_polyhedron(result.solid_polyhedron);
    const double cell_volume = box.volume();
    const double fraction_tolerance =
        256.0 * std::numeric_limits<double>::epsilon();
    double solid_volume = result.solid_geometry.volume;
    if (solid_volume / cell_volume <= fraction_tolerance) {
        solid_volume = 0.0;
        result.solid_geometry.volume = 0.0;
    } else if ((cell_volume - solid_volume) / cell_volume <= fraction_tolerance) {
        solid_volume = cell_volume;
        result.solid_geometry.volume = cell_volume;
        result.solid_geometry.centroid = box.center();
    }
    result.solid_volume_fraction = solid_volume / cell_volume;

    ConvexPolyhedron remaining = make_box_polyhedron(box);
    for (std::size_t plane_index = 0; plane_index < planes_.size(); ++plane_index) {
        const auto& plane = planes_[plane_index];
        const OrientedHalfSpace outside_half_space(
            plane.point(), plane.outward_normal() * -1.0,
            static_cast<std::uint64_t>(plane_index));
        auto outside_piece = clip_convex_polyhedron(
            remaining, outside_half_space, length_tolerance_,
            PolyhedronFaceKind::internal_partition);
        auto outside_geometry = measure_polyhedron(outside_piece);
        if (outside_geometry.positive_volume &&
            outside_geometry.volume / cell_volume > fraction_tolerance) {
            result.fluid_decomposition_polyhedra.push_back(
                std::move(outside_piece));
            result.fluid_decomposition_geometries.push_back(
                std::move(outside_geometry));
        }
        const OrientedHalfSpace inside_half_space(
            plane.point(), plane.outward_normal(),
            static_cast<std::uint64_t>(plane_index));
        remaining = clip_convex_polyhedron(
            remaining, inside_half_space, length_tolerance_,
            PolyhedronFaceKind::internal_partition);
        if (remaining.empty()) {
            break;
        }
    }

    const Vec3 cell_centroid = box.center();
    Vec3 relative_fluid_first_moment{};
    for (const auto& geometry : result.fluid_decomposition_geometries) {
        result.fluid_volume += geometry.volume;
        relative_fluid_first_moment =
            relative_fluid_first_moment +
            (geometry.centroid - cell_centroid) * geometry.volume;
    }
    if (result.fluid_volume / cell_volume <= fraction_tolerance) {
        result.fluid_volume = 0.0;
        result.fluid_decomposition_polyhedra.clear();
        result.fluid_decomposition_geometries.clear();
    } else if ((cell_volume - result.fluid_volume) / cell_volume <=
               fraction_tolerance) {
        result.fluid_volume = cell_volume;
    }
    result.fluid_volume_fraction = result.fluid_volume / cell_volume;
    result.volume_conservation_residual =
        solid_volume + result.fluid_volume - cell_volume;
    if (result.fluid_volume > 0.0) {
        result.fluid_centroid =
            cell_centroid + relative_fluid_first_moment / result.fluid_volume;
    }
    for (std::size_t face = 0; face < result.solid_polyhedron.faces.size(); ++face) {
        if (result.solid_polyhedron.faces[face].kind ==
            PolyhedronFaceKind::embedded_boundary) {
            result.embedded_boundary_area += result.solid_geometry.faces[face].area;
        }
    }
    result.cut = result.solid_volume_fraction > fraction_tolerance &&
                 result.fluid_volume_fraction > fraction_tolerance &&
                 result.embedded_boundary_area > 0.0;

    DisjointSets components(result.fluid_decomposition_polyhedra.size());
    for (std::size_t first = 0;
         first < result.fluid_decomposition_polyhedra.size(); ++first) {
        const auto& first_polyhedron = result.fluid_decomposition_polyhedra[first];
        const auto& first_geometry = result.fluid_decomposition_geometries[first];
        for (std::size_t second = first + 1;
             second < result.fluid_decomposition_polyhedra.size(); ++second) {
            const auto& second_polyhedron =
                result.fluid_decomposition_polyhedra[second];
            const auto& second_geometry =
                result.fluid_decomposition_geometries[second];
            bool connected = false;
            for (std::size_t first_face = 0;
                 first_face < first_polyhedron.faces.size() && !connected;
                 ++first_face) {
                const auto& first_face_record = first_polyhedron.faces[first_face];
                if (first_face_record.kind !=
                    PolyhedronFaceKind::internal_partition) {
                    continue;
                }
                for (std::size_t second_face = 0;
                     second_face < second_polyhedron.faces.size(); ++second_face) {
                    const auto& second_face_record =
                        second_polyhedron.faces[second_face];
                    if (second_face_record.kind ==
                            PolyhedronFaceKind::internal_partition &&
                        first_face_record.source_id ==
                            second_face_record.source_id &&
                        dot(first_geometry.faces[first_face].outward_normal,
                            second_geometry.faces[second_face].outward_normal) <
                            -1.0 + 512.0 *
                                       std::numeric_limits<double>::epsilon()) {
                        connected = true;
                        break;
                    }
                }
            }
            if (connected) {
                components.unite(first, second);
            }
        }
    }
    std::vector<std::size_t> roots;
    result.fluid_piece_component_ids.reserve(
        result.fluid_decomposition_polyhedra.size());
    for (std::size_t piece = 0;
         piece < result.fluid_decomposition_polyhedra.size(); ++piece) {
        const std::size_t root = components.find(piece);
        auto found = std::find(roots.begin(), roots.end(), root);
        if (found == roots.end()) {
            roots.push_back(root);
            found = roots.end() - 1;
        }
        result.fluid_piece_component_ids.push_back(
            static_cast<std::size_t>(found - roots.begin()));
    }
    result.fluid_component_count = roots.size();
    return result;
}

} // 命名空间 cartmesh
