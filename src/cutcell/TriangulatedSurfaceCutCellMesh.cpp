#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"

#include "cartmesh/cutcell/LocalTriangulatedCutCell.hpp"

#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/geometry/TriangleBoxIntersection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

constexpr std::array<std::uint8_t, 6> opposite_face = {1, 0, 3, 2, 5, 4};

struct BoxFace {
    double area{};
    Vec3 centroid{};
    Vec3 outward_normal{};
    std::vector<Vec3> vertices;
};

struct ArrangementPlane {
    Vec3 point{};
    Vec3 normal{};
    double offset{};
};

struct FluidComponentAnalysis {
    std::size_t component_count{};
    std::size_t convex_region_count{};
    bool resolved{true};
    std::uint64_t discarded_piece_count{};
    double discarded_piece_volume{};
    std::vector<FluidPolyhedronPiece> pieces;
};

[[nodiscard]] PointClassification classify_region_interior(
    const ConvexPolyhedron& region, const PolyhedronGeometry& geometry,
    const SurfaceClassifier& classifier, const Vec3& coordinate_origin) {
    const auto primary =
        classifier.classify(geometry.centroid + coordinate_origin).classification;
    if (primary == PointClassification::inside ||
        primary == PointClassification::outside) {
        return primary;
    }
    bool inside_seen = false;
    bool outside_seen = false;
    for (const auto& vertex : region.vertices) {
        // 凸区域质心到顶点的开线段仍在区域内。只在质心恰落在
        // 输入表面时取内部点，不修改几何或分类表面。
        const Vec3 sample = geometry.centroid * 0.875 + vertex * 0.125;
        const auto classification =
            classifier.classify(sample + coordinate_origin).classification;
        inside_seen = inside_seen || classification == PointClassification::inside;
        outside_seen = outside_seen || classification == PointClassification::outside;
        if (inside_seen && outside_seen) return PointClassification::conflict;
    }
    if (inside_seen) return PointClassification::inside;
    if (outside_seen) return PointClassification::outside;
    return PointClassification::conflict;
}

class DisjointSets {
  public:
    explicit DisjointSets(std::size_t size) : parent_(size) {
        for (std::size_t index = 0; index < size; ++index) parent_[index] = index;
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
        if (first != second) parent_[second] = first;
    }
  private:
    std::vector<std::size_t> parent_;
};

[[nodiscard]] std::array<BoxFace, 6> box_faces(const AABB& box) {
    const Vec3 minimum = box.minimum();
    const Vec3 maximum = box.maximum();
    const Vec3 center = box.center();
    const Vec3 extent = box.extent();
    return {{{extent.y * extent.z, {minimum.x, center.y, center.z}, {-1.0, 0.0, 0.0},
              {{minimum.x, minimum.y, minimum.z}, {minimum.x, minimum.y, maximum.z},
               {minimum.x, maximum.y, maximum.z}, {minimum.x, maximum.y, minimum.z}}},
             {extent.y * extent.z, {maximum.x, center.y, center.z}, {1.0, 0.0, 0.0},
              {{maximum.x, minimum.y, minimum.z}, {maximum.x, maximum.y, minimum.z},
               {maximum.x, maximum.y, maximum.z}, {maximum.x, minimum.y, maximum.z}}},
             {extent.x * extent.z, {center.x, minimum.y, center.z}, {0.0, -1.0, 0.0},
              {{minimum.x, minimum.y, minimum.z}, {maximum.x, minimum.y, minimum.z},
               {maximum.x, minimum.y, maximum.z}, {minimum.x, minimum.y, maximum.z}}},
             {extent.x * extent.z, {center.x, maximum.y, center.z}, {0.0, 1.0, 0.0},
              {{minimum.x, maximum.y, minimum.z}, {minimum.x, maximum.y, maximum.z},
               {maximum.x, maximum.y, maximum.z}, {maximum.x, maximum.y, minimum.z}}},
             {extent.x * extent.y, {center.x, center.y, minimum.z}, {0.0, 0.0, -1.0},
              {{minimum.x, minimum.y, minimum.z}, {minimum.x, maximum.y, minimum.z},
               {maximum.x, maximum.y, minimum.z}, {maximum.x, minimum.y, minimum.z}}},
             {extent.x * extent.y, {center.x, center.y, maximum.z}, {0.0, 0.0, 1.0},
              {{minimum.x, minimum.y, maximum.z}, {maximum.x, minimum.y, maximum.z},
               {maximum.x, maximum.y, maximum.z}, {minimum.x, maximum.y, maximum.z}}}}};
}

[[nodiscard]] bool positive_neighbor(const UniformCartesianGrid& grid,
                                     const CellKey& key, std::size_t axis,
                                     CellKey& neighbor) noexcept {
    neighbor = key;
    if (axis == 0) {
        if (key.i + 1U >= grid.nx()) return false;
        ++neighbor.i;
    } else if (axis == 1) {
        if (key.j + 1U >= grid.ny()) return false;
        ++neighbor.j;
    } else {
        if (key.k + 1U >= grid.nz()) return false;
        ++neighbor.k;
    }
    return true;
}

[[nodiscard]] AABB expanded_box(const AABB& box, double tolerance) {
    const Vec3 padding{tolerance, tolerance, tolerance};
    return AABB(box.minimum() - padding, box.maximum() + padding);
}

[[nodiscard]] double coordinate_ulp(double value) noexcept {
    const double magnitude = std::abs(value);
    return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) -
           magnitude;
}

[[nodiscard]] double bounds_coordinate_resolution(const AABB& bounds) noexcept {
    return std::max({coordinate_ulp(bounds.minimum().x),
                     coordinate_ulp(bounds.minimum().y),
                     coordinate_ulp(bounds.minimum().z),
                     coordinate_ulp(bounds.maximum().x),
                     coordinate_ulp(bounds.maximum().y),
                     coordinate_ulp(bounds.maximum().z)});
}

[[nodiscard]] bool canonical_plane_sign(const Vec3& normal) noexcept {
    if (normal.x != 0.0) return normal.x > 0.0;
    if (normal.y != 0.0) return normal.y > 0.0;
    return normal.z > 0.0;
}

[[nodiscard]] bool plane_less(const ArrangementPlane& first,
                              const ArrangementPlane& second) noexcept {
    if (first.normal.x != second.normal.x) return first.normal.x < second.normal.x;
    if (first.normal.y != second.normal.y) return first.normal.y < second.normal.y;
    if (first.normal.z != second.normal.z) return first.normal.z < second.normal.z;
    return first.offset < second.offset;
}

[[nodiscard]] std::vector<ArrangementPlane> intersecting_planes(
    const AABB& box, const TriangulatedSurfaceCutter& cutter,
    const Vec3& coordinate_origin, double length_tolerance) {
    std::vector<ArrangementPlane> planes;
    const AABB candidates = expanded_box(box, length_tolerance);
    for (const auto triangle_id : cutter.bvh().query(candidates)) {
        const auto& triangle = cutter.oriented_surface().triangles()[
            static_cast<std::size_t>(triangle_id)];
        if (!triangle_intersects_aabb(triangle, candidates)) continue;
        Vec3 normal = triangle.area_vector();
        const double magnitude = norm(normal);
        if (!(magnitude > 0.0)) continue;
        normal = normal / magnitude;
        if (!canonical_plane_sign(normal)) normal = normal * -1.0;
        const Vec3 point = triangle.vertices().front() - coordinate_origin;
        planes.push_back({point, normal, dot(normal, point)});
    }
    std::sort(planes.begin(), planes.end(), plane_less);
    const double angular_tolerance =
        512.0 * std::numeric_limits<double>::epsilon();
    std::vector<ArrangementPlane> unique;
    for (const auto& plane : planes) {
        if (unique.empty() ||
            norm(unique.back().normal - plane.normal) > angular_tolerance ||
            std::abs(unique.back().offset - plane.offset) > length_tolerance) {
            unique.push_back(plane);
        }
    }
    return unique;
}

[[nodiscard]] bool same_vertex_set(const ConvexPolyhedron& first,
                                   const PolyhedronFace& first_face,
                                   const ConvexPolyhedron& second,
                                   const PolyhedronFace& second_face,
                                   double tolerance) {
    if (first_face.vertex_indices.size() != second_face.vertex_indices.size()) return false;
    for (const auto first_vertex : first_face.vertex_indices) {
        const Vec3 point = first.vertices[first_vertex];
        const bool found = std::any_of(
            second_face.vertex_indices.begin(), second_face.vertex_indices.end(),
            [&](const std::uint32_t second_vertex) {
                return norm(point - second.vertices[second_vertex]) <= tolerance;
            });
        if (!found) return false;
    }
    return true;
}

[[nodiscard]] FluidComponentAnalysis analyze_fluid_components(
    const AABB& box, const TriangulatedSurfaceCutter& cutter,
    const SurfaceClassifier& classifier, double length_tolerance,
    double area_tolerance, bool retain_small_positive_pieces = false) {
    const Vec3 coordinate_origin = box.center();
    const AABB local_box(box.minimum() - coordinate_origin,
                         box.maximum() - coordinate_origin);
    const auto planes = intersecting_planes(
        box, cutter, coordinate_origin, length_tolerance);
    std::vector<ConvexPolyhedron> regions{make_box_polyhedron(local_box)};
    for (std::size_t plane_id = 0; plane_id < planes.size(); ++plane_id) {
        std::vector<ConvexPolyhedron> split_regions;
        split_regions.reserve(regions.size() * 2U);
        for (const auto& region : regions) {
            auto negative = clip_convex_polyhedron(
                region,
                OrientedHalfSpace(planes[plane_id].point,
                                  planes[plane_id].normal,
                                  static_cast<std::uint64_t>(plane_id)),
                length_tolerance, PolyhedronFaceKind::internal_partition);
            auto positive = clip_convex_polyhedron(
                region,
                OrientedHalfSpace(planes[plane_id].point,
                                  planes[plane_id].normal * -1.0,
                                  static_cast<std::uint64_t>(plane_id)),
                length_tolerance, PolyhedronFaceKind::internal_partition);
            if (!negative.empty()) split_regions.push_back(std::move(negative));
            if (!positive.empty()) split_regions.push_back(std::move(positive));
        }
        regions = std::move(split_regions);
    }

    std::vector<PolyhedronGeometry> geometries;
    std::vector<bool> fluid;
    geometries.reserve(regions.size());
    fluid.reserve(regions.size());
    FluidComponentAnalysis result;
    const double volume_tolerance =
        area_tolerance * std::max({local_box.extent().x, local_box.extent().y,
                                   local_box.extent().z});
    for (const auto& region : regions) {
        geometries.push_back(measure_polyhedron(region));
        if (!geometries.back().positive_volume ||
            (!(geometries.back().volume > 0.0)) ||
            (!retain_small_positive_pieces &&
             geometries.back().volume <= volume_tolerance)) {
            fluid.push_back(false);
            ++result.discarded_piece_count;
            result.discarded_piece_volume += geometries.back().volume;
            continue;
        }
        const auto classification = classify_region_interior(
            region, geometries.back(), classifier, coordinate_origin);
        if (classification == PointClassification::outside) {
            fluid.push_back(true);
            ++result.convex_region_count;
        } else if (classification == PointClassification::inside) {
            fluid.push_back(false);
        } else {
            fluid.push_back(false);
            result.resolved = false;
        }
    }
    DisjointSets components(regions.size());
    const double normal_tolerance =
        1024.0 * std::numeric_limits<double>::epsilon();
    for (std::size_t first = 0; first < regions.size(); ++first) {
        if (!fluid[first]) continue;
        for (std::size_t second = first + 1; second < regions.size(); ++second) {
            if (!fluid[second]) continue;
            bool adjacent = false;
            for (std::size_t first_face = 0;
                 first_face < regions[first].faces.size() && !adjacent;
                 ++first_face) {
                const auto& first_record = regions[first].faces[first_face];
                if (first_record.kind != PolyhedronFaceKind::internal_partition ||
                    geometries[first].faces[first_face].area <= area_tolerance) continue;
                for (std::size_t second_face = 0;
                     second_face < regions[second].faces.size(); ++second_face) {
                    const auto& second_record = regions[second].faces[second_face];
                    if (second_record.kind != PolyhedronFaceKind::internal_partition ||
                        first_record.source_id != second_record.source_id ||
                        geometries[second].faces[second_face].area <= area_tolerance) continue;
                    if (dot(geometries[first].faces[first_face].outward_normal,
                            geometries[second].faces[second_face].outward_normal) >
                        -1.0 + normal_tolerance) continue;
                    if (std::abs(geometries[first].faces[first_face].area -
                                 geometries[second].faces[second_face].area) >
                        area_tolerance) continue;
                    if (norm(geometries[first].faces[first_face].centroid -
                             geometries[second].faces[second_face].centroid) >
                        length_tolerance) continue;
                    if (same_vertex_set(regions[first], first_record, regions[second],
                                        second_record, length_tolerance)) {
                        adjacent = true;
                        break;
                    }
                }
            }
            if (adjacent) components.unite(first, second);
        }
    }
    std::vector<std::size_t> roots;
    for (std::size_t region = 0; region < regions.size(); ++region) {
        if (!fluid[region]) continue;
        const std::size_t root = components.find(region);
        if (std::find(roots.begin(), roots.end(), root) == roots.end()) roots.push_back(root);
    }
    result.component_count = roots.size();
    result.pieces.reserve(result.convex_region_count);
    for (std::size_t region = 0; region < regions.size(); ++region) {
        if (!fluid[region]) continue;
        const std::size_t root = components.find(region);
        const auto found = std::find(roots.begin(), roots.end(), root);
        for (auto& vertex : regions[region].vertices) {
            vertex = vertex + coordinate_origin;
        }
        geometries[region].centroid =
            geometries[region].centroid + coordinate_origin;
        for (auto& face : geometries[region].faces) {
            face.centroid = face.centroid + coordinate_origin;
        }
        result.pieces.push_back(
            {std::move(regions[region]), std::move(geometries[region]),
             static_cast<std::size_t>(found - roots.begin())});
    }
    return result;
}

[[nodiscard]] FluidCellGeometry full_fluid_cell(std::uint64_t cell_id,
                                                 const AABB& box) {
    FluidCellGeometry result;
    result.background_cell_id = cell_id;
    result.volume = box.volume();
    result.volume_fraction = 1.0;
    result.centroid = box.center();
    result.fluid_piece_count = 1;
    result.fluid_component_count = 1;
    const auto faces = box_faces(box);
    for (std::size_t face = 0; face < 6; ++face) {
        result.cartesian_faces[face] =
            {faces[face].area, 1.0, faces[face].centroid,
             faces[face].outward_normal, {faces[face].vertices}};
    }
    return result;
}

void recompute_closure(FluidCellGeometry& cell, double length_tolerance) {
    cell.oriented_area_vector_sum = {};
    for (const auto& face : cell.cartesian_faces) {
        cell.oriented_area_vector_sum =
            cell.oriented_area_vector_sum + face.outward_normal * face.area;
    }
    for (const auto& face : cell.embedded_boundary_faces) {
        cell.oriented_area_vector_sum =
            cell.oriented_area_vector_sum + face.outward_normal * face.area;
    }
    cell.area_vector_closure_residual = norm(cell.oriented_area_vector_sum);
    const auto edge_closure =
        analyze_boundary_edge_closure(cell, length_tolerance);
    cell.boundary_edge_imbalance_count = edge_closure.imbalanced_edge_count;
    cell.boundary_edge_closed = edge_closure.closed;
}

void rebuild_cartesian_apertures_from_fluid_pieces(
    FluidCellGeometry& cell, const AABB& box, double length_tolerance,
    double area_tolerance) {
    if (!cell.cut || cell.fluid_polyhedron_pieces.empty()) return;
    const auto full_faces = box_faces(box);
    std::array<Vec3, 6> first_moments{};
    for (std::size_t local_face = 0; local_face < 6; ++local_face) {
        auto& aperture = cell.cartesian_faces[local_face];
        aperture.area = 0.0;
        aperture.area_fraction = 0.0;
        aperture.centroid = full_faces[local_face].centroid;
        aperture.outward_normal = full_faces[local_face].outward_normal;
        aperture.oriented_boundary_loops.clear();
    }
    for (const auto& piece : cell.fluid_polyhedron_pieces) {
        if (piece.polyhedron.faces.size() != piece.geometry.faces.size()) {
            throw std::runtime_error(
                "Cut-cell 凸片面记录与几何测量数量不一致");
        }
        for (std::size_t face_index = 0;
             face_index < piece.polyhedron.faces.size(); ++face_index) {
            const auto& face = piece.polyhedron.faces[face_index];
            if (face.kind != PolyhedronFaceKind::cartesian) continue;
            if (face.source_id >= 6U) {
                throw std::runtime_error("Cut-cell Cartesian 面编号超出 0..5");
            }
            const auto local_face = static_cast<std::size_t>(face.source_id);
            const auto& geometry = piece.geometry.faces[face_index];
            if (geometry.area <= area_tolerance) continue;
            if (dot(geometry.outward_normal,
                    full_faces[local_face].outward_normal) <= 0.0) {
                throw std::runtime_error(
                    "Cut-cell 凸片的 Cartesian 面法向与背景盒不一致");
            }
            auto& aperture = cell.cartesian_faces[local_face];
            aperture.area += geometry.area;
            first_moments[local_face] =
                first_moments[local_face] + geometry.centroid * geometry.area;
            std::vector<Vec3> loop;
            loop.reserve(face.vertex_indices.size());
            for (const auto vertex : face.vertex_indices) {
                loop.push_back(piece.polyhedron.vertices[vertex]);
            }
            aperture.oriented_boundary_loops.push_back(std::move(loop));
        }
    }
    const auto point_on_face = [&](const Vec3& point, std::size_t local_face) {
        if (local_face == 0) {
            return std::abs(point.x - box.minimum().x) <= length_tolerance;
        }
        if (local_face == 1) {
            return std::abs(point.x - box.maximum().x) <= length_tolerance;
        }
        if (local_face == 2) {
            return std::abs(point.y - box.minimum().y) <= length_tolerance;
        }
        if (local_face == 3) {
            return std::abs(point.y - box.maximum().y) <= length_tolerance;
        }
        if (local_face == 4) {
            return std::abs(point.z - box.minimum().z) <= length_tolerance;
        }
        return std::abs(point.z - box.maximum().z) <= length_tolerance;
    };
    // 当实体边界恰好位于当前 Cartesian 面且流体在本单元内时，凸片的
    // box face 仍覆盖该区域。该区域是 wall 而不是跨单元开口，必须从
    // aperture 中扣除；反向环与随后写入的 embedded face 共同保持边链闭合。
    for (const auto& embedded : cell.embedded_boundary_faces) {
        for (std::size_t local_face = 0; local_face < 6; ++local_face) {
            if (dot(embedded.outward_normal,
                    full_faces[local_face].outward_normal) <= 0.5) {
                continue;
            }
            if (embedded.vertices.empty() ||
                !std::all_of(embedded.vertices.begin(), embedded.vertices.end(),
                             [&](const Vec3& point) {
                                 return point_on_face(point, local_face);
                             })) {
                continue;
            }
            auto& aperture = cell.cartesian_faces[local_face];
            aperture.area -= embedded.area;
            first_moments[local_face] =
                first_moments[local_face] - embedded.centroid * embedded.area;
            auto loop = embedded.vertices;
            std::reverse(loop.begin(), loop.end());
            aperture.oriented_boundary_loops.push_back(std::move(loop));
            break;
        }
    }
    for (std::size_t local_face = 0; local_face < 6; ++local_face) {
        auto& aperture = cell.cartesian_faces[local_face];
        const double full_area = full_faces[local_face].area;
        if (aperture.area > full_area + area_tolerance) {
            throw std::runtime_error(
                "Cut-cell 凸片在 Cartesian 面上的总开口面积超出背景面");
        }
        if (aperture.area < -area_tolerance) {
            throw std::runtime_error(
                "Cut-cell 共面 wall 面积大于凸片 Cartesian 面面积");
        }
        if (std::abs(aperture.area) <= area_tolerance) aperture.area = 0.0;
        aperture.area = std::min(aperture.area, full_area);
        aperture.area_fraction = aperture.area / full_area;
        if (aperture.area > 0.0) {
            aperture.centroid = first_moments[local_face] / aperture.area;
        }
    }
}

} // 匿名命名空间

LocalTriangulatedCutCellBuilder::LocalTriangulatedCutCellBuilder(
    const TriangulatedSurfaceCutter& cutter, const AABB& domain,
    const Vec3& representative_spacing, double geometric_tolerance,
    bool retain_small_positive_pieces)
    : cutter_(cutter), classifier_(cutter.bvh(), cutter.length_tolerance()),
      retain_small_positive_pieces_(retain_small_positive_pieces) {
    if (geometric_tolerance < 0.0 || !std::isfinite(geometric_tolerance)) {
        throw std::invalid_argument("\u5c40\u90e8 Cut-cell \u5bb9\u5dee\u5fc5\u987b\u662f\u975e\u8d1f\u6709\u9650\u6570");
    }
    if (!(representative_spacing.x > 0.0) ||
        !(representative_spacing.y > 0.0) ||
        !(representative_spacing.z > 0.0)) {
        throw std::invalid_argument("\u5c40\u90e8 Cut-cell \u4ee3\u8868\u5355\u5143\u5c3a\u5bf8\u5fc5\u987b\u4e3a\u6b63");
    }
    length_tolerance_ = std::max(
        {geometric_tolerance, cutter.length_tolerance(),
         256.0 * std::numeric_limits<double>::epsilon() *
             norm(domain.extent())});
    const double face_area_scale = std::max(
        {representative_spacing.x * representative_spacing.y,
         representative_spacing.x * representative_spacing.z,
         representative_spacing.y * representative_spacing.z});
    area_tolerance_ = std::max(
        length_tolerance_ * length_tolerance_,
        4096.0 * std::numeric_limits<double>::epsilon() * face_area_scale);
    topology_length_tolerance_ = std::max(
        length_tolerance_, 4.0 * bounds_coordinate_resolution(domain));
    closure_tolerance_ = std::max(
        area_tolerance_, topology_length_tolerance_ *
                             std::max({representative_spacing.x,
                                       representative_spacing.y,
                                       representative_spacing.z}));
}

LocalTriangulatedCutCellResult LocalTriangulatedCutCellBuilder::build(
    std::uint64_t background_cell_id, const AABB& box) const {
    if (!box.has_positive_volume()) {
        throw std::invalid_argument("\u5c40\u90e8 Cut-cell \u80cc\u666f\u76d2\u5fc5\u987b\u5177\u6709\u6b63\u4f53\u79ef");
    }
    auto components = analyze_fluid_components(
        box, cutter_, classifier_, length_tolerance_, area_tolerance_,
        retain_small_positive_pieces_);
    LocalTriangulatedCutCellResult result;
    result.discarded_numerical_piece_count = components.discarded_piece_count;
    result.discarded_numerical_piece_volume = components.discarded_piece_volume;
    result.classification_conflict = !components.resolved;
    result.component_analysis_pending =
        !components.resolved ||
        (components.pieces.empty() && components.discarded_piece_count > 0);

    double fluid_volume = 0.0;
    Vec3 first_moment{};
    for (const auto& piece : components.pieces) {
        fluid_volume += piece.geometry.volume;
        first_moment = first_moment +
                       piece.geometry.centroid * piece.geometry.volume;
    }
    const double volume_tolerance = area_tolerance_ * std::max(
        {box.extent().x, box.extent().y, box.extent().z});
    if (fluid_volume > box.volume() + volume_tolerance) {
        throw std::runtime_error(
            "\u5c40\u90e8 arrangement \u6d41\u4f53\u4f53\u79ef\u8d85\u51fa\u80cc\u666f\u76d2");
    }
    if (fluid_volume <= volume_tolerance) return result;
    if (std::abs(fluid_volume - box.volume()) <= volume_tolerance) {
        fluid_volume = box.volume();
    }

    const auto boundary = cutter_.clip_boundary_faces(box);
    const bool volume_cut =
        fluid_volume < box.volume() - volume_tolerance;
    result.has_fluid = true;
    if (!volume_cut && boundary.faces.empty()) {
        result.cell = full_fluid_cell(background_cell_id, box);
        return result;
    }

    auto& cell = result.cell;
    cell.background_cell_id = background_cell_id;
    cell.volume = fluid_volume;
    cell.volume_fraction = fluid_volume / box.volume();
    cell.centroid = first_moment / fluid_volume;
    cell.embedded_boundary_faces = boundary.faces;
    cell.fluid_piece_count = components.convex_region_count;
    cell.fluid_component_count = components.component_count;
    cell.fluid_polyhedron_pieces = std::move(components.pieces);
    cell.discarded_numerical_piece_count = components.discarded_piece_count;
    cell.discarded_numerical_piece_volume = components.discarded_piece_volume;
    cell.component_analysis_resolved = components.resolved;
    // 即使 wall 恰好与 Cartesian 面共面而当前单元体积仍为
    // 全流体，也要先用显式凸片重建 aperture 并扣除 wall。
    // cut 计数继续严格表示体积部分切割，共面 wall 另由
    // boundary_cell_count 记录。
    cell.cut = true;
    rebuild_cartesian_apertures_from_fluid_pieces(
        cell, box, topology_length_tolerance_, area_tolerance_);
    cell.cut = volume_cut;
    if (!volume_cut) {
        cell.fluid_polyhedron_pieces.clear();
        cell.fluid_piece_count = 1;
        cell.fluid_component_count = 1;
    }
    recompute_closure(cell, topology_length_tolerance_);
    result.embedded_boundary_area = boundary.total_area;
    return result;
}

ConvexCutCellMesh build_triangulated_cut_cell_mesh(
    const UniformCartesianGrid& grid,
    const TriangulatedSurfaceCutter& cutter,
    double geometric_tolerance) {
    if (geometric_tolerance < 0.0 || !std::isfinite(geometric_tolerance)) {
        throw std::invalid_argument("通用 Cut-cell 网格容差必须为非负有限数");
    }
    ConvexCutCellMesh result;
    const double length_tolerance = std::max(
        {geometric_tolerance, cutter.length_tolerance(),
         256.0 * std::numeric_limits<double>::epsilon() *
             norm(grid.domain().extent())});
    const double face_area_scale = std::max(
        {grid.spacing().x * grid.spacing().y,
         grid.spacing().x * grid.spacing().z,
         grid.spacing().y * grid.spacing().z});
    const double area_tolerance = std::max(
        length_tolerance * length_tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() * face_area_scale);
    const double topology_length_tolerance = std::max(
        length_tolerance, 4.0 * bounds_coordinate_resolution(grid.domain()));
    const double closure_tolerance = std::max(
        area_tolerance,
        topology_length_tolerance *
            std::max({grid.spacing().x, grid.spacing().y, grid.spacing().z}));
    SurfaceClassifier classifier(cutter.bvh(), cutter.length_tolerance());
    const LocalTriangulatedCutCellBuilder local_builder(
        cutter, grid.domain(), grid.spacing(), geometric_tolerance);
    std::vector<std::size_t> active_index(
        static_cast<std::size_t>(grid.cell_count()),
        std::numeric_limits<std::size_t>::max());
    result.fluid_cells.reserve(static_cast<std::size_t>(grid.cell_count()));

    for (std::uint64_t cell_id = 0; cell_id < grid.cell_count(); ++cell_id) {
        const CellKey key = grid.cell_key(cell_id);
        const AABB box = grid.cell_bounds(key);
        if (!cutter.bvh().intersects_surface(expanded_box(box, length_tolerance))) {
            const auto point = classifier.classify(box.center()).classification;
            if (point == PointClassification::inside) {
                ++result.full_solid_cell_count;
                continue;
            }
            if (point != PointClassification::outside) {
                ++result.classification_conflict_count;
            }
            active_index[static_cast<std::size_t>(cell_id)] = result.fluid_cells.size();
            result.fluid_cells.push_back(full_fluid_cell(cell_id, box));
            result.total_fluid_volume += box.volume();
            ++result.full_fluid_cell_count;
            continue;
        }
        auto local = local_builder.build(cell_id, box);
        if (!local.has_fluid) {
            ++result.full_solid_cell_count;
            if (local.classification_conflict) {
                ++result.classification_conflict_count;
            }
            result.discarded_numerical_piece_count +=
                local.discarded_numerical_piece_count;
            result.discarded_numerical_piece_volume +=
                local.discarded_numerical_piece_volume;
            if (local.component_analysis_pending) {
                ++result.component_analysis_pending_cell_count;
            }
            continue;
        }
        active_index[static_cast<std::size_t>(cell_id)] = result.fluid_cells.size();
        result.fluid_cells.push_back(std::move(local.cell));
        auto& cell = result.fluid_cells.back();
        result.discarded_numerical_piece_count +=
            local.discarded_numerical_piece_count;
        result.discarded_numerical_piece_volume +=
            local.discarded_numerical_piece_volume;
        if (local.component_analysis_pending) {
            ++result.component_analysis_pending_cell_count;
        }
        if (cell.cut) {
            if (cell.component_analysis_resolved &&
                cell.fluid_component_count > 1) {
                ++result.disconnected_fluid_cell_count;
            }
        }
        recompute_closure(cell, topology_length_tolerance);
        result.total_fluid_volume += cell.volume;
        result.total_embedded_boundary_area += local.embedded_boundary_area;
        if (!cell.embedded_boundary_faces.empty()) ++result.boundary_cell_count;
        if (cell.cut) {
            ++result.cut_cell_count;
        } else {
            ++result.full_fluid_cell_count;
        }
    }

    for (auto& cell : result.fluid_cells) {
        recompute_closure(cell, topology_length_tolerance);
        result.maximum_cell_area_closure_residual =
            std::max(result.maximum_cell_area_closure_residual,
                     cell.area_vector_closure_residual);
        if (cell.area_vector_closure_residual > closure_tolerance ||
            !cell.boundary_edge_closed) {
            ++result.nonclosed_cell_count;
        }
        if (!(cell.volume > 0.0)) {
            ++result.negative_volume_cell_count;
        }
    }

    for (std::uint64_t cell_id = 0; cell_id < grid.cell_count(); ++cell_id) {
        const CellKey key = grid.cell_key(cell_id);
        const std::size_t first_index = active_index[static_cast<std::size_t>(cell_id)];
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CellKey neighbor_key{};
            if (!positive_neighbor(grid, key, axis, neighbor_key)) continue;
            const std::uint64_t neighbor_id = grid.linear_id(neighbor_key);
            const std::size_t second_index =
                active_index[static_cast<std::size_t>(neighbor_id)];
            const std::uint8_t first_face = static_cast<std::uint8_t>(2 * axis + 1);
            const std::uint8_t second_face = opposite_face[first_face];
            const double first_area =
                first_index == std::numeric_limits<std::size_t>::max()
                    ? 0.0
                    : result.fluid_cells[first_index].cartesian_faces[first_face].area;
            const double second_area =
                second_index == std::numeric_limits<std::size_t>::max()
                    ? 0.0
                    : result.fluid_cells[second_index].cartesian_faces[second_face].area;
            const double area_mismatch = std::abs(first_area - second_area);
            result.maximum_shared_face_area_mismatch =
                std::max(result.maximum_shared_face_area_mismatch, area_mismatch);
            if (area_mismatch > closure_tolerance) {
                ++result.shared_face_mismatch_count;
            }
            if (first_index == std::numeric_limits<std::size_t>::max() ||
                second_index == std::numeric_limits<std::size_t>::max() ||
                std::max(first_area, second_area) <= area_tolerance) {
                continue;
            }
            const auto& first =
                result.fluid_cells[first_index].cartesian_faces[first_face];
            const auto& second =
                result.fluid_cells[second_index].cartesian_faces[second_face];
            const double centroid_mismatch = norm(first.centroid - second.centroid);
            result.maximum_shared_face_centroid_mismatch =
                std::max(result.maximum_shared_face_centroid_mismatch,
                         centroid_mismatch);
            const AABB first_box = grid.cell_bounds(key);
            const Vec3 reference = box_faces(first_box)[first_face].centroid;
            const double first_moment_mismatch = norm(
                (first.centroid - reference) * first_area -
                (second.centroid - reference) * second_area);
            result.maximum_shared_face_first_moment_mismatch = std::max(
                result.maximum_shared_face_first_moment_mismatch,
                first_moment_mismatch);
            const double first_moment_tolerance =
                area_tolerance * std::max({first_box.extent().x,
                                           first_box.extent().y,
                                           first_box.extent().z}) +
                topology_length_tolerance * std::max(first_area, second_area);
            if (first_moment_mismatch > first_moment_tolerance) {
                ++result.shared_face_mismatch_count;
            }
            result.internal_faces.push_back(
                {cell_id, neighbor_id, first_index, second_index,
                 first_face, second_face, 0.5 * (first_area + second_area),
                 (first.centroid + second.centroid) * 0.5,
                 first.outward_normal, area_mismatch, centroid_mismatch,
                 first_moment_mismatch});
        }
    }
    assign_global_fluid_regions(result, closure_tolerance);
    return result;
}

namespace {

struct AdaptiveReuseContext {
    const LinearOctree& old_tree;
    const ConvexCutCellMesh& old_mesh;
    const std::optional<AABB>& affected_bounds;
    IncrementalCutCellBuildStatistics& statistics;
    std::vector<std::uint8_t>& rebuilt_leaf_mask;
};

[[nodiscard]] bool geometry_can_be_reused(
    const AABB& box, double length_tolerance,
    const std::optional<AABB>& affected_bounds) {
    return !affected_bounds ||
           !affected_bounds->intersects(expanded_box(box, length_tolerance));
}

ConvexCutCellMesh build_adaptive_triangulated_cut_cell_mesh(
    const LinearOctree& tree, const TriangulatedSurfaceCutter& cutter,
    double geometric_tolerance, AdaptiveReuseContext* reuse) {
    if (geometric_tolerance < 0.0 || !std::isfinite(geometric_tolerance)) {
        throw std::invalid_argument("自适应 Cut-cell 网格容差必须为非负有限数");
    }
    ConvexCutCellMesh result;
    const double length_tolerance = std::max(
        {geometric_tolerance, cutter.length_tolerance(),
         256.0 * std::numeric_limits<double>::epsilon() *
             norm(tree.domain().extent())});
    const double minimum_scale = std::ldexp(
        std::max({tree.domain().extent().x, tree.domain().extent().y,
                  tree.domain().extent().z}),
        -static_cast<int>(tree.maximum_level()));
    const double face_area_scale = minimum_scale * minimum_scale;
    const double area_tolerance = std::max(
        length_tolerance * length_tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() * face_area_scale);
    const double topology_length_tolerance = std::max(
        length_tolerance, 4.0 * bounds_coordinate_resolution(tree.domain()));
    const double maximum_cell_scale = std::ldexp(
        std::max({tree.domain().extent().x, tree.domain().extent().y,
                  tree.domain().extent().z}),
        -static_cast<int>(tree.base_level()));
    const double closure_tolerance = std::max(
        area_tolerance, topology_length_tolerance * maximum_cell_scale);
    SurfaceClassifier classifier(cutter.bvh(), cutter.length_tolerance());
    const LocalTriangulatedCutCellBuilder local_builder(
        cutter, tree.domain(),
        Vec3{minimum_scale, minimum_scale, minimum_scale},
        geometric_tolerance);
    std::vector<std::size_t> active_index(
        static_cast<std::size_t>(tree.leaf_count()),
        std::numeric_limits<std::size_t>::max());
    result.fluid_cells.reserve(static_cast<std::size_t>(tree.leaf_count()));

    std::vector<const FluidCellGeometry*> old_fluid_cells;
    if (reuse) {
        reuse->statistics.old_leaf_count = reuse->old_tree.leaf_count();
        reuse->statistics.new_leaf_count = tree.leaf_count();
        reuse->rebuilt_leaf_mask.assign(static_cast<std::size_t>(tree.leaf_count()),
                                        std::uint8_t{1});
        old_fluid_cells.assign(
            static_cast<std::size_t>(reuse->old_tree.leaf_count()), nullptr);
        for (const auto& cell : reuse->old_mesh.fluid_cells) {
            if (cell.background_cell_id >= reuse->old_tree.leaf_count()) {
                throw std::invalid_argument(
                    "旧 Cut-cell 网格的 background ID 超出旧八叉树叶范围");
            }
            old_fluid_cells[static_cast<std::size_t>(cell.background_cell_id)] =
                &cell;
        }
    }

    for (std::uint64_t leaf_id = 0; leaf_id < tree.leaf_count(); ++leaf_id) {
        const OctreeNodeCode code = tree.leaf_code(leaf_id);
        const AABB box = tree.cell_bounds(code);
        if (reuse) {
            const auto old_leaf_id = reuse->old_tree.find_leaf(code);
            if (old_leaf_id && geometry_can_be_reused(
                                   box, length_tolerance,
                                   reuse->affected_bounds)) {
                ++reuse->statistics.reused_leaf_count;
                reuse->rebuilt_leaf_mask[static_cast<std::size_t>(leaf_id)] = 0;
                const auto* old_cell = old_fluid_cells[
                    static_cast<std::size_t>(*old_leaf_id)];
                if (!old_cell) {
                    ++result.full_solid_cell_count;
                    ++reuse->statistics.reused_solid_cell_count;
                    continue;
                }
                active_index[static_cast<std::size_t>(leaf_id)] =
                    result.fluid_cells.size();
                result.fluid_cells.push_back(*old_cell);
                auto& cell = result.fluid_cells.back();
                cell.background_cell_id = leaf_id;
                result.total_fluid_volume += cell.volume;
                for (const auto& face : cell.embedded_boundary_faces) {
                    result.total_embedded_boundary_area += face.area;
                }
                result.discarded_numerical_piece_count +=
                    cell.discarded_numerical_piece_count;
                result.discarded_numerical_piece_volume +=
                    cell.discarded_numerical_piece_volume;
                if (!cell.component_analysis_resolved) {
                    ++result.component_analysis_pending_cell_count;
                } else if (cell.fluid_component_count > 1) {
                    ++result.disconnected_fluid_cell_count;
                }
                if (!cell.embedded_boundary_faces.empty()) {
                    ++result.boundary_cell_count;
                }
                if (cell.cut) ++result.cut_cell_count;
                else ++result.full_fluid_cell_count;
                ++reuse->statistics.reused_fluid_cell_count;
                continue;
            }
            ++reuse->statistics.rebuilt_leaf_count;
        }
        if (!cutter.bvh().intersects_surface(expanded_box(box, length_tolerance))) {
            const auto point = classifier.classify(box.center()).classification;
            if (point == PointClassification::inside) {
                ++result.full_solid_cell_count;
                continue;
            }
            if (point != PointClassification::outside) {
                ++result.classification_conflict_count;
            }
            active_index[static_cast<std::size_t>(leaf_id)] =
                result.fluid_cells.size();
            result.fluid_cells.push_back(full_fluid_cell(leaf_id, box));
            result.total_fluid_volume += box.volume();
            ++result.full_fluid_cell_count;
            continue;
        }
        auto local = local_builder.build(leaf_id, box);
        if (!local.has_fluid) {
            ++result.full_solid_cell_count;
            if (local.classification_conflict) {
                ++result.classification_conflict_count;
            }
            result.discarded_numerical_piece_count +=
                local.discarded_numerical_piece_count;
            result.discarded_numerical_piece_volume +=
                local.discarded_numerical_piece_volume;
            if (local.component_analysis_pending) {
                ++result.component_analysis_pending_cell_count;
            }
            continue;
        }
        active_index[static_cast<std::size_t>(leaf_id)] =
            result.fluid_cells.size();
        result.fluid_cells.push_back(std::move(local.cell));
        auto& cell = result.fluid_cells.back();
        result.discarded_numerical_piece_count +=
            local.discarded_numerical_piece_count;
        result.discarded_numerical_piece_volume +=
            local.discarded_numerical_piece_volume;
        if (local.component_analysis_pending) {
            ++result.component_analysis_pending_cell_count;
        }
        if (cell.cut) {
            if (cell.component_analysis_resolved &&
                cell.fluid_component_count > 1) {
                ++result.disconnected_fluid_cell_count;
            }
        }
        recompute_closure(cell, topology_length_tolerance);
        result.total_fluid_volume += cell.volume;
        result.total_embedded_boundary_area += local.embedded_boundary_area;
        if (!cell.embedded_boundary_faces.empty()) ++result.boundary_cell_count;
        if (cell.cut) {
            ++result.cut_cell_count;
            if (reuse) ++reuse->statistics.rebuilt_cut_cell_count;
        }
        else ++result.full_fluid_cell_count;
    }

    for (auto& cell : result.fluid_cells) {
        recompute_closure(cell, topology_length_tolerance);
        result.maximum_cell_area_closure_residual =
            std::max(result.maximum_cell_area_closure_residual,
                     cell.area_vector_closure_residual);
        if (cell.area_vector_closure_residual > closure_tolerance ||
            !cell.boundary_edge_closed) {
            ++result.nonclosed_cell_count;
        }
        if (!(cell.volume > 0.0)) ++result.negative_volume_cell_count;
    }

    constexpr std::array<FaceDirection, 6> directions = {
        FaceDirection::negative_x, FaceDirection::positive_x,
        FaceDirection::negative_y, FaceDirection::positive_y,
        FaceDirection::negative_z, FaceDirection::positive_z};
    for (std::uint64_t leaf_id = 0; leaf_id < tree.leaf_count(); ++leaf_id) {
        const OctreeNodeCode code = tree.leaf_code(leaf_id);
        const auto node = decode_octree_node(code);
        const std::size_t first_index = active_index[static_cast<std::size_t>(leaf_id)];
        for (std::size_t local_face = 0; local_face < 6; ++local_face) {
            const auto neighbors = tree.face_neighbors(code, directions[local_face]);
            if (neighbors.empty()) continue;
            bool has_coarser_neighbor = false;
            bool has_finer_neighbor = false;
            std::vector<std::uint64_t> neighbor_ids;
            neighbor_ids.reserve(neighbors.size());
            for (const auto neighbor_code : neighbors) {
                const auto neighbor_node = decode_octree_node(neighbor_code);
                has_coarser_neighbor = has_coarser_neighbor || neighbor_node.level < node.level;
                has_finer_neighbor = has_finer_neighbor || neighbor_node.level > node.level;
                const auto neighbor_id = tree.find_leaf(neighbor_code);
                if (!neighbor_id) {
                    throw std::runtime_error("八叉树面邻居不在叶数组中");
                }
                neighbor_ids.push_back(*neighbor_id);
            }
            if (has_coarser_neighbor) continue;
            if (!has_finer_neighbor && leaf_id > neighbor_ids.front()) continue;

            const std::uint8_t opposite = opposite_face[local_face];
            const AABB first_box = tree.cell_bounds(code);
            const Vec3 reference = box_faces(first_box)[local_face].centroid;
            const double first_area =
                first_index == std::numeric_limits<std::size_t>::max()
                    ? 0.0
                    : result.fluid_cells[first_index].cartesian_faces[local_face].area;
            double neighbor_area_sum = 0.0;
            Vec3 neighbor_first_moment{};
            for (const auto neighbor_id : neighbor_ids) {
                const std::size_t second_index =
                    active_index[static_cast<std::size_t>(neighbor_id)];
                if (second_index == std::numeric_limits<std::size_t>::max()) continue;
                const auto& aperture =
                    result.fluid_cells[second_index].cartesian_faces[opposite];
                neighbor_area_sum += aperture.area;
                neighbor_first_moment =
                    neighbor_first_moment +
                    (aperture.centroid - reference) * aperture.area;
                if (first_index != std::numeric_limits<std::size_t>::max() &&
                    aperture.area > area_tolerance) {
                    const auto& first = result.fluid_cells[first_index];
                    result.internal_faces.push_back(
                        {leaf_id, neighbor_id, first_index, second_index,
                         static_cast<std::uint8_t>(local_face), opposite,
                         aperture.area, aperture.centroid,
                         first.cartesian_faces[local_face].outward_normal,
                         0.0, 0.0, 0.0});
                }
            }
            const double area_mismatch = std::abs(first_area - neighbor_area_sum);
            result.maximum_shared_face_area_mismatch =
                std::max(result.maximum_shared_face_area_mismatch, area_mismatch);
            if (area_mismatch > closure_tolerance) {
                ++result.shared_face_mismatch_count;
                continue;
            }
            if (first_area > area_tolerance && neighbor_area_sum > area_tolerance &&
                first_index != std::numeric_limits<std::size_t>::max()) {
                const Vec3 neighbor_centroid =
                    reference + neighbor_first_moment / neighbor_area_sum;
                const auto& first_aperture =
                    result.fluid_cells[first_index].cartesian_faces[local_face];
                const double centroid_mismatch = norm(
                    first_aperture.centroid - neighbor_centroid);
                result.maximum_shared_face_centroid_mismatch =
                    std::max(result.maximum_shared_face_centroid_mismatch,
                             centroid_mismatch);
                const double first_moment_mismatch = norm(
                    (first_aperture.centroid - reference) * first_area -
                    neighbor_first_moment);
                result.maximum_shared_face_first_moment_mismatch = std::max(
                    result.maximum_shared_face_first_moment_mismatch,
                    first_moment_mismatch);
                const double first_moment_tolerance =
                    area_tolerance * std::max({first_box.extent().x,
                                               first_box.extent().y,
                                               first_box.extent().z}) +
                    topology_length_tolerance *
                        std::max(first_area, neighbor_area_sum);
                if (first_moment_mismatch > first_moment_tolerance) {
                    ++result.shared_face_mismatch_count;
                }
            }
        }
    }
    assign_global_fluid_regions(result, closure_tolerance);
    return result;
}

} // namespace

ConvexCutCellMesh build_triangulated_cut_cell_mesh(
    const LinearOctree& tree, const TriangulatedSurfaceCutter& cutter,
    double geometric_tolerance) {
    return build_adaptive_triangulated_cut_cell_mesh(
        tree, cutter, geometric_tolerance, nullptr);
}

IncrementalCutCellBuildResult build_incremental_triangulated_cut_cell_mesh(
    const LinearOctree& old_tree, const ConvexCutCellMesh& old_mesh,
    const LinearOctree& new_tree,
    const TriangulatedSurfaceCutter& new_cutter,
    const std::optional<AABB>& affected_bounds,
    double geometric_tolerance) {
    if (old_tree.domain().minimum().x != new_tree.domain().minimum().x ||
        old_tree.domain().minimum().y != new_tree.domain().minimum().y ||
        old_tree.domain().minimum().z != new_tree.domain().minimum().z ||
        old_tree.domain().maximum().x != new_tree.domain().maximum().x ||
        old_tree.domain().maximum().y != new_tree.domain().maximum().y ||
        old_tree.domain().maximum().z != new_tree.domain().maximum().z ||
        old_tree.base_level() != new_tree.base_level() ||
        old_tree.maximum_level() != new_tree.maximum_level()) {
        throw std::invalid_argument(
            "增量 Cut-cell 重构要求旧、新八叉树使用相同固定域和层级范围");
    }
    IncrementalCutCellBuildResult result;
    AdaptiveReuseContext reuse{old_tree, old_mesh, affected_bounds,
                               result.statistics,
                               result.rebuilt_leaf_mask};
    result.mesh = build_adaptive_triangulated_cut_cell_mesh(
        new_tree, new_cutter, geometric_tolerance, &reuse);
    result.statistics.geometry_reuse_fraction =
        new_tree.leaf_count() == 0
            ? 0.0
            : static_cast<double>(result.statistics.reused_leaf_count) /
                  static_cast<double>(new_tree.leaf_count());
    return result;
}

} // 命名空间 cartmesh
