#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"

#include "cartmesh/geometry/SurfaceDiagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cartmesh {
namespace {

[[nodiscard]] bool vertex_less(const Vec3& first, const Vec3& second) noexcept {
    if (first.x != second.x) return first.x < second.x;
    if (first.y != second.y) return first.y < second.y;
    return first.z < second.z;
}

[[nodiscard]] std::array<Vec3, 3> canonical_vertices(
    const Triangle& triangle) noexcept {
    const auto& vertices = triangle.vertices();
    std::size_t first = 0;
    if (vertex_less(vertices[1], vertices[first])) first = 1;
    if (vertex_less(vertices[2], vertices[first])) first = 2;
    return {vertices[first], vertices[(first + 1U) % 3U],
            vertices[(first + 2U) % 3U]};
}

[[nodiscard]] bool triangle_less(const Triangle& first,
                                 const Triangle& second) noexcept {
    const auto first_vertices = canonical_vertices(first);
    const auto second_vertices = canonical_vertices(second);
    for (std::size_t vertex = 0; vertex < 3; ++vertex) {
        if (vertex_less(first_vertices[vertex], second_vertices[vertex])) return true;
        if (vertex_less(second_vertices[vertex], first_vertices[vertex])) return false;
    }
    return false;
}

[[nodiscard]] SurfaceMesh make_oriented_surface(
    const SurfaceMesh& surface,
    std::vector<std::uint64_t>& triangle_boundary_ids) {
    if (triangle_boundary_ids.size() != surface.triangles().size()) {
        throw std::invalid_argument(
            "三角片 boundary ID 数量必须与 STL 三角片数一致");
    }
    const auto diagnostics = diagnose_surface(surface);
    if (!diagnostics.valid_for_stage1_classification()) {
        throw std::invalid_argument(
            "通用 Cut-cell 需要封闭、定向一致且流形的三角表面");
    }
    const bool reverse = diagnostics.component_orientation_mismatch_count ==
                         diagnostics.connected_component_count;
    if (diagnostics.component_orientation_mismatch_count != 0 && !reverse) {
        throw std::invalid_argument(
            "多壳层 Cut-cell 要求所有分量符合嵌套方向，或仅整体统一反向");
    }
    struct TriangleWithBoundary {
        Triangle triangle;
        std::uint64_t boundary_id{};
    };
    std::vector<TriangleWithBoundary> records;
    records.reserve(surface.triangles().size());
    for (std::size_t index = 0; index < surface.triangles().size(); ++index) {
        const auto& triangle = surface.triangles()[index];
        const auto& vertex = triangle.vertices();
        records.push_back(
            {Triangle(vertex[0], vertex[reverse ? 2 : 1],
                      vertex[reverse ? 1 : 2]),
             triangle_boundary_ids[index]});
    }
    std::stable_sort(records.begin(), records.end(),
                     [](const auto& first, const auto& second) {
                         return triangle_less(first.triangle, second.triangle);
                     });
    std::vector<Triangle> triangles;
    triangles.reserve(records.size());
    for (std::size_t index = 0; index < records.size(); ++index) {
        const auto vertices = canonical_vertices(records[index].triangle);
        triangles.emplace_back(vertices[0], vertices[1], vertices[2]);
        triangle_boundary_ids[index] = records[index].boundary_id;
    }
    return SurfaceMesh(std::move(triangles), surface.format(), surface.name());
}

[[nodiscard]] AABB tetrahedron_bounds(const Vec3& a, const Vec3& b,
                                      const Vec3& c, const Vec3& d) {
    const Vec3 minimum{std::min({a.x, b.x, c.x, d.x}),
                       std::min({a.y, b.y, c.y, d.y}),
                       std::min({a.z, b.z, c.z, d.z})};
    const Vec3 maximum{std::max({a.x, b.x, c.x, d.x}),
                       std::max({a.y, b.y, c.y, d.y}),
                       std::max({a.z, b.z, c.z, d.z})};
    return AABB(minimum, maximum);
}

[[nodiscard]] std::array<Vec3, 6> box_face_normals() noexcept {
    return {{{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
             {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0},
             {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}};
}

[[nodiscard]] AABB expanded_box(const AABB& box, double tolerance) {
    const Vec3 padding{tolerance, tolerance, tolerance};
    return AABB(box.minimum() - padding, box.maximum() + padding);
}

[[nodiscard]] std::array<double, 6> box_face_coordinates(const AABB& box) noexcept {
    return {box.minimum().x, box.maximum().x, box.minimum().y,
            box.maximum().y, box.minimum().z, box.maximum().z};
}

[[nodiscard]] std::size_t face_axis(std::size_t face) noexcept { return face / 2U; }

[[nodiscard]] double component(const Vec3& value, std::size_t axis) noexcept {
    return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

[[nodiscard]] bool polygon_on_box_face(const std::vector<Vec3>& vertices,
                                       std::size_t face, const AABB& box,
                                       double tolerance) noexcept {
    const auto coordinates = box_face_coordinates(box);
    const std::size_t axis = face_axis(face);
    return !vertices.empty() &&
           std::all_of(vertices.begin(), vertices.end(), [&](const Vec3& vertex) {
               return std::abs(component(vertex, axis) - coordinates[face]) <=
                      tolerance;
           });
}

[[nodiscard]] PolygonGeometry measure_polygon(const std::vector<Vec3>& vertices) {
    PolygonGeometry result;
    if (vertices.size() < 3) {
        return result;
    }
    const Vec3 reference = vertices.front();
    Vec3 weighted_centroid{};
    for (std::size_t index = 1; index + 1 < vertices.size(); ++index) {
        const Vec3 area_vector =
            cross(vertices[index] - reference,
                  vertices[index + 1] - reference) *
            0.5;
        const double area = norm(area_vector);
        result.area_vector = result.area_vector + area_vector;
        result.area += area;
        weighted_centroid =
            weighted_centroid +
            (reference + vertices[index] + vertices[index + 1]) *
                (area / 3.0);
    }
    if (result.area > 0.0) {
        result.centroid = weighted_centroid / result.area;
    }
    const double vector_area = norm(result.area_vector);
    if (vector_area > 0.0) {
        result.outward_normal = result.area_vector / vector_area;
    }
    return result;
}

[[nodiscard]] std::vector<Vec3> clip_polygon_plane(
    const std::vector<Vec3>& input, std::size_t axis, double coordinate,
    bool keep_greater, double tolerance) {
    std::vector<Vec3> output;
    if (input.empty()) {
        return output;
    }
    for (std::size_t edge = 0; edge < input.size(); ++edge) {
        const Vec3 first = input[edge];
        const Vec3 second = input[(edge + 1U) % input.size()];
        const double first_distance =
            (component(first, axis) - coordinate) * (keep_greater ? 1.0 : -1.0);
        const double second_distance =
            (component(second, axis) - coordinate) * (keep_greater ? 1.0 : -1.0);
        const bool first_inside = first_distance >= -tolerance;
        const bool second_inside = second_distance >= -tolerance;
        if (first_inside) {
            output.push_back(first);
        }
        if (first_inside != second_inside) {
            const double fraction = first_distance /
                                    (first_distance - second_distance);
            output.push_back(first + (second - first) * fraction);
        }
    }
    std::vector<Vec3> deduplicated;
    for (const auto& vertex : output) {
        if (deduplicated.empty() ||
            norm(deduplicated.back() - vertex) > tolerance) {
            deduplicated.push_back(vertex);
        }
    }
    if (deduplicated.size() > 1 &&
        norm(deduplicated.front() - deduplicated.back()) <= tolerance) {
        deduplicated.pop_back();
    }
    return deduplicated;
}

[[nodiscard]] std::vector<Vec3> clip_triangle_to_box(const Triangle& triangle,
                                                     const AABB& box,
                                                     double tolerance) {
    std::vector<Vec3> polygon(triangle.vertices().begin(),
                              triangle.vertices().end());
    polygon = clip_polygon_plane(polygon, 0, box.minimum().x, true, tolerance);
    polygon = clip_polygon_plane(polygon, 0, box.maximum().x, false, tolerance);
    polygon = clip_polygon_plane(polygon, 1, box.minimum().y, true, tolerance);
    polygon = clip_polygon_plane(polygon, 1, box.maximum().y, false, tolerance);
    polygon = clip_polygon_plane(polygon, 2, box.minimum().z, true, tolerance);
    polygon = clip_polygon_plane(polygon, 2, box.maximum().z, false, tolerance);
    return polygon;
}

[[nodiscard]] ConvexPolyhedron clip_tetrahedron_to_box(
    ConvexPolyhedron tetrahedron, const AABB& box, double tolerance) {
    const std::array<OrientedHalfSpace, 6> half_spaces = {
        OrientedHalfSpace({box.minimum().x, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 0),
        OrientedHalfSpace({box.maximum().x, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1),
        OrientedHalfSpace({0.0, box.minimum().y, 0.0}, {0.0, -1.0, 0.0}, 2),
        OrientedHalfSpace({0.0, box.maximum().y, 0.0}, {0.0, 1.0, 0.0}, 3),
        OrientedHalfSpace({0.0, 0.0, box.minimum().z}, {0.0, 0.0, -1.0}, 4),
        OrientedHalfSpace({0.0, 0.0, box.maximum().z}, {0.0, 0.0, 1.0}, 5)};
    for (const auto& half_space : half_spaces) {
        tetrahedron = clip_convex_polyhedron(
            tetrahedron, half_space, tolerance,
            PolyhedronFaceKind::cartesian);
        if (tetrahedron.empty()) {
            break;
        }
    }
    return tetrahedron;
}

} // 匿名命名空间

TriangulatedSurfaceCutter::TriangulatedSurfaceCutter(
    const SurfaceMesh& surface, std::uint64_t boundary_id,
    double length_tolerance)
    : TriangulatedSurfaceCutter(
          surface,
          std::vector<std::uint64_t>(surface.triangles().size(), boundary_id),
          length_tolerance) {
    boundary_id_ = boundary_id;
}

TriangulatedSurfaceCutter::TriangulatedSurfaceCutter(
    const SurfaceMesh& surface,
    std::vector<std::uint64_t> triangle_boundary_ids,
    double length_tolerance)
    : oriented_surface_(make_oriented_surface(surface, triangle_boundary_ids)),
      bvh_(oriented_surface_),
      triangle_boundary_ids_(std::move(triangle_boundary_ids)),
      reference_(oriented_surface_.bounds().center()) {
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("通用 Cut-cell 长度容差必须为非负有限数");
    }
    const auto diagnostics = diagnose_surface(oriented_surface_);
    if (triangle_boundary_ids_.size() != oriented_surface_.triangles().size()) {
        throw std::invalid_argument("三角片 boundary ID 数量必须与 STL 三角片数一致");
    }
    if (!triangle_boundary_ids_.empty()) {
        boundary_id_ = triangle_boundary_ids_.front();
    }
    length_tolerance_ =
        std::max(length_tolerance, diagnostics.suggested_length_tolerance);
    tetrahedra_.reserve(oriented_surface_.triangles().size());
    for (const auto& triangle : oriented_surface_.triangles()) {
        const auto& vertex = triangle.vertices();
        const double determinant =
            dot(vertex[0] - reference_,
                cross(vertex[1] - reference_, vertex[2] - reference_));
        if (determinant == 0.0) {
            continue;
        }
        tetrahedra_.push_back(
            {make_tetrahedron_polyhedron(reference_, vertex[0], vertex[1],
                                         vertex[2]),
             tetrahedron_bounds(reference_, vertex[0], vertex[1], vertex[2]),
             determinant > 0.0 ? 1.0 : -1.0});
    }
}

TriangulatedSurfaceBoundaryClip
TriangulatedSurfaceCutter::clip_boundary_faces(const AABB& box) const {
    if (!box.has_positive_volume()) {
        throw std::invalid_argument("STL wall 裁剪的背景盒必须具有正体积");
    }
    TriangulatedSurfaceBoundaryClip result;
    const Vec3 cell_centroid = box.center();
    const AABB local_box(box.minimum() - cell_centroid,
                         box.maximum() - cell_centroid);
    const double local_scale = norm(box.extent());
    const double tolerance = std::max(
        length_tolerance_,
        256.0 * std::numeric_limits<double>::epsilon() * local_scale);
    const double area_tolerance = tolerance * tolerance;
    const auto face_normals = box_face_normals();

    for (const auto triangle_id : bvh_.query(expanded_box(box, tolerance))) {
        const auto& source_triangle = oriented_surface_.triangles()[
            static_cast<std::size_t>(triangle_id)];
        const auto& source_vertex = source_triangle.vertices();
        const Triangle triangle(source_vertex[0] - cell_centroid,
                                source_vertex[1] - cell_centroid,
                                source_vertex[2] - cell_centroid);
        auto polygon = clip_triangle_to_box(triangle, local_box, tolerance);
        const auto geometry = measure_polygon(polygon);
        if (geometry.area <= area_tolerance) continue;

        bool coplanar_with_box = false;
        bool fluid_side_inside_cell = true;
        for (std::size_t face = 0; face < 6; ++face) {
            if (!polygon_on_box_face(polygon, face, local_box, tolerance)) {
                continue;
            }
            coplanar_with_box = true;
            fluid_side_inside_cell =
                dot(geometry.outward_normal, face_normals[face]) < 0.0;
            break;
        }
        if (coplanar_with_box && !fluid_side_inside_cell) continue;

        EmbeddedBoundaryFaceGeometry embedded{
            triangle_boundary_ids_[static_cast<std::size_t>(triangle_id)],
            geometry.area,
            geometry.centroid + cell_centroid,
            geometry.outward_normal * -1.0,
            polygon};
        std::reverse(embedded.vertices.begin(), embedded.vertices.end());
        for (auto& vertex : embedded.vertices) {
            vertex = vertex + cell_centroid;
        }
        result.total_area += embedded.area;
        result.faces.push_back(std::move(embedded));
    }
    return result;
}

TriangulatedSurfaceCellCut TriangulatedSurfaceCutter::cut_box(
    const AABB& box) const {
    if (!box.has_positive_volume()) {
        throw std::invalid_argument("通用 Cut-cell 背景盒必须具有正体积");
    }
    TriangulatedSurfaceCellCut result;
    const double cell_volume = box.volume();
    const Vec3 cell_centroid = box.center();
    const AABB local_box(box.minimum() - cell_centroid,
                         box.maximum() - cell_centroid);
    const double local_scale = norm(box.extent());
    const double tolerance = std::max(
        length_tolerance_,
        256.0 * std::numeric_limits<double>::epsilon() * local_scale);
    const double area_tolerance = tolerance * tolerance;
    const auto face_normals = box_face_normals();
    std::array<Vec3, 6> solid_face_first_moment{};
    Vec3 relative_solid_first_moment{};

    for (const auto& contribution : tetrahedra_) {
        if (!contribution.bounds.intersects(box)) {
            continue;
        }
        auto local_tetrahedron = contribution.tetrahedron;
        for (auto& vertex : local_tetrahedron.vertices) {
            vertex = vertex - cell_centroid;
        }
        const auto clipped = clip_tetrahedron_to_box(
            std::move(local_tetrahedron), local_box, tolerance);
        const auto geometry = measure_polyhedron(clipped);
        if (!geometry.positive_volume) {
            continue;
        }
        const double signed_volume = contribution.coefficient * geometry.volume;
        result.solid_volume += signed_volume;
        relative_solid_first_moment =
            relative_solid_first_moment +
            geometry.centroid * signed_volume;
        for (std::size_t polyhedron_face = 0;
             polyhedron_face < clipped.faces.size(); ++polyhedron_face) {
            const auto& face = clipped.faces[polyhedron_face];
            std::vector<Vec3> vertices;
            vertices.reserve(face.vertex_indices.size());
            for (const auto vertex : face.vertex_indices) {
                vertices.push_back(clipped.vertices[vertex]);
            }
            for (std::size_t box_face = 0; box_face < 6; ++box_face) {
                if (!polygon_on_box_face(vertices, box_face, local_box,
                                         tolerance)) {
                    continue;
                }
                const double signed_area = contribution.coefficient *
                    dot(geometry.faces[polyhedron_face].area_vector,
                        face_normals[box_face]);
                if (std::abs(signed_area) <= area_tolerance) {
                    break;
                }
                result.solid_cartesian_faces[box_face].area += signed_area;
                solid_face_first_moment[box_face] =
                    solid_face_first_moment[box_face] +
                    geometry.faces[polyhedron_face].centroid * signed_area;
                if (signed_area < 0.0) {
                    std::reverse(vertices.begin(), vertices.end());
                }
                result.solid_cartesian_faces[box_face]
                    .oriented_boundary_loops.push_back(std::move(vertices));
                break;
            }
        }
    }

    std::array<double, 6> coplanar_surface_area{};
    std::array<Vec3, 6> coplanar_surface_first_moment{};
    std::array<std::vector<std::vector<Vec3>>, 6> coplanar_surface_loops{};
    for (const auto triangle_id : bvh_.query(expanded_box(box, tolerance))) {
        const auto& source_triangle = oriented_surface_.triangles()[
            static_cast<std::size_t>(triangle_id)];
        const auto& source_vertex = source_triangle.vertices();
        const Triangle triangle(source_vertex[0] - cell_centroid,
                                source_vertex[1] - cell_centroid,
                                source_vertex[2] - cell_centroid);
        auto polygon = clip_triangle_to_box(triangle, local_box, tolerance);
        const auto geometry = measure_polygon(polygon);
        if (geometry.area <= area_tolerance) {
            continue;
        }
        bool coplanar_with_box = false;
        bool fluid_side_inside_cell = true;
        for (std::size_t face = 0; face < 6; ++face) {
            if (!polygon_on_box_face(polygon, face, local_box, tolerance)) {
                continue;
            }
            coplanar_with_box = true;
            fluid_side_inside_cell =
                dot(geometry.outward_normal, face_normals[face]) < 0.0;
            coplanar_surface_area[face] += geometry.area;
            coplanar_surface_first_moment[face] =
                coplanar_surface_first_moment[face] +
                geometry.centroid * geometry.area;
            auto solid_loop = polygon;
            if (dot(measure_polygon(solid_loop).area_vector,
                    face_normals[face]) < 0.0) {
                std::reverse(solid_loop.begin(), solid_loop.end());
            }
            coplanar_surface_loops[face].push_back(std::move(solid_loop));
            break;
        }
        if (!coplanar_with_box || fluid_side_inside_cell) {
            EmbeddedBoundaryFaceGeometry embedded{
                triangle_boundary_ids_[static_cast<std::size_t>(triangle_id)],
                geometry.area, geometry.centroid + cell_centroid,
                geometry.outward_normal * -1.0, polygon};
            std::reverse(embedded.vertices.begin(), embedded.vertices.end());
            for (auto& vertex : embedded.vertices) {
                vertex = vertex + cell_centroid;
            }
            result.embedded_boundary_area += embedded.area;
            result.embedded_boundary_faces.push_back(std::move(embedded));
        }
    }

    const double fraction_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    if (result.solid_volume / cell_volume <= fraction_tolerance) {
        result.solid_volume = 0.0;
        relative_solid_first_moment = {};
    } else if ((cell_volume - result.solid_volume) / cell_volume <=
               fraction_tolerance) {
        result.solid_volume = cell_volume;
        relative_solid_first_moment = {};
    }
    if (result.solid_volume < 0.0 || result.solid_volume > cell_volume) {
        throw std::runtime_error("有向四面体链给出的单元固体体积超出背景盒");
    }
    if (result.solid_volume > 0.0) {
        result.solid_centroid =
            cell_centroid + relative_solid_first_moment / result.solid_volume;
    }
    result.fluid_volume = cell_volume - result.solid_volume;
    result.fluid_volume_fraction = result.fluid_volume / cell_volume;
    if (result.fluid_volume > 0.0) {
        result.fluid_centroid =
            cell_centroid - relative_solid_first_moment / result.fluid_volume;
    }
    result.volume_conservation_residual =
        result.solid_volume + result.fluid_volume - cell_volume;
    result.cut = result.solid_volume > fraction_tolerance * cell_volume &&
                 result.fluid_volume > fraction_tolerance * cell_volume;

    for (std::size_t face = 0; face < 6; ++face) {
        auto& region = result.solid_cartesian_faces[face];
        if (coplanar_surface_area[face] > region.area + area_tolerance) {
            region.area = coplanar_surface_area[face];
            solid_face_first_moment[face] =
                coplanar_surface_first_moment[face];
            region.oriented_boundary_loops =
                std::move(coplanar_surface_loops[face]);
        }
        if (std::abs(region.area) <= area_tolerance) {
            region.area = 0.0;
        }
        if (region.area > 0.0) {
            region.centroid = cell_centroid +
                              solid_face_first_moment[face] / region.area;
        }
        for (auto& loop : region.oriented_boundary_loops) {
            for (auto& vertex : loop) vertex = vertex + cell_centroid;
        }
    }
    return result;
}

} // 命名空间 cartmesh
