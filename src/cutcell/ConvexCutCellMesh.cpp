#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

constexpr std::array<std::uint8_t, 6> opposite_face = {1, 0, 3, 2, 5, 4};

struct BoxFaceGeometry {
    double area{};
    Vec3 centroid{};
    Vec3 outward_normal{};
    std::vector<Vec3> vertices;
};

struct SolidFacePatch {
    double area{};
    Vec3 centroid{};
    Vec3 outward_normal{};
    std::vector<Vec3> vertices;
};

using SolidFaceOccupancy = std::array<std::vector<SolidFacePatch>, 6>;

[[nodiscard]] std::array<BoxFaceGeometry, 6> box_faces(const AABB& box) noexcept {
    const Vec3 minimum = box.minimum();
    const Vec3 maximum = box.maximum();
    const Vec3 center = box.center();
    const Vec3 extent = box.extent();
    return {{{extent.y * extent.z,
              {minimum.x, center.y, center.z},
              {-1.0, 0.0, 0.0},
              {{minimum.x, minimum.y, minimum.z},
               {minimum.x, minimum.y, maximum.z},
               {minimum.x, maximum.y, maximum.z},
               {minimum.x, maximum.y, minimum.z}}},
             {extent.y * extent.z,
              {maximum.x, center.y, center.z},
              {1.0, 0.0, 0.0},
              {{maximum.x, minimum.y, minimum.z},
               {maximum.x, maximum.y, minimum.z},
               {maximum.x, maximum.y, maximum.z},
               {maximum.x, minimum.y, maximum.z}}},
             {extent.x * extent.z,
              {center.x, minimum.y, center.z},
              {0.0, -1.0, 0.0},
              {{minimum.x, minimum.y, minimum.z},
               {maximum.x, minimum.y, minimum.z},
               {maximum.x, minimum.y, maximum.z},
               {minimum.x, minimum.y, maximum.z}}},
             {extent.x * extent.z,
              {center.x, maximum.y, center.z},
              {0.0, 1.0, 0.0},
              {{minimum.x, maximum.y, minimum.z},
               {minimum.x, maximum.y, maximum.z},
               {maximum.x, maximum.y, maximum.z},
               {maximum.x, maximum.y, minimum.z}}},
             {extent.x * extent.y,
              {center.x, center.y, minimum.z},
              {0.0, 0.0, -1.0},
              {{minimum.x, minimum.y, minimum.z},
               {minimum.x, maximum.y, minimum.z},
               {maximum.x, maximum.y, minimum.z},
               {maximum.x, minimum.y, minimum.z}}},
             {extent.x * extent.y,
              {center.x, center.y, maximum.z},
              {0.0, 0.0, 1.0},
              {{minimum.x, minimum.y, maximum.z},
               {maximum.x, minimum.y, maximum.z},
               {maximum.x, maximum.y, maximum.z},
               {minimum.x, maximum.y, maximum.z}}}}};
}

[[nodiscard]] double default_length_tolerance(const UniformCartesianGrid& grid,
                                              double requested) noexcept {
    const double scale = norm(grid.domain().extent());
    return std::max(requested,
                    256.0 * std::numeric_limits<double>::epsilon() * scale);
}

[[nodiscard]] SolidFaceOccupancy extract_solid_face_occupancy(
    const ConvexSurfaceCutResult& cut) {
    SolidFaceOccupancy result;
    for (std::size_t face_index = 0;
         face_index < cut.solid_polyhedron.faces.size(); ++face_index) {
        const auto& face = cut.solid_polyhedron.faces[face_index];
        if (face.kind != PolyhedronFaceKind::cartesian) {
            continue;
        }
        if (face.source_id >= result.size()) {
            throw std::runtime_error("Cut-cell Cartesian 面编号超出 0..5");
        }
        const auto& geometry = cut.solid_geometry.faces[face_index];
        SolidFacePatch patch{geometry.area, geometry.centroid,
                             geometry.outward_normal, {}};
        patch.vertices.reserve(face.vertex_indices.size());
        for (const auto vertex : face.vertex_indices) {
            patch.vertices.push_back(cut.solid_polyhedron.vertices[vertex]);
        }
        result[static_cast<std::size_t>(face.source_id)].push_back(
            std::move(patch));
    }
    return result;
}

[[nodiscard]] double occupied_area(const std::vector<SolidFacePatch>& patches) {
    double result = 0.0;
    for (const auto& patch : patches) {
        result += patch.area;
    }
    return result;
}

void set_cartesian_aperture(CartesianFaceAperture& aperture,
                            const BoxFaceGeometry& full,
                            const std::vector<SolidFacePatch>& solid_patches,
                            double area_tolerance) {
    const double solid_area = occupied_area(solid_patches);
    Vec3 solid_first_moment{};
    for (const auto& patch : solid_patches) {
        solid_first_moment =
            solid_first_moment + patch.centroid * patch.area;
    }
    double fluid_area = full.area - solid_area;
    if (std::abs(fluid_area) <= area_tolerance) {
        fluid_area = 0.0;
    }
    if (fluid_area < 0.0 || fluid_area > full.area + area_tolerance) {
        throw std::runtime_error("Cut-cell 面开口面积位于 Cartesian 面范围之外");
    }
    fluid_area = std::min(fluid_area, full.area);
    aperture.area = fluid_area;
    aperture.area_fraction = fluid_area / full.area;
    aperture.outward_normal = full.outward_normal;
    aperture.oriented_boundary_loops.clear();
    aperture.oriented_boundary_loops.push_back(full.vertices);
    for (const auto& patch : solid_patches) {
        auto loop = patch.vertices;
        // fluid aperture = full Cartesian face - solid occupancy。来自当前
        // 单元的 solid face 与 full normal 同向，需要反转；共面接触从
        // 相邻固体单元取得时法向已经反向，不再重复反转。
        if (dot(patch.outward_normal, full.outward_normal) > 0.0) {
            std::reverse(loop.begin(), loop.end());
        }
        aperture.oriented_boundary_loops.push_back(std::move(loop));
    }
    aperture.centroid = fluid_area > 0.0
                            ? (full.centroid * full.area - solid_first_moment) /
                                  fluid_area
                            : full.centroid;
}

[[nodiscard]] FluidCellGeometry make_fluid_cell(
    std::uint64_t cell_id, const AABB& box, const ConvexSurfaceCutResult& cut,
    const SolidFaceOccupancy& solid_face_occupancy,
    double length_tolerance) {
    FluidCellGeometry result;
    result.background_cell_id = cell_id;
    result.volume = cut.fluid_volume;
    result.volume_fraction = cut.fluid_volume_fraction;
    result.centroid = cut.fluid_centroid;
    result.cut = cut.cut;
    result.fluid_piece_count = cut.fluid_decomposition_polyhedra.size();
    result.fluid_component_count = cut.fluid_component_count;
    if (cut.cut) {
        result.fluid_polyhedron_pieces.reserve(
            cut.fluid_decomposition_polyhedra.size());
        for (std::size_t piece = 0;
             piece < cut.fluid_decomposition_polyhedra.size(); ++piece) {
            result.fluid_polyhedron_pieces.push_back(
                {cut.fluid_decomposition_polyhedra[piece],
                 cut.fluid_decomposition_geometries[piece],
                 cut.fluid_piece_component_ids[piece]});
        }
    }

    const auto full_faces = box_faces(box);
    for (std::size_t face_index = 0;
         face_index < cut.solid_polyhedron.faces.size(); ++face_index) {
        const auto& face = cut.solid_polyhedron.faces[face_index];
        const auto& geometry = cut.solid_geometry.faces[face_index];
        if (face.kind == PolyhedronFaceKind::embedded_boundary) {
            EmbeddedBoundaryFaceGeometry embedded{
                face.source_id, geometry.area, geometry.centroid,
                geometry.outward_normal * -1.0, {}};
            embedded.vertices.reserve(face.vertex_indices.size());
            for (auto vertex = face.vertex_indices.rbegin();
                 vertex != face.vertex_indices.rend(); ++vertex) {
                embedded.vertices.push_back(cut.solid_polyhedron.vertices[*vertex]);
            }
            result.embedded_boundary_faces.push_back(std::move(embedded));
        }
    }

    const double area_tolerance = length_tolerance * length_tolerance;
    for (std::size_t face = 0; face < result.cartesian_faces.size(); ++face) {
        auto& aperture = result.cartesian_faces[face];
        set_cartesian_aperture(aperture, full_faces[face],
                               solid_face_occupancy[face], area_tolerance);
    }
    return result;
}

void add_coplanar_embedded_boundaries(
    FluidCellGeometry& cell, const std::vector<SolidFacePatch>& target,
    const std::vector<SolidFacePatch>& current, std::uint64_t boundary_id,
    double area_tolerance) {
    if (occupied_area(current) > area_tolerance) {
        return;
    }
    for (const auto& patch : target) {
        EmbeddedBoundaryFaceGeometry embedded{
            boundary_id, patch.area, patch.centroid,
            patch.outward_normal * -1.0, patch.vertices};
        std::reverse(embedded.vertices.begin(), embedded.vertices.end());
        cell.embedded_boundary_faces.push_back(std::move(embedded));
    }
}

void recompute_cell_closure(FluidCellGeometry& cell, double length_tolerance) {
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

[[nodiscard]] bool neighbor_in_positive_direction(const UniformCartesianGrid& grid,
                                                  const CellKey& key,
                                                  std::size_t axis,
                                                  CellKey& neighbor) noexcept {
    neighbor = key;
    if (axis == 0) {
        if (key.i + 1U >= grid.nx()) {
            return false;
        }
        ++neighbor.i;
    } else if (axis == 1) {
        if (key.j + 1U >= grid.ny()) {
            return false;
        }
        ++neighbor.j;
    } else {
        if (key.k + 1U >= grid.nz()) {
            return false;
        }
        ++neighbor.k;
    }
    return true;
}

} // 匿名命名空间

ConvexCutCellMesh build_convex_cut_cell_mesh(const UniformCartesianGrid& grid,
                                             const ConvexSurfaceCutter& cutter,
                                             double geometric_tolerance) {
    if (geometric_tolerance < 0.0 || !std::isfinite(geometric_tolerance)) {
        throw std::invalid_argument("Cut-cell 网格几何容差必须为非负有限数");
    }
    ConvexCutCellMesh result;
    const double length_tolerance =
        default_length_tolerance(grid, geometric_tolerance);
    const double area_tolerance = length_tolerance * length_tolerance;
    const double cell_face_area_scale = std::max(
        {grid.spacing().x * grid.spacing().y,
         grid.spacing().x * grid.spacing().z,
         grid.spacing().y * grid.spacing().z});
    const double closure_tolerance =
        std::max(area_tolerance, 2048.0 * std::numeric_limits<double>::epsilon() *
                                     cell_face_area_scale);

    std::vector<std::size_t> active_index(
        static_cast<std::size_t>(grid.cell_count()),
        std::numeric_limits<std::size_t>::max());
    std::vector<SolidFaceOccupancy> solid_face_occupancies(
        static_cast<std::size_t>(grid.cell_count()));
    result.fluid_cells.reserve(static_cast<std::size_t>(grid.cell_count()));
    for (std::uint64_t cell_id = 0; cell_id < grid.cell_count(); ++cell_id) {
        const CellKey key = grid.cell_key(cell_id);
        const AABB box = grid.cell_bounds(key);
        const auto cut = cutter.cut_box(box);
        auto occupancy = extract_solid_face_occupancy(cut);
        solid_face_occupancies[static_cast<std::size_t>(cell_id)] = occupancy;
        if (cut.fluid_volume <= 0.0) {
            ++result.full_solid_cell_count;
            continue;
        }
        active_index[static_cast<std::size_t>(cell_id)] = result.fluid_cells.size();
        result.fluid_cells.push_back(
            make_fluid_cell(cell_id, box, cut, occupancy, length_tolerance));
        const auto& fluid_cell = result.fluid_cells.back();
        result.total_fluid_volume += fluid_cell.volume;
        if (fluid_cell.fluid_component_count > 1) {
            ++result.disconnected_fluid_cell_count;
        }
        if (fluid_cell.cut) {
            ++result.cut_cell_count;
        } else {
            ++result.full_fluid_cell_count;
        }
    }

    // 当 STL 面与内部 Cartesian 面共面时，零体积接触侧不会自行产生
    // solid intersection。相邻两侧取闭集固体占据的并集，并把新增占据区
    // 作为流体侧的共面 embedded boundary，保证共享面与控制体闭合一致。
    for (std::uint64_t cell_id = 0; cell_id < grid.cell_count(); ++cell_id) {
        const CellKey key = grid.cell_key(cell_id);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CellKey neighbor_key{};
            if (!neighbor_in_positive_direction(grid, key, axis, neighbor_key)) {
                continue;
            }
            const std::uint64_t neighbor_id = grid.linear_id(neighbor_key);
            const std::uint8_t first_face = static_cast<std::uint8_t>(2 * axis + 1);
            const std::uint8_t second_face = opposite_face[first_face];
            const auto& first_solid =
                solid_face_occupancies[static_cast<std::size_t>(cell_id)][first_face];
            const auto& second_solid =
                solid_face_occupancies[static_cast<std::size_t>(neighbor_id)][second_face];
            const auto& target = occupied_area(first_solid) >= occupied_area(second_solid)
                                     ? first_solid
                                     : second_solid;
            const std::size_t first_index =
                active_index[static_cast<std::size_t>(cell_id)];
            const std::size_t second_index =
                active_index[static_cast<std::size_t>(neighbor_id)];
            if (first_index != std::numeric_limits<std::size_t>::max() &&
                occupied_area(target) > occupied_area(first_solid) + area_tolerance) {
                auto& cell = result.fluid_cells[first_index];
                set_cartesian_aperture(
                    cell.cartesian_faces[first_face],
                    box_faces(grid.cell_bounds(key))[first_face], target,
                    area_tolerance);
                add_coplanar_embedded_boundaries(
                    cell, target, first_solid, cutter.boundary_id(), area_tolerance);
            }
            if (second_index != std::numeric_limits<std::size_t>::max() &&
                occupied_area(target) > occupied_area(second_solid) + area_tolerance) {
                auto& cell = result.fluid_cells[second_index];
                set_cartesian_aperture(
                    cell.cartesian_faces[second_face],
                    box_faces(grid.cell_bounds(neighbor_key))[second_face], target,
                    area_tolerance);
                add_coplanar_embedded_boundaries(
                    cell, target, second_solid, cutter.boundary_id(), area_tolerance);
            }
        }
    }

    for (auto& fluid_cell : result.fluid_cells) {
        recompute_cell_closure(fluid_cell, length_tolerance);
        if (!fluid_cell.embedded_boundary_faces.empty()) {
            ++result.boundary_cell_count;
        }
        for (const auto& face : fluid_cell.embedded_boundary_faces) {
            result.total_embedded_boundary_area += face.area;
        }
        result.maximum_cell_area_closure_residual =
            std::max(result.maximum_cell_area_closure_residual,
                     fluid_cell.area_vector_closure_residual);
        if (fluid_cell.area_vector_closure_residual > closure_tolerance ||
            !fluid_cell.boundary_edge_closed) {
            ++result.nonclosed_cell_count;
        }
        if (!(fluid_cell.volume > 0.0)) {
            ++result.negative_volume_cell_count;
        }
    }

    for (std::uint64_t cell_id = 0; cell_id < grid.cell_count(); ++cell_id) {
        const std::size_t first_index = active_index[static_cast<std::size_t>(cell_id)];
        const CellKey key = grid.cell_key(cell_id);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            CellKey neighbor_key{};
            if (!neighbor_in_positive_direction(grid, key, axis, neighbor_key)) {
                continue;
            }
            const std::uint64_t neighbor_id = grid.linear_id(neighbor_key);
            const std::size_t second_index =
                active_index[static_cast<std::size_t>(neighbor_id)];
            const std::uint8_t first_face = static_cast<std::uint8_t>(2 * axis + 1);
            const std::uint8_t second_face = opposite_face[first_face];
            const double first_area = first_index == std::numeric_limits<std::size_t>::max()
                                          ? 0.0
                                          : result.fluid_cells[first_index]
                                                .cartesian_faces[first_face]
                                                .area;
            const double second_area =
                second_index == std::numeric_limits<std::size_t>::max()
                    ? 0.0
                    : result.fluid_cells[second_index]
                          .cartesian_faces[second_face]
                          .area;
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
            const auto& first_aperture =
                result.fluid_cells[first_index].cartesian_faces[first_face];
            const auto& second_aperture =
                result.fluid_cells[second_index].cartesian_faces[second_face];
            const double centroid_mismatch =
                norm(first_aperture.centroid - second_aperture.centroid);
            result.maximum_shared_face_centroid_mismatch =
                std::max(result.maximum_shared_face_centroid_mismatch,
                         centroid_mismatch);
            const AABB first_box = grid.cell_bounds(key);
            const Vec3 reference = box_faces(first_box)[first_face].centroid;
            const double first_moment_mismatch = norm(
                (first_aperture.centroid - reference) * first_area -
                (second_aperture.centroid - reference) * second_area);
            result.maximum_shared_face_first_moment_mismatch = std::max(
                result.maximum_shared_face_first_moment_mismatch,
                first_moment_mismatch);
            const double first_moment_tolerance =
                area_tolerance * std::max({first_box.extent().x,
                                           first_box.extent().y,
                                           first_box.extent().z}) +
                length_tolerance * std::max(first_area, second_area);
            if (first_moment_mismatch > first_moment_tolerance) {
                ++result.shared_face_mismatch_count;
            }
            result.internal_faces.push_back(
                {cell_id,
                 neighbor_id,
                 first_index,
                 second_index,
                 first_face,
                 second_face,
                 0.5 * (first_area + second_area),
                 (first_aperture.centroid + second_aperture.centroid) * 0.5,
                 first_aperture.outward_normal,
                 area_mismatch,
                 centroid_mismatch,
                 first_moment_mismatch});
        }
    }
    assign_global_fluid_regions(result, area_tolerance);
    return result;
}

SmallCutCellStatistics analyze_small_cut_cells(const ConvexCutCellMesh& mesh,
                                               double threshold) {
    if (!std::isfinite(threshold) || threshold < 0.0) {
        throw std::invalid_argument("小 Cut-cell 阈值必须是非负有限数");
    }
    SmallCutCellStatistics result;
    result.threshold = threshold;
    for (const auto& cell : mesh.fluid_cells) {
        if (!cell.cut) continue;
        result.minimum_cut_cell_volume_fraction =
            std::min(result.minimum_cut_cell_volume_fraction,
                     cell.volume_fraction);
        if (cell.volume_fraction >= threshold) continue;
        SmallCutCellEntry entry;
        entry.background_cell_id = cell.background_cell_id;
        entry.volume_fraction = cell.volume_fraction;
        entry.centroid = cell.centroid;
        for (const auto& face : cell.embedded_boundary_faces) {
            entry.boundary_ids.push_back(face.boundary_id);
        }
        std::sort(entry.boundary_ids.begin(), entry.boundary_ids.end());
        entry.boundary_ids.erase(
            std::unique(entry.boundary_ids.begin(), entry.boundary_ids.end()),
            entry.boundary_ids.end());
        result.cells.push_back(std::move(entry));
    }
    return result;
}

} // 命名空间 cartmesh
