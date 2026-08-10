#include "cartmesh/io/CutCellJsonWriter.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cartmesh {
namespace {

void write_vec3(std::ostream& output, const Vec3& value) {
    output << '[' << value.x << ',' << value.y << ',' << value.z << ']';
}

void write_loop(std::ostream& output, const std::vector<Vec3>& loop) {
    output << '[';
    for (std::size_t index = 0; index < loop.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        write_vec3(output, loop[index]);
    }
    output << ']';
}

[[nodiscard]] const char* face_kind_name(PolyhedronFaceKind kind) noexcept {
    switch (kind) {
    case PolyhedronFaceKind::cartesian: return "cartesian";
    case PolyhedronFaceKind::embedded_boundary: return "embedded_boundary";
    case PolyhedronFaceKind::internal_partition: return "internal_partition";
    }
    return "internal_partition";
}

} // 匿名命名空间

void write_cut_cell_geometry_json_impl(
    const std::filesystem::path& path, const AABB& domain,
    const std::string& background_description,
    const ConvexCutCellMesh& mesh, bool solver_ready_cut_cell_mesh) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开 Cut-cell 几何 JSON：" + path.string());
    }
    output << std::setprecision(17)
           << "{\n  \"schema\":\"cartmesh-cutcell-geometry-v1\",\n"
           << "  \"solverReadyCutCellMesh\":"
           << (solver_ready_cut_cell_mesh ? "true" : "false") << ",\n"
           << "  \"domain\":{\"minimum\":";
    write_vec3(output, domain.minimum());
    output << ",\"maximum\":";
    write_vec3(output, domain.maximum());
    output << "},\n" << background_description << "  \"fluidCells\":[\n";
    for (std::size_t cell_index = 0; cell_index < mesh.fluid_cells.size(); ++cell_index) {
        const auto& cell = mesh.fluid_cells[cell_index];
        output << "    {\"backgroundCellId\":" << cell.background_cell_id
               << ",\"volume\":" << cell.volume
               << ",\"volumeFraction\":" << cell.volume_fraction
               << ",\"centroid\":";
        write_vec3(output, cell.centroid);
        output << ",\"cut\":" << (cell.cut ? "true" : "false")
               << ",\"fluidPieceCount\":" << cell.fluid_piece_count
               << ",\"fluidComponentCount\":" << cell.fluid_component_count
               << ",\"fluidComponentRegionIds\":[";
        for (std::size_t component = 0;
             component < cell.fluid_component_region_ids.size(); ++component) {
            if (component != 0) output << ',';
            output << cell.fluid_component_region_ids[component];
        }
        output << ']'
               << ",\"boundaryEdgeClosed\":"
               << (cell.boundary_edge_closed ? "true" : "false")
               << ",\"boundaryEdgeImbalanceCount\":"
               << cell.boundary_edge_imbalance_count
               << ",\"areaVectorClosureResidual\":"
               << cell.area_vector_closure_residual << ",\"cartesianFaces\":[";
        for (std::size_t face_index = 0;
             face_index < cell.cartesian_faces.size(); ++face_index) {
            if (face_index != 0) {
                output << ',';
            }
            const auto& face = cell.cartesian_faces[face_index];
            output << "{\"localFace\":" << face_index << ",\"area\":"
                   << face.area << ",\"areaFraction\":" << face.area_fraction
                   << ",\"centroid\":";
            write_vec3(output, face.centroid);
            output << ",\"outwardNormal\":";
            write_vec3(output, face.outward_normal);
            output << ",\"orientedBoundaryLoops\":[";
            for (std::size_t loop = 0; loop < face.oriented_boundary_loops.size();
                 ++loop) {
                if (loop != 0) {
                    output << ',';
                }
                write_loop(output, face.oriented_boundary_loops[loop]);
            }
            output << "]}";
        }
        output << "],\"embeddedBoundaryFaces\":[";
        for (std::size_t face_index = 0;
             face_index < cell.embedded_boundary_faces.size(); ++face_index) {
            if (face_index != 0) {
                output << ',';
            }
            const auto& face = cell.embedded_boundary_faces[face_index];
            output << "{\"boundaryId\":" << face.boundary_id << ",\"area\":"
                   << face.area << ",\"centroid\":";
            write_vec3(output, face.centroid);
            output << ",\"outwardNormal\":";
            write_vec3(output, face.outward_normal);
            output << ",\"vertices\":";
            write_loop(output, face.vertices);
            output << '}';
        }
        output << "],\"fluidPolyhedronPieces\":[";
        for (std::size_t piece_index = 0;
             piece_index < cell.fluid_polyhedron_pieces.size(); ++piece_index) {
            if (piece_index != 0) output << ',';
            const auto& piece = cell.fluid_polyhedron_pieces[piece_index];
            output << "{\"componentId\":" << piece.component_id
                   << ",\"globalRegionId\":" << piece.global_region_id
                   << ",\"volume\":" << piece.geometry.volume
                   << ",\"centroid\":";
            write_vec3(output, piece.geometry.centroid);
            output << ",\"closed\":" << (piece.geometry.closed ? "true" : "false")
                   << ",\"positiveVolume\":"
                   << (piece.geometry.positive_volume ? "true" : "false")
                   << ",\"vertices\":";
            write_loop(output, piece.polyhedron.vertices);
            output << ",\"faces\":[";
            for (std::size_t face_index = 0;
                 face_index < piece.polyhedron.faces.size(); ++face_index) {
                if (face_index != 0) output << ',';
                const auto& face = piece.polyhedron.faces[face_index];
                const auto& geometry = piece.geometry.faces[face_index];
                output << "{\"kind\":\"" << face_kind_name(face.kind)
                       << "\",\"sourceId\":" << face.source_id
                       << ",\"vertexIndices\":[";
                for (std::size_t vertex = 0; vertex < face.vertex_indices.size();
                     ++vertex) {
                    if (vertex != 0) output << ',';
                    output << face.vertex_indices[vertex];
                }
                output << "],\"area\":" << geometry.area
                       << ",\"centroid\":";
                write_vec3(output, geometry.centroid);
                output << ",\"outwardNormal\":";
                write_vec3(output, geometry.outward_normal);
                output << '}';
            }
            output << "]}";
        }
        output << "]}" << (cell_index + 1U == mesh.fluid_cells.size() ? '\n' : ',')
               << '\n';
    }
    output << "  ],\n  \"internalFaces\":[\n";
    for (std::size_t face_index = 0; face_index < mesh.internal_faces.size();
         ++face_index) {
        const auto& face = mesh.internal_faces[face_index];
        output << "    {\"firstBackgroundCellId\":"
               << face.first_background_cell_id << ",\"secondBackgroundCellId\":"
               << face.second_background_cell_id << ",\"firstFluidCellIndex\":"
               << face.first_fluid_cell_index << ",\"secondFluidCellIndex\":"
               << face.second_fluid_cell_index << ",\"firstLocalFace\":"
               << static_cast<unsigned>(face.first_local_face)
               << ",\"secondLocalFace\":"
               << static_cast<unsigned>(face.second_local_face) << ",\"area\":"
               << face.area << ",\"centroid\":";
        write_vec3(output, face.centroid);
        output << ",\"normal\":";
        write_vec3(output, face.normal);
        output << ",\"areaMismatch\":" << face.area_mismatch
               << ",\"centroidMismatch\":" << face.centroid_mismatch
               << ",\"firstMomentMismatch\":"
               << face.first_moment_mismatch << '}'
               << (face_index + 1U == mesh.internal_faces.size() ? '\n' : ',')
               << '\n';
    }
    output << "  ],\n  \"componentInternalFaces\":[\n";
    for (std::size_t face_index = 0;
         face_index < mesh.component_internal_faces.size(); ++face_index) {
        const auto& face = mesh.component_internal_faces[face_index];
        output << "    {\"firstBackgroundCellId\":"
               << face.first_background_cell_id
               << ",\"secondBackgroundCellId\":"
               << face.second_background_cell_id
               << ",\"firstFluidCellIndex\":" << face.first_fluid_cell_index
               << ",\"secondFluidCellIndex\":" << face.second_fluid_cell_index
               << ",\"firstComponentId\":" << face.first_component_id
               << ",\"secondComponentId\":" << face.second_component_id
               << ",\"globalRegionId\":" << face.global_region_id
               << ",\"area\":" << face.area << ",\"centroid\":";
        write_vec3(output, face.centroid);
        output << ",\"normal\":";
        write_vec3(output, face.normal);
        output << '}'
               << (face_index + 1U == mesh.component_internal_faces.size()
                       ? '\n'
                       : ',')
               << '\n';
    }
    output << "  ],\n  \"globalFluidRegions\":[";
    for (std::size_t region = 0;
         region < mesh.global_fluid_region_volumes.size(); ++region) {
        if (region != 0) output << ',';
        output << "{\"regionId\":" << region << ",\"volume\":"
               << mesh.global_fluid_region_volumes[region] << '}';
    }
    output << "]\n}\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入 Cut-cell 几何 JSON 失败：" + path.string());
    }
}

void write_cut_cell_geometry_json(const std::filesystem::path& path,
                                  const UniformCartesianGrid& grid,
                                  const ConvexCutCellMesh& mesh,
                                  bool solver_ready_cut_cell_mesh) {
    std::ostringstream background;
    background << "  \"meshKind\":\"uniform_cartesian\",\n"
               << "  \"dimensions\":[" << grid.nx() << ',' << grid.ny()
               << ',' << grid.nz() << "],\n";
    write_cut_cell_geometry_json_impl(path, grid.domain(), background.str(),
                                      mesh, solver_ready_cut_cell_mesh);
}

void write_cut_cell_geometry_json(const std::filesystem::path& path,
                                  const LinearOctree& tree,
                                  const ConvexCutCellMesh& mesh,
                                  bool solver_ready_cut_cell_mesh) {
    std::ostringstream background;
    background << "  \"meshKind\":\"adaptive_linear_octree\",\n"
               << "  \"baseLevel\":"
               << static_cast<unsigned>(tree.base_level()) << ",\n"
               << "  \"maximumLevel\":"
               << static_cast<unsigned>(tree.maximum_level()) << ",\n"
               << "  \"leafNodeCodes\":[";
    for (std::uint64_t leaf = 0; leaf < tree.leaf_count(); ++leaf) {
        if (leaf != 0) background << ',';
        background << tree.leaf_code(leaf);
    }
    background << "],\n";
    write_cut_cell_geometry_json_impl(path, tree.domain(), background.str(),
                                      mesh, solver_ready_cut_cell_mesh);
}

} // 命名空间 cartmesh
