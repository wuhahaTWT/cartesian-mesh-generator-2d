#include "cartmesh/io/CutCellVtkWriter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace cartmesh {
namespace {

struct ProjectedVertex {
    double u{};
    double v{};
    std::uint32_t id{};
};

[[nodiscard]] std::vector<std::uint32_t> convex_face_vertices(
    const ConvexPolyhedron& polyhedron, const PolyhedronFace& face,
    const Vec3& outward_normal) {
    if (face.vertex_indices.size() <= 3U) return face.vertex_indices;
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& first, const Vec3& second) {
            return std::abs(dot(first, outward_normal)) <
                   std::abs(dot(second, outward_normal));
        });
    const Vec3 tangent_u = cross(reference_axis, outward_normal) /
                           norm(cross(reference_axis, outward_normal));
    const Vec3 tangent_v = cross(outward_normal, tangent_u);
    const Vec3 origin = polyhedron.vertices[face.vertex_indices.front()];
    std::vector<ProjectedVertex> points;
    points.reserve(face.vertex_indices.size());
    double scale = 0.0;
    for (const auto id : face.vertex_indices) {
        const Vec3 delta = polyhedron.vertices[id] - origin;
        const double u = dot(delta, tangent_u);
        const double v = dot(delta, tangent_v);
        points.push_back({u, v, id});
        scale = std::max({scale, std::abs(u), std::abs(v)});
    }
    std::sort(points.begin(), points.end(), [](const auto& first,
                                               const auto& second) {
        if (first.u != second.u) return first.u < second.u;
        if (first.v != second.v) return first.v < second.v;
        return first.id < second.id;
    });
    points.erase(std::unique(points.begin(), points.end(), [](const auto& first,
                                                              const auto& second) {
        return first.id == second.id;
    }), points.end());
    if (points.size() <= 3U) {
        std::vector<std::uint32_t> result;
        for (const auto& point : points) result.push_back(point.id);
        return result;
    }
    const double length_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(scale, std::numeric_limits<double>::min());
    const double cross_tolerance = length_tolerance * std::max(scale, length_tolerance);
    const auto turn = [](const ProjectedVertex& first,
                         const ProjectedVertex& second,
                         const ProjectedVertex& third) {
        return (second.u - first.u) * (third.v - first.v) -
               (second.v - first.v) * (third.u - first.u);
    };
    std::vector<ProjectedVertex> corners;
    corners.reserve(points.size() * 2U);
    for (const auto& point : points) {
        while (corners.size() >= 2U &&
               turn(corners[corners.size() - 2U], corners.back(), point) <=
                   cross_tolerance) {
            corners.pop_back();
        }
        corners.push_back(point);
    }
    const std::size_t lower_size = corners.size();
    for (std::size_t index = points.size() - 1U; index-- > 0U;) {
        const auto& point = points[index];
        while (corners.size() > lower_size &&
               turn(corners[corners.size() - 2U], corners.back(), point) <=
                   cross_tolerance) {
            corners.pop_back();
        }
        corners.push_back(point);
    }
    corners.pop_back();
    if (corners.size() < 3U) return face.vertex_indices;

    std::vector<std::uint32_t> result;
    result.reserve(points.size());
    for (std::size_t edge = 0; edge < corners.size(); ++edge) {
        const auto& first = corners[edge];
        const auto& second = corners[(edge + 1U) % corners.size()];
        const double du = second.u - first.u;
        const double dv = second.v - first.v;
        const double squared_length = du * du + dv * dv;
        const double edge_length = std::sqrt(squared_length);
        std::vector<std::pair<double, std::uint32_t>> on_edge;
        for (const auto& point : points) {
            const double local_u = point.u - first.u;
            const double local_v = point.v - first.v;
            const double parameter =
                (local_u * du + local_v * dv) / squared_length;
            const double distance =
                std::abs(du * local_v - dv * local_u) / edge_length;
            if (parameter >= -length_tolerance / edge_length &&
                parameter < 1.0 - length_tolerance / edge_length &&
                distance <= length_tolerance) {
                on_edge.emplace_back(parameter, point.id);
            }
        }
        std::sort(on_edge.begin(), on_edge.end());
        for (const auto [parameter, id] : on_edge) {
            static_cast<void>(parameter);
            if (result.empty() || result.back() != id) result.push_back(id);
        }
    }
    if (result.size() < 3U) return face.vertex_indices;
    Vec3 area_vector{};
    const Vec3 first_point = polyhedron.vertices[result.front()];
    for (std::size_t index = 1; index + 1U < result.size(); ++index) {
        area_vector = area_vector +
                      cross(polyhedron.vertices[result[index]] - first_point,
                            polyhedron.vertices[result[index + 1U]] - first_point);
    }
    if (dot(area_vector, outward_normal) < 0.0) {
        std::reverse(result.begin(), result.end());
    }
    return result;
}

} // 匿名命名空间

void write_embedded_boundary_vtp(const std::filesystem::path& path,
                                 const ConvexCutCellMesh& mesh) {
    std::uint64_t point_count = 0;
    std::uint64_t polygon_count = 0;
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            point_count += static_cast<std::uint64_t>(face.vertices.size());
            ++polygon_count;
        }
    }
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开 Cut-cell VTP 输出文件：" + path.string());
    }
    output << std::setprecision(17)
           << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"PolyData\" version=\"1.0\" "
              "byte_order=\"LittleEndian\">\n"
           << "  <PolyData>\n"
           << "    <Piece NumberOfPoints=\"" << point_count
           << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" "
              "NumberOfPolys=\""
           << polygon_count << "\">\n"
           << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            for (const auto& vertex : face.vertices) {
                output << vertex.x << ' ' << vertex.y << ' ' << vertex.z << '\n';
            }
        }
    }
    output << "        </DataArray>\n"
           << "      </Points>\n"
           << "      <Polys>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
              "format=\"ascii\">\n";
    std::uint64_t point_id = 0;
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex) {
                output << point_id++
                       << (vertex + 1U == face.vertices.size() ? '\n' : ' ');
            }
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" "
              "format=\"ascii\">\n";
    std::uint64_t offset = 0;
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            offset += static_cast<std::uint64_t>(face.vertices.size());
            output << offset << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "      </Polys>\n"
           << "      <CellData>\n"
           << "        <DataArray type=\"UInt64\" Name=\"background_cell_id\" "
              "format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (std::size_t face = 0; face < cell.embedded_boundary_faces.size(); ++face) {
            output << cell.background_cell_id << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt64\" Name=\"boundary_id\" "
              "format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            output << face.boundary_id << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Float64\" Name=\"area\" "
              "format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            output << face.area << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Float64\" Name=\"fluid_outward_normal\" "
              "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& face : cell.embedded_boundary_faces) {
            output << face.outward_normal.x << ' ' << face.outward_normal.y << ' '
                   << face.outward_normal.z << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "      </CellData>\n"
           << "    </Piece>\n"
           << "  </PolyData>\n"
           << "</VTKFile>\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入 Cut-cell VTP 输出文件失败：" + path.string());
    }
}

void write_fluid_polyhedra_vtu(const std::filesystem::path& path,
                               const ConvexCutCellMesh& mesh) {
    struct OutputPiece {
        const FluidCellGeometry* cell{};
        const FluidPolyhedronPiece* piece{};
        std::vector<std::size_t> used_vertices;
        std::vector<std::size_t> compact_vertex;
        std::vector<std::vector<std::uint32_t>> faces;
    };
    std::vector<OutputPiece> pieces;
    std::uint64_t point_count = 0;
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            if (!piece.geometry.closed || !piece.geometry.positive_volume) {
                throw std::runtime_error(
                    "不得写出非闭合或非正体积的 Cut-cell 多面体片：backgroundCell=" +
                    std::to_string(cell.background_cell_id) +
                    " closed=" + std::to_string(piece.geometry.closed) +
                    " positiveVolume=" +
                    std::to_string(piece.geometry.positive_volume) +
                    " volume=" + std::to_string(piece.geometry.volume) +
                    " vertices=" +
                    std::to_string(piece.polyhedron.vertices.size()) +
                    " faces=" + std::to_string(piece.polyhedron.faces.size()));
            }
            OutputPiece output_piece;
            output_piece.cell = &cell;
            output_piece.piece = &piece;
            output_piece.compact_vertex.assign(
                piece.polyhedron.vertices.size(),
                std::numeric_limits<std::size_t>::max());
            std::vector<bool> used(piece.polyhedron.vertices.size(), false);
            for (std::size_t face_index = 0;
                 face_index < piece.polyhedron.faces.size(); ++face_index) {
                output_piece.faces.push_back(convex_face_vertices(
                    piece.polyhedron, piece.polyhedron.faces[face_index],
                    piece.geometry.faces[face_index].outward_normal));
                for (const auto vertex : output_piece.faces.back()) {
                    if (vertex >= used.size()) {
                        throw std::runtime_error(
                            "Cut-cell 多面体面引用了越界顶点");
                    }
                    used[vertex] = true;
                }
            }
            for (std::size_t vertex = 0; vertex < used.size(); ++vertex) {
                if (!used[vertex]) continue;
                output_piece.compact_vertex[vertex] =
                    output_piece.used_vertices.size();
                output_piece.used_vertices.push_back(vertex);
            }
            point_count +=
                static_cast<std::uint64_t>(output_piece.used_vertices.size());
            pieces.push_back(std::move(output_piece));
        }
    }
    // meshio 5.x 将 polyhedron 按顶点数分块，并假设块按顶点数递增；稳定排序
    // 同时保留同类多面体原有的确定性背景单元顺序。
    std::stable_sort(
        pieces.begin(), pieces.end(),
        [](const OutputPiece& first, const OutputPiece& second) {
            return first.used_vertices.size() < second.used_vertices.size();
        });
    const std::uint64_t cell_count = static_cast<std::uint64_t>(pieces.size());
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开 Cut-cell polyhedron VTU：" + path.string());
    }
    output << std::setprecision(17)
           << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\">\n"
           << "  <UnstructuredGrid>\n"
           << "    <Piece NumberOfPoints=\"" << point_count
           << "\" NumberOfCells=\"" << cell_count << "\">\n"
           << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    for (const auto& output_piece : pieces) {
        for (const auto vertex_id : output_piece.used_vertices) {
            const auto& vertex =
                output_piece.piece->polyhedron.vertices[vertex_id];
                output << vertex.x << ' ' << vertex.y << ' ' << vertex.z << '\n';
        }
    }
    output << "        </DataArray>\n"
           << "      </Points>\n"
           << "      <Cells>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
              "format=\"ascii\">\n";
    std::uint64_t point_base = 0;
    for (const auto& output_piece : pieces) {
        for (std::size_t vertex = 0;
             vertex < output_piece.used_vertices.size(); ++vertex) {
            output << point_base + static_cast<std::uint64_t>(vertex)
                   << (vertex + 1U == output_piece.used_vertices.size() ? '\n' : ' ');
        }
        point_base += static_cast<std::uint64_t>(output_piece.used_vertices.size());
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" "
              "format=\"ascii\">\n";
    std::uint64_t connectivity_offset = 0;
    for (const auto& output_piece : pieces) {
        connectivity_offset +=
            static_cast<std::uint64_t>(output_piece.used_vertices.size());
        output << connectivity_offset << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt8\" Name=\"types\" "
              "format=\"ascii\">\n";
    for (std::uint64_t cell = 0; cell < cell_count; ++cell) output << "42\n";
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"faces\" "
              "format=\"ascii\">\n";
    point_base = 0;
    for (const auto& output_piece : pieces) {
            output << output_piece.faces.size();
            for (const auto& face : output_piece.faces) {
                output << ' ' << face.size();
                for (const auto vertex : face) {
                    output << ' ' << point_base +
                        static_cast<std::uint64_t>(
                            output_piece.compact_vertex[vertex]);
                }
            }
            output << '\n';
            point_base +=
                static_cast<std::uint64_t>(output_piece.used_vertices.size());
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"faceoffsets\" "
              "format=\"ascii\">\n";
    std::uint64_t face_offset = 0;
    for (const auto& output_piece : pieces) {
            ++face_offset;
            for (const auto& face : output_piece.faces) {
                face_offset += 1U +
                               static_cast<std::uint64_t>(face.size());
            }
            output << face_offset << '\n';
    }
    output << "        </DataArray>\n"
           << "      </Cells>\n"
           << "      <CellData>\n"
           << "        <DataArray type=\"UInt64\" Name=\"background_cell_id\" "
              "format=\"ascii\">\n";
    for (const auto& output_piece : pieces) {
        output << output_piece.cell->background_cell_id << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt64\" Name=\"component_id\" "
              "format=\"ascii\">\n";
    for (const auto& output_piece : pieces) {
        output << output_piece.piece->component_id << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt64\" Name=\"global_region_id\" "
              "format=\"ascii\">\n";
    for (const auto& output_piece : pieces) {
        output << output_piece.piece->global_region_id << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Float64\" Name=\"piece_volume\" "
              "format=\"ascii\">\n";
    for (const auto& output_piece : pieces) {
        output << output_piece.piece->geometry.volume << '\n';
    }
    output << "        </DataArray>\n"
           << "      </CellData>\n"
           << "    </Piece>\n"
           << "  </UnstructuredGrid>\n"
           << "</VTKFile>\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入 Cut-cell polyhedron VTU 失败：" + path.string());
    }
}

void write_fluid_tetrahedra_vtu(const std::filesystem::path& path,
                                const ConvexCutCellMesh& mesh) {
    const auto for_each_tetrahedron = [&](const auto& operation) {
        std::uint64_t piece_id = 0;
        for (const auto& cell : mesh.fluid_cells) {
            for (const auto& piece : cell.fluid_polyhedron_pieces) {
                const Vec3 interior = piece.geometry.centroid;
                for (const auto& face : piece.polyhedron.faces) {
                    if (face.vertex_indices.size() < 3U) continue;
                    const Vec3 first =
                        piece.polyhedron.vertices[face.vertex_indices.front()];
                    for (std::size_t vertex = 1U;
                         vertex + 1U < face.vertex_indices.size(); ++vertex) {
                        const Vec3 second = piece.polyhedron.vertices[
                            face.vertex_indices[vertex]];
                        const Vec3 third = piece.polyhedron.vertices[
                            face.vertex_indices[vertex + 1U]];
                        const double volume = std::abs(dot(
                            first - interior,
                            cross(second - interior, third - interior))) /
                                              6.0;
                        if (!(volume > 0.0)) continue;
                        operation(cell, piece, piece_id, interior, first, second,
                                  third, volume);
                    }
                }
                ++piece_id;
            }
        }
    };
    std::uint64_t tetrahedron_count = 0;
    std::uint64_t point_count = 0;
    for_each_tetrahedron([&](const auto&, const auto&, std::uint64_t,
                             const Vec3&, const Vec3&, const Vec3&, const Vec3&,
                             double) { ++tetrahedron_count; });
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            point_count += 1U +
                           static_cast<std::uint64_t>(
                               piece.polyhedron.vertices.size());
        }
    }
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开 Cut-cell tetra VTU：" + path.string());
    }
    output << std::setprecision(17)
           << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\">\n"
           << "  <UnstructuredGrid>\n"
           << "    <Piece NumberOfPoints=\"" << point_count
           << "\" NumberOfCells=\"" << tetrahedron_count << "\">\n"
           << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            output << piece.geometry.centroid.x << ' '
                   << piece.geometry.centroid.y << ' '
                   << piece.geometry.centroid.z << '\n';
            for (const auto& point : piece.polyhedron.vertices) {
                output << point.x << ' ' << point.y << ' ' << point.z << '\n';
            }
        }
    }
    output << "        </DataArray>\n"
           << "      </Points>\n"
           << "      <Cells>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
              "format=\"ascii\">\n";
    std::uint64_t point_base = 0;
    for (const auto& cell : mesh.fluid_cells) {
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            for (const auto& face : piece.polyhedron.faces) {
                if (face.vertex_indices.size() < 3U) continue;
                for (std::size_t vertex = 1U;
                     vertex + 1U < face.vertex_indices.size(); ++vertex) {
                    const Vec3 first = piece.polyhedron.vertices[
                        face.vertex_indices.front()];
                    const Vec3 second =
                        piece.polyhedron.vertices[face.vertex_indices[vertex]];
                    const Vec3 third = piece.polyhedron.vertices[
                        face.vertex_indices[vertex + 1U]];
                    const double volume = std::abs(dot(
                        first - piece.geometry.centroid,
                        cross(second - piece.geometry.centroid,
                              third - piece.geometry.centroid))) /
                                          6.0;
                    if (!(volume > 0.0)) continue;
                    output << point_base << ' '
                           << point_base + 1U + face.vertex_indices.front() << ' '
                           << point_base + 1U + face.vertex_indices[vertex] << ' '
                           << point_base + 1U + face.vertex_indices[vertex + 1U]
                           << '\n';
                }
            }
            point_base += 1U +
                          static_cast<std::uint64_t>(
                              piece.polyhedron.vertices.size());
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" "
              "format=\"ascii\">\n";
    for (std::uint64_t cell = 1; cell <= tetrahedron_count; ++cell) {
        output << 4U * cell << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt8\" Name=\"types\" "
              "format=\"ascii\">\n";
    for (std::uint64_t cell = 0; cell < tetrahedron_count; ++cell) {
        output << "10\n";
    }
    output << "        </DataArray>\n"
           << "      </Cells>\n"
           << "      <CellData>\n";
    const auto write_scalar = [&](const char* type, const char* name,
                                  const auto& value) {
        output << "        <DataArray type=\"" << type << "\" Name=\"" << name
               << "\" format=\"ascii\">\n";
        for_each_tetrahedron(
            [&](const auto& cell, const auto& piece, std::uint64_t piece_id,
                const Vec3&, const Vec3&, const Vec3&, const Vec3&, double volume) {
                output << value(cell, piece, piece_id, volume) << '\n';
            });
        output << "        </DataArray>\n";
    };
    write_scalar("UInt64", "background_cell_id",
                 [](const auto& cell, const auto&, std::uint64_t, double) {
                     return cell.background_cell_id;
                 });
    write_scalar("UInt64", "source_piece_id",
                 [](const auto&, const auto&, std::uint64_t piece_id, double) {
                     return piece_id;
                 });
    write_scalar("UInt64", "component_id",
                 [](const auto&, const auto& piece, std::uint64_t, double) {
                     return piece.component_id;
                 });
    write_scalar("UInt64", "global_region_id",
                 [](const auto&, const auto& piece, std::uint64_t, double) {
                     return piece.global_region_id;
                 });
    write_scalar("Float64", "tetra_volume",
                 [](const auto&, const auto&, std::uint64_t, double volume) {
                     return volume;
                 });
    output << "      </CellData>\n"
           << "    </Piece>\n"
           << "  </UnstructuredGrid>\n"
           << "</VTKFile>\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入 Cut-cell tetra VTU 失败：" + path.string());
    }
}

} // 命名空间 cartmesh
