#include "cartmesh/io/DiagnosticVtkWriter.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace cartmesh {
namespace {

struct Marker {
    Vec3 position{};
    std::uint8_t issue_code{};
};

[[nodiscard]] Vec3 edge_midpoint(const SurfaceDiagnosticLocation& location) noexcept {
    return (location.edge_start + location.edge_end) * 0.5;
}

} // 匿名命名空间

void write_surface_diagnostic_vtp(const std::filesystem::path& path, const SurfaceMesh& surface,
                                  const SurfaceDiagnostics& diagnostics) {
    std::vector<Marker> markers;
    const auto add_triangle_markers = [&](const std::vector<std::uint64_t>& triangle_ids,
                                          std::uint8_t issue_code) {
        for (const auto triangle_id : triangle_ids) {
            if (triangle_id < surface.triangles().size()) {
                markers.push_back(
                    {surface.triangles()[static_cast<std::size_t>(triangle_id)].centroid(),
                     issue_code});
            }
        }
    };
    add_triangle_markers(diagnostics.degenerate_triangle_examples, 1);
    add_triangle_markers(diagnostics.duplicate_triangle_examples, 2);
    for (const auto& location : diagnostics.boundary_edge_examples) {
        markers.push_back({edge_midpoint(location), 3});
    }
    for (const auto& location : diagnostics.non_manifold_edge_examples) {
        markers.push_back({edge_midpoint(location), 4});
    }
    for (const auto& location : diagnostics.orientation_conflict_examples) {
        markers.push_back({edge_midpoint(location), 5});
    }
    for (const auto& vertex : diagnostics.non_manifold_vertex_examples) {
        markers.push_back({vertex.position, 6});
    }
    for (const auto& component : diagnostics.components) {
        if (diagnostics.component_orientation_mismatch_count > 0 &&
            !component.orientation_matches_nesting) {
            markers.push_back({component.sample_position, 7});
        }
        if (component.bounding_diagonal <=
            diagnostics.small_component_diagonal_threshold) {
            markers.push_back({component.bounds.center(), 11});
        }
    }
    const auto add_pair_markers = [&](const auto& pairs, std::uint8_t issue_code) {
        for (const auto& pair : pairs) {
            markers.push_back({pair.first_position, issue_code});
            markers.push_back({pair.second_position, issue_code});
        }
    };
    add_pair_markers(diagnostics.overlapping_triangle_examples, 8);
    add_pair_markers(diagnostics.self_intersection_examples, 9);
    add_pair_markers(diagnostics.non_adjacent_contact_examples, 10);

    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开表面诊断 VTP：" + path.string());
    }
    output << std::setprecision(17);
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
           << "  <PolyData>\n"
           << "    <Piece NumberOfPoints=\"" << markers.size() << "\" NumberOfVerts=\""
           << markers.size()
           << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
           << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& marker : markers) {
        output << marker.position.x << ' ' << marker.position.y << ' ' << marker.position.z
               << '\n';
    }
    output << "        </DataArray>\n"
           << "      </Points>\n"
           << "      <Verts>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t index = 0; index < markers.size(); ++index) {
        output << index << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t index = 0; index < markers.size(); ++index) {
        output << index + 1 << '\n';
    }
    output << "        </DataArray>\n"
           << "      </Verts>\n"
           << "      <CellData>\n"
           << "        <DataArray type=\"UInt8\" Name=\"issue_code\" format=\"ascii\">\n";
    for (const auto& marker : markers) {
        output << static_cast<unsigned>(marker.issue_code) << '\n';
    }
    output << "        </DataArray>\n"
           << "      </CellData>\n"
           << "    </Piece>\n"
           << "  </PolyData>\n"
           << "</VTKFile>\n";
    if (!output) {
        throw std::runtime_error("写入表面诊断 VTP 失败：" + path.string());
    }
}

} // 命名空间 cartmesh
