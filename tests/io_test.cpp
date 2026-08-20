#include "cartmesh2d/io/MeshIO2D.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {
int failures = 0;
void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

CutCell2D fullCell(std::size_t id, std::uint64_t key,
                   const AABB2D& bounds) {
    CutCell2D cell;
    cell.sourceId = id;
    cell.sourceKey = key;
    cell.backgroundBounds = bounds;
    cell.kind = CutCellKind::Full;
    cell.fluidPolygon = {{
        bounds.min,
        {bounds.max.x, bounds.min.y},
        bounds.max,
        {bounds.min.x, bounds.max.y}
    }};
    cell.area = cell.fluidPolygon.area();
    cell.areaFraction = 1.0;
    cell.centroid = cell.fluidPolygon.centroid();
    return cell;
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream in(path);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}
} // namespace

int main() {
    const Domain2D domain{{{0.0,0.0},{2.0,1.0}}};
    BoundaryLoop boundary({{0.0,0.0},{2.0,0.0},{2.0,1.0},{0.0,1.0}});
    std::vector<CutCell2D> cells;
    cells.push_back(fullCell(10, 3, {{0.0,0.0},{1.0,1.0}}));
    cells.push_back(fullCell(20, 4, {{1.0,0.0},{2.0,1.0}}));
    const auto topology = buildGlobalTopology(cells, domain, boundary);
    check(topology.valid(), "IO fixture topology is valid");

    const auto base = std::filesystem::temp_directory_path() / "cartmesh2d_stage6_io";
    std::filesystem::create_directories(base);
    const auto vtkPath = base / "mesh.vtk";
    const auto cm2dPath = base / "mesh.cm2d";
    const auto cm2dPath2 = base / "mesh_repeat.cm2d";

    std::string error;
    check(writeLegacyVtk2D(topology, vtkPath, &error),
          "legacy VTK export succeeds: " + error);
    error.clear();
    check(writeCm2dTopology(topology, cm2dPath, &error),
          "CM2D export succeeds: " + error);
    error.clear();
    check(writeCm2dTopology(topology, cm2dPath2, &error),
          "repeated CM2D export succeeds: " + error);

    const std::string vtk = readText(vtkPath);
    check(vtk.find("DATASET UNSTRUCTURED_GRID") != std::string::npos,
          "VTK export declares unstructured grid");
    check(vtk.find("POINTS 6 double") != std::string::npos,
          "VTK export preserves point count");
    check(vtk.find("CELLS 2 10") != std::string::npos,
          "VTK export preserves cell count/connectivity size");

    const auto readback = readCm2dTopology(cm2dPath);
    check(readback.valid(), "CM2D independent read-back succeeds: " + readback.error);
    if (readback.valid()) {
        check(readback.topology.vertices.size() == topology.vertices.size(),
              "read-back vertex count matches memory topology");
        check(readback.topology.edges.size() == topology.edges.size(),
              "read-back edge count matches memory topology");
        check(readback.topology.cells.size() == topology.cells.size(),
              "read-back cell count matches memory topology");
        check(readback.topology.edges[1].owner == topology.edges[1].owner &&
              readback.topology.edges[1].neighbour == topology.edges[1].neighbour,
              "read-back owner/neighbour connectivity matches");
        check(readback.topology.cells[0].vertices == topology.cells[0].vertices &&
              readback.topology.cells[0].edges == topology.cells[0].edges,
              "read-back cell loops match exactly");
    }

    check(readText(cm2dPath) == readText(cm2dPath2),
          "same topology exports byte-identical CM2D text");

    auto invalid = topology;
    invalid.audit.nonManifoldEdges = 1;
    error.clear();
    check(!writeCm2dTopology(invalid, base / "invalid.cm2d", &error),
          "invalid topology is rejected before export");

    const auto corruptPath = base / "corrupt.cm2d";
    {
        std::ofstream out(corruptPath);
        out << "CM2D 99\n";
    }
    const auto corrupt = readCm2dTopology(corruptPath);
    check(!corrupt.valid() && !corrupt.error.empty(),
          "unsupported/corrupt CM2D input fails explicitly");

    std::filesystem::remove_all(base);

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-6B IO tests passed\n";
    return EXIT_SUCCESS;
}
