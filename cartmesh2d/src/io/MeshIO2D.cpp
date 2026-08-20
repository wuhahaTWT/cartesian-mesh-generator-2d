#include "cartmesh2d/io/MeshIO2D.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace cartmesh2d {
namespace {

void setError(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

[[nodiscard]] bool patchFromInt(int value, BoundaryPatch2D& patch) noexcept {
    switch (value) {
    case 0: patch = BoundaryPatch2D::None; return true;
    case 1: patch = BoundaryPatch2D::EmbeddedBoundary; return true;
    case 2: patch = BoundaryPatch2D::DomainBoundary; return true;
    case 3: patch = BoundaryPatch2D::Unclassified; return true;
    default: return false;
    }
}

[[nodiscard]] bool expectToken(std::istream& in, const char* expected,
                               std::string& error) {
    std::string token;
    if (!(in >> token) || token != expected) {
        error = std::string("expected token '") + expected + "'";
        return false;
    }
    return true;
}

} // namespace

bool writeLegacyVtk2D(const TopologyMesh2D& topology,
                      const std::filesystem::path& path,
                      std::string* error) {
    if (!topology.valid()) {
        setError(error, "cannot export invalid topology");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open VTK output file");
        return false;
    }
    std::size_t cellListSize = 0;
    for (const auto& cell : topology.cells) {
        if (cell.vertices.size() < 3) {
            setError(error, "VTK export encountered a cell with fewer than three vertices");
            return false;
        }
        cellListSize += cell.vertices.size() + 1;
    }

    out << std::setprecision(17);
    out << "# vtk DataFile Version 3.0\n";
    out << "cartmesh2d deterministic topology export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << topology.vertices.size() << " double\n";
    for (const auto& vertex : topology.vertices) {
        out << vertex.point.x << ' ' << vertex.point.y << " 0\n";
    }
    out << "CELLS " << topology.cells.size() << ' ' << cellListSize << "\n";
    for (const auto& cell : topology.cells) {
        out << cell.vertices.size();
        for (const std::size_t vertexId : cell.vertices) out << ' ' << vertexId;
        out << '\n';
    }
    out << "CELL_TYPES " << topology.cells.size() << "\n";
    for (std::size_t i = 0; i < topology.cells.size(); ++i) out << "7\n";
    out << "CELL_DATA " << topology.cells.size() << "\n";
    out << "SCALARS geometry_area double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& cell : topology.cells) out << cell.geometryArea << '\n';
    out << "SCALARS source_id unsigned_long 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& cell : topology.cells) out << cell.sourceId << '\n';

    if (!out.good()) {
        setError(error, "failed while writing VTK output");
        return false;
    }
    return true;
}

bool writeCm2dTopology(const TopologyMesh2D& topology,
                       const std::filesystem::path& path,
                       std::string* error) {
    if (!topology.valid()) {
        setError(error, "cannot export invalid topology");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open CM2D output file");
        return false;
    }
    out << std::setprecision(17);
    out << "CM2D 1\n";
    out << "VERTICES " << topology.vertices.size() << '\n';
    for (const auto& vertex : topology.vertices) {
        out << vertex.id << ' ' << vertex.point.x << ' ' << vertex.point.y << '\n';
    }
    out << "EDGES " << topology.edges.size() << '\n';
    for (const auto& edge : topology.edges) {
        const long long neighbour = edge.neighbour
            ? static_cast<long long>(*edge.neighbour) : -1LL;
        out << edge.id << ' ' << edge.v0 << ' ' << edge.v1 << ' '
            << edge.owner << ' ' << neighbour << ' '
            << static_cast<int>(edge.patch) << '\n';
    }
    out << "CELLS " << topology.cells.size() << '\n';
    for (const auto& cell : topology.cells) {
        out << cell.id << ' ' << cell.sourceId << ' ' << cell.sourceKey << ' '
            << cell.geometryArea << ' ' << cell.vertices.size();
        for (const std::size_t vertexId : cell.vertices) out << ' ' << vertexId;
        out << ' ' << cell.edges.size();
        for (const std::size_t edgeId : cell.edges) out << ' ' << edgeId;
        out << '\n';
    }
    out << "AUDIT "
        << topology.audit.duplicateVertices << ' '
        << topology.audit.duplicateEdges << ' '
        << topology.audit.orphanInternalEdges << ' '
        << topology.audit.nonManifoldEdges << ' '
        << topology.audit.unclassifiedBoundaryEdges << ' '
        << topology.audit.openCellLoops << ' '
        << topology.audit.areaMismatches << '\n';
    out << "END\n";
    if (!out.good()) {
        setError(error, "failed while writing CM2D output");
        return false;
    }
    return true;
}

MeshReadback2D readCm2dTopology(const std::filesystem::path& path) {
    MeshReadback2D result;
    std::ifstream in(path);
    if (!in) {
        result.error = "failed to open CM2D input file";
        return result;
    }
    std::string magic;
    int version = 0;
    if (!(in >> magic >> version) || magic != "CM2D" || version != 1) {
        result.error = "unsupported CM2D header/version";
        return result;
    }

    std::size_t vertexCount = 0;
    if (!expectToken(in, "VERTICES", result.error) || !(in >> vertexCount)) return result;
    result.topology.vertices.reserve(vertexCount);
    for (std::size_t i = 0; i < vertexCount; ++i) {
        Vertex2D vertex;
        if (!(in >> vertex.id >> vertex.point.x >> vertex.point.y)) {
            result.error = "failed to parse CM2D vertex";
            return result;
        }
        if (vertex.id != i) {
            result.error = "CM2D vertex ids must be contiguous and deterministic";
            return result;
        }
        result.topology.vertices.push_back(vertex);
    }

    std::size_t edgeCount = 0;
    if (!expectToken(in, "EDGES", result.error) || !(in >> edgeCount)) return result;
    result.topology.edges.reserve(edgeCount);
    for (std::size_t i = 0; i < edgeCount; ++i) {
        Edge2D edge;
        long long neighbour = -1;
        int patchValue = 0;
        if (!(in >> edge.id >> edge.v0 >> edge.v1 >> edge.owner >> neighbour >> patchValue)) {
            result.error = "failed to parse CM2D edge";
            return result;
        }
        if (edge.id != i || edge.v0 >= vertexCount || edge.v1 >= vertexCount) {
            result.error = "CM2D edge contains invalid deterministic ids or vertex references";
            return result;
        }
        if (neighbour >= 0) edge.neighbour = static_cast<std::size_t>(neighbour);
        if (!patchFromInt(patchValue, edge.patch)) {
            result.error = "CM2D edge contains invalid boundary patch value";
            return result;
        }
        result.topology.edges.push_back(edge);
    }

    std::size_t cellCount = 0;
    if (!expectToken(in, "CELLS", result.error) || !(in >> cellCount)) return result;
    result.topology.cells.reserve(cellCount);
    for (std::size_t i = 0; i < cellCount; ++i) {
        TopologyCell2D cell;
        std::size_t vertexLoopSize = 0;
        if (!(in >> cell.id >> cell.sourceId >> cell.sourceKey >> cell.geometryArea >> vertexLoopSize)) {
            result.error = "failed to parse CM2D cell header";
            return result;
        }
        if (cell.id != i || vertexLoopSize < 3) {
            result.error = "CM2D cell contains invalid deterministic id or vertex-loop size";
            return result;
        }
        cell.vertices.resize(vertexLoopSize);
        for (auto& vertexId : cell.vertices) {
            if (!(in >> vertexId) || vertexId >= vertexCount) {
                result.error = "CM2D cell contains invalid vertex reference";
                return result;
            }
        }
        std::size_t edgeLoopSize = 0;
        if (!(in >> edgeLoopSize) || edgeLoopSize != vertexLoopSize) {
            result.error = "CM2D cell edge loop must match its vertex-loop length";
            return result;
        }
        cell.edges.resize(edgeLoopSize);
        for (auto& edgeId : cell.edges) {
            if (!(in >> edgeId) || edgeId >= edgeCount) {
                result.error = "CM2D cell contains invalid edge reference";
                return result;
            }
        }
        result.topology.cells.push_back(std::move(cell));
    }

    if (!expectToken(in, "AUDIT", result.error)) return result;
    auto& audit = result.topology.audit;
    if (!(in >> audit.duplicateVertices >> audit.duplicateEdges >>
          audit.orphanInternalEdges >> audit.nonManifoldEdges >>
          audit.unclassifiedBoundaryEdges >> audit.openCellLoops >>
          audit.areaMismatches)) {
        result.error = "failed to parse CM2D topology audit";
        return result;
    }
    if (!expectToken(in, "END", result.error)) return result;

    for (const auto& edge : result.topology.edges) {
        if (edge.owner >= cellCount || (edge.neighbour && *edge.neighbour >= cellCount)) {
            result.error = "CM2D edge contains invalid owner/neighbour reference";
            return result;
        }
    }
    if (!result.topology.audit.pass()) {
        result.error = "CM2D read-back topology audit is not clean";
        return result;
    }
    return result;
}

} // namespace cartmesh2d
