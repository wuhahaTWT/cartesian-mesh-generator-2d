#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <string>

namespace cartmesh2d {
namespace {

void setError(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

void writeJsonString(std::ostream& out, const std::string& value) {
    out << '"';
    for (const char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << character; break;
        }
    }
    out << '"';
}

} // namespace

const char* boundaryLayerFailureReasonName(
    BoundaryLayerFailureReason2D reason) noexcept {
    switch (reason) {
    case BoundaryLayerFailureReason2D::None: return "none";
    case BoundaryLayerFailureReason2D::InvalidParameters: return "invalid_parameters";
    case BoundaryLayerFailureReason2D::InvalidChain: return "invalid_chain";
    case BoundaryLayerFailureReason2D::DegenerateSegment: return "degenerate_segment";
    case BoundaryLayerFailureReason2D::ConcaveCorner: return "concave_corner";
    case BoundaryLayerFailureReason2D::SharpCorner: return "sharp_corner";
    case BoundaryLayerFailureReason2D::ConflictingHalfPlanes: return "conflicting_half_planes";
    case BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit: return "thickness_exceeds_safe_limit";
    case BoundaryLayerFailureReason2D::EnvelopeSelfIntersection: return "envelope_self_intersection";
    case BoundaryLayerFailureReason2D::EnvelopeWallIntersection: return "envelope_wall_intersection";
    case BoundaryLayerFailureReason2D::ChainCollision: return "chain_collision";
    case BoundaryLayerFailureReason2D::NegativeAreaCell: return "negative_area_cell";
    case BoundaryLayerFailureReason2D::NonConvexCell: return "non_convex_cell";
    case BoundaryLayerFailureReason2D::SelfIntersectingCell: return "self_intersecting_cell";
    case BoundaryLayerFailureReason2D::OverlappingCells: return "overlapping_cells";
    case BoundaryLayerFailureReason2D::DuplicateGeometricVertex: return "duplicate_geometric_vertex";
    case BoundaryLayerFailureReason2D::TopologyConflict: return "topology_conflict";
    case BoundaryLayerFailureReason2D::IoFailure: return "io_failure";
    }
    return "unknown";
}

const char* wallVertexKindName(WallVertexKind2D kind) noexcept {
    switch (kind) {
    case WallVertexKind2D::Smooth: return "smooth";
    case WallVertexKind2D::MildConvex: return "mild_convex";
    case WallVertexKind2D::Concave: return "concave";
    case WallVertexKind2D::Sharp: return "sharp";
    case WallVertexKind2D::Degenerate: return "degenerate";
    }
    return "unknown";
}

bool writeBoundaryLayerLegacyVtk2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error) {
    if (!result.success() || result.strips.empty()) {
        setError(error, "cannot write VTK for a failed or empty boundary-layer candidate");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open boundary-layer VTK output");
        return false;
    }
    std::size_t pointCount = 0U;
    std::size_t cellCount = 0U;
    for (const auto& strip : result.strips) {
        pointCount += strip.vertices.size();
        cellCount += strip.cells.size();
    }
    out << std::setprecision(17);
    out << "# vtk DataFile Version 3.0\n";
    out << "cartmesh2d H4-1 boundary-layer strip\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << pointCount << " double\n";
    for (const auto& strip : result.strips) {
        for (const auto& vertex : strip.vertices) {
            out << vertex.point.x << ' ' << vertex.point.y << " 0\n";
        }
    }
    out << "CELLS " << cellCount << ' ' << cellCount * 5U << "\n";
    std::size_t pointOffset = 0U;
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) {
            out << "4";
            for (const auto id : cell.vertices) out << ' ' << pointOffset + id;
            out << '\n';
        }
        pointOffset += strip.vertices.size();
    }
    out << "CELL_TYPES " << cellCount << "\n";
    for (std::size_t i = 0; i < cellCount; ++i) out << "9\n";
    out << "CELL_DATA " << cellCount << "\n";
    out << "SCALARS layer_index int 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.layer << '\n';
    }
    out << "SCALARS wall_segment int 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.wallSegment << '\n';
    }
    out << "SCALARS cell_area double 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.area << '\n';
    }
    if (!out.good()) {
        setError(error, "failed while writing boundary-layer VTK output");
        return false;
    }
    return true;
}

bool writeBoundaryLayerReportJson2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error) {
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open boundary-layer JSON report");
        return false;
    }
    out << std::setprecision(17);
    out << "{\n  \"layer_status\": \""
        << (result.success() ? "success" : "failed") << "\",\n";
    if (!result.success()) {
        out << "  \"failure_reason\": \""
            << boundaryLayerFailureReasonName(result.failure.reason) << "\",\n";
        out << "  \"message\": ";
        writeJsonString(out, result.failure.message);
        out << ",\n  \"chain_id\": " << result.failure.chainId << ",\n";
        out << "  \"vertex_id\": ";
        if (result.failure.vertexId) out << *result.failure.vertexId; else out << "null";
        out << ",\n  \"edge_id\": ";
        if (result.failure.edgeId) out << *result.failure.edgeId; else out << "null";
        out << ",\n  \"cell_id\": ";
        if (result.failure.cellId) out << *result.failure.cellId; else out << "null";
        out << ",\n  \"requested_thickness\": " << result.failure.requestedThickness;
        out << ",\n  \"safe_thickness\": ";
        if (result.failure.safeThickness && std::isfinite(*result.failure.safeThickness)) {
            out << *result.failure.safeThickness;
        } else out << "null";
        out << ",\n  \"n_layers\": " << result.failure.nLayers;
        out << ",\n  \"growth_ratio\": " << result.failure.growthRatio << "\n}\n";
    } else {
        std::size_t totalVertices = 0U;
        std::size_t totalCells = 0U;
        double minArea = std::numeric_limits<double>::infinity();
        double maxArea = 0.0;
        double minThickness = std::numeric_limits<double>::infinity();
        double maxThickness = 0.0;
        for (const auto& strip : result.strips) {
            totalVertices += strip.metrics.vertexCount;
            totalCells += strip.metrics.cellCount;
            minArea = std::min(minArea, strip.metrics.minCellArea);
            maxArea = std::max(maxArea, strip.metrics.maxCellArea);
            minThickness = std::min(minThickness, strip.metrics.minLayerThickness);
            maxThickness = std::max(maxThickness, strip.metrics.maxLayerThickness);
        }
        out << "  \"failure_reason\": \"none\",\n";
        out << "  \"strip_count\": " << result.strips.size() << ",\n";
        out << "  \"layer_vertex_count\": " << totalVertices << ",\n";
        out << "  \"layer_cell_count\": " << totalCells << ",\n";
        out << "  \"min_cell_area\": " << minArea << ",\n";
        out << "  \"max_cell_area\": " << maxArea << ",\n";
        out << "  \"min_layer_thickness\": " << minThickness << ",\n";
        out << "  \"max_layer_thickness\": " << maxThickness << ",\n";
        out << "  \"strips\": [\n";
        for (std::size_t i = 0; i < result.strips.size(); ++i) {
            const auto& strip = result.strips[i];
            out << "    {\"chain_id\": " << strip.wallChain.id
                << ", \"patch\": ";
            writeJsonString(out, strip.wallChain.patchIdentity);
            out << ", \"closed\": " << (strip.wallChain.closed ? "true" : "false")
                << ", \"fluid_side\": \""
                << (strip.wallChain.fluidSide == FluidSide2D::Right ? "right" : "left")
                << "\", \"wall_vertex_count\": " << strip.wallChain.vertices.size()
                << ", \"layer_vertex_count\": " << strip.metrics.vertexCount
                << ", \"layer_cell_count\": " << strip.metrics.cellCount
                << ", \"n_layers\": " << strip.parameters.nLayers
                << ", \"first_layer_thickness\": "
                << strip.parameters.firstLayerThickness
                << ", \"total_thickness\": " << strip.parameters.totalThickness
                << ", \"growth_ratio\": " << strip.parameters.growthRatio
                << ", \"safe_thickness_limit\": ";
            if (std::isfinite(strip.metrics.safeThicknessLimit)) {
                out << strip.metrics.safeThicknessLimit;
            } else out << "null";
            out << '}';
            if (i + 1U != result.strips.size()) out << ',';
            out << '\n';
        }
        out << "  ]\n}\n";
    }
    if (!out.good()) {
        setError(error, "failed while writing boundary-layer JSON report");
        return false;
    }
    return true;
}

} // namespace cartmesh2d
