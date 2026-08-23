#include "cartmesh2d/quality/Quality2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace cartmesh2d {
namespace {

[[nodiscard]] double distance(const Point2D& a, const Point2D& b) noexcept {
    return std::hypot(b.x - a.x, b.y - a.y);
}

[[nodiscard]] Polygon2D polygonFromTopologyCell(const TopologyCell2D& cell,
                                                const TopologyMesh2D& topology,
                                                bool& ok) {
    Polygon2D polygon;
    polygon.vertices.reserve(cell.vertices.size());
    for (const std::size_t vertexId : cell.vertices) {
        if (vertexId >= topology.vertices.size()) {
            ok = false;
            return polygon;
        }
        polygon.vertices.push_back(topology.vertices[vertexId].point);
    }
    ok = polygon.vertices.size() >= 3;
    return polygon;
}

[[nodiscard]] unsigned sourceLevel(std::uint64_t key) noexcept {
    return static_cast<unsigned>(key & 0x3fU);
}

[[nodiscard]] std::string indent(int count) {
    return std::string(static_cast<std::size_t>(std::max(0, count)), ' ');
}

} // namespace

MeshQualityReport2D evaluateMeshQuality(const TopologyMesh2D& topology,
                                        const std::vector<CutCell2D>& sourceCutCells,
                                        const SmallCellReport2D* smallCells,
                                        const TolerancePolicy& tol) {
    MeshQualityReport2D report;
    report.vertexCount = topology.vertices.size();
    report.edgeCount = topology.edges.size();
    report.cellCount = topology.cells.size();
    report.topologyAudit = topology.audit;

    if (!topology.valid()) {
        report.issues.push_back({QualityIssueCode2D::InvalidTopology, 0,
                                 "quality evaluation requires a valid Stage 2D-4 topology"});
    }

    double minArea = std::numeric_limits<double>::infinity();
    double minEdge = std::numeric_limits<double>::infinity();
    double maxAspect = 0.0;
    double maxSkewness = 0.0;

    for (const auto& edge : topology.edges) {
        if (edge.v0 >= topology.vertices.size() || edge.v1 >= topology.vertices.size() ||
            edge.v0 == edge.v1) {
            report.issues.push_back({QualityIssueCode2D::InvalidEdgeGeometry, edge.id,
                                     "edge references invalid or identical vertices"});
            continue;
        }
        const double length = distance(topology.vertices[edge.v0].point,
                                       topology.vertices[edge.v1].point);
        const double edgeEps = tol.scale(std::max(1.0, length));
        if (!(length > edgeEps)) {
            std::ostringstream detail;
            detail << std::setprecision(17) << "edge length=" << length
                   << " is at or below quality epsilon=" << edgeEps
                   << " vertices=" << edge.v0 << ',' << edge.v1;
            report.issues.push_back({QualityIssueCode2D::InvalidEdgeGeometry, edge.id,
                                     detail.str()});
            continue;
        }
        minEdge = std::min(minEdge, length);
        if (edge.neighbour) {
            ++report.internalEdgeCount;
        } else {
            ++report.boundaryEdgeCount;
        }
    }

    for (const auto& cell : topology.cells) {
        bool polygonOk = false;
        const Polygon2D polygon = polygonFromTopologyCell(cell, topology, polygonOk);
        if (!polygonOk) {
            report.issues.push_back({QualityIssueCode2D::InvalidCellGeometry, cell.id,
                                     "cell does not contain a valid vertex loop"});
            continue;
        }

        const double area = polygon.area();
        const double areaEps = tol.scale(std::max(1.0, area));
        if (!(area > areaEps)) {
            std::ostringstream detail;
            detail << std::setprecision(17) << "cell area=" << area
                   << " is at or below quality epsilon=" << areaEps;
            report.issues.push_back({QualityIssueCode2D::InvalidCellGeometry, cell.id,
                                     detail.str()});
            continue;
        }
        minArea = std::min(minArea, area);

        double cellMinEdge = std::numeric_limits<double>::infinity();
        double cellMaxEdge = 0.0;
        for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
            const double length = distance(polygon.vertices[i],
                                           polygon.vertices[(i + 1) % polygon.vertices.size()]);
            cellMinEdge = std::min(cellMinEdge, length);
            cellMaxEdge = std::max(cellMaxEdge, length);
        }
        const double cellEdgeEps = tol.scale(std::max(1.0, cellMaxEdge));
        if (!(cellMinEdge > cellEdgeEps)) {
            std::ostringstream detail;
            detail << std::setprecision(17) << "cell minimum edge=" << cellMinEdge
                   << " is at or below quality epsilon=" << cellEdgeEps
                   << " maximum edge=" << cellMaxEdge;
            report.issues.push_back({QualityIssueCode2D::InvalidCellGeometry, cell.id,
                                     detail.str()});
            continue;
        }
        maxAspect = std::max(maxAspect, cellMaxEdge / cellMinEdge);

        const auto centroid = polygon.centroid();
        if (!centroid) {
            report.issues.push_back({QualityIssueCode2D::InvalidCellGeometry, cell.id,
                                     "cell polygon centroid could not be evaluated"});
            continue;
        }
        Point2D vertexMean{};
        for (const auto& p : polygon.vertices) {
            vertexMean.x += p.x;
            vertexMean.y += p.y;
        }
        const double invCount = 1.0 / static_cast<double>(polygon.vertices.size());
        vertexMean.x *= invCount;
        vertexMean.y *= invCount;
        const double skewness = distance(*centroid, vertexMean) / std::sqrt(area);
        maxSkewness = std::max(maxSkewness, skewness);
    }

    report.minCellArea = std::isfinite(minArea) ? minArea : 0.0;
    report.minEdgeLength = std::isfinite(minEdge) ? minEdge : 0.0;
    report.maxEdgeAspectRatio = maxAspect;
    report.maxCentroidSkewness = maxSkewness;

    bool hasCutFraction = false;
    for (const auto& cut : sourceCutCells) {
        if (cut.kind == CutCellKind::Empty || cut.kind == CutCellKind::Unsupported) continue;
        const unsigned level = sourceLevel(cut.sourceKey);
        ++report.levelDistribution[level];
        if (cut.kind == CutCellKind::Cut) {
            ++report.sourceCutCellCount;
            const double eps = tol.scale(std::max(1.0, std::abs(cut.areaFraction)));
            if (!(cut.areaFraction > eps && cut.areaFraction <= 1.0 + eps)) {
                report.issues.push_back({QualityIssueCode2D::InvalidSourceCutCell, cut.sourceId,
                                         "source Cut-cell alpha lies outside (0,1]"});
                continue;
            }
            report.minCutCellAreaFraction =
                std::min(report.minCutCellAreaFraction, cut.areaFraction);
            hasCutFraction = true;
        } else if (cut.kind == CutCellKind::Full) {
            ++report.sourceFullCellCount;
        }
    }
    if (!hasCutFraction) report.minCutCellAreaFraction = 1.0;

    if (smallCells != nullptr) {
        report.sourceSmallCellCount = smallCells->smallCellCount;
    }

    return report;
}

std::string qualityReportToJson(const MeshQualityReport2D& report, int indentSpaces) {
    const int step = std::max(0, indentSpaces);
    const std::string i1 = indent(step);
    const std::string i2 = indent(step * 2);
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n";
    out << i1 << "\"valid\": " << (report.valid() ? "true" : "false") << ",\n";
    out << i1 << "\"counts\": {\n";
    out << i2 << "\"vertices\": " << report.vertexCount << ",\n";
    out << i2 << "\"edges\": " << report.edgeCount << ",\n";
    out << i2 << "\"cells\": " << report.cellCount << ",\n";
    out << i2 << "\"internal_edges\": " << report.internalEdgeCount << ",\n";
    out << i2 << "\"boundary_edges\": " << report.boundaryEdgeCount << ",\n";
    out << i2 << "\"source_cut_cells\": " << report.sourceCutCellCount << ",\n";
    out << i2 << "\"source_full_cells\": " << report.sourceFullCellCount << ",\n";
    out << i2 << "\"source_small_cells\": " << report.sourceSmallCellCount << "\n";
    out << i1 << "},\n";
    out << i1 << "\"quality\": {\n";
    out << i2 << "\"min_cell_area\": " << report.minCellArea << ",\n";
    out << i2 << "\"min_edge_length\": " << report.minEdgeLength << ",\n";
    out << i2 << "\"max_edge_aspect_ratio\": " << report.maxEdgeAspectRatio << ",\n";
    out << i2 << "\"max_centroid_skewness\": " << report.maxCentroidSkewness << ",\n";
    out << i2 << "\"min_cut_cell_area_fraction\": " << report.minCutCellAreaFraction << "\n";
    out << i1 << "},\n";
    out << i1 << "\"level_distribution\": {";
    bool firstLevel = true;
    for (const auto& [level, count] : report.levelDistribution) {
        if (!firstLevel) out << ',';
        out << "\n" << i2 << "\"" << level << "\": " << count;
        firstLevel = false;
    }
    if (!report.levelDistribution.empty()) out << '\n' << i1;
    out << "},\n";
    out << i1 << "\"topology_audit\": {\n";
    out << i2 << "\"duplicate_vertices\": " << report.topologyAudit.duplicateVertices << ",\n";
    out << i2 << "\"duplicate_edges\": " << report.topologyAudit.duplicateEdges << ",\n";
    out << i2 << "\"orphan_internal_edges\": " << report.topologyAudit.orphanInternalEdges << ",\n";
    out << i2 << "\"non_manifold_edges\": " << report.topologyAudit.nonManifoldEdges << ",\n";
    out << i2 << "\"unclassified_boundary_edges\": " << report.topologyAudit.unclassifiedBoundaryEdges << ",\n";
    out << i2 << "\"open_cell_loops\": " << report.topologyAudit.openCellLoops << ",\n";
    out << i2 << "\"area_mismatches\": " << report.topologyAudit.areaMismatches << "\n";
    out << i1 << "},\n";
    out << i1 << "\"issue_count\": " << report.issues.size() << "\n";
    out << "}\n";
    return out.str();
}

} // namespace cartmesh2d
