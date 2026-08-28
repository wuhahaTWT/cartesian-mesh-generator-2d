#pragma once

#include "cartmesh2d/stabilization/SmallCell2D.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class QualityIssueCode2D {
    InvalidTopology,
    InvalidCellGeometry,
    InvalidEdgeGeometry,
    InvalidSourceCutCell
};

struct QualityIssue2D {
    QualityIssueCode2D code;
    std::size_t objectId = 0;
    std::string message;
};

struct MeshQualityReport2D {
    std::size_t vertexCount = 0;
    std::size_t edgeCount = 0;
    std::size_t cellCount = 0;
    std::size_t internalEdgeCount = 0;
    std::size_t boundaryEdgeCount = 0;

    std::size_t sourceCutCellCount = 0;
    std::size_t sourceFullCellCount = 0;
    std::size_t sourceSmallCellCount = 0;

    double minCellArea = 0.0;
    double minEdgeLength = 0.0;
    // Construction-only polygon diagnostics. These are deliberately not
    // solver or OpenFOAM aspect/skewness metrics.
    double maxCellEdgeLengthRatio = 0.0;
    double maxCentroidVertexMeanOffsetNormalized = 0.0;
    double minCutCellAreaFraction = 1.0;

    std::map<unsigned, std::size_t> levelDistribution;
    TopologyAudit2D topologyAudit;
    std::vector<QualityIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topologyAudit.pass();
    }
};

[[nodiscard]] MeshQualityReport2D evaluateMeshQuality(
    const TopologyMesh2D& topology,
    const std::vector<CutCell2D>& sourceCutCells = {},
    const SmallCellReport2D* smallCells = nullptr,
    const TolerancePolicy& tol = {});

[[nodiscard]] std::string qualityReportToJson(
    const MeshQualityReport2D& report,
    int indentSpaces = 2);

} // namespace cartmesh2d
