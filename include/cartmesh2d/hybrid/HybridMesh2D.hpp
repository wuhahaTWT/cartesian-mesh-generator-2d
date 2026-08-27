#pragma once

#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class HybridCellKind2D {
    BoundaryLayer,
    RemainderCut,
    RemainderCartesian,
    Transition
};

struct HybridCellRecord2D {
    std::size_t topologyCellId = 0;
    std::size_t sourceId = 0;
    HybridCellKind2D kind = HybridCellKind2D::RemainderCartesian;
    std::optional<std::size_t> layerIndex;
    std::optional<std::size_t> wallSegment;
};

struct HybridSourceCell2D {
    std::size_t id = 0;
    HybridCellKind2D kind = HybridCellKind2D::RemainderCartesian;
    Polygon2D polygon;
    double area = 0.0;
    std::optional<std::uint64_t> quadtreeSourceKey;
    std::optional<std::size_t> layerIndex;
    std::optional<std::size_t> wallSegment;
    std::vector<Segment2D> embeddedBoundary;
};

struct HybridInterfaceAudit2D {
    std::size_t interfaceEdgeCount = 0;
    std::size_t interfaceVertexCount = 0;
    std::size_t singleOwnerInterfaceEdges = 0;
    std::size_t wrongCellPairInterfaceEdges = 0;
    std::size_t nonTwoValentInterfaceVertices = 0;
    double expectedLength = 0.0;
    double actualLength = 0.0;
    double lengthError = 0.0;

    [[nodiscard]] bool pass(double tolerance) const noexcept {
        return interfaceEdgeCount > 0U && interfaceVertexCount > 0U &&
               singleOwnerInterfaceEdges == 0U &&
               wrongCellPairInterfaceEdges == 0U &&
               nonTwoValentInterfaceVertices == 0U &&
               std::abs(lengthError) <= tolerance;
    }
};

struct HybridMeshMetrics2D {
    std::size_t quadtreeLeafCount = 0;
    std::size_t remainderCartesianCellCount = 0;
    std::size_t remainderCutCellCount = 0;
    std::size_t boundaryLayerCellCount = 0;
    std::size_t transitionPolygonCount = 0;
    std::size_t remainderSmallCellCount = 0;
    std::size_t remainderAgglomeratedCellCount = 0;
    std::size_t solverCellCount = 0;
    std::size_t solverQualityAgglomerations = 0;
    std::size_t solverQualityRepartitions = 0;
    std::size_t unifiedVertexCount = 0;
    std::size_t unifiedEdgeCount = 0;
    std::size_t unifiedCellCount = 0;
    double solidArea = 0.0;
    double outerEnvelopeArea = 0.0;
    double layerArea = 0.0;
    double remainderArea = 0.0;
    double expectedFluidArea = 0.0;
    double actualFluidArea = 0.0;
    double areaError = 0.0;
    double selectedTransitionWidthMultiplier = 0.0;
};

enum class HybridMeshFailureReason2D {
    None,
    InvalidInput,
    UnsupportedWallSemantics,
    OuterEnvelopeOutsideDomain,
    InvalidOuterEnvelope,
    RemainderRefinementFailed,
    RemainderCutCellFailed,
    RemainderStabilizationFailed,
    LayerConversionFailed,
    UnifiedTopologyFailed,
    NonConformalInterface,
    AreaConservationFailed,
    RegionClassificationConflict,
    QualityFailed,
    SolverTopologyFailed,
    SolverQualityFailed,
    IoFailure
};

struct HybridMeshFailure2D {
    HybridMeshFailureReason2D reason = HybridMeshFailureReason2D::None;
    std::string message;
    std::optional<std::size_t> stripId;
    std::optional<std::size_t> leafId;
    std::optional<std::size_t> cellId;
    std::optional<std::size_t> edgeId;
};

struct HybridMeshPolicy2D {
    TolerancePolicy tolerance{};
    double areaToleranceMultiplier = 256.0;
    double interfaceToleranceMultiplier = 128.0;
    double transitionCellWidthMultiplier = 1.8;
    std::size_t transitionRingCount = 3U;
};

enum class HybridMeshStatus2D { Success, Failed };

struct HybridMeshBuildResult2D {
    HybridMeshStatus2D status = HybridMeshStatus2D::Failed;
    std::vector<BoundaryLayerStrip2D> strips;
    BoundaryRegion2D outerEnvelopeRegion{std::vector<BoundaryLoop>{}};
    std::vector<CutCell2D> remainderSourceCells;
    std::vector<HybridSourceCell2D> sourceCells;
    TopologyMesh2D topology;
    TopologyMesh2D solverTopology;
    std::vector<HybridCellRecord2D> cellRecords;
    HybridInterfaceAudit2D interfaceAudit;
    HybridInterfaceAudit2D solverInterfaceAudit;
    HybridMeshMetrics2D metrics;
    MeshQualityReport2D meshQuality;
    SolverQualityReport2D solverQuality;
    SmallCellReport2D remainderSmallCells;
    AgglomerationResult2D remainderStabilization;
    SolverTopologyResult2D solverTopologyReport;
    QuadtreeBalanceReport2D balance;
    HybridMeshFailure2D failure;

    [[nodiscard]] bool success() const noexcept {
        return status == HybridMeshStatus2D::Success;
    }
};

// H4-2 transaction: the H4-1 strips and all inputs are immutable. The outer
// envelopes are used as the solid boundary for a fresh remainder quadtree and
// Cut-cell pass. Layer and remainder polygons then enter one global topology
// build, where interface fragments share canonical vertex and edge IDs.
[[nodiscard]] HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    const HybridMeshPolicy2D& policy = {});

[[nodiscard]] const char* hybridMeshFailureReasonName(
    HybridMeshFailureReason2D reason) noexcept;
[[nodiscard]] const char* hybridCellKindName(HybridCellKind2D kind) noexcept;

[[nodiscard]] bool writeHybridLegacyVtk2D(
    const HybridMeshBuildResult2D& result,
    const std::filesystem::path& path,
    std::string* error = nullptr);

[[nodiscard]] bool writeHybridReportJson2D(
    const HybridMeshBuildResult2D& result,
    const std::filesystem::path& path,
    std::string* error = nullptr);

} // namespace cartmesh2d
