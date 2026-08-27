#pragma once

#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class HybridCellKind2D {
    BoundaryLayer,
    RemainderCut,
    RemainderCartesian
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

struct HybridTransitionPlan2D {
    std::size_t ringCount = 0U;
    std::size_t finalTangentialSubdivision = 1U;
    double targetCellSize = 0.0;
    double maxOuterEdgeLength = 0.0;
    double maxLastLayerSpacing = 0.0;
    double ringThickness = 0.0;
    double totalThickness = 0.0;
};

struct HybridMeshMetrics2D {
    std::size_t quadtreeLeafCount = 0;
    std::size_t remainderCartesianCellCount = 0;
    std::size_t remainderCutCellCount = 0;
    std::size_t boundaryLayerCellCount = 0;
    std::size_t transitionPolygonCount = 0;
    std::size_t transitionRingCount = 0;
    std::size_t transitionFinalTangentialSubdivision = 1;
    std::size_t remainderSmallCellCount = 0;
    std::size_t remainderAgglomeratedCellCount = 0;
    std::size_t solverCellCount = 0;
    std::size_t solverQualityAgglomerations = 0;
    std::size_t solverQualityRepartitions = 0;
    std::size_t unifiedVertexCount = 0;
    std::size_t unifiedEdgeCount = 0;
    std::size_t unifiedCellCount = 0;
    double transitionTargetCellSize = 0.0;
    double transitionMaxOuterEdgeLength = 0.0;
    double transitionMaxLastLayerSpacing = 0.0;
    double transitionRingThickness = 0.0;
    double transitionTotalThickness = 0.0;
    double solidArea = 0.0;
    double outerEnvelopeArea = 0.0;
    double layerArea = 0.0;
    double remainderArea = 0.0;
    double expectedFluidArea = 0.0;
    double actualFluidArea = 0.0;
    double areaError = 0.0;
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

    // Resolved transition controls used by the six-argument implementation.
    // Product callers should normally use the five-argument overload below,
    // which derives these values from interface length scale and remainder h.
    double transitionCellWidthMultiplier = 1.2;
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

// Expert/internal entry point with an already resolved transition policy.
[[nodiscard]] HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    const HybridMeshPolicy2D& policy);

// Resolve a conservative progressive-fan plan from geometry instead of
// case-specific tuning. The longest fixed outer edge controls tangential
// subdivision toward the remainder target h. At least three rows are retained
// so the interface never jumps directly from a fixed H4-1 edge to the final
// split density. Radial extent is driven by h and the last H4-1 normal spacing,
// avoiding geometry-specific width fitting.
[[nodiscard]] inline std::optional<HybridTransitionPlan2D>
resolveAutomaticHybridTransitionPlan2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const QuadtreeRefinementPolicy2D& remainderRefinement) noexcept {
    if (!boundaryLayers.success() || boundaryLayers.strips.empty() ||
        remainderRefinement.boundaryLevel >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return std::nullopt;
    }

    const double domainScale = std::max(domain.width(), domain.height());
    const double targetCellSize = std::ldexp(
        domainScale, -static_cast<int>(remainderRefinement.boundaryLevel));
    if (!std::isfinite(targetCellSize) || !(targetCellSize > 0.0)) {
        return std::nullopt;
    }

    double maxOuterEdgeLength = 0.0;
    double maxLastLayerSpacing = 0.0;
    for (const auto& strip : boundaryLayers.strips) {
        const auto outer = strip.outerEnvelope();
        if (outer.size() < 3U ||
            strip.parameters.cumulativeNormalDistances.empty()) {
            return std::nullopt;
        }
        for (std::size_t i = 0; i < outer.size(); ++i) {
            const auto& a = outer[i];
            const auto& b = outer[(i + 1U) % outer.size()];
            const double length = std::hypot(b.x - a.x, b.y - a.y);
            if (!std::isfinite(length) || !(length > 0.0)) return std::nullopt;
            maxOuterEdgeLength = std::max(maxOuterEdgeLength, length);
        }
        const auto& distances = strip.parameters.cumulativeNormalDistances;
        const double lastSpacing = distances.size() == 1U
            ? distances.front()
            : distances.back() - distances[distances.size() - 2U];
        if (!std::isfinite(lastSpacing) || !(lastSpacing > 0.0)) {
            return std::nullopt;
        }
        maxLastLayerSpacing = std::max(maxLastLayerSpacing, lastSpacing);
    }

    constexpr std::size_t minimumRingCount = 3U;
    constexpr std::size_t maximumRingCount = 8U;
    constexpr double targetTangentialMultiplier = 2.0;
    constexpr double minimumTotalWidthCells = 5.4;

    std::size_t ringCount = 1U;
    std::size_t finalSubdivision = 1U;
    const double targetTangentialLength =
        targetTangentialMultiplier * targetCellSize;
    while (maxOuterEdgeLength / static_cast<double>(finalSubdivision) >
           targetTangentialLength) {
        if (ringCount >= maximumRingCount) return std::nullopt;
        finalSubdivision *= 2U;
        ++ringCount;
    }
    while (ringCount < minimumRingCount) {
        finalSubdivision *= 2U;
        ++ringCount;
    }

    const double totalThickness = std::max(
        minimumTotalWidthCells * targetCellSize,
        static_cast<double>(ringCount) * maxLastLayerSpacing);
    const double ringThickness = totalThickness / static_cast<double>(ringCount);
    if (!std::isfinite(totalThickness) || !std::isfinite(ringThickness) ||
        !(ringThickness > 0.0)) {
        return std::nullopt;
    }

    HybridTransitionPlan2D plan;
    plan.ringCount = ringCount;
    plan.finalTangentialSubdivision = finalSubdivision;
    plan.targetCellSize = targetCellSize;
    plan.maxOuterEdgeLength = maxOuterEdgeLength;
    plan.maxLastLayerSpacing = maxLastLayerSpacing;
    plan.ringThickness = ringThickness;
    plan.totalThickness = totalThickness;
    return plan;
}

// Product entry point. H4-1 strips and the original wall are immutable; only
// the transition/remainder is adapted. The transition plan is deterministic
// for identical geometry and refinement input.
[[nodiscard]] inline HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement) {
    const auto plan = resolveAutomaticHybridTransitionPlan2D(
        boundaryLayers, domain, remainderRefinement);
    if (!plan) {
        HybridMeshBuildResult2D result;
        result.failure.reason = HybridMeshFailureReason2D::InvalidInput;
        result.failure.message =
            "automatic H4 transition plan could not be resolved safely";
        return result;
    }

    HybridMeshPolicy2D resolvedPolicy;
    resolvedPolicy.transitionRingCount = plan->ringCount;
    resolvedPolicy.transitionCellWidthMultiplier =
        plan->ringThickness / plan->targetCellSize;
    auto result = buildConformalHybridMesh2D(
        boundaryLayers, domain, originalWalls, remainderMaxLevel,
        remainderRefinement, resolvedPolicy);
    result.metrics.transitionRingCount = plan->ringCount;
    result.metrics.transitionFinalTangentialSubdivision =
        plan->finalTangentialSubdivision;
    result.metrics.transitionTargetCellSize = plan->targetCellSize;
    result.metrics.transitionMaxOuterEdgeLength = plan->maxOuterEdgeLength;
    result.metrics.transitionMaxLastLayerSpacing = plan->maxLastLayerSpacing;
    result.metrics.transitionRingThickness = plan->ringThickness;
    result.metrics.transitionTotalThickness = plan->totalThickness;
    return result;
}

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
