#pragma once

#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"
#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/QualityContract2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/topology/EdgeIncidence2D.hpp"

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
    RemainderCartesian,
    Termination,
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
    std::vector<std::uint64_t> refinementLineageKeys;
    std::optional<std::size_t> layerIndex;
    std::optional<std::size_t> wallSegment;
    std::optional<std::size_t> stripId;
    double localBackgroundH = 0.0;
    std::vector<Segment2D> embeddedBoundary;
    // True only for geometric cells whose exact one-cell identity is part of
    // the committed layer/termination construction. Remainder cells stay
    // repairable even when they carry the Termination reporting kind.
    bool solverImmutable = false;
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
    // Q5-1: radial rows the outermost ring is cut into so no transition row is
    // radially thicker than the remainder background cell it borders. 1 keeps
    // the historical single-row outer fan.
    std::size_t outerRingRadialSubdivision = 1U;
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
    std::size_t requestedBoundaryLayerCellCount = 0;
    std::size_t zeroLayerColumnCount = 0;
    std::size_t terminationCellCount = 0;
    std::size_t terminationEdgeCount = 0;
    std::size_t transitionPolygonCount = 0;
    std::size_t transitionRingCount = 0;
    std::size_t transitionFinalTangentialSubdivision = 1;
    std::size_t remainderSmallCellCount = 0;
    std::size_t remainderAgglomeratedCellCount = 0;
    std::size_t solverCellCount = 0;
    std::size_t solverQualityAgglomerations = 0;
    std::size_t solverQualityRepartitions = 0;
    std::size_t r1ShortFaceCandidates = 0;
    std::size_t r1LocalCandidates = 0;
    std::size_t r1LocalQualityEvaluations = 0;
    std::size_t r1AcceptedTransactions = 0;
    std::size_t r1CandidateGlobalTopologyBuilds = 0;
    std::size_t r1CandidateFullGlobalQualityEvaluations = 0;
    std::size_t r1GlobalOracleBuilds = 0;
    std::size_t r1MaximumWinnerGlobalOracleBuilds = 0;
    std::size_t r1AuthoritativeFullQualityEvaluations = 0;
    std::size_t q3TerminationCandidates = 0;
    std::size_t q3LocalCandidates = 0;
    std::size_t q3LocalQualityEvaluations = 0;
    std::size_t q3AcceptedTransactions = 0;
    std::size_t q3CandidateGlobalTopologyBuilds = 0;
    std::size_t q3CandidateFullGlobalQualityEvaluations = 0;
    std::size_t q3GlobalOracleBuilds = 0;
    std::size_t q3MaximumWinnerGlobalOracleBuilds = 0;
    std::size_t q3AuthoritativeFullQualityEvaluations = 0;
    std::size_t q32RepartitionCandidates = 0;
    std::size_t q32LocalCandidates = 0;
    std::size_t q32LocalQualityEvaluations = 0;
    std::size_t q32AcceptedTransactions = 0;
    std::size_t q32CandidateGlobalTopologyBuilds = 0;
    std::size_t q32CandidateFullGlobalQualityEvaluations = 0;
    std::size_t q32GlobalOracleBuilds = 0;
    std::size_t q32MaximumWinnerGlobalOracleBuilds = 0;
    std::size_t q32AuthoritativeFullQualityEvaluations = 0;
    std::size_t q33GroupedCandidates = 0;
    std::size_t q33LocalCandidates = 0;
    std::size_t q33LocalQualityEvaluations = 0;
    std::size_t q33AcceptedTransactions = 0;
    std::size_t q33CandidateGlobalTopologyBuilds = 0;
    std::size_t q33CandidateFullGlobalQualityEvaluations = 0;
    std::size_t q33GlobalOracleBuilds = 0;
    std::size_t q33MaximumWinnerGlobalOracleBuilds = 0;
    std::size_t q33AuthoritativeFullQualityEvaluations = 0;
    std::size_t q41ConstructionCandidates = 0;
    std::size_t q41LocalCandidates = 0;
    std::size_t q41LocalQualityEvaluations = 0;
    std::size_t q41AcceptedConstructions = 0;
    std::size_t q41CandidateGlobalTopologyBuilds = 0;
    std::size_t q41CandidateFullGlobalQualityEvaluations = 0;
    std::size_t q41GlobalOracleBuilds = 0;
    std::size_t q41MaximumWinnerGlobalOracleBuilds = 0;
    std::size_t q41AuthoritativeFullQualityEvaluations = 0;
    std::size_t q3HardVolumeRatioBefore = 0;
    std::size_t q3HardVolumeRatioAfter = 0;
    std::size_t q3HardFaceWeightBefore = 0;
    std::size_t q3HardFaceWeightAfter = 0;
    std::size_t q3HardShortFaceBefore = 0;
    std::size_t q3HardShortFaceAfter = 0;
    std::size_t q32HardVolumeRatioBefore = 0;
    std::size_t q32HardVolumeRatioAfter = 0;
    std::size_t q32HardFaceWeightBefore = 0;
    std::size_t q32HardFaceWeightAfter = 0;
    std::size_t q32HardShortFaceBefore = 0;
    std::size_t q32HardShortFaceAfter = 0;
    std::size_t q33HardVolumeRatioBefore = 0;
    std::size_t q33HardVolumeRatioAfter = 0;
    std::size_t q33HardFaceWeightBefore = 0;
    std::size_t q33HardFaceWeightAfter = 0;
    std::size_t q33HardShortFaceBefore = 0;
    std::size_t q33HardShortFaceAfter = 0;
    std::size_t q41HardVolumeRatioBefore = 0;
    std::size_t q41HardVolumeRatioAfter = 0;
    std::size_t q41HardFaceWeightBefore = 0;
    std::size_t q41HardFaceWeightAfter = 0;
    std::size_t q41HardShortFaceBefore = 0;
    std::size_t q41HardShortFaceAfter = 0;
    // Whole-build measured totals, so the R1 counts above can be read against
    // the cost of the rest of the pipeline instead of in isolation.
    std::size_t buildGlobalTopologyCalls = 0;
    std::size_t buildGlobalTopologyInputCells = 0;
    double buildGlobalTopologySeconds = 0.0;
    std::size_t fullSolverQualityEvaluations = 0;
    double solverRepairSeconds = 0.0;
    double r1MinimumFaceOverLocalHBefore = 0.0;
    double r1MinimumFaceOverLocalHAfter = 0.0;
    double r1RepairSeconds = 0.0;
    double q3MaximumVolumeRatioSeverityBefore = 0.0;
    double q3MaximumVolumeRatioSeverityAfter = 0.0;
    double q3MaximumFaceWeightSeverityBefore = 0.0;
    double q3MaximumFaceWeightSeverityAfter = 0.0;
    double q3RepairSeconds = 0.0;
    double q32RepairSeconds = 0.0;
    double q33RepairSeconds = 0.0;
    double q41ConstructionSelectionSeconds = 0.0;
    bool r1PatchOutsideStableIdsUnchanged = true;
    bool r1LocalDeltaMatchesGlobalOracle = true;
    bool r1LocalWinnerMatchesGlobalAuthority = true;
    bool q3PatchOutsideStableIdsUnchanged = true;
    bool q3LocalDeltaMatchesGlobalOracle = true;
    bool q3LocalWinnerMatchesGlobalAuthority = true;
    bool q3TransactionBoundReached = false;
    bool q32PatchOutsideStableIdsUnchanged = true;
    bool q32LocalDeltaMatchesGlobalOracle = true;
    bool q32LocalWinnerMatchesGlobalAuthority = true;
    bool q32TransactionBoundReached = false;
    bool q33PatchOutsideStableIdsUnchanged = true;
    bool q33LocalDeltaMatchesGlobalOracle = true;
    bool q33LocalWinnerMatchesGlobalAuthority = true;
    bool q33TransactionBoundReached = false;
    bool q41PatchOutsideStableIdsUnchanged = true;
    bool q41LocalDeltaMatchesGlobalOracle = true;
    bool q41LocalWinnerMatchesGlobalAuthority = true;
    bool q41ConstructionBoundReached = false;
    // True when Q4-1 was requested but every attempt with it enabled failed a
    // final gate, so the committed mesh is the unchanged non-Q4 hybrid. All
    // q41* counters above then describe no committed construction.
    bool q41ConstructionSelectionDeclined = false;
    // Q5-1: radial rows actually committed in the outermost transition ring and
    // whether the requested subdivision had to be declined to keep the mesh.
    std::size_t q51OuterTransitionRadialSubdivision = 1;
    bool q51OuterTransitionRadialDeclined = false;
    std::size_t unifiedVertexCount = 0;
    std::size_t unifiedEdgeCount = 0;
    std::size_t unifiedCellCount = 0;
    double transitionTargetCellSize = 0.0;
    double transitionMaxOuterEdgeLength = 0.0;
    double transitionMaxLastLayerSpacing = 0.0;
    double transitionRingThickness = 0.0;
    double transitionTotalThickness = 0.0;
    double terminationGrowthRatio = 0.0;
    double solidArea = 0.0;
    double outerEnvelopeArea = 0.0;
    double layerArea = 0.0;
    double remainderArea = 0.0;
    double expectedFluidArea = 0.0;
    double actualFluidArea = 0.0;
    double areaError = 0.0;
};

struct SourceLineageAudit2D {
    std::size_t solverCellCount = 0;
    std::size_t lineageCandidateChecks = 0;
    std::size_t oracleCandidateChecks = 0;
    std::size_t mismatchedCells = 0;
    bool oracleVerified = false;

    [[nodiscard]] bool pass() const noexcept { return mismatchedCells == 0; }
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
    bool sharedIntersectionConstruction = true;
    // Debug/acceptance oracle. Production uses propagated lineage directly;
    // the opt-in oracle rescans all source polygons and fails on disagreement.
    bool verifySourceLineageOracle = false;
    // Q3-1 is opt-in while validation is intentionally limited to narrow-gap.
    bool enableTerminationQualityOptimization = false;
    bool enableTerminationQualityRepartition = false;
    bool enableTerminationGroupedRepartition = false;
    // Q4-1 moves one bounded typed convex termination-patch construction
    // choice to the H4 remainder boundary, before unified topology commit.
    bool enableTerminationConstructionQualitySelection = false;
    // Q5-1 cuts the outermost transition row radially so it is not thicker
    // than the remainder background cell it borders. Opt-in while validation
    // is limited to the two no-termination geometries that motivated it.
    bool enableOuterTransitionRadialMatching = false;
    // Radial rows are added only while the outer row stays thicker than this
    // multiple of the remainder background cell size. A value of 1 means the
    // outer transition row may be at most one background cell thick.
    double outerTransitionRadialTargetCells = 1.0;
    std::size_t maximumOuterTransitionRadialSubdivision = 4U;
    TolerancePolicy tolerance{};
    double areaToleranceMultiplier = 256.0;
    double interfaceToleranceMultiplier = 128.0;

    // Resolved transition controls used by the six-argument implementation.
    // Product callers should normally use the five-argument overload below,
    // which derives these values from interface length scale and remainder h.
    double transitionCellWidthMultiplier = 1.2;
    std::size_t transitionRingCount = 3U;
    // Resolved Q5-1 value. 1 reproduces the historical single-row outer fan.
    std::size_t transitionOuterRingRadialSubdivision = 1U;
    double terminationGrowthRatio = 1.45;
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
    EdgeIncidenceStore2D constructionIncidence;
    EdgeIncidenceStore2D solverIncidence;
    std::vector<HybridCellRecord2D> cellRecords;
    HybridInterfaceAudit2D interfaceAudit;
    HybridInterfaceAudit2D solverInterfaceAudit;
    HybridMeshMetrics2D metrics;
    MeshQualityReport2D meshQuality;
    SolverQualityReport2D solverQuality;
    QualityContractReport2D qualityContract;
    SmallCellReport2D remainderSmallCells;
    AgglomerationResult2D remainderStabilization;
    SolverTopologyResult2D solverTopologyReport;
    SourceLineageAudit2D sourceLineageAudit;
    std::vector<CanonicalizedIntersection2D> canonicalizedIntersections;
    std::vector<ConstructionRecoveryRequest2D> constructionRecoveryRequests;
    std::vector<QuadtreeLocalRefinementReport2D> constructionRecoveryRefinements;
    QuadtreeBalanceReport2D balance;
    HybridMeshFailure2D failure;

    [[nodiscard]] bool success() const noexcept {
        return status == HybridMeshStatus2D::Success;
    }
};

enum class H4MeshMode2D { Failed, Hybrid, PureCutCellFallback };
enum class H4FallbackStage2D { None, RequestedLayers, LocalReduction, HybridCandidate };

struct PureCutCellFallback2D {
    std::vector<CutCell2D> sourceCells;
    TopologyMesh2D topology;
    TopologyMesh2D solverTopology;
    SmallCellReport2D smallCells;
    AgglomerationResult2D stabilization;
    MeshQualityReport2D meshQuality;
    SolverTopologyResult2D solverTopologyReport;
    SolverQualityReport2D solverQuality;
    QuadtreeBalanceReport2D balance;
    std::size_t quadtreeLeafCount = 0;
    double expectedFluidArea = 0.0;
    double actualFluidArea = 0.0;
    double areaError = 0.0;
    std::string failureMessage;

    [[nodiscard]] bool valid() const noexcept {
        return failureMessage.empty() && topology.valid() && solverTopology.valid() &&
               meshQuality.valid() && solverQuality.valid();
    }
};

struct RobustH4BuildResult2D {
    H4MeshMode2D mode = H4MeshMode2D::Failed;
    H4FallbackStage2D fallbackStage = H4FallbackStage2D::None;
    BoundaryLayerBuildResult2D requestedLayerCandidate;
    BoundaryLayerBuildResult2D localLayerCandidate;
    HybridMeshBuildResult2D hybridCandidate;
    PureCutCellFallback2D fallback;

    [[nodiscard]] bool success() const noexcept {
        return mode==H4MeshMode2D::Hybrid ||
               (mode==H4MeshMode2D::PureCutCellFallback && fallback.valid());
    }
};

// Transactional H4-3 product path. Pure Cut-cell is attempted only after the
// requested strip, local reduction/termination and conformal hybrid candidates
// have all failed their unchanged topology and solver-quality gates.
[[nodiscard]] RobustH4BuildResult2D buildRobustH4Mesh2D(
    const std::vector<WallChain2D>& wallChains,
    const LayerParameters2D& layerParameters,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t maxLevel,
    const QuadtreeRefinementPolicy2D& refinement,
    const BoundaryLayerPolicy2D& layerPolicy = {},
    const HybridMeshPolicy2D& hybridPolicy = {});

[[nodiscard]] const char* h4MeshModeName(H4MeshMode2D mode) noexcept;
[[nodiscard]] const char* h4FallbackStageName(H4FallbackStage2D stage) noexcept;

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
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    double radialTargetCells = 1.0,
    std::size_t maximumRadialSubdivision = 4U) noexcept {
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
    // Q5-1 radial matching. The fan is graded tangentially but every row keeps
    // the full ring thickness, so the outermost row can be several background
    // cells deep while the Cut cells it borders are a fraction of one. Report
    // the row count that brings the outer row to at most one background cell;
    // whether it is applied is a policy decision in the builder.
    plan.outerRingRadialSubdivision = 1U;
    if (radialTargetCells > 0.0 && maximumRadialSubdivision > 1U) {
        while (plan.outerRingRadialSubdivision < maximumRadialSubdivision &&
               ringThickness / static_cast<double>(plan.outerRingRadialSubdivision) >
                   radialTargetCells * targetCellSize) {
            ++plan.outerRingRadialSubdivision;
        }
    }
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
[[nodiscard]] inline HybridMeshBuildResult2D buildAutomaticHybridWithConstruction2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    const HybridMeshPolicy2D& basePolicy) {
    const auto plan = resolveAutomaticHybridTransitionPlan2D(
        boundaryLayers, domain, remainderRefinement,
        basePolicy.outerTransitionRadialTargetCells,
        basePolicy.maximumOuterTransitionRadialSubdivision);
    if (!plan) {
        HybridMeshBuildResult2D result;
        result.failure.reason = HybridMeshFailureReason2D::InvalidInput;
        result.failure.message =
            "automatic H4 transition plan could not be resolved safely";
        return result;
    }

    HybridMeshPolicy2D resolvedPolicy=basePolicy;
    resolvedPolicy.transitionRingCount = plan->ringCount;
    resolvedPolicy.transitionCellWidthMultiplier =
        plan->ringCount==0U?1.0:
        plan->ringThickness / plan->targetCellSize;
    resolvedPolicy.transitionOuterRingRadialSubdivision =
        basePolicy.enableOuterTransitionRadialMatching
            ?plan->outerRingRadialSubdivision:1U;
    HybridMeshBuildResult2D result;
    // Optional construction-quality selections are improvements, never
    // preconditions for producing a hybrid mesh. Relax them one at a time --
    // newest and most invasive first -- before the caller is allowed to
    // consider pure Cut-cell fallback: losing the whole boundary layer is a far
    // larger regression than declining a bounded construction choice.
    std::vector<HybridMeshPolicy2D> ladder{resolvedPolicy};
    if (resolvedPolicy.enableTerminationConstructionQualitySelection) {
        auto step=ladder.back();
        step.enableTerminationConstructionQualitySelection=false;
        ladder.push_back(step);
    }
    if (ladder.back().transitionOuterRingRadialSubdivision>1U) {
        auto step=ladder.back();
        step.transitionOuterRingRadialSubdivision=1U;
        ladder.push_back(step);
    }
    for (std::size_t pass=0;pass<ladder.size();++pass) {
        resolvedPolicy=ladder[pass];
        if (boundaryLayers.localReductionApplied) {
            // Different valid quadtree phases can place a graded termination front
            // arbitrarily close to a nested Cartesian line. Try a short, fixed and
            // deterministic family of geometric growth ratios; commit only a fully
            // solver-valid candidate. Every attempt is transactional.
            constexpr double candidates[]{1.45,1.55,1.50};
            for (const double growth:candidates) {
                resolvedPolicy.terminationGrowthRatio=growth;
                auto attempt=buildConformalHybridMesh2D(
                    boundaryLayers,domain,originalWalls,remainderMaxLevel,
                    remainderRefinement,resolvedPolicy);
                result=std::move(attempt);
                if (result.success()) break;
            }
        } else {
            result = buildConformalHybridMesh2D(
                boundaryLayers, domain, originalWalls, remainderMaxLevel,
                remainderRefinement, resolvedPolicy);
        }
        if (result.success()) break;
    }
    result.metrics.q41ConstructionSelectionDeclined=
        ladder.front().enableTerminationConstructionQualitySelection &&
        result.success() &&
        !resolvedPolicy.enableTerminationConstructionQualitySelection;
    result.metrics.q51OuterTransitionRadialSubdivision=
        result.success()?resolvedPolicy.transitionOuterRingRadialSubdivision:1U;
    result.metrics.q51OuterTransitionRadialDeclined=
        ladder.front().transitionOuterRingRadialSubdivision>1U &&
        result.success() &&
        resolvedPolicy.transitionOuterRingRadialSubdivision<=1U;
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

[[nodiscard]] inline HybridMeshBuildResult2D buildAutomaticHybridWithConstruction2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    bool sharedIntersectionConstruction) {
    HybridMeshPolicy2D policy;
    policy.sharedIntersectionConstruction=sharedIntersectionConstruction;
    return buildAutomaticHybridWithConstruction2D(
        boundaryLayers,domain,originalWalls,remainderMaxLevel,
        remainderRefinement,policy);
}

// Keep the original five-argument API and the six-argument policy overload
// unambiguous, including callers that explicitly pass {} as their policy.
[[nodiscard]] inline HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement) {
    return buildAutomaticHybridWithConstruction2D(boundaryLayers,domain,originalWalls,
        remainderMaxLevel,remainderRefinement,true);
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

// Wall-time measurements are written separately from the report because the
// report is compared byte-for-byte by the determinism regressions and timings
// are not reproducible.
[[nodiscard]] bool writeHybridProfileJson2D(
    const HybridMeshBuildResult2D& result,
    const std::filesystem::path& path,
    std::string* error = nullptr);

} // namespace cartmesh2d
