#pragma once

#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/topology/EdgeIncidence2D.hpp"

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace cartmesh2d {

// One cell of a patch-local evaluation scope.  `loop` and `polygon` are
// parallel: entry i of the loop is the stable identity of polygon vertex i.
struct PatchLocalCell2D {
    std::vector<StableVertexId2D> loop;
    Polygon2D polygon;
    // Parallel to `loop`: entry i marks the face from loop[i] to loop[i+1] as a
    // physical mesh boundary that has no neighbour cell anywhere in the mesh.
    // The producer sets it from topology; it is never inferred from geometry.
    std::vector<bool> physicalBoundaryFace;
    double localBackgroundH = 0.0;
    // Whether the dimensionless face/local_h termination contract rates this
    // cell.  Mirrors the rated-cell flag used by the global repair scan.
    bool ratedForFaceLength = false;
    // Whether this candidate replaces the cell's geometry.  Only in-patch
    // cells contribute cell metrics, and only faces incident to an in-patch
    // cell contribute face metrics.  Everything else is provably unchanged by
    // a boundary-locked patch transaction, so it is excluded from both the
    // base and the candidate scope and cancels in the comparison.
    bool inPatch = false;
};

// Authoritative solver metrics aggregated over one patch-local scope.  Every
// value comes from the same kernels as evaluateSolverQuality2D; no metric is
// approximated and no threshold is relaxed.
struct PatchLocalQuality2D {
    SolverQualityPolicy2D policy;
    std::size_t ratedCellCount = 0;
    std::size_t ratedFaceCount = 0;
    std::size_t ratedInternalFaceCount = 0;

    // Dimensionless termination contract (face length over local background h).
    std::size_t hardShortFaceCount = 0;
    double minimumFaceOverLocalH = 0.0;
    double maximumShortFaceSeverity = 0.0;
    double totalShortFaceSeverity = 0.0;

    // Q1 hard-contract counts kept separate from the legacy solver-safety
    // policy.  Q3 ranks these without changing SolverQualityPolicy2D.
    std::size_t hardVolumeRatioCount = 0;
    double maximumVolumeRatioSeverity = 0.0;
    double totalVolumeRatioSeverity = 0.0;
    std::size_t hardFaceWeightCount = 0;
    double maximumFaceWeightSeverity = 0.0;
    double totalFaceWeightSeverity = 0.0;
    std::size_t hardMinimumInteriorAngleCount = 0;
    std::size_t hardNonOrthogonalityCount = 0;
    std::size_t hardSkewnessCount = 0;
    std::size_t hardAspectCount = 0;

    double maxNonOrthogonalityDeg = 0.0;
    double maxInternalSkewness = 0.0;
    double maxBoundarySkewness = 0.0;
    double maxConcavityDeg = 0.0;
    double maxCellAspect = 0.0;
    double minInteriorAngleDeg = 180.0;
    double minFaceLength = 0.0;
    double minFaceWeight = 1.0;
    double minVolumeRatio = 1.0;
    double minCompactness = 1.0;

    // Solver-policy violations restricted to the scope, scored exactly as the
    // global candidate ranking scores them.
    std::size_t issueCount = 0;
    double maximumIssueSeverity = 0.0;
    double totalIssueSeverity = 0.0;

    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

struct PatchLocalHardLimits2D {
    double minimumFaceWeight = 0.10;
    double minimumVolumeRatio = 0.05;
    double minimumInteriorAngleDeg = 10.0;
    double maximumNonOrthogonalityDeg = 65.0;
    double maximumSkewness = 3.0;
    double maximumAspect = 50.0;
};

// Evaluates one scope.  Fails closed when a rated face does not carry its
// complete incidence, because an incomplete halo cannot bound the metrics of
// the faces it touches.
[[nodiscard]] PatchLocalQuality2D evaluatePatchLocalQuality2D(
    const std::vector<PatchLocalCell2D>& scope, double minimumFaceOverLocalH,
    const SolverQualityPolicy2D& policy = {}, const TolerancePolicy& tol = {},
    const PatchLocalHardLimits2D& hardLimits = {});

namespace detail {

// Sufficient condition for the local short-face lexicographic score to remain
// a strict improvement after any unchanged outside-of-scope faces are restored.
// When hard count is unchanged, maximum severity may not increase and the
// additive total severity must strictly decrease. A local maximum-only decrease
// is insufficient because an unchanged worse face outside the patch can hide it.
[[nodiscard]] constexpr bool patchLocalShortFaceImprovementImpliesGlobal(
    std::size_t candidateHard,double candidateMaximum,double candidateTotal,
    std::size_t baseHard,double baseMaximum,double baseTotal) noexcept {
    if (candidateHard<baseHard) return true;
    if (candidateHard!=baseHard) return false;
    return candidateMaximum<=baseMaximum && candidateTotal<baseTotal;
}

static_assert(patchLocalShortFaceImprovementImpliesGlobal(
    0U,100.0,100.0,1U,1.0,1.0));
static_assert(patchLocalShortFaceImprovementImpliesGlobal(
    1U,4.0,4.0,1U,5.0,5.0));
// Counterexample closed by R1 closeout: a lower local maximum accompanied by a
// higher total can become globally worse when an unchanged outside face owns
// the global maximum.
static_assert(!patchLocalShortFaceImprovementImpliesGlobal(
    1U,4.0,6.0,1U,5.0,5.0));
static_assert(!patchLocalShortFaceImprovementImpliesGlobal(
    1U,6.0,4.0,1U,5.0,5.0));

} // namespace detail

// Patch-local commit gate. Besides the authoritative solver no-worse metrics,
// this enforces the sufficient short-face condition above, so a local pass
// implies that restoring the unchanged remainder cannot reverse the global
// short-face lexicographic improvement.
[[nodiscard]] bool patchLocalQualityNoWorse2D(const PatchLocalQuality2D& candidate,
                                              const PatchLocalQuality2D& base) noexcept;

// Q3 gate: targeted Q1 face-weight/volume-ratio hard counts may improve, while
// the dimensionless short-face score and every non-target solver metric remain
// no worse.  Each targeted count is independently non-increasing.
[[nodiscard]] bool patchLocalTerminationQualityNoWorse2D(
    const PatchLocalQuality2D& candidate,
    const PatchLocalQuality2D& base) noexcept;

struct PatchLocalTerminationRank2D {
    std::ptrdiff_t hardViolationDelta = 0;
    std::ptrdiff_t hardVolumeRatioDelta = 0;
    double maximumVolumeRatioSeverityDelta = 0.0;
    double totalVolumeRatioSeverityDelta = 0.0;
    std::ptrdiff_t hardFaceWeightDelta = 0;
    double maximumFaceWeightSeverityDelta = 0.0;
    double totalFaceWeightSeverityDelta = 0.0;
    std::ptrdiff_t hardShortFaceDelta = 0;
    double maximumShortFaceSeverityDelta = 0.0;
    double totalShortFaceSeverityDelta = 0.0;
    double minInteriorAngleDelta = 0.0;
    double maxNonOrthogonalityDelta = 0.0;
    double maxInternalSkewnessDelta = 0.0;
    double maxCellAspectDelta = 0.0;
    std::size_t firstCellId = 0;
    std::size_t secondCellId = 0;
};

[[nodiscard]] PatchLocalTerminationRank2D patchLocalTerminationRank2D(
    const PatchLocalQuality2D& base,
    const PatchLocalQuality2D& candidate,
    std::size_t firstCellId,
    std::size_t secondCellId) noexcept;

[[nodiscard]] bool patchLocalTerminationRankBetter2D(
    const PatchLocalTerminationRank2D& candidate,
    const PatchLocalTerminationRank2D& current) noexcept;

// Deterministic total order used to pick the single winner among candidates
// whose patches, and therefore whose evaluation scopes, differ. Absolute local
// aggregates are not comparable across different scopes, so the ordering is on
// each candidate's improvement over its own base scope first, and falls back to
// absolute values only as a tiebreak.
struct PatchLocalRank2D {
    // Signed: negative means the candidate removed hard short faces.
    std::ptrdiff_t hardShortFaceDelta = 0;
    double maximumShortFaceSeverityDelta = 0.0;
    double totalShortFaceSeverityDelta = 0.0;
    std::ptrdiff_t issueDelta = 0;
    double maximumIssueSeverityDelta = 0.0;
    double totalIssueSeverityDelta = 0.0;
    // Absolute worst values that remain in the candidate scope.
    double maximumShortFaceSeverity = 0.0;
    double maxNonOrthogonalityDeg = 0.0;
    double maxInternalSkewness = 0.0;
    double maxCellAspect = 0.0;
    // Final tiebreak so the order is total regardless of evaluation order.
    std::size_t firstCellId = 0;
    std::size_t secondCellId = 0;
};

namespace detail {

[[nodiscard]] constexpr bool patchLocalRankKeyLess(
    const PatchLocalRank2D& candidate,const PatchLocalRank2D& current) noexcept {
    return std::tuple{
               candidate.hardShortFaceDelta,
               candidate.maximumShortFaceSeverityDelta,
               candidate.totalShortFaceSeverityDelta,
               candidate.issueDelta,
               candidate.maximumIssueSeverityDelta,
               candidate.totalIssueSeverityDelta,
               candidate.maximumShortFaceSeverity,
               candidate.maxNonOrthogonalityDeg,
               candidate.maxInternalSkewness,
               candidate.maxCellAspect,
               candidate.firstCellId,
               candidate.secondCellId}<
           std::tuple{
               current.hardShortFaceDelta,
               current.maximumShortFaceSeverityDelta,
               current.totalShortFaceSeverityDelta,
               current.issueDelta,
               current.maximumIssueSeverityDelta,
               current.totalIssueSeverityDelta,
               current.maximumShortFaceSeverity,
               current.maxNonOrthogonalityDeg,
               current.maxInternalSkewness,
               current.maxCellAspect,
               current.firstCellId,
               current.secondCellId};
}

constexpr PatchLocalRank2D relativeImprovementRank=[] {
    PatchLocalRank2D rank;
    rank.hardShortFaceDelta=-1;
    rank.maximumShortFaceSeverityDelta=-2.0;
    rank.totalShortFaceSeverityDelta=-8.0;
    rank.maximumShortFaceSeverity=8.0;
    rank.firstCellId=10U;
    rank.secondCellId=11U;
    return rank;
}();
constexpr PatchLocalRank2D lowerAbsoluteButSmallerImprovementRank=[] {
    PatchLocalRank2D rank;
    rank.hardShortFaceDelta=-1;
    rank.maximumShortFaceSeverityDelta=-2.0;
    rank.totalShortFaceSeverityDelta=-3.0;
    rank.maximumShortFaceSeverity=0.0;
    rank.firstCellId=1U;
    rank.secondCellId=2U;
    return rank;
}();
static_assert(patchLocalRankKeyLess(
    relativeImprovementRank,lowerAbsoluteButSmallerImprovementRank));

constexpr PatchLocalRank2D lowIdTieRank=[] {
    PatchLocalRank2D rank;
    rank.hardShortFaceDelta=-1;
    rank.totalShortFaceSeverityDelta=-1.0;
    rank.firstCellId=2U;
    rank.secondCellId=3U;
    return rank;
}();
constexpr PatchLocalRank2D highIdTieRank=[] {
    PatchLocalRank2D rank=lowIdTieRank;
    rank.firstCellId=7U;
    rank.secondCellId=8U;
    return rank;
}();
static_assert(patchLocalRankKeyLess(lowIdTieRank,highIdTieRank));

} // namespace detail

[[nodiscard]] PatchLocalRank2D patchLocalRank2D(const PatchLocalQuality2D& base,
                                                const PatchLocalQuality2D& candidate,
                                                std::size_t firstCellId,
                                                std::size_t secondCellId) noexcept;

[[nodiscard]] bool patchLocalRankBetter2D(const PatchLocalRank2D& candidate,
                                          const PatchLocalRank2D& current) noexcept;

struct PatchLocalScope2D {
    std::vector<PatchLocalCell2D> cells;
    // Input-order topology cell ids of the scope entries, so a candidate scope
    // can reuse the same halo without recomputing adjacency.
    std::vector<std::size_t> topologyCellIds;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

// Builds the base scope for a patch: the selected cells plus their one-ring
// neighbour halo, with stable identities from the incidence store and physical
// boundary flags taken from topology rather than inferred from geometry.
[[nodiscard]] PatchLocalScope2D buildPatchLocalScope2D(
    const TopologyMesh2D& topology, const EdgeIncidenceStore2D& incidence,
    const std::vector<std::size_t>& sortedPatchCellIds,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells);

// Exact comparison of two evaluations of the same scope. Counts and extremal
// aggregates must agree bit-for-bit; only the summed severities admit a
// relative tolerance, because summation order is the only permitted difference.
[[nodiscard]] bool patchLocalQualityMatches2D(const PatchLocalQuality2D& lhs,
                                              const PatchLocalQuality2D& rhs) noexcept;

} // namespace cartmesh2d
