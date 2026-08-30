#pragma once

#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

struct SolverTopologyProfile2D {
    double initialPartitionSeconds = 0.0;
    double buildGlobalTopologySeconds = 0.0;
    double fullQualitySeconds = 0.0;
    double candidateGenerationSeconds = 0.0;
    double candidatePolygonWorkSeconds = 0.0;
    double candidateGlobalRebuildSeconds = 0.0;
    double candidateQualitySeconds = 0.0;
    double sourceRepairSeconds = 0.0;
    double finalRepartitionSeconds = 0.0;
    std::size_t globalTopologyRebuildCalls = 0;
    std::size_t fullQualityCalls = 0;
    std::size_t candidateTopologyCount = 0;
    std::size_t candidateQualityEvaluationCount = 0;
    std::size_t sourceCandidatePairs = 0;
    std::size_t repartitionCandidatePairs = 0;
    std::size_t candidateSplits = 0;
    std::size_t repairPatchCount = 0;
    std::size_t acceptedTopologyCommitCount = 0;
    std::size_t maximumQualityIssueCount = 0;
    std::size_t acceptedSourceRepairs = 0;
    std::size_t acceptedRepartitions = 0;
    std::size_t sourceRepairIterations = 0;
    std::size_t repartitionIterations = 0;
};

struct SolverTopologyResult2D {
    TopologyMesh2D topology;
    std::size_t inputCellCount = 0;
    std::size_t outputCellCount = 0;
    std::size_t partitionedCellCount = 0;
    std::size_t qualityAgglomeratedSourceCellCount = 0;
    std::size_t qualityRepartitionCount = 0;
    std::vector<bool> immutableOutputCells;
    SolverTopologyProfile2D profile;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid() && outputCellCount>0;
    }
};

// Optional solver-repair constraints. An immutable input cell is retained as
// one cell, is never merged or repartitioned, and may acquire only collinear
// common-partition vertices on its edges. H4 uses this for fixed layer quads.
struct SolverTopologyConstraints2D {
    std::vector<bool> immutableInputCells;
    // Retained through the initial convex-partition pass to avoid collinear
    // common-partition triangulation. Unlike immutable cells, these may be
    // merged or repartitioned by later quality repair.
    std::vector<bool> preserveInputCells;
    // Optional unsplit source geometry for preserved/immutable cells. Global
    // topology may contain common-partition vertices that are not geometric
    // corners; retaining the original polygon prevents repair from treating
    // those collinear interface points as mandatory partition vertices.
    std::vector<std::optional<Polygon2D>> inputPolygonOverrides;
};

struct SolverLocalRepartitionResult2D {
    TopologyMesh2D topology;
    std::vector<bool> immutableCells;
    std::size_t repartitionCount = 0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid();
    }
};

struct SolverShortFaceRepairResult2D {
    TopologyMesh2D topology;
    std::vector<bool> immutableCells;
    bool applicable = false;
    bool accepted = false;
    std::size_t candidateCount = 0;
    // Candidates whose patch-local topology delta was valid and therefore
    // reached patch-local quality evaluation.
    std::size_t localCandidateCount = 0;
    std::size_t localQualityEvaluationCount = 0;
    // Must remain zero: candidate selection may not build global topology.
    std::size_t candidateGlobalTopologyBuildCount = 0;
    // Must remain zero: candidate selection may not run the full global
    // solver-quality evaluator.
    std::size_t candidateFullGlobalQualityEvaluationCount = 0;
    // At most one per transaction, for the selected winner only.
    std::size_t globalOracleBuildCount = 0;
    // Authoritative full-mesh quality evaluations: one for the pre-transaction
    // baseline and one for the winner. Never inside candidate selection.
    std::size_t authoritativeFullQualityEvaluationCount = 0;
    std::size_t hardFaceCountBefore = 0;
    std::size_t hardFaceCountAfter = 0;
    double minimumFaceOverLocalHBefore = 0.0;
    double minimumFaceOverLocalHAfter = 0.0;
    bool patchOutsideStableIdsUnchanged = false;
    bool localDeltaMatchesGlobalOracle = false;
    // True when the patch-local accept verdict for the winner was confirmed by
    // the authoritative global evaluation. False means the transaction was
    // rejected fail-closed without trying further candidates.
    bool localWinnerMatchesGlobalAuthority = false;
    double repairSeconds = 0.0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid();
    }
};

// Returns the input-order indices of a deterministic greedy independent set.
// Each halo must contain sorted unique cell/source IDs. A shared ID is a
// repair conflict and prevents both patches from being committed together.
[[nodiscard]] std::vector<std::size_t> selectIndependentSolverRepairPatches2D(
    const std::vector<std::vector<std::size_t>>& halos);

// Rebuilds and evaluates a closed cell patch with the same authoritative
// quality evaluator used by the full topology. Intended for local/full
// consistency tests and repair diagnostics.
[[nodiscard]] SolverQualityReport2D evaluateSolverQualityPatch2D(
    const TopologyMesh2D& topology,
    const std::vector<std::size_t>& sortedCellIds,
    const TolerancePolicy& tol = {});

// Replaces two adjacent solver cells by an area-identical alternative convex
// two-piece partition only when the complete solver-quality score improves.
[[nodiscard]] SolverLocalRepartitionResult2D repartitionSolverTopologyByQuality2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

// Retained H2 exhaustive reference for small H3 equivalence tests only.
[[nodiscard]] SolverLocalRepartitionResult2D
repartitionSolverTopologyByQualitySequentialReference2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

// Performs at most one deterministic short-face transaction. Every candidate is
// scored only by patch-local topology and patch-local quality over the affected
// patch plus its one-ring halo: candidate selection calls neither
// buildGlobalTopology() nor evaluateSolverQuality2D(). The single ranked winner
// then performs one global oracle/materialization and one authoritative full
// quality evaluation. A disagreement between the local verdict and that
// authoritative result fails the transaction closed; no further candidate is
// attempted.
[[nodiscard]] SolverShortFaceRepairResult2D repairSolverShortFaces2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::vector<bool>& immutableCells,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells,
    double minimumFaceOverLocalH,
    const TolerancePolicy& tol = {});

// Retains strictly convex cells and deterministically ear-clips cells with
// reflex or collinear vertices.  The latter occur at conformal 2:1 Cartesian
// transitions and are reported as concave polyhedra by OpenFOAM when extruded.
[[nodiscard]] SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const SolverTopologyConstraints2D& constraints,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
