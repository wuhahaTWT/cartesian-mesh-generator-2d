#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cartmesh2d {

struct SolverQualityReport2D;

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
    SolverTopologyProfile2D profile;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid() && outputCellCount>0;
    }
};

struct SolverLocalRepartitionResult2D {
    TopologyMesh2D topology;
    std::size_t repartitionCount = 0;
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

// Retains strictly convex cells and deterministically ear-clips cells with
// reflex or collinear vertices.  The latter occur at conformal 2:1 Cartesian
// transitions and are reported as concave polyhedra by OpenFOAM when extruded.
[[nodiscard]] SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
