#pragma once

#include "cartmesh2d/topology/EdgeIncidence2D.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace cartmesh2d {

struct PatchBoundaryLock2D {
    StableEdgeKey2D stableEdge;
    Point2D a;
    Point2D b;
    BoundaryPatch2D patch = BoundaryPatch2D::None;
    bool physicalBoundary = false;
};

struct TopologyPatchTransaction2D {
    std::uint64_t baseRevision = 0;
    std::vector<std::size_t> selectedCellIds;
    std::vector<PatchBoundaryLock2D> boundaryLocks;
    double originalPatchArea = 0.0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

struct TopologyPatchCommitGate2D {
    // Area and boundary-feature compatibility are evaluated internally.
    // The construction caller supplies the authoritative typed quality result.
    bool hardQualityPass = false;
    double qualityRank = 0.0;
};

struct TopologyPatchCommitResult2D {
    TopologyMesh2D topology;
    EdgeIncidenceStore2D incidence;
    bool accepted = false;
    std::uint64_t revision = 0;
    double originalPatchArea = 0.0;
    double candidatePatchArea = 0.0;
    std::size_t lockedBoundaryEdgeCount = 0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid() && incidence.valid();
    }
};

[[nodiscard]] TopologyPatchTransaction2D prepareTopologyPatchTransaction2D(
    const TopologyMesh2D& topology, const EdgeIncidenceStore2D& incidence,
    std::vector<std::size_t> selectedCellIds);

// Current R1-D oracle path: rebuilds the candidate globally, then proves the
// patch boundary lock, exact source replacement, area conservation and hard
// quality gate before advancing the revision. Rejection returns the unchanged
// base topology and incidence. A later incremental implementation must remain
// byte/topology equivalent to this function.
[[nodiscard]] TopologyPatchCommitResult2D evaluateTopologyPatchTransactionOracle2D(
    const TopologyMesh2D& baseTopology,
    const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<CutCell2D>& baseSourceCells,
    const std::vector<CutCell2D>& replacementCells,
    const Domain2D& domain, const BoundaryRegion2D& boundary,
    const TopologyPatchCommitGate2D& gate,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
