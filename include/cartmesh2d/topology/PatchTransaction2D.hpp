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

struct TopologyDeltaVertex2D {
    StableVertexId2D stableId = 0;
    Point2D point;
    bool existedAtBaseRevision = false;
};

struct TopologyDeltaEdge2D {
    StableEdgeKey2D stableEdge;
    std::size_t incidenceCount = 0;
    bool boundaryLocked = false;
};

struct TopologyDelta2D {
    std::uint64_t baseRevision = 0;
    std::uint64_t candidateRevision = 0;
    std::vector<std::pair<std::uint64_t,std::size_t>> removedSources;
    std::vector<std::pair<std::uint64_t,std::size_t>> addedSources;
    std::vector<TopologyDeltaVertex2D> vertices;
    std::vector<TopologyDeltaEdge2D> edges;
    std::size_t internalEdgeCount = 0;
    std::size_t lockedBoundaryEdgeCount = 0;
    double area = 0.0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

struct TopologyPatchCommitResult2D {
    TopologyMesh2D topology;
    EdgeIncidenceStore2D incidence;
    bool accepted = false;
    std::uint64_t revision = 0;
    double originalPatchArea = 0.0;
    double candidatePatchArea = 0.0;
    std::size_t lockedBoundaryEdgeCount = 0;
    TopologyDelta2D delta;
    std::size_t globalOracleBuildCount = 0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid() && incidence.valid();
    }
};

[[nodiscard]] TopologyPatchTransaction2D prepareTopologyPatchTransaction2D(
    const TopologyMesh2D& topology, const EdgeIncidenceStore2D& incidence,
    std::vector<std::size_t> selectedCellIds);

// Builds and audits only the replacement patch. One-incidence edges must match
// an exact boundary lock; all other edges must have two opposite incidences.
// New interior vertices receive deterministic IDs after the maximum base ID.
[[nodiscard]] TopologyDelta2D buildTopologyDelta2D(
    const TopologyMesh2D& baseTopology,
    const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<CutCell2D>& replacementCells,
    const TolerancePolicy& tol = {});

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
