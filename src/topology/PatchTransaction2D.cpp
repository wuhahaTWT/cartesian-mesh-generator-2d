#include "cartmesh2d/topology/PatchTransaction2D.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>

namespace cartmesh2d {
namespace {

using SourceIdentity = std::pair<std::uint64_t,std::size_t>;

[[nodiscard]] SourceIdentity sourceIdentity(const TopologyCell2D& cell) {
    return {cell.sourceKey,cell.sourceId};
}
[[nodiscard]] SourceIdentity sourceIdentity(const CutCell2D& cell) {
    return {cell.sourceKey,cell.sourceId};
}
[[nodiscard]] bool samePoint(const Point2D& a,const Point2D& b,
                             const TolerancePolicy& tol) {
    const double scale=std::max({1.0,std::abs(a.x),std::abs(a.y),
                                std::abs(b.x),std::abs(b.y)});
    return squaredNorm(a-b)<=tol.scale(scale)*tol.scale(scale);
}
[[nodiscard]] bool sameSegment(const Point2D& a,const Point2D& b,
                               const Point2D& c,const Point2D& d,
                               const TolerancePolicy& tol) {
    return (samePoint(a,c,tol) && samePoint(b,d,tol)) ||
           (samePoint(a,d,tol) && samePoint(b,c,tol));
}

} // namespace

TopologyPatchTransaction2D prepareTopologyPatchTransaction2D(
    const TopologyMesh2D& topology,const EdgeIncidenceStore2D& incidence,
    std::vector<std::size_t> selectedCellIds) {
    TopologyPatchTransaction2D result;
    result.baseRevision=incidence.revision;
    if (!topology.valid() || !incidence.valid() ||
        incidence.cellHalfEdges.size()!=topology.cells.size()) {
        result.issues.push_back("cannot prepare a patch transaction from invalid topology/incidence");
        return result;
    }
    std::sort(selectedCellIds.begin(),selectedCellIds.end());
    selectedCellIds.erase(std::unique(selectedCellIds.begin(),selectedCellIds.end()),
                          selectedCellIds.end());
    if (selectedCellIds.empty() || selectedCellIds.back()>=topology.cells.size()) {
        result.issues.push_back("patch transaction requires in-range selected cells");
        return result;
    }
    result.selectedCellIds=std::move(selectedCellIds);
    const std::set<std::size_t> selected(result.selectedCellIds.begin(),
                                         result.selectedCellIds.end());
    for (const auto cell:result.selectedCellIds)
        result.originalPatchArea+=topology.cells[cell].geometryArea;

    for (const auto& edge:topology.edges) {
        const bool ownerSelected=selected.contains(edge.owner);
        const bool neighbourSelected=edge.neighbour && selected.contains(*edge.neighbour);
        const bool interface=ownerSelected!=neighbourSelected;
        const bool physical=!edge.neighbour && ownerSelected;
        if (!interface && !physical) continue;
        result.boundaryLocks.push_back({incidence.stableEdgeKeys[edge.id],
            topology.vertices[edge.v0].point,topology.vertices[edge.v1].point,
            edge.patch,physical});
    }
    std::sort(result.boundaryLocks.begin(),result.boundaryLocks.end(),
              [](const auto& lhs,const auto& rhs) {
        return std::tie(lhs.stableEdge,lhs.physicalBoundary,lhs.patch)<
               std::tie(rhs.stableEdge,rhs.physicalBoundary,rhs.patch);
    });
    if (result.boundaryLocks.empty())
        result.issues.push_back("patch transaction has no locked boundary");
    return result;
}

TopologyPatchCommitResult2D evaluateTopologyPatchTransactionOracle2D(
    const TopologyMesh2D& baseTopology,const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<CutCell2D>& baseSourceCells,
    const std::vector<CutCell2D>& replacementCells,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TopologyPatchCommitGate2D& gate,const TolerancePolicy& tol) {
    TopologyPatchCommitResult2D result;
    result.topology=baseTopology;
    result.incidence=baseIncidence;
    result.revision=baseIncidence.revision;
    result.originalPatchArea=transaction.originalPatchArea;
    result.lockedBoundaryEdgeCount=transaction.boundaryLocks.size();
    if (!transaction.valid() || !baseTopology.valid() || !baseIncidence.valid() ||
        transaction.baseRevision!=baseIncidence.revision ||
        !std::isfinite(gate.qualityRank)) {
        result.issues.push_back("stale or invalid topology patch transaction");
        return result;
    }
    if (!gate.hardQualityPass) {
        result.issues.push_back("topology patch candidate failed the authoritative hard-quality gate");
        return result;
    }

    std::set<SourceIdentity> removed;
    for (const auto cellId:transaction.selectedCellIds)
        removed.insert(sourceIdentity(baseTopology.cells[cellId]));
    std::set<SourceIdentity> foundRemoved;
    std::vector<CutCell2D> candidateSources;
    candidateSources.reserve(baseSourceCells.size()-
        std::min(baseSourceCells.size(),removed.size())+replacementCells.size());
    for (const auto& cell:baseSourceCells) {
        if (removed.contains(sourceIdentity(cell))) foundRemoved.insert(sourceIdentity(cell));
        else candidateSources.push_back(cell);
    }
    if (foundRemoved!=removed) {
        result.issues.push_back("base source cells do not exactly cover selected patch identities");
        return result;
    }
    std::set<SourceIdentity> replacementIdentities;
    for (const auto& cell:replacementCells) {
        if (!cell.valid() || cell.kind==CutCellKind::Empty ||
            cell.kind==CutCellKind::Unsupported ||
            !replacementIdentities.insert(sourceIdentity(cell)).second) {
            result.issues.push_back("replacement patch contains invalid or duplicate source cells");
            return result;
        }
        result.candidatePatchArea+=cell.area;
        candidateSources.push_back(cell);
    }
    const double areaScale=std::max({1.0,result.originalPatchArea,
                                    result.candidatePatchArea});
    if (std::abs(result.candidatePatchArea-result.originalPatchArea)>tol.scale(areaScale)) {
        result.issues.push_back("replacement patch does not conserve source area");
        return result;
    }

    auto candidate=buildGlobalTopology(candidateSources,domain,boundary,tol);
    if (!candidate.valid()) {
        result.issues.push_back("replacement patch failed the global topology oracle");
        return result;
    }

    std::set<std::size_t> candidatePatchCells;
    for (const auto& cell:candidate.cells)
        if (replacementIdentities.contains(sourceIdentity(cell)))
            candidatePatchCells.insert(cell.id);
    if (candidatePatchCells.size()!=replacementCells.size()) {
        result.issues.push_back("replacement identities were not retained one-to-one");
        return result;
    }
    for (const auto& lock:transaction.boundaryLocks) {
        bool matched=false;
        for (const auto& edge:candidate.edges) {
            const auto& a=candidate.vertices[edge.v0].point;
            const auto& b=candidate.vertices[edge.v1].point;
            if (!sameSegment(lock.a,lock.b,a,b,tol)) continue;
            const bool ownerSelected=candidatePatchCells.contains(edge.owner);
            const bool neighbourSelected=edge.neighbour &&
                candidatePatchCells.contains(*edge.neighbour);
            if (lock.physicalBoundary) {
                matched=!edge.neighbour && ownerSelected && edge.patch==lock.patch;
            } else matched=edge.neighbour && ownerSelected!=neighbourSelected;
            if (matched) break;
        }
        if (!matched) {
            result.issues.push_back("replacement patch changed a locked boundary edge");
            return result;
        }
    }

    auto candidateIncidence=buildEdgeIncidenceStore2D(candidate,baseIncidence.revision+1U);
    if (!candidateIncidence.valid()) {
        result.issues.push_back("replacement patch failed edge-incidence audit");
        return result;
    }
    result.topology=std::move(candidate);
    result.incidence=std::move(candidateIncidence);
    result.revision=baseIncidence.revision+1U;
    result.accepted=true;
    result.issues.clear();
    return result;
}

} // namespace cartmesh2d
