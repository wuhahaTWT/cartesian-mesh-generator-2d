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

struct DirectedDeltaUse {
    StableVertexId2D from = 0;
    StableVertexId2D to = 0;
};

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

TopologyDelta2D buildTopologyDelta2D(
    const TopologyMesh2D& baseTopology,const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<CutCell2D>& replacementCells,const TolerancePolicy& tol) {
    TopologyDelta2D result;
    result.baseRevision=transaction.baseRevision;
    result.candidateRevision=transaction.baseRevision+1U;
    if (!baseTopology.valid() || !baseIncidence.valid() || !transaction.valid() ||
        transaction.baseRevision!=baseIncidence.revision) {
        result.issues.push_back("cannot build topology delta from invalid or stale base");
        return result;
    }
    for (const auto cell:transaction.selectedCellIds)
        result.removedSources.push_back(sourceIdentity(baseTopology.cells[cell]));
    std::sort(result.removedSources.begin(),result.removedSources.end());

    std::vector<Point2D> newPoints;
    for (const auto& cell:replacementCells) {
        if (!cell.valid() || cell.kind==CutCellKind::Empty ||
            cell.kind==CutCellKind::Unsupported) {
            result.issues.push_back("topology delta replacement cell is invalid");
            return result;
        }
        result.addedSources.push_back(sourceIdentity(cell));
        result.area+=cell.area;
        for (const auto& point:cell.fluidPolygon.vertices) {
            bool found=false;
            for (const auto& base:baseTopology.vertices)
                if (samePoint(point,base.point,tol)) { found=true;break; }
            if (!found) newPoints.push_back(point);
        }
    }
    std::sort(result.addedSources.begin(),result.addedSources.end());
    if (std::adjacent_find(result.addedSources.begin(),result.addedSources.end())!=
        result.addedSources.end()) {
        result.issues.push_back("topology delta has duplicate replacement source identity");
        return result;
    }
    std::sort(newPoints.begin(),newPoints.end(),[](const auto& lhs,const auto& rhs) {
        return std::tie(lhs.x,lhs.y)<std::tie(rhs.x,rhs.y);
    });
    newPoints.erase(std::unique(newPoints.begin(),newPoints.end(),[&](const auto& lhs,
                                                                     const auto& rhs) {
        return samePoint(lhs,rhs,tol);
    }),newPoints.end());
    StableVertexId2D nextId=0;
    for (const auto id:baseIncidence.stableVertexIds)
        nextId=std::max(nextId,id+1U);
    for (std::size_t ordinal=0;ordinal<newPoints.size();++ordinal) {
        result.vertices.push_back({nextId+ordinal,newPoints[ordinal],false});
    }

    const auto stableVertex=[&](const Point2D& point)->std::optional<StableVertexId2D> {
        for (std::size_t vertex=0;vertex<baseTopology.vertices.size();++vertex)
            if (samePoint(point,baseTopology.vertices[vertex].point,tol))
                return baseIncidence.stableVertexIds[vertex];
        for (std::size_t vertex=0;vertex<result.vertices.size();++vertex)
            if (samePoint(point,result.vertices[vertex].point,tol))
                return result.vertices[vertex].stableId;
        return std::nullopt;
    };
    for (std::size_t vertex=0;vertex<baseTopology.vertices.size();++vertex) {
        bool used=false;
        for (const auto& cell:replacementCells) for (const auto& point:cell.fluidPolygon.vertices)
            if (samePoint(point,baseTopology.vertices[vertex].point,tol)) used=true;
        if (used) result.vertices.push_back({baseIncidence.stableVertexIds[vertex],
                                             baseTopology.vertices[vertex].point,true});
    }
    std::sort(result.vertices.begin(),result.vertices.end(),[](const auto& lhs,const auto& rhs) {
        return lhs.stableId<rhs.stableId;
    });

    std::map<StableEdgeKey2D,std::vector<DirectedDeltaUse>> uses;
    for (const auto& cell:replacementCells) {
        const auto& polygon=cell.fluidPolygon.vertices;
        for (std::size_t edge=0;edge<polygon.size();++edge) {
            const auto from=stableVertex(polygon[edge]);
            const auto to=stableVertex(polygon[(edge+1U)%polygon.size()]);
            if (!from || !to || *from==*to) {
                result.issues.push_back("topology delta could not assign distinct stable endpoints");
                return result;
            }
            const auto endpoints=std::minmax(*from,*to);
            uses[{endpoints.first,endpoints.second}].push_back({*from,*to});
        }
    }

    std::set<std::size_t> matchedLocks;
    for (const auto& [key,incidences]:uses) {
        TopologyDeltaEdge2D edge{key,incidences.size(),false};
        if (incidences.size()==1U) {
            const auto& use=incidences.front();
            const auto pointFor=[&](StableVertexId2D id)->std::optional<Point2D> {
                for (const auto& vertex:result.vertices)
                    if (vertex.stableId==id) return vertex.point;
                return std::nullopt;
            };
            const auto a=pointFor(use.from),b=pointFor(use.to);
            for (std::size_t lockId=0;lockId<transaction.boundaryLocks.size();++lockId) {
                const auto& lock=transaction.boundaryLocks[lockId];
                if (a && b && sameSegment(*a,*b,lock.a,lock.b,tol)) {
                    edge.boundaryLocked=true;
                    matchedLocks.insert(lockId);
                    break;
                }
            }
            if (!edge.boundaryLocked)
                result.issues.push_back("topology delta has an unlocked one-incidence edge");
            else ++result.lockedBoundaryEdgeCount;
        } else if (incidences.size()==2U) {
            if (incidences[0].from!=incidences[1].to ||
                incidences[0].to!=incidences[1].from)
                result.issues.push_back("topology delta internal incidences are not opposite");
            else ++result.internalEdgeCount;
        } else result.issues.push_back("topology delta edge is non-manifold");
        result.edges.push_back(edge);
    }
    if (matchedLocks.size()!=transaction.boundaryLocks.size())
        result.issues.push_back("topology delta does not preserve every locked boundary edge");
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

    result.delta=buildTopologyDelta2D(baseTopology,baseIncidence,transaction,
                                      replacementCells,tol);
    if (!result.delta.valid()) {
        result.issues=result.delta.issues;
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

    ++result.globalOracleBuildCount;
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
