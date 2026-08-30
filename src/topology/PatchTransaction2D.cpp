#include "cartmesh2d/topology/PatchTransaction2D.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>

namespace cartmesh2d {
namespace {

using SourceIdentity = TopologySourceIdentity2D;

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
struct DirectedDeltaUse {
    StableVertexId2D from = 0;
    StableVertexId2D to = 0;
};

[[nodiscard]] StableEdgeKey2D edgeKey(StableVertexId2D from,
                                      StableVertexId2D to) {
    const auto endpoints=std::minmax(from,to);
    return {endpoints.first,endpoints.second};
}

void sortEdgeUses(std::vector<RevisionedEdgeUse2D>& uses) {
    std::sort(uses.begin(),uses.end(),[](const auto& lhs,const auto& rhs) {
        return std::tie(lhs.source,lhs.from,lhs.to)<
               std::tie(rhs.source,rhs.from,rhs.to);
    });
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

TopologyDelta2D buildTopologyDelta2D(
    const TopologyMesh2D& baseTopology,const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<TopologyReplacementCell2D>& replacementCells,
    const TolerancePolicy& tol) {
    TopologyDelta2D result;
    (void)tol;
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
    std::set<StableEdgeKey2D> affectedStableEdges;
    for (const auto cell:transaction.selectedCellIds)
        for (const auto edge:baseTopology.cells[cell].edges)
            affectedStableEdges.insert(baseIncidence.stableEdgeKeys[edge]);

    std::map<StableVertexId2D,TopologyDeltaVertex2D> identifiedVertices;
    for (const auto& replacement:replacementCells) {
        const auto& cell=replacement.cell;
        if (!replacement.valid()) {
            result.issues.push_back("topology delta replacement cell is invalid");
            return result;
        }
        result.addedSources.push_back(sourceIdentity(cell));
        result.area+=cell.area;
        for (std::size_t i=0;i<replacement.vertices.size();++i) {
            const auto& vertex=replacement.vertices[i];
            const auto& point=cell.fluidPolygon.vertices[i];
            if (vertex.point.x!=point.x || vertex.point.y!=point.y) {
                result.issues.push_back("replacement stable identity point is not its polygon vertex");
                return result;
            }
            const auto found=identifiedVertices.find(vertex.stableId);
            if (found!=identifiedVertices.end() &&
                (found->second.point.x!=point.x || found->second.point.y!=point.y ||
                 found->second.constructionKey!=vertex.constructionKey ||
                 found->second.existedAtBaseRevision!=vertex.existedAtBaseRevision)) {
                result.issues.push_back("replacement stable identity is internally inconsistent");
                return result;
            }
            identifiedVertices.emplace(vertex.stableId,vertex);
        }
    }
    std::sort(result.addedSources.begin(),result.addedSources.end());
    if (std::adjacent_find(result.addedSources.begin(),result.addedSources.end())!=
        result.addedSources.end()) {
        result.issues.push_back("topology delta has duplicate replacement source identity");
        return result;
    }
    std::map<StableVertexId2D,std::size_t> denseBaseVertex;
    for (std::size_t vertex=0;vertex<baseIncidence.stableVertexIds.size();++vertex)
        denseBaseVertex.emplace(baseIncidence.stableVertexIds[vertex],vertex);
    for (const auto& [stableId,vertex]:identifiedVertices) {
        const auto base=denseBaseVertex.find(stableId);
        if (vertex.existedAtBaseRevision) {
            if (base==denseBaseVertex.end() ||
                baseTopology.vertices[base->second].point.x!=vertex.point.x ||
                baseTopology.vertices[base->second].point.y!=vertex.point.y) {
                result.issues.push_back("replacement references a stale existing stable vertex");
                return result;
            }
        } else if (base!=denseBaseVertex.end() ||
                   vertex.constructionKey.kind!=StableVertexKeyKind2D::PatchGenerated) {
            result.issues.push_back("new replacement vertex lacks a unique typed patch identity");
            return result;
        }
        result.vertices.push_back(vertex);
    }

    std::map<StableEdgeKey2D,std::vector<DirectedDeltaUse>> uses;
    for (const auto& replacement:replacementCells) {
        const auto& cell=replacement.cell;
        const auto& polygon=cell.fluidPolygon.vertices;
        TopologyDeltaCell2D deltaCell;
        deltaCell.source=sourceIdentity(cell);
        deltaCell.area=cell.area;
        for (std::size_t edge=0;edge<polygon.size();++edge) {
            const auto from=replacement.vertices[edge].stableId;
            const auto to=replacement.vertices[(edge+1U)%polygon.size()].stableId;
            if (from==to) {
                result.issues.push_back("topology delta could not assign distinct stable endpoints");
                return result;
            }
            deltaCell.vertices.push_back(from);
            uses[edgeKey(from,to)].push_back({from,to});
        }
        result.cells.push_back(std::move(deltaCell));
    }
    std::sort(result.cells.begin(),result.cells.end(),[](const auto& lhs,const auto& rhs) {
        return lhs.source<rhs.source;
    });

    std::set<std::size_t> matchedLocks;
    for (const auto& [key,incidences]:uses) {
        affectedStableEdges.insert(key);
        TopologyDeltaEdge2D edge{key,incidences.size(),false};
        if (incidences.size()==1U) {
            for (std::size_t lockId=0;lockId<transaction.boundaryLocks.size();++lockId) {
                const auto& lock=transaction.boundaryLocks[lockId];
                if (key==lock.stableEdge) {
                    edge.boundaryLocked=true;
                    edge.patch=lock.patch;
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
    result.affectedStableEdges.assign(affectedStableEdges.begin(),affectedStableEdges.end());
    return result;
}

RevisionedTopology2D buildRevisionedTopology2D(
    const TopologyMesh2D& topology,const EdgeIncidenceStore2D& incidence) {
    RevisionedTopology2D result;
    result.revision=incidence.revision;
    if (!topology.valid() || !incidence.valid() ||
        incidence.stableVertexIds.size()!=topology.vertices.size()) {
        result.issues.push_back("cannot initialize revisioned topology from invalid input");
        return result;
    }
    for (std::size_t vertex=0;vertex<topology.vertices.size();++vertex) {
        RevisionedVertex2D record;
        record.point=topology.vertices[vertex].point;
        record.constructionKey={StableVertexKeyKind2D::LegacyCanonical,
                                incidence.stableVertexIds[vertex],0U,0U,0U};
        if (!result.vertices.emplace(incidence.stableVertexIds[vertex],record).second) {
            result.issues.push_back("revisioned topology has duplicate stable vertex id");
            return result;
        }
    }
    for (const auto& cell:topology.cells) {
        TopologyDeltaCell2D record;
        record.source=sourceIdentity(cell);
        record.area=cell.geometryArea;
        for (const auto vertex:cell.vertices)
            record.vertices.push_back(incidence.stableVertexIds.at(vertex));
        if (!result.cells.emplace(record.source,record).second) {
            result.issues.push_back("revisioned topology has duplicate source identity");
            return result;
        }
        for (std::size_t edge=0;edge<record.vertices.size();++edge) {
            const auto from=record.vertices[edge];
            const auto to=record.vertices[(edge+1U)%record.vertices.size()];
            auto& vertex=result.vertices.at(from);
            ++vertex.referenceCount;
            vertex.state=RevisionedVertexState2D::Active;
            result.edgeIncidences[edgeKey(from,to)].push_back({record.source,from,to});
        }
    }
    for (const auto& edge:topology.edges)
        result.edgePatches[incidence.stableEdgeKeys[edge.id]]=edge.patch;
    for (auto& [key,uses]:result.edgeIncidences) {
        (void)key;
        sortEdgeUses(uses);
        if (uses.empty() || uses.size()>2U ||
            (uses.size()==2U && (uses[0].from!=uses[1].to ||
                                 uses[0].to!=uses[1].from))) {
            result.issues.push_back("revisioned topology edge incidence is invalid");
            return result;
        }
    }
    return result;
}

RevisionedTopology2D applyTopologyDelta2D(
    const RevisionedTopology2D& base,const TopologyDelta2D& delta) {
    RevisionedTopology2D result=base;
    if (!base.valid() || !delta.valid() || base.revision!=delta.baseRevision ||
        delta.candidateRevision!=delta.baseRevision+1U) {
        result.issues.push_back("cannot apply invalid or stale topology delta");
        return result;
    }
    std::set<StableEdgeKey2D> affectedEdges;
    for (const auto& source:delta.removedSources) {
        const auto found=result.cells.find(source);
        if (found==result.cells.end()) {
            result.issues.push_back("topology delta removes an unknown source");
            return result;
        }
        const auto& vertices=found->second.vertices;
        for (std::size_t edge=0;edge<vertices.size();++edge) {
            auto& vertex=result.vertices.at(vertices[edge]);
            if (vertex.referenceCount==0U) {
                result.issues.push_back("removed source underflowed stable vertex references");
                return result;
            }
            --vertex.referenceCount;
            const auto key=edgeKey(vertices[edge],vertices[(edge+1U)%vertices.size()]);
            affectedEdges.insert(key);
            auto edgeFound=result.edgeIncidences.find(key);
            if (edgeFound==result.edgeIncidences.end()) {
                result.issues.push_back("removed source references a missing stable edge");
                return result;
            }
            auto& uses=edgeFound->second;
            uses.erase(std::remove_if(uses.begin(),uses.end(),[&](const auto& use) {
                return use.source==source;
            }),uses.end());
        }
        result.cells.erase(found);
    }
    for (const auto& vertex:delta.vertices) {
        const auto found=result.vertices.find(vertex.stableId);
        if (vertex.existedAtBaseRevision) {
            if (found==result.vertices.end() || found->second.point.x!=vertex.point.x ||
                found->second.point.y!=vertex.point.y) {
                result.issues.push_back("topology delta changed an existing stable vertex");
                return result;
            }
        } else if (found!=result.vertices.end() ||
                   !result.vertices.emplace(vertex.stableId,RevisionedVertex2D{
                       vertex.point,vertex.constructionKey,0U,
                       RevisionedVertexState2D::Tombstone}).second) {
            result.issues.push_back("topology delta reused a stable id for a new vertex");
            return result;
        }
    }
    for (const auto& cell:delta.cells) {
        if (!result.cells.emplace(cell.source,cell).second) {
            result.issues.push_back("topology delta adds an existing source identity");
            return result;
        }
        for (std::size_t edge=0;edge<cell.vertices.size();++edge) {
            const auto from=cell.vertices[edge];
            const auto to=cell.vertices[(edge+1U)%cell.vertices.size()];
            if (!result.vertices.contains(from) || !result.vertices.contains(to) || from==to) {
                result.issues.push_back("topology delta cell references an invalid stable vertex");
                return result;
            }
            auto& vertex=result.vertices.at(from);
            ++vertex.referenceCount;
            vertex.state=RevisionedVertexState2D::Active;
            const auto key=edgeKey(from,to);
            affectedEdges.insert(key);
            result.edgeIncidences[key].push_back({cell.source,from,to});
        }
    }
    for (const auto& key:affectedEdges) {
        auto found=result.edgeIncidences.find(key);
        if (found==result.edgeIncidences.end()) continue;
        auto& uses=found->second;
        if (uses.empty()) {
            result.edgeIncidences.erase(found);
            result.edgePatches.erase(key);
            continue;
        }
        sortEdgeUses(uses);
        if (uses.size()>2U || (uses.size()==2U &&
            (uses[0].from!=uses[1].to || uses[0].to!=uses[1].from))) {
            result.issues.push_back("applied topology delta created invalid edge incidence");
            return result;
        }
        const auto deltaEdge=std::find_if(delta.edges.begin(),delta.edges.end(),
            [&](const auto& edge) { return edge.stableEdge==key; });
        if (deltaEdge!=delta.edges.end()) result.edgePatches[key]=deltaEdge->patch;
    }
    for (auto& [id,vertex]:result.vertices) {
        (void)id;
        vertex.state=vertex.referenceCount>0U
            ?RevisionedVertexState2D::Active:RevisionedVertexState2D::Tombstone;
    }
    result.revision=delta.candidateRevision;
    return result;
}

bool revisionedTopologyMatchesOracle2D(
    const RevisionedTopology2D& revisioned,const TopologyMesh2D& oracle,
    const TolerancePolicy& tol) {
    if (!revisioned.valid() || !oracle.valid() ||
        revisioned.cells.size()!=oracle.cells.size() ||
        revisioned.edgeIncidences.size()!=oracle.edges.size() ||
        revisioned.edgePatches.size()!=oracle.edges.size()) return false;
    std::vector<StableVertexId2D> active;
    for (const auto& [id,vertex]:revisioned.vertices)
        if (vertex.state==RevisionedVertexState2D::Active) {
            if (vertex.referenceCount==0U) return false;
            active.push_back(id);
        } else if (vertex.referenceCount!=0U) return false;
    if (active.size()!=oracle.vertices.size()) return false;
    std::vector<StableVertexId2D> oracleStable(oracle.vertices.size());
    std::set<StableVertexId2D> matched;
    for (std::size_t dense=0;dense<oracle.vertices.size();++dense) {
        std::optional<StableVertexId2D> id;
        for (const auto candidate:active) {
            if (matched.contains(candidate)) continue;
            if (samePoint(revisioned.vertices.at(candidate).point,
                          oracle.vertices[dense].point,tol)) {
                if (id) return false;
                id=candidate;
            }
        }
        if (!id) return false;
        oracleStable[dense]=*id;
        matched.insert(*id);
    }
    const auto cyclicIdsEqual=[](const std::vector<StableVertexId2D>& stable,
                                 const std::vector<StableVertexId2D>& dense) {
        if (stable.size()!=dense.size() || stable.empty()) return false;
        for (std::size_t offset=0;offset<dense.size();++offset) {
            bool same=true;
            for (std::size_t i=0;i<stable.size();++i)
                if (stable[i]!=dense[(i+offset)%dense.size()]) { same=false;break; }
            if (same) return true;
        }
        return false;
    };
    double revisionedArea=0.0,oracleArea=0.0;
    for (const auto& cell:oracle.cells) {
        const auto found=revisioned.cells.find(sourceIdentity(cell));
        std::vector<StableVertexId2D> loop;
        for (const auto vertex:cell.vertices) loop.push_back(oracleStable.at(vertex));
        if (found==revisioned.cells.end() ||
            std::abs(found->second.area-cell.geometryArea)>
                tol.scale(std::max({1.0,found->second.area,cell.geometryArea})) ||
            !cyclicIdsEqual(found->second.vertices,loop)) return false;
        revisionedArea+=found->second.area;
        oracleArea+=cell.geometryArea;
    }
    if (std::abs(revisionedArea-oracleArea)>
        tol.scale(std::max({1.0,revisionedArea,oracleArea}))) return false;
    for (const auto& edge:oracle.edges) {
        const auto key=edgeKey(oracleStable.at(edge.v0),oracleStable.at(edge.v1));
        const auto uses=revisioned.edgeIncidences.find(key);
        const auto patch=revisioned.edgePatches.find(key);
        if (uses==revisioned.edgeIncidences.end() ||
            patch==revisioned.edgePatches.end() || patch->second!=edge.patch ||
            uses->second.size()!=(edge.neighbour?2U:1U)) return false;
        std::set<SourceIdentity> expected{sourceIdentity(oracle.cells.at(edge.owner))};
        if (edge.neighbour) expected.insert(sourceIdentity(oracle.cells.at(*edge.neighbour)));
        std::set<SourceIdentity> actual;
        for (const auto& use:uses->second) {
            if (edgeKey(use.from,use.to)!=key) return false;
            actual.insert(use.source);
        }
        if (actual!=expected) return false;
    }
    return true;
}

TopologyPatchCommitResult2D evaluateTopologyPatchTransactionOracle2D(
    const TopologyMesh2D& baseTopology,const EdgeIncidenceStore2D& baseIncidence,
    const TopologyPatchTransaction2D& transaction,
    const std::vector<CutCell2D>& baseSourceCells,
    const std::vector<TopologyReplacementCell2D>& replacementCells,
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
    result.profile.removedCellCount=result.delta.removedSources.size();
    result.profile.addedCellCount=result.delta.cells.size();
    result.profile.localEdgeCount=result.delta.edges.size();
    result.profile.affectedEdgeCount=result.delta.affectedStableEdges.size();
    result.profile.addedStableVertexCount=static_cast<std::size_t>(std::count_if(
        result.delta.vertices.begin(),result.delta.vertices.end(),
        [](const auto& vertex) { return !vertex.existedAtBaseRevision; }));
    result.profile.retainedStableVertexCount=result.delta.vertices.size()-
        result.profile.addedStableVertexCount;
    result.profile.cacheInvalidatedEdgeCount=result.profile.affectedEdgeCount;
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
    for (const auto& replacement:replacementCells) {
        const auto& cell=replacement.cell;
        if (!replacement.valid() ||
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
    ++result.profile.globalOracleBuildCount;
    auto candidate=buildGlobalTopology(candidateSources,domain,boundary,tol,
                                       baseTopology.constructionRegistry);
    if (!candidate.valid()) {
        result.issues.push_back("replacement patch failed the global topology oracle");
        return result;
    }
    const auto revisionedBase=buildRevisionedTopology2D(baseTopology,baseIncidence);
    result.revisionedTopology=applyTopologyDelta2D(revisionedBase,result.delta);
    result.deltaMatchesGlobalOracle=revisionedTopologyMatchesOracle2D(
        result.revisionedTopology,candidate,tol);
    if (!result.revisionedTopology.valid() || !result.deltaMatchesGlobalOracle) {
        result.issues.push_back("patch-local topology delta disagrees with global oracle");
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
        const auto uses=result.revisionedTopology.edgeIncidences.find(lock.stableEdge);
        const auto patch=result.revisionedTopology.edgePatches.find(lock.stableEdge);
        bool matched=uses!=result.revisionedTopology.edgeIncidences.end() &&
                     patch!=result.revisionedTopology.edgePatches.end() &&
                     patch->second==lock.patch;
        std::size_t selectedUses=0U;
        if (matched) for (const auto& use:uses->second)
            if (replacementIdentities.contains(use.source)) ++selectedUses;
        if (lock.physicalBoundary)
            matched=matched && uses->second.size()==1U && selectedUses==1U;
        else matched=matched && uses->second.size()==2U && selectedUses==1U;
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
