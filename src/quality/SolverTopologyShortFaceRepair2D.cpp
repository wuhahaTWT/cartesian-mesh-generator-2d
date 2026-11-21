#include "SolverTopologyDetail2D.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace cartmesh2d {
using namespace solver_topology_detail;

SolverShortFaceRepairResult2D repairSolverShortFaces2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells,double minimumFaceOverLocalH,
    const TolerancePolicy& tol) {
    const auto repairStart=ProfileClock::now();
    SolverShortFaceRepairResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        immutableCells.size()!=topology.cells.size() ||
        localBackgroundH.size()!=topology.cells.size() ||
        ratedCells.size()!=topology.cells.size() || !(minimumFaceOverLocalH>0.0)) {
        result.issues.push_back("revisioned short-face repair requires valid aligned inputs");
        return result;
    }
    struct FaceScore {
        std::size_t hardCount=0;
        double minimum=std::numeric_limits<double>::infinity();
        double maximumSeverity=0.0;
        double totalSeverity=0.0;
    };
    const auto scoreFaces=[&](const TopologyMesh2D& mesh,
                              const std::vector<double>& localH,
                              const std::vector<bool>& rated) {
        FaceScore score;
        for (const auto& edge:mesh.edges) {
            double h=0.0;
            if (rated.at(edge.owner)) h=std::max(h,localH.at(edge.owner));
            if (edge.neighbour && rated.at(*edge.neighbour))
                h=std::max(h,localH.at(*edge.neighbour));
            if (!(h>0.0)) continue;
            const double ratio=std::sqrt(squaredNorm(
                mesh.vertices[edge.v1].point-mesh.vertices[edge.v0].point))/h;
            score.minimum=std::min(score.minimum,ratio);
            if (ratio<minimumFaceOverLocalH) {
                ++score.hardCount;
                const double severity=minimumFaceOverLocalH/
                    (ratio+std::numeric_limits<double>::min());
                score.maximumSeverity=std::max(score.maximumSeverity,severity);
                score.totalSeverity+=severity;
            }
        }
        if (!std::isfinite(score.minimum)) score.minimum=0.0;
        return score;
    };
    const auto better=[](const FaceScore& lhs,const FaceScore& rhs) {
        return std::tie(lhs.hardCount,lhs.maximumSeverity,lhs.totalSeverity) <
               std::tie(rhs.hardCount,rhs.maximumSeverity,rhs.totalSeverity);
    };
    const auto legacyNoWorse=[](const SolverQualityReport2D& candidate,
                                const SolverQualityReport2D& current) {
        if (!candidate.valid()) return false;
        const auto upper=[](double after,double before) {
            return after-before<=1.0e-8*std::max(1.0,std::abs(before));
        };
        const auto lower=[](double after,double before) {
            return before-after<=1.0e-8*std::max(1.0,std::abs(before));
        };
        return upper(candidate.maxNonOrthogonalityDeg,current.maxNonOrthogonalityDeg) &&
               upper(candidate.maxInternalSkewness,current.maxInternalSkewness) &&
               upper(candidate.maxBoundarySkewness,current.maxBoundarySkewness) &&
               upper(candidate.maxConcavityDeg,current.maxConcavityDeg) &&
               upper(candidate.maxCellAspect,current.maxCellAspect) &&
               lower(candidate.minInteriorAngleDeg,current.minInteriorAngleDeg) &&
               lower(candidate.minFaceLength,current.minFaceLength) &&
               lower(candidate.minFaceWeight,current.minFaceWeight) &&
               lower(candidate.minVolumeRatio,current.minVolumeRatio) &&
               lower(candidate.minCompactness,current.minCompactness);
    };
    const auto currentQuality=evaluateSolverQuality2D(topology,{},tol);
    ++result.authoritativeFullQualityEvaluationCount;
    if (!currentQuality.valid()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const auto currentScore=scoreFaces(topology,localBackgroundH,ratedCells);
    result.hardFaceCountBefore=currentScore.hardCount;
    result.hardFaceCountAfter=currentScore.hardCount;
    result.minimumFaceOverLocalHBefore=currentScore.minimum;
    result.minimumFaceOverLocalHAfter=currentScore.minimum;
    if (currentScore.hardCount==0U) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }

    std::vector<std::tuple<double,std::size_t,std::size_t>> shortEdges;
    for (const auto& edge:topology.edges) {
        if (!edge.neighbour) continue;
        double h=0.0;
        if (ratedCells[edge.owner]) h=std::max(h,localBackgroundH[edge.owner]);
        if (ratedCells[*edge.neighbour]) h=std::max(h,localBackgroundH[*edge.neighbour]);
        if (!(h>0.0)) continue;
        const double ratio=std::sqrt(squaredNorm(
            topology.vertices[edge.v1].point-topology.vertices[edge.v0].point))/h;
        if (ratio<minimumFaceOverLocalH)
            shortEdges.emplace_back(ratio,std::min(edge.owner,*edge.neighbour),
                                    std::max(edge.owner,*edge.neighbour));
    }
    std::sort(shortEdges.begin(),shortEdges.end());
    if (shortEdges.empty()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const std::vector<std::size_t> affected{std::get<1>(shortEdges.front()),
                                            std::get<2>(shortEdges.front())};
    result.affectedCells=affected;
    // This R1 checkpoint migrates the narrow-gap Q2-B mode only: the failing
    // atomic face is a common-partition fragment on an immutable layer
    // support. Mutable/mutable sharp-tail patches still belong to the later
    // migration stage and must retain their pre-R1 behavior.
    if (!immutableCells[affected[0]] && !immutableCells[affected[1]]) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    result.applicable=true;
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    if (!immutableCells[affected[0]] && !immutableCells[affected[1]])
        pairs.emplace_back(affected[0],affected[1]);
    for (const auto cell:affected) {
        if (immutableCells[cell]) continue;
        for (const auto edgeId:topology.cells[cell].edges) {
            const auto& edge=topology.edges[edgeId];
            if (!edge.neighbour) continue;
            const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
            if (!immutableCells[other])
                pairs.emplace_back(std::min(cell,other),std::max(cell,other));
        }
    }
    std::sort(pairs.begin(),pairs.end());
    pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());

    std::vector<std::vector<std::size_t>> mutablePatches;
    for (const auto& [first,second]:pairs) mutablePatches.push_back({first,second});
    std::set<std::vector<std::size_t>> seenTriples;
    for (const auto& [first,second]:pairs) {
        for (const auto cell:{first,second}) {
            for (const auto edgeId:topology.cells[cell].edges) {
                const auto& edge=topology.edges[edgeId];
                if (!edge.neighbour) continue;
                const auto third=edge.owner==cell?*edge.neighbour:edge.owner;
                if (third==first || third==second || immutableCells[third]) continue;
                std::vector<std::size_t> key{first,second,third};
                std::sort(key.begin(),key.end());
                if (seenTriples.insert(key).second)
                    mutablePatches.push_back({first,second,third});
            }
        }
    }

    const auto incidence=buildEdgeIncidenceStore2D(topology,0U);
    if (!incidence.valid()) {
        result.issues.push_back("short-face base incidence is invalid");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    std::map<std::pair<double,double>,TopologyDeltaVertex2D> identityByPoint;
    std::map<StableVertexId2D,TopologyDeltaVertex2D> identityById;
    for (std::size_t dense=0;dense<topology.vertices.size();++dense) {
        const auto id=incidence.stableVertexIds[dense];
        const TopologyDeltaVertex2D identity{id,
            {StableVertexKeyKind2D::LegacyCanonical,id,0U,0U,0U},
            topology.vertices[dense].point,true};
        identityByPoint.emplace(
            std::pair{topology.vertices[dense].point.x,topology.vertices[dense].point.y},
            identity);
        identityById.emplace(id,identity);
    }
    const auto sourceCell=[&](std::size_t cellId,const Polygon2D& polygon) {
        CutCell2D cell=makeCell(topology.cells[cellId].sourceId,polygon,tol,
                                topology.cells[cellId].sourceLineage);
        cell.sourceKey=topology.cells[cellId].sourceKey;
        cell.refinementLineageKeys=topology.cells[cellId].refinementLineageKeys;
        if (topology.constructionRegistry) {
            for (const auto& point:polygon.vertices)
                cell.canonicalVertexIds.push_back(
                    identityByPoint.at({point.x,point.y}).stableId);
            cell.canonicalRegistry=topology.constructionRegistry.get();
        }
        return cell;
    };
    std::vector<CutCell2D> baseSources;
    baseSources.reserve(topology.cells.size());
    for (const auto& cell:topology.cells)
        baseSources.push_back(sourceCell(cell.id,topologyCellPolygon(topology,cell.id)));

    // Phase 1: every candidate is generated, validated and ranked using only
    // the affected patch and its one-ring halo. No global topology build and no
    // full global quality evaluation may happen in this loop.
    struct ShortFaceCandidate2D {
        std::size_t first = 0;
        std::size_t second = 0;
        std::vector<std::size_t> mutableCells;
        std::vector<std::size_t> selected;
        TopologyPatchTransaction2D transaction;
        std::vector<TopologyReplacementCell2D> replacements;
        PatchLocalRank2D rank;
    };
    std::optional<ShortFaceCandidate2D> winner;
    std::string localRejectReason;
    // Measured, not asserted: any accidental global rebuild or full-quality
    // evaluation inside the selection loop shows up in these deltas.
    const auto selectionStartGlobalBuilds=globalTopologyBuildCount2D();
    const auto selectionStartFullQuality=solverQualityEvaluationCount2D();
    for (std::size_t patchIndex=0;patchIndex<mutablePatches.size();++patchIndex) {
        // Preserve two-cell winners. A neighbouring third mutable cell is
        // considered only if every pair failed; no layer cell may be merged.
        if (patchIndex==pairs.size() && winner) break;
        const auto& mutablePatch=mutablePatches[patchIndex];
        const auto first=mutablePatch[0],second=mutablePatch[1];
        ++result.candidateCount;
        std::optional<Polygon2D> merged=topologyCellPolygon(topology,first);
        for (std::size_t i=1;i<mutablePatch.size() && merged;++i)
            merged=mergeAdjacentPolygonsSimple(*merged,
                topologyCellPolygon(topology,mutablePatch[i]),tol);
        if (!merged) { ++result.rejectedUnionCount; continue; }
        const auto repairUnion=removeArtificialCollinearVertices(*merged,tol,true);
        if (!strictlyConvex(repairUnion)) {
            ++result.rejectedConvexityCount;
            continue;
        }
        if (underDeterminedBoundaryCell(repairUnion,domain,boundary,tol)) {
            ++result.rejectedBoundaryCount;
            continue;
        }
        std::vector<std::size_t> selected=mutablePatch;
        std::map<std::size_t,Polygon2D> simplifiedImmutable;
        std::set<std::size_t> localImmutable;
        auto incidentCells=mutablePatch;
        incidentCells.insert(incidentCells.end(),affected.begin(),affected.end());
        for (const auto cell:incidentCells) {
            if (immutableCells[cell]) localImmutable.insert(cell);
            for (const auto edgeId:topology.cells[cell].edges) {
                const auto& edge=topology.edges[edgeId];
                if (!edge.neighbour) continue;
                const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
                if (immutableCells[other]) localImmutable.insert(other);
            }
        }
        for (const auto cell:localImmutable) {
            const auto original=topologyCellPolygon(topology,cell);
            const auto simplified=removeArtificialCollinearVertices(original,tol,true);
            if (simplified.vertices.size()!=original.vertices.size()) {
                selected.push_back(cell);
                simplifiedImmutable.emplace(cell,simplified);
            }
        }
        std::sort(selected.begin(),selected.end());
        selected.erase(std::unique(selected.begin(),selected.end()),selected.end());
        const auto transaction=prepareTopologyPatchTransaction2D(topology,incidence,selected);
        if (!transaction.valid()) continue;
        // Simplifying several replacement polygons can remove a fragment on
        // one side while retaining an immutable layer corner on the other.
        // Build one atomic partition for both sides, using existing identities
        // only. Points removed from every replacement are not reintroduced.
        std::set<StableVertexId2D> replacementEndpoints;
        const auto collectEndpoints=[&](const Polygon2D& polygon) {
            for (const auto& point:polygon.vertices)
                replacementEndpoints.insert(identityByPoint.at({point.x,point.y}).stableId);
        };
        collectEndpoints(repairUnion);
        for (const auto& [cell,polygon]:simplifiedImmutable) {
            (void)cell;
            collectEndpoints(polygon);
        }
        const auto identified=[&](std::size_t sourceId,Polygon2D polygon) {
            // Insert locked stable endpoints on any simplified outer patch
            // edge. Identity comes from the lock key; geometry is only used to
            // order already-identified endpoints along that local edge.
            Polygon2D locked;
            for (std::size_t i=0;i<polygon.vertices.size();++i) {
                const auto a=polygon.vertices[i];
                const auto b=polygon.vertices[(i+1U)%polygon.vertices.size()];
                std::vector<std::pair<double,TopologyDeltaVertex2D>> points;
                const auto add=[&](StableVertexId2D id) {
                    const auto found=identityById.find(id);
                    if (found==identityById.end() ||
                        !pointOnSegment(found->second.point,{a,b},tol)) return;
                    const auto direction=b-a;
                    const double denom=squaredNorm(direction);
                    const double t=denom>0.0?dot(found->second.point-a,direction)/denom:0.0;
                    points.emplace_back(t,found->second);
                };
                add(identityByPoint.at({a.x,a.y}).stableId);
                for (const auto id:replacementEndpoints) add(id);
                for (const auto& lock:transaction.boundaryLocks) {
                    add(lock.stableEdge.v0);
                    add(lock.stableEdge.v1);
                }
                std::sort(points.begin(),points.end(),[](const auto& lhs,const auto& rhs) {
                    return std::tie(lhs.first,lhs.second.stableId)<
                           std::tie(rhs.first,rhs.second.stableId);
                });
                points.erase(std::unique(points.begin(),points.end(),[](const auto& lhs,
                                                                        const auto& rhs) {
                    return lhs.second.stableId==rhs.second.stableId;
                }),points.end());
                for (const auto& [t,vertex]:points) {
                    (void)t;
                    if (locked.vertices.empty() || locked.vertices.back().x!=vertex.point.x ||
                        locked.vertices.back().y!=vertex.point.y)
                        locked.vertices.push_back(vertex.point);
                }
            }
            if (locked.vertices.size()>1U &&
                locked.vertices.front().x==locked.vertices.back().x &&
                locked.vertices.front().y==locked.vertices.back().y)
                locked.vertices.pop_back();
            TopologyReplacementCell2D replacement;
            replacement.cell=sourceCell(sourceId,locked);
            for (const auto& point:locked.vertices)
                replacement.vertices.push_back(identityByPoint.at({point.x,point.y}));
            return replacement;
        };
        std::vector<TopologyReplacementCell2D> replacements;
        std::vector<std::pair<double,bool>> replacementMetadata;
        replacements.push_back(identified(first,repairUnion));
        std::vector<std::size_t> lineage;
        double mergedH=0.0;
        bool mergedRated=false;
        for (const auto cell:mutablePatch) {
            lineage=mergedLineage(lineage,topology.cells[cell].sourceLineage);
            mergedH=std::max(mergedH,localBackgroundH[cell]);
            mergedRated=mergedRated || ratedCells[cell];
        }
        replacements.back().cell.sourceLineage=std::move(lineage);
        replacementMetadata.emplace_back(mergedH,mergedRated);
        for (const auto& [cell,polygon]:simplifiedImmutable) {
            replacements.push_back(identified(cell,polygon));
            replacementMetadata.emplace_back(localBackgroundH[cell],ratedCells[cell]);
        }
        const auto localDelta=buildTopologyDelta2D(
            topology,incidence,transaction,replacements,tol);
        if (!localDelta.valid()) {
            if (!localDelta.issues.empty()) localRejectReason=localDelta.issues.front();
            continue;
        }
        ++result.localCandidateCount;

        const auto baseScope=buildPatchLocalScope2D(
            topology,incidence,selected,localBackgroundH,ratedCells);
        if (!baseScope.valid()) {
            localRejectReason=baseScope.issues.front();
            continue;
        }
        std::set<StableEdgeKey2D> physicalLocks;
        for (const auto& lock:transaction.boundaryLocks)
            if (lock.physicalBoundary) physicalLocks.insert(lock.stableEdge);
        std::vector<PatchLocalCell2D> candidateScope;
        for (std::size_t i=0;i<baseScope.cells.size();++i)
            if (!baseScope.cells[i].inPatch) candidateScope.push_back(baseScope.cells[i]);
        for (std::size_t i=0;i<replacements.size();++i) {
            const auto& replacement=replacements[i];
            PatchLocalCell2D entry;
            entry.inPatch=true;
            entry.localBackgroundH=replacementMetadata[i].first;
            entry.ratedForFaceLength=replacementMetadata[i].second;
            entry.polygon=replacement.cell.fluidPolygon;
            for (std::size_t local=0;local<replacement.vertices.size();++local) {
                const auto from=replacement.vertices[local].stableId;
                const auto to=replacement.vertices[
                    (local+1U)%replacement.vertices.size()].stableId;
                const auto endpoints=std::minmax(from,to);
                entry.loop.push_back(from);
                entry.physicalBoundaryFace.push_back(
                    physicalLocks.contains(StableEdgeKey2D{endpoints.first,endpoints.second}));
            }
            candidateScope.push_back(std::move(entry));
        }
        const auto baseLocal=evaluatePatchLocalQuality2D(
            baseScope.cells,minimumFaceOverLocalH,{},tol);
        const auto candidateLocal=evaluatePatchLocalQuality2D(
            candidateScope,minimumFaceOverLocalH,{},tol);
        result.localQualityEvaluationCount+=2U;
        if (!baseLocal.valid() || !candidateLocal.valid()) {
            localRejectReason=baseLocal.valid()?candidateLocal.issues.front()
                                               :baseLocal.issues.front();
            continue;
        }
        // Outside cells and faces are unchanged by the boundary-locked delta.
        // The touched metrics must stay within the measured global baseline,
        // and the local short-face improvement must survive restoring outside
        // faces. The winner still passes the original full-mesh no-worse gate.
        if (!patchLocalQualityWithinGlobalBaseline2D(candidateLocal,baseLocal,currentQuality)) {
            localRejectReason="patch-local candidate exceeds the global quality baseline";
            continue;
        }
        if (!(std::tie(candidateLocal.hardShortFaceCount,
                       candidateLocal.maximumShortFaceSeverity,
                       candidateLocal.totalShortFaceSeverity)<
              std::tie(baseLocal.hardShortFaceCount,
                       baseLocal.maximumShortFaceSeverity,
                       baseLocal.totalShortFaceSeverity))) {
            localRejectReason="patch-local candidate does not reduce the local short-face score";
            continue;
        }
        ShortFaceCandidate2D accepted{first,second,mutablePatch,selected,transaction,replacements,
            patchLocalRank2D(baseLocal,candidateLocal,first,second)};
        if (!winner || patchLocalRankBetter2D(accepted.rank,winner->rank))
            winner=std::move(accepted);
    }
    if (!winner) {
        result.candidateGlobalTopologyBuildCount=
            globalTopologyBuildCount2D()-selectionStartGlobalBuilds;
        result.candidateFullGlobalQualityEvaluationCount=
            solverQualityEvaluationCount2D()-selectionStartFullQuality;
        result.issues.push_back(localRejectReason.empty()
            ?"no patch-local short-face candidate passed the local quality gate"
            :localRejectReason);
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    result.candidateGlobalTopologyBuildCount=
        globalTopologyBuildCount2D()-selectionStartGlobalBuilds;
    result.candidateFullGlobalQualityEvaluationCount=
        solverQualityEvaluationCount2D()-selectionStartFullQuality;

    // Phase 2: exactly one global oracle/materialization and one authoritative
    // full-quality evaluation, for the selected winner only. A disagreement
    // fails the transaction closed rather than retrying other candidates.
    const auto& first=winner->first;
    const auto& selected=winner->selected;
    const auto committed=evaluateTopologyPatchTransactionOracle2D(
        topology,incidence,winner->transaction,baseSources,winner->replacements,
        domain,boundary,{true,0.0},tol);
    result.globalOracleBuildCount+=committed.globalOracleBuildCount;
    if (!committed.accepted || !committed.deltaMatchesGlobalOracle) {
        result.issues.push_back(committed.issues.empty()
            ?"winner patch transaction failed the global oracle"
            :committed.issues.front());
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    std::map<TopologySourceIdentity2D,std::pair<double,bool>> metadata;
    std::set<TopologySourceIdentity2D> immutableIdentities;
    for (const auto& cell:topology.cells) {
        const TopologySourceIdentity2D identity{cell.sourceKey,cell.sourceId};
        metadata[identity]={localBackgroundH[cell.id],ratedCells[cell.id]};
        if (immutableCells[cell.id]) immutableIdentities.insert(identity);
    }
    double mergedH=0.0;
    bool mergedRated=false;
    for (const auto cell:winner->mutableCells) {
        mergedH=std::max(mergedH,localBackgroundH[cell]);
        mergedRated=mergedRated || ratedCells[cell];
    }
    metadata[{topology.cells[first].sourceKey,topology.cells[first].sourceId}]=
        {mergedH,mergedRated};
    std::vector<double> candidateH;
    std::vector<bool> candidateRated,candidateImmutable;
    for (const auto& cell:committed.topology.cells) {
        const auto identity=TopologySourceIdentity2D{cell.sourceKey,cell.sourceId};
        const auto item=metadata.at(identity);
        candidateH.push_back(item.first);
        candidateRated.push_back(item.second);
        candidateImmutable.push_back(immutableIdentities.contains(identity));
    }
    const auto candidateQuality=evaluateSolverQuality2D(committed.topology,{},tol);
    ++result.authoritativeFullQualityEvaluationCount;
    const auto candidateScore=scoreFaces(committed.topology,candidateH,candidateRated);
    result.localWinnerMatchesGlobalAuthority=
        legacyNoWorse(candidateQuality,currentQuality) &&
        better(candidateScore,currentScore);
    if (!result.localWinnerMatchesGlobalAuthority) {
        result.issues.push_back(
            "patch-local winner disagrees with the authoritative global quality result");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const std::set<std::size_t> selectedSet(selected.begin(),selected.end());
    bool outsideStable=true;
    const auto revisionedBase=buildRevisionedTopology2D(topology,incidence);
    for (const auto& cell:topology.cells) if (!selectedSet.contains(cell.id)) {
        const TopologySourceIdentity2D identity{cell.sourceKey,cell.sourceId};
        outsideStable=outsideStable &&
            committed.revisionedTopology.cells.at(identity).vertices==
            revisionedBase.cells.at(identity).vertices;
    }
    result.topology=committed.topology;
    result.immutableCells=std::move(candidateImmutable);
    result.accepted=true;
    result.hardFaceCountAfter=candidateScore.hardCount;
    result.minimumFaceOverLocalHAfter=candidateScore.minimum;
    result.patchOutsideStableIdsUnchanged=outsideStable;
    result.localDeltaMatchesGlobalOracle=committed.deltaMatchesGlobalOracle;
    result.repairSeconds=profileSeconds(repairStart);
    return result;
}

} // namespace cartmesh2d
