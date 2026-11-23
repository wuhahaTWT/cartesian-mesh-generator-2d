#include "SolverTopologyDetail2D.hpp"

#include <cmath>

namespace cartmesh2d {
using namespace solver_topology_detail;

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return buildSolverTopology2D(
        topology,domain,boundary,SolverTopologyConstraints2D{},tol);
}

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const SolverTopologyConstraints2D& constraints,
    const TolerancePolicy& tol) {
    SolverTopologyResult2D result;
    result.inputCellCount=topology.cells.size();
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        (!constraints.immutableInputCells.empty() &&
         constraints.immutableInputCells.size()!=topology.cells.size()) ||
        (!constraints.preserveInputCells.empty() &&
         constraints.preserveInputCells.size()!=topology.cells.size()) ||
        (!constraints.inputPolygonOverrides.empty() &&
         constraints.inputPolygonOverrides.size()!=topology.cells.size())) {
        result.issues.push_back("solver topology repair requires valid inputs");
        return result;
    }
    std::vector<Polygon2D> sourcePolygons;
    std::vector<std::vector<std::size_t>> sourceLineages;
    sourcePolygons.reserve(topology.cells.size());
    sourceLineages.reserve(topology.cells.size());
    for (std::size_t cellId=0;cellId<topology.cells.size();++cellId) {
        const auto& cell=topology.cells[cellId];
        sourceLineages.push_back(cell.sourceLineage.empty()
            ?std::vector<std::size_t>{cell.sourceId}:cell.sourceLineage);
        if (!constraints.inputPolygonOverrides.empty() &&
            constraints.inputPolygonOverrides[cellId].has_value()) {
            const auto& polygon=*constraints.inputPolygonOverrides[cellId];
            const BoundaryLoop loop(polygon.vertices);
            const double areaScale=std::max(std::abs(polygon.area()),
                                            std::abs(cell.geometryArea));
            if (!loop.diagnose(tol).valid() || !(polygon.signedArea()>0.0) ||
                std::abs(polygon.area()-cell.geometryArea)>
                    tol.absolute*tol.absolute+tol.relative*areaScale) {
                result.issues.push_back(
                    "solver topology polygon override does not match input cell");
                return result;
            }
            sourcePolygons.push_back(polygon);
            continue;
        }
        Polygon2D polygon;
        polygon.vertices.reserve(cell.vertices.size());
        for (const auto vertex:cell.vertices) {
            if (vertex>=topology.vertices.size()) {
                result.issues.push_back("solver topology cell has invalid vertex index");
                return result;
            }
            polygon.vertices.push_back(topology.vertices[vertex].point);
        }
        sourcePolygons.push_back(std::move(polygon));
    }

    const auto initialPartitionStart=ProfileClock::now();
    PartitionAttempt2D partition=partitionSourcePolygons(
        sourcePolygons,domain,boundary,tol,&result.profile,false,
        constraints.immutableInputCells,constraints.preserveInputCells,
        sourceLineages,topology.constructionRegistry,
        topology.embeddedBoundaryFragments);
    std::vector<bool> immutableSources=constraints.immutableInputCells;
    std::vector<bool> preserveSources=constraints.preserveInputCells;
    result.profile.initialPartitionSeconds=profileSeconds(initialPartitionStart);
    if (!partition.error.empty()) {
        result.issues.push_back(partition.error);
        return result;
    }
    const auto sourceRepairStart=ProfileClock::now();
    for (std::size_t iteration=0;iteration<32;++iteration) {
        const auto quality=timedFullQuality(partition.topology,tol,&result.profile,false);
        if (quality.valid()) break;
        ++result.profile.sourceRepairIterations;
        const auto generationStart=ProfileClock::now();
        const auto pairs=sourceMergePairs(partition,quality,immutableSources);
        result.profile.candidateGenerationSeconds+=profileSeconds(generationStart);
        result.profile.sourceCandidatePairs+=pairs.size();
        const auto adjacency=sourceAdjacency(partition,sourcePolygons.size());
        std::vector<std::vector<std::size_t>> cellsForSource(sourcePolygons.size());
        for (std::size_t cell=0;cell<partition.sourceForCell.size();++cell)
            cellsForSource[partition.sourceForCell[cell]].push_back(cell);
        std::vector<SourceMergeProposal2D> proposals;
        std::vector<std::pair<std::size_t,std::size_t>> fallbackPairs;
        const auto locallyImproves=[](const LocalQualityRank2D& rank) noexcept {
            return rank.issueDelta<0 ||
                   (rank.issueDelta==0 &&
                    (rank.maximumSeverityDelta<-1.0e-12 ||
                     rank.totalSeverityDelta<-1.0e-12));
        };
        for (const auto& [first,second]:pairs) {
            if (second>=sourcePolygons.size()) continue;
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                sourcePolygons[first],sourcePolygons[second],tol);
            if (merged) {
                const auto pieces=partitionOnePolygon(*merged,domain,boundary,tol);
                if (pieces) {
                    const auto sourceHalo=sourcePairHalo(first,second,adjacency);
                    std::vector<std::size_t> cellHalo;
                    for (const auto source:sourceHalo)
                        cellHalo.insert(cellHalo.end(),cellsForSource[source].begin(),
                                        cellsForSource[source].end());
                    std::sort(cellHalo.begin(),cellHalo.end());
                    auto removed=cellsForSource[first];
                    removed.insert(removed.end(),cellsForSource[second].begin(),
                                   cellsForSource[second].end());
                    std::sort(removed.begin(),removed.end());
                    // A merge changes centres on interfaces with its neighbours.
                    // Rank those real interfaces as well as the merged interior;
                    // isolated-polygon scoring misses their face weight and skew.
                    // Because every face and cell outside cellHalo is unchanged,
                    // a pair that does not improve within cellHalo cannot improve
                    // the global QualityScore2D; pre-screening both batch and
                    // fallback candidates with one LocalPatchQualityScope avoids
                    // rebuilding the full topology on provably non-improving pairs.
                    bool screenedOut=false;
                    std::optional<LocalQualityRank2D> rank;
                    const auto scope=makeLocalPatchQualityScope(
                        partition.topology,cellHalo,removed,tol);
                    if (scope.valid()) {
                        rank=scope.evaluate(*pieces);
                        if (rank) {
                            std::vector<Polygon2D> originalPieces;
                            originalPieces.reserve(removed.size());
                            for (const auto cell:removed)
                                originalPieces.push_back(topologyCellPolygon(partition.topology,cell));
                            const auto baseRank=scope.evaluate(originalPieces);
                            if (baseRank) {
                                rankRelativeToBase(*rank,*baseRank);
                                if (!locallyImproves(*rank)) screenedOut=true;
                            }
                        }
                    }
                    if (!screenedOut) {
                        fallbackPairs.emplace_back(first,second);
                        if (!rank) rank=localPartitionQualityRank(*merged,*pieces,tol);
                        if (rank) {
                            SourceMergeProposal2D proposal;
                            proposal.first=first;
                            proposal.second=second;
                            proposal.merged=*merged;
                            proposal.rank=*rank;
                            proposal.halo=sourceHalo;
                            proposals.push_back(std::move(proposal));
                        }
                    }
                }
            }
            result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
        }
        auto selected=selectIndependentProposals(std::move(proposals));
        result.profile.repairPatchCount+=selected.size();
        bool acceptedBatch=false;
        for (std::size_t batchSize=selected.size();batchSize>0;) {
            const std::vector<SourceMergeProposal2D> batch(
                selected.begin(),selected.begin()+static_cast<std::ptrdiff_t>(batchSize));
            auto candidatePolygons=applySourceMergeBatch(sourcePolygons,batch);
            auto candidateImmutable=applySourceMergeProtection(immutableSources,batch);
            auto candidatePreserve=applySourceMergePreservation(preserveSources,batch);
            auto candidateLineages=applySourceMergeLineage(sourceLineages,batch);
            auto candidate=partitionSourcePolygons(
                candidatePolygons,domain,boundary,tol,&result.profile,true,
                candidateImmutable,candidatePreserve,candidateLineages,
                topology.constructionRegistry,
                topology.embeddedBoundaryFragments);
            if (candidate.error.empty()) {
                const auto candidateQuality=timedFullQuality(
                    candidate.topology,tol,&result.profile,true);
                if (betterQualityScore(qualityScore(candidateQuality),qualityScore(quality))) {
                    partition=std::move(candidate);
                    sourcePolygons=std::move(candidatePolygons);
                    immutableSources=std::move(candidateImmutable);
                    preserveSources=std::move(candidatePreserve);
                    sourceLineages=std::move(candidateLineages);
                    result.qualityAgglomeratedSourceCellCount+=batchSize;
                    result.profile.acceptedSourceRepairs+=batchSize;
                    ++result.profile.acceptedTopologyCommitCount;
                    acceptedBatch=true;
                    break;
                }
            }
            if (batchSize==1) break;
            batchSize=(batchSize+1)/2;
        }
        if (acceptedBatch) continue;

        // Preserve the H2 exhaustive path as a correctness fallback when the
        // independent batch is invalid or does not improve the global score.
        std::optional<PartitionAttempt2D> bestPartition;
        std::optional<std::vector<Polygon2D>> bestPolygons;
        std::optional<std::vector<bool>> bestImmutable;
        std::optional<std::vector<bool>> bestPreserve;
        std::optional<std::vector<std::vector<std::size_t>>> bestLineages;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:fallbackPairs) {
            // No candidate can improve on a zero-issue authoritative score.
            if (bestScore.issueCount==0U) break;
            if (second>=sourcePolygons.size()) continue;
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                sourcePolygons[first],sourcePolygons[second],tol);
            if (!merged) {
                result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                continue;
            }
            auto candidatePolygons=sourcePolygons;
            auto candidateImmutable=immutableSources;
            auto candidatePreserve=preserveSources;
            auto candidateLineages=sourceLineages;
            candidatePolygons[first]=*merged;
            candidatePolygons.erase(
                candidatePolygons.begin()+static_cast<std::ptrdiff_t>(second));
            if (!candidateImmutable.empty()) {
                candidateImmutable.erase(
                    candidateImmutable.begin()+static_cast<std::ptrdiff_t>(second));
            }
            if (!candidatePreserve.empty()) {
                candidatePreserve[first]=false;
                candidatePreserve.erase(
                    candidatePreserve.begin()+static_cast<std::ptrdiff_t>(second));
            }
            candidateLineages[first]=mergedLineage(
                candidateLineages[first],candidateLineages[second]);
            candidateLineages.erase(
                candidateLineages.begin()+static_cast<std::ptrdiff_t>(second));
            result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
            auto candidate=partitionSourcePolygons(
                candidatePolygons,domain,boundary,tol,&result.profile,true,
                candidateImmutable,candidatePreserve,candidateLineages,
                topology.constructionRegistry,
                topology.embeddedBoundaryFragments);
            if (!candidate.error.empty()) continue;
            const auto candidateQuality=timedFullQuality(
                candidate.topology,tol,&result.profile,true);
            const auto candidateScore=qualityScore(candidateQuality);
            if (betterQualityScore(candidateScore,bestScore)) {
                bestScore=candidateScore;
                bestPartition=std::move(candidate);
                bestPolygons=std::move(candidatePolygons);
                bestImmutable=std::move(candidateImmutable);
                bestPreserve=std::move(candidatePreserve);
                bestLineages=std::move(candidateLineages);
            }
        }
        if (!bestPartition) break;
        partition=std::move(*bestPartition);
        sourcePolygons=std::move(*bestPolygons);
        immutableSources=std::move(*bestImmutable);
        preserveSources=std::move(*bestPreserve);
        sourceLineages=std::move(*bestLineages);
        ++result.qualityAgglomeratedSourceCellCount;
        ++result.profile.acceptedSourceRepairs;
        ++result.profile.acceptedTopologyCommitCount;
    }
    result.profile.sourceRepairSeconds=profileSeconds(sourceRepairStart);
    const auto repartitionStart=ProfileClock::now();
    auto repartitioned=repartitionSolverTopologyByQualityImpl(
        partition.topology,domain,boundary,tol,&result.profile,true,
        partition.immutableForCell);
    result.profile.finalRepartitionSeconds=profileSeconds(repartitionStart);
    if (!repartitioned.valid()) {
        result.issues.insert(result.issues.end(),repartitioned.issues.begin(),
                             repartitioned.issues.end());
        return result;
    }
    partition.topology=std::move(repartitioned.topology);
    result.qualityRepartitionCount=repartitioned.repartitionCount;
    result.topology=std::move(partition.topology);
    result.immutableOutputCells=std::move(repartitioned.immutableCells);
    result.partitionedCellCount=partition.partitionedSourceCellCount;
    result.outputCellCount=result.topology.cells.size();
    return result;
}

} // namespace cartmesh2d
