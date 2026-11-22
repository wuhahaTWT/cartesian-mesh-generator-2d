#include "SolverTopologyDetail2D.hpp"

#include <cmath>
#include <limits>
#include <set>

namespace cartmesh2d {
namespace solver_topology_detail {

SolverLocalRepartitionResult2D repartitionSolverTopologyByQualityImpl(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool useBatch,
    const std::vector<bool>& initialImmutableCells,bool allowCollinear,
    const SolverQualityPolicy2D& policy) {
    SolverLocalRepartitionResult2D result;
    result.topology=topology;
    result.immutableCells=initialImmutableCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back("local solver repartition requires valid inputs");
        return result;
    }
    for (std::size_t iteration=0;iteration<128;++iteration) {
        const auto quality=timedFullQuality(result.topology,tol,profile,false,policy);
        if (quality.valid()) break;
        if (profile) ++profile->repartitionIterations;
        const auto generationStart=ProfileClock::now();
        const auto pairs=solverRepartitionPairs(
            result.topology,quality,result.immutableCells);
        const auto needsWallCuts=[&](std::size_t first,std::size_t second) {
            return std::any_of(quality.issues.begin(),quality.issues.end(),[&](const auto& issue) {
                return issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness &&
                    (issue.cellId==first || issue.cellId==second);
            });
        };
        if (profile) {
            profile->candidateGenerationSeconds+=profileSeconds(generationStart);
            profile->repartitionCandidatePairs+=pairs.size();
        }
        bool acceptedBatch=false;
        if (useBatch) {
            std::vector<RepartitionProposal2D> proposals;
            for (const auto& [first,second]:pairs) {
                const auto polygonStart=ProfileClock::now();
                const auto merged=mergeAdjacentPolygonsSimple(
                    topologyCellPolygon(result.topology,first),
                    topologyCellPolygon(result.topology,second),tol);
                if (!merged) {
                    if (profile) {
                        profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                    }
                    continue;
                }
                const auto splits=convexTwoPieceSplits(*merged,domain,boundary,tol,
                                                      needsWallCuts(first,second),allowCollinear);
                if (profile) profile->candidateSplits+=splits.size();
                std::optional<RepartitionProposal2D> best;
                const auto halo=cellPairHalo(result.topology,first,second);
                const std::vector<std::size_t> removed{first,second};
                const auto baseRank=localReplacementQualityRank(result.topology,halo,removed,
                    {topologyCellPolygon(result.topology,first),
                     topologyCellPolygon(result.topology,second)},tol,policy);
                const auto consider=[&](std::vector<Polygon2D> pieces) {
                    auto rank=localReplacementQualityRank(
                        result.topology,halo,removed,pieces,tol,policy);
                    if (!rank) return;
                    if (baseRank) rankRelativeToBase(*rank,*baseRank);
                    RepartitionProposal2D proposal;
                    proposal.first=first;
                    proposal.second=second;
                    proposal.pieces=std::move(pieces);
                    proposal.rank=*rank;
                    proposal.halo=halo;
                    if (!best || betterLocalQualityRank(proposal.rank,best->rank))
                        best=std::move(proposal);
                };
                // Exact-union agglomeration already exists in the exhaustive
                // repair path below. Include it in the same independent batch
                // machinery so many small bad interfaces need one global gate,
                // instead of one full rebuild for every union/split candidate.
                const bool convexUnion=strictlyConvex(allowCollinear?
                    removeArtificialCollinearVertices(*merged,tol):*merged);
                bool admissibleUnion=convexUnion;
                if (!convexUnion && allowCollinear) {
                    // The exhaustive alternate below already permits exact
                    // concave unions under these same cell-policy bounds.
                    // Score them with their neighbours before rebuilding the
                    // whole mesh; the authoritative gate is still unchanged.
                    const auto metrics=evaluateSolverCellMetrics2D(*merged,tol);
                    admissibleUnion=metrics.valid &&
                        metrics.maxConcavityDeg<=quality.policy.maxConcavityDeg &&
                        metrics.minInteriorAngleDeg>=quality.policy.minInteriorAngleDeg &&
                        metrics.hydraulicAspect<=quality.policy.maxCellAspect;
                }
                if (admissibleUnion &&
                    !underDeterminedBoundaryCell(*merged,domain,boundary,tol))
                    consider({*merged});
                for (const auto& [firstPiece,secondPiece]:splits)
                    consider({firstPiece,secondPiece});
                if (best) proposals.push_back(std::move(*best));
                if (profile) {
                    profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                }
            }
            // The exhaustive fallback also splits a single skewed wall cell.
            // Rank those same candidates alongside pair repairs, so disjoint
            // wall defects can share one authoritative global rebuild.
            std::set<std::size_t> wallCells;
            for (const auto& issue:quality.issues)
                if (issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness &&
                    issue.cellId<result.topology.cells.size() &&
                    (result.immutableCells.empty() || !result.immutableCells[issue.cellId]))
                    wallCells.insert(issue.cellId);
            for (const auto cell:wallCells) {
                const auto polygonStart=ProfileClock::now();
                const auto polygon=topologyCellPolygon(result.topology,cell);
                const auto splits=convexTwoPieceSplits(polygon,domain,boundary,tol,true,allowCollinear);
                if (profile) profile->candidateSplits+=splits.size();
                const auto halo=cellPairHalo(result.topology,cell,cell);
                const std::vector<std::size_t> removed{cell};
                const auto baseRank=localReplacementQualityRank(
                    result.topology,halo,removed,{polygon},tol,policy);
                std::optional<RepartitionProposal2D> best;
                for (const auto& [firstPiece,secondPiece]:splits) {
                    std::vector<Polygon2D> pieces{firstPiece,secondPiece};
                    auto rank=localReplacementQualityRank(result.topology,halo,removed,pieces,tol,policy);
                    if (!rank) continue;
                    if (baseRank) rankRelativeToBase(*rank,*baseRank);
                    RepartitionProposal2D proposal{cell,cell,std::move(pieces),*rank,halo};
                    if (!best || betterLocalQualityRank(proposal.rank,best->rank)) best=std::move(proposal);
                }
                if (best) proposals.push_back(std::move(*best));
                if (profile) profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
            }
            auto selected=selectIndependentProposals(std::move(proposals));
            if (profile) profile->repairPatchCount+=selected.size();
            for (std::size_t batchSize=selected.size();batchSize>0;) {
                const std::vector<RepartitionProposal2D> batch(
                    selected.begin(),selected.begin()+static_cast<std::ptrdiff_t>(batchSize));
                auto candidate=applyRepartitionBatch(
                    result.topology,batch,domain,boundary,tol,profile,
                    result.immutableCells);
                if (candidate.topology.valid()) {
                    const auto candidateQuality=timedFullQuality(
                        candidate.topology,tol,profile,true,policy);
                    if (betterQualityScore(qualityScore(candidateQuality),qualityScore(quality))) {
                        result.topology=std::move(candidate.topology);
                        result.immutableCells=std::move(candidate.immutableCells);
                        result.repartitionCount+=batchSize;
                        if (profile) {
                            profile->acceptedRepartitions+=batchSize;
                            ++profile->acceptedTopologyCommitCount;
                        }
                        acceptedBatch=true;
                        break;
                    }
                }
                if (batchSize==1) break;
                batchSize=(batchSize+1)/2;
            }
        }
        if (acceptedBatch) continue;

        // Fail closed to the H2 exhaustive global check if the locally ranked
        // independent batch does not strictly improve the authoritative score.
        std::optional<RepartitionBatch2D> bestTopology;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:pairs) {
            // Zero issues has score (0,0,0), the absolute minimum. The first
            // such candidate already wins the exhaustive strict ordering;
            // later candidates cannot replace it. Keep the sequential oracle
            // exhaustive so equivalence remains independently testable.
            if (useBatch && bestScore.issueCount==0U) break;
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                topologyCellPolygon(result.topology,first),
                topologyCellPolygon(result.topology,second),tol);
            if (!merged) {
                if (profile) profile->candidatePolygonWorkSeconds+=
                    profileSeconds(polygonStart);
                continue;
            }
            // A deterministic convex partition may create a tiny triangle
            // beside a regular cell. If their exact union is already convex,
            // retaining two cells cannot improve interpolation weight or
            // volume ratio. Allow a true solver-cell agglomeration, subject to
            // the same full-topology validity and strict quality-score gate.
            if (strictlyConvex(allowCollinear?removeArtificialCollinearVertices(*merged,tol):*merged) &&
                !underDeterminedBoundaryCell(*merged,domain,boundary,tol)) {
                auto candidate=agglomerateCellPair(
                    result.topology,first,second,*merged,domain,boundary,tol,
                    profile,result.immutableCells);
                if (candidate.topology.valid()) {
                    const auto candidateQuality=timedFullQuality(
                        candidate.topology,tol,profile,true,policy);
                    const auto candidateScore=qualityScore(candidateQuality);
                    if (betterQualityScore(candidateScore,bestScore)) {
                        bestScore=candidateScore;
                        bestTopology=std::move(candidate);
                    }
                }
            }
            const auto splits=convexTwoPieceSplits(*merged,domain,boundary,tol,
                                                  needsWallCuts(first,second),allowCollinear);
            if (profile) {
                profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
            }
            for (const auto& [firstPiece,secondPiece]:splits) {
                if (useBatch && bestScore.issueCount==0U) break;
                auto candidate=repartitionPair(result.topology,first,second,
                                               firstPiece,secondPiece,
                                               domain,boundary,tol,profile,
                                               result.immutableCells);
                if (!candidate.topology.valid()) continue;
                const auto candidateQuality=timedFullQuality(
                    candidate.topology,tol,profile,true,policy);
                const auto candidateScore=qualityScore(candidateQuality);
                if (betterQualityScore(candidateScore,bestScore)) {
                    bestScore=candidateScore;
                    bestTopology=std::move(candidate);
                }
            }
        }
        if (!bestTopology) {
            std::set<std::size_t> wallCells;
            for (const auto& issue:quality.issues)
                if (issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness &&
                    issue.cellId<result.topology.cells.size() &&
                    (result.immutableCells.empty() || !result.immutableCells[issue.cellId]))
                    wallCells.insert(issue.cellId);
            for (const auto cell:wallCells) {
                if (useBatch && bestScore.issueCount==0U) break;
                const auto splits=convexTwoPieceSplits(
                    topologyCellPolygon(result.topology,cell),domain,boundary,tol,true,allowCollinear);
                for (const auto& [firstPiece,secondPiece]:splits) {
                    if (useBatch && bestScore.issueCount==0U) break;
                    auto candidate=repartitionPair(result.topology,cell,cell,
                        firstPiece,secondPiece,domain,boundary,tol,profile,result.immutableCells);
                    if (!candidate.topology.valid()) continue;
                    const auto candidateScore=qualityScore(timedFullQuality(
                        candidate.topology,tol,profile,true,policy));
                    if (betterQualityScore(candidateScore,bestScore)) {
                        bestScore=candidateScore;
                        bestTopology=std::move(candidate);
                    }
                }
            }
        }
        if (!bestTopology && allowCollinear) {
            // A curved physical wall can make the exact union slightly concave.
            // Convex-only splitting then leaves a sliver with a poor face weight
            // even though its unsplit polygon satisfies the existing Solver
            // policy. Try that exact union only after the convex repairs fail.
            // No boundary point is removed, no immutable layer cell is touched,
            // and both the cell metrics and authoritative global score gate it.
            for (const auto& [first,second]:pairs) {
                if (useBatch && bestScore.issueCount==0U) break;
                const auto merged=mergeAdjacentPolygonsSimple(
                    topologyCellPolygon(result.topology,first),
                    topologyCellPolygon(result.topology,second),tol);
                if (!merged ||
                    strictlyConvex(removeArtificialCollinearVertices(*merged,tol)) ||
                    underDeterminedBoundaryCell(*merged,domain,boundary,tol)) continue;
                const auto metrics=evaluateSolverCellMetrics2D(*merged,tol);
                if (!metrics.valid ||
                    metrics.maxConcavityDeg>quality.policy.maxConcavityDeg ||
                    metrics.minInteriorAngleDeg<quality.policy.minInteriorAngleDeg ||
                    metrics.hydraulicAspect>quality.policy.maxCellAspect) continue;
                auto candidate=agglomerateCellPair(
                    result.topology,first,second,*merged,domain,boundary,tol,
                    profile,result.immutableCells);
                if (!candidate.topology.valid()) continue;
                const auto candidateScore=qualityScore(timedFullQuality(
                    candidate.topology,tol,profile,true,policy));
                if (betterQualityScore(candidateScore,bestScore)) {
                    bestScore=candidateScore;
                    bestTopology=std::move(candidate);
                }
            }
        }
        if (!bestTopology) break;
        result.topology=std::move(bestTopology->topology);
        result.immutableCells=std::move(bestTopology->immutableCells);
        ++result.repartitionCount;
        if (profile) ++profile->acceptedRepartitions;
        if (profile) ++profile->acceptedTopologyCommitCount;
    }
    if (!allowCollinear && !evaluateSolverQuality2D(result.topology,policy,tol).valid()) {
        auto alternate=repartitionSolverTopologyByQualityImpl(
            topology,domain,boundary,tol,profile,useBatch,initialImmutableCells,true,policy);
        if (alternate.valid() && evaluateSolverQuality2D(alternate.topology,policy,tol).valid())
            return alternate;
    }
    return result;
}

} // namespace solver_topology_detail

using namespace solver_topology_detail;

std::vector<std::size_t> selectIndependentSolverRepairPatches2D(
    const std::vector<std::vector<std::size_t>>& halos) {
    std::vector<std::size_t> selected;
    for (std::size_t candidate=0;candidate<halos.size();++candidate) {
        bool conflict=false;
        for (const auto accepted:selected) {
            if (sortedSetsIntersect(halos[candidate],halos[accepted])) {
                conflict=true;
                break;
            }
        }
        if (!conflict) selected.push_back(candidate);
    }
    return selected;
}

SolverQualityReport2D evaluateSolverQualityPatch2D(
    const TopologyMesh2D& topology,
    const std::vector<std::size_t>& sortedCellIds,
    const TolerancePolicy& tol) {
    if (!topology.valid() || sortedCellIds.empty() ||
        !std::is_sorted(sortedCellIds.begin(),sortedCellIds.end()) ||
        std::adjacent_find(sortedCellIds.begin(),sortedCellIds.end())!=sortedCellIds.end() ||
        sortedCellIds.back()>=topology.cells.size()) {
        return evaluateSolverQuality2D(TopologyMesh2D{}, {}, tol);
    }
    const auto boundary=patchBoundaryRegion(topology,sortedCellIds,tol);
    if (!boundary) return evaluateSolverQuality2D(TopologyMesh2D{}, {}, tol);
    const Domain2D domain{boundary->bounds()};
    std::vector<CutCell2D> cells;
    cells.reserve(sortedCellIds.size());
    for (const auto cell:sortedCellIds) {
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol));
    }
    const auto patch=buildGlobalTopology(cells,domain,*boundary,tol);
    return evaluateSolverQuality2D(patch,{},tol);
}

SolverLocalRepartitionResult2D repartitionSolverTopologyByQuality2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return repartitionSolverTopologyByQualityImpl(
        topology,domain,boundary,tol,nullptr,true);
}

SolverLocalRepartitionResult2D
repartitionSolverTopologyByQualitySequentialReference2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return repartitionSolverTopologyByQualityImpl(
        topology,domain,boundary,tol,nullptr,false);
}

SolverLocalRepartitionResult2D improveSolverForTargetPolicy2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const SolverQualityPolicy2D& policy,const TolerancePolicy& tol) {
    if (!immutableCells.empty() && immutableCells.size()!=topology.cells.size()) {
        SolverLocalRepartitionResult2D invalid;
        invalid.topology=topology;
        invalid.issues.push_back("target solver repair requires one immutable flag per cell");
        return invalid;
    }
    return repartitionSolverTopologyByQualityImpl(
        topology,domain,boundary,tol,nullptr,true,immutableCells,true,policy);
}

double extrudedCellDeterminant2D(const Polygon2D& polygon,double thickness) {
    if (!(thickness>0.0) || !std::isfinite(thickness) || polygon.vertices.size()<3U)
        return std::numeric_limits<double>::quiet_NaN();
    double scale=thickness;
    std::vector<Vector2D> edges;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        const auto edge=polygon.vertices[(i+1U)%polygon.vertices.size()]-polygon.vertices[i];
        const double length=std::hypot(edge.x,edge.y);
        if (!(length>0.0)) return std::numeric_limits<double>::quiet_NaN();
        scale=std::max(scale,length);
        edges.push_back(edge);
    }
    std::vector<double> lengths;
    double perimeter=0.0;
    for (auto& edge:edges) {
        edge=edge*(1.0/scale);
        lengths.push_back(std::hypot(edge.x,edge.y));
        perimeter+=lengths.back();
    }
    double minorSum=0.0;
    for (std::size_t i=0;i<edges.size();++i)
        for (std::size_t j=i+1U;j<edges.size();++j) {
            const double minor=cross(edges[i],edges[j]);
            minorSum+=minor*minor/(lengths[i]*lengths[j]);
        }
    // Stable pairwise-minor form of the in-plane area tensor determinant.
    // Empty front/back faces contribute 2*A to the third diagonal entry.
    const double area=polygon.area()/scale/scale,t=thickness/scale;
    const double areaSum=t*perimeter+2.0*area;
    return 27.0*minorSum*t*t*(2.0*area)/(areaSum*areaSum*areaSum);
}

SolverLocalRepartitionResult2D improveSolverExtrudedDeterminant2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    double thickness,const SolverQualityPolicy2D& policy,const TolerancePolicy& tol) {
    SolverLocalRepartitionResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        !(thickness>0.0) || !std::isfinite(thickness) ||
        (!immutableCells.empty() && immutableCells.size()!=topology.cells.size())) {
        result.issues.push_back("extruded determinant repair requires valid inputs");
        return result;
    }
    const auto failures=[&](const TopologyMesh2D& mesh) {
        std::vector<std::size_t> ids;
        for (const auto& cell:mesh.cells) {
            const double value=extrudedCellDeterminant2D(topologyCellPolygon(mesh,cell.id),thickness);
            if (!std::isfinite(value) || value<0.001) ids.push_back(cell.id);
        }
        return ids;
    };
    auto failed=failures(result.topology);
    for (std::size_t iteration=0;iteration<32U && !failed.empty();++iteration) {
        const auto locked=[&](std::size_t cell) {
            return !result.immutableCells.empty() && result.immutableCells[cell];
        };
        const auto beforeQuality=evaluateSolverQuality2D(result.topology,policy,tol);
        const auto beforeDirectional=evaluateDirectionalConnectivity2D(result.topology);
        std::set<std::pair<std::size_t,std::size_t>> pairs;
        for (const auto cell:failed) {
            if (locked(cell)) continue;
            for (const auto edgeId:result.topology.cells[cell].edges) {
                const auto& edge=result.topology.edges[edgeId];
                if (edge.neighbour && !locked(edge.owner) && !locked(*edge.neighbour))
                    pairs.emplace(std::min(edge.owner,*edge.neighbour),std::max(edge.owner,*edge.neighbour));
            }
        }
        bool accepted=false;
        for (const auto& [first,second]:pairs) {
            const auto merged=mergeAdjacentPolygonsSimple(topologyCellPolygon(result.topology,first),
                topologyCellPolygon(result.topology,second),tol);
            if (!merged || extrudedCellDeterminant2D(*merged,thickness)<0.001) continue;
            const auto metrics=evaluateSolverCellMetrics2D(*merged,tol);
            if (!metrics.valid || metrics.maxConcavityDeg>std::max(1e-5,beforeQuality.maxConcavityDeg)) continue;
            auto candidate=agglomerateCellPair(result.topology,first,second,*merged,
                domain,boundary,tol,nullptr,result.immutableCells);
            if (!candidate.topology.valid()) continue;
            const auto quality=evaluateSolverQuality2D(candidate.topology,policy,tol);
            const auto directional=evaluateDirectionalConnectivity2D(candidate.topology);
            if (!quality.valid() || !directional.issues.empty() ||
                directional.failedCells.size()>beforeDirectional.failedCells.size() ||
                (beforeDirectional.minimumMeasured && directional.minimumMeasured &&
                 *directional.minimumMeasured<std::min(*beforeDirectional.minimumMeasured,minimumDirectionalDeterminant2D))) continue;
            auto candidateFailed=failures(candidate.topology);
            if (candidateFailed.size()>=failed.size()) continue;
            result.topology=std::move(candidate.topology);
            result.immutableCells=std::move(candidate.immutableCells);
            failed=std::move(candidateFailed);
            ++result.repartitionCount;
            accepted=true;
            break;
        }
        if (!accepted) break;
    }
    return result;
}

SolverDirectionalRepairResult2D improveSolverDirectionalConnectivity2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const TolerancePolicy& tol) {
    SolverDirectionalRepairResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    result.before=evaluateDirectionalConnectivity2D(topology);
    result.after=result.before;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        (!immutableCells.empty() && immutableCells.size()!=topology.cells.size()) ||
        !result.before.issues.empty()) {
        result.issues.push_back("directional repair requires valid inputs");
        return result;
    }
    auto currentQuality=evaluateSolverQuality2D(topology,{},tol);
    if (!currentQuality.valid()) {
        result.issues.push_back("directional repair requires the original Solver policy to pass");
        return result;
    }
    const auto noWorse=[](const SolverQualityReport2D& after,const SolverQualityReport2D& before) {
        const auto upper=[](double a,double b) { return a-b<=1e-8*std::max(1.0,std::abs(b)); };
        const auto lower=[](double a,double b) { return b-a<=1e-8*std::max(1.0,std::abs(b)); };
        return after.valid() &&
            upper(after.maxNonOrthogonalityDeg,before.maxNonOrthogonalityDeg) &&
            upper(after.maxInternalSkewness,before.maxInternalSkewness) &&
            upper(after.maxBoundarySkewness,before.maxBoundarySkewness) &&
            upper(after.maxConcavityDeg,before.maxConcavityDeg) &&
            upper(after.maxCellAspect,before.maxCellAspect) &&
            lower(after.minInteriorAngleDeg,before.minInteriorAngleDeg) &&
            lower(after.minFaceLength,before.minFaceLength) &&
            lower(after.minFaceWeight,before.minFaceWeight) &&
            lower(after.minVolumeRatio,before.minVolumeRatio) &&
            lower(after.minCompactness,before.minCompactness);
    };
    for (std::size_t iteration=0;iteration<32U && !result.after.failedCells.empty();++iteration) {
        const auto immutable=[&](std::size_t cell) {
            return !result.immutableCells.empty() && result.immutableCells[cell];
        };
        std::set<std::pair<std::size_t,std::size_t>> pairs;
        for (const auto cell:result.after.failedCells) {
            if (immutable(cell)) continue;
            for (const auto edgeId:result.topology.cells[cell].edges) {
                const auto& edge=result.topology.edges[edgeId];
                if (!edge.neighbour || immutable(edge.owner) || immutable(*edge.neighbour)) continue;
                pairs.emplace(std::min(edge.owner,*edge.neighbour),std::max(edge.owner,*edge.neighbour));
            }
        }
        bool accepted=false;
        for (const auto& [first,second]:pairs) {
            ++result.candidateCount;
            const auto inPair=[&](std::size_t cell) { return cell==first || cell==second; };
            // All outer atomic faces are retained by an exact union. Removing
            // their common internal face is the only change to this tensor.
            std::set<std::size_t> outerInternalEdges;
            for (const auto cell:{first,second}) {
                for (const auto edgeId:result.topology.cells[cell].edges) {
                    const auto& edge=result.topology.edges[edgeId];
                    if (edge.neighbour && inPair(edge.owner)!=inPair(*edge.neighbour))
                        outerInternalEdges.insert(edgeId);
                }
            }
            std::vector<Vector2D> directions;
            for (const auto id:outerInternalEdges) {
                const auto& edge=result.topology.edges[id];
                directions.push_back(result.topology.vertices[edge.v1].point-
                                     result.topology.vertices[edge.v0].point);
            }
            const double mergedDeterminant=directionalDeterminant2D(directions);
            if (!std::isfinite(mergedDeterminant) ||
                mergedDeterminant<minimumDirectionalDeterminant2D) continue;
            const auto merged=mergeAdjacentPolygonsSimple(
                topologyCellPolygon(result.topology,first),topologyCellPolygon(result.topology,second),tol);
            if (!merged) continue;
            const auto metrics=evaluateSolverCellMetrics2D(*merged,tol);
            if (!metrics.valid || metrics.maxConcavityDeg>currentQuality.policy.maxConcavityDeg ||
                metrics.minInteriorAngleDeg<currentQuality.policy.minInteriorAngleDeg ||
                metrics.hydraulicAspect>currentQuality.policy.maxCellAspect) continue;
            auto candidate=agglomerateCellPair(result.topology,first,second,*merged,
                domain,boundary,tol,nullptr,result.immutableCells);
            if (!candidate.topology.valid()) continue;
            auto candidateQuality=evaluateSolverQuality2D(candidate.topology,{},tol);
            if (!noWorse(candidateQuality,currentQuality)) continue;
            auto connectivity=evaluateDirectionalConnectivity2D(candidate.topology);
            if (!connectivity.issues.empty() || !connectivity.minimumMeasured ||
                connectivity.failedCells.size()>=result.after.failedCells.size() ||
                *connectivity.minimumMeasured<*result.after.minimumMeasured) continue;
            result.topology=std::move(candidate.topology);
            result.immutableCells=std::move(candidate.immutableCells);
            result.after=std::move(connectivity);
            currentQuality=std::move(candidateQuality);
            ++result.acceptedCount;
            accepted=true;
            break;
        }
        if (!accepted) break;
    }
    return result;
}

} // namespace cartmesh2d
