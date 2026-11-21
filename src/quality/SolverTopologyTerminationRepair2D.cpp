#include "SolverTopologyDetail2D.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace cartmesh2d {
using namespace solver_topology_detail;

SolverTerminationQualityRepairResult2D repairSolverTerminationQuality2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const std::vector<bool>& terminationCells,
    const std::vector<bool>& cartesianCells,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells,double minimumFaceOverLocalH,
    double minimumFaceWeight,double minimumVolumeRatio,
    TerminationQualityCandidateMode2D candidateMode,
    const TolerancePolicy& tol) {
    const auto repairStart=ProfileClock::now();
    SolverTerminationQualityRepairResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    result.terminationCells=terminationCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        immutableCells.size()!=topology.cells.size() ||
        terminationCells.size()!=topology.cells.size() ||
        cartesianCells.size()!=topology.cells.size() ||
        localBackgroundH.size()!=topology.cells.size() ||
        ratedCells.size()!=topology.cells.size() || !(minimumFaceOverLocalH>0.0) ||
        !(minimumFaceWeight>0.0) || !(minimumVolumeRatio>0.0)) {
        result.issues.push_back(
            "Q3 termination repair requires valid aligned inputs and hard limits");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    struct TargetScore {
        std::size_t volumeCount=0U;
        std::size_t faceWeightCount=0U;
        double maximumVolumeSeverity=0.0;
        double totalVolumeSeverity=0.0;
        double maximumFaceWeightSeverity=0.0;
        double totalFaceWeightSeverity=0.0;
        std::size_t minimumAngleCount=0U;
        std::size_t nonOrthogonalityCount=0U;
        std::size_t skewnessCount=0U;
        std::size_t aspectCount=0U;
    };
    struct ShortScore {
        std::size_t count=0U;
        double maximumSeverity=0.0;
        double totalSeverity=0.0;
    };
    const auto cellMetrics=[&](const TopologyMesh2D& mesh) {
        std::vector<SolverCellMetrics2D> metrics(mesh.cells.size());
        for (const auto& cell:mesh.cells)
            metrics[cell.id]=evaluateSolverCellMetrics2D(
                topologyCellPolygon(mesh,cell.id),tol);
        return metrics;
    };
    const auto scoreTarget=[&](const TopologyMesh2D& mesh,
                               const std::vector<bool>& rated) {
        TargetScore score;
        const auto metrics=cellMetrics(mesh);
        for (std::size_t cell=0;cell<metrics.size();++cell) {
            if (!rated[cell]) continue;
            if (metrics[cell].minInteriorAngleDeg<10.0) ++score.minimumAngleCount;
            if (metrics[cell].hydraulicAspect>50.0) ++score.aspectCount;
        }
        for (const auto& edge:mesh.edges) {
            if (!edge.neighbour ||
                (!rated[edge.owner] && !rated[*edge.neighbour])) continue;
            const auto face=evaluateSolverInternalFaceMetrics2D(
                mesh.vertices[edge.v0].point,mesh.vertices[edge.v1].point,
                metrics[edge.owner].centroid,metrics[*edge.neighbour].centroid,
                metrics[edge.owner].area,metrics[*edge.neighbour].area,tol);
            if (!face.orientationValid) continue;
            if (face.nonOrthogonalityDeg>65.0) ++score.nonOrthogonalityCount;
            if (face.skewness>3.0) ++score.skewnessCount;
            if (face.volumeRatio<minimumVolumeRatio) {
                ++score.volumeCount;
                const double severity=minimumVolumeRatio/
                    (face.volumeRatio+std::numeric_limits<double>::min());
                score.maximumVolumeSeverity=
                    std::max(score.maximumVolumeSeverity,severity);
                score.totalVolumeSeverity+=severity;
            }
            if (face.faceWeight<minimumFaceWeight) {
                ++score.faceWeightCount;
                const double severity=minimumFaceWeight/
                    (face.faceWeight+std::numeric_limits<double>::min());
                score.maximumFaceWeightSeverity=
                    std::max(score.maximumFaceWeightSeverity,severity);
                score.totalFaceWeightSeverity+=severity;
            }
        }
        return score;
    };
    const auto scoreShort=[&](const TopologyMesh2D& mesh,
                              const std::vector<double>& localH,
                              const std::vector<bool>& rated) {
        ShortScore score;
        for (const auto& edge:mesh.edges) {
            double h=0.0;
            if (rated[edge.owner]) h=std::max(h,localH[edge.owner]);
            if (edge.neighbour && rated[*edge.neighbour])
                h=std::max(h,localH[*edge.neighbour]);
            if (!(h>0.0)) continue;
            const double ratio=std::sqrt(squaredNorm(
                mesh.vertices[edge.v1].point-mesh.vertices[edge.v0].point))/h;
            if (ratio<minimumFaceOverLocalH) {
                ++score.count;
                const double severity=minimumFaceOverLocalH/
                    (ratio+std::numeric_limits<double>::min());
                score.maximumSeverity=std::max(score.maximumSeverity,severity);
                score.totalSeverity+=severity;
            }
        }
        return score;
    };
    const auto targetBetter=[](const TargetScore& candidate,const TargetScore& base) {
        const std::size_t candidateCount=
            candidate.volumeCount+candidate.faceWeightCount;
        const std::size_t baseCount=base.volumeCount+base.faceWeightCount;
        return candidate.volumeCount<=base.volumeCount &&
               candidate.faceWeightCount<=base.faceWeightCount &&
               candidate.minimumAngleCount<=base.minimumAngleCount &&
               candidate.nonOrthogonalityCount<=base.nonOrthogonalityCount &&
               candidate.skewnessCount<=base.skewnessCount &&
               candidate.aspectCount<=base.aspectCount &&
               candidateCount<baseCount;
    };
    const auto shortNoWorse=[](const ShortScore& candidate,const ShortScore& base) {
        return candidate.count<=base.count &&
               candidate.maximumSeverity<=base.maximumSeverity &&
               candidate.totalSeverity<=base.totalSeverity;
    };
    const auto solverNoWorse=[](const SolverQualityReport2D& candidate,
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
    const auto constructionContextNoWorse=[&](
        const SolverQualityReport2D& candidate,
        const SolverQualityReport2D& current) {
        const auto upper=[](double after,double before) {
            return after-before<=1.0e-8*std::max(1.0,std::abs(before));
        };
        const auto lower=[](double after,double before) {
            return before-after<=1.0e-8*std::max(1.0,std::abs(before));
        };
        return candidate.issues.size()<=current.issues.size() &&
               upper(candidate.maxNonOrthogonalityDeg,
                     current.maxNonOrthogonalityDeg) &&
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
    const bool constructionSelection=
        candidateMode==
            TerminationQualityCandidateMode2D::ConstructionAgglomeration;
    if (!constructionSelection && !currentQuality.valid()) {
        result.issues.push_back(
            "Q3 termination repair requires a solver-safe input topology");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const auto currentTarget=scoreTarget(topology,ratedCells);
    const auto currentShort=scoreShort(topology,localBackgroundH,ratedCells);
    result.hardVolumeRatioCountBefore=currentTarget.volumeCount;
    result.hardVolumeRatioCountAfter=currentTarget.volumeCount;
    result.hardFaceWeightCountBefore=currentTarget.faceWeightCount;
    result.hardFaceWeightCountAfter=currentTarget.faceWeightCount;
    result.maximumVolumeRatioSeverityBefore=currentTarget.maximumVolumeSeverity;
    result.maximumVolumeRatioSeverityAfter=currentTarget.maximumVolumeSeverity;
    result.maximumFaceWeightSeverityBefore=currentTarget.maximumFaceWeightSeverity;
    result.maximumFaceWeightSeverityAfter=currentTarget.maximumFaceWeightSeverity;
    result.hardShortFaceCountBefore=currentShort.count;
    result.hardShortFaceCountAfter=currentShort.count;

    struct TargetEdge {
        std::size_t hardCount=0U;
        double volumeSeverity=0.0;
        double faceWeightSeverity=0.0;
        std::size_t first=0U;
        std::size_t second=0U;
    };
    const auto metrics=cellMetrics(topology);
    std::vector<TargetEdge> targetEdges;
    for (const auto& edge:topology.edges) {
        if (!edge.neighbour ||
            (!terminationCells[edge.owner] && !terminationCells[*edge.neighbour]) ||
            (!ratedCells[edge.owner] && !ratedCells[*edge.neighbour])) continue;
        if (candidateMode!=TerminationQualityCandidateMode2D::Agglomeration &&
            candidateMode!=
                TerminationQualityCandidateMode2D::ConstructionAgglomeration) {
            const bool terminationCartesian=
                (terminationCells[edge.owner] && cartesianCells[*edge.neighbour]) ||
                (cartesianCells[edge.owner] && terminationCells[*edge.neighbour]);
            if (!terminationCartesian) continue;
        }
        const auto face=evaluateSolverInternalFaceMetrics2D(
            topology.vertices[edge.v0].point,topology.vertices[edge.v1].point,
            metrics[edge.owner].centroid,metrics[*edge.neighbour].centroid,
            metrics[edge.owner].area,metrics[*edge.neighbour].area,tol);
        if (!face.orientationValid) continue;
        TargetEdge target;
        target.first=std::min(edge.owner,*edge.neighbour);
        target.second=std::max(edge.owner,*edge.neighbour);
        if (face.volumeRatio<minimumVolumeRatio) {
            ++target.hardCount;
            target.volumeSeverity=minimumVolumeRatio/
                (face.volumeRatio+std::numeric_limits<double>::min());
        }
        if (face.faceWeight<minimumFaceWeight) {
            ++target.hardCount;
            target.faceWeightSeverity=minimumFaceWeight/
                (face.faceWeight+std::numeric_limits<double>::min());
        }
        if (target.hardCount>0U) targetEdges.push_back(target);
    }
    std::sort(targetEdges.begin(),targetEdges.end(),[](const auto& lhs,const auto& rhs) {
        return std::tuple{-static_cast<std::ptrdiff_t>(lhs.hardCount),
                          -lhs.volumeSeverity,-lhs.faceWeightSeverity,
                          lhs.first,lhs.second}<
               std::tuple{-static_cast<std::ptrdiff_t>(rhs.hardCount),
                          -rhs.volumeSeverity,-rhs.faceWeightSeverity,
                          rhs.first,rhs.second};
    });
    if (targetEdges.empty()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    result.applicable=true;
    std::vector<std::vector<std::size_t>> patches;
    if (candidateMode==TerminationQualityCandidateMode2D::Agglomeration ||
        candidateMode==
            TerminationQualityCandidateMode2D::ConstructionAgglomeration) {
        const std::array<std::size_t,2> affected{
            targetEdges.front().first,targetEdges.front().second};
        if (!immutableCells[affected[0]] && !immutableCells[affected[1]])
            patches.push_back({affected[0],affected[1]});
        for (const auto cell:affected) {
            if (immutableCells[cell]) continue;
            for (const auto edgeId:topology.cells[cell].edges) {
                const auto& edge=topology.edges[edgeId];
                if (!edge.neighbour) continue;
                const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
                if (!immutableCells[other])
                    patches.push_back({std::min(cell,other),std::max(cell,other)});
            }
        }
    } else if (candidateMode==TerminationQualityCandidateMode2D::Repartition) {
        // Q3-2 scans only a small deterministic prefix of the remaining
        // termination-Cartesian hard faces. Each patch is the face's direct
        // two-cell union; no one-ring pair expansion or global search occurs.
        constexpr std::size_t patchBound=4U;
        for (const auto& target:targetEdges) {
            if (patches.size()==patchBound) break;
            if (immutableCells[target.first] || immutableCells[target.second]) continue;
            patches.push_back({target.first,target.second});
        }
    } else {
        struct GroupTarget {
            std::size_t hardCount=0U;
            double maximumVolumeSeverity=0.0;
            double maximumFaceWeightSeverity=0.0;
            std::vector<std::size_t> cells;
        };
        std::map<std::size_t,std::vector<TargetEdge>> byCartesian;
        for (const auto& target:targetEdges) {
            const std::size_t cartesian=cartesianCells[target.first]
                ?target.first:target.second;
            byCartesian[cartesian].push_back(target);
        }
        std::vector<GroupTarget> groups;
        for (auto& [cartesian,edges]:byCartesian) {
            std::sort(edges.begin(),edges.end(),[](const auto& lhs,const auto& rhs) {
                return std::tuple{-static_cast<std::ptrdiff_t>(lhs.hardCount),
                                  -lhs.volumeSeverity,-lhs.faceWeightSeverity,
                                  lhs.first,lhs.second}<
                       std::tuple{-static_cast<std::ptrdiff_t>(rhs.hardCount),
                                  -rhs.volumeSeverity,-rhs.faceWeightSeverity,
                                  rhs.first,rhs.second};
            });
            std::vector<std::size_t> hardTerminationNeighbours;
            for (const auto& edge:edges) {
                const std::size_t termination=terminationCells[edge.first]
                    ?edge.first:edge.second;
                if (std::find(hardTerminationNeighbours.begin(),
                              hardTerminationNeighbours.end(),termination)==
                    hardTerminationNeighbours.end())
                    hardTerminationNeighbours.push_back(termination);
            }
            std::vector<std::size_t> companionNeighbours;
            for (const auto edgeId:topology.cells[cartesian].edges) {
                const auto& edge=topology.edges[edgeId];
                if (!edge.neighbour) continue;
                const std::size_t other=edge.owner==cartesian
                    ?*edge.neighbour:edge.owner;
                if (!terminationCells[other] ||
                    std::find(hardTerminationNeighbours.begin(),
                              hardTerminationNeighbours.end(),other)!=
                        hardTerminationNeighbours.end()) continue;
                companionNeighbours.push_back(other);
            }
            std::sort(companionNeighbours.begin(),companionNeighbours.end(),
                      [&](std::size_t lhs,std::size_t rhs) {
                return std::tie(metrics[lhs].area,lhs)<
                       std::tie(metrics[rhs].area,rhs);
            });
            companionNeighbours.erase(std::unique(companionNeighbours.begin(),
                                                   companionNeighbours.end()),
                                      companionNeighbours.end());
            if (hardTerminationNeighbours.size()<2U &&
                companionNeighbours.empty()) continue;
            // One fixed template: the Cartesian cell plus the two highest
            // priority termination neighbours. The second may be a non-hard
            // thin companion on the same Cartesian fan; larger fans are not
            // searched.
            const std::size_t firstTermination=hardTerminationNeighbours[0];
            const std::size_t secondTermination=hardTerminationNeighbours.size()>1U
                ?hardTerminationNeighbours[1]:companionNeighbours[0];
            GroupTarget group;
            group.cells={cartesian,firstTermination,secondTermination};
            for (const auto& edge:edges) {
                const std::size_t termination=terminationCells[edge.first]
                    ?edge.first:edge.second;
                if (termination!=firstTermination &&
                    termination!=secondTermination) continue;
                group.hardCount+=edge.hardCount;
                group.maximumVolumeSeverity=std::max(
                    group.maximumVolumeSeverity,edge.volumeSeverity);
                group.maximumFaceWeightSeverity=std::max(
                    group.maximumFaceWeightSeverity,edge.faceWeightSeverity);
            }
            groups.push_back(std::move(group));
        }
        std::sort(groups.begin(),groups.end(),[](const auto& lhs,const auto& rhs) {
            return std::tuple{-static_cast<std::ptrdiff_t>(lhs.hardCount),
                              -lhs.maximumVolumeSeverity,
                              -lhs.maximumFaceWeightSeverity,lhs.cells}<
                   std::tuple{-static_cast<std::ptrdiff_t>(rhs.hardCount),
                              -rhs.maximumVolumeSeverity,
                              -rhs.maximumFaceWeightSeverity,rhs.cells};
        });
        constexpr std::size_t groupedPatchBound=4U;
        for (const auto& group:groups) {
            if (patches.size()==groupedPatchBound) break;
            if (std::any_of(group.cells.begin(),group.cells.end(),
                            [&](std::size_t cell) { return immutableCells[cell]; }))
                continue;
            patches.push_back(group.cells);
        }
    }
    std::sort(patches.begin(),patches.end());
    patches.erase(std::unique(patches.begin(),patches.end()),patches.end());
    constexpr std::size_t candidateBound=16U;
    if (patches.size()>candidateBound) patches.resize(candidateBound);
    if (patches.empty()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }

    const auto incidence=buildEdgeIncidenceStore2D(topology,0U);
    if (!incidence.valid()) {
        result.issues.push_back("Q3 termination base incidence is invalid");
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

    struct Candidate {
        std::size_t first=0U;
        std::size_t second=0U;
        std::vector<std::size_t> selected;
        TopologyPatchTransaction2D transaction;
        std::vector<TopologyReplacementCell2D> replacements;
        PatchLocalTerminationRank2D rank;
    };
    std::optional<Candidate> winner;
    const auto selectionStartGlobalBuilds=globalTopologyBuildCount2D();
    const auto selectionStartFullQuality=solverQualityEvaluationCount2D();
    const auto lockPolygon=[&](const Polygon2D& polygon,
                               const TopologyPatchTransaction2D& transaction) {
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
                const double denominator=squaredNorm(direction);
                const double position=denominator>0.0
                    ?dot(found->second.point-a,direction)/denominator:0.0;
                points.emplace_back(position,found->second);
            };
            add(identityByPoint.at({a.x,a.y}).stableId);
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
            for (const auto& [position,vertex]:points) {
                (void)position;
                if (locked.vertices.empty() ||
                    locked.vertices.back().x!=vertex.point.x ||
                    locked.vertices.back().y!=vertex.point.y)
                    locked.vertices.push_back(vertex.point);
            }
        }
        if (locked.vertices.size()>1U &&
            locked.vertices.front().x==locked.vertices.back().x &&
            locked.vertices.front().y==locked.vertices.back().y)
            locked.vertices.pop_back();
        return locked;
    };
    for (const auto& selected:patches) {
        const std::size_t first=selected[0];
        const std::size_t second=selected[1];
        const bool agglomerationMode=
            candidateMode==TerminationQualityCandidateMode2D::Agglomeration ||
            candidateMode==
                TerminationQualityCandidateMode2D::ConstructionAgglomeration;
        if (agglomerationMode)
            ++result.candidateCount;
        std::optional<Polygon2D> merged=topologyCellPolygon(topology,selected[0]);
        for (std::size_t i=1U;i<selected.size() && merged;++i)
            merged=mergeAdjacentPolygonsSimple(
                *merged,topologyCellPolygon(topology,selected[i]),tol);
        if (!merged) continue;
        const auto repairUnion=removeArtificialCollinearVertices(*merged,tol);
        const auto transaction=prepareTopologyPatchTransaction2D(
            topology,incidence,selected);
        if (!transaction.valid()) continue;
        const auto baseScope=buildPatchLocalScope2D(
            topology,incidence,selected,localBackgroundH,ratedCells);
        if (!baseScope.valid()) continue;
        std::set<StableEdgeKey2D> physicalLocks;
        for (const auto& lock:transaction.boundaryLocks)
            if (lock.physicalBoundary) physicalLocks.insert(lock.stableEdge);
        const PatchLocalHardLimits2D hardLimits{
            minimumFaceWeight,minimumVolumeRatio};
        double patchLocalH=0.0;
        bool patchRated=false;
        for (const auto cell:selected) {
            patchLocalH=std::max(patchLocalH,localBackgroundH[cell]);
            patchRated=patchRated || ratedCells[cell];
        }
        std::optional<PatchLocalQuality2D> baseLocal;
        const auto evaluateReplacements=[&](
            std::vector<TopologyReplacementCell2D> replacements,
            bool countCandidate,bool eligibleForWinner) {
            if (countCandidate) ++result.candidateCount;
            const auto localDelta=buildTopologyDelta2D(
                topology,incidence,transaction,replacements,tol);
            if (!localDelta.valid()) return false;
            ++result.localCandidateCount;
            if (!baseLocal) {
                baseLocal=evaluatePatchLocalQuality2D(
                    baseScope.cells,minimumFaceOverLocalH,{},tol,hardLimits);
                ++result.localQualityEvaluationCount;
            }
            if (!baseLocal->valid()) return false;
            std::vector<PatchLocalCell2D> candidateScope;
            for (const auto& entry:baseScope.cells)
                if (!entry.inPatch) candidateScope.push_back(entry);
            for (const auto& replacement:replacements) {
                PatchLocalCell2D entry;
                entry.inPatch=true;
                entry.localBackgroundH=patchLocalH;
                entry.ratedForFaceLength=patchRated;
                entry.polygon=replacement.cell.fluidPolygon;
                for (std::size_t i=0;i<replacement.vertices.size();++i) {
                    const auto from=replacement.vertices[i].stableId;
                    const auto to=replacement.vertices[
                        (i+1U)%replacement.vertices.size()].stableId;
                    const auto endpoints=std::minmax(from,to);
                    entry.loop.push_back(from);
                    entry.physicalBoundaryFace.push_back(
                        physicalLocks.contains({endpoints.first,endpoints.second}));
                }
                candidateScope.push_back(std::move(entry));
            }
            const auto candidateLocal=evaluatePatchLocalQuality2D(
                candidateScope,minimumFaceOverLocalH,{},tol,hardLimits);
            ++result.localQualityEvaluationCount;
            if (!patchLocalTerminationQualityNoWorse2D(candidateLocal,*baseLocal))
                return false;
            const std::size_t baseHard=baseLocal->hardVolumeRatioCount+
                                       baseLocal->hardFaceWeightCount;
            const std::size_t candidateHard=candidateLocal.hardVolumeRatioCount+
                                            candidateLocal.hardFaceWeightCount;
            if (candidateHard>=baseHard) return false;
            Candidate accepted{first,second,selected,transaction,
                std::move(replacements),patchLocalTerminationRank2D(
                    *baseLocal,candidateLocal,first,second)};
            if (eligibleForWinner &&
                (!winner || patchLocalTerminationRankBetter2D(
                    accepted.rank,winner->rank)))
                winner=std::move(accepted);
            return true;
        };
        bool agglomerationCanImprove=false;
        if (selected.size()==2U &&
            strictlyConvex(repairUnion) &&
            !underDeterminedBoundaryCell(repairUnion,domain,boundary,tol)) {
            const auto locked=lockPolygon(repairUnion,transaction);
            TopologyReplacementCell2D replacement;
            replacement.cell=sourceCell(first,locked);
            replacement.cell.sourceLineage=mergedLineage(
                topology.cells[first].sourceLineage,
                topology.cells[second].sourceLineage);
            for (const auto& point:locked.vertices)
                replacement.vertices.push_back(identityByPoint.at({point.x,point.y}));
            agglomerationCanImprove=evaluateReplacements(
                {std::move(replacement)},
                false,
                agglomerationMode);
        }
        if (agglomerationMode ||
            agglomerationCanImprove) continue;

        std::vector<std::size_t> lineage;
        for (const auto cell:selected)
            lineage=mergedLineage(lineage,topology.cells[cell].sourceLineage);
        if (candidateMode==TerminationQualityCandidateMode2D::Repartition) {
            auto splits=convexTwoPieceSplits(repairUnion,domain,boundary,tol);
            std::sort(splits.begin(),splits.end(),[&](const auto& lhs,const auto& rhs) {
                const auto splitKey=[&](const auto& split) {
                    const double balance=std::abs(std::log(
                        split.first.area()/split.second.area()));
                    const auto a=identityByPoint.at({split.first.vertices.front().x,
                                                     split.first.vertices.front().y}).stableId;
                    const auto b=identityByPoint.at({split.first.vertices.back().x,
                                                     split.first.vertices.back().y}).stableId;
                    const auto endpoints=std::minmax(a,b);
                    return std::tuple{balance,endpoints.first,endpoints.second};
                };
                return splitKey(lhs)<splitKey(rhs);
            });
            constexpr std::size_t candidatesPerPatch=4U;
            if (splits.size()>candidatesPerPatch) splits.resize(candidatesPerPatch);
            for (const auto& [rawFirst,rawSecond]:splits) {
                const std::array<Polygon2D,2> pieces{
                    lockPolygon(rawFirst,transaction),
                    lockPolygon(rawSecond,transaction)};
                std::vector<TopologyReplacementCell2D> replacements(2U);
                for (std::size_t i=0;i<2U;++i) {
                    replacements[i].cell=sourceCell(selected[i],pieces[i]);
                    replacements[i].cell.sourceLineage=lineage;
                    for (const auto& point:pieces[i].vertices)
                        replacements[i].vertices.push_back(
                            identityByPoint.at({point.x,point.y}));
                }
                evaluateReplacements(std::move(replacements),true,true);
            }
            continue;
        }

        struct ThreePieceSplit {
            std::array<Polygon2D,3> pieces;
            double balance=0.0;
            std::array<StableVertexId2D,4> diagonals{};
        };
        std::vector<ThreePieceSplit> threePieceSplits;
        for (const auto& pieces:convexThreePieceFanSplits(
                 repairUnion,domain,boundary,tol)) {
            ThreePieceSplit candidate{pieces};
            const std::array<double,3> areas{
                pieces[0].area(),pieces[1].area(),pieces[2].area()};
            const auto [minimum,maximum]=std::minmax_element(
                areas.begin(),areas.end());
            candidate.balance=*maximum/(*minimum+
                std::numeric_limits<double>::min());
            const auto diagonalKey=[&](const Polygon2D& piece) {
                const auto a=identityByPoint.at({piece.vertices.front().x,
                                                 piece.vertices.front().y}).stableId;
                const auto b=identityByPoint.at({piece.vertices.back().x,
                                                 piece.vertices.back().y}).stableId;
                // minmax returns references; the key must outlive these locals.
                return std::pair<StableVertexId2D,StableVertexId2D>{
                    std::min(a,b),std::max(a,b)};
            };
            const auto firstDiagonal=diagonalKey(pieces[0]);
            const auto secondDiagonal=diagonalKey(pieces[2]);
            candidate.diagonals={firstDiagonal.first,firstDiagonal.second,
                                 secondDiagonal.first,secondDiagonal.second};
            threePieceSplits.push_back(std::move(candidate));
        }
        std::sort(threePieceSplits.begin(),threePieceSplits.end(),
                  [](const auto& lhs,const auto& rhs) {
            return std::tie(lhs.balance,lhs.diagonals)<
                   std::tie(rhs.balance,rhs.diagonals);
        });
        constexpr std::size_t groupedCandidatesPerPatch=16U;
        if (threePieceSplits.size()>groupedCandidatesPerPatch)
            threePieceSplits.resize(groupedCandidatesPerPatch);
        for (const auto& split:threePieceSplits) {
            std::vector<TopologyReplacementCell2D> replacements(3U);
            for (std::size_t i=0;i<3U;++i) {
                const auto piece=lockPolygon(split.pieces[i],transaction);
                replacements[i].cell=sourceCell(selected[i],piece);
                replacements[i].cell.sourceLineage=lineage;
                for (const auto& point:piece.vertices)
                    replacements[i].vertices.push_back(
                        identityByPoint.at({point.x,point.y}));
            }
            evaluateReplacements(std::move(replacements),true,true);
        }
    }
    result.candidateGlobalTopologyBuildCount=
        globalTopologyBuildCount2D()-selectionStartGlobalBuilds;
    result.candidateFullGlobalQualityEvaluationCount=
        solverQualityEvaluationCount2D()-selectionStartFullQuality;
    if (!winner) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }

    const auto committed=evaluateTopologyPatchTransactionOracle2D(
        topology,incidence,winner->transaction,baseSources,winner->replacements,
        domain,boundary,{true,0.0},tol);
    result.globalOracleBuildCount+=committed.globalOracleBuildCount;
    if (!committed.accepted || !committed.deltaMatchesGlobalOracle) {
        result.issues.push_back(committed.issues.empty()
            ?"Q3 winner patch transaction failed the global oracle"
            :committed.issues.front());
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    std::map<TopologySourceIdentity2D,std::tuple<double,bool,bool,bool>> metadata;
    for (const auto& cell:topology.cells) {
        metadata[{cell.sourceKey,cell.sourceId}]={
            localBackgroundH[cell.id],ratedCells[cell.id],immutableCells[cell.id],
            terminationCells[cell.id]};
    }
    double winnerLocalH=0.0;
    bool winnerRated=false,winnerTermination=false;
    for (const auto cell:winner->selected) {
        winnerLocalH=std::max(winnerLocalH,localBackgroundH[cell]);
        winnerRated=winnerRated || ratedCells[cell];
        winnerTermination=winnerTermination || terminationCells[cell];
    }
    const std::size_t retainedSources=
        (candidateMode==TerminationQualityCandidateMode2D::Agglomeration ||
         candidateMode==
             TerminationQualityCandidateMode2D::ConstructionAgglomeration)
        ?1U:winner->selected.size();
    for (std::size_t i=0;i<retainedSources;++i) {
        const auto cell=winner->selected[i];
        metadata[{topology.cells[cell].sourceKey,topology.cells[cell].sourceId}]={
            winnerLocalH,winnerRated,false,winnerTermination};
    }
    std::vector<double> candidateH;
    std::vector<bool> candidateRated,candidateImmutable,candidateTermination;
    for (const auto& cell:committed.topology.cells) {
        const auto item=metadata.at({cell.sourceKey,cell.sourceId});
        candidateH.push_back(std::get<0>(item));
        candidateRated.push_back(std::get<1>(item));
        candidateImmutable.push_back(std::get<2>(item));
        candidateTermination.push_back(std::get<3>(item));
    }
    const auto candidateQuality=evaluateSolverQuality2D(committed.topology,{},tol);
    ++result.authoritativeFullQualityEvaluationCount;
    const auto candidateTarget=scoreTarget(committed.topology,candidateRated);
    const auto candidateShort=scoreShort(
        committed.topology,candidateH,candidateRated);
    result.localWinnerMatchesGlobalAuthority=
        (constructionSelection
             ?constructionContextNoWorse(candidateQuality,currentQuality)
             :solverNoWorse(candidateQuality,currentQuality)) &&
        targetBetter(candidateTarget,currentTarget) &&
        shortNoWorse(candidateShort,currentShort);
    if (!result.localWinnerMatchesGlobalAuthority) {
        // Winner-only oracle: fail this transaction closed without retrying a
        // lower-ranked candidate. The already committed input topology remains
        // the valid Q3 result for the outer convergence loop.
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const std::set<std::size_t> selectedSet(
        winner->selected.begin(),winner->selected.end());
    const auto revisionedBase=buildRevisionedTopology2D(topology,incidence);
    bool outsideStable=true;
    for (const auto& cell:topology.cells) if (!selectedSet.contains(cell.id)) {
        const TopologySourceIdentity2D identity{cell.sourceKey,cell.sourceId};
        outsideStable=outsideStable &&
            committed.revisionedTopology.cells.at(identity).vertices==
            revisionedBase.cells.at(identity).vertices;
    }
    result.topology=committed.topology;
    result.immutableCells=std::move(candidateImmutable);
    result.terminationCells=std::move(candidateTermination);
    result.accepted=true;
    result.hardVolumeRatioCountAfter=candidateTarget.volumeCount;
    result.hardFaceWeightCountAfter=candidateTarget.faceWeightCount;
    result.maximumVolumeRatioSeverityAfter=candidateTarget.maximumVolumeSeverity;
    result.maximumFaceWeightSeverityAfter=candidateTarget.maximumFaceWeightSeverity;
    result.hardShortFaceCountAfter=candidateShort.count;
    result.patchOutsideStableIdsUnchanged=outsideStable;
    result.localDeltaMatchesGlobalOracle=committed.deltaMatchesGlobalOracle;
    result.repairSeconds=profileSeconds(repairStart);
    return result;
}

} // namespace cartmesh2d
