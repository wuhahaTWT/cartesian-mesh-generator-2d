#include "SolverTopologyDetail2D.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>

namespace cartmesh2d {
using namespace solver_topology_detail;

namespace {

// Lightweight copy of an unchanged topology cell into a CutCell2D ready for
// buildGlobalTopology. Preserves exact canonical vertex handles when present so
// unchanged cells skip polygon centroid/area recomputation and registry interning.
[[nodiscard]] CutCell2D makeUnchangedCell(const TopologyMesh2D& mesh, std::size_t id) {
    const auto& old = mesh.cells[id];
    CutCell2D cell;
    cell.sourceId = old.sourceId;
    cell.sourceKey = old.sourceKey;
    cell.sourceLineage = old.sourceLineage;
    cell.refinementLineageKeys = old.refinementLineageKeys;
    cell.kind = CutCellKind::Cut;
    cell.area = old.geometryArea;
    cell.areaFraction = 1.0;
    cell.fluidPolygon.vertices.reserve(old.vertices.size());
    const bool hasHandles = mesh.constructionRegistry &&
                            mesh.canonicalVertexIds.size() == mesh.vertices.size();
    if (hasHandles) {
        cell.canonicalVertexIds.reserve(old.vertices.size());
        cell.canonicalRegistry = mesh.constructionRegistry.get();
    }
    AABB2D box{{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
               {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}};
    for (const auto vid : old.vertices) {
        const auto& pt = mesh.vertices[vid].point;
        cell.fluidPolygon.vertices.push_back(pt);
        if (hasHandles) cell.canonicalVertexIds.push_back(mesh.canonicalVertexIds[vid]);
        box.min.x = std::min(box.min.x, pt.x);
        box.min.y = std::min(box.min.y, pt.y);
        box.max.x = std::max(box.max.x, pt.x);
        box.max.y = std::max(box.max.y, pt.y);
    }
    cell.backgroundBounds = box;
    return cell;
}

template<class SimplifyUnion>
void repartitionAgglomeratedCells(SolverShapeAgglomerationResult2D& result,
    std::vector<double>& referenceArea,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const SolverShapeAgglomerationPolicy2D& policy,
    const TolerancePolicy& tol,const SimplifyUnion& simplifyUnion) {
    // Restore input-grid resolution after merging, and repartition sharp-tip
    // cut cells whose centroid cannot support all their faces. Only the cell
    // and its direct neighbours participate; no physical boundary is moved.
    std::optional<SolverQualityReport2D> cachedLegacy;
    for (std::size_t pass=0;pass<policy.maximumPasses;++pass) {
        const auto& mesh=result.topology;
        struct SplitProposal {
            std::vector<std::size_t> members,halo;
            std::vector<Polygon2D> pieces;
            LocalQualityRank2D rank;
            double referenceArea;
            bool improvesTarget=false;
        };
        std::set<std::size_t> issueCells;
        for (const auto& issue:result.after.issues) {
            if (issue.cellId<mesh.cells.size()) issueCells.insert(issue.cellId);
            const bool faceIssue=internalInterfaceIssue(issue.code) ||
                issue.code==SolverQualityIssueCode2D::InvalidFaceGeometry ||
                issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness ||
                issue.code==SolverQualityIssueCode2D::ShortFace;
            if (faceIssue && issue.edgeId<mesh.edges.size()) {
                const auto& edge=mesh.edges[issue.edgeId];
                if (edge.owner<mesh.cells.size()) issueCells.insert(edge.owner);
                if (edge.neighbour && *edge.neighbour<mesh.cells.size())
                    issueCells.insert(*edge.neighbour);
            }
        }
        const auto usable=[&](const Polygon2D& polygon) {
            const auto centre=polygon.centroid(tol);
            if (!centre) return false;
            for (std::size_t i=0;i<polygon.vertices.size();++i)
                if (!evaluateSolverFaceGeometry2D(polygon.vertices[i],
                    polygon.vertices[(i+1)%polygon.vertices.size()],*centre).valid) return false;
            return true;
        };

        const auto splitChoices=[&](const Polygon2D& polygon,
                                    const std::vector<Point2D>& anchorVertices) {
            std::vector<std::vector<Polygon2D>> choices;
            const auto n=polygon.vertices.size();
            const auto add=[&](std::size_t start,std::size_t edge,const Point2D& q,bool vertexEnd) {
                Polygon2D first,second;
                for (std::size_t k=start;;k=(k+1)%n) {
                    first.vertices.push_back(polygon.vertices[k]);if (k==edge) break;
                }
                if (!vertexEnd) first.vertices.push_back(q);
                second.vertices.push_back(q);
                for (std::size_t k=(edge+1)%n;k!=start;k=(k+1)%n)
                    second.vertices.push_back(polygon.vertices[k]);
                second.vertices.push_back(polygon.vertices[start]);
                if (first.vertices.size()<3 || second.vertices.size()<3 ||
                    !(first.signedArea()>0) || !(second.signedArea()>0) ||
                    !BoundaryLoop(first.vertices).diagnose(tol).valid() ||
                    !BoundaryLoop(second.vertices).diagnose(tol).valid()) return;
                if (std::abs(first.area()+second.area()-polygon.area())>
                    tol.absolute*tol.absolute+tol.relative*polygon.area()) return;
                choices.push_back({std::move(first),std::move(second)});
            };
            for (std::size_t start=0;start<n;++start) {
                for (std::size_t end=start+2;end<n;++end) {
                    if (start==0 && end+1==n) continue;
                    add(start,end,polygon.vertices[end],true);
                }
                const auto p=polygon.vertices[start];
                auto incoming=p-polygon.vertices[(start+n-1)%n];
                auto outgoing=polygon.vertices[(start+1)%n]-p;
                const double inNorm=std::sqrt(squaredNorm(incoming));
                const double outNorm=std::sqrt(squaredNorm(outgoing));
                if (!(inNorm>0.0) || !(outNorm>0.0)) continue;
                incoming=incoming*(1.0/inNorm);
                outgoing=outgoing*(1.0/outNorm);
                const Vector2D firstNormal{-incoming.y,incoming.x};
                const Vector2D secondNormal{-outgoing.y,outgoing.x};
                const Vector2D bisector{firstNormal.x+secondNormal.x,firstNormal.y+secondNormal.y};
                const std::array<Vector2D,6> directions{
                    firstNormal,secondNormal,bisector,bisector*(-1.0),
                    incoming,outgoing*(-1.0)};
                for (const auto rawDirection:directions) {
                    const double dirNorm=std::sqrt(squaredNorm(rawDirection));
                    if (!(dirNorm>1.0e-12)) continue;
                    const Vector2D direction=rawDirection*(1.0/dirNorm);
                    std::optional<std::pair<std::size_t,Point2D>> hit;
                    double closest=std::numeric_limits<double>::infinity();
                    for (std::size_t edge=0;edge<n;++edge) {
                        if (edge==start || (edge+1)%n==start) continue;
                        const auto a=polygon.vertices[edge];
                        const auto b=polygon.vertices[(edge+1)%n];
                        const auto tangent=b-a;
                        const double edgeLength=std::sqrt(squaredNorm(tangent));
                        const double determinant=cross(direction,tangent);
                        if (std::abs(determinant)<=1.0e-15*edgeLength) continue;
                        const double distance=cross(a-p,tangent)/determinant;
                        const double fraction=cross(a-p,direction)/determinant;
                        if (distance<=0.0 || fraction<=1.0e-9 || fraction>=1.0-1.0e-9 || distance>=closest) continue;
                        Point2D point=a+tangent*fraction;
                        // Snap ray hits within roundoff of an existing shared topology
                        // vertex (which simplifyUnion may have dropped as collinear) so
                        // global edge matching never sees ULP-mismatched T-junctions.
                        const double snapLimit=std::max(tol.absolute,tol.relative*edgeLength);
                        const double snapLimitSq=snapLimit*snapLimit;
                        double bestDistSq=snapLimitSq;
                        for (const auto& anchor:anchorVertices) {
                            const double distSq=squaredNorm(anchor-point);
                            if (distSq<=bestDistSq) {
                                bestDistSq=distSq;
                                point=anchor;
                            }
                        }
                        if (squaredNorm(point-a)<=snapLimitSq || squaredNorm(point-b)<=snapLimitSq) continue;
                        closest=distance;hit={{edge,point}};
                    }
                    if (hit) add(start,hit->first,hit->second,false);
                }
            }
            return choices;
        };

        const auto propose=[&](const std::vector<std::size_t>& members,const Polygon2D& polygon)
            ->std::optional<SplitProposal> {
            std::vector<std::size_t> halo;
            std::vector<Polygon2D> originals;
            double resolutionArea=0;
            bool badGeometry=false;
            for (const auto id:members) {
                if (!result.immutableCells.empty() && result.immutableCells[id]) return std::nullopt;
                const auto neighbours=cellPairHalo(mesh,id,id);
                halo.insert(halo.end(),neighbours.begin(),neighbours.end());
                originals.push_back(topologyCellPolygon(mesh,id));
                badGeometry=badGeometry || !usable(originals.back());
                resolutionArea=std::max(resolutionArea,referenceArea[id]);
            }
            std::sort(halo.begin(),halo.end());halo.erase(std::unique(halo.begin(),halo.end()),halo.end());
            const double areaTolerance=tol.absolute*tol.absolute+tol.relative*resolutionArea;
            const bool oversized=members.size()==1 && polygon.area()>resolutionArea+areaTolerance;
            const auto scope=makeLocalPatchQualityScope(mesh,halo,members,tol);
            if (!scope.valid()) return std::nullopt;
            const auto beforePatch=scope.buildPatch(originals);
            const auto before=beforePatch?scope.rankFromPatch(*beforePatch,policy.target):std::nullopt;
            const auto beforeLegacy=(!badGeometry && oversized && beforePatch)
                ? scope.rankFromPatch(*beforePatch,{})
                : std::nullopt;
            std::optional<LocalQualityRank2D> bestRank;
            std::vector<Polygon2D> bestPieces;
            double bestLargest=std::numeric_limits<double>::infinity();
            const auto consider=[&](std::vector<Polygon2D> pieces) {
                if (pieces.size()<std::max<std::size_t>(2,members.size())) return;
                double largest=0.0,sum=0.0;
                for (const auto& piece:pieces) {
                    if (!usable(piece)) return;
                    largest=std::max(largest,piece.area());sum+=piece.area();
                }
                if (std::abs(sum-polygon.area())>tol.absolute*tol.absolute+
                    tol.relative*polygon.area()) return;
                if (members.size()>1 && largest>resolutionArea+areaTolerance) return;
                if (oversized && largest>=polygon.area()-areaTolerance) return;
                const double excess=std::max(0.0,largest-resolutionArea);
                const double bestExcess=std::max(0.0,bestLargest-resolutionArea);
                if (oversized && bestRank && excess>bestExcess) return;
                const auto candidatePatch=scope.buildPatch(pieces);
                if (!candidatePatch) return;
                auto rank=scope.rankFromPatch(*candidatePatch,policy.target);
                if (!rank) return;
                if (!badGeometry && !oversized && (!before ||
                    !betterQualityScore(rank->issues,before->issues))) return;
                // Size restoration must not create the thin, poorly connected
                // triangles that agglomeration was introduced to remove.
                if (!badGeometry && before && betterQualityScore(before->issues,rank->issues)) return;
                if (!badGeometry && oversized) {
                    const auto candidateLegacy=scope.rankFromPatch(*candidatePatch,{});
                    if (!candidateLegacy ||
                        (beforeLegacy && betterQualityScore(beforeLegacy->issues,candidateLegacy->issues))) return;
                }
                if (!bestRank || (oversized && excess<bestExcess) ||
                    ((!oversized || excess==bestExcess) && betterLocalQualityRank(*rank,*bestRank))) {
                    bestRank=*rank;bestPieces=std::move(pieces);bestLargest=largest;
                }
            };
            std::vector<Point2D> anchorVertices;
            for (const auto& orig:originals)
                anchorVertices.insert(anchorVertices.end(),orig.vertices.begin(),orig.vertices.end());
            for (const auto& pieces:splitChoices(polygon,anchorVertices)) consider(pieces);
            if (badGeometry || bestPieces.empty() ||
                (bestRank && bestRank->issues.issueCount>0U && polygon.vertices.size()<=12U)) {
                for (const auto& pieces:convexThreePieceFanSplits(polygon,domain,boundary,tol))
                    consider({pieces.begin(),pieces.end()});
                const auto pieces=partitionOnePolygon(polygon,domain,boundary,tol);
                if (pieces) consider(*pieces);
            }
            if (!bestRank) return std::nullopt;
            const bool improvesTarget=before && betterQualityScore(bestRank->issues,before->issues);
            if (before) rankRelativeToBase(*bestRank,*before);
            return SplitProposal{members,halo,std::move(bestPieces),*bestRank,resolutionArea,improvesTarget};
        };
        const auto commitProposals=[&](std::vector<SplitProposal>& proposals,
                                       bool fixingBadGeometry,
                                       bool restoringSizeOnly) -> bool {
            if (proposals.empty()) return false;
            std::sort(proposals.begin(),proposals.end(),[](const SplitProposal& a,const SplitProposal& b) {
                if (betterLocalQualityRank(a.rank,b.rank)) return true;
                if (betterLocalQualityRank(b.rank,a.rank)) return false;
                return a.members<b.members;
            });
            std::vector<std::vector<std::size_t>> halos;
            for (const auto& proposal:proposals) halos.push_back(proposal.halo);
            const auto selected=selectIndependentSolverRepairPatches2D(halos);
            if (restoringSizeOnly && !cachedLegacy) cachedLegacy=evaluateSolverQuality2D(mesh,{},tol);
            const auto targetScoreBefore=qualityScore(result.after);
            const auto legacyScoreBefore=restoringSizeOnly?qualityScore(*cachedLegacy):QualityScore2D{};
            const auto invalidFaceCount=[](const SolverQualityReport2D& report) {
                std::size_t count=0;
                for (const auto& issue:report.issues)
                    if (issue.code==SolverQualityIssueCode2D::InvalidFaceGeometry) ++count;
                return count;
            };
            const std::size_t badFacesBefore=invalidFaceCount(result.after);
            const auto trySubset=[&](const std::vector<std::size_t>& indices) -> bool {
                if (indices.empty()) return false;
                std::map<std::size_t,const SplitProposal*> replacements;
                std::set<std::size_t> removed;
                std::size_t changed=0,added=0;
                for (const auto index:indices) {
                    const auto& proposal=proposals[index];
                    replacements.emplace(proposal.members.front(),&proposal);
                    removed.insert(proposal.members.begin()+1,proposal.members.end());
                    changed+=proposal.members.size();
                    added+=proposal.pieces.size()-proposal.members.size();
                }
                std::size_t nextSource=0;
                for (const auto& cell:mesh.cells) nextSource=std::max(nextSource,cell.sourceId+1);
                std::vector<CutCell2D> cells;
                cells.reserve(mesh.cells.size()+added);
                std::map<std::pair<std::uint64_t,std::size_t>,std::pair<double,bool>> metadata;
                for (std::size_t id=0;id<mesh.cells.size();++id) {
                    if (removed.contains(id)) continue;
                    const auto& old=mesh.cells[id];
                    const auto replacement=replacements.find(id);
                    const bool isLocked=!result.immutableCells.empty() && result.immutableCells[id];
                    if (replacement==replacements.end()) {
                        auto cell=makeUnchangedCell(mesh,id);
                        metadata[{cell.sourceKey,cell.sourceId}]={referenceArea[id],isLocked};
                        cells.push_back(std::move(cell));
                        continue;
                    }
                    const auto& pieces=replacement->second->pieces;
                    auto lineage=old.sourceLineage;
                    auto refinementKeys=old.refinementLineageKeys;
                    const double resolutionArea=replacement->second->referenceArea;
                    for (const auto member:replacement->second->members) {
                        lineage=mergedLineage(lineage,mesh.cells[member].sourceLineage);
                        const auto& keys=mesh.cells[member].refinementLineageKeys;
                        refinementKeys.insert(refinementKeys.end(),keys.begin(),keys.end());
                    }
                    std::sort(refinementKeys.begin(),refinementKeys.end());
                    refinementKeys.erase(std::unique(refinementKeys.begin(),refinementKeys.end()),refinementKeys.end());
                    for (std::size_t piece=0;piece<pieces.size();++piece) {
                        auto cell=makeCell(piece==0?old.sourceId:nextSource++,pieces[piece],tol,lineage);
                        cell.sourceKey=old.sourceKey;cell.refinementLineageKeys=refinementKeys;
                        metadata[{cell.sourceKey,cell.sourceId}]={resolutionArea,isLocked};
                        cells.push_back(std::move(cell));
                    }
                }
                auto candidate=buildGlobalTopology(cells,domain,boundary,tol,
                    mesh.constructionRegistry,mesh.embeddedBoundaryFragments);
                if (!candidate.valid()) return false;
                const auto afterTarget=evaluateSolverQuality2D(candidate,policy.target,tol);
                std::optional<SolverQualityReport2D> afterLegacy;
                const std::size_t badFacesAfter=invalidFaceCount(afterTarget);
                if (fixingBadGeometry) {
                    if (badFacesAfter>=badFacesBefore) return false;
                } else {
                    const bool targetImproved=betterQualityScore(qualityScore(afterTarget),targetScoreBefore);
                    const bool targetWorse=betterQualityScore(targetScoreBefore,qualityScore(afterTarget));
                    if (badFacesAfter>badFacesBefore || targetWorse ||
                        (!restoringSizeOnly && !targetImproved)) return false;
                    if (restoringSizeOnly) {
                        afterLegacy=evaluateSolverQuality2D(candidate,{},tol);
                        if (betterQualityScore(legacyScoreBefore,qualityScore(*afterLegacy))) return false;
                    }
                }
                std::vector<double> newAreas;
                std::vector<bool> locked;
                newAreas.reserve(candidate.cells.size());
                locked.reserve(candidate.cells.size());
                for (const auto& cell:candidate.cells) {
                    const auto& entry=metadata.at({cell.sourceKey,cell.sourceId});
                    newAreas.push_back(entry.first);locked.push_back(entry.second);
                }
                result.topology=std::move(candidate);referenceArea=std::move(newAreas);
                result.immutableCells=std::move(locked);
                result.repartitionedCellCount+=changed;result.addedCellCount+=added;
                result.after=afterTarget;
                cachedLegacy=std::move(afterLegacy);
                return true;
            };
            // Try full independent batch first, then halved prefixes; if prefix 1
            // fails due to a global halo interaction, scan individual candidates
            // across all sorted proposals so a rejected patch never starves
            // overlapping or unrelated valid repairs.
            for (std::size_t count=selected.size();count>0;count=count==1?0:(count+1)/2) {
                std::vector<std::size_t> subset(selected.begin(),
                    selected.begin()+static_cast<std::ptrdiff_t>(count));
                if (trySubset(subset)) return true;
            }
            const std::size_t singleLimit=std::min<std::size_t>(proposals.size(),12U);
            for (std::size_t k=0;k<singleLimit;++k) {
                if (!selected.empty() && k==selected.front()) continue;
                if (trySubset({k})) return true;
            }
            return false;
        };

        // Strict three-tier priority prevents geometry recovery from admitting
        // unrelated size-restoration splits that worsen global quality.
        std::vector<SplitProposal> badProposals,issueProposals,sizeProposals;
        for (std::size_t id=0;id<mesh.cells.size();++id) {
            if (!result.immutableCells.empty() && result.immutableCells[id]) continue;
            const bool hasIssue=issueCells.contains(id);
            const double areaMargin=tol.absolute*tol.absolute+tol.relative*referenceArea[id];
            // evaluateSolverQuality2D reports every cell with an invalid centroid
            // or nonpositive face normal distance in result.after.issues, so
            // any cell outside issueCells that is within its reference area can
            // be skipped in O(1) without constructing its Polygon2D.
            if (!hasIssue && mesh.cells[id].geometryArea<=referenceArea[id]+areaMargin) continue;
            const auto polygon=topologyCellPolygon(mesh,id);
            const bool badGeometry=hasIssue && !usable(polygon);
            const bool oversized=polygon.area()>referenceArea[id]+areaMargin;
            if (!badGeometry && !hasIssue && !oversized) continue;
            auto best=propose({id},polygon);
            if (badGeometry || hasIssue) {
                std::set<std::size_t> neighbours;
                for (const auto edgeId:mesh.cells[id].edges) {
                    const auto& edge=mesh.edges[edgeId];
                    if (edge.neighbour) neighbours.insert(edge.owner==id?*edge.neighbour:edge.owner);
                }
                for (const auto neighbour:neighbours) {
                    const auto merged=mergeAdjacentPolygonsSimple(
                        polygon,topologyCellPolygon(mesh,neighbour),tol);
                    if (!merged) continue;
                    const std::vector<std::size_t> members{std::min(id,neighbour),std::max(id,neighbour)};
                    auto candidate=propose(members,simplifyUnion(*merged));
                    if (candidate && (!best || betterLocalQualityRank(candidate->rank,best->rank)))
                        best=std::move(candidate);
                }
            }
            if (!best) continue;
            if (badGeometry) badProposals.push_back(std::move(*best));
            else if (hasIssue && best->improvesTarget) issueProposals.push_back(std::move(*best));
            else if (oversized) sizeProposals.push_back(std::move(*best));
        }
        bool committed=false;
        if (!badProposals.empty()) committed=commitProposals(badProposals,true,false);
        if (!committed && !issueProposals.empty()) committed=commitProposals(issueProposals,false,false);
        if (!committed && !sizeProposals.empty()) committed=commitProposals(sizeProposals,false,true);
        if (!committed) break;
    }

}

} // namespace

SolverShapeAgglomerationResult2D agglomerateSolverByShape2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const SolverShapeAgglomerationPolicy2D& policy,const TolerancePolicy& tol) {
    SolverShapeAgglomerationResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    const SolverQualityPolicy2D legacyPolicy;
    const auto& target=policy.target;
    bool finitePolicy=true;
    for (const double value:{target.maxNonOrthogonalityDeg,target.maxInternalSkewness,
        target.maxBoundarySkewness,target.maxConcavityDeg,target.maxCellAspect,
        target.minInteriorAngleDeg,target.minFaceLength,target.minFaceWeight,target.minVolumeRatio})
        finitePolicy=finitePolicy && std::isfinite(value);
    if (!topology.valid() || topology.cells.empty() || !domain.valid(tol) ||
        !boundary.diagnose(tol).valid() ||
        (!immutableCells.empty() && immutableCells.size()!=topology.cells.size()) ||
        !std::isfinite(policy.maximumAreaRatio) || policy.maximumAreaRatio<1 ||
        policy.maximumPasses==0 || policy.maximumPasses>1024 ||
        policy.maximumPatchCells<2 || policy.maximumPatchCells>8 || !finitePolicy ||
        target.maxNonOrthogonalityDeg<=0 || target.maxNonOrthogonalityDeg>legacyPolicy.maxNonOrthogonalityDeg ||
        target.maxInternalSkewness<=0 || target.maxInternalSkewness>legacyPolicy.maxInternalSkewness ||
        target.maxBoundarySkewness<=0 || target.maxBoundarySkewness>legacyPolicy.maxBoundarySkewness ||
        target.maxConcavityDeg<0 || target.maxConcavityDeg>legacyPolicy.maxConcavityDeg ||
        target.maxCellAspect<=0 || target.maxCellAspect>legacyPolicy.maxCellAspect ||
        target.minInteriorAngleDeg<legacyPolicy.minInteriorAngleDeg || target.minInteriorAngleDeg>=180 ||
        target.minFaceLength<legacyPolicy.minFaceLength ||
        target.minFaceWeight<legacyPolicy.minFaceWeight || target.minFaceWeight>.5 ||
        target.minVolumeRatio<legacyPolicy.minVolumeRatio || target.minVolumeRatio>1) {
        result.issues.push_back("shape agglomeration requires valid inputs and finite budgets");
        return result;
    }
    std::vector<double> referenceArea;
    for (auto& cell:result.topology.cells) {
        if (cell.sourceLineage.empty()) cell.sourceLineage={cell.sourceId};
        const auto box=topologyCellPolygon(topology,cell.id).bounds();
        referenceArea.push_back((box.max.x-box.min.x)*(box.max.y-box.min.y));
    }
    result.before=evaluateSolverQuality2D(result.topology,policy.target,tol);
    result.after=result.before;
    std::set<std::pair<double,double>> physicalVertices;
    for (const auto& loop:boundary.loops()) for (const auto& p:loop.vertices())
        physicalVertices.emplace(p.x,p.y);
    if (topology.constructionRegistry) for (const auto& record:topology.constructionRegistry->records())
        if (physicalVertices.contains({record.originalPoint.x,record.originalPoint.y}))
            physicalVertices.emplace(record.canonicalVertex.point.x,record.canonicalVertex.point.y);
    const auto simplifyUnion=[&](Polygon2D polygon) {
        const auto original=polygon;
        // A grid/wall intersection ceases to be a face endpoint after both
        // incident cells are united. Remove only redundant construction
        // vertices, never an input wall vertex. Common partitions with cells
        // outside the union are restored by the authoritative topology build.
        for (std::size_t i=0;i<polygon.vertices.size() && polygon.vertices.size()>3;) {
            const auto n=polygon.vertices.size();
            const auto a=polygon.vertices[(i+n-1)%n],b=polygon.vertices[i],c=polygon.vertices[(i+1)%n];
            const auto u=b-a,v=c-b;
            const double scale=std::sqrt(squaredNorm(u)*squaredNorm(v));
            const double magnitude=std::max({std::abs(a.x),std::abs(a.y),std::abs(b.x),
                std::abs(b.y),std::abs(c.x),std::abs(c.y)});
            const double roundoff=tol.constructionRoundoffScale(magnitude)*
                (std::sqrt(squaredNorm(u))+std::sqrt(squaredNorm(v)));
            const double budget=std::max(64*std::numeric_limits<double>::epsilon()*scale,
                std::min(tol.relative*scale,roundoff));
            if (!physicalVertices.contains({b.x,b.y}) && scale>0 && dot(u,v)>0 &&
                std::abs(cross(u,v))<=budget) {
                polygon.vertices.erase(polygon.vertices.begin()+static_cast<std::ptrdiff_t>(i));
                i=0;
            } else ++i;
        }
        return std::abs(polygon.area()-original.area())<=
            tol.absolute*tol.absolute+tol.relative*original.area()?polygon:original;
    };
    std::optional<SolverQualityReport2D> cachedLegacy;
    for (;result.passes<policy.maximumPasses && !result.after.valid();++result.passes) {
        const auto& mesh=result.topology;
        if (!cachedLegacy) cachedLegacy=evaluateSolverQuality2D(mesh,{},tol);
        const auto& legacy=*cachedLegacy;
        std::set<std::pair<std::size_t,std::size_t>> pairs;
        for (const auto& issue:result.after.issues) {
            std::vector<std::size_t> affected{issue.cellId};
            const bool faceIssue=internalInterfaceIssue(issue.code) ||
                issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness ||
                issue.code==SolverQualityIssueCode2D::ShortFace;
            if (faceIssue && issue.edgeId<mesh.edges.size() && mesh.edges[issue.edgeId].neighbour) {
                affected.push_back(mesh.edges[issue.edgeId].owner);
                affected.push_back(*mesh.edges[issue.edgeId].neighbour);
            }
            for (const auto cell:affected) {
                if (cell>=mesh.cells.size()) continue;
                for (const auto edgeId:mesh.cells[cell].edges) {
                    const auto& edge=mesh.edges[edgeId];
                    if (!edge.neighbour) continue;
                    const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
                    if (!result.immutableCells.empty() &&
                        (result.immutableCells[cell] || result.immutableCells[other])) continue;
                    pairs.emplace(std::min(cell,other),std::max(cell,other));
                }
            }
        }
        struct ShapeProposal {
            std::vector<std::size_t> members,halo;
            Polygon2D polygon;
            LocalQualityRank2D rank;
        };
        std::vector<ShapeProposal> proposals;
        const auto neighbours=[&](const std::vector<std::size_t>& members) {
            std::set<std::size_t> ids(members.begin(),members.end());
            for (const auto cell:members)
                for (const auto edge:mesh.cells[cell].edges) {
                    const auto& e=mesh.edges[edge];
                    if (e.neighbour) {ids.insert(e.owner);ids.insert(*e.neighbour);}
                }
            return std::vector<std::size_t>(ids.begin(),ids.end());
        };
        std::map<std::vector<std::size_t>,Polygon2D> frontier;
        for (const auto& [first,second]:pairs) {
            const auto merged=mergeAdjacentPolygonsSimple(
                topologyCellPolygon(mesh,first),topologyCellPolygon(mesh,second),tol);
            if (merged) frontier.emplace(std::vector<std::size_t>{first,second},*merged);
        }
        // A pair can be stuck although a three/four-cell union improves the
        // whole patch. Grow only when no smaller union improves local quality.
        for (std::size_t width=2;width<=policy.maximumPatchCells && proposals.empty();++width) {
          for (const auto& [members,rawMerged]:frontier) {
            const auto merged=simplifyUnion(rawMerged);
            double areaScale=0;
            for (const auto cell:members) areaScale=std::max(areaScale,referenceArea[cell]);
            if (merged.area()>policy.maximumAreaRatio*areaScale) continue;
            const auto centroid=merged.centroid(tol);
            if (!centroid) continue;
            bool admissible=true;
            for (std::size_t i=0;i<merged.vertices.size();++i)
                admissible=admissible && evaluateSolverFaceGeometry2D(merged.vertices[i],
                    merged.vertices[(i+1)%merged.vertices.size()],*centroid).valid;
            if (!admissible) continue;
            ++result.candidateCount;
            const auto halo=neighbours(members);
            std::vector<Polygon2D> originals;
            originals.reserve(members.size());
            for (const auto cell:members) originals.push_back(topologyCellPolygon(mesh,cell));
            const auto scope=makeLocalPatchQualityScope(mesh,halo,members,tol);
            if (!scope.valid()) continue;
            const auto beforePatch=scope.buildPatch(originals);
            const auto afterPatch=scope.buildPatch({merged});
            if (!beforePatch || !afterPatch) continue;
            const auto before=scope.rankFromPatch(*beforePatch,policy.target);
            auto after=scope.rankFromPatch(*afterPatch,policy.target);
            if (!before || !after || !betterQualityScore(after->issues,before->issues)) continue;
            const auto oldLegacy=scope.rankFromPatch(*beforePatch,{});
            const auto newLegacy=scope.rankFromPatch(*afterPatch,{});
            if (!oldLegacy || !newLegacy || betterQualityScore(oldLegacy->issues,newLegacy->issues)) continue;
            rankRelativeToBase(*after,*before);
            proposals.push_back({members,halo,merged,*after});
          }
          if (proposals.empty() && width<policy.maximumPatchCells) {
            std::map<std::vector<std::size_t>,Polygon2D> next;
            for (const auto& [members,rawMerged]:frontier) {
              double areaScale=0;
              for (const auto cell:members) areaScale=std::max(areaScale,referenceArea[cell]);
              if (rawMerged.area()>policy.maximumAreaRatio*areaScale) continue;
              const auto halo=neighbours(members);
              for (const auto other:halo) {
                if (std::binary_search(members.begin(),members.end(),other) ||
                    (!result.immutableCells.empty() && result.immutableCells[other])) continue;
                auto grown=members;grown.push_back(other);std::sort(grown.begin(),grown.end());
                if (next.contains(grown)) continue;
                const auto polygon=mergeAdjacentPolygonsSimple(rawMerged,topologyCellPolygon(mesh,other),tol);
                if (polygon) next.emplace(std::move(grown),*polygon);
              }
            }
            frontier=std::move(next);
          }
        }
        std::sort(proposals.begin(),proposals.end(),[](const auto& a,const auto& b) {
            if (betterLocalQualityRank(a.rank,b.rank)) return true;
            if (betterLocalQualityRank(b.rank,a.rank)) return false;
            return a.members<b.members;
        });
        std::vector<std::vector<std::size_t>> halos;
        for (const auto& p:proposals) halos.push_back(p.halo);
        const auto selected=selectIndependentSolverRepairPatches2D(halos);
        bool committed=false;
        for (std::size_t count=selected.size();count>0;count=count==1?0:(count+1)/2) {
            std::vector<bool> removed(mesh.cells.size(),false);
            std::vector<const ShapeProposal*> replacement(mesh.cells.size(),nullptr);
            auto areas=referenceArea;
            std::size_t mergedCount=0;
            for (std::size_t k=0;k<count;++k) {
                const auto& p=proposals[selected[k]];
                replacement[p.members.front()]=&p;
                for (const auto cell:p.members) {
                    removed[cell]=cell!=p.members.front();
                    areas[p.members.front()]=std::max(areas[p.members.front()],areas[cell]);
                }
                mergedCount+=p.members.size()-1;
            }
            std::vector<CutCell2D> cells;
            std::vector<bool> locked;
            std::vector<double> newAreas;
            cells.reserve(mesh.cells.size()-mergedCount);
            locked.reserve(mesh.cells.size()-mergedCount);
            newAreas.reserve(mesh.cells.size()-mergedCount);
            for (std::size_t i=0;i<mesh.cells.size();++i) if (!removed[i]) {
                if (!replacement[i]) {
                    cells.push_back(makeUnchangedCell(mesh,i));
                } else {
                    const auto& old=mesh.cells[i];
                    auto lineage=old.sourceLineage;
                    auto refinementKeys=old.refinementLineageKeys;
                    for (const auto member:replacement[i]->members) {
                        lineage=mergedLineage(lineage,mesh.cells[member].sourceLineage);
                        const auto& keys=mesh.cells[member].refinementLineageKeys;
                        refinementKeys.insert(refinementKeys.end(),keys.begin(),keys.end());
                    }
                    std::sort(refinementKeys.begin(),refinementKeys.end());
                    refinementKeys.erase(std::unique(refinementKeys.begin(),refinementKeys.end()),refinementKeys.end());
                    auto cell=makeCell(old.sourceId,replacement[i]->polygon,tol,std::move(lineage));
                    cell.sourceKey=old.sourceKey;
                    cell.refinementLineageKeys=std::move(refinementKeys);
                    cells.push_back(std::move(cell));
                }
                locked.push_back(!result.immutableCells.empty() && result.immutableCells[i]);
                newAreas.push_back(areas[i]);
            }
            auto candidate=buildGlobalTopology(cells,domain,boundary,tol,
                mesh.constructionRegistry,mesh.embeddedBoundaryFragments);
            if (!candidate.valid()) continue;
            const auto after=evaluateSolverQuality2D(candidate,policy.target,tol);
            const auto newLegacy=evaluateSolverQuality2D(candidate,{},tol);
            if (!betterQualityScore(qualityScore(after),qualityScore(result.after)) ||
                betterQualityScore(qualityScore(legacy),qualityScore(newLegacy))) continue;
            referenceArea=std::move(newAreas);
            result.topology=std::move(candidate);
            result.immutableCells=std::move(locked);
            result.after=after;
            cachedLegacy=newLegacy;
            result.mergedCellCount+=mergedCount;
            for (std::size_t k=0;k<count;++k) result.maximumUnionCellCount=
                std::max(result.maximumUnionCellCount,proposals[selected[k]].members.size());
            committed=true;
            break;
        }
        if (!committed) {result.stopReason="no-improving-union";break;}
    }


    repartitionAgglomeratedCells(result,referenceArea,domain,boundary,policy,tol,simplifyUnion);

    if (result.after.valid()) result.stopReason="target-met";
    else if (result.stopReason.empty()) result.stopReason="pass-limit";
    return result;
}

} // namespace cartmesh2d
