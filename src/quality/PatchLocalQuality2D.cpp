#include "cartmesh2d/quality/PatchLocalQuality2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>

namespace cartmesh2d {
namespace {

struct ScopeUse {
    std::size_t cell = 0;
    StableVertexId2D from = 0;
    StableVertexId2D to = 0;
    Point2D a;
    Point2D b;
    bool declaredPhysicalBoundary = false;
};

[[nodiscard]] StableEdgeKey2D edgeKey(StableVertexId2D from,
                                      StableVertexId2D to) {
    const auto endpoints=std::minmax(from,to);
    return {endpoints.first,endpoints.second};
}

[[nodiscard]] double lowerBoundSeverity(double measured,double limit) noexcept {
    return limit/(measured+std::numeric_limits<double>::min());
}

[[nodiscard]] double upperBoundSeverity(double measured,double limit) noexcept {
    return measured/(limit+std::numeric_limits<double>::min());
}

[[nodiscard]] bool noWorseUpper(double after,double before) noexcept {
    return after-before<=1.0e-8*std::max(1.0,std::abs(before));
}

[[nodiscard]] bool noWorseLower(double after,double before) noexcept {
    return before-after<=1.0e-8*std::max(1.0,std::abs(before));
}

} // namespace

PatchLocalQuality2D evaluatePatchLocalQuality2D(
    const std::vector<PatchLocalCell2D>& scope,double minimumFaceOverLocalH,
    const SolverQualityPolicy2D& policy,const TolerancePolicy& tol) {
    PatchLocalQuality2D result;
    result.policy=policy;
    if (scope.empty() || !(minimumFaceOverLocalH>0.0)) {
        result.issues.push_back("patch-local quality requires a non-empty scope and a positive face ratio");
        return result;
    }
    std::vector<SolverCellMetrics2D> metrics(scope.size());
    for (std::size_t index=0;index<scope.size();++index) {
        const auto& cell=scope[index];
        if (cell.loop.size()!=cell.polygon.vertices.size() || cell.loop.size()<3U ||
            cell.physicalBoundaryFace.size()!=cell.loop.size()) {
            result.issues.push_back("patch-local scope cell has inconsistent identity, geometry or boundary flags");
            return result;
        }
        metrics[index]=evaluateSolverCellMetrics2D(cell.polygon,tol);
        if (!metrics[index].valid) {
            result.issues.push_back("patch-local scope cell polygon or centroid is invalid");
            return result;
        }
        if (!cell.inPatch) continue;
        ++result.ratedCellCount;
        result.maxCellAspect=std::max(result.maxCellAspect,metrics[index].hydraulicAspect);
        result.maxConcavityDeg=std::max(result.maxConcavityDeg,metrics[index].maxConcavityDeg);
        result.minInteriorAngleDeg=std::min(result.minInteriorAngleDeg,
                                            metrics[index].minInteriorAngleDeg);
        result.minCompactness=std::min(result.minCompactness,metrics[index].compactness);
        const auto addIssue=[&](double severity) {
            ++result.issueCount;
            result.maximumIssueSeverity=std::max(result.maximumIssueSeverity,severity);
            result.totalIssueSeverity+=severity;
        };
        if (metrics[index].hydraulicAspect>policy.maxCellAspect)
            addIssue(upperBoundSeverity(metrics[index].hydraulicAspect,policy.maxCellAspect));
        if (metrics[index].maxConcavityDeg>policy.maxConcavityDeg)
            addIssue(upperBoundSeverity(metrics[index].maxConcavityDeg,policy.maxConcavityDeg));
        if (metrics[index].minInteriorAngleDeg<policy.minInteriorAngleDeg)
            addIssue(lowerBoundSeverity(metrics[index].minInteriorAngleDeg,
                                        policy.minInteriorAngleDeg));
    }

    std::map<StableEdgeKey2D,std::vector<ScopeUse>> uses;
    for (std::size_t index=0;index<scope.size();++index) {
        const auto& cell=scope[index];
        for (std::size_t i=0;i<cell.loop.size();++i) {
            const auto from=cell.loop[i];
            const auto to=cell.loop[(i+1U)%cell.loop.size()];
            if (from==to) {
                result.issues.push_back("patch-local scope face has identical endpoints");
                return result;
            }
            uses[edgeKey(from,to)].push_back(
                {index,from,to,cell.polygon.vertices[i],
                 cell.polygon.vertices[(i+1U)%cell.loop.size()],
                 cell.physicalBoundaryFace[i]});
        }
    }

    double minimumFaceLength=std::numeric_limits<double>::infinity();
    double minimumRatio=std::numeric_limits<double>::infinity();
    for (const auto& [key,incidences]:uses) {
        (void)key;
        if (incidences.size()>2U) {
            result.issues.push_back("patch-local scope face is non-manifold");
            return result;
        }
        const bool physical=incidences.size()==1U && incidences.front().declaredPhysicalBoundary;
        if (incidences.size()==1U && !physical) {
            // The real neighbour lies outside the scope.  That is only legal
            // when the truncated side is a halo cell; a patch face without its
            // neighbour cannot be evaluated and must fail closed.
            if (scope[incidences.front().cell].inPatch) {
                result.issues.push_back("patch-local scope truncates a patch face without its neighbour");
                return result;
            }
            continue;
        }
        if (incidences.size()==2U &&
            (incidences[0].from!=incidences[1].to || incidences[0].to!=incidences[1].from)) {
            result.issues.push_back("patch-local scope internal face incidences are not opposite");
            return result;
        }
        const bool rated=std::any_of(incidences.begin(),incidences.end(),
            [&](const ScopeUse& use) { return scope[use.cell].inPatch; });
        if (!rated) continue;

        const auto& owner=incidences.front();
        if (incidences.size()==2U &&
            (owner.a.x!=incidences.back().b.x || owner.a.y!=incidences.back().b.y ||
             owner.b.x!=incidences.back().a.x || owner.b.y!=incidences.back().a.y)) {
            result.issues.push_back("patch-local scope internal face endpoints disagree between its two cells");
            return result;
        }
        const Point2D a=owner.a;
        const Point2D b=owner.b;
        const double faceLength=std::hypot(b.x-a.x,b.y-a.y);
        ++result.ratedFaceCount;
        minimumFaceLength=std::min(minimumFaceLength,faceLength);
        const auto addIssue=[&](double severity) {
            ++result.issueCount;
            result.maximumIssueSeverity=std::max(result.maximumIssueSeverity,severity);
            result.totalIssueSeverity+=severity;
        };
        if (faceLength<=policy.minFaceLength)
            addIssue(lowerBoundSeverity(faceLength,policy.minFaceLength));

        double localH=0.0;
        for (const auto& use:incidences) {
            if (scope[use.cell].ratedForFaceLength)
                localH=std::max(localH,scope[use.cell].localBackgroundH);
        }
        if (localH>0.0) {
            const double ratio=faceLength/localH;
            minimumRatio=std::min(minimumRatio,ratio);
            if (ratio<minimumFaceOverLocalH) {
                ++result.hardShortFaceCount;
                const double severity=minimumFaceOverLocalH/
                    (ratio+std::numeric_limits<double>::min());
                result.maximumShortFaceSeverity=
                    std::max(result.maximumShortFaceSeverity,severity);
                result.totalShortFaceSeverity+=severity;
            }
        }

        if (physical) {
            const double skewness=evaluateSolverBoundaryFaceSkewness2D(
                a,b,metrics[owner.cell].centroid);
            result.maxBoundarySkewness=std::max(result.maxBoundarySkewness,skewness);
            if (skewness>policy.maxBoundarySkewness)
                addIssue(upperBoundSeverity(skewness,policy.maxBoundarySkewness));
            continue;
        }
        const auto& neighbourMetrics=metrics[incidences.back().cell];
        const auto face=evaluateSolverInternalFaceMetrics2D(
            a,b,metrics[owner.cell].centroid,neighbourMetrics.centroid,
            metrics[owner.cell].area,neighbourMetrics.area,tol);
        if (!face.orientationValid) continue;
        ++result.ratedInternalFaceCount;
        result.maxNonOrthogonalityDeg=std::max(result.maxNonOrthogonalityDeg,
                                               face.nonOrthogonalityDeg);
        result.maxInternalSkewness=std::max(result.maxInternalSkewness,face.skewness);
        result.minFaceWeight=std::min(result.minFaceWeight,face.faceWeight);
        result.minVolumeRatio=std::min(result.minVolumeRatio,face.volumeRatio);
        if (face.nonOrthogonalityDeg>policy.maxNonOrthogonalityDeg)
            addIssue(upperBoundSeverity(face.nonOrthogonalityDeg,
                                        policy.maxNonOrthogonalityDeg));
        if (face.skewness>policy.maxInternalSkewness)
            addIssue(upperBoundSeverity(face.skewness,policy.maxInternalSkewness));
        if (face.faceWeight<policy.minFaceWeight)
            addIssue(lowerBoundSeverity(face.faceWeight,policy.minFaceWeight));
        if (face.volumeRatio<policy.minVolumeRatio)
            addIssue(lowerBoundSeverity(face.volumeRatio,policy.minVolumeRatio));
    }
    result.minFaceLength=std::isfinite(minimumFaceLength)?minimumFaceLength:0.0;
    result.minimumFaceOverLocalH=std::isfinite(minimumRatio)?minimumRatio:0.0;
    if (result.ratedCellCount==0U)
        result.issues.push_back("patch-local scope contains no in-patch cell");
    return result;
}

bool patchLocalQualityNoWorse2D(const PatchLocalQuality2D& candidate,
                               const PatchLocalQuality2D& base) noexcept {
    if (!candidate.valid() || !base.valid()) return false;
    const bool solverNoWorse=
        noWorseUpper(candidate.maxNonOrthogonalityDeg,base.maxNonOrthogonalityDeg) &&
        noWorseUpper(candidate.maxInternalSkewness,base.maxInternalSkewness) &&
        noWorseUpper(candidate.maxBoundarySkewness,base.maxBoundarySkewness) &&
        noWorseUpper(candidate.maxConcavityDeg,base.maxConcavityDeg) &&
        noWorseUpper(candidate.maxCellAspect,base.maxCellAspect) &&
        noWorseLower(candidate.minInteriorAngleDeg,base.minInteriorAngleDeg) &&
        noWorseLower(candidate.minFaceLength,base.minFaceLength) &&
        noWorseLower(candidate.minFaceWeight,base.minFaceWeight) &&
        noWorseLower(candidate.minVolumeRatio,base.minVolumeRatio) &&
        noWorseLower(candidate.minCompactness,base.minCompactness);
    return solverNoWorse && detail::patchLocalShortFaceImprovementImpliesGlobal(
        candidate.hardShortFaceCount,candidate.maximumShortFaceSeverity,
        candidate.totalShortFaceSeverity,base.hardShortFaceCount,
        base.maximumShortFaceSeverity,base.totalShortFaceSeverity);
}

PatchLocalRank2D patchLocalRank2D(const PatchLocalQuality2D& base,
                                 const PatchLocalQuality2D& candidate,
                                 std::size_t firstCellId,
                                 std::size_t secondCellId) noexcept {
    PatchLocalRank2D rank;
    const auto signedDelta=[](std::size_t after,std::size_t before) {
        return static_cast<std::ptrdiff_t>(after)-static_cast<std::ptrdiff_t>(before);
    };
    rank.hardShortFaceDelta=signedDelta(candidate.hardShortFaceCount,
                                        base.hardShortFaceCount);
    rank.maximumShortFaceSeverityDelta=
        candidate.maximumShortFaceSeverity-base.maximumShortFaceSeverity;
    rank.totalShortFaceSeverityDelta=
        candidate.totalShortFaceSeverity-base.totalShortFaceSeverity;
    rank.issueDelta=signedDelta(candidate.issueCount,base.issueCount);
    rank.maximumIssueSeverityDelta=
        candidate.maximumIssueSeverity-base.maximumIssueSeverity;
    rank.totalIssueSeverityDelta=
        candidate.totalIssueSeverity-base.totalIssueSeverity;
    rank.maximumShortFaceSeverity=candidate.maximumShortFaceSeverity;
    rank.maxNonOrthogonalityDeg=candidate.maxNonOrthogonalityDeg;
    rank.maxInternalSkewness=candidate.maxInternalSkewness;
    rank.maxCellAspect=candidate.maxCellAspect;
    rank.firstCellId=firstCellId;
    rank.secondCellId=secondCellId;
    return rank;
}

bool patchLocalRankBetter2D(const PatchLocalRank2D& candidate,
                           const PatchLocalRank2D& current) noexcept {
    return detail::patchLocalRankKeyLess(candidate,current);
}

PatchLocalScope2D buildPatchLocalScope2D(
    const TopologyMesh2D& topology,const EdgeIncidenceStore2D& incidence,
    const std::vector<std::size_t>& sortedPatchCellIds,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells) {
    PatchLocalScope2D result;
    if (!topology.valid() || !incidence.valid() ||
        incidence.stableVertexIds.size()!=topology.vertices.size() ||
        localBackgroundH.size()!=topology.cells.size() ||
        ratedCells.size()!=topology.cells.size() || sortedPatchCellIds.empty() ||
        !std::is_sorted(sortedPatchCellIds.begin(),sortedPatchCellIds.end()) ||
        std::adjacent_find(sortedPatchCellIds.begin(),sortedPatchCellIds.end())!=
            sortedPatchCellIds.end() ||
        sortedPatchCellIds.back()>=topology.cells.size()) {
        result.issues.push_back("patch-local scope requires valid aligned inputs");
        return result;
    }
    const auto inPatch=[&](std::size_t cell) {
        return std::binary_search(sortedPatchCellIds.begin(),sortedPatchCellIds.end(),cell);
    };
    std::vector<std::size_t> scopeCells=sortedPatchCellIds;
    for (const auto cell:sortedPatchCellIds) {
        for (const auto edgeId:topology.cells[cell].edges) {
            const auto& edge=topology.edges[edgeId];
            if (!edge.neighbour) continue;
            const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
            if (!inPatch(other)) scopeCells.push_back(other);
        }
    }
    std::sort(scopeCells.begin(),scopeCells.end());
    scopeCells.erase(std::unique(scopeCells.begin(),scopeCells.end()),scopeCells.end());
    result.topologyCellIds=scopeCells;
    result.cells.reserve(scopeCells.size());
    for (const auto cellId:scopeCells) {
        const auto& cell=topology.cells[cellId];
        PatchLocalCell2D entry;
        entry.localBackgroundH=localBackgroundH[cellId];
        entry.ratedForFaceLength=ratedCells[cellId];
        entry.inPatch=inPatch(cellId);
        entry.loop.reserve(cell.vertices.size());
        entry.polygon.vertices.reserve(cell.vertices.size());
        entry.physicalBoundaryFace.reserve(cell.edges.size());
        for (std::size_t local=0;local<cell.vertices.size();++local) {
            entry.loop.push_back(incidence.stableVertexIds[cell.vertices[local]]);
            entry.polygon.vertices.push_back(topology.vertices[cell.vertices[local]].point);
            entry.physicalBoundaryFace.push_back(
                !topology.edges[cell.edges[local]].neighbour.has_value());
        }
        result.cells.push_back(std::move(entry));
    }
    return result;
}

bool patchLocalQualityMatches2D(const PatchLocalQuality2D& lhs,
                               const PatchLocalQuality2D& rhs) noexcept {
    if (lhs.valid()!=rhs.valid()) return false;
    const auto sameSum=[](double a,double b) {
        return std::abs(a-b)<=1.0e-12*std::max({1.0,std::abs(a),std::abs(b)});
    };
    return lhs.ratedCellCount==rhs.ratedCellCount &&
           lhs.ratedFaceCount==rhs.ratedFaceCount &&
           lhs.ratedInternalFaceCount==rhs.ratedInternalFaceCount &&
           lhs.hardShortFaceCount==rhs.hardShortFaceCount &&
           lhs.issueCount==rhs.issueCount &&
           lhs.minimumFaceOverLocalH==rhs.minimumFaceOverLocalH &&
           lhs.maximumShortFaceSeverity==rhs.maximumShortFaceSeverity &&
           lhs.maxNonOrthogonalityDeg==rhs.maxNonOrthogonalityDeg &&
           lhs.maxInternalSkewness==rhs.maxInternalSkewness &&
           lhs.maxBoundarySkewness==rhs.maxBoundarySkewness &&
           lhs.maxConcavityDeg==rhs.maxConcavityDeg &&
           lhs.maxCellAspect==rhs.maxCellAspect &&
           lhs.minInteriorAngleDeg==rhs.minInteriorAngleDeg &&
           lhs.minFaceLength==rhs.minFaceLength &&
           lhs.minFaceWeight==rhs.minFaceWeight &&
           lhs.minVolumeRatio==rhs.minVolumeRatio &&
           lhs.minCompactness==rhs.minCompactness &&
           lhs.maximumIssueSeverity==rhs.maximumIssueSeverity &&
           sameSum(lhs.totalShortFaceSeverity,rhs.totalShortFaceSeverity) &&
           sameSum(lhs.totalIssueSeverity,rhs.totalIssueSeverity);
}

} // namespace cartmesh2d
