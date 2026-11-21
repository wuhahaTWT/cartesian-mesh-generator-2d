#pragma once

#include "cartmesh2d/quality/PatchLocalQuality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/topology/PatchTransaction2D.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace cartmesh2d::solver_topology_detail {

using ProfileClock = std::chrono::steady_clock;

struct PartitionAttempt2D {
    TopologyMesh2D topology;
    std::vector<std::size_t> sourceForCell;
    std::vector<bool> immutableForCell;
    std::size_t partitionedSourceCellCount = 0;
    std::string error;
};

struct QualityScore2D {
    std::size_t issueCount = 0;
    double maximumSeverity = 0.0;
    double totalSeverity = 0.0;
};

struct LocalQualityRank2D {
    // Source proposals span different neighbourhoods; compare improvement
    // relative to each unchanged base scope before absolute tie-breaks.
    std::ptrdiff_t issueDelta = 0;
    double maximumSeverityDelta = 0.0;
    double totalSeverityDelta = 0.0;
    QualityScore2D issues;
    double maxNonOrthogonality = 0.0;
    double maxInternalSkewness = 0.0;
    double maxBoundarySkewness = 0.0;
    double maxCellAspect = 0.0;
    double negativeMinInteriorAngle = 0.0;
    double negativeMinFaceWeight = 0.0;
    double negativeMinVolumeRatio = 0.0;
};

struct RepartitionBatch2D {
    TopologyMesh2D topology;
    std::vector<bool> immutableCells;
};

struct SourceMergeProposal2D {
    std::size_t first = 0;
    std::size_t second = 0;
    Polygon2D merged;
    LocalQualityRank2D rank;
    std::vector<std::size_t> halo;
};

struct RepartitionProposal2D {
    std::size_t first = 0;
    std::size_t second = 0;
    std::vector<Polygon2D> pieces;
    LocalQualityRank2D rank;
    std::vector<std::size_t> halo;
};

[[nodiscard]] double profileSeconds(ProfileClock::time_point start) noexcept;
[[nodiscard]] bool pointInOrOnTriangle(const Point2D& point,
                                       const Point2D& a,const Point2D& b,
                                       const Point2D& c) noexcept;
[[nodiscard]] bool strictlyConvex(const Polygon2D& polygon) noexcept;
[[nodiscard]] bool cartesianRectangle(const Polygon2D& polygon) noexcept;
[[nodiscard]] bool boundaryEdge(const Point2D& a,const Point2D& b,
                                const Domain2D& domain,
                                const BoundaryRegion2D& boundary,
                                const TolerancePolicy& tol);
[[nodiscard]] unsigned boundaryEdgeCount(const Polygon2D& polygon,
                                         const Domain2D& domain,
                                         const BoundaryRegion2D& boundary,
                                         const TolerancePolicy& tol);
[[nodiscard]] bool underDeterminedBoundaryCell(const Polygon2D& polygon,
                                               const Domain2D& domain,
                                               const BoundaryRegion2D& boundary,
                                               const TolerancePolicy& tol);
[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygons(
    const Polygon2D& first,const Polygon2D& second);
[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygonsSimple(
    const Polygon2D& first,const Polygon2D& second,const TolerancePolicy& tol);
[[nodiscard]] double minimumInteriorAngle(const Polygon2D& polygon);
[[nodiscard]] std::optional<std::vector<Polygon2D>> triangulate(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol);
[[nodiscard]] std::vector<std::size_t> mergedLineage(
    const std::vector<std::size_t>& first,
    const std::vector<std::size_t>& second);
[[nodiscard]] Polygon2D removeArtificialCollinearVertices(
    Polygon2D polygon,const TolerancePolicy& tol,bool constructionRoundoff=false);
[[nodiscard]] CutCell2D makeCell(std::size_t id,const Polygon2D& polygon,
                                 const TolerancePolicy& tol,
                                 std::vector<std::size_t> sourceLineage={});

[[nodiscard]] QualityScore2D qualityScore(const SolverQualityReport2D& quality);
[[nodiscard]] bool betterQualityScore(const QualityScore2D& candidate,
                                      const QualityScore2D& current) noexcept;
[[nodiscard]] bool betterLocalQualityRank(const LocalQualityRank2D& candidate,
                                          const LocalQualityRank2D& current) noexcept;
void rankRelativeToBase(LocalQualityRank2D& rank,const LocalQualityRank2D& base);

[[nodiscard]] std::optional<std::vector<Polygon2D>> partitionOnePolygon(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol);
[[nodiscard]] std::optional<LocalQualityRank2D> localPartitionQualityRank(
    const Polygon2D& outer,const std::vector<Polygon2D>& pieces,
    const TolerancePolicy& tol);
[[nodiscard]] bool internalInterfaceIssue(SolverQualityIssueCode2D code) noexcept;
[[nodiscard]] PartitionAttempt2D partitionSourcePolygons(
    const std::vector<Polygon2D>& sourcePolygons,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile=nullptr,bool candidate=false,
    const std::vector<bool>& immutableSources={},
    const std::vector<bool>& preserveSources={},
    const std::vector<std::vector<std::size_t>>& sourceLineages={},
    std::shared_ptr<IntersectionRegistry2D> registry={},
    const std::vector<Segment2D>& knownEmbeddedBoundary={});

[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> sourceMergePairs(
    const PartitionAttempt2D& partition,const SolverQualityReport2D& quality,
    const std::vector<bool>& immutableSources);
[[nodiscard]] Polygon2D topologyCellPolygon(const TopologyMesh2D& topology,
                                            std::size_t cellId);
[[nodiscard]] std::vector<std::pair<Polygon2D,Polygon2D>> convexTwoPieceSplits(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    bool wallNormalCuts=false,bool allowCollinear=false);
[[nodiscard]] std::vector<std::array<Polygon2D,3>> convexThreePieceFanSplits(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol);
[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> solverRepartitionPairs(
    const TopologyMesh2D& topology,const SolverQualityReport2D& quality,
    const std::vector<bool>& immutableCells={});
[[nodiscard]] RepartitionBatch2D repartitionPair(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol,SolverTopologyProfile2D* profile,
    const std::vector<bool>& immutableCells);
[[nodiscard]] RepartitionBatch2D agglomerateCellPair(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second,
    const Polygon2D& merged,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,const std::vector<bool>& immutableCells);

[[nodiscard]] bool sortedSetsIntersect(const std::vector<std::size_t>& first,
                                       const std::vector<std::size_t>& second) noexcept;
[[nodiscard]] std::vector<std::vector<std::size_t>> sourceAdjacency(
    const PartitionAttempt2D& partition,std::size_t sourceCount);
[[nodiscard]] std::vector<std::size_t> sourcePairHalo(
    std::size_t first,std::size_t second,
    const std::vector<std::vector<std::size_t>>& adjacency);
[[nodiscard]] std::vector<std::size_t> cellPairHalo(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second);
[[nodiscard]] std::optional<BoundaryRegion2D> patchBoundaryRegion(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& cells,
    const TolerancePolicy& tol);
[[nodiscard]] std::optional<LocalQualityRank2D> localReplacementQualityRank(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& halo,
    const std::vector<std::size_t>& removedCells,
    const std::vector<Polygon2D>& replacementPieces,
    const TolerancePolicy& tol,const SolverQualityPolicy2D& policy = {});

template<class Proposal>
[[nodiscard]] std::vector<Proposal> selectIndependentProposals(
    std::vector<Proposal> proposals) {
    std::sort(proposals.begin(),proposals.end(),[](const Proposal& first,const Proposal& second) {
        if (betterLocalQualityRank(first.rank,second.rank)) return true;
        if (betterLocalQualityRank(second.rank,first.rank)) return false;
        return std::tie(first.first,first.second)<std::tie(second.first,second.second);
    });
    std::vector<std::vector<std::size_t>> halos;
    halos.reserve(proposals.size());
    for (const auto& proposal:proposals) halos.push_back(proposal.halo);
    const auto indices=selectIndependentSolverRepairPatches2D(halos);
    std::vector<Proposal> selected;
    selected.reserve(indices.size());
    for (const auto index:indices) selected.push_back(std::move(proposals[index]));
    return selected;
}

[[nodiscard]] std::vector<Polygon2D> applySourceMergeBatch(
    const std::vector<Polygon2D>& sourcePolygons,
    const std::vector<SourceMergeProposal2D>& selected);
[[nodiscard]] std::vector<bool> applySourceMergeProtection(
    const std::vector<bool>& immutableSources,
    const std::vector<SourceMergeProposal2D>& selected);
[[nodiscard]] std::vector<std::vector<std::size_t>> applySourceMergeLineage(
    const std::vector<std::vector<std::size_t>>& sourceLineages,
    const std::vector<SourceMergeProposal2D>& selected);
[[nodiscard]] std::vector<bool> applySourceMergePreservation(
    const std::vector<bool>& preserveSources,
    const std::vector<SourceMergeProposal2D>& selected);
[[nodiscard]] RepartitionBatch2D applyRepartitionBatch(
    const TopologyMesh2D& topology,
    const std::vector<RepartitionProposal2D>& selected,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol,SolverTopologyProfile2D* profile,
    const std::vector<bool>& immutableCells);
[[nodiscard]] SolverQualityReport2D timedFullQuality(
    const TopologyMesh2D& topology,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool candidate,
    const SolverQualityPolicy2D& policy = {});
[[nodiscard]] SolverLocalRepartitionResult2D repartitionSolverTopologyByQualityImpl(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool useBatch,
    const std::vector<bool>& initialImmutableCells={},bool allowCollinear=false,
    const SolverQualityPolicy2D& policy = {});

} // namespace cartmesh2d::solver_topology_detail
