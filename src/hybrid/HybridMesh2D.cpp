#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include "cartmesh2d/hybrid/TransitionCanonicalization2D.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

using H4ProfileClock = std::chrono::steady_clock;

[[nodiscard]] double h4ProfileSeconds(H4ProfileClock::time_point start) noexcept {
    return std::chrono::duration<double>(H4ProfileClock::now()-start).count();
}

// Process-wide instrumentation counter, single-threaded by construction like
// globalTopologyBuilds in Topology2D.cpp.
std::size_t conformalHybridBuilds=0U;

// Local candidate-search targets retained from the termination repair. They rank
// bounded exact-union/repartition candidates only; they are not an acceptance
// contract and are never reported as a product verdict.
constexpr double minimumRepairFaceOverLocalH=0.01;
constexpr double minimumRepairFaceWeight=0.10;
constexpr double minimumRepairVolumeRatio=0.05;


[[nodiscard]] double polygonBoundsArea(const Polygon2D& polygon) noexcept {
    const auto bounds = polygon.bounds();
    return (bounds.max.x - bounds.min.x) * (bounds.max.y - bounds.min.y);
}

[[nodiscard]] double segmentLength(const Point2D& a, const Point2D& b) noexcept {
    return std::sqrt(squaredNorm(b - a));
}

[[nodiscard]] Polygon2D topologyCellPolygon(
    const TopologyMesh2D& topology,std::size_t cellId) {
    Polygon2D polygon;
    if (cellId>=topology.cells.size()) return polygon;
    polygon.vertices.reserve(topology.cells[cellId].vertices.size());
    for (const auto vertexId:topology.cells[cellId].vertices) {
        if (vertexId>=topology.vertices.size()) return {};
        polygon.vertices.push_back(topology.vertices[vertexId].point);
    }
    return polygon;
}

[[nodiscard]] bool axisAlignedCartesianPolygon(const Polygon2D& polygon) noexcept {
    if (polygon.vertices.size()!=4U) return false;
    const auto bounds=polygon.bounds();
    return std::all_of(polygon.vertices.begin(),polygon.vertices.end(),
                       [&](const Point2D& point) {
        return (point.x==bounds.min.x || point.x==bounds.max.x) &&
               (point.y==bounds.min.y || point.y==bounds.max.y);
    });
}

[[nodiscard]] Polygon2D removeConstructionCollinearVertices(
    Polygon2D polygon,const TolerancePolicy& tol) {
    bool changed=true;
    while (changed && polygon.vertices.size()>3U) {
        changed=false;
        for (std::size_t i=0;i<polygon.vertices.size();++i) {
            const auto n=polygon.vertices.size();
            const auto& previous=polygon.vertices[(i+n-1U)%n];
            const auto& current=polygon.vertices[i];
            const auto& next=polygon.vertices[(i+1U)%n];
            const Vector2D incoming=current-previous;
            const Vector2D outgoing=next-current;
            const double scale=std::sqrt(squaredNorm(incoming)*squaredNorm(outgoing));
            if (!(scale>0.0) || dot(incoming,outgoing)<0.0 ||
                std::abs(cross(incoming,outgoing))>
                    64.0*std::numeric_limits<double>::epsilon()*scale) continue;
            auto candidate=polygon;
            candidate.vertices.erase(candidate.vertices.begin()+
                static_cast<std::ptrdiff_t>(i));
            const double areaScale=std::max(polygon.area(),candidate.area());
            if (std::abs(candidate.area()-polygon.area())>
                tol.absolute*tol.absolute+tol.relative*areaScale) continue;
            polygon=std::move(candidate);
            changed=true;
            break;
        }
    }
    return polygon;
}


[[nodiscard]] bool convexPolygon(const Polygon2D& polygon,
                                 const TolerancePolicy& tol) noexcept {
    if (polygon.vertices.size()<3U || !(polygon.signedArea()>0.0)) return false;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        const auto& a=polygon.vertices[(i+polygon.vertices.size()-1U)%
                                       polygon.vertices.size()];
        const auto& b=polygon.vertices[i];
        const auto& c=polygon.vertices[(i+1U)%polygon.vertices.size()];
        const double scale=std::sqrt(squaredNorm(b-a)*squaredNorm(c-b));
        if (!(scale>0.0) || cross(b-a,c-b)<=tol.areaScale(std::sqrt(scale))) return false;
    }
    return true;
}

[[nodiscard]] bool pointInOrOnTriangle(const Point2D& point,
                                       const Point2D& a,const Point2D& b,
                                       const Point2D& c,
                                       const TolerancePolicy& tol) noexcept {
    const double scale=std::max({segmentLength(a,b),segmentLength(b,c),
                                 segmentLength(c,a)});
    const double epsilon=tol.areaScale(scale);
    return cross(b-a,point-a)>=-epsilon &&
           cross(c-b,point-b)>=-epsilon &&
           cross(a-c,point-c)>=-epsilon;
}

[[nodiscard]] std::optional<std::vector<Polygon2D>> triangulateTermination(
    const Polygon2D& polygon,const TolerancePolicy& tol) {
    if (convexPolygon(polygon,tol)) return std::vector<Polygon2D>{polygon};
    if (polygon.vertices.size()<3U || !(polygon.signedArea()>0.0)) {
        return std::nullopt;
    }
    std::vector<std::size_t> remaining(polygon.vertices.size());
    std::iota(remaining.begin(),remaining.end(),0U);
    std::vector<Polygon2D> pieces;
    while (remaining.size()>3U) {
        std::optional<std::size_t> selected;
        double selectedArea=-1.0;
        for (std::size_t position=0;position<remaining.size();++position) {
            const auto previous=remaining[(position+remaining.size()-1U)%remaining.size()];
            const auto current=remaining[position];
            const auto next=remaining[(position+1U)%remaining.size()];
            const auto& a=polygon.vertices[previous];
            const auto& b=polygon.vertices[current];
            const auto& c=polygon.vertices[next];
            const double turn=cross(b-a,c-b);
            const double edgeProduct=segmentLength(a,b)*segmentLength(b,c);
            if (!(turn>tol.areaScale(std::sqrt(edgeProduct)))) continue;
            bool contains=false;
            for (const auto vertex:remaining) {
                if (vertex==previous || vertex==current || vertex==next) continue;
                if (pointInOrOnTriangle(polygon.vertices[vertex],a,b,c,tol)) {
                    contains=true;
                    break;
                }
            }
            const double area=Polygon2D{{a,b,c}}.area();
            if (!contains && area>selectedArea) {
                selected=position;
                selectedArea=area;
            }
        }
        if (!selected) return std::nullopt;
        const auto position=*selected;
        const auto previous=remaining[(position+remaining.size()-1U)%remaining.size()];
        const auto current=remaining[position];
        const auto next=remaining[(position+1U)%remaining.size()];
        pieces.push_back(Polygon2D{{polygon.vertices[previous],
                                    polygon.vertices[current],
                                    polygon.vertices[next]}});
        remaining.erase(remaining.begin()+static_cast<std::ptrdiff_t>(position));
    }
    pieces.push_back(Polygon2D{{polygon.vertices[remaining[0]],
                                polygon.vertices[remaining[1]],
                                polygon.vertices[remaining[2]]}});
    for (const auto& piece:pieces) {
        if (!convexPolygon(piece,tol)) return std::nullopt;
    }
    const double area=std::accumulate(
        pieces.begin(),pieces.end(),0.0,
        [](double sum,const Polygon2D& piece) { return sum+piece.area(); });
    const double tolerance=tol.absolute*tol.absolute+
                           tol.relative*polygon.area();
    if (std::abs(area-polygon.area())>tolerance) return std::nullopt;
    return pieces;
}

[[nodiscard]] bool pointOnRegionBoundary(const Point2D& point,
                                         const BoundaryRegion2D& region,
                                         const TolerancePolicy& tol) noexcept {
    for (const auto& loop : region.loops()) {
        const auto& vertices = loop.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (pointOnSegment(point,
                               {vertices[i], vertices[(i + 1U) % vertices.size()]},
                               tol)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool edgeOnRegionBoundary(const Point2D& a, const Point2D& b,
                                        const BoundaryRegion2D& region,
                                        const TolerancePolicy& tol) noexcept {
    const Point2D midpoint{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
    return pointOnRegionBoundary(a, region, tol) &&
           pointOnRegionBoundary(b, region, tol) &&
           pointOnRegionBoundary(midpoint, region, tol);
}

[[nodiscard]] bool policyValid(const HybridMeshPolicy2D& policy) noexcept {
    return std::isfinite(policy.remainderSmallCellAreaFraction) &&
           policy.remainderSmallCellAreaFraction>0.0 &&
           policy.remainderSmallCellAreaFraction<1.0 &&
           std::isfinite(policy.tolerance.absolute) &&
           std::isfinite(policy.tolerance.relative) &&
           policy.tolerance.absolute >= 0.0 &&
           policy.tolerance.relative >= 0.0 &&
           std::isfinite(policy.areaToleranceMultiplier) &&
           policy.areaToleranceMultiplier >= 1.0 &&
           std::isfinite(policy.interfaceToleranceMultiplier) &&
           policy.interfaceToleranceMultiplier >= 1.0 &&
           std::isfinite(policy.transitionCellWidthMultiplier) &&
           policy.transitionCellWidthMultiplier > 0.0;
}

[[nodiscard]] HybridMeshBuildResult2D failed(
    HybridMeshFailureReason2D reason, std::string message,
    std::optional<std::size_t> stripId = std::nullopt,
    std::optional<std::size_t> leafId = std::nullopt,
    std::optional<std::size_t> cellId = std::nullopt,
    std::optional<std::size_t> edgeId = std::nullopt) {
    HybridMeshBuildResult2D result;
    result.failure.reason = reason;
    result.failure.message = std::move(message);
    result.failure.stripId = stripId;
    result.failure.leafId = leafId;
    result.failure.cellId = cellId;
    result.failure.edgeId = edgeId;
    return result;
}

[[nodiscard]] bool sameTopology(const TopologyMesh2D& lhs,
                                const TopologyMesh2D& rhs) noexcept {
    if (lhs.vertices.size() != rhs.vertices.size() ||
        lhs.edges.size() != rhs.edges.size() ||
        lhs.cells.size() != rhs.cells.size()) return false;
    for (std::size_t i = 0; i < lhs.vertices.size(); ++i) {
        if (lhs.vertices[i].id != rhs.vertices[i].id ||
            lhs.vertices[i].point.x != rhs.vertices[i].point.x ||
            lhs.vertices[i].point.y != rhs.vertices[i].point.y) return false;
    }
    for (std::size_t i = 0; i < lhs.edges.size(); ++i) {
        const auto& a = lhs.edges[i];
        const auto& b = rhs.edges[i];
        if (a.id != b.id || a.v0 != b.v0 || a.v1 != b.v1 ||
            a.owner != b.owner || a.neighbour != b.neighbour ||
            a.patch != b.patch) return false;
    }
    for (std::size_t i = 0; i < lhs.cells.size(); ++i) {
        const auto& a = lhs.cells[i];
        const auto& b = rhs.cells[i];
        if (a.id != b.id || a.sourceId != b.sourceId ||
            a.sourceKey != b.sourceKey || a.geometryArea != b.geometryArea ||
            a.sourceLineage != b.sourceLineage ||
            a.refinementLineageKeys != b.refinementLineageKeys ||
            a.vertices != b.vertices || a.edges != b.edges) return false;
    }
    return true;
}

[[nodiscard]] CutCell2D topologyAdapter(const HybridSourceCell2D& source,
                                        const TolerancePolicy& tol) {
    CutCell2D adapter;
    adapter.sourceId = source.id;
    // This is an ordering key only. Quadtree provenance remains explicit in
    // HybridSourceCell2D::quadtreeSourceKey and is never decoded for layers.
    adapter.sourceKey = source.id;
    adapter.sourceLineage = {source.id};
    adapter.refinementLineageKeys = source.refinementLineageKeys;
    adapter.backgroundBounds = source.polygon.bounds();
    adapter.kind = source.kind == HybridCellKind2D::RemainderCartesian
        ? CutCellKind::Full : CutCellKind::Cut;
    adapter.fluidPolygon = source.polygon;
    adapter.area = source.area;
    const double boundsArea = polygonBoundsArea(source.polygon);
    adapter.areaFraction = boundsArea > 0.0
        ? std::clamp(source.area / boundsArea, 0.0, 1.0) : 0.0;
    adapter.centroid = source.polygon.centroid(tol);
    adapter.embeddedBoundary = source.embeddedBoundary;
    return adapter;
}

[[nodiscard]] HybridInterfaceAudit2D auditInterface(
    const TopologyMesh2D& topology,const BoundaryRegion2D& outerRegion,
    const BoundaryRegion2D& originalWalls,
    const std::vector<bool>& layerCells,const TolerancePolicy& tol) {
    HybridInterfaceAudit2D audit;
    for (const auto& loop:outerRegion.loops()) {
        const auto& vertices=loop.vertices();
        for (std::size_t edge=0;edge<vertices.size();++edge) {
            const auto& a=vertices[edge];
            const auto& b=vertices[(edge+1U)%vertices.size()];
            if (!edgeOnRegionBoundary(a,b,originalWalls,tol)) {
                audit.expectedLength+=segmentLength(a,b);
            }
        }
    }
    std::map<std::size_t,std::size_t> degree;
    for (const auto& edge:topology.edges) {
        const auto& a=topology.vertices[edge.v0].point;
        const auto& b=topology.vertices[edge.v1].point;
        if (!edgeOnRegionBoundary(a,b,outerRegion,tol)) continue;
        if (edgeOnRegionBoundary(a,b,originalWalls,tol)) continue;
        ++audit.interfaceEdgeCount;
        audit.actualLength+=segmentLength(a,b);
        ++degree[edge.v0];
        ++degree[edge.v1];
        if (!edge.neighbour) {
            ++audit.singleOwnerInterfaceEdges;
            continue;
        }
        const bool ownerLayer=edge.owner<layerCells.size() && layerCells[edge.owner];
        const bool neighbourLayer=*edge.neighbour<layerCells.size() &&
                                  layerCells[*edge.neighbour];
        if (ownerLayer==neighbourLayer) ++audit.wrongCellPairInterfaceEdges;
    }
    audit.interfaceVertexCount=degree.size();
    for (const auto& [vertex,value]:degree) {
        (void)vertex;
        const bool wallEndpoint=value==1U && pointOnRegionBoundary(
            topology.vertices[vertex].point,originalWalls,tol);
        if (value!=2U && !wallEndpoint) ++audit.nonTwoValentInterfaceVertices;
    }
    audit.lengthError=audit.actualLength-audit.expectedLength;
    return audit;
}

struct SolverSourceMetadata2D {
    bool resolved=false;
    HybridCellKind2D kind=HybridCellKind2D::RemainderCartesian;
    double localBackgroundH=0.0;
};

[[nodiscard]] std::vector<SolverSourceMetadata2D> solverSourceMetadata(
    const TopologyMesh2D& topology,const std::vector<HybridSourceCell2D>& sources,
    std::vector<std::optional<std::size_t>>& sourceForCell,
    SourceLineageAudit2D& lineageAudit,bool verifyOracle,
    const TolerancePolicy& tol) {
    std::vector<SolverSourceMetadata2D> metadata(topology.cells.size());
    sourceForCell.assign(topology.cells.size(),std::nullopt);
    lineageAudit={};
    lineageAudit.solverCellCount=topology.cells.size();
    lineageAudit.oracleVerified=verifyOracle;
    for (const auto& cell:topology.cells) {
        Polygon2D polygon;
        for (const auto vertex:cell.vertices) {
            if (vertex<topology.vertices.size()) {
                polygon.vertices.push_back(topology.vertices[vertex].point);
            }
        }
        const auto centroid=polygon.centroid(tol);
        if (!centroid) continue;
        const auto locate=[&](const std::vector<std::size_t>& candidates,
                              std::size_t& checkCount) {
            std::optional<std::size_t> found;
            std::optional<std::size_t> boundaryCandidate;
            for (const auto sourceId:candidates) {
                ++checkCount;
                if (sourceId>=sources.size()) continue;
                const auto& source=sources[sourceId];
                if (!source.polygon.bounds().contains(*centroid,tol)) continue;
                const auto state=classifyPointInPolygon(*centroid,source.polygon,tol);
                if (state==PointInPolygon::Inside) {
                    found=sourceId;
                    break;
                }
                if (state==PointInPolygon::Boundary && !boundaryCandidate) {
                    boundaryCandidate=sourceId;
                }
            }
            return found?found:boundaryCandidate;
        };
        auto lineage=cell.sourceLineage;
        if (lineage.empty()) lineage.push_back(cell.sourceId);
        std::sort(lineage.begin(),lineage.end());
        lineage.erase(std::unique(lineage.begin(),lineage.end()),lineage.end());
        const auto indexed=locate(lineage,lineageAudit.lineageCandidateChecks);
        if (verifyOracle) {
            std::vector<std::size_t> allSources(sources.size());
            std::iota(allSources.begin(),allSources.end(),0U);
            const auto oracle=locate(allSources,lineageAudit.oracleCandidateChecks);
            if (indexed!=oracle) {
                ++lineageAudit.mismatchedCells;
            }
        }
        sourceForCell[cell.id]=indexed;
        if (!sourceForCell[cell.id]) continue;
        const auto& source=sources[*sourceForCell[cell.id]];
        metadata[cell.id].resolved=true;
        metadata[cell.id].kind=source.kind;
        metadata[cell.id].localBackgroundH=source.localBackgroundH;
    }
    return metadata;
}

} // namespace

std::optional<std::vector<Polygon2D>> detail::partitionTerminationPolygon2D(
    const Polygon2D& polygon,const TolerancePolicy& tol) {
    if (!BoundaryLoop(polygon.vertices).diagnose(tol).valid()) return std::nullopt;
    return triangulateTermination(polygon,tol);
}

HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    const HybridMeshPolicy2D& policy) {
    ++conformalHybridBuilds;
    const auto buildStartGlobalTopologies=globalTopologyBuildCount2D();
    const auto buildStartGlobalTopologyCells=globalTopologyBuildInputCells2D();
    const auto buildStartGlobalTopologySeconds=globalTopologyBuildSeconds2D();
    const auto buildStartFullQuality=solverQualityEvaluationCount2D();
    if (!policyValid(policy) || !domain.valid(policy.tolerance) ||
        !originalWalls.diagnose(policy.tolerance).valid() ||
        !boundaryLayers.success() || boundaryLayers.strips.empty()) {
        return failed(HybridMeshFailureReason2D::InvalidInput,
                      "H4-2 requires valid domain, original walls and successful H4-1 strips");
    }
    if (boundaryLayers.strips.size() != originalWalls.loops().size()) {
        return failed(HybridMeshFailureReason2D::InvalidInput,
                      "H4-1 strip count must match original wall-loop count");
    }
    if (remainderRefinement.minimumLevel > remainderMaxLevel ||
        remainderRefinement.boundaryLevel > remainderMaxLevel) {
        return failed(HybridMeshFailureReason2D::InvalidInput,
                      "remainder refinement level exceeds remainderMaxLevel");
    }

    std::vector<BoundaryLoop> outerLoops;
    outerLoops.reserve(boundaryLayers.strips.size());
    bool localTermination=boundaryLayers.localReductionApplied;
    for (std::size_t stripId = 0; stripId < boundaryLayers.strips.size(); ++stripId) {
        const auto& strip = boundaryLayers.strips[stripId];
        localTermination=localTermination || std::any_of(
            strip.actualLayerCounts.begin(),strip.actualLayerCounts.end(),
            [&](std::size_t count) { return count<strip.parameters.nLayers; });
        if (!strip.wallChain.closed || strip.wallChain.fluidSide != FluidSide2D::Right ||
            strip.wallChain.orientation != (policy.fluidRegion == FluidRegion2D::Interior
                ? WallChainOrientation2D::Clockwise : WallChainOrientation2D::CounterClockwise)) {
            return failed(HybridMeshFailureReason2D::UnsupportedWallSemantics,
                          "closed wall strips must march toward the selected fluid region",
                          stripId);
        }
        const auto outer = strip.outerEnvelope();
        BoundaryLoop outerLoop(outer);
        const auto diagnostics = outerLoop.diagnose(policy.tolerance);
        if (!diagnostics.valid() ||
            diagnostics.orientation != (policy.fluidRegion == FluidRegion2D::Interior
                ? LoopOrientation::Clockwise : LoopOrientation::CounterClockwise)) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "H4-1 envelope orientation disagrees with the selected fluid region",
                          stripId);
        }
        const double domainScale = std::max(domain.width(), domain.height());
        const double epsilon = policy.tolerance.scale(domainScale);
        for (const auto& point : outer) {
            if (!domain.bounds.contains(point, policy.tolerance) ||
                point.x <= domain.bounds.min.x + epsilon ||
                point.x >= domain.bounds.max.x - epsilon ||
                point.y <= domain.bounds.min.y + epsilon ||
                point.y >= domain.bounds.max.y - epsilon) {
                return failed(HybridMeshFailureReason2D::OuterEnvelopeOutsideDomain,
                              "outer envelope must lie strictly inside the Cartesian domain",
                              stripId);
            }
        }
        outerLoops.push_back(std::move(outerLoop));
    }

    BoundaryRegion2D outerRegion(outerLoops);
    if (!outerRegion.diagnose(policy.tolerance).valid()) {
        return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                      "outer envelope loops intersect, touch or have invalid nesting");
    }
    const auto outerDepths = outerRegion.nestingDepths(policy.tolerance);
    if (std::any_of(outerDepths.begin(), outerDepths.end(),
                    [](std::size_t depth) { return depth != 0U; })) {
        return failed(HybridMeshFailureReason2D::UnsupportedWallSemantics,
                      "nested outer envelopes are outside the H4-2 fixed-strip scope");
    }

    // Add a graded solver-repairable fan outside the immutable H4-1 strip.
    // The first row preserves every outer-envelope edge as one shared face;
    // later rows double the tangential resolution gradually. This prevents a
    // long layer face from meeting many tiny Cartesian fragments in one jump.
    const double transitionRingThickness=policy.transitionCellWidthMultiplier*std::ldexp(
        std::max(domain.width(),domain.height()),
        -static_cast<int>(remainderRefinement.boundaryLevel));
    std::vector<BoundaryLoop> remainderBoundaryLoops;
    std::vector<Polygon2D> transitionPolygons;
    std::size_t terminationBufferRows=0U;
    double terminationBufferGrowthRatio=0.0;
    double terminationBufferOuterRowThickness=0.0;
    double terminationBufferRowCap=0.0;
    bool terminationBufferRowCapReachable=true;
    remainderBoundaryLoops.reserve(boundaryLayers.strips.size());
    if (policy.fluidRegion == FluidRegion2D::Interior || policy.directLayerConnection) {
        // The inward layer front bounds the Cartesian remainder directly.
        // A common edge partition connects both without extending offsets
        // across the throat or adding thin subdivided transition columns.
        remainderBoundaryLoops=outerLoops;
    } else if (localTermination) {
        double lastLayerSpacing=0.0;
        for (const auto& strip:boundaryLayers.strips) {
            const auto& cumulative=strip.parameters.cumulativeNormalDistances;
            if (!cumulative.empty()) {
                lastLayerSpacing=std::max(
                    lastLayerSpacing,cumulative.size()==1U?cumulative.front():
                    cumulative.back()-cumulative[cumulative.size()-2U]);
            }
        }
        if (!(lastLayerSpacing>0.0)) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "local termination requires positive final layer spacing");
        }
        std::vector<WallChain2D> terminationChains;
        terminationChains.reserve(outerLoops.size());
        for (std::size_t loopId=0;loopId<outerLoops.size();++loopId) {
            const auto chain=makeClosedWallChain2D(
                outerLoops[loopId],loopId,"termination_"+std::to_string(loopId),
                policy.fluidRegion == FluidRegion2D::Interior
                    ? WallFluidRegion2D::Interior : WallFluidRegion2D::Exterior);
            if (!chain.success()) {
                return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                              "local termination front could not form a wall chain",loopId);
            }
            terminationChains.push_back(*chain.chain);
        }
        LayerParameters2D terminationParameters;
        terminationParameters.nLayers=std::max<std::size_t>(3U,policy.transitionRingCount);
        terminationParameters.thicknessMode=LayerThicknessMode2D::FirstLayerThickness;
        terminationParameters.thickness=lastLayerSpacing;
        // A non-dyadic ratio prevents the outer front from repeatedly landing
        // almost on nested quadtree lines while retaining monotone grading.
        terminationParameters.growthRatio=policy.terminationGrowthRatio;
        // Q5-2 EXPERIMENT: cap every buffer row at one remainder background cell
        // while holding the total buffer thickness -- and therefore the outer
        // front position -- at the value the ungated march would have produced.
        if (policy.enableTerminationBufferRadialMatching) {
            const double remainderCellSize=std::ldexp(
                std::max(domain.width(),domain.height()),
                -static_cast<int>(remainderRefinement.boundaryLevel));
            const double cap=policy.outerTransitionRadialTargetCells*remainderCellSize;
            const auto rowThickness=[](double total,std::size_t rows,double ratio) {
                if (rows==0U) return 0.0;
                if (std::abs(ratio-1.0)<1.0e-12)
                    return total/static_cast<double>(rows);
                const double power=std::pow(ratio,static_cast<double>(rows-1U));
                return total*power*(ratio-1.0)/
                       (std::pow(ratio,static_cast<double>(rows))-1.0);
            };
            const auto baseline=resolveLayerParameters2D(terminationParameters);
            if (baseline.success() && cap>0.0 &&
                rowThickness(baseline.parameters->totalThickness,
                             terminationParameters.nLayers,
                             terminationParameters.growthRatio)>cap) {
                const double total=baseline.parameters->totalThickness;
                std::size_t rows=terminationParameters.nLayers;
                const std::size_t maximumRows=
                    policy.maximumTerminationBufferRows;
                while (rows<maximumRows && total/static_cast<double>(rows)>cap)
                    ++rows;
                // The bounded row count may run out before the cap is reachable.
                // A uniform march is the thinnest outer row this row count can
                // produce, so if that still exceeds the cap the rule cannot be
                // satisfied here: leave the historical march untouched instead of
                // reporting a matched buffer that still violates its own cap.
                // This is also what makes ratio 1 a valid bisection lower bound,
                // because rowThickness(total,rows,1) is exactly total/rows.
                if (total/static_cast<double>(rows)<=cap) {
                    // Keep as much near-wall grading as the cap allows: a larger
                    // ratio makes the first row thinner at fixed total thickness.
                    double low=1.0,high=terminationParameters.growthRatio;
                    if (rowThickness(total,rows,high)>cap) {
                        for (int step=0;step<60;++step) {
                            const double middle=0.5*(low+high);
                            if (rowThickness(total,rows,middle)>cap) high=middle;
                            else low=middle;
                        }
                    } else {
                        low=high;
                    }
                    terminationParameters.nLayers=rows;
                    terminationParameters.thicknessMode=
                        LayerThicknessMode2D::TotalThickness;
                    terminationParameters.thickness=total;
                    terminationParameters.growthRatio=low;
                } else {
                    terminationBufferRowCapReachable=false;
                }
            }
            terminationBufferRowCap=cap;
        }
        BoundaryLayerPolicy2D terminationPolicy;
        terminationPolicy.tolerance=policy.tolerance;
        terminationPolicy.permitConcaveTerminationMarching=true;
        terminationBufferRows=terminationParameters.nLayers;
        terminationBufferGrowthRatio=terminationParameters.growthRatio;
        if (const auto resolved=resolveLayerParameters2D(terminationParameters);
            resolved.success() &&
            resolved.parameters->cumulativeNormalDistances.size()>=2U) {
            const auto& d=resolved.parameters->cumulativeNormalDistances;
            terminationBufferOuterRowThickness=d.back()-d[d.size()-2U];
        } else if (const auto single=resolveLayerParameters2D(terminationParameters);
                   single.success()) {
            terminationBufferOuterRowThickness=
                single.parameters->cumulativeNormalDistances.empty()
                    ?0.0:single.parameters->cumulativeNormalDistances.front();
        }
        const auto terminationLayers=buildLocallyReducedBoundaryLayerStrips2D(
            terminationChains,terminationParameters,terminationPolicy);
        if (!terminationLayers.success()) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          std::string("graded local termination buffer failed: ")+
                          boundaryLayerFailureReasonName(terminationLayers.failure.reason)+
                          " "+terminationLayers.failure.message);
        }
        for (const auto& strip:terminationLayers.strips) {
            BoundaryLoop loop(strip.outerEnvelope());
            if (!loop.diagnose(policy.tolerance).valid()) {
                return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                              "graded local termination outer loop is invalid");
            }
            remainderBoundaryLoops.push_back(std::move(loop));
            for (const auto& cell:strip.cells) {
                Polygon2D polygon;
                for (const auto id:cell.vertices) {
                    polygon.vertices.push_back(strip.vertices[id].point);
                }
                transitionPolygons.push_back(std::move(polygon));
            }
        }
    } else {
    if (policy.transitionRingCount==0U || !(transitionRingThickness>0.0)) {
        return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                      "H4 transition fan requires positive width and ring count");
    }
    for (std::size_t stripId=0;stripId<boundaryLayers.strips.size();++stripId) {
        const auto& strip=boundaryLayers.strips[stripId];
        const std::size_t outerRing=strip.ringVertexIds.size()-1U;
        const std::size_t previousRing=outerRing-1U;
        const auto& outerIds=strip.ringVertexIds[outerRing];
        const auto& previousIds=strip.ringVertexIds[previousRing];
        const auto& cumulative=strip.parameters.cumulativeNormalDistances;
        const double lastNormalSpacing=cumulative.size()==1U
            ?cumulative.front()
            :cumulative.back()-cumulative[cumulative.size()-2U];
        if (!(lastNormalSpacing>0.0)) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "H4 transition requires positive final layer spacing",stripId);
        }
        std::vector<Point2D> innerLoop;
        innerLoop.reserve(outerIds.size());
        for (const auto vertex:outerIds) innerLoop.push_back(strip.vertices[vertex].point);
        std::size_t innerSubdivision=1U;
        for (std::size_t ring=0;ring<policy.transitionRingCount;++ring) {
            // Q5-1: the fan is graded tangentially but every ring keeps the
            // full ring thickness, so the outermost ring can be several
            // background cells deep while the Cut cells it borders are a
            // fraction of one. Cut that last ring into radial rows. The last
            // row's offset is arithmetically identical to the single-row
            // offset, so the committed outer envelope -- and therefore the
            // whole remainder quadtree and its Cut cells -- is unchanged.
            const bool lastRing=ring+1U==policy.transitionRingCount;
            const std::size_t radialRows=lastRing
                ?std::max<std::size_t>(
                     1U,policy.transitionOuterRingRadialSubdivision)
                :1U;
            const std::size_t ringSubdivision=ring==0U?1U:innerSubdivision*2U;
            for (std::size_t row=0;row<radialRows;++row) {
            // Only the first row of a ring carries that ring's tangential
            // refinement; later radial rows keep it, so every added face is
            // interior to the transition strip.
            const std::size_t outerSubdivision=
                row==0U?ringSubdivision:innerSubdivision;
            std::vector<Point2D> transitionOuter;
            transitionOuter.reserve(outerIds.size()*outerSubdivision);
            const double factor=transitionRingThickness*
                                (static_cast<double>(ring)+
                                 static_cast<double>(row+1U)/
                                 static_cast<double>(radialRows))/
                                lastNormalSpacing;
            for (std::size_t segment=0;segment<outerIds.size();++segment) {
                const std::size_t next=(segment+1U)%outerIds.size();
                const auto& a=strip.vertices[outerIds[segment]].point;
                const auto& b=strip.vertices[outerIds[next]].point;
                const Vector2D da=strip.vertices[outerIds[segment]].point-
                                  strip.vertices[previousIds[segment]].point;
                const Vector2D db=strip.vertices[outerIds[next]].point-
                                  strip.vertices[previousIds[next]].point;
                for (std::size_t sub=0;sub<outerSubdivision;++sub) {
                    const double t=static_cast<double>(sub)/
                                   static_cast<double>(outerSubdivision);
                    const Point2D base{a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t};
                    const Vector2D delta{da.x+(db.x-da.x)*t,
                                         da.y+(db.y-da.y)*t};
                    transitionOuter.push_back(base+delta*factor);
                }
            }
            BoundaryLoop transitionLoop(transitionOuter);
            if (!transitionLoop.diagnose(policy.tolerance).valid()) {
                return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                              "layer-aware transition envelope is invalid",stripId);
            }
            for (const auto& point:transitionOuter) {
                if (!domain.bounds.contains(point,policy.tolerance)) {
                    return failed(HybridMeshFailureReason2D::OuterEnvelopeOutsideDomain,
                                  "layer-aware transition leaves the Cartesian domain",stripId);
                }
            }
            const std::size_t ratio=outerSubdivision/innerSubdivision;
            const std::size_t innerCount=outerIds.size()*innerSubdivision;
            for (std::size_t segment=0;segment<innerCount;++segment) {
                Polygon2D polygon;
                polygon.vertices.push_back(innerLoop[segment]);
                for (std::size_t sub=0;sub<=ratio;++sub) {
                    polygon.vertices.push_back(
                        transitionOuter[(segment*ratio+sub)%transitionOuter.size()]);
                }
                polygon.vertices.push_back(innerLoop[(segment+1U)%innerLoop.size()]);
                BoundaryLoop cellLoop(polygon.vertices);
                if (!(polygon.signedArea()>0.0) ||
                    !cellLoop.diagnose(policy.tolerance).valid()) {
                    return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                                  "layer-aware transition fan cell is invalid",
                                  stripId,std::nullopt,segment);
                }
                transitionPolygons.push_back(std::move(polygon));
            }
            innerLoop=std::move(transitionOuter);
            innerSubdivision=outerSubdivision;
            }
        }
        BoundaryLoop transitionLoop(innerLoop);
        if (!transitionLoop.diagnose(policy.tolerance).valid()) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "layer-aware transition envelope is invalid",stripId);
        }
        remainderBoundaryLoops.push_back(std::move(transitionLoop));
    }
    }
    IntersectionRegistry2D envelopeRegistry({TransitionCanonicalizationPolicy2D{}.minimumFaceFraction});
    std::string canonicalizationError;
    if (!canonicalizeTransitionEnvelope2D(
            remainderBoundaryLoops,transitionPolygons,outerRegion,domain,
            remainderRefinement.boundaryLevel,transitionRingThickness,
            localTermination,envelopeRegistry,policy.tolerance,canonicalizationError)) {
        return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,canonicalizationError);
    }
    BoundaryRegion2D remainderBoundaryRegion(std::move(remainderBoundaryLoops));
    if (!remainderBoundaryRegion.diagnose(policy.tolerance).valid()) {
        return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                      "layer-aware transition envelopes intersect or nest");
    }

    std::optional<Quadtree2D> remainderTree;
    try {
        remainderTree.emplace(domain, remainderMaxLevel, remainderBoundaryRegion,
                              policy.tolerance);
        remainderTree->refine(remainderBoundaryRegion, remainderRefinement,
                              policy.tolerance);
    } catch (const std::exception& exception) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      std::string("remainder refinement failed: ") + exception.what());
    }
    QuadtreeBalanceReport2D balance;
    try {
        balance = remainderTree->enforceTwoToOneBalance(
            remainderBoundaryRegion, policy.tolerance);
    } catch (const std::exception& exception) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      std::string("remainder 2:1 balance failed: ") + exception.what());
    }
    if (balance.violationsAfter != 0U || !remainderTree->deterministicOrderingValid()) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      "remainder quadtree is not deterministically ordered and 2:1 balanced");
    }

    std::vector<CutCell2D> remainderSourceCells;
    std::shared_ptr<IntersectionRegistry2D> constructionRegistry;
    std::vector<ConstructionRecoveryRequest2D> constructionRecoveryRequests;
    std::vector<QuadtreeLocalRefinementReport2D> constructionRecoveryRefinements;
    std::map<std::uint64_t,std::vector<std::uint64_t>> refinementAncestors;
    std::size_t nextSourceId = 0U;
    std::size_t remainderCartesianCount = 0U;
    std::size_t remainderCutCount = 0U;
    std::size_t transitionPolygonCount = 0U;
    double remainderArea = 0.0;
    for (std::size_t recoveryPass=0;;++recoveryPass) {
        if (recoveryPass>remainderMaxLevel+1U) {
            return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                          "typed construction recovery exceeded bounded refinement passes");
        }
        remainderSourceCells.clear();
        remainderSourceCells.reserve(remainderTree->leaves().size());
        nextSourceId=0U;
        remainderCartesianCount=0U;
        remainderCutCount=0U;
        transitionPolygonCount=0U;
        remainderArea=0.0;
        if (policy.sharedIntersectionConstruction) {
            constructionRegistry=std::make_shared<IntersectionRegistry2D>();
            constructionRegistry->configureGrid(domain.bounds,remainderMaxLevel);
            const double h=std::ldexp(std::min(
                domain.bounds.max.x-domain.bounds.min.x,
                domain.bounds.max.y-domain.bounds.min.y),
                -static_cast<int>(remainderMaxLevel));
            // Use H4's established feature classification. Physical wall
            // vertices remain immutable through every recovery pass.
            for (const auto& strip:boundaryLayers.strips) {
                for (std::size_t i=0;i<strip.wallChain.vertices.size();++i) {
                    const auto kind=strip.wallVertexKinds.at(i);
                    const auto feature=kind==WallVertexKind2D::Sharp
                        ?IntersectionFeature2D::WallSharpCorner:
                        kind==WallVertexKind2D::Concave
                        ?IntersectionFeature2D::WallConcaveCorner:
                        IntersectionFeature2D::Smooth;
                    constructionRegistry->registerGridCornerAnchor(
                        strip.wallChain.vertices[i],h,feature);
                }
            }
            for (const auto& loop:remainderBoundaryRegion.loops())
                for (const auto& vertex:loop.vertices())
                    constructionRegistry->registerGridCornerAnchor(
                        vertex,h,IntersectionFeature2D::TransitionVertex);
        } else constructionRegistry.reset();

        std::optional<std::uint64_t> recoveryLeafKey;
        for (const auto& leaf : remainderTree->leaves()) {
            auto components = constructionRegistry
                ?buildCutCellsShared(leaf,remainderBoundaryRegion,*constructionRegistry,
                                    IntersectionSource2D::TransitionEnvelopeCartesian,
                                    policy.fluidRegion,policy.tolerance)
                :buildCutCells(leaf,remainderBoundaryRegion,
                               policy.fluidRegion,policy.tolerance);
            bool requestedRecovery=false;
            for (auto& component : components) {
                if (!component.constructionRecoveryRequests.empty()) {
                    constructionRecoveryRequests.insert(
                        constructionRecoveryRequests.end(),
                        component.constructionRecoveryRequests.begin(),
                        component.constructionRecoveryRequests.end());
                    recoveryLeafKey=leaf.key;
                    requestedRecovery=true;
                    break;
                }
                if (!component.valid() || component.kind == CutCellKind::Unsupported) {
                    const std::string detail = component.issues.empty()
                        ? "unsupported remainder Cut-cell"
                        : component.issues.front().message;
                    return failed(HybridMeshFailureReason2D::RemainderCutCellFailed,
                                  detail, std::nullopt, leaf.id);
                }
                if (component.kind == CutCellKind::Empty) continue;
                component.sourceId = nextSourceId;
                component.sourceKey = leaf.key;
                if (const auto lineage=refinementAncestors.find(leaf.key);
                    lineage!=refinementAncestors.end()) {
                    component.refinementLineageKeys=lineage->second;
                }
                if (!(component.area > 0.0) || !component.centroid) {
                    return failed(HybridMeshFailureReason2D::RemainderCutCellFailed,
                                  "remainder cell has non-positive area or missing centroid",
                                  std::nullopt, leaf.id);
                }
                const auto centroidState = remainderBoundaryRegion.classifyPoint(
                    *component.centroid,policy.tolerance);
                if (centroidState == (policy.fluidRegion == FluidRegion2D::Interior
                        ? PointInPolygon::Outside : PointInPolygon::Inside) &&
                    convexPolygon(component.fluidPolygon,policy.tolerance)) {
                    std::ostringstream detail;
                    detail<<"remainder cell centroid lies outside the selected remainder region at ("
                          <<component.centroid->x<<','<<component.centroid->y
                          <<") polygon=";
                    for (const auto& point:component.fluidPolygon.vertices)
                        detail<<'('<<point.x<<','<<point.y<<')';
                    return failed(HybridMeshFailureReason2D::RegionClassificationConflict,
                                  detail.str(),std::nullopt,leaf.id);
                }
                const auto kind = component.kind == CutCellKind::Full
                    ? HybridCellKind2D::RemainderCartesian
                    : HybridCellKind2D::RemainderCut;
                if (kind == HybridCellKind2D::RemainderCartesian) {
                    ++remainderCartesianCount;
                } else {
                    ++remainderCutCount;
                    if (component.fluidPolygon.vertices.size() != 4U)
                        ++transitionPolygonCount;
                }
                remainderArea += component.area;
                remainderSourceCells.push_back(std::move(component));
                ++nextSourceId;
            }
            if (requestedRecovery) break;
        }
        if (!recoveryLeafKey) break;
        auto refinement=remainderTree->refineLeavesWithClosure(
            {*recoveryLeafKey},remainderBoundaryRegion,policy.tolerance);
        if (!refinement.pass() || refinement.requestedRefinedLeaves!=1U) {
            return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                          "typed construction conflict cannot be locally refined");
        }
        for (const auto& lineage:refinement.lineage) {
            auto ancestors=refinementAncestors[lineage.parentKey];
            ancestors.push_back(lineage.parentKey);
            for (const auto child:lineage.childKeys)
                refinementAncestors[child]=ancestors;
        }
        constructionRecoveryRefinements.push_back(std::move(refinement));
    }

    const auto remainderTopology=buildGlobalTopology(
        remainderSourceCells,domain,remainderBoundaryRegion,policy.tolerance,constructionRegistry);
    if (!remainderTopology.valid()) {
        std::ostringstream detail;
        detail<<"remainder topology failed before H4 stabilization";
        for (std::size_t i=0;i<std::min<std::size_t>(remainderTopology.issues.size(),4U);++i) {
            const auto& issue=remainderTopology.issues[i];
            detail<<" [code="<<static_cast<unsigned>(issue.code)<<" object="<<issue.objectId
                  <<" "<<issue.message<<']';
        }
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                      detail.str());
    }
    SmallCellPolicy2D smallPolicy;
    smallPolicy.areaFractionThreshold=policy.remainderSmallCellAreaFraction;
    auto remainderSmallCells=analyzeSmallCells(
        remainderSourceCells,remainderTopology,smallPolicy,policy.tolerance);
    if (!remainderSmallCells.valid()) {
        return failed(HybridMeshFailureReason2D::RemainderStabilizationFailed,
                      "hybrid remainder small-cell analysis failed");
    }
    auto remainderStabilization=agglomerateSmallCells(
        remainderSourceCells,remainderTopology,remainderSmallCells,
        domain,remainderBoundaryRegion,policy.tolerance);
    if (!remainderStabilization.valid() &&
        !remainderStabilization.issues.empty() &&
        remainderStabilization.issues.front().code ==
            AgglomerationIssueCode2D::SmallTargetUnsupported) {
        // Concave outer-envelope corners can pinch two small remainder cells
        // together before reaching a stable neighbour. Retry once with bounded
        // small-cell chaining while preserving exact polygon union & area audits.
        auto chainedSmallCells = remainderSmallCells;
        chainedSmallCells.policy.allowBoundedSmallChaining = true;
        auto chained = agglomerateSmallCells(
            remainderSourceCells, remainderTopology, chainedSmallCells,
            domain, remainderBoundaryRegion, policy.tolerance);
        if (chained.valid()) {
            remainderSmallCells = std::move(chainedSmallCells);
            remainderStabilization = std::move(chained);
        }
    }
    if (!remainderStabilization.valid()) {
        const std::string detail=remainderStabilization.issues.empty()
            ?"hybrid remainder agglomeration failed"
            :remainderStabilization.issues.front().message;
        return failed(HybridMeshFailureReason2D::RemainderStabilizationFailed,detail);
    }

    // Q4-1 construction choice: the dominant remaining narrow-gap family is a
    // termination remainder patch produced before the graded strip meets the
    // Cartesian remainder. Select retain-vs-strictly-convex union in the
    // worst face's bounded one-ring before either source is committed to the
    // unified H4 topology. Candidate scoring remains patch-local; only the
    // unique winner reaches the construction-context topology oracle.
    std::size_t q41ConstructionCandidates=0U,q41LocalCandidates=0U;
    std::size_t q41LocalQualityEvaluations=0U,q41AcceptedConstructions=0U;
    std::size_t q41CandidateGlobalTopologyBuilds=0U;
    std::size_t q41CandidateFullGlobalQualityEvaluations=0U;
    std::size_t q41GlobalOracleBuilds=0U,q41MaximumWinnerGlobalOracleBuilds=0U;
    std::size_t q41AuthoritativeFullQualityEvaluations=0U;
    std::size_t q41HardVolumeRatioBefore=0U,q41HardVolumeRatioAfter=0U;
    std::size_t q41HardFaceWeightBefore=0U,q41HardFaceWeightAfter=0U;
    std::size_t q41HardShortFaceBefore=0U,q41HardShortFaceAfter=0U;
    double q41ConstructionSelectionSeconds=0.0;
    bool q41PatchOutsideStableIdsUnchanged=true;
    bool q41LocalDeltaMatchesGlobalOracle=true;
    bool q41LocalWinnerMatchesGlobalAuthority=true;
    bool q41ConstructionBoundReached=false;
    std::vector<AgglomeratedCell2D> remainderConstructionCells=
        remainderStabilization.cells;
    if (policy.enableTerminationConstructionQualitySelection && localTermination) {
        // Include the already-built graded termination strip in the temporary
        // construction context. On the remainder alone its interface is an
        // artificial boundary and legacy boundary-skewness is not meaningful;
        // in this context that same interface has its real two-cell incidence.
        const std::size_t remainderConstructionSourceCount=
            remainderStabilization.cells.size();
        std::vector<CutCell2D> constructionContextSources;
        constructionContextSources.reserve(
            remainderConstructionSourceCount+transitionPolygons.size());
        for (const auto& cell:remainderStabilization.cells) {
            CutCell2D source;
            source.sourceId=constructionContextSources.size();
            source.sourceKey=source.sourceId;
            source.sourceLineage={source.sourceId};
            source.backgroundBounds=cell.polygon.bounds();
            source.kind=CutCellKind::Cut;
            source.fluidPolygon=cell.polygon;
            source.area=cell.area;
            source.centroid=cell.centroid;
            constructionContextSources.push_back(std::move(source));
        }
        for (const auto& polygon:transitionPolygons) {
            CutCell2D source;
            source.sourceId=constructionContextSources.size();
            source.sourceKey=source.sourceId;
            source.sourceLineage={source.sourceId};
            source.backgroundBounds=polygon.bounds();
            source.kind=CutCellKind::Cut;
            source.fluidPolygon=polygon;
            source.area=polygon.area();
            source.centroid=polygon.centroid(policy.tolerance);
            constructionContextSources.push_back(std::move(source));
        }
        auto constructionTopology=buildGlobalTopology(
            constructionContextSources,domain,outerRegion,policy.tolerance,
            constructionRegistry);
        if (!constructionTopology.valid()) {
            return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                          "Q4-1 construction context topology is invalid");
        }
        const auto characterize=[&](const TopologyMesh2D& candidate,
                                    std::vector<bool>& termination,
                                    std::vector<bool>& cartesian,
                                    std::vector<double>& localH,
                                    std::vector<bool>& rated,
                                    std::vector<bool>& immutable) {
            termination.clear();
            cartesian.clear();
            localH.clear();
            rated.assign(candidate.cells.size(),true);
            immutable.clear();
            termination.reserve(candidate.cells.size());
            cartesian.reserve(candidate.cells.size());
            localH.reserve(candidate.cells.size());
            immutable.reserve(candidate.cells.size());
            for (const auto& cell:candidate.cells) {
                const auto polygon=topologyCellPolygon(candidate,cell.id);
                auto groups=cell.sourceLineage;
                if (groups.empty()) groups.push_back(cell.sourceId);
                const bool transitionSource=std::any_of(
                    groups.begin(),groups.end(),[&](std::size_t groupId) {
                        return groupId>=remainderConstructionSourceCount;
                    });
                bool onTermination=false;
                for (std::size_t edge=0;edge<polygon.vertices.size();++edge) {
                    if (edgeOnRegionBoundary(
                            polygon.vertices[edge],
                            polygon.vertices[(edge+1U)%polygon.vertices.size()],
                            remainderBoundaryRegion,policy.tolerance)) {
                        onTermination=true;
                        break;
                    }
                }
                termination.push_back(onTermination || transitionSource);
                cartesian.push_back(!transitionSource &&
                                    axisAlignedCartesianPolygon(polygon));
                immutable.push_back(transitionSource);
                double h=std::numeric_limits<double>::infinity();
                for (const auto groupId:groups) {
                    if (groupId>=remainderConstructionSourceCount) {
                        h=std::min(h,transitionRingThickness);
                        continue;
                    }
                    for (const auto sourceId:
                         remainderStabilization.cells[groupId].memberSourceIds) {
                        if (sourceId>=remainderSourceCells.size()) continue;
                        const auto& bounds=remainderSourceCells[sourceId].backgroundBounds;
                        h=std::min(h,std::max(bounds.max.x-bounds.min.x,
                                              bounds.max.y-bounds.min.y));
                    }
                }
                if (!std::isfinite(h)) {
                    const auto bounds=polygon.bounds();
                    h=std::max(bounds.max.x-bounds.min.x,
                               bounds.max.y-bounds.min.y);
                }
                localH.push_back(h);
            }
        };
        constexpr std::size_t maximumQ41Constructions=8U;
        bool converged=false;
        for (std::size_t iteration=0;iteration<maximumQ41Constructions;++iteration) {
            std::vector<bool> termination,cartesian,rated,immutable;
            std::vector<double> localH;
            characterize(constructionTopology,termination,cartesian,localH,rated,
                         immutable);
            auto selection=repairSolverTerminationQuality2D(
                constructionTopology,domain,outerRegion,immutable,
                termination,cartesian,localH,rated,
                minimumRepairFaceOverLocalH,minimumRepairFaceWeight,
                minimumRepairVolumeRatio,
                TerminationQualityCandidateMode2D::ConstructionAgglomeration,
                policy.tolerance);
            if (iteration==0U) {
                q41HardVolumeRatioBefore=selection.hardVolumeRatioCountBefore;
                q41HardFaceWeightBefore=selection.hardFaceWeightCountBefore;
                q41HardShortFaceBefore=selection.hardShortFaceCountBefore;
            }
            q41HardVolumeRatioAfter=selection.hardVolumeRatioCountAfter;
            q41HardFaceWeightAfter=selection.hardFaceWeightCountAfter;
            q41HardShortFaceAfter=selection.hardShortFaceCountAfter;
            q41ConstructionCandidates+=selection.candidateCount;
            q41LocalCandidates+=selection.localCandidateCount;
            q41LocalQualityEvaluations+=selection.localQualityEvaluationCount;
            q41CandidateGlobalTopologyBuilds+=
                selection.candidateGlobalTopologyBuildCount;
            q41CandidateFullGlobalQualityEvaluations+=
                selection.candidateFullGlobalQualityEvaluationCount;
            q41GlobalOracleBuilds+=selection.globalOracleBuildCount;
            q41MaximumWinnerGlobalOracleBuilds=std::max(
                q41MaximumWinnerGlobalOracleBuilds,selection.globalOracleBuildCount);
            q41AuthoritativeFullQualityEvaluations+=
                selection.authoritativeFullQualityEvaluationCount;
            q41ConstructionSelectionSeconds+=selection.repairSeconds;
            if (!selection.issues.empty()) {
                return failed(HybridMeshFailureReason2D::SolverTopologyFailed,
                              std::string("Q4-1 construction selection failed: ")+
                              selection.issues.front());
            }
            if (!selection.applicable || !selection.accepted) {
                converged=true;
                break;
            }
            q41PatchOutsideStableIdsUnchanged=
                q41PatchOutsideStableIdsUnchanged &&
                selection.patchOutsideStableIdsUnchanged;
            q41LocalDeltaMatchesGlobalOracle=
                q41LocalDeltaMatchesGlobalOracle &&
                selection.localDeltaMatchesGlobalOracle;
            q41LocalWinnerMatchesGlobalAuthority=
                q41LocalWinnerMatchesGlobalAuthority &&
                selection.localWinnerMatchesGlobalAuthority;
            constructionTopology=std::move(selection.topology);
            ++q41AcceptedConstructions;
        }
        q41ConstructionBoundReached=!converged &&
            q41AcceptedConstructions==maximumQ41Constructions;

        remainderConstructionCells.clear();
        remainderConstructionCells.reserve(constructionTopology.cells.size());
        for (const auto& cell:constructionTopology.cells) {
            auto groups=cell.sourceLineage;
            if (groups.empty()) groups.push_back(cell.sourceId);
            if (std::any_of(groups.begin(),groups.end(),[&](std::size_t groupId) {
                    return groupId>=remainderConstructionSourceCount;
                })) continue;
            AgglomeratedCell2D converted;
            converted.id=remainderConstructionCells.size();
            converted.polygon=removeConstructionCollinearVertices(
                topologyCellPolygon(constructionTopology,cell.id),policy.tolerance);
            converted.area=converted.polygon.area();
            converted.centroid=converted.polygon.centroid(policy.tolerance);
            for (const auto groupId:groups) {
                if (groupId>=remainderStabilization.cells.size()) continue;
                converted.memberTopologyCellIds.insert(
                    converted.memberTopologyCellIds.end(),
                    remainderStabilization.cells[groupId].memberTopologyCellIds.begin(),
                    remainderStabilization.cells[groupId].memberTopologyCellIds.end());
                converted.memberSourceIds.insert(
                    converted.memberSourceIds.end(),
                    remainderStabilization.cells[groupId].memberSourceIds.begin(),
                    remainderStabilization.cells[groupId].memberSourceIds.end());
            }
            std::sort(converted.memberTopologyCellIds.begin(),
                      converted.memberTopologyCellIds.end());
            converted.memberTopologyCellIds.erase(
                std::unique(converted.memberTopologyCellIds.begin(),
                            converted.memberTopologyCellIds.end()),
                converted.memberTopologyCellIds.end());
            std::sort(converted.memberSourceIds.begin(),converted.memberSourceIds.end());
            converted.memberSourceIds.erase(
                std::unique(converted.memberSourceIds.begin(),
                            converted.memberSourceIds.end()),
                converted.memberSourceIds.end());
            remainderConstructionCells.push_back(std::move(converted));
        }
    }

    std::vector<HybridSourceCell2D> hybridSources;
    std::size_t requestedLayerCells=0U;
    for (const auto& strip:boundaryLayers.strips) requestedLayerCells+=strip.cells.size();
    std::size_t terminationCellCount=0U;
    hybridSources.reserve(remainderConstructionCells.size()+requestedLayerCells);
    for (const auto& cell:remainderConstructionCells) {
        const auto& polygon=cell.polygon;
        bool cartesian=polygon.vertices.size()==4U;
        const auto bounds=polygon.bounds();
        for (const auto& point:polygon.vertices) {
            const bool corner=(point.x==bounds.min.x || point.x==bounds.max.x) &&
                              (point.y==bounds.min.y || point.y==bounds.max.y);
            cartesian=cartesian && corner;
        }
        bool termination=false;
        if (localTermination) {
            for (std::size_t edge=0;edge<polygon.vertices.size();++edge) {
                if (edgeOnRegionBoundary(
                        polygon.vertices[edge],
                        polygon.vertices[(edge+1U)%polygon.vertices.size()],
                        remainderBoundaryRegion,policy.tolerance)) {
                    termination=true;
                    break;
                }
            }
        }
        const auto pieces=termination
            ?detail::partitionTerminationPolygon2D(polygon,policy.tolerance)
            :std::optional<std::vector<Polygon2D>>(
                std::vector<Polygon2D>{polygon});
        if (!pieces) {
            std::ostringstream detail;
            detail<<std::setprecision(17)
                  <<"termination Cut-cell could not be partitioned before topology area="
                  <<polygon.area()<<" polygon=[";
            for (const auto& point:polygon.vertices) detail<<'('<<point.x<<','<<point.y<<')';
            detail<<']';
            return failed(HybridMeshFailureReason2D::LayerConversionFailed,
                          detail.str());
        }
        for (const auto& piece:*pieces) {
            HybridSourceCell2D source;
            source.id=hybridSources.size();
            source.kind=termination?HybridCellKind2D::Termination:
                        (cartesian?HybridCellKind2D::RemainderCartesian
                                  :HybridCellKind2D::RemainderCut);
            if (termination) ++terminationCellCount;
            source.polygon=piece;
            source.area=piece.area();
            source.localBackgroundH=std::numeric_limits<double>::infinity();
            for (const auto memberSourceId:cell.memberSourceIds) {
                if (memberSourceId>=remainderSourceCells.size()) continue;
                const auto& sourceBounds=
                    remainderSourceCells[memberSourceId].backgroundBounds;
                source.localBackgroundH=std::min(
                    source.localBackgroundH,
                    std::max(sourceBounds.max.x-sourceBounds.min.x,
                             sourceBounds.max.y-sourceBounds.min.y));
                source.refinementLineageKeys.insert(
                    source.refinementLineageKeys.end(),
                    remainderSourceCells[memberSourceId].refinementLineageKeys.begin(),
                    remainderSourceCells[memberSourceId].refinementLineageKeys.end());
            }
            std::sort(source.refinementLineageKeys.begin(),
                      source.refinementLineageKeys.end());
            source.refinementLineageKeys.erase(
                std::unique(source.refinementLineageKeys.begin(),
                            source.refinementLineageKeys.end()),
                source.refinementLineageKeys.end());
            if (!std::isfinite(source.localBackgroundH)) {
                const auto pieceBounds=piece.bounds();
                source.localBackgroundH=std::max(
                    pieceBounds.max.x-pieceBounds.min.x,
                    pieceBounds.max.y-pieceBounds.min.y);
            }
            if (!termination && cell.memberSourceIds.size()==1U) {
                const auto sourceId=cell.memberSourceIds.front();
                if (sourceId<remainderSourceCells.size()) {
                    source.quadtreeSourceKey=remainderSourceCells[sourceId].sourceKey;
                }
            }
            hybridSources.push_back(std::move(source));
        }
    }

    double transitionArea=0.0;
    for (const auto& polygon:transitionPolygons) {
        HybridSourceCell2D source;
        source.id=hybridSources.size();
        source.kind=localTermination?HybridCellKind2D::Termination:
                                     HybridCellKind2D::Transition;
        source.polygon=polygon;
        source.area=polygon.area();
        source.localBackgroundH=transitionRingThickness;
        transitionArea+=source.area;
        if (localTermination) ++terminationCellCount;
        hybridSources.push_back(std::move(source));
    }

    double layerArea = 0.0;
    std::size_t layerCellCount = 0U;
    for (std::size_t stripId = 0; stripId < boundaryLayers.strips.size(); ++stripId) {
        const auto& strip = boundaryLayers.strips[stripId];
        for (const auto& layerCell : strip.cells) {
            HybridSourceCell2D source;
            source.id = hybridSources.size();
            source.kind = HybridCellKind2D::BoundaryLayer;
            source.solverImmutable=true;
            source.layerIndex=layerCell.layer;
            source.wallSegment=layerCell.wallSegment;
            source.stripId=stripId;
            source.localBackgroundH=segmentLength(
                strip.wallChain.segments[layerCell.wallSegment].a,
                strip.wallChain.segments[layerCell.wallSegment].b);
            for (const auto vertexId : layerCell.vertices) {
                source.polygon.vertices.push_back(strip.vertices[vertexId].point);
            }
            source.area = source.polygon.signedArea();
            const double boundsArea = polygonBoundsArea(source.polygon);
            const auto centroid=source.polygon.centroid(policy.tolerance);
            if (!(source.area > 0.0) || !centroid || !(boundsArea > 0.0)) {
                return failed(HybridMeshFailureReason2D::LayerConversionFailed,
                              "H4-1 layer quad could not become a positive topology source",
                              stripId, std::nullopt, layerCell.id);
            }
            const auto originalState = originalWalls.classifyPoint(*centroid,
                                                                    policy.tolerance);
            const auto outerState = outerRegion.classifyPoint(*centroid,
                                                               policy.tolerance);
            if (policy.fluidRegion == FluidRegion2D::Interior
                ? (originalState == PointInPolygon::Outside || outerState == PointInPolygon::Inside)
                : (originalState == PointInPolygon::Inside || outerState == PointInPolygon::Outside)) {
                return failed(HybridMeshFailureReason2D::RegionClassificationConflict,
                              "layer cell is not between original wall and outer envelope",
                              stripId, std::nullopt, layerCell.id);
            }
            if (layerCell.layer == 0U) {
                source.embeddedBoundary.push_back(
                    strip.wallChain.segments[layerCell.wallSegment]);
            }
            layerArea += source.area;
            hybridSources.push_back(std::move(source));
            ++layerCellCount;
        }
    }

    std::vector<CutCell2D> topologyAdapters;
    topologyAdapters.reserve(hybridSources.size());
    for (const auto& source:hybridSources) {
        topologyAdapters.push_back(topologyAdapter(source,policy.tolerance));
    }
    TopologyMesh2D topology = buildGlobalTopology(
        topologyAdapters, domain, originalWalls, policy.tolerance,constructionRegistry);
    if (!topology.valid()) {
        std::string detail = topology.issues.empty()
            ? "unified global topology audit failed"
            : topology.issues.front().message;
        if (!topology.issues.empty() &&
            topology.issues.front().code==TopologyIssueCode2D::UnclassifiedBoundaryEdge) {
            const auto& e=topology.edges[topology.issues.front().objectId];
            const auto a=topology.vertices[e.v0].point;
            const auto b=topology.vertices[e.v1].point;
            std::ostringstream diagnostics;
            diagnostics<<std::setprecision(17)<<" edge="<<e.id<<" owner="<<e.owner
                       <<" a=("<<a.x<<','<<a.y<<") b=("<<b.x<<','<<b.y<<')';
            detail+=diagnostics.str();
        }
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed, detail);
    }

    // Rebuild a second time from immutable sources. This is both a transaction
    // boundary check and a direct determinism guard for canonical interface IDs.
    const TopologyMesh2D repeatedTopology = buildGlobalTopology(
        topologyAdapters, domain, originalWalls, policy.tolerance,constructionRegistry);
    if (!repeatedTopology.valid() || !sameTopology(topology, repeatedTopology)) {
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                      "repeated unified topology build changed IDs or connectivity");
    }

    std::vector<HybridCellRecord2D> records;
    records.reserve(topology.cells.size());
    for (const auto& cell : topology.cells) {
        if (cell.sourceId>=hybridSources.size()) {
            return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                          "topology cell lost its hybrid source metadata",
                          std::nullopt, std::nullopt, cell.id);
        }
        const auto& source=hybridSources[cell.sourceId];
        records.push_back({cell.id, cell.sourceId, source.kind,
                           source.layerIndex, source.wallSegment});
    }

    std::vector<bool> layerTopologyCells;
    layerTopologyCells.reserve(records.size());
    for (const auto& record:records) {
        layerTopologyCells.push_back(record.kind==HybridCellKind2D::BoundaryLayer);
    }
    const auto interfaceAudit=auditInterface(
        topology,outerRegion,originalWalls,layerTopologyCells,policy.tolerance);
    const double interfaceTolerance = policy.interfaceToleranceMultiplier *
        policy.tolerance.scale(interfaceAudit.expectedLength);
    if (!interfaceAudit.pass(interfaceTolerance)) {
        return failed(HybridMeshFailureReason2D::NonConformalInterface,
                      "outer-envelope interface is not a closed two-owner common partition");
    }

    const double domainArea = domain.width() * domain.height();
    const double solidArea = originalWalls.area(policy.tolerance);
    const double outerArea = outerRegion.area(policy.tolerance);
    const bool interior = policy.fluidRegion == FluidRegion2D::Interior;
    const double expectedFluidArea = interior ? solidArea : domainArea - solidArea;
    const double remainderBoundaryArea = remainderBoundaryRegion.area(policy.tolerance);
    const double actualFluidArea = std::accumulate(
        hybridSources.begin(), hybridSources.end(), 0.0,
        [](double sum, const HybridSourceCell2D& cell) { return sum + cell.area; });
    const double areaError = actualFluidArea - expectedFluidArea;
    const double areaTolerance = policy.areaToleranceMultiplier *
        (policy.tolerance.absolute * policy.tolerance.absolute +
         policy.tolerance.relative * std::max(1.0, expectedFluidArea));
    if (std::abs(areaError) > areaTolerance ||
        std::abs(layerArea - (interior ? solidArea - outerArea : outerArea - solidArea)) > areaTolerance ||
        std::abs(remainderArea -
                 (interior ? remainderBoundaryArea : domainArea-remainderBoundaryArea))>
            areaTolerance ||
        std::abs(transitionArea-
                 (interior ? outerArea-remainderBoundaryArea : remainderBoundaryArea-outerArea))>
            areaTolerance) {
        return failed(HybridMeshFailureReason2D::AreaConservationFailed,
                      "hybrid layer/remainder areas do not close to the selected fluid region");
    }

    auto meshQuality = evaluateMeshQuality(
        topology,remainderSourceCells,&remainderSmallCells,policy.tolerance);
    if (!meshQuality.valid()) {
        const std::string detail = meshQuality.issues.empty()
            ? "hybrid mesh quality audit failed"
            : meshQuality.issues.front().message;
        return failed(HybridMeshFailureReason2D::QualityFailed, detail);
    }
    SolverTopologyConstraints2D solverConstraints;
    solverConstraints.immutableInputCells.reserve(topology.cells.size());
    solverConstraints.preserveInputCells.reserve(topology.cells.size());
    solverConstraints.inputPolygonOverrides.reserve(topology.cells.size());
    for (const auto& cell:topology.cells) {
        const auto& source=hybridSources[cell.sourceId];
        // Preserve the established solver-repair constraints: fixed layer cells
        // stay immutable and geometric transition sources retain their shape.
        const bool transition=source.kind==HybridCellKind2D::Transition ||
            (source.kind==HybridCellKind2D::RemainderCut &&
             policy.fluidRegion != FluidRegion2D::Interior &&
             !source.quadtreeSourceKey.has_value());
        solverConstraints.immutableInputCells.push_back(source.solverImmutable);
        solverConstraints.preserveInputCells.push_back(transition);
        // Global common-partition vertices are interface bookkeeping, not
        // geometric corners of the solver source. Preserve every hybrid
        // source polygon so repair cannot repartition on artificial collinear
        // points and create overlapping or non-manifold fragments.
        solverConstraints.inputPolygonOverrides.emplace_back(source.polygon);
    }
    auto solverTopologyReport=buildSolverTopology2D(
        topology,domain,originalWalls,solverConstraints,policy.tolerance);
    if (!solverTopologyReport.valid()) {
        const std::string detail=solverTopologyReport.issues.empty()
            ?"hybrid constrained solver topology failed"
            :solverTopologyReport.issues.front();
        return failed(HybridMeshFailureReason2D::SolverTopologyFailed,detail);
    }
    std::size_t r1ShortFaceCandidates=0U;
    std::size_t r1LocalCandidates=0U;
    std::size_t r1LocalQualityEvaluations=0U;
    std::size_t r1AcceptedTransactions=0U;
    std::size_t r1CandidateGlobalTopologyBuilds=0U;
    std::size_t r1CandidateFullGlobalQualityEvaluations=0U;
    std::size_t r1GlobalOracleBuilds=0U;
    std::size_t r1MaximumWinnerGlobalOracleBuilds=0U;
    std::size_t r1AuthoritativeFullQualityEvaluations=0U;
    double r1MinimumFaceOverLocalHBefore=0.0;
    double r1MinimumFaceOverLocalHAfter=0.0;
    double r1RepairSeconds=0.0;
    bool r1PatchOutsideStableIdsUnchanged=true;
    bool r1LocalDeltaMatchesGlobalOracle=true;
    bool r1LocalWinnerMatchesGlobalAuthority=true;
    std::size_t q3TerminationCandidates=0U;
    std::size_t q3LocalCandidates=0U;
    std::size_t q3LocalQualityEvaluations=0U;
    std::size_t q3AcceptedTransactions=0U;
    std::size_t q3CandidateGlobalTopologyBuilds=0U;
    std::size_t q3CandidateFullGlobalQualityEvaluations=0U;
    std::size_t q3GlobalOracleBuilds=0U;
    std::size_t q3MaximumWinnerGlobalOracleBuilds=0U;
    std::size_t q3AuthoritativeFullQualityEvaluations=0U;
    std::size_t q3HardVolumeRatioBefore=0U,q3HardVolumeRatioAfter=0U;
    std::size_t q3HardFaceWeightBefore=0U,q3HardFaceWeightAfter=0U;
    std::size_t q3HardShortFaceBefore=0U,q3HardShortFaceAfter=0U;
    double q3MaximumVolumeRatioSeverityBefore=0.0;
    double q3MaximumVolumeRatioSeverityAfter=0.0;
    double q3MaximumFaceWeightSeverityBefore=0.0;
    double q3MaximumFaceWeightSeverityAfter=0.0;
    double q3RepairSeconds=0.0;
    bool q3PatchOutsideStableIdsUnchanged=true;
    bool q3LocalDeltaMatchesGlobalOracle=true;
    bool q3LocalWinnerMatchesGlobalAuthority=true;
    bool q3TransactionBoundReached=false;
    std::size_t q32RepartitionCandidates=0U,q32LocalCandidates=0U;
    std::size_t q32LocalQualityEvaluations=0U,q32AcceptedTransactions=0U;
    std::size_t q32CandidateGlobalTopologyBuilds=0U;
    std::size_t q32CandidateFullGlobalQualityEvaluations=0U;
    std::size_t q32GlobalOracleBuilds=0U,q32MaximumWinnerGlobalOracleBuilds=0U;
    std::size_t q32AuthoritativeFullQualityEvaluations=0U;
    std::size_t q32HardVolumeRatioBefore=0U,q32HardVolumeRatioAfter=0U;
    std::size_t q32HardFaceWeightBefore=0U,q32HardFaceWeightAfter=0U;
    std::size_t q32HardShortFaceBefore=0U,q32HardShortFaceAfter=0U;
    double q32RepairSeconds=0.0;
    bool q32PatchOutsideStableIdsUnchanged=true;
    bool q32LocalDeltaMatchesGlobalOracle=true;
    bool q32LocalWinnerMatchesGlobalAuthority=true;
    bool q32TransactionBoundReached=false;
    std::size_t q33GroupedCandidates=0U,q33LocalCandidates=0U;
    std::size_t q33LocalQualityEvaluations=0U,q33AcceptedTransactions=0U;
    std::size_t q33CandidateGlobalTopologyBuilds=0U;
    std::size_t q33CandidateFullGlobalQualityEvaluations=0U;
    std::size_t q33GlobalOracleBuilds=0U,q33MaximumWinnerGlobalOracleBuilds=0U;
    std::size_t q33AuthoritativeFullQualityEvaluations=0U;
    std::size_t q33HardVolumeRatioBefore=0U,q33HardVolumeRatioAfter=0U;
    std::size_t q33HardFaceWeightBefore=0U,q33HardFaceWeightAfter=0U;
    std::size_t q33HardShortFaceBefore=0U,q33HardShortFaceAfter=0U;
    double q33RepairSeconds=0.0;
    bool q33PatchOutsideStableIdsUnchanged=true;
    bool q33LocalDeltaMatchesGlobalOracle=true;
    bool q33LocalWinnerMatchesGlobalAuthority=true;
    bool q33TransactionBoundReached=false;
    if (localTermination) {
        bool converged=false;
        const double minimumFaceFraction=minimumRepairFaceOverLocalH;
        for (std::size_t iteration=0;iteration<32U;++iteration) {
            std::vector<std::optional<std::size_t>> repairSources;
            SourceLineageAudit2D repairLineageAudit;
            const auto repairMetadata=solverSourceMetadata(
                solverTopologyReport.topology,hybridSources,repairSources,
                repairLineageAudit,false,policy.tolerance);
            std::vector<double> localH;
            std::vector<bool> rated;
            for (const auto& metadata:repairMetadata) {
                localH.push_back(metadata.localBackgroundH);
                rated.push_back(metadata.resolved &&
                                metadata.kind!=HybridCellKind2D::BoundaryLayer);
            }
            auto repair=repairSolverShortFaces2D(
                solverTopologyReport.topology,domain,originalWalls,
                solverTopologyReport.immutableOutputCells,localH,rated,
                minimumFaceFraction,policy.tolerance);
            if (iteration==0U)
                r1MinimumFaceOverLocalHBefore=repair.minimumFaceOverLocalHBefore;
            r1MinimumFaceOverLocalHAfter=repair.minimumFaceOverLocalHAfter;
            r1ShortFaceCandidates+=repair.candidateCount;
            r1LocalCandidates+=repair.localCandidateCount;
            r1LocalQualityEvaluations+=repair.localQualityEvaluationCount;
            r1CandidateGlobalTopologyBuilds+=repair.candidateGlobalTopologyBuildCount;
            r1CandidateFullGlobalQualityEvaluations+=
                repair.candidateFullGlobalQualityEvaluationCount;
            r1GlobalOracleBuilds+=repair.globalOracleBuildCount;
            r1MaximumWinnerGlobalOracleBuilds=std::max(
                r1MaximumWinnerGlobalOracleBuilds,repair.globalOracleBuildCount);
            r1AuthoritativeFullQualityEvaluations+=
                repair.authoritativeFullQualityEvaluationCount;
            r1RepairSeconds+=repair.repairSeconds;
            // Only meaningful once a winner reached the authoritative check.
            if (repair.globalOracleBuildCount>0U)
                r1LocalWinnerMatchesGlobalAuthority=
                    r1LocalWinnerMatchesGlobalAuthority &&
                    repair.localWinnerMatchesGlobalAuthority;
            if (!repair.valid()) {
                const std::string detail=repair.issues.empty()
                    ?"R1 patch-local short-face repair failed":repair.issues.front();
                auto rejected=failed(HybridMeshFailureReason2D::SolverTopologyFailed,detail);
                rejected.shortFaceFailure=HybridShortFaceFailure2D{
                    std::move(repair),std::move(localH),std::move(rated),minimumFaceFraction};
                return rejected;
            }
            if (!repair.accepted) {
                converged=!repair.applicable || repair.hardFaceCountBefore==0U;
                break;
            }
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            r1PatchOutsideStableIdsUnchanged=
                r1PatchOutsideStableIdsUnchanged && repair.patchOutsideStableIdsUnchanged;
            r1LocalDeltaMatchesGlobalOracle=
                r1LocalDeltaMatchesGlobalOracle && repair.localDeltaMatchesGlobalOracle;
            ++r1AcceptedTransactions;
            ++solverTopologyReport.qualityRepartitionCount;
            ++solverTopologyReport.profile.acceptedTopologyCommitCount;
        }
        if (!converged)
            return failed(HybridMeshFailureReason2D::SolverTopologyFailed,
                          "R1 patch-local short-face transactions did not converge");
        solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
    }
    if (localTermination && (policy.enableTerminationQualityOptimization ||
                             policy.enableTerminationQualityRepartition ||
                             policy.enableTerminationGroupedRepartition)) {
        bool converged=false;
        // Q3-1 is deliberately bounded: at most 32 winner-only transactions
        // in one build, even when additional eligible faces remain.
        constexpr std::size_t maximumQ3Transactions=32U;
        for (std::size_t iteration=0;iteration<maximumQ3Transactions;++iteration) {
            std::vector<std::optional<std::size_t>> repairSources;
            SourceLineageAudit2D repairLineageAudit;
            const auto repairMetadata=solverSourceMetadata(
                solverTopologyReport.topology,hybridSources,repairSources,
                repairLineageAudit,false,policy.tolerance);
            std::vector<double> localH;
            std::vector<bool> rated,termination,cartesian;
            localH.reserve(repairMetadata.size());
            rated.reserve(repairMetadata.size());
            termination.reserve(repairMetadata.size());
            cartesian.reserve(repairMetadata.size());
            for (const auto& metadata:repairMetadata) {
                localH.push_back(metadata.localBackgroundH);
                rated.push_back(metadata.resolved &&
                                metadata.kind!=HybridCellKind2D::BoundaryLayer);
                termination.push_back(metadata.resolved &&
                                      metadata.kind==HybridCellKind2D::Termination);
                cartesian.push_back(metadata.resolved &&
                                    metadata.kind==HybridCellKind2D::RemainderCartesian);
            }
            auto repair=repairSolverTerminationQuality2D(
                solverTopologyReport.topology,domain,originalWalls,
                solverTopologyReport.immutableOutputCells,termination,cartesian,
                localH,rated,
                minimumRepairFaceOverLocalH,minimumRepairFaceWeight,
                minimumRepairVolumeRatio,
                TerminationQualityCandidateMode2D::Agglomeration,policy.tolerance);
            if (iteration==0U) {
                q3HardVolumeRatioBefore=repair.hardVolumeRatioCountBefore;
                q3HardFaceWeightBefore=repair.hardFaceWeightCountBefore;
                q3HardShortFaceBefore=repair.hardShortFaceCountBefore;
                q3MaximumVolumeRatioSeverityBefore=
                    repair.maximumVolumeRatioSeverityBefore;
                q3MaximumFaceWeightSeverityBefore=
                    repair.maximumFaceWeightSeverityBefore;
            }
            q3HardVolumeRatioAfter=repair.hardVolumeRatioCountAfter;
            q3HardFaceWeightAfter=repair.hardFaceWeightCountAfter;
            q3HardShortFaceAfter=repair.hardShortFaceCountAfter;
            q3MaximumVolumeRatioSeverityAfter=
                repair.maximumVolumeRatioSeverityAfter;
            q3MaximumFaceWeightSeverityAfter=
                repair.maximumFaceWeightSeverityAfter;
            q3TerminationCandidates+=repair.candidateCount;
            q3LocalCandidates+=repair.localCandidateCount;
            q3LocalQualityEvaluations+=repair.localQualityEvaluationCount;
            q3CandidateGlobalTopologyBuilds+=
                repair.candidateGlobalTopologyBuildCount;
            q3CandidateFullGlobalQualityEvaluations+=
                repair.candidateFullGlobalQualityEvaluationCount;
            q3GlobalOracleBuilds+=repair.globalOracleBuildCount;
            q3MaximumWinnerGlobalOracleBuilds=std::max(
                q3MaximumWinnerGlobalOracleBuilds,repair.globalOracleBuildCount);
            q3AuthoritativeFullQualityEvaluations+=
                repair.authoritativeFullQualityEvaluationCount;
            q3RepairSeconds+=repair.repairSeconds;
            if (repair.globalOracleBuildCount>0U)
                q3LocalWinnerMatchesGlobalAuthority=
                    q3LocalWinnerMatchesGlobalAuthority &&
                    repair.localWinnerMatchesGlobalAuthority;
            if (!repair.valid()) {
                const std::string detail=repair.issues.empty()
                    ?"Q3 termination quality repair failed":repair.issues.front();
                return failed(HybridMeshFailureReason2D::SolverTopologyFailed,detail);
            }
            if (!repair.accepted) {
                converged=true;
                break;
            }
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            q3PatchOutsideStableIdsUnchanged=
                q3PatchOutsideStableIdsUnchanged &&
                repair.patchOutsideStableIdsUnchanged;
            q3LocalDeltaMatchesGlobalOracle=
                q3LocalDeltaMatchesGlobalOracle && repair.localDeltaMatchesGlobalOracle;
            ++q3AcceptedTransactions;
            ++solverTopologyReport.qualityRepartitionCount;
            ++solverTopologyReport.profile.acceptedTopologyCommitCount;
        }
        // Hitting the explicit bound is a successful partial Q3 pass, not a
        // claim of global convergence. The report exposes accepted/count data.
        q3TransactionBoundReached=!converged &&
            q3AcceptedTransactions==maximumQ3Transactions;
        solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
    }
    if (localTermination && (policy.enableTerminationQualityRepartition ||
                             policy.enableTerminationGroupedRepartition)) {
        bool converged=false;
        // Q3-2 keeps the search local and cell-count neutral: no more than 16
        // two-cell/two-cell winner-only repartition transactions per build.
        constexpr std::size_t maximumQ32Transactions=16U;
        for (std::size_t iteration=0;iteration<maximumQ32Transactions;++iteration) {
            std::vector<std::optional<std::size_t>> repairSources;
            SourceLineageAudit2D repairLineageAudit;
            const auto repairMetadata=solverSourceMetadata(
                solverTopologyReport.topology,hybridSources,repairSources,
                repairLineageAudit,false,policy.tolerance);
            std::vector<double> localH;
            std::vector<bool> rated,termination,cartesian;
            localH.reserve(repairMetadata.size());
            rated.reserve(repairMetadata.size());
            termination.reserve(repairMetadata.size());
            cartesian.reserve(repairMetadata.size());
            for (const auto& metadata:repairMetadata) {
                localH.push_back(metadata.localBackgroundH);
                rated.push_back(metadata.resolved &&
                                metadata.kind!=HybridCellKind2D::BoundaryLayer);
                termination.push_back(metadata.resolved &&
                                      metadata.kind==HybridCellKind2D::Termination);
                cartesian.push_back(metadata.resolved &&
                                    metadata.kind==HybridCellKind2D::RemainderCartesian);
            }
            auto repair=repairSolverTerminationQuality2D(
                solverTopologyReport.topology,domain,originalWalls,
                solverTopologyReport.immutableOutputCells,termination,cartesian,
                localH,rated,minimumRepairFaceOverLocalH,
                minimumRepairFaceWeight,minimumRepairVolumeRatio,
                TerminationQualityCandidateMode2D::Repartition,policy.tolerance);
            if (iteration==0U) {
                q32HardVolumeRatioBefore=repair.hardVolumeRatioCountBefore;
                q32HardFaceWeightBefore=repair.hardFaceWeightCountBefore;
                q32HardShortFaceBefore=repair.hardShortFaceCountBefore;
            }
            q32HardVolumeRatioAfter=repair.hardVolumeRatioCountAfter;
            q32HardFaceWeightAfter=repair.hardFaceWeightCountAfter;
            q32HardShortFaceAfter=repair.hardShortFaceCountAfter;
            q32RepartitionCandidates+=repair.candidateCount;
            q32LocalCandidates+=repair.localCandidateCount;
            q32LocalQualityEvaluations+=repair.localQualityEvaluationCount;
            q32CandidateGlobalTopologyBuilds+=
                repair.candidateGlobalTopologyBuildCount;
            q32CandidateFullGlobalQualityEvaluations+=
                repair.candidateFullGlobalQualityEvaluationCount;
            q32GlobalOracleBuilds+=repair.globalOracleBuildCount;
            q32MaximumWinnerGlobalOracleBuilds=std::max(
                q32MaximumWinnerGlobalOracleBuilds,repair.globalOracleBuildCount);
            q32AuthoritativeFullQualityEvaluations+=
                repair.authoritativeFullQualityEvaluationCount;
            q32RepairSeconds+=repair.repairSeconds;
            if (repair.globalOracleBuildCount>0U)
                q32LocalWinnerMatchesGlobalAuthority=
                    q32LocalWinnerMatchesGlobalAuthority &&
                    repair.localWinnerMatchesGlobalAuthority;
            if (!repair.valid()) {
                const std::string detail=repair.issues.empty()
                    ?"Q3-2 termination repartition failed":repair.issues.front();
                return failed(HybridMeshFailureReason2D::SolverTopologyFailed,detail);
            }
            if (!repair.accepted) {
                converged=true;
                break;
            }
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            q32PatchOutsideStableIdsUnchanged=
                q32PatchOutsideStableIdsUnchanged &&
                repair.patchOutsideStableIdsUnchanged;
            q32LocalDeltaMatchesGlobalOracle=
                q32LocalDeltaMatchesGlobalOracle &&
                repair.localDeltaMatchesGlobalOracle;
            ++q32AcceptedTransactions;
            ++solverTopologyReport.qualityRepartitionCount;
            ++solverTopologyReport.profile.acceptedTopologyCommitCount;
        }
        q32TransactionBoundReached=!converged &&
            q32AcceptedTransactions==maximumQ32Transactions;
        solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
    }
    if (localTermination && policy.enableTerminationGroupedRepartition) {
        bool converged=false;
        // Q3-3 handles only a fixed Cartesian-plus-two-termination template.
        // Three cells are replaced by three cells, with at most 8 commits.
        constexpr std::size_t maximumQ33Transactions=8U;
        for (std::size_t iteration=0;iteration<maximumQ33Transactions;++iteration) {
            std::vector<std::optional<std::size_t>> repairSources;
            SourceLineageAudit2D repairLineageAudit;
            const auto repairMetadata=solverSourceMetadata(
                solverTopologyReport.topology,hybridSources,repairSources,
                repairLineageAudit,false,policy.tolerance);
            std::vector<double> localH;
            std::vector<bool> rated,termination,cartesian;
            for (const auto& metadata:repairMetadata) {
                localH.push_back(metadata.localBackgroundH);
                rated.push_back(metadata.resolved &&
                                metadata.kind!=HybridCellKind2D::BoundaryLayer);
                termination.push_back(metadata.resolved &&
                                      metadata.kind==HybridCellKind2D::Termination);
                cartesian.push_back(metadata.resolved &&
                                    metadata.kind==HybridCellKind2D::RemainderCartesian);
            }
            auto repair=repairSolverTerminationQuality2D(
                solverTopologyReport.topology,domain,originalWalls,
                solverTopologyReport.immutableOutputCells,termination,cartesian,
                localH,rated,minimumRepairFaceOverLocalH,
                minimumRepairFaceWeight,minimumRepairVolumeRatio,
                TerminationQualityCandidateMode2D::GroupedRepartition,
                policy.tolerance);
            if (iteration==0U) {
                q33HardVolumeRatioBefore=repair.hardVolumeRatioCountBefore;
                q33HardFaceWeightBefore=repair.hardFaceWeightCountBefore;
                q33HardShortFaceBefore=repair.hardShortFaceCountBefore;
            }
            q33HardVolumeRatioAfter=repair.hardVolumeRatioCountAfter;
            q33HardFaceWeightAfter=repair.hardFaceWeightCountAfter;
            q33HardShortFaceAfter=repair.hardShortFaceCountAfter;
            q33GroupedCandidates+=repair.candidateCount;
            q33LocalCandidates+=repair.localCandidateCount;
            q33LocalQualityEvaluations+=repair.localQualityEvaluationCount;
            q33CandidateGlobalTopologyBuilds+=
                repair.candidateGlobalTopologyBuildCount;
            q33CandidateFullGlobalQualityEvaluations+=
                repair.candidateFullGlobalQualityEvaluationCount;
            q33GlobalOracleBuilds+=repair.globalOracleBuildCount;
            q33MaximumWinnerGlobalOracleBuilds=std::max(
                q33MaximumWinnerGlobalOracleBuilds,repair.globalOracleBuildCount);
            q33AuthoritativeFullQualityEvaluations+=
                repair.authoritativeFullQualityEvaluationCount;
            q33RepairSeconds+=repair.repairSeconds;
            if (repair.globalOracleBuildCount>0U)
                q33LocalWinnerMatchesGlobalAuthority=
                    q33LocalWinnerMatchesGlobalAuthority &&
                    repair.localWinnerMatchesGlobalAuthority;
            if (!repair.valid()) {
                const std::string detail=repair.issues.empty()
                    ?"Q3-3 grouped termination repartition failed":
                     repair.issues.front();
                return failed(HybridMeshFailureReason2D::SolverTopologyFailed,detail);
            }
            if (!repair.accepted) {
                converged=true;
                break;
            }
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            q33PatchOutsideStableIdsUnchanged=
                q33PatchOutsideStableIdsUnchanged &&
                repair.patchOutsideStableIdsUnchanged;
            q33LocalDeltaMatchesGlobalOracle=
                q33LocalDeltaMatchesGlobalOracle &&
                repair.localDeltaMatchesGlobalOracle;
            ++q33AcceptedTransactions;
            ++solverTopologyReport.qualityRepartitionCount;
            ++solverTopologyReport.profile.acceptedTopologyCommitCount;
        }
        q33TransactionBoundReached=!converged &&
            q33AcceptedTransactions==maximumQ33Transactions;
        solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
    }
    const auto directionalStart=H4ProfileClock::now();
    auto directionalQuality=evaluateDirectionalConnectivity2D(solverTopologyReport.topology);
    const auto directionalMinimumBefore=directionalQuality.minimumMeasured;
    std::size_t directionalCandidates=0U,directionalAccepted=0U;
    if (directionalQuality.issues.empty() && !directionalQuality.failedCells.empty() &&
        evaluateSolverQuality2D(solverTopologyReport.topology,{},policy.tolerance).valid()) {
        auto repair=improveSolverDirectionalConnectivity2D(
            solverTopologyReport.topology,domain,originalWalls,
            solverTopologyReport.immutableOutputCells,policy.tolerance);
        directionalCandidates=repair.candidateCount;
        if (!repair.valid()) {
            return failed(HybridMeshFailureReason2D::SolverTopologyFailed,
                          repair.issues.empty()?"directional repair failed":repair.issues.front());
        }
        if (repair.acceptedCount>0U) {
            directionalAccepted=repair.acceptedCount;
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
            solverTopologyReport.qualityRepartitionCount+=directionalAccepted;
            solverTopologyReport.profile.acceptedTopologyCommitCount+=directionalAccepted;
            directionalQuality=std::move(repair.after);
        }
    }
    const double directionalSeconds=h4ProfileSeconds(directionalStart);

    // OpenFOAM's packaged meshQualityDict uses 65 degrees for non-orthogonality.
    // Try to repair that concrete target on the final solver topology while the
    // fixed layer cells are still locked. Remaining target issues are reported,
    // but do not redefine the existing native Solver acceptance gate.
    SolverQualityPolicy2D targetPolicy;
    targetPolicy.maxNonOrthogonalityDeg=65.0;
    const auto targetBefore=evaluateSolverQuality2D(
        solverTopologyReport.topology,targetPolicy,policy.tolerance);
    std::size_t targetIssueCountBefore=targetBefore.issues.size();
    std::size_t targetIssueCountAfter=targetIssueCountBefore;
    std::size_t targetRepartitionCount=0U;
    bool targetRepairAttempted=!targetBefore.valid();
    bool targetRepairCommitted=false;
    bool targetRepairStructuralFailure=false;
    if (targetRepairAttempted) {
        auto repair=improveSolverForTargetPolicy2D(
            solverTopologyReport.topology,domain,originalWalls,
            solverTopologyReport.immutableOutputCells,targetPolicy,policy.tolerance);
        if (repair.valid()) {
            const auto targetAfter=evaluateSolverQuality2D(
                repair.topology,targetPolicy,policy.tolerance);
            const auto defaultAfter=evaluateSolverQuality2D(
                repair.topology,{},policy.tolerance);
            const auto directionalAfter=evaluateDirectionalConnectivity2D(
                repair.topology);
            targetIssueCountAfter=targetAfter.issues.size();
            if (repair.repartitionCount>0U &&
                targetIssueCountAfter<=targetIssueCountBefore &&
                defaultAfter.valid() && directionalAfter.valid()) {
                targetRepartitionCount=repair.repartitionCount;
                targetRepairCommitted=true;
                solverTopologyReport.topology=std::move(repair.topology);
                solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
                solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
                solverTopologyReport.qualityRepartitionCount+=targetRepartitionCount;
                solverTopologyReport.profile.acceptedTopologyCommitCount+=
                    targetRepartitionCount;
            } else {
                targetIssueCountAfter=targetIssueCountBefore;
            }
        } else {
            targetRepairStructuralFailure=true;
        }
    }
    std::size_t determinantRepartitionCount=0U;
    if (policy.targetExtrusionThickness) {
        auto repair=improveSolverExtrudedDeterminant2D(
            solverTopologyReport.topology,domain,originalWalls,
            solverTopologyReport.immutableOutputCells,
            *policy.targetExtrusionThickness,targetPolicy,policy.tolerance);
        const auto repairedTarget=repair.valid()
            ?evaluateSolverQuality2D(repair.topology,targetPolicy,policy.tolerance)
            :SolverQualityReport2D{};
        if (repair.valid() && repair.repartitionCount>0U &&
            repairedTarget.issues.size()<=targetIssueCountAfter &&
            evaluateSolverQuality2D(repair.topology,{},policy.tolerance).valid() &&
            evaluateDirectionalConnectivity2D(repair.topology).valid()) {
            determinantRepartitionCount=repair.repartitionCount;
            solverTopologyReport.topology=std::move(repair.topology);
            solverTopologyReport.immutableOutputCells=std::move(repair.immutableCells);
            solverTopologyReport.outputCellCount=solverTopologyReport.topology.cells.size();
            solverTopologyReport.qualityRepartitionCount+=determinantRepartitionCount;
            solverTopologyReport.profile.acceptedTopologyCommitCount+=
                determinantRepartitionCount;
        }
    }
    // Both target-solver repair passes may change the final polygon set.
    directionalQuality=evaluateDirectionalConnectivity2D(solverTopologyReport.topology);
    targetIssueCountAfter=evaluateSolverQuality2D(
        solverTopologyReport.topology,targetPolicy,policy.tolerance).issues.size();
    if (std::count(solverTopologyReport.immutableOutputCells.begin(),
        solverTopologyReport.immutableOutputCells.end(),true)!=
        static_cast<std::ptrdiff_t>(layerCellCount)) {
        return failed(HybridMeshFailureReason2D::SolverTopologyFailed,
                      "solver repair changed the number of fixed H4-1 layer cells");
    }
    std::vector<bool> solverLayerCells;
    solverLayerCells.reserve(solverTopologyReport.topology.cells.size());
    for (const auto& cell:solverTopologyReport.topology.cells) {
        Polygon2D polygon;
        for (const auto vertex:cell.vertices) {
            polygon.vertices.push_back(
                solverTopologyReport.topology.vertices[vertex].point);
        }
        const auto centroid=polygon.centroid(policy.tolerance);
        solverLayerCells.push_back(centroid &&
            outerRegion.classifyPoint(*centroid,policy.tolerance)==
                (interior ? PointInPolygon::Outside : PointInPolygon::Inside));
    }
    const auto solverInterfaceAudit=auditInterface(
        solverTopologyReport.topology,outerRegion,originalWalls,solverLayerCells,
        policy.tolerance);
    if (!solverInterfaceAudit.pass(interfaceTolerance)) {
        return failed(HybridMeshFailureReason2D::NonConformalInterface,
                      "solver repair changed the shared outer-envelope interface");
    }
    auto solverQuality=evaluateSolverQuality2D(
        solverTopologyReport.topology,{},policy.tolerance);
    if (!solverQuality.exportable()) {
        std::ostringstream detail;
        detail<<"solver quality remains invalid after constrained repair: issues="
              <<solverQuality.issues.size()<<" max_nonorthogonality="
              <<solverQuality.maxNonOrthogonalityDeg<<" min_face_weight="
              <<solverQuality.minFaceWeight<<" min_volume_ratio="
              <<solverQuality.minVolumeRatio;
        const std::size_t diagnosticCount=std::min<std::size_t>(4U,solverQuality.issues.size());
        for (std::size_t i=0;i<diagnosticCount;++i) {
            const auto& issue=solverQuality.issues[i];
            detail<<" issue["<<i<<"]=(code="<<static_cast<int>(issue.code)
                  <<",cell="<<issue.cellId<<",edge="<<issue.edgeId
                  <<",measured="<<issue.measured<<')';
            if (issue.cellId<solverTopologyReport.topology.cells.size()) {
                detail<<" vertices=";
                for (const auto vertex:solverTopologyReport.topology.cells[issue.cellId].vertices) {
                    const auto& point=solverTopologyReport.topology.vertices[vertex].point;
                    detail<<'('<<point.x<<','<<point.y<<')';
                }
            }
        }
        auto failure=failed(HybridMeshFailureReason2D::SolverQualityFailed,detail.str());
        failure.solverTopology=solverTopologyReport.topology;
        failure.solverQuality=solverQuality;
        failure.solverTopologyReport=std::move(solverTopologyReport);
        return failure;
    }
    std::vector<std::optional<std::size_t>> sourceForSolverCell;
    SourceLineageAudit2D sourceLineageAudit;
    (void)solverSourceMetadata(
        solverTopologyReport.topology,hybridSources,sourceForSolverCell,
        sourceLineageAudit,policy.verifySourceLineageOracle,policy.tolerance);
    if (!sourceLineageAudit.pass()) {
        return failed(HybridMeshFailureReason2D::RegionClassificationConflict,
                      "R1-A solver source lineage disagrees with full-scan oracle");
    }
    auto constructionIncidence=buildEdgeIncidenceStore2D(topology,0U);
    auto solverIncidence=buildEdgeIncidenceStore2D(solverTopologyReport.topology,0U);
    if (!constructionIncidence.valid() || !solverIncidence.valid()) {
        const auto& issues=!constructionIncidence.valid()
            ?constructionIncidence.issues:solverIncidence.issues;
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                      issues.empty()?"edge-incidence audit failed":issues.front().message);
    }

    HybridMeshBuildResult2D result;
    result.status = HybridMeshStatus2D::Success;
    result.strips = boundaryLayers.strips;
    result.outerEnvelopeRegion = std::move(outerRegion);
    result.remainderSourceCells=std::move(remainderSourceCells);
    result.sourceCells = std::move(hybridSources);
    result.topology = std::move(topology);
    result.solverTopology=solverTopologyReport.topology;
    result.constructionIncidence=std::move(constructionIncidence);
    result.solverIncidence=std::move(solverIncidence);
    result.cellRecords = std::move(records);
    result.interfaceAudit = interfaceAudit;
    result.solverInterfaceAudit=solverInterfaceAudit;
    result.meshQuality = std::move(meshQuality);
    result.solverQuality = std::move(solverQuality);
    result.metrics.directionalRepairCandidateCount=directionalCandidates;
    result.metrics.directionalRepairAcceptedCount=directionalAccepted;
    result.metrics.directionalRepairSeconds=directionalSeconds;
    result.metrics.directionalMinimumBefore=directionalMinimumBefore;
    result.metrics.directionalMinimumAfter=directionalQuality.minimumMeasured;
    result.metrics.directionalFailedCellCount=directionalQuality.failedCells.size();
    result.metrics.targetPolicyRepairIssueCountBefore=targetIssueCountBefore;
    result.metrics.targetPolicyRepairIssueCountAfter=targetIssueCountAfter;
    result.metrics.targetPolicyRepairRepartitionCount=targetRepartitionCount;
    result.metrics.targetPolicyRepairAttempted=targetRepairAttempted;
    result.metrics.targetPolicyRepairCommitted=targetRepairCommitted;
    result.metrics.targetPolicyRepairStructuralFailure=targetRepairStructuralFailure;
    result.metrics.extrudedDeterminantRepairRepartitionCount=
        determinantRepartitionCount;
    result.remainderSmallCells=std::move(remainderSmallCells);
    result.remainderStabilization=std::move(remainderStabilization);
    result.solverTopologyReport=std::move(solverTopologyReport);
    result.sourceLineageAudit=sourceLineageAudit;
    result.canonicalizedIntersections=envelopeRegistry.records();
    result.constructionRecoveryRequests=std::move(constructionRecoveryRequests);
    result.constructionRecoveryRefinements=std::move(constructionRecoveryRefinements);
    for (const auto& cell:result.remainderSourceCells) {
        for (auto record:cell.canonicalizedIntersections) {
            record.source=IntersectionSource2D::TransitionEnvelopeCartesian;
            record.sourceId=cell.sourceId;
            result.canonicalizedIntersections.push_back(record);
        }
    }
    for (std::size_t i=0U;i<result.canonicalizedIntersections.size();++i) {
        auto& record=result.canonicalizedIntersections[i];
        record.id=i;
        for (const auto& vertex:result.solverTopology.vertices) {
            if (segmentLength(vertex.point,record.canonicalVertex.point)<=
                policy.tolerance.relative*record.localH) {
                record.solverVertexId=vertex.id;
                break;
            }
        }
    }
    result.balance = balance;
    result.metrics.quadtreeLeafCount = remainderTree->leaves().size();
    result.metrics.remainderCartesianCellCount = remainderCartesianCount;
    result.metrics.remainderCutCellCount = remainderCutCount;
    result.metrics.boundaryLayerCellCount = layerCellCount;
    for (const auto& strip:boundaryLayers.strips) {
        result.metrics.requestedBoundaryLayerCellCount+=
            strip.metrics.requestedColumnCellCount;
        result.metrics.zeroLayerColumnCount+=strip.metrics.zeroLayerColumnCount;
        result.metrics.terminationEdgeCount+=strip.metrics.terminationEdgeCount;
    }
    result.metrics.terminationCellCount=terminationCellCount;
    result.metrics.transitionPolygonCount = transitionPolygonCount+
                                            transitionPolygons.size();
    result.metrics.remainderSmallCellCount=result.remainderSmallCells.smallCellCount;
    result.metrics.remainderAgglomeratedCellCount=
        result.remainderStabilization.mergedSmallCellCount;
    result.metrics.solverCellCount=result.solverTopology.cells.size();
    result.metrics.solverQualityAgglomerations=
        result.solverTopologyReport.qualityAgglomeratedSourceCellCount;
    result.metrics.solverQualityRepartitions=
        result.solverTopologyReport.qualityRepartitionCount;
    result.metrics.r1ShortFaceCandidates=r1ShortFaceCandidates;
    result.metrics.r1LocalCandidates=r1LocalCandidates;
    result.metrics.r1LocalQualityEvaluations=r1LocalQualityEvaluations;
    result.metrics.r1AcceptedTransactions=r1AcceptedTransactions;
    result.metrics.r1CandidateGlobalTopologyBuilds=r1CandidateGlobalTopologyBuilds;
    result.metrics.r1CandidateFullGlobalQualityEvaluations=
        r1CandidateFullGlobalQualityEvaluations;
    result.metrics.r1GlobalOracleBuilds=r1GlobalOracleBuilds;
    result.metrics.r1MaximumWinnerGlobalOracleBuilds=r1MaximumWinnerGlobalOracleBuilds;
    result.metrics.r1AuthoritativeFullQualityEvaluations=
        r1AuthoritativeFullQualityEvaluations;
    result.metrics.r1RepairSeconds=r1RepairSeconds;
    result.metrics.r1MinimumFaceOverLocalHBefore=r1MinimumFaceOverLocalHBefore;
    result.metrics.r1MinimumFaceOverLocalHAfter=r1MinimumFaceOverLocalHAfter;
    result.metrics.r1PatchOutsideStableIdsUnchanged=r1PatchOutsideStableIdsUnchanged;
    result.metrics.r1LocalDeltaMatchesGlobalOracle=r1LocalDeltaMatchesGlobalOracle;
    result.metrics.r1LocalWinnerMatchesGlobalAuthority=
        r1LocalWinnerMatchesGlobalAuthority;
    result.metrics.q3TerminationCandidates=q3TerminationCandidates;
    result.metrics.q3LocalCandidates=q3LocalCandidates;
    result.metrics.q3LocalQualityEvaluations=q3LocalQualityEvaluations;
    result.metrics.q3AcceptedTransactions=q3AcceptedTransactions;
    result.metrics.q3CandidateGlobalTopologyBuilds=q3CandidateGlobalTopologyBuilds;
    result.metrics.q3CandidateFullGlobalQualityEvaluations=
        q3CandidateFullGlobalQualityEvaluations;
    result.metrics.q3GlobalOracleBuilds=q3GlobalOracleBuilds;
    result.metrics.q3MaximumWinnerGlobalOracleBuilds=
        q3MaximumWinnerGlobalOracleBuilds;
    result.metrics.q3AuthoritativeFullQualityEvaluations=
        q3AuthoritativeFullQualityEvaluations;
    result.metrics.q3HardVolumeRatioBefore=q3HardVolumeRatioBefore;
    result.metrics.q3HardVolumeRatioAfter=q3HardVolumeRatioAfter;
    result.metrics.q3HardFaceWeightBefore=q3HardFaceWeightBefore;
    result.metrics.q3HardFaceWeightAfter=q3HardFaceWeightAfter;
    result.metrics.q3HardShortFaceBefore=q3HardShortFaceBefore;
    result.metrics.q3HardShortFaceAfter=q3HardShortFaceAfter;
    result.metrics.q3MaximumVolumeRatioSeverityBefore=
        q3MaximumVolumeRatioSeverityBefore;
    result.metrics.q3MaximumVolumeRatioSeverityAfter=
        q3MaximumVolumeRatioSeverityAfter;
    result.metrics.q3MaximumFaceWeightSeverityBefore=
        q3MaximumFaceWeightSeverityBefore;
    result.metrics.q3MaximumFaceWeightSeverityAfter=
        q3MaximumFaceWeightSeverityAfter;
    result.metrics.q3RepairSeconds=q3RepairSeconds;
    result.metrics.q3PatchOutsideStableIdsUnchanged=
        q3PatchOutsideStableIdsUnchanged;
    result.metrics.q3LocalDeltaMatchesGlobalOracle=
        q3LocalDeltaMatchesGlobalOracle;
    result.metrics.q3LocalWinnerMatchesGlobalAuthority=
        q3LocalWinnerMatchesGlobalAuthority;
    result.metrics.q3TransactionBoundReached=q3TransactionBoundReached;
    result.metrics.q32RepartitionCandidates=q32RepartitionCandidates;
    result.metrics.q32LocalCandidates=q32LocalCandidates;
    result.metrics.q32LocalQualityEvaluations=q32LocalQualityEvaluations;
    result.metrics.q32AcceptedTransactions=q32AcceptedTransactions;
    result.metrics.q32CandidateGlobalTopologyBuilds=
        q32CandidateGlobalTopologyBuilds;
    result.metrics.q32CandidateFullGlobalQualityEvaluations=
        q32CandidateFullGlobalQualityEvaluations;
    result.metrics.q32GlobalOracleBuilds=q32GlobalOracleBuilds;
    result.metrics.q32MaximumWinnerGlobalOracleBuilds=
        q32MaximumWinnerGlobalOracleBuilds;
    result.metrics.q32AuthoritativeFullQualityEvaluations=
        q32AuthoritativeFullQualityEvaluations;
    result.metrics.q32HardVolumeRatioBefore=q32HardVolumeRatioBefore;
    result.metrics.q32HardVolumeRatioAfter=q32HardVolumeRatioAfter;
    result.metrics.q32HardFaceWeightBefore=q32HardFaceWeightBefore;
    result.metrics.q32HardFaceWeightAfter=q32HardFaceWeightAfter;
    result.metrics.q32HardShortFaceBefore=q32HardShortFaceBefore;
    result.metrics.q32HardShortFaceAfter=q32HardShortFaceAfter;
    result.metrics.q32RepairSeconds=q32RepairSeconds;
    result.metrics.q32PatchOutsideStableIdsUnchanged=
        q32PatchOutsideStableIdsUnchanged;
    result.metrics.q32LocalDeltaMatchesGlobalOracle=
        q32LocalDeltaMatchesGlobalOracle;
    result.metrics.q32LocalWinnerMatchesGlobalAuthority=
        q32LocalWinnerMatchesGlobalAuthority;
    result.metrics.q32TransactionBoundReached=q32TransactionBoundReached;
    result.metrics.q33GroupedCandidates=q33GroupedCandidates;
    result.metrics.q33LocalCandidates=q33LocalCandidates;
    result.metrics.q33LocalQualityEvaluations=q33LocalQualityEvaluations;
    result.metrics.q33AcceptedTransactions=q33AcceptedTransactions;
    result.metrics.q33CandidateGlobalTopologyBuilds=
        q33CandidateGlobalTopologyBuilds;
    result.metrics.q33CandidateFullGlobalQualityEvaluations=
        q33CandidateFullGlobalQualityEvaluations;
    result.metrics.q33GlobalOracleBuilds=q33GlobalOracleBuilds;
    result.metrics.q33MaximumWinnerGlobalOracleBuilds=
        q33MaximumWinnerGlobalOracleBuilds;
    result.metrics.q33AuthoritativeFullQualityEvaluations=
        q33AuthoritativeFullQualityEvaluations;
    result.metrics.q33HardVolumeRatioBefore=q33HardVolumeRatioBefore;
    result.metrics.q33HardVolumeRatioAfter=q33HardVolumeRatioAfter;
    result.metrics.q33HardFaceWeightBefore=q33HardFaceWeightBefore;
    result.metrics.q33HardFaceWeightAfter=q33HardFaceWeightAfter;
    result.metrics.q33HardShortFaceBefore=q33HardShortFaceBefore;
    result.metrics.q33HardShortFaceAfter=q33HardShortFaceAfter;
    result.metrics.q33RepairSeconds=q33RepairSeconds;
    result.metrics.q33PatchOutsideStableIdsUnchanged=
        q33PatchOutsideStableIdsUnchanged;
    result.metrics.q33LocalDeltaMatchesGlobalOracle=
        q33LocalDeltaMatchesGlobalOracle;
    result.metrics.q33LocalWinnerMatchesGlobalAuthority=
        q33LocalWinnerMatchesGlobalAuthority;
    result.metrics.q33TransactionBoundReached=q33TransactionBoundReached;
    result.metrics.q41ConstructionCandidates=q41ConstructionCandidates;
    result.metrics.q41LocalCandidates=q41LocalCandidates;
    result.metrics.q41LocalQualityEvaluations=q41LocalQualityEvaluations;
    result.metrics.q41AcceptedConstructions=q41AcceptedConstructions;
    result.metrics.q41CandidateGlobalTopologyBuilds=
        q41CandidateGlobalTopologyBuilds;
    result.metrics.q41CandidateFullGlobalQualityEvaluations=
        q41CandidateFullGlobalQualityEvaluations;
    result.metrics.q41GlobalOracleBuilds=q41GlobalOracleBuilds;
    result.metrics.q41MaximumWinnerGlobalOracleBuilds=
        q41MaximumWinnerGlobalOracleBuilds;
    result.metrics.q41AuthoritativeFullQualityEvaluations=
        q41AuthoritativeFullQualityEvaluations;
    result.metrics.q41HardVolumeRatioBefore=q41HardVolumeRatioBefore;
    result.metrics.q41HardVolumeRatioAfter=q41HardVolumeRatioAfter;
    result.metrics.q41HardFaceWeightBefore=q41HardFaceWeightBefore;
    result.metrics.q41HardFaceWeightAfter=q41HardFaceWeightAfter;
    result.metrics.q41HardShortFaceBefore=q41HardShortFaceBefore;
    result.metrics.q41HardShortFaceAfter=q41HardShortFaceAfter;
    result.metrics.q41ConstructionSelectionSeconds=
        q41ConstructionSelectionSeconds;
    result.metrics.q41PatchOutsideStableIdsUnchanged=
        q41PatchOutsideStableIdsUnchanged;
    result.metrics.q41LocalDeltaMatchesGlobalOracle=
        q41LocalDeltaMatchesGlobalOracle;
    result.metrics.q41LocalWinnerMatchesGlobalAuthority=
        q41LocalWinnerMatchesGlobalAuthority;
    result.metrics.q41ConstructionBoundReached=q41ConstructionBoundReached;
    result.metrics.q52TerminationBufferRows=terminationBufferRows;
    result.metrics.q52TerminationBufferGrowthRatio=terminationBufferGrowthRatio;
    result.metrics.q52TerminationBufferOuterRowThickness=
        terminationBufferOuterRowThickness;
    result.metrics.q52TerminationBufferRowCap=terminationBufferRowCap;
    result.metrics.q52TerminationBufferRowCapReachable=
        terminationBufferRowCapReachable;
    result.metrics.buildGlobalTopologyCalls=
        globalTopologyBuildCount2D()-buildStartGlobalTopologies;
    result.metrics.buildGlobalTopologyInputCells=
        globalTopologyBuildInputCells2D()-buildStartGlobalTopologyCells;
    result.metrics.buildGlobalTopologySeconds=
        globalTopologyBuildSeconds2D()-buildStartGlobalTopologySeconds;
    result.metrics.fullSolverQualityEvaluations=
        solverQualityEvaluationCount2D()-buildStartFullQuality;
    result.metrics.solverRepairSeconds=
        result.solverTopologyReport.profile.initialPartitionSeconds+
        result.solverTopologyReport.profile.sourceRepairSeconds+
        result.solverTopologyReport.profile.finalRepartitionSeconds+
        r1RepairSeconds+q3RepairSeconds+q32RepairSeconds+q33RepairSeconds;
    result.metrics.unifiedVertexCount = result.topology.vertices.size();
    result.metrics.unifiedEdgeCount = result.topology.edges.size();
    result.metrics.unifiedCellCount = result.topology.cells.size();
    result.metrics.solidArea = solidArea;
    result.metrics.outerEnvelopeArea = outerArea;
    result.metrics.layerArea = layerArea;
    result.metrics.remainderArea = remainderArea+transitionArea;
    result.metrics.expectedFluidArea = expectedFluidArea;
    result.metrics.actualFluidArea = actualFluidArea;
    result.metrics.areaError = areaError;
    result.metrics.terminationGrowthRatio=policy.terminationGrowthRatio;
    return result;
}


std::size_t conformalHybridBuildCount2D() noexcept {
    return conformalHybridBuilds;
}


} // namespace cartmesh2d
