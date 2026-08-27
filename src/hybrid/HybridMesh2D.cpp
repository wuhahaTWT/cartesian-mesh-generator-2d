#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <algorithm>
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

[[nodiscard]] double polygonBoundsArea(const Polygon2D& polygon) noexcept {
    const auto bounds = polygon.bounds();
    return (bounds.max.x - bounds.min.x) * (bounds.max.y - bounds.min.y);
}

[[nodiscard]] double segmentLength(const Point2D& a, const Point2D& b) noexcept {
    return std::sqrt(squaredNorm(b - a));
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
        if (!(scale>0.0) || cross(b-a,c-b)<=tol.scale(scale)) return false;
    }
    return true;
}

[[nodiscard]] bool pointInOrOnTriangle(const Point2D& point,
                                       const Point2D& a,const Point2D& b,
                                       const Point2D& c,
                                       const TolerancePolicy& tol) noexcept {
    const double scale=std::max({segmentLength(a,b),segmentLength(b,c),
                                 segmentLength(c,a),1.0});
    const double epsilon=tol.scale(scale)*scale;
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
            if (!(turn>tol.scale(std::max(1.0,std::abs(turn))))) continue;
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
    const double area=std::accumulate(
        pieces.begin(),pieces.end(),0.0,
        [](double sum,const Polygon2D& piece) { return sum+piece.area(); });
    const double tolerance=tol.absolute*tol.absolute+
                           tol.relative*std::max(1.0,polygon.area());
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
    return std::isfinite(policy.tolerance.absolute) &&
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

void setError(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

void writeJsonString(std::ostream& out, const std::string& value) {
    out << '"';
    for (const char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << character; break;
        }
    }
    out << '"';
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

} // namespace

HybridMeshBuildResult2D buildConformalHybridMesh2D(
    const BoundaryLayerBuildResult2D& boundaryLayers,
    const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,
    std::size_t remainderMaxLevel,
    const QuadtreeRefinementPolicy2D& remainderRefinement,
    const HybridMeshPolicy2D& policy) {
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
            strip.wallChain.orientation != WallChainOrientation2D::CounterClockwise) {
            return failed(HybridMeshFailureReason2D::UnsupportedWallSemantics,
                          "H4-2 currently supports closed exterior wall strips only",
                          stripId);
        }
        const auto outer = strip.outerEnvelope();
        BoundaryLoop outerLoop(outer);
        const auto diagnostics = outerLoop.diagnose(policy.tolerance);
        if (!diagnostics.valid() ||
            diagnostics.orientation != LoopOrientation::CounterClockwise) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "H4-1 outer envelope is not a valid counter-clockwise loop",
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
    remainderBoundaryLoops.reserve(boundaryLayers.strips.size());
    if (localTermination) {
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
                outerLoops[loopId],loopId,"termination_"+std::to_string(loopId));
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
        BoundaryLayerPolicy2D terminationPolicy;
        terminationPolicy.tolerance=policy.tolerance;
        terminationPolicy.permitConcaveTerminationMarching=true;
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
            const std::size_t outerSubdivision=ring==0U?1U:innerSubdivision*2U;
            std::vector<Point2D> transitionOuter;
            transitionOuter.reserve(outerIds.size()*outerSubdivision);
            const double factor=transitionRingThickness*static_cast<double>(ring+1U)/
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
        BoundaryLoop transitionLoop(innerLoop);
        if (!transitionLoop.diagnose(policy.tolerance).valid()) {
            return failed(HybridMeshFailureReason2D::InvalidOuterEnvelope,
                          "layer-aware transition envelope is invalid",stripId);
        }
        remainderBoundaryLoops.push_back(std::move(transitionLoop));
    }
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
    remainderSourceCells.reserve(remainderTree->leaves().size());
    std::size_t nextSourceId = 0U;
    std::size_t remainderCartesianCount = 0U;
    std::size_t remainderCutCount = 0U;
    std::size_t transitionPolygonCount = 0U;
    double remainderArea = 0.0;
    for (const auto& leaf : remainderTree->leaves()) {
        auto components = buildCutCells(leaf, remainderBoundaryRegion,
                                        FluidRegion2D::Exterior,
                                        policy.tolerance);
        for (auto& component : components) {
            if (!component.valid() || component.kind == CutCellKind::Unsupported) {
                const std::string detail = component.issues.empty()
                    ? "unsupported remainder Cut-cell"
                    : component.issues.front().message;
                return failed(HybridMeshFailureReason2D::RemainderCutCellFailed,
                              detail, std::nullopt, leaf.id);
            }
            if (component.kind == CutCellKind::Empty) continue;
            if (!(component.area > 0.0) || !component.centroid) {
                return failed(HybridMeshFailureReason2D::RemainderCutCellFailed,
                              "remainder cell has non-positive area or missing centroid",
                              std::nullopt, leaf.id);
            }
            const auto centroidState = remainderBoundaryRegion.classifyPoint(*component.centroid,
                                                                  policy.tolerance);
            // A certified exterior component may be concave around a stepped
            // termination corner and its area centroid may lie inside the
            // excluded envelope. That is not a classification conflict: it is
            // partitioned into convex solver cells below. For a convex source,
            // however, an inside centroid proves the cutter chose the wrong face.
            if (centroidState == PointInPolygon::Inside &&
                convexPolygon(component.fluidPolygon,policy.tolerance)) {
                std::ostringstream detail;
                detail<<"remainder cell centroid lies inside the outer envelope at ("
                      <<component.centroid->x<<','<<component.centroid->y
                      <<") polygon=";
                for (const auto& point:component.fluidPolygon.vertices) {
                    detail<<'('<<point.x<<','<<point.y<<')';
                }
                return failed(HybridMeshFailureReason2D::RegionClassificationConflict,
                              detail.str(),
                              std::nullopt, leaf.id);
            }
            component.sourceId = nextSourceId;
            component.sourceKey = leaf.key;
            const auto kind = component.kind == CutCellKind::Full
                ? HybridCellKind2D::RemainderCartesian
                : HybridCellKind2D::RemainderCut;
            if (kind == HybridCellKind2D::RemainderCartesian) {
                ++remainderCartesianCount;
            } else {
                ++remainderCutCount;
                if (component.fluidPolygon.vertices.size() != 4U) {
                    ++transitionPolygonCount;
                }
            }
            remainderArea += component.area;
            remainderSourceCells.push_back(std::move(component));
            ++nextSourceId;
        }
    }

    const auto remainderTopology=buildGlobalTopology(
        remainderSourceCells,domain,remainderBoundaryRegion,policy.tolerance);
    if (!remainderTopology.valid()) {
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                      "remainder topology failed before H4 stabilization");
    }
    SmallCellPolicy2D smallPolicy;
    smallPolicy.areaFractionThreshold=0.10;
    auto remainderSmallCells=analyzeSmallCells(
        remainderSourceCells,remainderTopology,smallPolicy,policy.tolerance);
    if (!remainderSmallCells.valid()) {
        return failed(HybridMeshFailureReason2D::RemainderStabilizationFailed,
                      "hybrid remainder small-cell analysis failed");
    }
    auto remainderStabilization=agglomerateSmallCells(
        remainderSourceCells,remainderTopology,remainderSmallCells,
        domain,remainderBoundaryRegion,policy.tolerance);
    if (!remainderStabilization.valid()) {
        const std::string detail=remainderStabilization.issues.empty()
            ?"hybrid remainder agglomeration failed"
            :remainderStabilization.issues.front().message;
        return failed(HybridMeshFailureReason2D::RemainderStabilizationFailed,detail);
    }

    std::vector<HybridSourceCell2D> hybridSources;
    std::size_t requestedLayerCells=0U;
    for (const auto& strip:boundaryLayers.strips) requestedLayerCells+=strip.cells.size();
    std::size_t terminationCellCount=0U;
    hybridSources.reserve(remainderStabilization.cells.size()+requestedLayerCells);
    for (const auto& cell:remainderStabilization.cells) {
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
            ?triangulateTermination(polygon,policy.tolerance)
            :std::optional<std::vector<Polygon2D>>(
                std::vector<Polygon2D>{polygon});
        if (!pieces) {
            return failed(HybridMeshFailureReason2D::LayerConversionFailed,
                          "termination Cut-cell could not be partitioned before topology");
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
                                     HybridCellKind2D::RemainderCut;
        source.polygon=polygon;
        source.area=polygon.area();
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
            if (originalState == PointInPolygon::Inside ||
                outerState == PointInPolygon::Outside) {
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
        topologyAdapters, domain, originalWalls, policy.tolerance);
    if (!topology.valid()) {
        const std::string detail = topology.issues.empty()
            ? "unified global topology audit failed"
            : topology.issues.front().message;
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed, detail);
    }

    // Rebuild a second time from immutable sources. This is both a transaction
    // boundary check and a direct determinism guard for canonical interface IDs.
    const TopologyMesh2D repeatedTopology = buildGlobalTopology(
        topologyAdapters, domain, originalWalls, policy.tolerance);
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
    const double expectedFluidArea = domainArea - solidArea;
    const double actualFluidArea = std::accumulate(
        hybridSources.begin(), hybridSources.end(), 0.0,
        [](double sum, const HybridSourceCell2D& cell) { return sum + cell.area; });
    const double areaError = actualFluidArea - expectedFluidArea;
    const double areaTolerance = policy.areaToleranceMultiplier *
        (policy.tolerance.absolute * policy.tolerance.absolute +
         policy.tolerance.relative * std::max(1.0, expectedFluidArea));
    if (std::abs(areaError) > areaTolerance ||
        std::abs(layerArea - (outerArea - solidArea)) > areaTolerance ||
        std::abs(remainderArea -
                 (domainArea-remainderBoundaryRegion.area(policy.tolerance)))>
            areaTolerance ||
        std::abs(transitionArea-
                 (remainderBoundaryRegion.area(policy.tolerance)-outerArea))>
            areaTolerance) {
        return failed(HybridMeshFailureReason2D::AreaConservationFailed,
                      "hybrid layer/remainder areas do not close to domain minus solid");
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
        const bool transition=source.kind==HybridCellKind2D::RemainderCut &&
                              !source.quadtreeSourceKey.has_value();
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
            outerRegion.classifyPoint(*centroid,policy.tolerance)==PointInPolygon::Inside);
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
    if (!solverQuality.valid()) {
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
        }
        auto failure=failed(HybridMeshFailureReason2D::SolverQualityFailed,detail.str());
        failure.solverTopology=solverTopologyReport.topology;
        failure.solverQuality=solverQuality;
        failure.solverTopologyReport=std::move(solverTopologyReport);
        return failure;
    }

    HybridMeshBuildResult2D result;
    result.status = HybridMeshStatus2D::Success;
    result.strips = boundaryLayers.strips;
    result.outerEnvelopeRegion = std::move(outerRegion);
    result.remainderSourceCells=std::move(remainderSourceCells);
    result.sourceCells = std::move(hybridSources);
    result.topology = std::move(topology);
    result.solverTopology=solverTopologyReport.topology;
    result.cellRecords = std::move(records);
    result.interfaceAudit = interfaceAudit;
    result.solverInterfaceAudit=solverInterfaceAudit;
    result.meshQuality = std::move(meshQuality);
    result.solverQuality = std::move(solverQuality);
    result.remainderSmallCells=std::move(remainderSmallCells);
    result.remainderStabilization=std::move(remainderStabilization);
    result.solverTopologyReport=std::move(solverTopologyReport);
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

[[nodiscard]] PureCutCellFallback2D buildPureCutCellFallback2D(
    const Domain2D& domain,const BoundaryRegion2D& originalWalls,
    std::size_t maxLevel,const QuadtreeRefinementPolicy2D& refinement,
    const HybridMeshPolicy2D& policy) {
    PureCutCellFallback2D result;
    if (!domain.valid(policy.tolerance) ||
        !originalWalls.diagnose(policy.tolerance).valid() ||
        refinement.minimumLevel>maxLevel || refinement.boundaryLevel>maxLevel) {
        result.failureMessage="pure Cut-cell fallback requires valid inputs";
        return result;
    }
    std::optional<Quadtree2D> tree;
    try {
        tree.emplace(domain,maxLevel,originalWalls,policy.tolerance);
        tree->refine(originalWalls,refinement,policy.tolerance);
        result.balance=tree->enforceTwoToOneBalance(originalWalls,policy.tolerance);
    } catch (const std::exception& exception) {
        result.failureMessage=std::string("pure Cut-cell fallback refinement failed: ")+
                              exception.what();
        return result;
    }
    if (result.balance.violationsAfter!=0U || !tree->deterministicOrderingValid()) {
        result.failureMessage="pure Cut-cell fallback quadtree is not balanced/deterministic";
        return result;
    }
    result.quadtreeLeafCount=tree->leaves().size();
    std::size_t sourceId=0U;
    for (const auto& leaf:tree->leaves()) {
        auto components=buildCutCells(
            leaf,originalWalls,FluidRegion2D::Exterior,policy.tolerance);
        for (auto& component:components) {
            if (!component.valid() || component.kind==CutCellKind::Unsupported) {
                result.failureMessage=component.issues.empty()
                    ?"pure Cut-cell fallback encountered unsupported geometry"
                    :component.issues.front().message;
                return result;
            }
            component.sourceId=sourceId++;
            component.sourceKey=leaf.key;
            if (component.kind!=CutCellKind::Empty) result.actualFluidArea+=component.area;
            result.sourceCells.push_back(std::move(component));
        }
    }
    result.expectedFluidArea=domain.width()*domain.height()-
                             originalWalls.area(policy.tolerance);
    result.areaError=result.actualFluidArea-result.expectedFluidArea;
    const double areaTolerance=policy.areaToleranceMultiplier*
        (policy.tolerance.absolute*policy.tolerance.absolute+
         policy.tolerance.relative*std::max(1.0,result.expectedFluidArea));
    if (std::abs(result.areaError)>areaTolerance) {
        result.failureMessage="pure Cut-cell fallback failed fluid-area conservation";
        return result;
    }
    const auto sourceTopology=buildGlobalTopology(
        result.sourceCells,domain,originalWalls,policy.tolerance);
    if (!sourceTopology.valid()) {
        result.failureMessage=sourceTopology.issues.empty()
            ?"pure Cut-cell fallback source topology failed"
            :sourceTopology.issues.front().message;
        return result;
    }
    SmallCellPolicy2D smallPolicy;
    smallPolicy.areaFractionThreshold=0.10;
    result.smallCells=analyzeSmallCells(
        result.sourceCells,sourceTopology,smallPolicy,policy.tolerance);
    if (!result.smallCells.valid()) {
        result.failureMessage="pure Cut-cell fallback small-cell analysis failed";
        return result;
    }
    result.stabilization=agglomerateSmallCells(
        result.sourceCells,sourceTopology,result.smallCells,
        domain,originalWalls,policy.tolerance);
    if (!result.stabilization.valid()) {
        result.failureMessage=result.stabilization.issues.empty()
            ?"pure Cut-cell fallback agglomeration failed"
            :result.stabilization.issues.front().message;
        return result;
    }
    result.topology=result.stabilization.topology;
    result.meshQuality=evaluateMeshQuality(
        result.topology,result.sourceCells,&result.smallCells,policy.tolerance);
    if (!result.meshQuality.valid()) {
        result.failureMessage="pure Cut-cell fallback mesh-quality audit failed";
        return result;
    }
    result.solverTopologyReport=buildSolverTopology2D(
        result.topology,domain,originalWalls,policy.tolerance);
    if (!result.solverTopologyReport.valid()) {
        result.failureMessage=result.solverTopologyReport.issues.empty()
            ?"pure Cut-cell fallback solver topology failed"
            :result.solverTopologyReport.issues.front();
        return result;
    }
    result.solverTopology=result.solverTopologyReport.topology;
    result.solverQuality=evaluateSolverQuality2D(
        result.solverTopology,{},policy.tolerance);
    if (!result.solverQuality.valid()) {
        std::ostringstream message;
        message<<"pure Cut-cell fallback solver quality failed: issues="
               <<result.solverQuality.issues.size()
               <<" max_nonorthogonality="
               <<result.solverQuality.maxNonOrthogonalityDeg
               <<" min_face_weight="<<result.solverQuality.minFaceWeight
               <<" min_volume_ratio="<<result.solverQuality.minVolumeRatio;
        result.failureMessage=message.str();
    }
    return result;
}

RobustH4BuildResult2D buildRobustH4Mesh2D(
    const std::vector<WallChain2D>& wallChains,
    const LayerParameters2D& layerParameters,const Domain2D& domain,
    const BoundaryRegion2D& originalWalls,std::size_t maxLevel,
    const QuadtreeRefinementPolicy2D& refinement,
    const BoundaryLayerPolicy2D& layerPolicy,
    const HybridMeshPolicy2D& hybridPolicy) {
    RobustH4BuildResult2D result;
    result.requestedLayerCandidate=buildBoundaryLayerStrips2D(
        wallChains,layerParameters,layerPolicy);
    if (result.requestedLayerCandidate.success()) {
        result.hybridCandidate=buildConformalHybridMesh2D(
            result.requestedLayerCandidate,domain,originalWalls,maxLevel,refinement);
        if (result.hybridCandidate.success()) {
            result.mode=H4MeshMode2D::Hybrid;
            return result;
        }
        result.fallbackStage=H4FallbackStage2D::HybridCandidate;
    } else {
        result.fallbackStage=H4FallbackStage2D::RequestedLayers;
    }

    result.localLayerCandidate=buildLocallyReducedBoundaryLayerStrips2D(
        wallChains,layerParameters,layerPolicy);
    if (result.localLayerCandidate.success()) {
        result.hybridCandidate=buildConformalHybridMesh2D(
            result.localLayerCandidate,domain,originalWalls,maxLevel,refinement);
        if (result.hybridCandidate.success()) {
            result.mode=H4MeshMode2D::Hybrid;
            result.fallbackStage=H4FallbackStage2D::None;
            return result;
        }
        result.fallbackStage=H4FallbackStage2D::HybridCandidate;
    } else {
        result.fallbackStage=H4FallbackStage2D::LocalReduction;
    }

    result.fallback=buildPureCutCellFallback2D(
        domain,originalWalls,maxLevel,refinement,hybridPolicy);
    result.mode=result.fallback.valid()
        ?H4MeshMode2D::PureCutCellFallback:H4MeshMode2D::Failed;
    return result;
}

const char* h4MeshModeName(H4MeshMode2D mode) noexcept {
    switch (mode) {
    case H4MeshMode2D::Failed: return "failed";
    case H4MeshMode2D::Hybrid: return "hybrid";
    case H4MeshMode2D::PureCutCellFallback: return "pure_cutcell_fallback";
    }
    return "unknown";
}

const char* h4FallbackStageName(H4FallbackStage2D stage) noexcept {
    switch (stage) {
    case H4FallbackStage2D::None: return "none";
    case H4FallbackStage2D::RequestedLayers: return "requested_layers";
    case H4FallbackStage2D::LocalReduction: return "local_reduction";
    case H4FallbackStage2D::HybridCandidate: return "hybrid_candidate";
    }
    return "unknown";
}

const char* hybridMeshFailureReasonName(HybridMeshFailureReason2D reason) noexcept {
    switch (reason) {
    case HybridMeshFailureReason2D::None: return "none";
    case HybridMeshFailureReason2D::InvalidInput: return "invalid_input";
    case HybridMeshFailureReason2D::UnsupportedWallSemantics: return "unsupported_wall_semantics";
    case HybridMeshFailureReason2D::OuterEnvelopeOutsideDomain: return "outer_envelope_outside_domain";
    case HybridMeshFailureReason2D::InvalidOuterEnvelope: return "invalid_outer_envelope";
    case HybridMeshFailureReason2D::RemainderRefinementFailed: return "remainder_refinement_failed";
    case HybridMeshFailureReason2D::RemainderCutCellFailed: return "remainder_cutcell_failed";
    case HybridMeshFailureReason2D::RemainderStabilizationFailed: return "remainder_stabilization_failed";
    case HybridMeshFailureReason2D::LayerConversionFailed: return "layer_conversion_failed";
    case HybridMeshFailureReason2D::UnifiedTopologyFailed: return "unified_topology_failed";
    case HybridMeshFailureReason2D::NonConformalInterface: return "nonconformal_interface";
    case HybridMeshFailureReason2D::AreaConservationFailed: return "area_conservation_failed";
    case HybridMeshFailureReason2D::RegionClassificationConflict: return "region_classification_conflict";
    case HybridMeshFailureReason2D::QualityFailed: return "quality_failed";
    case HybridMeshFailureReason2D::SolverTopologyFailed: return "solver_topology_failed";
    case HybridMeshFailureReason2D::SolverQualityFailed: return "solver_quality_failed";
    case HybridMeshFailureReason2D::IoFailure: return "io_failure";
    }
    return "unknown";
}

const char* hybridCellKindName(HybridCellKind2D kind) noexcept {
    switch (kind) {
    case HybridCellKind2D::BoundaryLayer: return "boundary_layer";
    case HybridCellKind2D::RemainderCut: return "remainder_cut";
    case HybridCellKind2D::RemainderCartesian: return "remainder_cartesian";
    case HybridCellKind2D::Termination: return "termination";
    }
    return "unknown";
}

bool writeHybridLegacyVtk2D(const HybridMeshBuildResult2D& result,
                            const std::filesystem::path& path,
                            std::string* error) {
    if (!result.success() || !result.topology.valid() ||
        result.cellRecords.size() != result.topology.cells.size()) {
        setError(error, "cannot write failed or inconsistent H4-2 hybrid candidate");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open H4-2 VTK output");
        return false;
    }
    std::size_t cellListSize = 0U;
    for (const auto& cell : result.topology.cells) {
        cellListSize += cell.vertices.size() + 1U;
    }
    out << std::setprecision(17);
    out << "# vtk DataFile Version 3.0\n";
    out << "cartmesh2d H4-2 conformal hybrid mesh\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << result.topology.vertices.size() << " double\n";
    for (const auto& vertex : result.topology.vertices) {
        out << vertex.point.x << ' ' << vertex.point.y << " 0\n";
    }
    out << "CELLS " << result.topology.cells.size() << ' ' << cellListSize << "\n";
    for (const auto& cell : result.topology.cells) {
        out << cell.vertices.size();
        for (const auto id : cell.vertices) out << ' ' << id;
        out << '\n';
    }
    out << "CELL_TYPES " << result.topology.cells.size() << "\n";
    for (std::size_t i = 0; i < result.topology.cells.size(); ++i) out << "7\n";
    out << "CELL_DATA " << result.topology.cells.size() << "\n";
    out << "SCALARS hybrid_kind int 1\nLOOKUP_TABLE default\n";
    for (const auto& record : result.cellRecords) {
        out << static_cast<int>(record.kind) << '\n';
    }
    out << "SCALARS layer_index int 1\nLOOKUP_TABLE default\n";
    for (const auto& record : result.cellRecords) {
        out << (record.layerIndex ? static_cast<long long>(*record.layerIndex) : -1LL)
            << '\n';
    }
    out << "SCALARS geometry_area double 1\nLOOKUP_TABLE default\n";
    for (const auto& cell : result.topology.cells) out << cell.geometryArea << '\n';
    if (!out.good()) {
        setError(error, "failed while writing H4-2 VTK output");
        return false;
    }
    return true;
}

bool writeHybridReportJson2D(const HybridMeshBuildResult2D& result,
                             const std::filesystem::path& path,
                             std::string* error) {
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open H4-2 JSON report");
        return false;
    }
    out << std::setprecision(17);
    out << "{\n  \"hybrid_status\": \""
        << (result.success() ? "success" : "failed") << "\",\n";
    if (!result.success()) {
        out << "  \"failure_reason\": \""
            << hybridMeshFailureReasonName(result.failure.reason) << "\",\n";
        out << "  \"message\": ";
        writeJsonString(out, result.failure.message);
        out << ",\n  \"strip_id\": ";
        if (result.failure.stripId) out << *result.failure.stripId; else out << "null";
        out << ",\n  \"leaf_id\": ";
        if (result.failure.leafId) out << *result.failure.leafId; else out << "null";
        out << ",\n  \"cell_id\": ";
        if (result.failure.cellId) out << *result.failure.cellId; else out << "null";
        out << ",\n  \"edge_id\": ";
        if (result.failure.edgeId) out << *result.failure.edgeId; else out << "null";
        out << "\n}\n";
    } else {
        const auto& metrics = result.metrics;
        const auto& interface = result.interfaceAudit;
        out << "  \"failure_reason\": \"none\",\n";
        out << "  \"quadtree_leaf_count\": " << metrics.quadtreeLeafCount << ",\n";
        out << "  \"boundary_layer_cell_count\": " << metrics.boundaryLayerCellCount << ",\n";
        out << "  \"requested_boundary_layer_cell_count\": "
            << metrics.requestedBoundaryLayerCellCount << ",\n";
        out << "  \"zero_layer_column_count\": "
            << metrics.zeroLayerColumnCount << ",\n";
        out << "  \"termination_cell_count\": "
            << metrics.terminationCellCount << ",\n";
        out << "  \"termination_edge_count\": "
            << metrics.terminationEdgeCount << ",\n";
        out << "  \"remainder_cut_cell_count\": " << metrics.remainderCutCellCount << ",\n";
        out << "  \"remainder_cartesian_cell_count\": " << metrics.remainderCartesianCellCount << ",\n";
        out << "  \"transition_polygon_count\": " << metrics.transitionPolygonCount << ",\n";
        out << "  \"termination_growth_ratio\": "
            << metrics.terminationGrowthRatio << ",\n";
        out << "  \"vertex_count\": " << metrics.unifiedVertexCount << ",\n";
        out << "  \"edge_count\": " << metrics.unifiedEdgeCount << ",\n";
        out << "  \"cell_count\": " << metrics.unifiedCellCount << ",\n";
        out << "  \"solid_area\": " << metrics.solidArea << ",\n";
        out << "  \"outer_envelope_area\": " << metrics.outerEnvelopeArea << ",\n";
        out << "  \"layer_area\": " << metrics.layerArea << ",\n";
        out << "  \"remainder_area\": " << metrics.remainderArea << ",\n";
        out << "  \"expected_fluid_area\": " << metrics.expectedFluidArea << ",\n";
        out << "  \"actual_fluid_area\": " << metrics.actualFluidArea << ",\n";
        out << "  \"area_error\": " << metrics.areaError << ",\n";
        out << "  \"interface_edge_count\": " << interface.interfaceEdgeCount << ",\n";
        out << "  \"interface_vertex_count\": " << interface.interfaceVertexCount << ",\n";
        out << "  \"single_owner_interface_edges\": " << interface.singleOwnerInterfaceEdges << ",\n";
        out << "  \"wrong_cell_pair_interface_edges\": " << interface.wrongCellPairInterfaceEdges << ",\n";
        out << "  \"non_two_valent_interface_vertices\": " << interface.nonTwoValentInterfaceVertices << ",\n";
        out << "  \"expected_interface_length\": " << interface.expectedLength << ",\n";
        out << "  \"actual_interface_length\": " << interface.actualLength << ",\n";
        out << "  \"interface_length_error\": " << interface.lengthError << ",\n";
        out << "  \"topology_valid\": " << (result.topology.valid() ? "true" : "false") << ",\n";
        out << "  \"mesh_quality_valid\": " << (result.meshQuality.valid() ? "true" : "false") << ",\n";
        out << "  \"solver_quality_valid\": " << (result.solverQuality.valid() ? "true" : "false") << ",\n";
        out << "  \"solver_quality_issue_count\": " << result.solverQuality.issues.size() << ",\n";
        out << "  \"max_non_orthogonality_deg\": "
            << result.solverQuality.maxNonOrthogonalityDeg << ",\n";
        out << "  \"min_face_weight\": "
            << result.solverQuality.minFaceWeight << ",\n";
        out << "  \"min_volume_ratio\": "
            << result.solverQuality.minVolumeRatio << ",\n";
        out << "  \"balance_violations_after\": " << result.balance.violationsAfter << "\n}\n";
    }
    if (!out.good()) {
        setError(error, "failed while writing H4-2 JSON report");
        return false;
    }
    return true;
}

} // namespace cartmesh2d
