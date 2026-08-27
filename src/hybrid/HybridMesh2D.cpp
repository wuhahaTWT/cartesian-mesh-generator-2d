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

constexpr std::uint64_t kLayerSourceKeyBase = std::uint64_t{1} << 63U;

struct SourceMetadata2D {
    HybridCellKind2D kind = HybridCellKind2D::RemainderCartesian;
    std::optional<std::size_t> layerIndex;
    std::optional<std::size_t> wallSegment;
};

[[nodiscard]] double polygonBoundsArea(const Polygon2D& polygon) noexcept {
    const auto bounds = polygon.bounds();
    return (bounds.max.x - bounds.min.x) * (bounds.max.y - bounds.min.y);
}

[[nodiscard]] double segmentLength(const Point2D& a, const Point2D& b) noexcept {
    return std::sqrt(squaredNorm(b - a));
}

[[nodiscard]] double regionPerimeter(const BoundaryRegion2D& region) noexcept {
    double result = 0.0;
    for (const auto& loop : region.loops()) {
        const auto& vertices = loop.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            result += segmentLength(vertices[i], vertices[(i + 1U) % vertices.size()]);
        }
    }
    return result;
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
           policy.interfaceToleranceMultiplier >= 1.0;
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
    for (std::size_t stripId = 0; stripId < boundaryLayers.strips.size(); ++stripId) {
        const auto& strip = boundaryLayers.strips[stripId];
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

    std::optional<Quadtree2D> remainderTree;
    try {
        remainderTree.emplace(domain, remainderMaxLevel, outerRegion, policy.tolerance);
        remainderTree->refine(outerRegion, remainderRefinement, policy.tolerance);
    } catch (const std::exception& exception) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      std::string("remainder refinement failed: ") + exception.what());
    }
    QuadtreeBalanceReport2D balance;
    try {
        balance = remainderTree->enforceTwoToOneBalance(outerRegion, policy.tolerance);
    } catch (const std::exception& exception) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      std::string("remainder 2:1 balance failed: ") + exception.what());
    }
    if (balance.violationsAfter != 0U || !remainderTree->deterministicOrderingValid()) {
        return failed(HybridMeshFailureReason2D::RemainderRefinementFailed,
                      "remainder quadtree is not deterministically ordered and 2:1 balanced");
    }

    std::vector<CutCell2D> sourceCells;
    std::map<std::size_t, SourceMetadata2D> metadata;
    sourceCells.reserve(remainderTree->leaves().size());
    std::size_t nextSourceId = 0U;
    std::size_t remainderCartesianCount = 0U;
    std::size_t remainderCutCount = 0U;
    std::size_t transitionPolygonCount = 0U;
    double remainderArea = 0.0;
    for (const auto& leaf : remainderTree->leaves()) {
        auto components = buildCutCells(leaf, outerRegion, FluidRegion2D::Exterior,
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
            const auto centroidState = outerRegion.classifyPoint(*component.centroid,
                                                                  policy.tolerance);
            if (centroidState == PointInPolygon::Inside) {
                return failed(HybridMeshFailureReason2D::RegionClassificationConflict,
                              "remainder cell centroid lies inside the outer envelope",
                              std::nullopt, leaf.id);
            }
            component.sourceId = nextSourceId;
            component.sourceKey = leaf.key;
            const auto kind = component.kind == CutCellKind::Full
                ? HybridCellKind2D::RemainderCartesian
                : HybridCellKind2D::RemainderCut;
            metadata.emplace(nextSourceId, SourceMetadata2D{kind, std::nullopt,
                                                            std::nullopt});
            if (kind == HybridCellKind2D::RemainderCartesian) {
                ++remainderCartesianCount;
            } else {
                ++remainderCutCount;
                if (component.fluidPolygon.vertices.size() != 4U) {
                    ++transitionPolygonCount;
                }
            }
            remainderArea += component.area;
            sourceCells.push_back(std::move(component));
            ++nextSourceId;
        }
    }

    double layerArea = 0.0;
    std::size_t layerCellCount = 0U;
    std::size_t globalLayerId = 0U;
    for (std::size_t stripId = 0; stripId < boundaryLayers.strips.size(); ++stripId) {
        const auto& strip = boundaryLayers.strips[stripId];
        for (const auto& layerCell : strip.cells) {
            CutCell2D source;
            source.sourceId = nextSourceId;
            source.sourceKey = kLayerSourceKeyBase + globalLayerId;
            source.kind = CutCellKind::Cut;
            for (const auto vertexId : layerCell.vertices) {
                source.fluidPolygon.vertices.push_back(strip.vertices[vertexId].point);
            }
            source.backgroundBounds = source.fluidPolygon.bounds();
            source.area = source.fluidPolygon.signedArea();
            const double boundsArea = polygonBoundsArea(source.fluidPolygon);
            source.areaFraction = boundsArea > 0.0
                ? std::clamp(source.area / boundsArea, 0.0, 1.0) : 0.0;
            source.centroid = source.fluidPolygon.centroid(policy.tolerance);
            if (!(source.area > 0.0) || !source.centroid || !(boundsArea > 0.0)) {
                return failed(HybridMeshFailureReason2D::LayerConversionFailed,
                              "H4-1 layer quad could not become a positive topology source",
                              stripId, std::nullopt, layerCell.id);
            }
            const auto originalState = originalWalls.classifyPoint(*source.centroid,
                                                                    policy.tolerance);
            const auto outerState = outerRegion.classifyPoint(*source.centroid,
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
            metadata.emplace(nextSourceId,
                             SourceMetadata2D{HybridCellKind2D::BoundaryLayer,
                                              layerCell.layer,
                                              layerCell.wallSegment});
            layerArea += source.area;
            sourceCells.push_back(std::move(source));
            ++nextSourceId;
            ++globalLayerId;
            ++layerCellCount;
        }
    }

    TopologyMesh2D topology = buildGlobalTopology(sourceCells, domain, originalWalls,
                                                   policy.tolerance);
    if (!topology.valid()) {
        const std::string detail = topology.issues.empty()
            ? "unified global topology audit failed"
            : topology.issues.front().message;
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed, detail);
    }

    // Rebuild a second time from immutable sources. This is both a transaction
    // boundary check and a direct determinism guard for canonical interface IDs.
    const TopologyMesh2D repeatedTopology = buildGlobalTopology(
        sourceCells, domain, originalWalls, policy.tolerance);
    if (!repeatedTopology.valid() || !sameTopology(topology, repeatedTopology)) {
        return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                      "repeated unified topology build changed IDs or connectivity");
    }

    std::vector<HybridCellRecord2D> records;
    records.reserve(topology.cells.size());
    for (const auto& cell : topology.cells) {
        const auto found = metadata.find(cell.sourceId);
        if (found == metadata.end()) {
            return failed(HybridMeshFailureReason2D::UnifiedTopologyFailed,
                          "topology cell lost its hybrid source metadata",
                          std::nullopt, std::nullopt, cell.id);
        }
        records.push_back({cell.id, cell.sourceId, found->second.kind,
                           found->second.layerIndex, found->second.wallSegment});
    }

    HybridInterfaceAudit2D interfaceAudit;
    interfaceAudit.expectedLength = regionPerimeter(outerRegion);
    std::map<std::size_t, std::size_t> interfaceVertexDegree;
    for (const auto& edge : topology.edges) {
        const auto& a = topology.vertices[edge.v0].point;
        const auto& b = topology.vertices[edge.v1].point;
        if (!edgeOnRegionBoundary(a, b, outerRegion, policy.tolerance)) continue;
        ++interfaceAudit.interfaceEdgeCount;
        interfaceAudit.actualLength += segmentLength(a, b);
        ++interfaceVertexDegree[edge.v0];
        ++interfaceVertexDegree[edge.v1];
        if (!edge.neighbour) {
            ++interfaceAudit.singleOwnerInterfaceEdges;
            continue;
        }
        const auto ownerKind = records[edge.owner].kind;
        const auto neighbourKind = records[*edge.neighbour].kind;
        const bool ownerLayer = ownerKind == HybridCellKind2D::BoundaryLayer;
        const bool neighbourLayer = neighbourKind == HybridCellKind2D::BoundaryLayer;
        if (ownerLayer == neighbourLayer) {
            ++interfaceAudit.wrongCellPairInterfaceEdges;
        }
    }
    interfaceAudit.interfaceVertexCount = interfaceVertexDegree.size();
    for (const auto& [vertexId, degree] : interfaceVertexDegree) {
        (void)vertexId;
        if (degree != 2U) ++interfaceAudit.nonTwoValentInterfaceVertices;
    }
    interfaceAudit.lengthError = interfaceAudit.actualLength -
                                 interfaceAudit.expectedLength;
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
        sourceCells.begin(), sourceCells.end(), 0.0,
        [](double sum, const CutCell2D& cell) { return sum + cell.area; });
    const double areaError = actualFluidArea - expectedFluidArea;
    const double areaTolerance = policy.areaToleranceMultiplier *
        (policy.tolerance.absolute * policy.tolerance.absolute +
         policy.tolerance.relative * std::max(1.0, expectedFluidArea));
    if (std::abs(areaError) > areaTolerance ||
        std::abs(layerArea - (outerArea - solidArea)) > areaTolerance ||
        std::abs(remainderArea - (domainArea - outerArea)) > areaTolerance) {
        return failed(HybridMeshFailureReason2D::AreaConservationFailed,
                      "hybrid layer/remainder areas do not close to domain minus solid");
    }

    auto meshQuality = evaluateMeshQuality(topology, sourceCells, nullptr,
                                           policy.tolerance);
    if (!meshQuality.valid()) {
        const std::string detail = meshQuality.issues.empty()
            ? "hybrid mesh quality audit failed"
            : meshQuality.issues.front().message;
        return failed(HybridMeshFailureReason2D::QualityFailed, detail);
    }
    auto solverQuality = evaluateSolverQuality2D(topology, {}, policy.tolerance);

    HybridMeshBuildResult2D result;
    result.status = HybridMeshStatus2D::Success;
    result.strips = boundaryLayers.strips;
    result.outerEnvelopeRegion = std::move(outerRegion);
    result.sourceCells = std::move(sourceCells);
    result.topology = std::move(topology);
    result.cellRecords = std::move(records);
    result.interfaceAudit = interfaceAudit;
    result.meshQuality = std::move(meshQuality);
    result.solverQuality = std::move(solverQuality);
    result.balance = balance;
    result.metrics.quadtreeLeafCount = remainderTree->leaves().size();
    result.metrics.remainderCartesianCellCount = remainderCartesianCount;
    result.metrics.remainderCutCellCount = remainderCutCount;
    result.metrics.boundaryLayerCellCount = layerCellCount;
    result.metrics.transitionPolygonCount = transitionPolygonCount;
    result.metrics.unifiedVertexCount = result.topology.vertices.size();
    result.metrics.unifiedEdgeCount = result.topology.edges.size();
    result.metrics.unifiedCellCount = result.topology.cells.size();
    result.metrics.solidArea = solidArea;
    result.metrics.outerEnvelopeArea = outerArea;
    result.metrics.layerArea = layerArea;
    result.metrics.remainderArea = remainderArea;
    result.metrics.expectedFluidArea = expectedFluidArea;
    result.metrics.actualFluidArea = actualFluidArea;
    result.metrics.areaError = areaError;
    return result;
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
    case HybridMeshFailureReason2D::LayerConversionFailed: return "layer_conversion_failed";
    case HybridMeshFailureReason2D::UnifiedTopologyFailed: return "unified_topology_failed";
    case HybridMeshFailureReason2D::NonConformalInterface: return "nonconformal_interface";
    case HybridMeshFailureReason2D::AreaConservationFailed: return "area_conservation_failed";
    case HybridMeshFailureReason2D::RegionClassificationConflict: return "region_classification_conflict";
    case HybridMeshFailureReason2D::QualityFailed: return "quality_failed";
    case HybridMeshFailureReason2D::IoFailure: return "io_failure";
    }
    return "unknown";
}

const char* hybridCellKindName(HybridCellKind2D kind) noexcept {
    switch (kind) {
    case HybridCellKind2D::BoundaryLayer: return "boundary_layer";
    case HybridCellKind2D::RemainderCut: return "remainder_cut";
    case HybridCellKind2D::RemainderCartesian: return "remainder_cartesian";
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
        out << "  \"remainder_cut_cell_count\": " << metrics.remainderCutCellCount << ",\n";
        out << "  \"remainder_cartesian_cell_count\": " << metrics.remainderCartesianCellCount << ",\n";
        out << "  \"transition_polygon_count\": " << metrics.transitionPolygonCount << ",\n";
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
        out << "  \"balance_violations_after\": " << result.balance.violationsAfter << "\n}\n";
    }
    if (!out.good()) {
        setError(error, "failed while writing H4-2 JSON report");
        return false;
    }
    return true;
}

} // namespace cartmesh2d
