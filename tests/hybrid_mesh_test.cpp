#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

std::vector<Point2D> ellipse(std::size_t count, double a, double b,
                             Point2D center = {}) {
    std::vector<Point2D> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double angle = 2.0 * std::numbers::pi *
                             static_cast<double>(i) / static_cast<double>(count);
        points.push_back({center.x + a * std::cos(angle),
                          center.y + b * std::sin(angle)});
    }
    return points;
}

std::vector<Point2D> superellipse(std::size_t count, double a, double b) {
    std::vector<Point2D> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double angle = 2.0 * std::numbers::pi *
                             static_cast<double>(i) / static_cast<double>(count);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        points.push_back({a * std::copysign(std::sqrt(std::abs(cosine)), cosine),
                          b * std::copysign(std::sqrt(std::abs(sine)), sine)});
    }
    return points;
}

BoundaryLayerBuildResult2D layersFor(const BoundaryLoop& wall,
                                     std::size_t count,
                                     double first,
                                     double ratio) {
    const auto chain = makeClosedWallChain2D(wall, 0U, "wall_0");
    if (!chain.success()) return {};
    return buildBoundaryLayerStrip2D(
        *chain.chain,
        {count, LayerThicknessMode2D::FirstLayerThickness, first, ratio});
}

bool samePoints(const std::vector<BoundaryLayerVertex2D>& lhs,
                const std::vector<BoundaryLayerVertex2D>& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i].id != rhs[i].id || lhs[i].ring != rhs[i].ring ||
            lhs[i].chainVertex != rhs[i].chainVertex ||
            lhs[i].point.x != rhs[i].point.x ||
            lhs[i].point.y != rhs[i].point.y) return false;
    }
    return true;
}

void checkInterfacePairs(const HybridMeshBuildResult2D& result,
                         const std::string& label) {
    std::size_t checked = 0U;
    for (const auto& edge : result.topology.edges) {
        if (!edge.neighbour) continue;
        const auto ownerKind = result.cellRecords[edge.owner].kind;
        const auto neighbourKind = result.cellRecords[*edge.neighbour].kind;
        if ((ownerKind == HybridCellKind2D::BoundaryLayer) ==
            (neighbourKind == HybridCellKind2D::BoundaryLayer)) continue;
        ++checked;
        check(edge.v0 < result.topology.vertices.size() &&
              edge.v1 < result.topology.vertices.size(),
              label + " interface uses global vertex IDs");
    }
    check(checked == result.interfaceAudit.interfaceEdgeCount,
          label + " every geometric interface edge is one layer/remainder pair");
}

void printSolverDiagnostics(const HybridMeshBuildResult2D& result) {
    for (std::size_t i = 0;
         i < std::min<std::size_t>(4U, result.solverQuality.issues.size()); ++i) {
        const auto& issue = result.solverQuality.issues[i];
        if (issue.edgeId >= result.solverTopology.edges.size()) continue;
        const auto& edge = result.solverTopology.edges[issue.edgeId];
        const auto& a = result.solverTopology.vertices[edge.v0].point;
        const auto& b = result.solverTopology.vertices[edge.v1].point;
        std::cerr << " diagnostic edge=" << issue.edgeId << " (" << a.x << ',' << a.y
                  << ")-(" << b.x << ',' << b.y << ") owner=" << edge.owner
                  << " neighbour=";
        if (!edge.neighbour) {
            std::cerr << "none\n";
            continue;
        }
        std::cerr << *edge.neighbour << '\n';
        for (const auto cellId : {edge.owner, *edge.neighbour}) {
            std::cerr << "  cell " << cellId << ':';
            for (const auto vertex : result.solverTopology.cells[cellId].vertices) {
                const auto& point = result.solverTopology.vertices[vertex].point;
                std::cerr << " (" << point.x << ',' << point.y << ')';
            }
            std::cerr << '\n';
        }
    }
}

void checkAutomaticPlan(const BoundaryLayerBuildResult2D& layers,
                        const Domain2D& domain,
                        const QuadtreeRefinementPolicy2D& refinement,
                        const std::string& label) {
    const auto plan = resolveAutomaticHybridTransitionPlan2D(
        layers, domain, refinement);
    check(plan.has_value(), label + " automatic transition plan resolves");
    if (!plan) return;
    check(plan->ringCount >= 3U,
          label + " keeps a genuinely progressive transition fan");
    check(plan->finalTangentialSubdivision >= 4U,
          label + " has progressive tangential subdivision");
    check(plan->maxOuterEdgeLength /
              static_cast<double>(plan->finalTangentialSubdivision) <=
          2.0 * plan->targetCellSize * (1.0 + 1.0e-12),
          label + " final fan spacing is tied to remainder target h");
    check(plan->ringThickness >= plan->maxLastLayerSpacing,
          label + " transition radial spacing does not contract below last layer");
}

void checkSolverReady(const HybridMeshBuildResult2D& result,
                      const std::string& label) {
    check(result.success(), label + " hybrid build succeeds");
    if (!result.success()) {
        std::cerr << label << " failure: "
                  << hybridMeshFailureReasonName(result.failure.reason) << ' '
                  << result.failure.message << '\n';
        printSolverDiagnostics(result);
        return;
    }
    check(result.solverTopology.valid(), label + " solver topology is valid");
    check(result.solverQuality.valid(), label + " solver-quality gate passes");
    check(result.solverQuality.maxNonOrthogonalityDeg <= 70.0,
          label + " non-orthogonality stays within production gate");
    check(result.solverQuality.minFaceWeight >= 0.05,
          label + " face weight stays within production gate");
    check(result.solverQuality.minVolumeRatio >= 0.01,
          label + " neighbouring volume ratio stays within production gate");
    check(result.metrics.transitionRingCount >= 3U &&
          result.metrics.transitionRingThickness > 0.0,
          label + " records the automatic transition plan");
}

} // namespace

int main() {
    const Domain2D domain{{{-2.0, -2.0}, {2.0, 2.0}}};
    const BoundaryLoop circleWall(ellipse(32U, 1.0, 1.0, {0.07, 0.03}));
    const auto circleLayers = layersFor(circleWall, 4U, 0.02, 1.2);
    check(circleLayers.success(), "H4-1 circle prerequisite succeeds");
    const auto savedLayerVertices = circleLayers.success()
        ? circleLayers.strips.front().vertices
        : std::vector<BoundaryLayerVertex2D>{};
    QuadtreeRefinementPolicy2D refinement;
    refinement.minimumLevel = 3U;
    refinement.boundaryLevel = 6U;
    checkAutomaticPlan(circleLayers, domain, refinement, "circle");
    const auto circleHybrid = buildConformalHybridMesh2D(
        circleLayers, domain, BoundaryRegion2D(circleWall), 6U, refinement);
    checkSolverReady(circleHybrid, "circle");
    HybridMeshPolicy2D lineageOraclePolicy;
    lineageOraclePolicy.verifySourceLineageOracle=true;
    const auto circleLineageVerified=buildAutomaticHybridWithConstruction2D(
        circleLayers,domain,BoundaryRegion2D(circleWall),6U,refinement,
        lineageOraclePolicy);
    check(circleLineageVerified.success() &&
              circleLineageVerified.sourceLineageAudit.oracleVerified &&
              circleLineageVerified.sourceLineageAudit.mismatchedCells==0U &&
              circleLineageVerified.sourceLineageAudit.oracleCandidateChecks>
                  circleLineageVerified.sourceLineageAudit.lineageCandidateChecks,
          "R1-A opt-in source-lineage oracle matches the propagated production path");
    check(circleLayers.success() &&
          samePoints(savedLayerVertices, circleLayers.strips.front().vertices),
          "H4-2 transaction does not mutate the H4-1 strip");
    if (circleHybrid.success()) {
        check(circleHybrid.topology.valid() && circleHybrid.meshQuality.valid(),
              "circle unified topology and base quality pass");
        check(circleHybrid.metrics.boundaryLayerCellCount == 128U,
              "circle preserves all four H4-1 layers");
        check(circleHybrid.metrics.remainderCutCellCount > 0U &&
              circleHybrid.metrics.remainderCartesianCellCount > 0U,
              "circle contains true remainder Cut-cells and Cartesian cells");
        check(circleHybrid.metrics.transitionPolygonCount > 0U,
              "circle transition contains legal general polygons");
        check(circleHybrid.interfaceAudit.pass(1.0e-8) &&
              circleHybrid.interfaceAudit.interfaceEdgeCount ==
                  circleWall.vertices().size(),
              "circle preserves the fixed H4-1 outer-envelope partition as a closed interface");
        check(circleHybrid.solverInterfaceAudit.pass(1.0e-8),
              "circle solver repair preserves the H4-1 shared interface");
        check(std::abs(circleHybrid.metrics.areaError) < 1.0e-8,
              "circle hybrid conserves domain minus solid area");
        checkInterfacePairs(circleHybrid, "circle");
    }

    const Domain2D ellipseDomain{{{-3.0, -2.0}, {3.0, 2.0}}};
    const BoundaryLoop ellipseWall(ellipse(48U, 1.5, 0.75));
    const auto ellipseLayers = layersFor(ellipseWall, 3U, 0.025, 1.15);
    QuadtreeRefinementPolicy2D ellipseRefinement;
    ellipseRefinement.minimumLevel = 3U;
    ellipseRefinement.boundaryLevel = 6U;
    checkAutomaticPlan(ellipseLayers, ellipseDomain, ellipseRefinement, "ellipse");
    const auto ellipseHybrid = buildConformalHybridMesh2D(
        ellipseLayers, ellipseDomain, BoundaryRegion2D(ellipseWall), 6U,
        ellipseRefinement);
    checkSolverReady(ellipseHybrid, "ellipse");
    if (ellipseHybrid.success()) {
        check(ellipseHybrid.interfaceAudit.pass(1.0e-8),
              "ellipse interface is closed and two-owner conformal");
        check(ellipseHybrid.topology.audit.nonManifoldEdges == 0U &&
              ellipseHybrid.topology.audit.unclassifiedBoundaryEdges == 0U,
              "ellipse has no non-manifold or unclassified interface edge");
        check(std::abs(ellipseHybrid.metrics.areaError) < 1.0e-8,
              "ellipse hybrid area closes");
        checkInterfacePairs(ellipseHybrid, "ellipse");
    }

    const Domain2D superDomain{{{-2.8, -1.8}, {2.8, 1.8}}};
    const BoundaryLoop superWall(superellipse(24U, 1.8, 0.8));
    const auto superLayers = layersFor(superWall, 3U, 0.015, 1.15);
    QuadtreeRefinementPolicy2D superRefinement;
    superRefinement.minimumLevel = 3U;
    superRefinement.boundaryLevel = 6U;
    checkAutomaticPlan(superLayers, superDomain, superRefinement, "superellipse");
    const auto superHybrid = buildConformalHybridMesh2D(
        superLayers, superDomain, BoundaryRegion2D(superWall), 6U,
        superRefinement);
    checkSolverReady(superHybrid, "superellipse");
    if (superHybrid.success()) {
        const auto& micro=superHybrid.qualityContract.ordinaryMetrics.at("face_length_over_local_background_h");
        check(micro.worst>=QualityContract2D{}.transition.faceOverLocalBackgroundH.hard,
              "Q2 superellipse final solver face/local_h satisfies Q1 hard limit");
        check(!superHybrid.canonicalizedIntersections.empty(),
              "Q2 retains canonical intersection provenance");
        check(superHybrid.interfaceAudit.pass(1.0e-8) &&
              superHybrid.solverInterfaceAudit.pass(1.0e-8),
              "superellipse geometric and solver interfaces stay conformal");
        check(std::abs(superHybrid.metrics.areaError) < 1.0e-8,
              "superellipse hybrid area closes");
    }

    // Generalization check: change both boundary resolution and H4-1 layer
    // parameters. No transition width/ring count is supplied by this test.
    QuadtreeRefinementPolicy2D variedRefinement;
    variedRefinement.minimumLevel = 3U;
    variedRefinement.boundaryLevel = 7U;
    const auto variedLayers = layersFor(circleWall, 3U, 0.015, 1.10);
    checkAutomaticPlan(variedLayers, domain, variedRefinement, "varied circle");
    const auto variedHybrid = buildConformalHybridMesh2D(
        variedLayers, domain, BoundaryRegion2D(circleWall), 7U,
        variedRefinement);
    checkSolverReady(variedHybrid, "varied circle");
    if (circleHybrid.success() && variedHybrid.success()) {
        check(circleHybrid.metrics.transitionTargetCellSize !=
              variedHybrid.metrics.transitionTargetCellSize,
              "automatic fan responds to changed boundaryLevel");
    }

    const Domain2D tooSmall{{{-1.02, -1.02}, {1.02, 1.02}}};
    const auto rejected = buildConformalHybridMesh2D(
        circleLayers, tooSmall, BoundaryRegion2D(circleWall), 5U,
        QuadtreeRefinementPolicy2D{2U, 5U, {}, {}});
    check(!rejected.success() &&
          rejected.failure.reason ==
              HybridMeshFailureReason2D::OuterEnvelopeOutsideDomain,
          "outer envelope touching/exceeding domain fails before remainder mutation");
    check(circleLayers.success() &&
          samePoints(savedLayerVertices, circleLayers.strips.front().vertices),
          "failed hybrid candidate also leaves H4-1 input unchanged");

    if (failures == 0) {
        std::cout << "cartmesh2d H4 solver-ready hybrid tests: PASS\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " H4 solver-ready test(s) failed\n";
    return EXIT_FAILURE;
}
