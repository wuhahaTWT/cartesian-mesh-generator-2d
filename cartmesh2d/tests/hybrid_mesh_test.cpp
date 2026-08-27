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

std::vector<Point2D> superellipse(std::size_t count,double a,double b) {
    std::vector<Point2D> points;
    for (std::size_t i=0;i<count;++i) {
        const double angle=2.0*std::numbers::pi*static_cast<double>(i)/
                           static_cast<double>(count);
        const double cosine=std::cos(angle),sine=std::sin(angle);
        points.push_back({a*std::copysign(std::sqrt(std::abs(cosine)),cosine),
                          b*std::copysign(std::sqrt(std::abs(sine)),sine)});
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
    for (std::size_t i=0;i<std::min<std::size_t>(4U,result.solverQuality.issues.size());++i) {
        const auto& issue=result.solverQuality.issues[i];
        if (issue.edgeId>=result.solverTopology.edges.size()) continue;
        const auto& edge=result.solverTopology.edges[issue.edgeId];
        const auto& a=result.solverTopology.vertices[edge.v0].point;
        const auto& b=result.solverTopology.vertices[edge.v1].point;
        std::cerr<<" diagnostic edge="<<issue.edgeId<<" ("<<a.x<<','<<a.y
                 <<")-("<<b.x<<','<<b.y<<") owner="<<edge.owner<<" neighbour=";
        if (edge.neighbour) std::cerr<<*edge.neighbour;
        else std::cerr<<"none";
        std::cerr<<'\n';
        for (const auto cellId:{edge.owner,*edge.neighbour}) {
            std::cerr<<"  cell "<<cellId<<':' ;
            for (const auto vertex:result.solverTopology.cells[cellId].vertices) {
                const auto& point=result.solverTopology.vertices[vertex].point;
                std::cerr<<" ("<<point.x<<','<<point.y<<')';
            }
            std::cerr<<'\n';
        }
    }
}

} // namespace

int main() {
    if (std::getenv("H4_SWEEP")!=nullptr) {
        const BoundaryLoop wall(superellipse(24U,1.8,0.8));
        const auto layers=layersFor(wall,3U,0.015,1.15);
        QuadtreeRefinementPolicy2D refine; refine.minimumLevel=3U; refine.boundaryLevel=6U;
        const Domain2D sweepDomain{{{-2.8,-1.8},{2.8,1.8}}};
        for (const double factor:{1.8}) {
            HybridMeshPolicy2D policy; policy.transitionCellWidthMultiplier=factor;
            const auto value=buildConformalHybridMesh2D(
                layers,sweepDomain,BoundaryRegion2D(wall),6U,refine,policy);
            std::cerr<<"factor="<<factor<<" success="<<value.success();
            if (!value.success()) std::cerr<<" "<<value.failure.message;
            else std::cerr<<" nonorth="<<value.solverQuality.maxNonOrthogonalityDeg
                          <<" weight="<<value.solverQuality.minFaceWeight
                          <<" ratio="<<value.solverQuality.minVolumeRatio;
            std::cerr<<'\n';
        }
        return 0;
    }
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
    const auto circleHybrid = buildConformalHybridMesh2D(
        circleLayers, domain, BoundaryRegion2D(circleWall), 6U, refinement);
    if (!circleHybrid.success()) {
        std::cerr << "circle hybrid failure: "
                  << hybridMeshFailureReasonName(circleHybrid.failure.reason)
                  << " " << circleHybrid.failure.message << '\n';
        printSolverDiagnostics(circleHybrid);
    }
    check(circleHybrid.success(), "circle layer and remainder build one hybrid mesh");
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
              circleHybrid.interfaceAudit.interfaceEdgeCount > 32U,
              "circle outer envelope is split into a closed conformal interface");
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
    const auto ellipseHybrid = buildConformalHybridMesh2D(
        ellipseLayers, ellipseDomain, BoundaryRegion2D(ellipseWall), 6U,
        ellipseRefinement);
    if (!ellipseHybrid.success()) {
        std::cerr << "ellipse hybrid failure: "
                  << hybridMeshFailureReasonName(ellipseHybrid.failure.reason)
                  << " " << ellipseHybrid.failure.message << '\n';
        printSolverDiagnostics(ellipseHybrid);
    }
    check(ellipseHybrid.success(), "non-uniform-curvature ellipse hybrid succeeds");
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
        std::cout << "cartmesh2d H4-2 conformal hybrid tests: PASS\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " H4-2 test(s) failed\n";
    return EXIT_FAILURE;
}
