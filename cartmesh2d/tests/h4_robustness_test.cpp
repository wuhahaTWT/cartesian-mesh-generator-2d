#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool condition,const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr<<"FAIL: "<<message<<'\n';
    }
}

std::vector<Point2D> circle(std::size_t count,double radius) {
    std::vector<Point2D> points;
    points.reserve(count);
    for (std::size_t i=0;i<count;++i) {
        const double angle=2.0*std::numbers::pi*static_cast<double>(i)/
                           static_cast<double>(count);
        points.push_back({radius*std::cos(angle),radius*std::sin(angle)});
    }
    return points;
}

BoundaryLayerBuildResult2D buildLayers(
    const BoundaryRegion2D& walls,const LayerParameters2D& parameters) {
    std::vector<WallChain2D> chains;
    for (std::size_t i=0;i<walls.loops().size();++i) {
        const auto chain=makeClosedWallChain2D(
            walls.loops()[i],i,"wall_"+std::to_string(i));
        if (!chain.success()) return {};
        chains.push_back(*chain.chain);
    }
    return buildBoundaryLayerStrips2D(chains,parameters);
}

QuadtreeRefinementPolicy2D refinement(std::size_t level=6U) {
    QuadtreeRefinementPolicy2D policy;
    policy.minimumLevel=3U;
    policy.boundaryLevel=level;
    return policy;
}

void checkFallback(const RobustH4BuildResult2D& result,
                   H4FallbackStage2D stage,const std::string& label) {
    if (!result.success()) {
        std::cerr<<label<<" fallback failure: "
                 <<result.fallback.failureMessage<<'\n';
    }
    check(result.success(),label+" produces a usable fallback mesh");
    check(result.mode==H4MeshMode2D::PureCutCellFallback,
          label+" is explicitly marked pure Cut-cell fallback");
    check(result.fallbackStage==stage,label+" records the failing H4 stage");
    check(!result.layerEnabled(),label+" never reports boundary layers enabled");
    if (!result.fallback.valid()) return;
    check(result.fallback.topology.valid() &&
          result.fallback.solverTopology.valid(),
          label+" fallback topology is valid");
    check(result.fallback.meshQuality.valid() &&
          result.fallback.solverQuality.valid(),
          label+" fallback passes unchanged solver-quality thresholds");
    check(std::abs(result.fallback.areaError)<1.0e-8,
          label+" fallback conserves fluid area");
    check(result.fallback.topology.audit.nonManifoldEdges==0U,
          label+" fallback has no non-manifold edges");
}

} // namespace

int main() {
    const LayerParameters2D standardLayers{
        3U,LayerThicknessMode2D::FirstLayerThickness,0.02,1.1};

    const BoundaryRegion2D circleWalls(BoundaryLoop(circle(32U,1.0)));
    const Domain2D circleDomain{{{-2.0,-2.0},{2.0,2.0}}};
    const auto circleLayers=buildLayers(circleWalls,standardLayers);
    const auto successful=buildRobustH4Mesh2D(
        circleLayers,circleDomain,circleWalls,6U,refinement());
    check(successful.success() && successful.mode==H4MeshMode2D::Hybrid,
          "smooth convex wall commits the hybrid candidate");
    check(successful.layerEnabled() &&
          successful.fallbackStage==H4FallbackStage2D::None,
          "successful hybrid does not invoke fallback");

    const BoundaryRegion2D concaveWalls(BoundaryLoop({
        {0.0,0.0},{3.0,0.0},{3.0,3.0},{2.0,3.0},
        {2.0,1.0},{1.0,1.0},{1.0,3.0},{0.0,3.0}}));
    const auto concaveLayers=buildLayers(concaveWalls,standardLayers);
    check(!concaveLayers.success() &&
          concaveLayers.failure.reason==BoundaryLayerFailureReason2D::ConcaveCorner,
          "concave minimum case retains its structured layer failure");
    const auto concave=buildRobustH4Mesh2D(
        concaveLayers,Domain2D{{{-1.0,-1.0},{4.0,4.0}}},concaveWalls,
        6U,refinement());
    checkFallback(concave,H4FallbackStage2D::BoundaryLayer,"concave wall");
    check(concave.layerCandidate.failure.vertexId.has_value(),
          "concave failure records the first vertex");

    const BoundaryRegion2D sharpWalls(BoundaryLoop({
        {0.0,0.0},{3.0,0.0},{0.05,0.1}}));
    const LayerParameters2D sharpParameters{
        2U,LayerThicknessMode2D::FirstLayerThickness,0.01,1.0};
    const auto sharpLayers=buildLayers(sharpWalls,sharpParameters);
    check(!sharpLayers.success() &&
          sharpLayers.failure.reason==BoundaryLayerFailureReason2D::SharpCorner,
          "sharp minimum case retains its structured layer failure");
    const auto sharp=buildRobustH4Mesh2D(
        sharpLayers,Domain2D{{{-1.0,-1.0},{4.0,1.1}}},sharpWalls,
        7U,refinement(7U));
    checkFallback(sharp,H4FallbackStage2D::BoundaryLayer,"sharp wall");
    check(sharp.layerCandidate.failure.vertexId.has_value(),
          "sharp failure records the first vertex");
    const auto rejectedFallback=buildRobustH4Mesh2D(
        sharpLayers,Domain2D{{{-1.0,-1.0},{4.0,2.0}}},sharpWalls,
        7U,refinement(7U));
    check(!rejectedFallback.success() &&
          rejectedFallback.mode==H4MeshMode2D::Failed &&
          rejectedFallback.fallback.failureMessage.find("solver quality failed")!=
              std::string::npos,
          "fallback candidate that misses solver quality fails closed with metrics");

    const BoundaryRegion2D gapWalls({
        BoundaryLoop({{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}}),
        BoundaryLoop({{1.2,0.0},{2.2,0.0},{2.2,1.0},{1.2,1.0}})});
    const LayerParameters2D gapParameters{
        1U,LayerThicknessMode2D::FirstLayerThickness,0.11,1.0};
    const auto gapLayers=buildLayers(gapWalls,gapParameters);
    check(!gapLayers.success() &&
          gapLayers.failure.reason==
              BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit &&
          gapLayers.failure.safeThickness.has_value(),
          "narrow-gap case retains requested and safe thickness");
    const auto gap=buildRobustH4Mesh2D(
        gapLayers,Domain2D{{{-1.0,-1.0},{3.2,2.0}}},gapWalls,
        7U,refinement(7U));
    checkFallback(gap,H4FallbackStage2D::BoundaryLayer,"narrow gap");

    HybridMeshPolicy2D forcedHybridFailure;
    forcedHybridFailure.transitionCellWidthMultiplier=20.0;
    const auto rejectedHybrid=buildRobustH4Mesh2D(
        circleLayers,circleDomain,circleWalls,6U,refinement(),
        forcedHybridFailure);
    checkFallback(rejectedHybrid,H4FallbackStage2D::HybridCandidate,
                  "rejected hybrid candidate");
    check(rejectedHybrid.layerCandidate.success() &&
          !rejectedHybrid.hybridCandidate.success(),
          "H4-2 failure preserves the successful H4-1 candidate");

    const auto repeat=buildRobustH4Mesh2D(
        concaveLayers,Domain2D{{{-1.0,-1.0},{4.0,4.0}}},concaveWalls,
        6U,refinement());
    check(repeat.fallback.topology.cells.size()==
              concave.fallback.topology.cells.size() &&
          repeat.fallback.topology.edges.size()==
              concave.fallback.topology.edges.size(),
          "fallback construction is deterministic in size");

    const auto reportPath=std::filesystem::temp_directory_path()/
                          "cartmesh2d-h4-3-concave.json";
    std::string reportError;
    check(writeRobustH4ReportJson2D(concave,reportPath,&reportError),
          "structured H4-3 report is writable");
    std::ifstream reportStream(reportPath);
    std::ostringstream reportText;
    reportText<<reportStream.rdbuf();
    check(reportText.str().find("\"mesh_mode\": \"pure_cutcell_fallback\"")!=
              std::string::npos &&
          reportText.str().find("\"layer_status\": \"failed\"")!=
              std::string::npos &&
          reportText.str().find("\"layer_failure_reason\": \"concave_corner\"")!=
              std::string::npos,
          "report distinguishes fallback from a successful boundary layer");

    if (failures==0) {
        std::cout<<"cartmesh2d H4-3 robustness and rollback tests: PASS\n";
        return EXIT_SUCCESS;
    }
    std::cerr<<failures<<" H4-3 test(s) failed\n";
    return EXIT_FAILURE;
}
