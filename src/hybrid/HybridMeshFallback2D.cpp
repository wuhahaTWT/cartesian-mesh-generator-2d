#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include "cartmesh2d/cutcell/CutCell2D.hpp"
#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/stabilization/SmallCell2D.hpp"
#include "cartmesh2d/topology/Topology2D.hpp"

#include <chrono>
#include <cmath>
#include <optional>
#include <sstream>

namespace cartmesh2d {
namespace {

using H4ProfileClock = std::chrono::steady_clock;

[[nodiscard]] double h4ProfileSeconds(H4ProfileClock::time_point start) noexcept {
    return std::chrono::duration<double>(H4ProfileClock::now() - start).count();
}

} // namespace

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
            leaf,originalWalls,policy.fluidRegion,policy.tolerance);
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
    result.expectedFluidArea=policy.fluidRegion == FluidRegion2D::Interior
        ? originalWalls.area(policy.tolerance)
        : domain.width()*domain.height()-originalWalls.area(policy.tolerance);
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
    smallPolicy.areaFractionThreshold=policy.remainderSmallCellAreaFraction;
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
    if (!result.solverQuality.exportable()) {
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
    const auto totalStart=H4ProfileClock::now();
    const auto startConformalBuilds=conformalHybridBuildCount2D();
    const auto finish=[&]() -> RobustH4BuildResult2D& {
        result.profile.totalSeconds=h4ProfileSeconds(totalStart);
        result.profile.conformalHybridBuildCalls=
            conformalHybridBuildCount2D()-startConformalBuilds;
        return result;
    };

    auto stageStart=H4ProfileClock::now();
    auto inwardPolicy=layerPolicy;
    if (hybridPolicy.fluidRegion == FluidRegion2D::Interior)
        inwardPolicy.permitConcaveTerminationMarching=true;
    result.requestedLayerCandidate=hybridPolicy.fluidRegion == FluidRegion2D::Interior
        ?buildLocallyReducedBoundaryLayerStrips2D(wallChains,layerParameters,inwardPolicy)
        :buildBoundaryLayerStrips2D(wallChains,layerParameters,inwardPolicy);
    result.profile.requestedLayerSeconds=h4ProfileSeconds(stageStart);
    if (result.requestedLayerCandidate.success()) {
        stageStart=H4ProfileClock::now();
        ++result.profile.requestedHybridAttempts;
        result.hybridCandidate=buildAutomaticHybridWithConstruction2D(
            result.requestedLayerCandidate,domain,originalWalls,maxLevel,refinement,
            hybridPolicy);
        result.profile.requestedHybridSeconds=h4ProfileSeconds(stageStart);
        if (result.hybridCandidate.success()) {
            result.mode=H4MeshMode2D::Hybrid;
            return finish();
        }
        result.fallbackStage=H4FallbackStage2D::HybridCandidate;
    } else {
        result.fallbackStage=H4FallbackStage2D::RequestedLayers;
    }

    if (hybridPolicy.fluidRegion != FluidRegion2D::Interior) {
    stageStart=H4ProfileClock::now();
    inwardPolicy.permitConcaveTerminationMarching=true;
    result.localLayerCandidate=buildLocallyReducedBoundaryLayerStrips2D(
        wallChains,layerParameters,inwardPolicy);
    result.profile.localLayerSeconds=h4ProfileSeconds(stageStart);
    if (result.localLayerCandidate.success()) {
        stageStart=H4ProfileClock::now();
        ++result.profile.localHybridAttempts;
        result.hybridCandidate=buildAutomaticHybridWithConstruction2D(
            result.localLayerCandidate,domain,originalWalls,maxLevel,refinement,
            hybridPolicy);
        result.profile.localHybridSeconds=h4ProfileSeconds(stageStart);
        if (result.hybridCandidate.success()) {
            result.mode=H4MeshMode2D::Hybrid;
            result.fallbackStage=H4FallbackStage2D::None;
            return finish();
        }
        result.fallbackStage=H4FallbackStage2D::HybridCandidate;
    } else {
        result.fallbackStage=H4FallbackStage2D::LocalReduction;
    }

    }

    stageStart=H4ProfileClock::now();
    ++result.profile.pureCutCellFallbackAttempts;
    result.fallback=buildPureCutCellFallback2D(
        domain,originalWalls,maxLevel,refinement,hybridPolicy);
    result.profile.pureCutCellFallbackSeconds=h4ProfileSeconds(stageStart);
    result.mode=result.fallback.valid()
        ?H4MeshMode2D::PureCutCellFallback:H4MeshMode2D::Failed;
    return finish();
}


} // namespace cartmesh2d
