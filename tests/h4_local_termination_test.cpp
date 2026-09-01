#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures=0;

void check(bool condition,const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr<<"FAIL: "<<message<<'\n';
    }
}

bool sameTopology(const TopologyMesh2D& lhs,const TopologyMesh2D& rhs) {
    if (lhs.vertices.size()!=rhs.vertices.size() ||
        lhs.edges.size()!=rhs.edges.size() || lhs.cells.size()!=rhs.cells.size()) {
        return false;
    }
    for (std::size_t i=0;i<lhs.vertices.size();++i) {
        if (lhs.vertices[i].point.x!=rhs.vertices[i].point.x ||
            lhs.vertices[i].point.y!=rhs.vertices[i].point.y) return false;
    }
    for (std::size_t i=0;i<lhs.cells.size();++i) {
        if (lhs.cells[i].vertices!=rhs.cells[i].vertices ||
            lhs.cells[i].edges!=rhs.cells[i].edges) return false;
    }
    for (std::size_t i=0;i<lhs.edges.size();++i) {
        const auto& a=lhs.edges[i];
        const auto& b=rhs.edges[i];
        if (a.v0!=b.v0 || a.v1!=b.v1 || a.owner!=b.owner ||
            a.neighbour!=b.neighbour || a.patch!=b.patch) return false;
    }
    return true;
}

LayerParameters2D fourLayers() {
    return {4U,LayerThicknessMode2D::FirstLayerThickness,0.012,1.15};
}

QuadtreeRefinementPolicy2D refinement(std::size_t level=7U) {
    QuadtreeRefinementPolicy2D result;
    result.minimumLevel=3U;
    result.boundaryLevel=level;
    return result;
}

std::vector<Point2D> subdividedRectangle(
    double x0,double x1,double y0,double y1,std::size_t divisions) {
    std::vector<Point2D> points;
    for (std::size_t i=0;i<divisions;++i) {
        const double t=static_cast<double>(i)/static_cast<double>(divisions);
        points.push_back({x0+(x1-x0)*t,y0});
    }
    for (std::size_t i=0;i<divisions;++i) {
        const double t=static_cast<double>(i)/static_cast<double>(divisions);
        points.push_back({x1,y0+(y1-y0)*t});
    }
    for (std::size_t i=0;i<divisions;++i) {
        const double t=static_cast<double>(i)/static_cast<double>(divisions);
        points.push_back({x1-(x1-x0)*t,y1});
    }
    for (std::size_t i=0;i<divisions;++i) {
        const double t=static_cast<double>(i)/static_cast<double>(divisions);
        points.push_back({x0,y1-(y1-y0)*t});
    }
    return points;
}

void printCounts(const BoundaryLayerBuildResult2D& result,
                 const std::string& label) {
    std::cerr<<label<<" count histogram:";
    for (const auto& strip:result.strips) {
        std::vector<std::size_t> histogram(strip.parameters.nLayers+1U,0U);
        for (const auto count:strip.actualLayerCounts) ++histogram[count];
        std::cerr<<" [";
        for (std::size_t count=0;count<histogram.size();++count) {
            std::cerr<<count<<':'<<histogram[count]<<' ';
        }
        std::cerr<<']';
    }
    std::cerr<<'\n';
}

void checkLocalStrip(const BoundaryLayerBuildResult2D& layers,
                     const std::string& label,bool requireZeroLayer=true) {
    if (!layers.success()) {
        std::cerr<<label<<" layer failure: "
                 <<boundaryLayerFailureReasonName(layers.failure.reason)<<' '
                 <<layers.failure.message<<'\n';
    }
    check(layers.success(),label+" locally reduced strip succeeds");
    if (!layers.success()) return;
    check(layers.localReductionApplied,label+" records local reduction");
    std::size_t retained=0U,zero=0U,terminations=0U;
    for (const auto& strip:layers.strips) {
        retained+=strip.cells.size();
        zero+=strip.metrics.zeroLayerColumnCount;
        terminations+=strip.metrics.terminationEdgeCount;
        check(BoundaryLoop(strip.outerEnvelope()).diagnose().valid(),
              label+" has a valid stepped outer envelope");
        for (std::size_t segment=0;segment<strip.actualLayerCounts.size();++segment) {
            const auto next=(segment+1U)%strip.actualLayerCounts.size();
            const auto lhs=strip.actualLayerCounts[segment];
            const auto rhs=strip.actualLayerCounts[next];
            check(lhs<=rhs+1U && rhs<=lhs+1U,
                  label+" adjacent columns differ by at most one layer");
        }
    }
    check(retained>0U,label+" retains boundary-layer cells");
    if (requireZeroLayer) {
        check(zero>0U,label+" contains a true 0-layer patch");
    }
    check(terminations>0U,label+" contains real layer termination edges");
}

void checkHybrid(const HybridMeshBuildResult2D& hybrid,
                 const std::string& label,bool requireZeroLayer=true) {
    if (!hybrid.success()) {
        std::cerr<<label<<" hybrid failure: "
                 <<hybridMeshFailureReasonName(hybrid.failure.reason)<<' '
                 <<hybrid.failure.message<<'\n';
        const auto& topology=hybrid.solverTopology;
        const std::size_t count=std::min<std::size_t>(
            hybrid.solverQuality.issues.size(),12U);
        for (std::size_t issueId=0;issueId<count;++issueId) {
            const auto& issue=hybrid.solverQuality.issues[issueId];
            std::cerr<<"  issue "<<issueId<<" code="
                     <<static_cast<int>(issue.code)<<" cell="<<issue.cellId
                     <<" edge="<<issue.edgeId<<" measured="<<issue.measured;
            if (issue.cellId<topology.cells.size()) {
                std::cerr<<" polygon=";
                for (const auto vertex:topology.cells[issue.cellId].vertices) {
                    const auto& point=topology.vertices[vertex].point;
                    std::cerr<<'('<<point.x<<','<<point.y<<')';
                }
            }
            if (issue.edgeId<topology.edges.size()) {
                const auto& edge=topology.edges[issue.edgeId];
                std::cerr<<" face=("<<topology.vertices[edge.v0].point.x<<','
                         <<topology.vertices[edge.v0].point.y<<")-("
                         <<topology.vertices[edge.v1].point.x<<','
                         <<topology.vertices[edge.v1].point.y<<')';
                if (edge.neighbour) std::cerr<<" neighbour="<<*edge.neighbour;
                if (edge.neighbour && *edge.neighbour<topology.cells.size()) {
                    std::cerr<<" neighbour_polygon=";
                    for (const auto vertex:topology.cells[*edge.neighbour].vertices) {
                        const auto& point=topology.vertices[vertex].point;
                        std::cerr<<'('<<point.x<<','<<point.y<<')';
                    }
                }
            }
            std::cerr<<'\n';
        }
    }
    check(hybrid.success(),label+" commits a hybrid mesh instead of global fallback");
    if (!hybrid.success()) return;
    check(hybrid.metrics.boundaryLayerCellCount>0U &&
          hybrid.metrics.boundaryLayerCellCount<
              hybrid.metrics.requestedBoundaryLayerCellCount,
          label+" reports retained and locally removed layer cells");
    check((!requireZeroLayer || hybrid.metrics.zeroLayerColumnCount>0U) &&
          hybrid.metrics.terminationCellCount>0U &&
          hybrid.metrics.terminationEdgeCount>0U,
          label+" reports 0-layer columns and topology termination cells");
    check(hybrid.topology.valid() && hybrid.solverTopology.valid(),
          label+" base and solver topology pass");
    check(hybrid.interfaceAudit.pass(1.0e-8) &&
          hybrid.solverInterfaceAudit.pass(1.0e-8),
          label+" geometric and solver termination interfaces are conformal");
    check(std::abs(hybrid.metrics.areaError)<1.0e-8,
          label+" conserves domain-minus-solid area");
    check(hybrid.solverQuality.valid(),label+" passes unchanged solver quality");
}

} // namespace

int main() {
    const BoundaryLoop concaveWall({
        {0.0,0.0},{3.0,0.0},{3.0,3.0},{2.0,3.0},
        {2.0,1.0},{1.0,1.0},{1.0,3.0},{0.0,3.0}});
    const auto concaveChain=makeClosedWallChain2D(concaveWall,0U,"wall_0");
    const auto concaveLayers=concaveChain.success()
        ?buildLocallyReducedBoundaryLayerStrip2D(*concaveChain.chain,fourLayers())
        :BoundaryLayerBuildResult2D{};
    checkLocalStrip(concaveLayers,"concave L");
    if (concaveLayers.success()) printCounts(concaveLayers,"concave L");
    const auto concaveHybrid=buildConformalHybridMesh2D(
        concaveLayers,Domain2D{{{-1.0,-1.0},{4.0,4.0}}},
        BoundaryRegion2D(concaveWall),8U,refinement(8U));
    checkHybrid(concaveHybrid,"concave L");

    const BoundaryLoop sharpWall({
        {-1.0,0.0},{-0.5,-0.5},{0.5,-0.4},{1.5,-0.25},
        {2.3,-0.10},{2.7,-0.035},{3.0,0.0},
        {2.7,0.035},{2.3,0.10},{1.5,0.25},{0.5,0.4},{-0.5,0.5}});
    const auto sharpChain=makeClosedWallChain2D(sharpWall,0U,"wall_0");
    const auto sharpLayers=sharpChain.success()
        ?buildLocallyReducedBoundaryLayerStrip2D(*sharpChain.chain,fourLayers())
        :BoundaryLayerBuildResult2D{};
    checkLocalStrip(sharpLayers,"sharp trailing edge");
    if (sharpLayers.success()) {
        printCounts(sharpLayers,"sharp trailing edge");
        std::vector<std::size_t> unique=sharpLayers.strips.front().actualLayerCounts;
        std::sort(unique.begin(),unique.end());
        unique.erase(std::unique(unique.begin(),unique.end()),unique.end());
        check(unique==std::vector<std::size_t>({0U,1U,2U,3U,4U}),
              "sharp taper contains a real 4-to-3-to-2-to-1-to-0 sequence");
    }
    const auto sharpHybrid=buildConformalHybridMesh2D(
        sharpLayers,Domain2D{{{-2.0,-1.5},{4.0,1.5}}},
        BoundaryRegion2D(sharpWall),8U,refinement(8U));
    checkHybrid(sharpHybrid,"sharp trailing edge");
    const auto sharpRepeat=buildConformalHybridMesh2D(
        sharpLayers,Domain2D{{{-2.0,-1.5},{4.0,1.5}}},
        BoundaryRegion2D(sharpWall),8U,refinement(8U));
    check(sharpRepeat.success() &&
          sameTopology(sharpHybrid.topology,sharpRepeat.topology) &&
          sameTopology(sharpHybrid.solverTopology,sharpRepeat.solverTopology),
          "sharp termination connectivity/owner/neighbour is deterministic");

    const BoundaryLoop leftWall(subdividedRectangle(0.0,1.0,0.0,1.0,4U));
    const BoundaryLoop rightWall(subdividedRectangle(1.08,2.08,0.0,1.0,4U));
    const auto leftChain=makeClosedWallChain2D(leftWall,0U,"wall_0");
    const auto rightChain=makeClosedWallChain2D(rightWall,1U,"wall_1");
    const auto gapLayers=leftChain.success() && rightChain.success()
        ?buildLocallyReducedBoundaryLayerStrips2D(
             {*leftChain.chain,*rightChain.chain},fourLayers())
        :BoundaryLayerBuildResult2D{};
    checkLocalStrip(gapLayers,"narrow gap",false);
    if (gapLayers.success()) printCounts(gapLayers,"narrow gap");
    const auto gapHybrid=buildConformalHybridMesh2D(
        gapLayers,Domain2D{{{-1.0,-1.0},{3.08,2.0}}},
        BoundaryRegion2D({leftWall,rightWall}),8U,refinement(8U));
    checkHybrid(gapHybrid,"narrow gap",false);

    // Q4-1 construction selection must never cost the whole boundary layer.
    // On this taper the Q4 candidate fails a final gate, so the committed mesh
    // must be the unchanged non-Q4 hybrid rather than pure Cut-cell fallback.
    HybridMeshPolicy2D q4Policy;
    q4Policy.enableTerminationConstructionQualitySelection=true;
    const auto sharpQ4=sharpChain.success()
        ?buildRobustH4Mesh2D(
             {*sharpChain.chain},fourLayers(),Domain2D{{{-2.0,-1.5},{4.0,1.5}}},
             BoundaryRegion2D(sharpWall),8U,refinement(8U),{},q4Policy)
        :RobustH4BuildResult2D{};
    check(sharpQ4.success() && sharpQ4.mode==H4MeshMode2D::Hybrid,
          "sharp taper keeps a hybrid mesh when Q4-1 selection is requested");
    check(sharpQ4.hybridCandidate.metrics.q41ConstructionSelectionDeclined,
          "declined Q4-1 selection is reported instead of silently succeeding");
    check(sharpQ4.hybridCandidate.metrics.boundaryLayerCellCount==
              sharpHybrid.metrics.boundaryLayerCellCount,
          "declined Q4-1 selection retains every non-Q4 boundary-layer cell");

    // Q5-2 re-resolves the graded termination buffer with more, thinner rows.
    // The march is locally reduced, so that moves the stepped front and changes
    // the remainder: the outcome must be measured, and the re-resolved front may
    // only be committed when it strictly lowers the typed hard count.
    HybridMeshPolicy2D q52Policy;
    q52Policy.enableTerminationBufferRadialMatching=true;
    const auto sharpQ52=sharpChain.success()
        ?buildRobustH4Mesh2D(
             {*sharpChain.chain},fourLayers(),Domain2D{{{-2.0,-1.5},{4.0,1.5}}},
             BoundaryRegion2D(sharpWall),8U,refinement(8U),{},q52Policy)
        :RobustH4BuildResult2D{};
    check(sharpQ52.success() && sharpQ52.mode==H4MeshMode2D::Hybrid,
          "sharp taper keeps a hybrid mesh under a Q5-2 policy");
    if (sharpQ52.mode==H4MeshMode2D::Hybrid) {
        const auto& m=sharpQ52.hybridCandidate.metrics;
        check(m.q52TerminationBufferRadialCommitted==
                  (m.q52TerminationBufferHardWithMatching<
                   m.q52TerminationBufferHardWithHistoricalMarch),
              "Q5-2 commits the re-resolved buffer only on a strict hard-count win");
        check(m.q52TerminationBufferHardWithHistoricalMarch>0U,
              "Q5-2 reports the historical march hard count it was measured against");
        check(m.q52TerminationBufferRadialCommitted,
              "the sharp taper is a measured Q5-2 win");
    }

    // A bounded row/subdivision search may run out before its own thickness
    // target is reachable. Neither rule may then report a matched plan: doing so
    // would publish a cap the committed mesh does not satisfy.
    for (const double target:{1.0,0.1,0.01,0.001}) {
        const auto tight=resolveAutomaticHybridTransitionPlan2D(
            sharpLayers,Domain2D{{{-2.0,-1.5},{4.0,1.5}}},refinement(8U),target,4U);
        check(tight.has_value(),
              "Q5-1 bounded-search plan resolves for every scanned target");
        if (!tight) continue;
        const double cap=target*tight->targetCellSize;
        const double row=tight->ringThickness/
            static_cast<double>(tight->outerRingRadialSubdivision);
        check(tight->outerRingRadialTargetReachable
                  ?row<=cap*(1.0+1.0e-12)
                  :tight->outerRingRadialSubdivision==1U,
              "Q5-1 reports a matched outer row only when it really meets the cap");
    }
    HybridMeshPolicy2D q52Tight;
    q52Tight.enableTerminationBufferRadialMatching=true;
    q52Tight.outerTransitionRadialTargetCells=0.01;
    const auto sharpTight=sharpChain.success()
        ?buildRobustH4Mesh2D(
             {*sharpChain.chain},fourLayers(),Domain2D{{{-2.0,-1.5},{4.0,1.5}}},
             BoundaryRegion2D(sharpWall),8U,refinement(8U),{},q52Tight)
        :RobustH4BuildResult2D{};
    check(sharpTight.success() && sharpTight.mode==H4MeshMode2D::Hybrid,
          "an unreachable Q5-2 cap still leaves a hybrid mesh");
    if (sharpTight.mode==H4MeshMode2D::Hybrid) {
        const auto& m=sharpTight.hybridCandidate.metrics;
        check(!m.q52TerminationBufferRowCapReachable,
              "Q5-2 reports an unreachable row cap instead of pretending to match");
        check(!m.q52TerminationBufferRadialCommitted,
              "Q5-2 commits nothing when its row cap is unreachable");
        check(m.q52TerminationBufferRowCap<=0.0 ||
              m.q52TerminationBufferOuterRowThickness<=
                  m.q52TerminationBufferRowCap*(1.0+1.0e-12) ||
              !m.q52TerminationBufferRadialCommitted,
              "a committed Q5-2 buffer never exceeds the cap it reports");
    }

    const BoundaryLoop fallbackWall({{0.3,0.3},{0.7,0.3},{0.7,0.7},{0.3,0.7}});    const auto fallbackChain=makeClosedWallChain2D(fallbackWall,0U,"wall_0");
    LayerParameters2D invalidLayers=fourLayers();
    invalidLayers.nLayers=0U;
    const auto fallback=fallbackChain.success()
        ?buildRobustH4Mesh2D(
             {*fallbackChain.chain},invalidLayers,Domain2D{{{0.0,0.0},{1.0,1.0}}},
             BoundaryRegion2D(fallbackWall),6U,refinement(6U))
        :RobustH4BuildResult2D{};
    check(fallback.success() &&
          fallback.mode==H4MeshMode2D::PureCutCellFallback,
          "pure Cut-cell remains an explicit final fallback");
    check(!fallback.requestedLayerCandidate.success() &&
          !fallback.localLayerCandidate.success() &&
          fallback.fallbackStage==H4FallbackStage2D::LocalReduction,
          "fallback report proves requested and local layer attempts failed first");

    if (failures==0) {
        std::cout<<"cartmesh2d H4-3 local termination tests: PASS\n";
        return EXIT_SUCCESS;
    }
    std::cerr<<failures<<" H4-3 local termination test(s) failed\n";
    return EXIT_FAILURE;
}