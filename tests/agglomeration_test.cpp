#include "cartmesh2d/stabilization/SmallCell2D.hpp"

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
    if (!condition) { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

CutCell2D polygonCell(std::size_t id, std::uint64_t key,
                      const std::vector<Point2D>& polygon,
                      const AABB2D& background, CutCellKind kind) {
    CutCell2D cell;
    cell.sourceId=id; cell.sourceKey=key; cell.backgroundBounds=background; cell.kind=kind;
    cell.fluidPolygon={polygon}; cell.area=cell.fluidPolygon.area();
    const double bg=(background.max.x-background.min.x)*(background.max.y-background.min.y);
    cell.areaFraction=cell.area/bg; cell.centroid=cell.fluidPolygon.centroid();
    return cell;
}
}

int main() {
    const Domain2D domain{{{0.0,0.0},{2.0,1.0}}};
    BoundaryLoop boundary({{0.9,0.0},{2.0,0.0},{2.0,1.0},{0.9,1.0}});
    std::vector<CutCell2D> cells;
    cells.push_back(polygonCell(10,10,{{0.9,0.0},{1.0,0.0},{1.0,1.0},{0.9,1.0}},
                                {{0.0,0.0},{1.0,1.0}},CutCellKind::Cut));
    cells.push_back(polygonCell(20,20,{{1.0,0.0},{2.0,0.0},{2.0,1.0},{1.0,1.0}},
                                {{1.0,0.0},{2.0,1.0}},CutCellKind::Full));
    const auto topo=buildGlobalTopology(cells,domain,boundary);
    SmallCellPolicy2D policy; policy.areaFractionThreshold=0.2;
    const auto analysis=analyzeSmallCells(cells,topo,policy);
    const auto merged=agglomerateSmallCells(cells,topo,analysis,domain,boundary);
    check(merged.valid(), "manufactured sliver agglomeration succeeds");
    check(merged.inputCellCount==2 && merged.outputCellCount==1,
          "one small cell is absorbed into one stable cell");
    check(merged.mergedSmallCellCount==1, "exactly one small cell is merged");
    check(merged.cells.size()==1 && merged.cells[0].memberTopologyCellIds.size()==2,
          "agglomerated cell retains both source members");
    check(std::abs(merged.totalAreaBefore-1.1)<=1e-12 &&
          std::abs(merged.totalAreaAfter-1.1)<=1e-12 && merged.areaError<=1e-12,
          "manufactured merge conserves area");
    check(merged.topology.valid(), "rebuilt topology passes Stage 2D-4 audit");
    if (!merged.cells.empty()) {
        check(std::abs(merged.cells[0].area-1.1)<=1e-12, "merged polygon has exact union area");
        check(merged.cells[0].polygon.vertices.size()==4,
              "collinear shared interface disappears from rectangular union boundary");
    }

    SmallCellPolicy2D lowPolicy; lowPolicy.areaFractionThreshold=0.05;
    const auto stableAnalysis=analyzeSmallCells(cells,topo,lowPolicy);
    const auto stable=agglomerateSmallCells(cells,topo,stableAnalysis,domain,boundary);
    check(stable.valid(), "identity agglomeration succeeds when no cells are small");
    check(stable.inputCellCount==2 && stable.outputCellCount==2 && stable.mergedSmallCellCount==0,
          "no-small case preserves cell count");
    check(stable.areaError<=1e-12, "no-small case preserves area");

    auto badAnalysis=stableAnalysis;
    badAnalysis.records[0].status=SmallCellStatus2D::CandidateFound;
    badAnalysis.records[0].targetTopologyCellId=0;
    const auto bad=agglomerateSmallCells(cells,topo,badAnalysis,domain,boundary);
    check(!bad.valid(), "self-targeting candidate is rejected explicitly");

    // Preserve the historical small-cell agglomeration fixture as an explicit
    // interior-flow regression. Product/CLI defaults are tested separately as
    // exterior fluid around a solid obstacle.
    std::vector<Point2D> circle;
    constexpr std::size_t segments=64;
    constexpr double cx=0.07, cy=0.03;
    for (std::size_t i=0;i<segments;++i) {
        const double a=2.0*std::numbers::pi*static_cast<double>(i)/static_cast<double>(segments);
        circle.push_back({cx+std::cos(a),cy+std::sin(a)});
    }
    BoundaryLoop circleBoundary(circle);
    const Domain2D circleDomain{{{-2.0,-2.0},{2.0,2.0}}};
    Quadtree2D tree(circleDomain,4,circleBoundary);
    QuadtreeRefinementPolicy2D refine; refine.boundaryLevel=4;
    tree.refine(circleBoundary,refine);
    const auto balance=tree.enforceTwoToOneBalance(circleBoundary);
    check(balance.violationsAfter==0,"shifted-circle tree is 2:1 balanced");
    std::vector<CutCell2D> circleCells;
    for (const auto& leaf:tree.leaves()) {
        circleCells.push_back(buildCutCell(leaf,circleBoundary,FluidRegion2D::Interior));
    }
    const auto circleTopo=buildGlobalTopology(circleCells,circleDomain,circleBoundary);
    check(circleTopo.valid(),"shifted-circle source topology valid");
    const auto circleAnalysis=analyzeSmallCells(circleCells,circleTopo);
    check(circleAnalysis.valid() && circleAnalysis.smallCellCount==8 && circleAnalysis.unresolvedCount==0,
          "shifted-circle baseline has eight resolvable small cells");
    const auto circleMerged=agglomerateSmallCells(circleCells,circleTopo,circleAnalysis,
                                                   circleDomain,circleBoundary);
    check(circleMerged.valid(),"shifted-circle agglomeration succeeds");
    check(circleMerged.mergedSmallCellCount==8,"all eight shifted-circle small cells are merged");
    check(circleMerged.outputCellCount+8==circleMerged.inputCellCount,
          "each resolved small cell reduces output cell count by one");
    check(circleMerged.areaError<=1e-10,"shifted-circle total fluid area is conserved");
    check(circleMerged.topology.audit.duplicateVertices==0 &&
          circleMerged.topology.audit.duplicateEdges==0 &&
          circleMerged.topology.audit.orphanInternalEdges==0 &&
          circleMerged.topology.audit.nonManifoldEdges==0 &&
          circleMerged.topology.audit.unclassifiedBoundaryEdges==0 &&
          circleMerged.topology.audit.openCellLoops==0 &&
          circleMerged.topology.audit.areaMismatches==0,
          "shifted-circle rebuilt topology has zero Stage 2D-4 audit errors");
    for (const auto& c:circleMerged.cells) {
        check(c.area>0.0 && c.polygon.vertices.size()>=3,
              "every stabilized cell has positive simple polygon geometry");
    }

    if (failures) { std::cerr<<failures<<" test(s) failed\n"; return EXIT_FAILURE; }
    std::cout<<"cartmesh2d 2D-5B agglomeration tests passed\n";
    return EXIT_SUCCESS;
}
