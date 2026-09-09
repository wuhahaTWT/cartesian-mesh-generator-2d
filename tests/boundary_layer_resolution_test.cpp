#include "cartmesh2d/sizing/BoundaryLayerResolution2D.hpp"
#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include <cmath>
#include <iostream>
using namespace cartmesh2d;
namespace {
int failures=0;
void check(bool pass,const char* message) {if (!pass) {++failures;std::cerr<<message<<'\n';}}
HybridMeshBuildResult2D fixture(double scale) {
    HybridMeshBuildResult2D hybrid;
    BoundaryLayerStrip2D strip;
    strip.wallChain.fluidSide=FluidSide2D::Left;
    strip.wallChain.segments={{{0,0},{scale,0}}};
    strip.parameters.nLayers=2;
    hybrid.strips={strip};
    auto& mesh=hybrid.solverTopology;
    const std::vector<Point2D> points{{0,0},{1,0},{1,.2},{0,.2},{1,.5},{0,.5}};
    for (const auto p:points) mesh.vertices.push_back({mesh.vertices.size(),{p.x*scale,p.y*scale}});
    for (std::size_t i=0;i<2;++i) {
        TopologyCell2D cell;cell.id=i;cell.sourceId=i;cell.sourceLineage={i};
        cell.vertices=i==0?std::vector<std::size_t>{0,1,2,3}:std::vector<std::size_t>{3,2,4,5};
        cell.geometryArea=(i==0?.2:.3)*scale*scale;mesh.cells.push_back(cell);
        HybridSourceCell2D source;source.id=i;source.kind=HybridCellKind2D::BoundaryLayer;
        source.stripId=0;source.wallSegment=0;source.layerIndex=i;
        for (const auto v:cell.vertices) source.polygon.vertices.push_back(mesh.vertices[v].point);
        source.area=cell.geometryArea;hybrid.sourceCells.push_back(source);
    }
    mesh.edges={{0,0,1,0,std::nullopt,BoundaryPatch2D::EmbeddedBoundary},
                {1,3,2,0,1,BoundaryPatch2D::None}};
    return hybrid;
}
}
int main() {
    for (double scale:{.001,1.,1000.}) {
        const auto original=makeClosedWallChain2D(BoundaryLoop({{0,0},{scale,0},{scale,.5*scale},{0,.5*scale}}),0,"wall");
        const auto refined=refineWallChainToSize2D(*original.chain,.125*scale);
        check(refined.success() && refined.chain->segmentCount()==24,
              "tangential subdivision count is scale invariant");
        for (const auto p:original.chain->vertices) {
            bool retained=false;
            for (const auto q:refined.chain->vertices) retained=retained || (p.x==q.x && p.y==q.y);
            check(retained,"tangential sizing preserves exact input corners");
        }
        for (const auto segment:refined.chain->segments)
            check(std::hypot(segment.b.x-segment.a.x,segment.b.y-segment.a.y)/scale<=.125+1.e-12,
                  "tangential segments satisfy the requested relative size");
        auto h=fixture(scale);
        const auto full=measureBoundaryLayerResolution2D(h.solverTopology,scale,2,&h);
        check(full.status=="evaluated" && full.retainedCells==2 && full.requestedCells==2,
              "final two-layer column must retain both source polygons");
        check(full.columns[0].retainedLayers==2 && std::abs(full.fullLayerWallLength/full.wallLength-1.)<1.e-12,
              "full coverage requires the complete connected column");
        check(full.columns[0].firstLayerNormalHeightMin &&
              std::abs(*full.columns[0].firstLayerNormalHeightMin-.2)<1.e-12,
              "normal front height uses the reference length");
        const auto partial=measureBoundaryLayerResolution2D(h.solverTopology,scale,3,&h);
        check(partial.requestedCells==3 && partial.retainedCells==2 && partial.fullLayerWallLength==0 &&
              partial.firstLayerWallLength==partial.wallLength,"local stopping differs from first-layer coverage");
    }
    auto h=fixture(1.);
    const auto coarse=measureBoundaryLayerResolution2D(h.solverTopology,1.,2,&h,{},.1);
    check(coarse.firstLayerHeightExceededWallLength==coarse.wallLength,
          "preserved layers can still exceed the requested first-layer height");
    h.solverTopology.vertices[4].point.y=.6;
    auto report=measureBoundaryLayerResolution2D(h.solverTopology,1.,2,&h);
    check(report.status=="incomplete" && report.sourceMismatches==1 && report.retainedCells==1,
          "source labels and stale cell areas cannot hide changed final vertices");
    h=fixture(1.);h.solverTopology.cells[0].sourceLineage={1};
    report=measureBoundaryLayerResolution2D(h.solverTopology,1.,2,&h);
    check(report.retainedCells==0 && report.sourceMismatches==2,"misassigned lineage cannot certify retained layers");
    h=fixture(1.);h.solverTopology.edges[1].neighbour.reset();
    report=measureBoundaryLayerResolution2D(h.solverTopology,1.,2,&h);
    check(report.status=="incomplete" && report.continuityMismatches==1 &&
          report.columns[0].retainedLayers==1 && report.fullLayerWallLength==0,
          "disconnected upper cell cannot count as complete wall coverage");
    h=fixture(1.);
    report=measureBoundaryLayerResolution2D(h.solverTopology,1.,2);
    check(report.status=="no_layers_retained" && report.retainedCells==0,
          "pure fallback explicitly reports zero retained layers");
    return failures?1:0;
}
