#include "cartmesh2d/topology/SharedEdgePartition2D.hpp"
#include "cartmesh2d/topology/Topology2D.hpp"
#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace cartmesh2d;
static_assert(requires(const BoundaryLayerBuildResult2D& layers,const Domain2D& domain,
                       const BoundaryRegion2D& walls,const QuadtreeRefinementPolicy2D& refinement) {
    buildConformalHybridMesh2D(layers,domain,walls,4,refinement);
    buildConformalHybridMesh2D(layers,domain,walls,4,refinement,{});
});
namespace {
int failures=0;
void check(bool ok,const char* message) {
    if (!ok) {++failures;std::cerr<<"FAIL: "<<message<<'\n';}
}
template<class F> void rejects(F f,const char* message) {
    bool threw=false;try {f();} catch (const std::exception&) {threw=true;}
    check(threw,message);
}
CutCell2D rectangle(std::size_t id,AABB2D box) {
    CutCell2D c;c.sourceId=id;c.sourceKey=id;c.backgroundBounds=box;
    c.kind=CutCellKind::Full;c.area=(box.max.x-box.min.x)*(box.max.y-box.min.y);
    c.areaFraction=1;c.fluidPolygon={{{box.min.x,box.min.y},{box.max.x,box.min.y},
                                    {box.max.x,box.max.y},{box.min.x,box.max.y}}};
    c.centroid=c.fluidPolygon.centroid();return c;
}
bool identical(const TopologyMesh2D& a,const TopologyMesh2D& b) {
    if (!a.valid() || !b.valid() || a.vertices.size()!=b.vertices.size() ||
        a.cells.size()!=b.cells.size() || a.edges.size()!=b.edges.size()) return false;
    for (std::size_t i=0;i<a.vertices.size();++i)
        if (a.vertices[i].point.x!=b.vertices[i].point.x || a.vertices[i].point.y!=b.vertices[i].point.y) return false;
    for (std::size_t i=0;i<a.cells.size();++i)
        if (a.cells[i].vertices!=b.cells[i].vertices || a.cells[i].edges!=b.cells[i].edges ||
            a.cells[i].geometryArea!=b.cells[i].geometryArea) return false;
    for (std::size_t i=0;i<a.edges.size();++i) {
        const auto& x=a.edges[i];const auto& y=b.edges[i];
        if (x.v0!=y.v0 || x.v1!=y.v1 || x.owner!=y.owner || x.neighbour!=y.neighbour || x.patch!=y.patch) return false;
    }
    return true;
}
}

int main(int argc,char** argv) {
    // Optional isolated topology benchmark; generated geometry is identical
    if (argc==2) {
        const auto n=static_cast<std::size_t>(std::stoul(argv[1]));
        const Domain2D domain{{{0,0},{1,1}}};
        const BoundaryRegion2D wall(BoundaryLoop({{2,2},{3,2},{3,3},{2,3}}));
        std::vector<CutCell2D> cells;
        for (std::size_t y=0;y<n;++y) for (std::size_t x=0;x<n;++x)
            cells.push_back(rectangle(y*n+x,{{double(x)/double(n),double(y)/double(n)},
                                           {double(x+1)/double(n),double(y+1)/double(n)}}));
        const auto start=std::chrono::steady_clock::now();
        const auto legacy=buildGlobalTopology(cells,domain,wall);
        const auto middle=std::chrono::steady_clock::now();
        const auto shared=buildGlobalTopology(cells,domain,wall,{},std::make_shared<IntersectionRegistry2D>());
        const auto end=std::chrono::steady_clock::now();
        check(identical(legacy,shared),"benchmark topology and patch identity unchanged");
        std::cout<<"{\"cells\":"<<cells.size()<<",\"legacy_seconds\":"
            <<std::chrono::duration<double>(middle-start).count()<<",\"shared_seconds\":"
            <<std::chrono::duration<double>(end-middle).count()<<",\"partition_cache_hits\":"
            <<shared.sharedPartitionCacheHits<<",\"identical\":"<<(failures?"false":"true")<<"}\n";
        return failures?1:0;
    }
    for (double scale:{1.e-6,1.,1.e6}) {
        IntersectionRegistry2D registry;
        registry.configureGrid({{0,0},{scale,scale}},4);
        const Segment2D s{{.125*scale,.25*scale},{.875*scale,.75*scale}};
        const auto support=registry.registerSegment(s,scale/4,IntersectionSource2D::WallCartesian);
        check(registry.registerSegment({s.b,s.a},scale/8,IntersectionSource2D::WallCartesian)==support,
              "reversed support has one identity");
        const auto x=registry.gridLine(0,.5*scale);
        const auto id=registry.intersectGridLine(support,x,scale/4);
        check(id==registry.intersectGridLine(support,x,scale/8),"coarse/fine calls reuse one intersection");
        const auto y=registry.gridLine(1,.5*scale);
        check(id==registry.intersectGridLine(support,y,scale/8),"exact grid corner shares one vertex");
        const auto& exactStore=registry.shadowVertexStore();
        check(exactStore.resolveExactKey(
                  {StableVertexKeyKind2D::SourceVertex,
                   static_cast<std::uint64_t>(support),0U,0U,0U}) ==
                  static_cast<StableVertexId2D>(registry.internVertex(
                      s.a,scale/8,IntersectionFeature2D::Smooth)),
              "source endpoint has a typed exact alias");
        check(exactStore.resolveExactKey(
                  {StableVertexKeyKind2D::WallGridIntersection,
                   static_cast<std::uint64_t>(support),0U,8U,0U}) == id &&
              exactStore.resolveExactKey(
                  {StableVertexKeyKind2D::GridVertex,0U,0U,8U,8U}) == id,
              "wall-grid event and grid corner resolve to one stable id");
        check(registry.events().size()==2 && registry.intersectionCacheHits()==1,"one solve per support/grid line");
        check(registry.events()[0].localH==scale/8,"event records smallest incident local_h");
        const auto sharp=registry.internVertex(s.a,scale/8,IntersectionFeature2D::WallSharpCorner);
        const auto endpoint=registry.intersectGridLine(support,registry.gridLine(0,s.a.x),scale/8);
        check(endpoint==sharp && registry.vertices()[endpoint].feature==IntersectionFeature2D::WallSharpCorner,
              "incident sharp endpoint is immutable and keeps classification");
        const auto gap=registry.registerSegment({{s.a.x,s.a.y+scale*1.e-8},{s.b.x,s.b.y+scale*1.e-8}},
                                                scale/8,IntersectionSource2D::WallCartesian);
        check(registry.intersectGridLine(gap,x,scale/8)!=id,"nearby nonincident wall support cannot snap across gap");
        const auto report=intersectionConstructionToJson(registry,{id},0,0);
        check(report==intersectionConstructionToJson(registry,{id},0,0) &&
              report.find("source_segment")!=std::string::npos && report.find("wall_sharp_corner")!=std::string::npos,
              "deterministic provenance includes source geometry and feature");
        const auto& stable=registry.shadowVertexStore().records();
        check(stable.size()==registry.vertices().size() &&
              !stable[id].sourceRefs.empty(),
              "shared construction events retain stable source lineage");
        rejects([&]{(void)registry.intersectGridLine(support,x,0);},"invalid local_h rejected on cached event");
        rejects([&]{(void)registry.intersectGridLine(support,registry.gridLine(0,0),scale);},"out-of-segment intersection rejected");
        rejects([&]{(void)registry.gridLine(0,.123*scale);},"non-dyadic coordinate is not silently rounded");
        rejects([&]{registry.configureGrid({{0,0},{scale,scale}},4);},"context cannot be reconfigured");
        const auto parallel=registry.registerSegment({{0,scale/4},{scale,scale/4}},scale,
                                                     IntersectionSource2D::WallCartesian);
        rejects([&]{(void)registry.intersectGridLine(parallel,registry.gridLine(1,scale/4),scale);},
                "collinear support is not an isolated intersection");
    }
    {
        IntersectionRegistry2D r;
        rejects([&]{r.configureGrid({{0,0},{std::numeric_limits<double>::infinity(),1}},4);},"nonfinite grid rejected");
        rejects([&]{(void)r.internVertex({0,0},-1);},"invalid vertex scale rejected");
    }
    {
        IntersectionRegistry2D r;r.configureGrid({{0,0},{1,1}},4);
        const auto a=r.registerSegment({{.1,.5},{.9,.5}},.25,IntersectionSource2D::WallCartesian);
        const auto b=r.registerSegment({{.1,.5+2.e-16},{.9,.5+2.e-16}},.25,IntersectionSource2D::WallCartesian);
        (void)r.intersectGridLine(a,r.gridLine(0,.5),.25);
        rejects([&]{(void)r.intersectGridLine(b,r.gridLine(0,.5),.25);},
                "even sub-roundoff nonincident walls must not be silently welded");
    }
    {
        IntersectionRegistry2D r;r.configureGrid({{0,0},{1,1}},4);
        (void)r.internVertex({.5,.5},.25,IntersectionFeature2D::WallSharpCorner);
        const auto support=r.registerSegment({{.1,.5+2.e-16},{.9,.5+2.e-16}},.25,
                                             IntersectionSource2D::WallCartesian);
        rejects([&]{(void)r.intersectGridLine(support,r.gridLine(0,.5),.25);},
                "pre-registered nonincident sharp feature is protected even without prior events");
        const auto small=r.registerSegment({{.5-4.e-16,.4},{.5+4.e-16,.4}},1,
                                           IntersectionSource2D::WallCartesian);
        const auto a=r.internVertex({.5-4.e-16,.4},1), b=r.internVertex({.5+4.e-16,.4},1);
        const auto middle=r.intersectGridLine(small,r.gridLine(0,.5),1);
        check(middle!=a && middle!=b,"coarse local_h cannot collapse a resolved short support");
    }
    {
        IntersectionRegistry2D r;r.configureGrid({{0,0},{1,1}},4);
        const auto support=r.registerSegment({{.1,.5+2.e-16},{.9,.5+2.e-16}},.25,
                                             IntersectionSource2D::WallCartesian);
        (void)r.intersectGridLine(support,r.gridLine(0,.5),.25);
        rejects([&]{(void)r.internVertex({.5,.5},.25,IntersectionFeature2D::WallSharpCorner);},
                "late feature registration cannot silently legitimize a nonincident snap");
    }
    {
        const Domain2D domain{{{0,0},{2,1}}};
        const BoundaryRegion2D wall(BoundaryLoop({{3,2},{4,2},{4,3},{3,3}}));
        auto cells=std::vector{rectangle(0,{{0,0},{1,1}}),rectangle(1,{{1,0},{2,.5}}),rectangle(2,{{1,.5},{2,1}})};
        auto registry=std::make_shared<IntersectionRegistry2D>();
        const auto original=buildGlobalTopology(cells,domain,wall);
        const auto shared=buildGlobalTopology(cells,domain,wall,{},registry);
        check(identical(original,shared),"coarse/fine partition, area, ownership and patches match legacy");
        check(shared.sharedPartitionCacheHits>0 && shared.canonicalVertexIds.size()==shared.vertices.size(),
              "actual topology consumes common edge partitions and retains canonical handles");
        (void)registry->internVertex({1,.25},.5); // rejected candidate
        std::reverse(cells.begin(),cells.end());
        const auto again=buildGlobalTopology(cells,domain,wall,{},registry);
        check(identical(shared,again) && shared.canonicalVertexIds==again.canonicalVertexIds,
              "reversed input and unused candidate do not alter topology or handles");
        auto bad=cells;
        bad[0].canonicalVertexIds={999};bad[0].canonicalRegistry=registry.get();
        check(!buildGlobalTopology(bad,domain,wall,{},registry).valid(),"stale handles fail explicitly");
        bad=cells;bad[0].canonicalVertexIds={0,1,2,3};
        check(!buildGlobalTopology(bad,domain,wall,{},registry).valid(),"foreign context fails explicitly");
    }
    {
        const Domain2D domain{{{0,0},{1,1}}};
        const BoundaryRegion2D wall(BoundaryLoop({{.25,.2},{.75,.35},{.75,.7},{.25,.7}}));
        Quadtree2D tree(domain,3,wall);QuadtreeRefinementPolicy2D policy;policy.minimumLevel=2;policy.boundaryLevel=3;
        tree.refine(wall,policy);(void)tree.enforceTwoToOneBalance(wall);
        auto registry=std::make_shared<IntersectionRegistry2D>();registry->configureGrid(domain.bounds,3);
        std::vector<CutCell2D> cells;
        for (const auto& leaf:tree.leaves()) {
            auto cut=buildCutCellsShared(leaf,wall,*registry);
            for (const auto& c:cut) check(c.valid() && c.canonicalVertexIds.size()==c.fluidPolygon.vertices.size(),
                                         "real Cut-cell returns registry handles");
            cells.insert(cells.end(),cut.begin(),cut.end());
        }
        const auto mesh=buildGlobalTopology(cells,domain,wall,{},registry);
        check(mesh.valid() && registry->intersectionCacheHits()>0,"adjacent real cuts reuse intersections and form conformal mesh");
        double area=0;for (const auto& c:mesh.cells) area+=c.geometryArea;
        check(std::abs(area-(1.-Polygon2D{wall.loops()[0].vertices()}.area()))<1.e-12,"exterior fluid area conserved without deleting cells");
        auto foreign=std::make_shared<IntersectionRegistry2D>();
        check(!buildGlobalTopology(cells,domain,wall,{},foreign).valid(),"Cut-cell handles cannot cross contexts");
    }
    return failures?1:0;
}
