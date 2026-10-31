#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <iomanip>

using namespace cartmesh2d;
namespace {
int failures=0;
void check(bool pass,const char* message) {
    if (!pass) {++failures;std::cerr<<message<<'\n';}
}
CutCell2D cell(std::size_t id,std::vector<Point2D> points) {
    CutCell2D c;c.sourceId=id;c.sourceKey=id;c.sourceLineage={id};
    c.refinementLineageKeys={100+id};c.kind=CutCellKind::Cut;
    c.fluidPolygon={std::move(points)};c.area=c.fluidPolygon.area();
    c.areaFraction=1;c.centroid=c.fluidPolygon.centroid();
    c.backgroundBounds=c.fluidPolygon.bounds();return c;
}
double area(const TopologyMesh2D& mesh) {
    double sum=0;for (const auto& c:mesh.cells) sum+=c.geometryArea;return sum;
}
std::string signature(const TopologyMesh2D& mesh) {
    std::ostringstream out;out<<std::setprecision(17);
    for (const auto& v:mesh.vertices) out<<v.point.x<<' '<<v.point.y<<'\n';
    for (const auto& c:mesh.cells) {
        out<<c.sourceKey<<' '<<c.geometryArea;
        for (const auto v:c.vertices) out<<' '<<v;
        out<<'\n';
    }
    for (const auto& e:mesh.edges) out<<e.v0<<' '<<e.v1<<' '<<e.owner<<' '
        <<e.neighbour.value_or(std::numeric_limits<std::size_t>::max())<<' '<<int(e.patch)<<'\n';
    return out.str();
}
void sliver(double scale,double rotation) {
    const auto transform=[&](Point2D p) {
        return Point2D{scale*(std::cos(rotation)*p.x-std::sin(rotation)*p.y),
            scale*(std::sin(rotation)*p.x+std::cos(rotation)*p.y)};
    };
    std::vector<Point2D> outline{{0,0},{2,0},{2,1},{0,1}};
    for (auto& p:outline) p=transform(p);
    BoundaryRegion2D boundary{BoundaryLoop(outline)};
    const Domain2D domain{{{-3*scale,-3*scale},{3*scale,3*scale}}};
    std::vector<CutCell2D> cells{
        cell(10,{{0,0},{1,0},{1,.01}}),
        cell(11,{{0,0},{1,.01},{1,1},{0,1}}),
        cell(12,{{1,0},{2,0},{2,1},{1,1}})};
    for (auto& c:cells) {
        auto points=c.fluidPolygon.vertices;for (auto& p:points) p=transform(p);
        c=cell(c.sourceId,std::move(points));
    }
    const auto raw=buildGlobalTopology(cells,domain,boundary);
    const auto result=agglomerateSolverByShape2D(raw,domain,boundary,{false,false,true});
    check(raw.valid() && !result.before.valid(),"sliver must trigger shape repair");
    check(result.targetMet() && result.mergedCellCount==1 && result.topology.cells.size()==2,
        "sliver exact union must meet the target without adding diagonals");
    check(std::abs(area(result.topology)-area(raw))<=TolerancePolicy{}.relative*area(raw),
        "union preserves real fluid area under scaling and rotation");
    bool locked=false,lineage=false;
    for (std::size_t i=0;i<result.topology.cells.size();++i) {
        const auto& c=result.topology.cells[i];
        if (c.sourceId==12) locked=result.immutableCells[i] && c.sourceKey==12;
        if (c.sourceLineage==std::vector<std::size_t>{10,11})
            lineage=c.refinementLineageKeys==std::vector<std::uint64_t>{110,111};
    }
    check(locked && lineage,"immutable identity and merged source/refinement provenance survive");
    std::set<std::pair<double,double>> vertices;
    for (const auto& v:result.topology.vertices) vertices.emplace(v.point.x,v.point.y);
    for (const auto& p:outline) check(vertices.contains({p.x,p.y}),"every input physical corner is retained");
    const auto again=agglomerateSolverByShape2D(raw,domain,boundary,{false,false,true});
    check(signature(result.topology)==signature(again.topology),
        "same input has deterministic topology and serialization");
    const auto stopped=agglomerateSolverByShape2D(raw,domain,boundary,{true,true,true});
    check(stopped.valid() && !stopped.targetMet() && stopped.stopReason=="no-improving-union" &&
        stopped.topology.cells.size()==raw.cells.size(),"unrepairable immutable sliver is explicit, never a pass");
    SolverShapeAgglomerationPolicy2D smallBudget;smallBudget.maximumAreaRatio=1;
    const auto bounded=agglomerateSolverByShape2D(raw,domain,boundary,{},smallBudget);
    double largestBox=0;
    for (const auto& c:cells) {
        const auto b=c.fluidPolygon.bounds();
        largestBox=std::max(largestBox,(b.max.x-b.min.x)*(b.max.y-b.min.y));
    }
    bool withinBudget=bounded.valid();
    for (const auto& c:bounded.topology.cells)
        withinBudget=withinBudget && c.geometryArea<=largestBox*(1+TolerancePolicy{}.relative);
    check(withinBudget,"resolution budget cannot be bypassed by repeated unions");
    SolverShapeAgglomerationPolicy2D invalid;invalid.target.minInteriorAngleDeg=std::numeric_limits<double>::quiet_NaN();
    check(!agglomerateSolverByShape2D(raw,domain,boundary,{},invalid).valid(),
        "nonfinite target is rejected rather than comparing as a pass");
}
}
int main() {
    for (const double s:{.001,1.,1000.}) for (const double a:{0.,.31}) sliver(s,a);
    const auto c=cell(0,{{0,0},{1,0},{1,1},{.5,.95},{0,1}});
    const BoundaryRegion2D wall{BoundaryLoop(c.fluidPolygon.vertices)};
    const Domain2D domain{{{-1,-1},{2,2}}};
    const auto raw=buildGlobalTopology({c},domain,wall);
    const auto result=agglomerateSolverByShape2D(raw,domain,wall);
    check(result.targetMet() && result.topology.cells.size()==1 && result.after.maxConcavityDeg>0,
        "admissible physical reflex wall is retained as one cell");
    const auto bad=cell(0,{{0,0},{2,0},{2,2},{1.9,2},{1.9,.1},{0,.1}});
    const BoundaryRegion2D badWall{BoundaryLoop(bad.fluidPolygon.vertices)};
    const Domain2D badDomain{{{-1,-1},{3,3}}};
    const auto rejected=agglomerateSolverByShape2D(buildGlobalTopology({bad},badDomain,badWall),badDomain,badWall);
    check(rejected.valid() && !rejected.targetMet(),"centroid outside a face half-plane remains rejected");
    return failures?1:0;
}
