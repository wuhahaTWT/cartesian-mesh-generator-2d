#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;
void require(bool b,const char* text){if(!b)throw std::runtime_error(text);}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::runtime_error&){caught=true;}require(caught,"invalid grid accepted");}
int main(){try{
    const std::vector<double> x={0,.2,.5,1},y={0,.01,.05,.2,1};
    const auto t=makeRectilinearMesh2D(x,y),repeat=makeRectilinearMesh2D(x,y);
    const auto m=fv::makeFvMesh2D(t);
    require(m.cells.size()==12&&m.faces.size()==31&&t.vertices.size()==20,"wrong grid cardinality");
    double area=0;std::size_t interior=0;
    for(const auto& c:m.cells)area+=c.area;
    for(const auto& f:m.faces)if(f.neighbour)++interior;
    require(std::abs(area-1)<1e-14&&interior==17,"partition area/incidence differs");
    for(std::size_t i=0;i<t.edges.size();++i) {
        const auto& a=t.edges[i];const auto& b=repeat.edges[i];
        require(a.v0==b.v0&&a.v1==b.v1&&a.owner==b.owner&&a.neighbour==b.neighbour,"nondeterministic shared topology");
        require(std::abs(m.faces[i].correction.x)<1e-14&&std::abs(m.faces[i].correction.y)<1e-14,"rectilinear face is not orthogonal");
    }
    std::vector<double> field,bc;
    for(const auto& c:m.cells)field.push_back(2*c.centre.x-3*c.centre.y+1);
    for(const auto& f:m.faces)bc.push_back(2*f.centre.x-3*f.centre.y+1);
    for(const auto& g:fv::reconstructGradient(m,field,bc))
        require(std::abs(g.x-2)<1e-12&&std::abs(g.y+3)<1e-12,"nonuniform mesh lost affine gradient");
    for(const auto& bad:std::vector<std::vector<double>>{{},{0},{0,0},{1,0},{0,std::numeric_limits<double>::infinity()}})
        rejects([&]{(void)makeRectilinearMesh2D(bad,y);});
    rejects([&]{(void)makeRectilinearMesh2D({0,1},{0,1e-5,1});});
    rejects([&]{(void)makeRectilinearMesh2D(x,{0,std::numeric_limits<double>::quiet_NaN()});});
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
