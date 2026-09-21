#include "cartmesh2d/fv/Diffusion2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
using namespace cartmesh2d;
namespace {
int failures=0;
void check(bool test,const std::string& message) { if(!test) {++failures;std::cerr<<"FAIL: "<<message<<'\n';} }
template<class F> void rejects(F run,const std::string& message) {
    try {run();check(false,message);} catch(const std::exception&) {}
}
TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D t; std::map<std::pair<double,double>,std::size_t> vertices;
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for(const auto& p:polygons) {
        TopologyCell2D c; c.id=t.cells.size();c.geometryArea=p.area();
        for(auto v:p.vertices) {
            auto [it,inserted]=vertices.emplace(std::make_pair(v.x,v.y),t.vertices.size());
            if(inserted) t.vertices.push_back({it->second,v});
            c.vertices.push_back(it->second);
        }
        for(std::size_t i=0;i<c.vertices.size();++i) {
            const auto a=c.vertices[i],b=c.vertices[(i+1)%c.vertices.size()];
            auto [it,inserted]=edges.emplace(std::minmax(a,b),t.edges.size());
            if(inserted) t.edges.push_back({it->second,a,b,c.id,{},BoundaryPatch2D::DomainBoundary});
            else {auto& e=t.edges[it->second];e.neighbour=c.id;e.patch=BoundaryPatch2D::None;}
            c.edges.push_back(it->second);
        }
        t.cells.push_back(c);
    }
    return t;
}
TopologyMesh2D skewGrid(int n) {
    const auto point=[&](int i,int j) {const double y=static_cast<double>(j)/n;return Point2D{static_cast<double>(i)/n+.25*y,y};};
    std::vector<Polygon2D> polygons;
    for(int j=0;j<n;++j) for(int i=0;i<n;++i) polygons.push_back({{point(i,j),point(i+1,j),point(i+1,j+1),point(i,j+1)}});
    return fromPolygons(polygons);
}
void linearTest(const TopologyMesh2D& t) {
    const auto mesh=fv::makeFvMesh2D(t);
    const auto exact=[](Point2D p) {return 1+p.x+2*p.y;};
    std::vector<double> u,bc(mesh.faces.size());
    for(const auto& c:mesh.cells) u.push_back(exact(c.centre));
    for(std::size_t i=0;i<bc.size();++i) bc[i]=exact(mesh.faces[i].centre);
    const auto g=fv::reconstructGradient(mesh,u,bc);
    for(auto v:g) check(std::abs(v.x-1)<1e-11 && std::abs(v.y-2)<1e-11,"linear gradient exact on final polygons");
    fv::DiffusionProblem2D problem{1,[](Point2D){return 0.;},[&](Point2D p,BoundaryPatch2D){return exact(p);}};
    const auto r=fv::solveDiffusion2D(mesh,problem);
    check(r.converged,"linear full corrected equation converges");
    for(std::size_t i=0;i<u.size();++i) check(std::abs(r.values[i]-u[i])<2e-8,"linear solution recovered");
    auto residual=r.sourceIntegrals; for(auto&v:residual)v=-v;
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id]; residual[f.owner]+=r.fluxes[id];
        if(f.neighbour) residual[*f.neighbour]-=r.fluxes[id];
    }
    for(double v:residual)check(std::abs(v)<1e-8,"independent signed cell flux balance");
    check(mesh.maxClosureError<1e-14,"normal closure");
    auto limited=fv::DiffusionControls2D{};limited.maxCorrections=1;
    check(!fv::solveDiffusion2D(mesh,problem,limited).converged,"outer limit is not false success");
    limited.maxCorrections=400;limited.maxLinearIterations=1;
    rejects([&]{(void)fv::solveDiffusion2D(mesh,problem,limited);},"linear solver limit must fail");
    problem.boundaryValue=[](Point2D,BoundaryPatch2D){return 0.;};
    const auto zero=fv::solveDiffusion2D(mesh,problem);
    check(zero.converged && zero.globalBalance==0,"zero source/zero boundary well defined");
}
double sineError(int n) {
    const auto mesh=fv::makeFvMesh2D(skewGrid(n));
    const auto exact=[](Point2D p){return std::sin(std::numbers::pi*p.x)*std::sin(std::numbers::pi*p.y);};
    fv::DiffusionProblem2D problem{1,[&](Point2D p){return 2*std::numbers::pi*std::numbers::pi*exact(p);},
        [&](Point2D p,BoundaryPatch2D){return exact(p);}};
    const auto r=fv::solveDiffusion2D(mesh,problem);check(r.converged,"sine convergence");
    double e=0,a=0;for(std::size_t i=0;i<mesh.cells.size();++i) {const auto& c=mesh.cells[i];e+=c.area*std::pow(r.values[i]-exact(c.centre),2);a+=c.area;}
    return std::sqrt(e/a);
}
}
int main() {
 try {
    const auto mesh=skewGrid(5);linearTest(mesh);
    linearTest(fromPolygons({{{{0,0},{1,0},{1,.5},{1,1},{0,1}}},{{{1,0},{2,0},{2,.5},{1,.5}}},{{{1,.5},{2,.5},{2,1},{1,1}}}}));
    const double coarse=sineError(8),fine=sineError(16);
    check(fine<coarse/2.5,"smooth manufactured solution improves on skew refinement");
    auto invalid=mesh;invalid.cells[0].geometryArea*=1.1;
    rejects([&]{(void)fv::makeFvMesh2D(invalid);},"forged clean audit does not hide wrong area");
    invalid=mesh;invalid.cells[0].edges[0]=invalid.cells[0].edges[1];
    rejects([&]{(void)fv::makeFvMesh2D(invalid);},"duplicate/mismatched edge rejected");
    invalid=mesh;invalid.edges[0].owner=mesh.cells.size()+1;
    rejects([&]{(void)fv::makeFvMesh2D(invalid);},"out-of-range owner rejected");
    invalid=mesh;invalid.edges[0].patch=BoundaryPatch2D::Unclassified;
    rejects([&]{(void)fv::makeFvMesh2D(invalid);},"unclassified boundary rejected");
    auto cached=fv::makeFvMesh2D(mesh);
    fv::DiffusionProblem2D extreme{1,[](Point2D){return 0.;},[](Point2D,BoundaryPatch2D){return 1e154;}};
    rejects([&]{(void)fv::solveDiffusion2D(cached,extreme);},"overflow cannot be reported as convergence");
    cached.faces[0].owner=cached.cells.size();
    rejects([&]{(void)fv::solveDiffusion2D(cached,extreme);},"public API rejects corrupt cached owner before indexing");
    std::cout<<"skew-grid sine L2: "<<coarse<<" -> "<<fine<<"; failures: "<<failures<<'\n';
 }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
 return failures?1:0;
}
