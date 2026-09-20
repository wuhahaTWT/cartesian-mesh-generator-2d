#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d {
namespace {
void check(bool ok,const char* message) {
    if(!ok)throw std::runtime_error(message);
}
std::size_t product(std::size_t a,std::size_t b) {
    check(!b || a<=std::numeric_limits<std::size_t>::max()/b,"Rectilinear mesh size overflow");
    return a*b;
}
void coordinates(const std::vector<double>& values) {
    check(values.size()>=2,"Rectilinear mesh needs at least two coordinates per axis");
    for(std::size_t i=0;i<values.size();++i)
        check(std::isfinite(values[i]) && (!i || values[i]>values[i-1]),
              "Rectilinear coordinates must be finite and strictly increasing");
}
}
TopologyMesh2D makeRectilinearMesh2D(const std::vector<double>& x,const std::vector<double>& y) {
    coordinates(x);coordinates(y);
    const auto nx=x.size()-1,ny=y.size()-1,nc=product(nx,ny);
    const auto nh=product(nx,ny+1),nv=product(nx+1,ny);
    check(nh<=std::numeric_limits<std::size_t>::max()-nv,"Rectilinear face count overflow");
    TopologyMesh2D t;t.vertices.reserve(product(nx+1,ny+1));t.edges.reserve(nh+nv);t.cells.reserve(nc);
    const auto vertex=[&](std::size_t i,std::size_t j){return j*(nx+1)+i;};
    const auto horizontal=[&](std::size_t i,std::size_t j){return j*nx+i;};
    const auto vertical=[&](std::size_t i,std::size_t j){return nh+j*(nx+1)+i;};
    for(std::size_t j=0;j<=ny;++j)for(std::size_t i=0;i<=nx;++i)
        t.vertices.push_back({vertex(i,j),{x[i],y[j]}});
    // Shared edges are stored once, oriented counterclockwise for their owner.
    for(std::size_t j=0;j<=ny;++j)for(std::size_t i=0;i<nx;++i) {
        const auto a=vertex(i,j),b=vertex(i+1,j);
        Edge2D e{horizontal(i,j),j?b:a,j?a:b,j?(j-1)*nx+i:i,{},BoundaryPatch2D::DomainBoundary};
        if(j&&j<ny){e.neighbour=j*nx+i;e.patch=BoundaryPatch2D::None;}
        t.edges.push_back(e);
    }
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<=nx;++i) {
        const auto a=vertex(i,j),b=vertex(i,j+1);
        Edge2D e{vertical(i,j),i?a:b,i?b:a,j*nx+(i?i-1:0),{},BoundaryPatch2D::DomainBoundary};
        if(i&&i<nx){e.neighbour=j*nx+i;e.patch=BoundaryPatch2D::None;}
        t.edges.push_back(e);
    }
    for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i) {
        TopologyCell2D c;c.id=j*nx+i;c.sourceId=c.id;
        c.geometryArea=(x[i+1]-x[i])*(y[j+1]-y[j]);
        check(std::isfinite(c.geometryArea)&&c.geometryArea>0,"Rectilinear cell area out of range");
        c.vertices={vertex(i,j),vertex(i+1,j),vertex(i+1,j+1),vertex(i,j+1)};
        c.edges={horizontal(i,j),vertical(i+1,j),horizontal(i,j+1),vertical(i,j)};
        t.cells.push_back(std::move(c));
    }
    const auto quality=evaluateSolverQuality2D(t);
    check(quality.valid(),"Rectilinear mesh fails the default Solver quality gate");
    return t;
}
}
