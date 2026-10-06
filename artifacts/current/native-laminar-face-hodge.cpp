// Research-only face-flux Hodge prototype. Not a Navier--Stokes replacement.
// q_f = integral(U.n ds), globally oriented by the product owner normal.
// For each untouched polygon let N_f=S_f, R_f=C_f-C_i. N^T R=V I.
// H_i=R R^T/V + alpha*(I-N*(N^T N)^(-1)*N^T), alpha=tr(R R^T/V)/2.
// Then H_i N=R, H_i is SPD, and q^T H_i q=V |U|^2 for constant U.
// B is exact signed incidence: pressure force is -B^T p, divergence B q.
// This pairs pressure work without discarding unequal distance or skewness.
// Native assembly/geometry/manufactured data; Python only sparse linear algebra.
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <fstream>
#include <iomanip>
#include <limits>
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cartmesh2d;
using cartmesh2d::fv::FvMesh2D;
TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D topology;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = topology.cells.size();
        cell.geometryArea = polygon.area();
        for (const auto point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            const auto [it, inserted] = vertices.emplace(key, topology.vertices.size());
            if (inserted) topology.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            const auto key = std::minmax(a, b);
            const auto [it, inserted] = edges.emplace(key, topology.edges.size());
            if (inserted) {
                topology.edges.push_back({it->second, a, b, cell.id, {},
                                          BoundaryPatch2D::DomainBoundary});
            } else {
                auto& edge = topology.edges[it->second];
                edge.neighbour = cell.id;
                edge.patch = BoundaryPatch2D::None;
            }
            cell.edges.push_back(it->second);
        }
        topology.cells.push_back(std::move(cell));
    }
    return topology;
}

FvMesh2D rectangularMesh(std::size_t nx, std::size_t ny, double shear = 0.0) {
    std::vector<Polygon2D> polygons;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const auto point = [shear](std::size_t x, std::size_t y) {
                const double yy = static_cast<double>(y);
                return Point2D{static_cast<double>(x) + shear * yy, yy};
            };
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1),
                                 point(i, j + 1)}});
        }
    }
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons(polygons));
}

FvMesh2D normalizedSkewMesh(std::size_t n) {
    std::vector<Polygon2D> polygons;
    const double inverse = 1.0 / static_cast<double>(n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const auto point = [inverse](std::size_t x, std::size_t y) {
                return Point2D{inverse * (static_cast<double>(x) +
                                          0.3 * static_cast<double>(y)),
                               inverse * static_cast<double>(y)};
            };
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1),
                                 point(i, j + 1)}});
        }
    }
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons(polygons));
}



int main(int argc,char**argv) {
    if(argc!=3)throw std::runtime_error("usage: face-hodge mesh|skew:N prefix");
    const std::string input=argv[1],prefix=argv[2];
    FvMesh2D m;
    if(input.starts_with("skew:"))m=normalizedSkewMesh(std::stoul(input.substr(5)));
    else {const auto t=readCm2dTopology(input);if(!t.valid())throw std::runtime_error(t.error);m=cartmesh2d::fv::makeFvMesh2D(t.topology);}
    const auto n=m.cells.size(),nf=m.faces.size();
    std::ofstream matrix(prefix+".hodge.csv"),cells(prefix+".cells.csv"),faces(prefix+".faces.csv");
    for(auto*out:{&matrix,&cells,&faces}){if(!*out)throw std::runtime_error("output open");*out<<std::setprecision(17);}
    matrix<<"row,col,value\n";
    double momentMax=0,constantMax=0,minPivot=std::numeric_limits<double>::max(),area=0;
    std::size_t nonzeros=0;
    for(const auto&cell:m.cells) {
        const auto i=static_cast<std::size_t>(&cell-m.cells.data()),k=cell.faces.size();
        std::vector<Vector2D> normal(k),r(k);std::vector<double>sign(k),h(k*k);
        double xx=0,xy=0,yy=0,rr=0,identity[4]{};
        for(std::size_t a=0;a<k;++a){const auto&f=m.faces[cell.faces[a]];sign[a]=f.owner==i?1:-1;
            normal[a]={sign[a]*f.areaVector.x,sign[a]*f.areaVector.y};r[a]=f.centre-cell.centre;
            xx+=normal[a].x*normal[a].x;xy+=normal[a].x*normal[a].y;yy+=normal[a].y*normal[a].y;
            rr+=dot(r[a],r[a]);identity[0]+=normal[a].x*r[a].x;identity[1]+=normal[a].x*r[a].y;
            identity[2]+=normal[a].y*r[a].x;identity[3]+=normal[a].y*r[a].y;}
        for(int a=0;a<4;++a)momentMax=std::max(momentMax,std::abs(identity[a]/cell.area-((a==0||a==3)?1:0)));
        const double det=xx*yy-xy*xy,alpha=rr/(2*cell.area);
        if(!(det>0 && alpha>0))throw std::runtime_error("invalid geometric rank");
        for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<k;++b){
            const double projector=(normal[a].x*(yy*normal[b].x-xy*normal[b].y)+normal[a].y*(xx*normal[b].y-xy*normal[b].x))/det;
            h[a*k+b]=dot(r[a],r[b])/cell.area+alpha*((a==b?1:0)-projector);
            matrix<<cell.faces[a]<<','<<cell.faces[b]<<','<<sign[a]*sign[b]*h[a*k+b]<<'\n';++nonzeros;}
        for(std::size_t a=0;a<k;++a)for(int axis=0;axis<2;++axis){double hn=0;
            for(std::size_t b=0;b<k;++b)hn+=h[a*k+b]*(axis?normal[b].y:normal[b].x);
            constantMax=std::max(constantMax,std::abs(hn-(axis?r[a].y:r[a].x)));}
        // Verify local SPD numerically, with no tuned spectral floor.
        auto chol=h;
        for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<=a;++b){double v=chol[a*k+b];
            for(std::size_t c=0;c<b;++c)v-=chol[a*k+c]*chol[b*k+c];
            if(a==b){minPivot=std::min(minPivot,v);if(!(v>0))throw std::runtime_error("Hodge Cholesky failure");chol[a*k+b]=std::sqrt(v);}
            else chol[a*k+b]=v/chol[b*k+b];}
        area+=cell.area;
    }
    double xmin=1e300,xmax=-1e300,ymin=1e300,ymax=-1e300;
    for(const auto&f:m.faces){xmin=std::min(xmin,f.centre.x);xmax=std::max(xmax,f.centre.x);ymin=std::min(ymin,f.centre.y);ymax=std::max(ymax,f.centre.y);}
    if(input.starts_with("skew:")){xmin=0;xmax=1.3;ymin=0;ymax=1;}
    const double pi=std::acos(-1.),kx=pi/(xmax-xmin),ky=pi/(ymax-ymin);
    const auto nonlinear=[&](Point2D p){return std::sin(kx*(p.x-xmin))*std::sin(ky*(p.y-ymin));};
    cells<<"cell,x,y,area,affinePressure,nonlinearPressure,sourceIntegral\n";
    for(std::size_t i=0;i<n;++i){const auto&c=m.cells[i];const double p=nonlinear(c.centre);
        cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<1+c.centre.x+.3*c.centre.y<<','<<p<<','<<(kx*kx+ky*ky)*p*c.area<<'\n';}
    faces<<"face,owner,neighbour,x,y,sx,sy,affinePressure,nonlinearPressure,affineFlux,nonlinearFlux,trialFlux\n";
    for(std::size_t j=0;j<nf;++j){const auto&f=m.faces[j];const auto p=f.centre;
        const double gx=kx*std::cos(kx*(p.x-xmin))*std::sin(ky*(p.y-ymin));
        const double gy=ky*std::sin(kx*(p.x-xmin))*std::cos(ky*(p.y-ymin));
        const double trial=f.neighbour?(std::sin(3*p.x+2*p.y)*f.areaVector.x+std::cos(p.x-4*p.y)*f.areaVector.y):0;
        faces<<j<<','<<f.owner<<','<<(f.neighbour?static_cast<long>(*f.neighbour):-1)<<','<<p.x<<','<<p.y<<','<<f.areaVector.x<<','<<f.areaVector.y
            <<','<<1+p.x+.3*p.y<<','<<nonlinear(p)<<','<<-f.areaVector.x-.3*f.areaVector.y<<','<<-gx*f.areaVector.x-gy*f.areaVector.y<<','<<trial<<'\n';}
    for(auto*out:{&matrix,&cells,&faces}){out->close();if(!*out)throw std::runtime_error("output write");}
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<n<<",\"faces\":"<<nf<<",\"area\":"<<area<<",\"triplets\":"<<nonzeros
      <<",\"relativeMomentIdentityMax\":"<<momentMax<<",\"constantFluxConsistencyMax\":"<<constantMax
      <<",\"minimumLocalCholeskyPivot\":"<<minPivot<<"}\n";
}
