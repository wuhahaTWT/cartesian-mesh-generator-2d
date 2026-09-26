#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <map>
namespace fv_test {
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
inline FvMesh2D rectangle(int nx=24, int ny=8, double length=4, bool warped=false) {
    TopologyMesh2D t;
    for (int j=0;j<=ny;++j) for (int i=0;i<=nx;++i) {
        double x=length*i/nx,y=double(j)/ny;
        if(warped && i>0 && i<nx && j>0 && j<ny) {
            const double pi=std::acos(-1.),amplitude=.2*std::min(length/nx,1./ny);
            x+=amplitude*std::sin(2*pi*length*i/nx/length)*std::sin(pi*double(j)/ny);
            y+=amplitude*std::sin(pi*length*i/nx/length)*std::sin(2*pi*double(j)/ny);
        }
        t.vertices.push_back({t.vertices.size(),{x,y}});
    }
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for (int j=0;j<ny;++j) for (int i=0;i<nx;++i) {
        TopologyCell2D cell;
        cell.id=t.cells.size();cell.geometryArea=length/nx/ny;
        const auto a=static_cast<std::size_t>(j*(nx+1)+i);
        cell.vertices={a,a+1,a+static_cast<std::size_t>(nx)+2,a+static_cast<std::size_t>(nx)+1};
        if(warped) {
            cell.geometryArea=0;
            for(std::size_t k=0;k<4;++k) {
                const auto p=t.vertices[cell.vertices[k]].point,q=t.vertices[cell.vertices[(k+1)%4]].point;
                cell.geometryArea+=(p.x*q.y-p.y*q.x)/2;
            }
        }
        for (std::size_t k=0;k<4;++k) {
            const auto x=cell.vertices[k],y=cell.vertices[(k+1)%4];
            const auto [it,inserted]=edges.emplace(std::minmax(x,y),t.edges.size());
            if (inserted) t.edges.push_back({it->second,x,y,cell.id,{},BoundaryPatch2D::DomainBoundary});
            else {t.edges[it->second].neighbour=cell.id;t.edges[it->second].patch=BoundaryPatch2D::None;}
            cell.edges.push_back(it->second);
        }
        t.cells.push_back(cell);
    }
    return makeFvMesh2D(t);
}

}
