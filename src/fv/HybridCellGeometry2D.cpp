#include "cartmesh2d/fv/HybridCellGeometry2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace cartmesh2d::fv {
namespace {
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("hybrid geometry: ")+message);}
constexpr double epsBudget=4096*std::numeric_limits<double>::epsilon();
}
HybridCellGeometry2D hybridCellGeometry2D(const FvMesh2D& mesh,std::size_t cell) {
    check(cell<mesh.cells.size(),"invalid cell index");const auto& c=mesh.cells[cell];
    HybridCellGeometry2D out;out.faces=c.faces;out.area=c.area;check(std::isfinite(c.area)&&c.area>0,"nonpositive area");
    const auto n=out.faces.size();check(n>=3&&std::isfinite(c.centre.x)&&std::isfinite(c.centre.y),"invalid cell geometry");
    auto& s=out.areaVectors;auto& d=out.offsets;auto& g=out.gradientColumns;auto& w=out.stabilizationWeights;auto& r=out.remainder;
    s.resize(n);d.resize(n);g.resize(n);w.resize(n);r.resize(n*n);out.minimumNormalDistance=std::numeric_limits<double>::infinity();
        double sx=0,sy=0,perimeter=0;double moment[2][2]{},scale[2][2]{};
        for(std::size_t i=0;i<n;++i){
            check(out.faces[i]<mesh.faces.size(),"invalid face index");
            check(std::find(out.faces.begin(),out.faces.begin()+static_cast<std::ptrdiff_t>(i),out.faces[i])==out.faces.begin()+static_cast<std::ptrdiff_t>(i),"repeated cell face");
            const auto& f=mesh.faces[out.faces[i]];check(f.owner==cell||(f.neighbour&&*f.neighbour==cell),"unrelated face");const double sign=f.owner==cell?1:-1;
            s[i]={sign*f.areaVector.x,sign*f.areaVector.y};d[i]={f.centre.x-c.centre.x,f.centre.y-c.centre.y};
            const double length=std::hypot(s[i].x,s[i].y),normalDistance=dot(s[i],d[i])/length;
            check(std::isfinite(normalDistance)&&normalDistance>0,"cell centre outside strict face kernel");
            out.minimumNormalDistance=std::min(out.minimumNormalDistance,normalDistance);
            g[i]={s[i].x/c.area,s[i].y/c.area};w[i]=length/normalDistance;
            check(std::isfinite(g[i].x)&&std::isfinite(g[i].y)&&std::isfinite(w[i])&&w[i]>0,"nonfinite local geometry coefficient");
            sx+=s[i].x;sy+=s[i].y;perimeter+=length;
            const double ss[]={s[i].x,s[i].y},dd[]={d[i].x,d[i].y};
            for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b){moment[a][b]+=ss[a]*dd[b];scale[a][b]+=std::abs(ss[a]*dd[b]);}
        }
        const double closure=std::hypot(sx,sy)/(epsBudget*perimeter);out.closureRoundoffRatio=std::max(out.closureRoundoffRatio,closure);check(closure<=1,"face closure exceeds roundoff budget");
        for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b){double e=std::abs(moment[a][b]-(a==b?c.area:0))/(epsBudget*(c.area+scale[a][b]));out.momentRoundoffRatio=std::max(out.momentRoundoffRatio,e);check(e<=1,"first geometric moment exceeds roundoff budget");}
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){r[i*n+j]=(i==j?1.:0.)-dot(d[i],g[j]);check(std::isfinite(r[i*n+j]),"nonfinite geometric remainder");}
    return out;
}
}
