#include "cartmesh2d/fv/WallGradient2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
constexpr std::size_t terms=6;
using Row=std::array<double,terms>;
struct PointSample {std::size_t index;bool boundary;Vector2D offset;};
struct ImageCell {std::size_t cell;Vector2D offset;};

// A single polynomial supplies the value and both derivatives at the face.
// Distance decay localizes the fit; column scaling and pivoted Householder QR
// handle dimensional units/aspect ratio without forming normal equations.
bool fit(const std::vector<PointSample>& points,Vector2D normal,Vector2D target,
         double h,WallGradientStencil2D& result) {
    const Vector2D tangent{-normal.y,normal.x};
    double hn=0,ht=0;
    for(const auto& p:points){hn=std::max(hn,std::abs(dot(p.offset,normal)));ht=std::max(ht,std::abs(dot(p.offset,tangent)));}
    if(points.size()<terms||!(hn>0&&ht>0&&h>0))return false;
    std::vector<Row> a(points.size()),q(points.size());
    std::vector<double> weight(points.size());
    Row columnNorm{};std::array<std::size_t,terms> permutation{0,1,2,3,4,5};
    for(std::size_t i=0;i<points.size();++i) {
        const double x=dot(points[i].offset,normal)/hn,y=dot(points[i].offset,tangent)/ht;
        const double r=std::hypot(points[i].offset.x,points[i].offset.y)/h;
        weight[i]=r<.5?1:std::pow(2*r,-5);
        a[i]={1,x,y,.5*x*x,x*y,.5*y*y};
        for(std::size_t j=0;j<terms;++j){a[i][j]*=weight[i];columnNorm[j]=std::hypot(columnNorm[j],a[i][j]);}
    }
    for(double n:columnNorm)if(!(n>0&&std::isfinite(n)))return false;
    for(auto& row:a)for(std::size_t j=0;j<terms;++j)row[j]/=columnNorm[j];
    std::array<std::vector<double>,terms> reflectors;
    double largest=0,smallest=1;
    for(std::size_t k=0;k<terms;++k) {
        std::size_t pivot=k;double best=-1;
        for(std::size_t j=k;j<terms;++j){double norm=0;for(std::size_t i=k;i<a.size();++i)norm=std::hypot(norm,a[i][j]);if(norm>best){best=norm;pivot=j;}}
        // Same dimensionless rank criterion as the existing wall-only QR.
        if(!(best>std::sqrt(std::numeric_limits<double>::epsilon())))return false;
        for(auto& row:a)std::swap(row[k],row[pivot]);std::swap(permutation[k],permutation[pivot]);
        auto& v=reflectors[k];v.resize(a.size()-k);
        for(std::size_t i=k;i<a.size();++i)v[i-k]=a[i][k];
        v[0]+=std::copysign(best,v[0]);double length=0;for(double x:v)length=std::hypot(length,x);for(auto& x:v)x/=length;
        for(std::size_t j=k;j<terms;++j){double projection=0;for(std::size_t i=k;i<a.size();++i)projection+=v[i-k]*a[i][j];for(std::size_t i=k;i<a.size();++i)a[i][j]-=2*v[i-k]*projection;}
        largest=std::max(largest,std::abs(a[k][k]));smallest=std::min(smallest,std::abs(a[k][k]));
    }
    for(std::size_t j=0;j<terms;++j)q[j][j]=1;
    for(std::size_t k=terms;k-->0;)for(std::size_t j=0;j<terms;++j){const auto& v=reflectors[k];double projection=0;for(std::size_t i=k;i<q.size();++i)projection+=v[i-k]*q[i][j];for(std::size_t i=k;i<q.size();++i)q[i][j]-=2*v[i-k]*projection;}
    const double tx=dot(target,normal)/hn,ty=dot(target,tangent)/ht;
    result.samples.clear();result.pivotRatio=smallest/largest;
    for(std::size_t i=0;i<points.size();++i) {
        Row x{},c{};
        for(std::size_t j=terms;j-->0;){double rhs=q[i][j];for(std::size_t k=j+1;k<terms;++k)rhs-=a[j][k]*x[k];x[j]=rhs/a[j][j];}
        for(std::size_t j=0;j<terms;++j)c[permutation[j]]=x[j]*weight[i]/columnNorm[permutation[j]];
        const double dn=(c[1]+c[3]*tx+c[4]*ty)/hn,dt=(c[2]+c[4]*tx+c[5]*ty)/ht;
        const Vector2D w{dn*normal.x+dt*tangent.x,dn*normal.y+dt*tangent.y};
        const double value=c[0]+c[1]*tx+c[2]*ty+.5*c[3]*tx*tx+c[4]*tx*ty+.5*c[5]*ty*ty;
        if(!std::isfinite(w.x)||!std::isfinite(w.y)||!std::isfinite(value))return false;
        result.samples.push_back({points[i].index,points[i].boundary,w,value});
    }
    return true;
}
}

WallGradientStencil2D quadraticFaceGradient2D(const FvMesh2D& mesh,std::size_t target,
    const std::vector<bool>& prescribed,const std::vector<std::optional<std::size_t>>& partners) {
    if(target>=mesh.faces.size()||prescribed.size()!=mesh.faces.size()||partners.size()!=mesh.faces.size())
        throw std::runtime_error("quadratic face gradient: invalid target or boundary map");
    const auto& face=mesh.faces[target];
    if(!face.neighbour&&!partners[target]&&!prescribed[target])throw std::runtime_error("quadratic face gradient: target has no prescribed value or neighbor");
    const double length=std::hypot(face.areaVector.x,face.areaVector.y);
    const Vector2D normal{face.areaVector.x/length,face.areaVector.y/length};
    const auto centre=mesh.cells[face.owner].centre;
    // A boundary fit is centered on its fluid owner, not anchored to a very
    // short wall segment. Internal fits are centered on the shared face.
    const auto origin=(!face.neighbour&&!partners[target])?centre:face.centre;
    const Vector2D evaluation{face.centre.x-origin.x,face.centre.y-origin.y};
    // Use geometric extent, not sqrt(cut area), which vanishes for slivers.
    double h=length;
    auto includeExtent=[&](std::size_t cell){const auto c=mesh.cells[cell].centre;for(auto id:mesh.cells[cell].faces){const auto& f=mesh.faces[id];h=std::max({h,std::hypot(f.areaVector.x,f.areaVector.y),2*std::hypot(f.centre.x-c.x,f.centre.y-c.y)});}};
    includeExtent(face.owner);
    if(face.neighbour)includeExtent(*face.neighbour);
    else if(partners[target]) {
        if(*partners[target]>=mesh.faces.size())throw std::runtime_error("quadratic face gradient: invalid periodic pair");
        includeExtent(mesh.faces[*partners[target]].owner);
    }
    const double epsilon=128*std::numeric_limits<double>::epsilon()*h;
    auto same=[&](Vector2D a,Vector2D b){return std::hypot(a.x-b.x,a.y-b.y)<=epsilon;};
    std::vector<ImageCell> cells{{face.owner,{centre.x-origin.x,centre.y-origin.y}}};
    auto addCell=[&](ImageCell candidate){for(const auto& old:cells)if(old.cell==candidate.cell&&same(old.offset,candidate.offset))return;cells.push_back(candidate);};
    auto otherCell=[&](const ImageCell& image,std::size_t id)->std::optional<ImageCell> {
        const auto& f=mesh.faces[id];const auto c=mesh.cells[image.cell].centre;
        std::optional<std::size_t> next;Vector2D delta{};
        if(f.neighbour){next=f.owner==image.cell?*f.neighbour:f.owner;const auto p=mesh.cells[*next].centre;delta={p.x-c.x,p.y-c.y};}
        else if(partners[id]) {
            const auto other=*partners[id];
            if(other>=mesh.faces.size()||mesh.faces[other].neighbour||partners[other]!=id)throw std::runtime_error("quadratic face gradient: invalid periodic pair");
            const auto& partner=mesh.faces[other];next=partner.owner;const auto p=mesh.cells[*next].centre;
            delta={(f.centre.x-c.x)+(p.x-partner.centre.x),(f.centre.y-c.y)+(p.y-partner.centre.y)};
        }
        if(!next)return {};
        return ImageCell{*next,{image.offset.x+delta.x,image.offset.y+delta.y}};
    };
    if(const auto other=otherCell(cells.front(),target))addCell(*other);
    std::vector<PointSample> points;
    auto addPoint=[&](PointSample p){for(const auto& old:points)if(old.index==p.index&&old.boundary==p.boundary&&same(old.offset,p.offset))return;points.push_back(p);};
    std::size_t begin=0,end=cells.size();WallGradientStencil2D result;
    for(unsigned ring=0;ring<=5;++ring) {
        for(std::size_t k=begin;k<end;++k) {
            const auto image=cells[k];const auto& cell=mesh.cells[image.cell];addPoint({image.cell,false,image.offset});
            for(auto id:cell.faces)if(prescribed[id]){const auto& f=mesh.faces[id];addPoint({id,true,{image.offset.x+f.centre.x-cell.centre.x,image.offset.y+f.centre.y-cell.centre.y}});}
        }
        if(ring>=2&&fit(points,normal,evaluation,h,result)){result.rings=ring;return result;}
        if(ring==5)break;
        for(std::size_t k=begin;k<end;++k)for(auto id:mesh.cells[cells[k].cell].faces) {
            if(const auto other=otherCell(cells[k],id))addCell(*other);
            if(cells.size()>512)throw std::runtime_error("quadratic face gradient: excessive stencil at face "+std::to_string(target));
        }
        begin=end;end=cells.size();
        if(begin==end){if(fit(points,normal,evaluation,h,result)){result.rings=ring;return result;}break;}
    }
    throw std::runtime_error("quadratic face gradient: rank-deficient stencil at face "+std::to_string(target)+"; refine the neighborhood or select linear gradients");
}
} // namespace cartmesh2d::fv
