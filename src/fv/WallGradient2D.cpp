#include "cartmesh2d/fv/WallGradient2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
constexpr std::size_t terms=5;
using Row=std::array<double,terms>;
struct PointSample {std::size_t index;bool boundary;Vector2D offset;};
struct ImageCell {std::size_t cell;Vector2D offset;};
// Householder QR with column pivoting, in scaled wall-normal/tangent axes.
// No normal equations: squaring the condition number is particularly harmful
// in stretched layers and very small cut cells. Only two derivative rows of
// the left inverse are retained; all factorization work happens once.
bool fit(const std::vector<PointSample>& points,Vector2D normal,WallGradientStencil2D& result) {
    const Vector2D tangent{-normal.y,normal.x};
    double hn=0,ht=0;
    for(const auto& p:points){hn=std::max(hn,std::abs(dot(p.offset,normal)));ht=std::max(ht,std::abs(dot(p.offset,tangent)));}
    if(points.size()<terms||hn==0||ht==0)return false;
    std::vector<Row> a(points.size()),q(points.size());std::vector<double> weight(points.size());
    Row columnNorm{};std::array<std::size_t,terms> permutation{0,1,2,3,4};
    for(std::size_t i=0;i<points.size();++i) {
        const double x=dot(points[i].offset,normal)/hn,y=dot(points[i].offset,tangent)/ht;
        weight[i]=1/std::hypot(x,y);a[i]={x,y,.5*x*x,x*y,.5*y*y};
        for(std::size_t j=0;j<terms;++j){a[i][j]*=weight[i];columnNorm[j]=std::hypot(columnNorm[j],a[i][j]);}
    }
    for(double n:columnNorm)if(!(n>0&&std::isfinite(n)))return false;
    for(auto& row:a)for(std::size_t j=0;j<terms;++j)row[j]/=columnNorm[j];
    std::array<std::vector<double>,terms> reflectors;
    double largest=0,smallest=1;
    for(std::size_t k=0;k<terms;++k) {
        std::size_t pivot=k;double best=-1;
        for(std::size_t j=k;j<terms;++j){double norm=0;for(std::size_t i=k;i<a.size();++i)norm=std::hypot(norm,a[i][j]);if(norm>best){best=norm;pivot=j;}}
        // A dimensionless rank/roundoff criterion, not a physical-error gate.
        // sqrt(eps) leaves about half the double precision digits available.
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
    result.samples.clear();result.pivotRatio=smallest/largest;
    for(std::size_t i=0;i<points.size();++i) {
        Row x{},coefficient{};
        for(std::size_t j=terms;j-->0;){double rhs=q[i][j];for(std::size_t k=j+1;k<terms;++k)rhs-=a[j][k]*x[k];x[j]=rhs/a[j][j];}
        for(std::size_t j=0;j<terms;++j)coefficient[permutation[j]]=x[j]*weight[i]/columnNorm[permutation[j]];
        const Vector2D w{coefficient[0]*normal.x/hn+coefficient[1]*tangent.x/ht,coefficient[0]*normal.y/hn+coefficient[1]*tangent.y/ht};
        if(!std::isfinite(w.x)||!std::isfinite(w.y))return false;
        result.samples.push_back({points[i].index,points[i].boundary,w});
    }
    return true;
}
}
WallGradientStencil2D quadraticWallGradient2D(const FvMesh2D& mesh,std::size_t target,
    const std::vector<bool>& prescribed,const std::vector<std::optional<std::size_t>>& partners) {
    if(target>=mesh.faces.size()||prescribed.size()!=mesh.faces.size()||partners.size()!=mesh.faces.size()||
       mesh.faces[target].neighbour||!prescribed[target]||partners[target])throw std::runtime_error("quadratic wall gradient: invalid target or boundary map");
    const auto& wall=mesh.faces[target];const double length=std::hypot(wall.areaVector.x,wall.areaVector.y);
    const Vector2D normal{wall.areaVector.x/length,wall.areaVector.y/length};
    const auto centre=mesh.cells[wall.owner].centre;
    std::vector<ImageCell> cells{{wall.owner,{centre.x-wall.centre.x,centre.y-wall.centre.y}}};
    std::vector<PointSample> points;
    // Offset comparisons are local and dimensional; no global coordinate or
    // unit-sized tolerance causes images to collapse under SI rescaling.
    double scale=length;
    for(auto id:mesh.cells[wall.owner].faces){const auto& f=mesh.faces[id];scale=std::max(scale,std::hypot(f.centre.x-centre.x,f.centre.y-centre.y));}
    const double epsilon=128*std::numeric_limits<double>::epsilon()*scale;
    auto same=[&](Vector2D a,Vector2D b){return std::hypot(a.x-b.x,a.y-b.y)<=epsilon;};
    auto addPoint=[&](PointSample p){if(std::hypot(p.offset.x,p.offset.y)<=epsilon)return;for(const auto& old:points)if(old.index==p.index&&old.boundary==p.boundary&&same(old.offset,p.offset))return;points.push_back(p);};
    std::size_t begin=0,end=1;WallGradientStencil2D result;
    for(unsigned ring=0;ring<=5;++ring) {
        for(std::size_t k=begin;k<end;++k) {
            const auto image=cells[k];const auto& cell=mesh.cells[image.cell];addPoint({image.cell,false,image.offset});
            for(auto id:cell.faces) {
                const auto& face=mesh.faces[id];
                if(prescribed[id])addPoint({id,true,{image.offset.x+face.centre.x-cell.centre.x,image.offset.y+face.centre.y-cell.centre.y}});
            }
        }
        if(ring>=2&&fit(points,normal,result)){result.rings=ring;return result;}
        if(ring==5)break;
        for(std::size_t k=begin;k<end;++k) {
            const auto image=cells[k];const auto& cell=mesh.cells[image.cell];
            for(auto id:cell.faces) {
                const auto& face=mesh.faces[id];std::optional<std::size_t> next;Vector2D delta{};
                if(face.neighbour){next=face.owner==image.cell?*face.neighbour:face.owner;const auto c=mesh.cells[*next].centre;delta={c.x-cell.centre.x,c.y-cell.centre.y};}
                else if(partners[id]) {
                    const auto other=*partners[id];
                    if(other>=mesh.faces.size()||mesh.faces[other].neighbour||partners[other]!=id)throw std::runtime_error("quadratic wall gradient: invalid periodic pair");
                    const auto& partner=mesh.faces[other];next=partner.owner;const auto c=mesh.cells[*next].centre;
                    delta={(face.centre.x-cell.centre.x)+(c.x-partner.centre.x),(face.centre.y-cell.centre.y)+(c.y-partner.centre.y)};
                }
                if(!next)continue;
                const ImageCell candidate{*next,{image.offset.x+delta.x,image.offset.y+delta.y}};
                bool seen=false;for(const auto& old:cells)if(old.cell==candidate.cell&&same(old.offset,candidate.offset)){seen=true;break;}
                if(!seen)cells.push_back(candidate);
                if(cells.size()>512)throw std::runtime_error("quadratic wall gradient: excessive stencil at face "+std::to_string(target));
            }
        }
        begin=end;end=cells.size();if(begin==end){if(fit(points,normal,result)){result.rings=ring;return result;}break;}
    }
    throw std::runtime_error("quadratic wall gradient: rank-deficient stencil at face "+std::to_string(target)+"; refine the wall neighbourhood or select linear wall gradients");
}
} // namespace cartmesh2d::fv
