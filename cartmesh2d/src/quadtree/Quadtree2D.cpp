#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
namespace cartmesh2d {
namespace {
struct LatticeFace { std::size_t leafIndex = 0; std::uint64_t begin = 0; std::uint64_t end = 0; };
using FaceMap = std::map<std::uint64_t, std::vector<LatticeFace>>;
std::uint64_t pow2(std::size_t exponent) { if (exponent >= 63) throw std::invalid_argument("quadtree level exceeds lattice capacity"); return std::uint64_t{1} << exponent; }
bool intervalOverlapPositive(std::uint64_t a0, std::uint64_t a1, std::uint64_t b0, std::uint64_t b1) noexcept { return std::max(a0,b0) < std::min(a1,b1); }
void appendMatches(const FaceMap& positiveFaces, const FaceMap& negativeFaces, std::vector<FaceNeighborPair2D>& result) {
    for (const auto& [coordinate, positives] : positiveFaces) {
        const auto it = negativeFaces.find(coordinate); if (it == negativeFaces.end()) continue;
        for (const auto& a : positives) for (const auto& b : it->second) {
            if (a.leafIndex == b.leafIndex || !intervalOverlapPositive(a.begin,a.end,b.begin,b.end)) continue;
            result.push_back({std::min(a.leafIndex,b.leafIndex), std::max(a.leafIndex,b.leafIndex)});
        }
    }
}
}
Point2D QuadtreeLeaf2D::center() const noexcept { return {(bounds.min.x+bounds.max.x)*0.5,(bounds.min.y+bounds.max.y)*0.5}; }
double QuadtreeLeaf2D::area() const noexcept { return (bounds.max.x-bounds.min.x)*(bounds.max.y-bounds.min.y); }
double distancePointToSegment(const Point2D& point, const Segment2D& segment) noexcept {
    const Vector2D ab = segment.b-segment.a; const double denom = squaredNorm(ab);
    if (denom <= 0.0) return std::sqrt(squaredNorm(point-segment.a));
    const double t = std::clamp(dot(point-segment.a,ab)/denom,0.0,1.0); const Point2D closest = segment.a + ab*t;
    return std::sqrt(squaredNorm(point-closest));
}
double distanceAABBToBoundary(const AABB2D& box, const BoundaryLoop& boundary, const TolerancePolicy& tol) {
    const auto& vertices = boundary.vertices(); if (vertices.size()<2) return std::numeric_limits<double>::infinity();
    const Point2D corners[4]={{box.min.x,box.min.y},{box.max.x,box.min.y},{box.max.x,box.max.y},{box.min.x,box.max.y}};
    double best=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<vertices.size();++i){ const Segment2D edge{vertices[i],vertices[(i+1)%vertices.size()]};
        if(segmentIntersectsClosedAABB(edge,box,tol)) return 0.0;
        for(const auto& corner:corners) best=std::min(best,distancePointToSegment(corner,edge));
        for(const auto& endpoint:{edge.a,edge.b}){ const double dx=endpoint.x<box.min.x?box.min.x-endpoint.x:endpoint.x>box.max.x?endpoint.x-box.max.x:0.0; const double dy=endpoint.y<box.min.y?box.min.y-endpoint.y:endpoint.y>box.max.y?endpoint.y-box.max.y:0.0; best=std::min(best,std::hypot(dx,dy)); }
    } return best;
}
Quadtree2D::Quadtree2D(Domain2D domain,std::size_t maxLevel,const BoundaryLoop& boundary,const TolerancePolicy& tol):domain_(domain),maxLevel_(maxLevel){
    if (!domain_.valid(tol)) throw std::invalid_argument("invalid quadtree domain");
    if (maxLevel_ > 28) throw std::invalid_argument("maxLevel > 28 unsupported");
    if (!boundary.diagnose(tol).valid()) throw std::invalid_argument("invalid boundary loop");
    leaves_.push_back(makeLeaf(0, 0, 0, boundary, tol));
    sortAndAssignIds();
}
std::uint64_t Quadtree2D::mortonPath(std::size_t level,std::uint64_t ix,std::uint64_t iy) noexcept { std::uint64_t code=0; for(std::size_t bit=0;bit<level;++bit){ code|=((ix>>bit)&1ULL)<<(2*bit); code|=((iy>>bit)&1ULL)<<(2*bit+1);} return code; }
std::uint64_t Quadtree2D::makeKey(std::size_t level,std::uint64_t ix,std::uint64_t iy) noexcept { return (mortonPath(level,ix,iy)<<6)|static_cast<std::uint64_t>(level); }
AABB2D Quadtree2D::boundsFor(std::size_t level,std::uint64_t ix,std::uint64_t iy) const noexcept { const double n=static_cast<double>(std::uint64_t{1}<<level); const double dx=domain_.width()/n,dy=domain_.height()/n; const double x0=domain_.bounds.min.x+static_cast<double>(ix)*dx,y0=domain_.bounds.min.y+static_cast<double>(iy)*dy; return {{x0,y0},{x0+dx,y0+dy}}; }
QuadtreeLeaf2D Quadtree2D::makeLeaf(std::size_t level,std::uint64_t ix,std::uint64_t iy,const BoundaryLoop& boundary,const TolerancePolicy& tol) const { QuadtreeLeaf2D leaf; leaf.key=makeKey(level,ix,iy); leaf.level=level; leaf.ix=ix; leaf.iy=iy; leaf.bounds=boundsFor(level,ix,iy); const CartesianCell2D cell{0,0,0,leaf.bounds,CellClass::Outside}; leaf.classification=classifyCartesianCell(cell,boundary,tol); return leaf; }
bool Quadtree2D::splitLeafAt(std::size_t index,const BoundaryLoop& boundary,const TolerancePolicy& tol){ if(index>=leaves_.size()) return false; const auto parent=leaves_[index]; if(parent.level>=maxLevel_) return false; const std::size_t l=parent.level+1; const std::uint64_t x=parent.ix*2,y=parent.iy*2; QuadtreeLeaf2D children[4]={makeLeaf(l,x,y,boundary,tol),makeLeaf(l,x+1,y,boundary,tol),makeLeaf(l,x,y+1,boundary,tol),makeLeaf(l,x+1,y+1,boundary,tol)}; leaves_.erase(leaves_.begin()+static_cast<std::ptrdiff_t>(index)); leaves_.insert(leaves_.end(),std::begin(children),std::end(children)); sortAndAssignIds(); return true; }
bool Quadtree2D::refineLeafByKey(std::uint64_t key,const BoundaryLoop& boundary,const TolerancePolicy& tol){ const auto it=std::find_if(leaves_.begin(),leaves_.end(),[key](const auto& leaf){return leaf.key==key;}); if(it==leaves_.end()) return false; return splitLeafAt(static_cast<std::size_t>(std::distance(leaves_.begin(),it)),boundary,tol); }
void Quadtree2D::refine(const BoundaryLoop& boundary,const QuadtreeRefinementPolicy2D& policy,const TolerancePolicy& tol){ if(!boundary.diagnose(tol).valid()) throw std::invalid_argument("invalid boundary loop"); if(policy.boundaryLevel>maxLevel_) throw std::invalid_argument("boundaryLevel exceeds maxLevel"); for(const auto& band:policy.distanceBands){ if(band.distance<0.0) throw std::invalid_argument("negative distance"); if(band.targetLevel>maxLevel_) throw std::invalid_argument("distance level exceeds maxLevel"); }
    bool changed=true; while(changed){ changed=false; for(std::size_t i=0;i<leaves_.size();++i){ const auto leaf=leaves_[i]; std::size_t requested=leaf.classification==CellClass::Intersected?policy.boundaryLevel:0; const double distance=distanceAABBToBoundary(leaf.bounds,boundary,tol); for(const auto& band:policy.distanceBands) if(distance<=band.distance+tol.scale(std::max(1.0,band.distance))) requested=std::max(requested,band.targetLevel); if (leaf.level < requested) { changed = splitLeafAt(i, boundary, tol); if (changed) break; } } }
}
std::vector<FaceNeighborPair2D> Quadtree2D::faceNeighbors() const { FaceMap left,right,bottom,top; for(std::size_t index=0;index<leaves_.size();++index){ const auto& leaf=leaves_[index]; const std::uint64_t scale=pow2(maxLevel_-leaf.level),x0=leaf.ix*scale,x1=(leaf.ix+1)*scale,y0=leaf.iy*scale,y1=(leaf.iy+1)*scale; left[x0].push_back({index,y0,y1}); right[x1].push_back({index,y0,y1}); bottom[y0].push_back({index,x0,x1}); top[y1].push_back({index,x0,x1}); }
    std::vector<FaceNeighborPair2D> result; appendMatches(right,left,result); appendMatches(top,bottom,result); std::sort(result.begin(),result.end(),[](const auto&a,const auto&b){return std::tie(a.first,a.second)<std::tie(b.first,b.second);}); result.erase(std::unique(result.begin(),result.end(),[](const auto&a,const auto&b){return a.first==b.first&&a.second==b.second;}),result.end()); return result; }
std::size_t Quadtree2D::countBalanceViolations() const { std::size_t count=0; for(const auto& p:faceNeighbors()){ const auto a=leaves_[p.first].level,b=leaves_[p.second].level; if((a>b?a-b:b-a)>1) ++count; } return count; }
QuadtreeBalanceReport2D Quadtree2D::enforceTwoToOneBalance(const BoundaryLoop& boundary,const TolerancePolicy& tol){ QuadtreeBalanceReport2D report; report.violationsBefore=countBalanceViolations(); while(true){ std::vector<std::uint64_t> keys; for(const auto& p:faceNeighbors()){ const auto&a=leaves_[p.first]; const auto&b=leaves_[p.second]; const auto diff=a.level>b.level?a.level-b.level:b.level-a.level; if(diff<=1) continue; keys.push_back((a.level<b.level?a:b).key); } if(keys.empty()) break; std::sort(keys.begin(),keys.end()); keys.erase(std::unique(keys.begin(),keys.end()),keys.end()); for(auto key:keys) if(refineLeafByKey(key,boundary,tol)) ++report.refinedLeaves; ++report.iterations; if(report.iterations>maxLevel_+1) throw std::runtime_error("2:1 balance failed to converge"); } report.violationsAfter=countBalanceViolations(); return report; }
double Quadtree2D::totalLeafArea() const noexcept { double total=0.0; for(const auto& leaf:leaves_) total+=leaf.area(); return total; }
bool Quadtree2D::deterministicOrderingValid() const noexcept { for(std::size_t i=0;i<leaves_.size();++i) if(leaves_[i].id!=i) return false; return true; }
void Quadtree2D::sortAndAssignIds(){ std::sort(leaves_.begin(),leaves_.end(),[this](const auto&a,const auto&b){ const std::uint64_t as=std::uint64_t{1}<<(2*(maxLevel_-a.level)),bs=std::uint64_t{1}<<(2*(maxLevel_-b.level)); const auto aa=mortonPath(a.level,a.ix,a.iy)*as,ba=mortonPath(b.level,b.ix,b.iy)*bs; if(aa!=ba) return aa<ba; if(a.level!=b.level) return a.level<b.level; return a.key<b.key; }); for(std::size_t i=0;i<leaves_.size();++i) leaves_[i].id=i; }
} // namespace cartmesh2d
