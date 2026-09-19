#include "cartmesh2d/fv/WallDistance2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* error) {if(!ok)throw std::runtime_error(error);}
struct Segment { Point2D a,b; std::size_t face; };
struct Box { double x0,y0,x1,y1; };
struct Node { Box box{}; std::size_t begin=0,end=0,left=0,right=0; bool leaf=false; };
Box bounds(const Segment& s) {
    return {std::min(s.a.x,s.b.x),std::min(s.a.y,s.b.y),std::max(s.a.x,s.b.x),std::max(s.a.y,s.b.y)};
}
Box merged(Box a,Box b) {
    return {std::min(a.x0,b.x0),std::min(a.y0,b.y0),std::max(a.x1,b.x1),std::max(a.y1,b.y1)};
}
double boxDistance(Point2D p,Box b) {
    return std::hypot(std::max({b.x0-p.x,0.,p.x-b.x1}),std::max({b.y0-p.y,0.,p.y-b.y1}));
}
double segmentDistance(Point2D p,const Segment& s) {
    const auto d=s.b-s.a;const double length=std::hypot(d.x,d.y);
    require(std::isfinite(length)&&length>0,"Wall distance degenerate segment");
    const Vector2D tangent{d.x/length,d.y/length};
    const double projection=dot(p-s.a,tangent);
    require(std::isfinite(projection),"Wall distance numerical range exceeded");
    const auto closest=s.a+tangent*std::clamp(projection,0.,length);
    const double result=std::hypot(p.x-closest.x,p.y-closest.y);
    require(std::isfinite(result),"Wall distance numerical range exceeded");return result;
}
class Index {
    std::vector<Segment> segments;
    std::vector<std::size_t> order;
    std::vector<Node> nodes;
    std::size_t build(std::size_t begin,std::size_t end) {
        const auto id=nodes.size();nodes.emplace_back();
        auto box=bounds(segments[order[begin]]);
        for(auto i=begin+1;i<end;++i)box=merged(box,bounds(segments[order[i]]));
        nodes[id].box=box;nodes[id].begin=begin;nodes[id].end=end;
        if(end-begin<=8) {nodes[id].leaf=true;return id;}
        const bool x=(box.x1-box.x0)>=(box.y1-box.y0);
        std::stable_sort(order.begin()+static_cast<std::ptrdiff_t>(begin),
                         order.begin()+static_cast<std::ptrdiff_t>(end),[&](auto a,auto b){
            const auto aa=bounds(segments[a]),bb=bounds(segments[b]);
            const double ca=x?std::midpoint(aa.x0,aa.x1):std::midpoint(aa.y0,aa.y1);
            const double cb=x?std::midpoint(bb.x0,bb.x1):std::midpoint(bb.y0,bb.y1);
            return ca!=cb?ca<cb:segments[a].face<segments[b].face;
        });
        const auto mid=begin+(end-begin)/2;
        nodes[id].left=build(begin,mid);nodes[id].right=build(mid,end);return id;
    }
    void search(std::size_t id,Point2D p,double& best,std::size_t& face,std::size_t& tests) const {
        const auto& node=nodes[id];
        if(boxDistance(p,node.box)>best)return;
        if(node.leaf) {
            for(auto i=node.begin;i<node.end;++i) {
                const auto& s=segments[order[i]];++tests;const double d=segmentDistance(p,s);
                if(d<best||(d==best&&s.face<face)){best=d;face=s.face;}
            }
            return;
        }
        const bool leftFirst=boxDistance(p,nodes[node.left].box)<=boxDistance(p,nodes[node.right].box);
        search(leftFirst?node.left:node.right,p,best,face,tests);
        search(leftFirst?node.right:node.left,p,best,face,tests);
    }
public:
    explicit Index(std::vector<Segment> input):segments(std::move(input)) {
        order.resize(segments.size());std::iota(order.begin(),order.end(),0);
        nodes.reserve(2*segments.size());build(0,segments.size());
    }
    void nearest(Point2D p,double& distance,std::size_t& face,std::size_t& tests) const {
        distance=std::numeric_limits<double>::infinity();face=std::numeric_limits<std::size_t>::max();
        search(0,p,distance,face,tests);
    }
};
}
WallDistanceResult2D computeWallDistance2D(const FvMesh2D& mesh,const std::vector<bool>& walls) {
    validateFvMesh2D(mesh);
    require(walls.size()==mesh.faces.size(),"Wall distance mask dimensions mismatch");
    std::vector<Segment> segments;
    for(std::size_t id=0;id<walls.size();++id)if(walls[id]) {
        const auto& f=mesh.faces[id];
        require(!f.neighbour,"Wall distance wall mask includes internal face");
        const Vector2D halfTangent{-.5*f.areaVector.y,.5*f.areaVector.x};
        const auto a=f.centre+halfTangent*(-1),b=f.centre+halfTangent;
        require(std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(b.x)&&std::isfinite(b.y),
            "Wall distance invalid segment endpoints");
        segments.push_back({a,b,id});
    }
    require(!segments.empty(),"Wall distance requires at least one wall face");
    Index index(std::move(segments));WallDistanceResult2D result;
    result.distance.resize(mesh.cells.size());result.nearestFace.resize(mesh.cells.size());
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        index.nearest(mesh.cells[i].centre,result.distance[i],result.nearestFace[i],result.segmentTests);
        require(std::isfinite(result.distance[i])&&result.distance[i]>0,
            "Wall distance nonpositive cell-to-wall distance");
    }
    return result;
}
}
