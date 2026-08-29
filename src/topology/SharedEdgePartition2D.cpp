#include "cartmesh2d/topology/SharedEdgePartition2D.hpp"
#include <algorithm>
#include <cmath>
#include <tuple>

namespace cartmesh2d {
SharedEdgePartition2D::SharedEdgePartition2D(const IntersectionRegistry2D& registry,
                                           std::vector<std::size_t> active)
    : handles_(std::move(active)) {
    const auto& vertices=registry.vertices();
    std::sort(handles_.begin(),handles_.end(),[&](auto a,auto b) {
        const auto& p=vertices.at(a).point;const auto& q=vertices.at(b).point;
        return std::tie(p.x,p.y,a)<std::tie(q.x,q.y,b);
    });
    handles_.erase(std::unique(handles_.begin(),handles_.end()),handles_.end());
    for (std::size_t i=0;i<handles_.size();++i) {
        const auto p=vertices.at(handles_[i]).point;
        points_.push_back(p);denseIds_.emplace(handles_[i],i);
        columns_[p.x].emplace(p.y,i);
    }
}

const std::vector<std::pair<double,std::size_t>>& SharedEdgePartition2D::partition(
    std::size_t a,std::size_t b,double eps,const TolerancePolicy& tol) {
    // Canonical direction is dense-ID order, independent of polygon winding.
    if (a>b) std::swap(a,b);
    const auto key=std::pair(a,b);
    if (const auto it=partitions_.find(key);it!=partitions_.end()) {++cacheHits_;return it->second;}
    const auto p=points_.at(a),q=points_.at(b);const auto d=q-p;
    std::vector<std::pair<double,std::size_t>> split{{0,a},{1,b}};
    const double tEps=eps/std::max(std::sqrt(squaredNorm(d)),eps);
    for (auto x=columns_.lower_bound(std::min(p.x,q.x)-eps);
         x!=columns_.end() && x->first<=std::max(p.x,q.x)+eps;++x) {
        for (auto y=x->second.lower_bound(std::min(p.y,q.y)-eps);
             y!=x->second.end() && y->first<=std::max(p.y,q.y)+eps;++y) {
            const auto id=y->second;
            if (id==a || id==b || !pointOnSegment(points_[id],{p,q},tol)) continue;
            const double t=dot(points_[id]-p,d)/squaredNorm(d);
            if (t>tEps && t<1-tEps) split.push_back({t,id});
        }
    }
    std::sort(split.begin(),split.end());
    return partitions_.emplace(key,std::move(split)).first->second;
}
}
