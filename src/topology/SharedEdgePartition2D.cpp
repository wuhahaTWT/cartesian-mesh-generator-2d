#include "cartmesh2d/topology/SharedEdgePartition2D.hpp"
#include <algorithm>
#include <cmath>
#include <tuple>

namespace cartmesh2d {
SharedEdgePartition2D::SharedEdgePartition2D(const IntersectionRegistry2D& registry,
                                           std::vector<std::size_t> active)
    : handles_(std::move(active)) {
    const auto& vertices=registry.vertices();
    // Each polygon repeats its incident handles. Deduplicate cheap integer
    // keys before the coordinate sort, retaining the same final (x,y,id) order.
    std::sort(handles_.begin(),handles_.end());
    handles_.erase(std::unique(handles_.begin(),handles_.end()),handles_.end());
    std::sort(handles_.begin(),handles_.end(),[&](auto a,auto b) {
        const auto& p=vertices.at(a).point;const auto& q=vertices.at(b).point;
        return std::tie(p.x,p.y,a)<std::tie(q.x,q.y,b);
    });
    denseIds_.assign(vertices.size(),std::numeric_limits<std::size_t>::max());
    points_.reserve(handles_.size());
    for (std::size_t i=0;i<handles_.size();++i) {
        const auto p=vertices.at(handles_[i]).point;
        points_.push_back(p);denseIds_[handles_[i]]=i;
        // Coordinate order already supplies sorted columns and sorted rows.
        // Match map::emplace's first-id rule for coincident registry entries.
        if (columns_.empty() || columns_.back().x!=p.x)
            columns_.push_back({p.x,{}});
        auto& rows=columns_.back().byY;
        if (rows.empty() || rows.back().first!=p.y) rows.emplace_back(p.y,i);
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
    for (auto x=std::lower_bound(columns_.begin(),columns_.end(),std::min(p.x,q.x)-eps,
                                [](const Column& column,double value) { return column.x<value; });
         x!=columns_.end() && x->x<=std::max(p.x,q.x)+eps;++x) {
        for (auto y=std::lower_bound(x->byY.begin(),x->byY.end(),std::min(p.y,q.y)-eps,
                                    [](const auto& row,double value) { return row.first<value; });
             y!=x->byY.end() && y->first<=std::max(p.y,q.y)+eps;++y) {
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
