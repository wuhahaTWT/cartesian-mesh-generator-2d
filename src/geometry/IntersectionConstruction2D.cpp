#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace cartmesh2d {
namespace {
// Shared construction only absorbs arithmetic roundoff. The legacy sampling
// policy may be larger, but must not enlarge this construction budget.
constexpr double arithmeticFractionOfLocalH =
    IntersectionRegistryPolicy2D{}.snapFractionOfLocalH;
}

void IntersectionRegistry2D::configureGrid(const AABB2D& bounds, std::size_t level) {
    if (gridConfigured_ || level>30U || !std::isfinite(bounds.min.x) ||
        !std::isfinite(bounds.min.y) || !std::isfinite(bounds.max.x-bounds.min.x) ||
        !std::isfinite(bounds.max.y-bounds.min.y) || !(bounds.max.x>bounds.min.x) ||
        !(bounds.max.y>bounds.min.y)) throw std::invalid_argument("invalid or repeated construction grid");
    gridBounds_=bounds;
    gridLevel_=level;
    gridConfigured_=true;
}

GridLineIdentity2D IntersectionRegistry2D::gridLine(unsigned axis,double value) const {
    if (!gridConfigured_ || axis>1U || !std::isfinite(value))
        throw std::invalid_argument("invalid construction grid line");
    const double lo=axis==0U?gridBounds_.min.x:gridBounds_.min.y;
    const double hi=axis==0U?gridBounds_.max.x:gridBounds_.max.y;
    const auto count=std::uint64_t{1}<<gridLevel_;
    const double logical=(value-lo)/(hi-lo)*static_cast<double>(count);
    const double rounded=std::round(logical);
    // This validates an already known dyadic grid side, not a geometric snap.
    if (rounded<0 || rounded>static_cast<double>(count) ||
        std::abs(logical-rounded)>arithmeticFractionOfLocalH*static_cast<double>(count))
        throw std::invalid_argument("coordinate is not a construction grid side");
    return {axis,static_cast<std::uint64_t>(rounded)};
}

double IntersectionRegistry2D::gridCoordinate(GridLineIdentity2D line) const {
    if (!gridConfigured_ || line.axis>1U || line.coordinate>(std::uint64_t{1}<<gridLevel_))
        throw std::invalid_argument("invalid construction grid key");
    const double lo=line.axis==0U?gridBounds_.min.x:gridBounds_.min.y;
    const double hi=line.axis==0U?gridBounds_.max.x:gridBounds_.max.y;
    if (line.coordinate==(std::uint64_t{1}<<gridLevel_)) return hi;
    return lo+static_cast<double>(line.coordinate)*std::ldexp(hi-lo,-static_cast<int>(gridLevel_));
}

std::size_t IntersectionRegistry2D::internVertex(const Point2D& p,double h,IntersectionFeature2D feature) {
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !(h>0) || !std::isfinite(h))
        throw std::invalid_argument("invalid canonical vertex");
    const auto key=std::pair(p.x,p.y);
    if (const auto it=exactVertices_.find(key);it!=exactVertices_.end()) {
        auto& vertex=vertices_[it->second];
        const bool inputAnchor=feature==IntersectionFeature2D::Smooth ||
            feature==IntersectionFeature2D::TransitionVertex ||
            feature==IntersectionFeature2D::WallSharpCorner ||
            feature==IntersectionFeature2D::WallConcaveCorner;
        if (inputAnchor) {
            const auto found=vertexEvents_.find(vertex.id);
            if (found!=vertexEvents_.end()) for (const auto eventId:found->second) {
                const auto& event=events_[eventId];const auto& support=supports_[event.supportId];
                if (event.displacement>0 && vertex.id!=support.a && vertex.id!=support.b)
                    throw std::runtime_error("late nonincident feature conflicts with grid-corner snap; local support resampling required");
            }
        }
        vertex.localH=std::min(vertex.localH,h);
        // Metadata upgrades never move an existing vertex. Prefer the most
        // constrained classification, independent of registration order.
        const auto rank=[](IntersectionFeature2D f) {
            switch (f) {
            case IntersectionFeature2D::WallSharpCorner: return 0;
            case IntersectionFeature2D::WallConcaveCorner: return 1;
            case IntersectionFeature2D::TransitionVertex: return 2;
            case IntersectionFeature2D::Smooth: return 3;
            case IntersectionFeature2D::CartesianGridVertex: return 4;
            case IntersectionFeature2D::CartesianGridLine: return 5;
            case IntersectionFeature2D::None: return 6;
            }
            return 6;
        };
        if (rank(feature)<rank(vertex.feature)) vertex.feature=feature;
        return it->second;
    }
    const auto id=addCanonicalVertex(p,h,feature);
    exactVertices_.emplace(key,id);
    return id;
}

std::size_t IntersectionRegistry2D::registerSegment(const Segment2D& s,double h,IntersectionSource2D source) {
    const auto feature=source==IntersectionSource2D::TransitionEnvelopeCartesian
        ?IntersectionFeature2D::TransitionVertex:IntersectionFeature2D::Smooth;
    const auto a=internVertex(s.a,h,feature), b=internVertex(s.b,h,feature);
    if (a==b) throw std::invalid_argument("zero length construction support");
    const auto key=std::tuple(std::min(a,b),std::max(a,b),source);
    if (const auto it=supportKeys_.find(key);it!=supportKeys_.end()) return it->second;
    const auto id=supports_.size();
    // Store one direction for all callers, including a reversed segment.
    const bool forward=std::tie(s.a.x,s.a.y)<std::tie(s.b.x,s.b.y);
    supports_.push_back({forward?s:Segment2D{s.b,s.a},forward?a:b,forward?b:a,source});
    supportKeys_.emplace(key,id);
    return id;
}

std::string intersectionConstructionToJson(const IntersectionRegistry2D& registry,
    const std::vector<std::size_t>& handles,std::size_t partitions,std::size_t hits) {
    std::ostringstream out;out<<std::setprecision(17);
    out<<"{\"format\":\"cartmesh2d-shared-construction-v1\",\"intersection_evaluations\":"
       <<registry.events().size()<<",\"intersection_cache_hits\":"<<registry.intersectionCacheHits()
       <<",\"solver_partition_count\":"<<partitions<<",\"solver_partition_cache_hits\":"<<hits
       <<",\"solver_vertex_handles\":[";
    for (std::size_t i=0;i<handles.size();++i) out<<(i?",":"")<<handles[i];
    out<<"],\"events\":[";
    for (std::size_t i=0;i<registry.events().size();++i) {
        const auto& e=registry.events()[i];const auto& v=registry.vertices()[e.canonicalVertex];
        out<<(i?",":"")<<"{\"support_id\":"<<e.supportId<<",\"grid_axis\":"<<e.gridLine.axis
           <<",\"grid_coordinate\":"<<e.gridLine.coordinate<<",\"source\":\""<<intersectionSourceName(e.source)
           <<"\",\"original_position\":["<<e.originalPoint.x<<','<<e.originalPoint.y
           <<"],\"canonical_handle\":"<<v.id<<",\"canonical_position\":["<<v.point.x<<','<<v.point.y
           <<"],\"displacement\":"<<e.displacement<<",\"local_h\":"<<e.localH
           <<",\"displacement_over_local_h\":"<<e.displacement/e.localH
           <<",\"source_segment\":[["<<e.sourceSegment.a.x<<','<<e.sourceSegment.a.y
           <<"],["<<e.sourceSegment.b.x<<','<<e.sourceSegment.b.y<<"]]"
           <<",\"feature_classification\":\""<<intersectionFeatureName(v.feature)<<"\"}";
    }
    out<<"]}";return out.str();
}

std::size_t IntersectionRegistry2D::intersectGridLine(std::size_t support,GridLineIdentity2D line,double h) {
    if (!(h>0) || !std::isfinite(h)) throw std::invalid_argument("intersection requires positive local_h");
    const auto key=std::pair(support,line);
    if (const auto it=eventKeys_.find(key);it!=eventKeys_.end()) {
        ++cacheHits_;
        events_[it->second].localH=std::min(events_[it->second].localH,h);
        auto& vertex=vertices_[events_[it->second].canonicalVertex];
        vertex.localH=std::min(vertex.localH,h);
        return events_[it->second].canonicalVertex;
    }
    const auto& s=supports_.at(support);
    const auto d=s.segment.b-s.segment.a;
    const double delta=line.axis==0U?d.x:d.y;
    const double origin=line.axis==0U?s.segment.a.x:s.segment.a.y;
    const double target=gridCoordinate(line);
    if (delta==0) throw std::invalid_argument("parallel support has no isolated intersection");
    const double t=(target-origin)/delta;
    if (t<0 || t>1) throw std::invalid_argument("grid line lies outside construction support");
    const Point2D raw=s.segment.a+d*t;
    Point2D p=raw;
    if (line.axis==0U) p.x=target; else p.y=target;
    // Only incident endpoints or this support's arithmetic grid corner may
    // absorb roundoff. No nearest-feature search or Q1-sized movement occurs.
    const double gridH=std::ldexp(std::min(gridBounds_.max.x-gridBounds_.min.x,
                                         gridBounds_.max.y-gridBounds_.min.y),
                                  -static_cast<int>(gridLevel_));
    const double eps=arithmeticFractionOfLocalH*
        std::min({h,gridH,vertices_[s.a].localH,vertices_[s.b].localH,
                  std::sqrt(squaredNorm(d))});
    std::size_t id;
    if (std::sqrt(squaredNorm(p-s.segment.a))<=eps) id=s.a;
    else if (std::sqrt(squaredNorm(p-s.segment.b))<=eps) id=s.b;
    else {
        auto feature=IntersectionFeature2D::CartesianGridLine;
        const unsigned other=1U-line.axis;
        const double value=other==0U?p.x:p.y;
        const double lo=other==0U?gridBounds_.min.x:gridBounds_.min.y;
        const double hi=other==0U?gridBounds_.max.x:gridBounds_.max.y;
        const auto count=std::uint64_t{1}<<gridLevel_;
        const double logical=std::round((value-lo)/(hi-lo)*static_cast<double>(count));
        if (logical>=0 && logical<=static_cast<double>(count)) {
            const double targetOther=gridCoordinate({other,static_cast<std::uint64_t>(logical)});
            if (std::abs(targetOther-value)<=eps) {
                if (other==0U) p.x=targetOther;else p.y=targetOther;
                feature=IntersectionFeature2D::CartesianGridVertex;
            }
        }
        id=internVertex(p,h,feature);
        const auto anchorFeature=vertices_[id].feature;
        if (squaredNorm(p-raw)>0 &&
            (anchorFeature==IntersectionFeature2D::Smooth ||
             anchorFeature==IntersectionFeature2D::TransitionVertex ||
             anchorFeature==IntersectionFeature2D::WallSharpCorner ||
             anchorFeature==IntersectionFeature2D::WallConcaveCorner))
            throw std::runtime_error("unsafe snap onto nonincident feature vertex; local support resampling required");
    }
    const auto canonical=vertices_[id].point;
    // Check endpoints too: registration order cannot make an unrelated wall
    // incident to a previously moved intersection. Exact incidences are legal.
    for (const auto eventId:vertexEvents_[id]) {
        const auto& event=events_[eventId];
        if (event.supportId==support) continue;
        const auto& prior=supports_[event.supportId];
        if (s.a!=prior.a && s.a!=prior.b && s.b!=prior.a && s.b!=prior.b &&
            (event.displacement>0 || squaredNorm(canonical-raw)>0))
            throw std::runtime_error("unsafe nonincident grid-corner snap; local support resampling required");
    }
    eventKeys_.emplace(key,events_.size());
    vertexEvents_[id].push_back(events_.size());
    events_.push_back({support,line,raw,id,std::sqrt(squaredNorm(canonical-raw)),h,s.source,s.segment});
    return id;
}

} // namespace cartmesh2d
