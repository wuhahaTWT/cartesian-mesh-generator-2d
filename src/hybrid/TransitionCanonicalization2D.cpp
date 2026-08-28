#include "cartmesh2d/hybrid/TransitionCanonicalization2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace cartmesh2d {
namespace {
double segmentLength(const Point2D& a,const Point2D& b) {
    return std::sqrt(squaredNorm(b-a));
}
bool pointOnRegionBoundary(const Point2D& p,const BoundaryRegion2D& region,
                           const TolerancePolicy& tol) {
    for (const auto& loop:region.loops()) {
        const auto& v=loop.vertices();
        for (std::size_t i=0;i<v.size();++i)
            if (pointOnSegment(p,{v[i],v[(i+1U)%v.size()]},tol)) return true;
    }
    return false;
}
bool edgeOnRegionBoundary(const Point2D& a,const Point2D& b,
                          const BoundaryRegion2D& region,const TolerancePolicy& tol) {
    return pointOnRegionBoundary(a,region,tol) && pointOnRegionBoundary(b,region,tol) &&
           pointOnRegionBoundary(a+(b-a)*0.5,region,tol);
}
double distanceToSegment(const Point2D& p,const Point2D& a,const Point2D& b) {
    const auto d=b-a;
    const double t=squaredNorm(d)>0.0?std::clamp(dot(p-a,d)/squaredNorm(d),0.0,1.0):0.0;
    return segmentLength(p,a+d*t);
}
// A displacement must stay inside the clearance to nonincident supports;
// proximity alone must never bridge a narrow gap or cross a fixed H4 layer.
double displacementBudget(const Point2D& p,
    const std::vector<BoundaryLoop>& loops,std::size_t loopId,std::size_t edgeId,
    bool atVertex,const BoundaryRegion2D& fixedInterface) {
    double clearance=std::numeric_limits<double>::infinity();
    for (std::size_t l=0;l<loops.size();++l) {
        const auto& v=loops[l].vertices();
        for (std::size_t e=0;e<v.size();++e) {
            if (l==loopId && (e==edgeId || (atVertex && (e+1U)%v.size()==edgeId))) continue;
            clearance=std::min(clearance,distanceToSegment(p,v[e],v[(e+1U)%v.size()]));
        }
    }
    for (const auto& loop:fixedInterface.loops()) {
        const auto& v=loop.vertices();
        for (std::size_t e=0;e<v.size();++e)
            clearance=std::min(clearance,distanceToSegment(p,v[e],v[(e+1U)%v.size()]));
    }
    return TransitionCanonicalizationPolicy2D{}.clearanceFraction*clearance;
}
[[nodiscard]] IntersectionFeature2D transitionVertexFeature(
    const std::vector<Point2D>& vertices, std::size_t index) noexcept {
    if (vertices.size()<3U) return IntersectionFeature2D::TransitionVertex;
    const auto& previous=vertices[(index+vertices.size()-1U)%vertices.size()];
    const auto& current=vertices[index];
    const auto& next=vertices[(index+1U)%vertices.size()];
    const Vector2D incoming=current-previous;
    const Vector2D outgoing=next-current;
    const double lengths=std::sqrt(squaredNorm(incoming)*squaredNorm(outgoing));
    if (!(lengths>0.0)) return IntersectionFeature2D::TransitionVertex;
    const double turn=std::atan2(std::abs(cross(incoming,outgoing)),
                                 dot(incoming,outgoing));
    // Polygonal sampling of a smooth curve has small turns. A 30 degree or
    // larger turn is retained as a real feature and is never moved to a grid
    // anchor by the registry.
    if (turn<TransitionCanonicalizationPolicy2D{}.featureTurnRadians)
        return IntersectionFeature2D::TransitionVertex;
    return cross(incoming,outgoing)>0.0
        ?IntersectionFeature2D::WallSharpCorner
        :IntersectionFeature2D::WallConcaveCorner;
}


} // namespace

bool canonicalizeTransitionEnvelope2D(
    std::vector<BoundaryLoop>& remainderBoundaryLoops,
    std::vector<Polygon2D>& transitionPolygons,
    const BoundaryRegion2D& outerRegion, const Domain2D& domain,
    std::size_t boundaryLevel, double transitionRingThickness,
    bool localTermination, IntersectionRegistry2D& envelopeRegistry,
    const TolerancePolicy& tolerance, std::string& error) {
    // Stepped fronts contain immutable H4 termination vertices. Local grid
    // refinement/front deformation experiments failed the unchanged solver
    // gate; retain those sources and their visible Q1 failures. Q2 is partial
    // until that route is safely resolved (see the Q2 validation document).
    if (localTermination) return true;
    const TransitionCanonicalizationPolicy2D canonicalPolicy;
    const double minimumFaceFraction=canonicalPolicy.minimumFaceFraction;
    // Align only smooth outer transition samples that would otherwise land
    // within one Q1 hard-face fraction of a Cartesian grid line. This changes
    // the transition sampling before cutting, so both sides construct the
    // same point; sharp/concave feature vertices are explicitly protected.
    const double gridHx=std::ldexp(domain.width(),
                                  -static_cast<int>(boundaryLevel));
    const double gridHy=std::ldexp(domain.height(),
                                  -static_cast<int>(boundaryLevel));
    std::size_t envelopeFeatureId=0U;
    std::size_t supportId=0U;
    for (auto& polygon:transitionPolygons) {
        std::vector<Point2D> partitioned;
        for (std::size_t e=0;e<polygon.vertices.size();++e) {
            const auto a=polygon.vertices[e];
            const auto b=polygon.vertices[(e+1U)%polygon.vertices.size()];
            partitioned.push_back(a);
            std::vector<std::pair<double,Point2D>> onEdge;
            const auto d=b-a;
            for (const auto& loop:remainderBoundaryLoops) for (const auto& p:loop.vertices()) {
                const double t=dot(p-a,d)/squaredNorm(d);
                if (t>0.0 && t<1.0 && !nearlyEqual(p,a,tolerance) &&
                    !nearlyEqual(p,b,tolerance) && pointOnSegment(p,{a,b},tolerance)) {
                    onEdge.push_back({t,p});
                }
            }
            std::sort(onEdge.begin(),onEdge.end(),[](const auto& lhs,const auto& rhs) {return lhs.first<rhs.first;});
            for (const auto& item:onEdge) partitioned.push_back(item.second);
        }
        polygon.vertices=std::move(partitioned);
    }
    for (std::size_t loopId=0U;loopId<remainderBoundaryLoops.size();++loopId) {
        const auto originalVertices=remainderBoundaryLoops[loopId].vertices();
        std::vector<Point2D> adjusted=originalVertices;
        for (std::size_t i=0U;i<adjusted.size();++i) {
            ++supportId;
            if (pointOnRegionBoundary(originalVertices[i],outerRegion,tolerance)) continue;
            const auto feature=transitionVertexFeature(originalVertices,i);
            const bool protectedVertex=
                feature==IntersectionFeature2D::WallSharpCorner ||
                feature==IntersectionFeature2D::WallConcaveCorner;
            const auto featureIdentity=protectedVertex
                ?std::optional<std::size_t>(envelopeFeatureId++)
                :std::nullopt;
            const auto& original=originalVertices[i];
            const double gx=domain.bounds.min.x+
                std::round((original.x-domain.bounds.min.x)/gridHx)*gridHx;
            const double gy=domain.bounds.min.y+
                std::round((original.y-domain.bounds.min.y)/gridHy)*gridHy;
            const Point2D projections[]{{gx,original.y},{original.x,gy},{gx,gy}};
            const double budget=displacementBudget(original,remainderBoundaryLoops,
                loopId,i,true,outerRegion);
            for (const auto& projection:projections) {
                if (segmentLength(projection,original)>budget) continue;
                // A small perpendicular distance to an almost parallel grid
                // line is not a short face. Require a genuinely nearby
                // intersection along an incident envelope segment as well.
                bool shortIncident=false;
                const auto previous=originalVertices[(i+originalVertices.size()-1U)%originalVertices.size()];
                const auto next=originalVertices[(i+1U)%originalVertices.size()];
                for (const auto& other:{previous,next}) {
                    const auto d=other-original;
                    for (int axis=0;axis<2;++axis) {
                        const double delta=axis==0?d.x:d.y;
                        const double target=axis==0?projection.x:projection.y;
                        const double origin=axis==0?original.x:original.y;
                        if (delta==0.0 || target==origin) continue;
                        const double t=(target-origin)/delta;
                        if (t>=0.0 && t<=1.0 && t*segmentLength(original,other)<
                            minimumFaceFraction*transitionRingThickness) shortIncident=true;
                    }
                }
                if (!shortIncident) continue;
                (void)envelopeRegistry.addCanonicalVertex(
                    projection,transitionRingThickness,
                    IntersectionFeature2D::CartesianGridLine,std::nullopt,supportId);
            }
            adjusted[i]=envelopeRegistry.canonicalize(
                original,transitionRingThickness,
                IntersectionSource2D::TransitionEnvelopeCartesian,
                i,feature,featureIdentity,supportId);
            if (adjusted[i].x==original.x && adjusted[i].y==original.y) continue;
            for (auto& polygon:transitionPolygons) {
                for (auto& point:polygon.vertices) {
                    if (point.x==original.x && point.y==original.y) point=adjusted[i];
                }
            }
        }
        // A segment may also graze a Cartesian corner without either sample
        // being near a grid line. Route that local envelope sample through the
        // corner on BOTH sides before cutting, rather than snapping one of
        // two independently built Cut-cell fragments afterwards.
        std::vector<Point2D> routed;
        for (std::size_t i=0;i<adjusted.size();++i) {
            const auto a=adjusted[i];
            const auto b=adjusted[(i+1U)%adjusted.size()];
            routed.push_back(a);
            if (edgeOnRegionBoundary(a,b,outerRegion,tolerance)) continue;
            const auto d=b-a;
            std::vector<std::pair<double,Point2D>> corners;
            const auto xmin=static_cast<long long>(std::ceil((std::min(a.x,b.x)-domain.bounds.min.x)/gridHx));
            const auto xmax=static_cast<long long>(std::floor((std::max(a.x,b.x)-domain.bounds.min.x)/gridHx));
            if (d.x!=0.0) for (auto ix=xmin;ix<=xmax;++ix) {
                const double x=domain.bounds.min.x+static_cast<double>(ix)*gridHx;
                const double t=(x-a.x)/d.x;
                if (!(t>0.0 && t<1.0)) continue;
                const double y=a.y+t*d.y;
                const double gy=domain.bounds.min.y+std::round((y-domain.bounds.min.y)/gridHy)*gridHy;
                if (d.y==0.0) continue;
                const double otherT=(gy-a.y)/d.y;
                if (!(otherT>0.0 && otherT<1.0) ||
                    std::min(std::abs((t-otherT)*d.x),std::abs(y-gy))>=minimumFaceFraction*std::max(gridHx,gridHy)) continue;
                const Point2D corner{x,gy};
                if (segmentLength(a,corner)<minimumFaceFraction*transitionRingThickness ||
                    segmentLength(b,corner)<minimumFaceFraction*transitionRingThickness) continue;
                const Point2D foot=a+d*(dot(corner-a,d)/squaredNorm(d));
                if (segmentLength(foot,corner)>displacementBudget(foot,
                    remainderBoundaryLoops,loopId,i,false,outerRegion)) continue;
                ++supportId;
                (void)envelopeRegistry.addCanonicalVertex(corner,transitionRingThickness,
                    IntersectionFeature2D::CartesianGridVertex,std::nullopt,supportId);
                const auto canonical=envelopeRegistry.canonicalize(foot,transitionRingThickness,
                    IntersectionSource2D::TransitionEnvelopeCartesian,i,
                    IntersectionFeature2D::Smooth,std::nullopt,supportId);
                if (canonical.x!=corner.x || canonical.y!=corner.y) continue;
                corners.push_back({t,corner});
            }
            std::sort(corners.begin(),corners.end(),[](const auto& lhs,const auto& rhs) {
                return lhs.first<rhs.first;
            });
            for (const auto& [t,corner]:corners) {
                (void)t;
                routed.push_back(corner);
            }
            if (!corners.empty()) for (auto& polygon:transitionPolygons) {
                std::vector<Point2D> replacement;
                for (std::size_t j=0;j<polygon.vertices.size();++j) {
                    const auto p=polygon.vertices[j];
                    const auto q=polygon.vertices[(j+1U)%polygon.vertices.size()];
                    replacement.push_back(p);
                    if (p.x==a.x && p.y==a.y && q.x==b.x && q.y==b.y) {
                        for (const auto& item:corners) replacement.push_back(item.second);
                    } else if (p.x==b.x && p.y==b.y && q.x==a.x && q.y==a.y) {
                        for (auto it=corners.rbegin();it!=corners.rend();++it) replacement.push_back(it->second);
                    }
                }
                polygon.vertices=std::move(replacement);
            }
        }
        BoundaryLoop adjustedLoop(std::move(routed));
        if (!adjustedLoop.diagnose(tolerance).valid()) {
            error="feature-safe transition/grid alignment produced an invalid envelope";
            return false;
        }
        remainderBoundaryLoops[loopId]=std::move(adjustedLoop);
    }

    for (const auto& polygon:transitionPolygons) {
        if (!(polygon.signedArea()>0.0) || !BoundaryLoop(polygon.vertices).diagnose(tolerance).valid()) {
            error="canonical transition sampling produced an invalid source polygon";
            return false;
        }
    }
    return true;
}

} // namespace cartmesh2d
