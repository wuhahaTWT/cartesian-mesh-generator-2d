#include "cartmesh2d/geometry/BoundarySimplification2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>

namespace cartmesh2d {
namespace {

[[nodiscard]] double pointSegmentDistance(const Point2D& point,
                                           const Point2D& a,
                                           const Point2D& b) noexcept {
    const Vector2D edge=b-a;
    const double denominator=squaredNorm(edge);
    if (!(denominator>0.0)) return std::sqrt(squaredNorm(point-a));
    const double parameter=std::clamp(dot(point-a,edge)/denominator,0.0,1.0);
    const Point2D projection=a+edge*parameter;
    return std::sqrt(squaredNorm(point-projection));
}

[[nodiscard]] std::vector<Point2D> simplifyOpenChain(
    const std::vector<Point2D>& chain,double maximumDeviation) {
    if (chain.size()<=2) return chain;
    std::vector<bool> keep(chain.size(),false);
    keep.front()=true;
    keep.back()=true;
    std::vector<std::pair<std::size_t,std::size_t>> pending{{0,chain.size()-1}};
    while (!pending.empty()) {
        const auto [first,last]=pending.back();
        pending.pop_back();
        std::optional<std::size_t> split;
        double largest=maximumDeviation;
        for (std::size_t i=first+1;i<last;++i) {
            const double distance=pointSegmentDistance(chain[i],chain[first],chain[last]);
            if (distance>largest) {
                largest=distance;
                split=i;
            }
        }
        if (!split) continue;
        keep[*split]=true;
        pending.emplace_back(*split,last);
        pending.emplace_back(first,*split);
    }
    std::vector<Point2D> result;
    for (std::size_t i=0;i<chain.size();++i) {
        if (keep[i]) result.push_back(chain[i]);
    }
    return result;
}

[[nodiscard]] std::vector<Point2D> forwardChain(
    const std::vector<Point2D>& vertices,std::size_t first,std::size_t last) {
    std::vector<Point2D> chain;
    for (std::size_t i=first;;i=(i+1)%vertices.size()) {
        chain.push_back(vertices[i]);
        if (i==last) break;
    }
    return chain;
}

[[nodiscard]] std::vector<Point2D> simplifyClosedLoop(
    const std::vector<Point2D>& vertices,double maximumDeviation) {
    if (vertices.size()<=3 || !(maximumDeviation>0.0)) return vertices;
    const auto anchorIt=std::min_element(vertices.begin(),vertices.end(),
        [](const Point2D& lhs,const Point2D& rhs) {
            return std::tie(lhs.x,lhs.y)<std::tie(rhs.x,rhs.y);
        });
    const std::size_t first=static_cast<std::size_t>(anchorIt-vertices.begin());
    std::size_t second=first;
    double farthest=-1.0;
    for (std::size_t i=0;i<vertices.size();++i) {
        const double distance=squaredNorm(vertices[i]-vertices[first]);
        if (distance>farthest ||
            (distance==farthest &&
             std::tie(vertices[i].x,vertices[i].y)<
             std::tie(vertices[second].x,vertices[second].y))) {
            farthest=distance;
            second=i;
        }
    }
    if (second==first) return vertices;
    const auto firstHalf=simplifyOpenChain(
        forwardChain(vertices,first,second),maximumDeviation);
    const auto secondHalf=simplifyOpenChain(
        forwardChain(vertices,second,first),maximumDeviation);
    std::vector<Point2D> result=firstHalf;
    if (secondHalf.size()>2) {
        result.insert(result.end(),secondHalf.begin()+1,secondHalf.end()-1);
    }
    return result.size()>=3?result:vertices;
}

[[nodiscard]] std::size_t vertexCount(const BoundaryRegion2D& boundary) {
    return std::accumulate(boundary.loops().begin(),boundary.loops().end(),std::size_t{0},
        [](std::size_t total,const BoundaryLoop& loop) {
            return total+loop.vertices().size();
        });
}

[[nodiscard]] double measuredDeviation(const BoundaryRegion2D& source,
                                       const BoundaryRegion2D& simplified) {
    double maximum=0.0;
    for (std::size_t loop=0;loop<source.loops().size();++loop) {
        const auto& target=simplified.loops()[loop].vertices();
        for (const auto& point:source.loops()[loop].vertices()) {
            double nearest=std::numeric_limits<double>::infinity();
            for (std::size_t i=0;i<target.size();++i) {
                nearest=std::min(nearest,pointSegmentDistance(
                    point,target[i],target[(i+1)%target.size()]));
            }
            maximum=std::max(maximum,nearest);
        }
    }
    return maximum;
}

} // namespace

BoundarySimplificationResult2D simplifyBoundaryRegion2D(
    const BoundaryRegion2D& source,double maximumDeviation,
    const TolerancePolicy& tol) {
    BoundarySimplificationResult2D result;
    result.report.requestedDeviation=maximumDeviation;
    result.report.originalArea=source.area(tol);
    result.report.originalVertexCount=vertexCount(source);
    if (!source.diagnose(tol).valid() || !std::isfinite(maximumDeviation) ||
        maximumDeviation<0.0) {
        result.error="boundary simplification requires a valid source and finite non-negative deviation";
        return result;
    }

    double applied=maximumDeviation;
    const auto sourceDepths=source.nestingDepths(tol);
    for (std::size_t attempt=1;attempt<=32;++attempt) {
        std::vector<BoundaryLoop> loops;
        loops.reserve(source.loops().size());
        for (const auto& loop:source.loops()) {
            loops.emplace_back(simplifyClosedLoop(loop.vertices(),applied));
        }
        BoundaryRegion2D candidate(std::move(loops));
        if (candidate.diagnose(tol).valid() &&
            candidate.nestingDepths(tol)==sourceDepths &&
            candidate.normalizeAlternating(tol)) {
            result.report.appliedDeviation=applied;
            result.report.attempts=attempt;
            result.report.simplifiedArea=candidate.area(tol);
            result.report.simplifiedVertexCount=vertexCount(candidate);
            result.report.measuredMaxDeviation=measuredDeviation(source,candidate);
            if (result.report.measuredMaxDeviation>
                applied+tol.scale(std::max(1.0,applied))) {
                result.error="boundary simplification exceeded its geometric deviation bound";
                return result;
            }
            result.boundary=std::move(candidate);
            return result;
        }
        applied*=0.5;
    }
    result.error="boundary simplification could not preserve loop/region topology";
    return result;
}

} // namespace cartmesh2d
