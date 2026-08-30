#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/quality/PatchLocalQuality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/topology/PatchTransaction2D.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <numbers>
#include <optional>
#include <set>
#include <tuple>

namespace cartmesh2d {
namespace {

using ProfileClock = std::chrono::steady_clock;

[[nodiscard]] double profileSeconds(ProfileClock::time_point start) noexcept {
    return std::chrono::duration<double>(ProfileClock::now()-start).count();
}

[[nodiscard]] bool pointInOrOnTriangle(const Point2D& point,
                                       const Point2D& a,const Point2D& b,
                                       const Point2D& c) noexcept {
    return orientationSign(a,b,point)>=0 && orientationSign(b,c,point)>=0 &&
           orientationSign(c,a,point)>=0;
}

[[nodiscard]] bool strictlyConvex(const Polygon2D& polygon) noexcept {
    if (polygon.vertices.size()<3) return false;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        const auto n=polygon.vertices.size();
        const Point2D& previous=polygon.vertices[(i+n-1)%n];
        const Point2D& current=polygon.vertices[i];
        const Point2D& next=polygon.vertices[(i+1)%n];
        if (orientationSign(previous,current,next)<=0) return false;
        const Vector2D incoming=current-previous;
        const Vector2D outgoing=next-current;
        const double scale=std::sqrt(squaredNorm(incoming)*squaredNorm(outgoing));
        // OpenFOAM's face-plane concavity check treats numerically collinear
        // prism corners as concave even when an exact binary predicate sees a
        // tiny positive turn. Keep such transition vertices out of a single
        // solver cell and partition them explicitly.
        if (!(scale>0.0) || cross(incoming,outgoing)<=1.0e-10*scale) return false;
    }
    return true;
}

[[nodiscard]] bool boundaryEdge(const Point2D& a,const Point2D& b,
                                const Domain2D& domain,
                                const BoundaryRegion2D& boundary,
                                const TolerancePolicy& tol) {
    const auto near=[&](double lhs,double rhs) {
        return std::abs(lhs-rhs)<=tol.scale(std::max(std::abs(lhs),std::abs(rhs)));
    };
    if ((near(a.x,domain.bounds.min.x) && near(b.x,domain.bounds.min.x)) ||
        (near(a.x,domain.bounds.max.x) && near(b.x,domain.bounds.max.x)) ||
        (near(a.y,domain.bounds.min.y) && near(b.y,domain.bounds.min.y)) ||
        (near(a.y,domain.bounds.max.y) && near(b.y,domain.bounds.max.y))) return true;
    for (const auto& loop:boundary.loops()) {
        const auto& vertices=loop.vertices();
        for (std::size_t i=0;i<vertices.size();++i) {
            const Segment2D segment{vertices[i],vertices[(i+1)%vertices.size()]};
            if (pointOnSegment(a,segment,tol) && pointOnSegment(b,segment,tol)) return true;
        }
    }
    return false;
}

[[nodiscard]] unsigned boundaryEdgeCount(const Polygon2D& polygon,
                                         const Domain2D& domain,
                                         const BoundaryRegion2D& boundary,
                                         const TolerancePolicy& tol) {
    unsigned count=0;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        count+=static_cast<unsigned>(boundaryEdge(
            polygon.vertices[i],polygon.vertices[(i+1)%polygon.vertices.size()],
            domain,boundary,tol));
    }
    return count;
}

[[nodiscard]] bool underDeterminedBoundaryCell(const Polygon2D& polygon,
                                               const Domain2D& domain,
                                               const BoundaryRegion2D& boundary,
                                               const TolerancePolicy& tol) {
    return polygon.vertices.size()==3 &&
           boundaryEdgeCount(polygon,domain,boundary,tol)>=2;
}

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygons(
    const Polygon2D& first,const Polygon2D& second) {
    struct DirectedEdge { Point2D a; Point2D b; };
    std::vector<DirectedEdge> remaining;
    const auto add=[&](const Polygon2D& polygon) {
        for (std::size_t i=0;i<polygon.vertices.size();++i) {
            DirectedEdge edge{polygon.vertices[i],polygon.vertices[(i+1)%polygon.vertices.size()]};
            const auto reverse=std::find_if(remaining.begin(),remaining.end(),
                [&](const DirectedEdge& other) {
                    return other.a.x==edge.b.x && other.a.y==edge.b.y &&
                           other.b.x==edge.a.x && other.b.y==edge.a.y;
                });
            if (reverse==remaining.end()) remaining.push_back(edge);
            else remaining.erase(reverse);
        }
    };
    add(first);
    add(second);
    if (remaining.size()<3 ||
        remaining.size()!=first.vertices.size()+second.vertices.size()-2) return std::nullopt;

    const auto firstEdge=std::min_element(remaining.begin(),remaining.end(),
        [](const DirectedEdge& lhs,const DirectedEdge& rhs) {
            return std::tie(lhs.a.x,lhs.a.y,lhs.b.x,lhs.b.y)<
                   std::tie(rhs.a.x,rhs.a.y,rhs.b.x,rhs.b.y);
        });
    const Point2D start=firstEdge->a;
    Point2D next=firstEdge->b;
    std::vector<Point2D> ordered{start};
    remaining.erase(firstEdge);
    while (!remaining.empty()) {
        ordered.push_back(next);
        const auto following=std::find_if(remaining.begin(),remaining.end(),
            [&](const DirectedEdge& edge) {
                return edge.a.x==next.x && edge.a.y==next.y;
            });
        if (following==remaining.end()) return std::nullopt;
        next=following->b;
        remaining.erase(following);
    }
    if (next.x!=start.x || next.y!=start.y) return std::nullopt;
    Polygon2D merged{ordered};
    const double sourceArea=first.area()+second.area();
    const double areaTolerance=1.0e-12+1.0e-10*sourceArea;
    return strictlyConvex(merged) && std::abs(merged.area()-sourceArea)<=areaTolerance
        ?std::optional<Polygon2D>(merged):std::nullopt;
}

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygonsSimple(
    const Polygon2D& first,const Polygon2D& second,const TolerancePolicy& tol) {
    struct DirectedEdge { Point2D a; Point2D b; };
    std::vector<DirectedEdge> remaining;
    const auto same=[](const Point2D& a,const Point2D& b) {
        return a.x==b.x && a.y==b.y;
    };
    const auto add=[&](const Polygon2D& polygon) {
        for (std::size_t i=0;i<polygon.vertices.size();++i) {
            DirectedEdge edge{polygon.vertices[i],
                              polygon.vertices[(i+1)%polygon.vertices.size()]};
            const auto reverse=std::find_if(remaining.begin(),remaining.end(),
                [&](const DirectedEdge& other) {
                    return same(other.a,edge.b) && same(other.b,edge.a);
                });
            if (reverse==remaining.end()) remaining.push_back(edge);
            else remaining.erase(reverse);
        }
    };
    add(first);
    add(second);
    if (remaining.size()<3) return std::nullopt;
    const auto firstEdge=std::min_element(remaining.begin(),remaining.end(),
        [](const DirectedEdge& lhs,const DirectedEdge& rhs) {
            return std::tie(lhs.a.x,lhs.a.y,lhs.b.x,lhs.b.y)<
                   std::tie(rhs.a.x,rhs.a.y,rhs.b.x,rhs.b.y);
        });
    const Point2D start=firstEdge->a;
    Point2D next=firstEdge->b;
    std::vector<Point2D> ordered{start};
    remaining.erase(firstEdge);
    while (!remaining.empty()) {
        ordered.push_back(next);
        const auto following=std::find_if(remaining.begin(),remaining.end(),
            [&](const DirectedEdge& edge) { return same(edge.a,next); });
        if (following==remaining.end()) return std::nullopt;
        next=following->b;
        remaining.erase(following);
    }
    if (!same(next,start)) return std::nullopt;
    Polygon2D merged{ordered};
    const double sourceArea=first.area()+second.area();
    const double areaTolerance=tol.absolute*tol.absolute+tol.relative*sourceArea;
    BoundaryLoop loop(ordered);
    if (!loop.diagnose(tol).valid() || !(merged.signedArea()>0.0) ||
        std::abs(merged.area()-sourceArea)>areaTolerance) return std::nullopt;
    return merged;
}

[[nodiscard]] double minimumInteriorAngle(const Polygon2D& polygon) {
    double minimum=360.0;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        const auto n=polygon.vertices.size();
        const Vector2D a=polygon.vertices[(i+n-1)%n]-polygon.vertices[i];
        const Vector2D b=polygon.vertices[(i+1)%n]-polygon.vertices[i];
        const double denominator=std::sqrt(squaredNorm(a)*squaredNorm(b));
        if (!(denominator>0.0)) return 0.0;
        minimum=std::min(minimum,
            std::acos(std::clamp(dot(a,b)/denominator,-1.0,1.0))*180.0/
            std::numbers::pi);
    }
    return minimum;
}

[[nodiscard]] std::optional<std::vector<Polygon2D>> triangulate(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    if (polygon.vertices.size()<3 || !(polygon.signedArea()>0.0)) return std::nullopt;
    const std::size_t polygonSize=polygon.vertices.size();
    // A polygon with a reflex or collinear transition corner often admits
    // more than one two-piece convex partition. Prefer the most area-balanced
    // valid diagonal; this keeps every conformal boundary vertex and
    // prevents a geometrically legal but solver-hostile sliver beside a much
    // larger cell.  The area identity also rejects diagonals outside the
    // polygon without relying on a floating-point point-in-polygon sample.
    std::optional<std::vector<Polygon2D>> bestConvexSplit;
    unsigned bestConvexBoundaryPenalty=std::numeric_limits<unsigned>::max();
    double bestAreaBalance=-1.0;
    double bestConvexAngle=-1.0;
    for (std::size_t i=0;i<polygonSize;++i) {
        for (std::size_t j=i+1;j<polygonSize;++j) {
            if (j==i+1 || (i==0 && j+1==polygonSize)) continue;
            Polygon2D first,second;
            for (std::size_t k=i;;k=(k+1)%polygonSize) {
                first.vertices.push_back(polygon.vertices[k]);
                if (k==j) break;
            }
            for (std::size_t k=j;;k=(k+1)%polygonSize) {
                second.vertices.push_back(polygon.vertices[k]);
                if (k==i) break;
            }
            if (!strictlyConvex(first) || !strictlyConvex(second)) continue;
            const double firstArea=first.area();
            const double secondArea=second.area();
            const double areaError=std::abs(firstArea+secondArea-polygon.area());
            if (areaError>tol.absolute*tol.absolute+tol.relative*polygon.area()) continue;
            const unsigned penalty=
                static_cast<unsigned>(underDeterminedBoundaryCell(first,domain,boundary,tol))+
                static_cast<unsigned>(underDeterminedBoundaryCell(second,domain,boundary,tol));
            const double balance=std::min(firstArea,secondArea)/polygon.area();
            const double angle=std::min(minimumInteriorAngle(first),
                                        minimumInteriorAngle(second));
            if (!bestConvexSplit || penalty<bestConvexBoundaryPenalty ||
                (penalty==bestConvexBoundaryPenalty && balance>bestAreaBalance) ||
                (penalty==bestConvexBoundaryPenalty && balance==bestAreaBalance &&
                 angle>bestConvexAngle)) {
                bestConvexSplit=std::vector<Polygon2D>{first,second};
                bestConvexBoundaryPenalty=penalty;
                bestAreaBalance=balance;
                bestConvexAngle=angle;
            }
        }
    }
    if (bestConvexSplit) return bestConvexSplit;

    std::vector<std::size_t> remaining(polygon.vertices.size());
    std::iota(remaining.begin(),remaining.end(),0);
    std::vector<Polygon2D> triangles;
    triangles.reserve(polygon.vertices.size()-2);
    while (remaining.size()>3) {
        std::optional<std::size_t> selected;
        unsigned selectedBoundaryPenalty=2;
        double selectedArea=-1.0;
        double selectedAngle=-1.0;
        for (std::size_t pos=0;pos<remaining.size();++pos) {
            const std::size_t prev=remaining[(pos+remaining.size()-1)%remaining.size()];
            const std::size_t current=remaining[pos];
            const std::size_t next=remaining[(pos+1)%remaining.size()];
            const Point2D& a=polygon.vertices[prev];
            const Point2D& b=polygon.vertices[current];
            const Point2D& c=polygon.vertices[next];
            if (orientationSign(a,b,c)<=0) continue;
            bool containsOther=false;
            for (const auto vertex:remaining) {
                if (vertex==prev || vertex==current || vertex==next) continue;
                if (pointInOrOnTriangle(polygon.vertices[vertex],a,b,c)) {
                    containsOther=true;
                    break;
                }
            }
            if (containsOther) continue;
            const unsigned boundaryEdges=
                static_cast<unsigned>(boundaryEdge(a,b,domain,boundary,tol))+
                static_cast<unsigned>(boundaryEdge(b,c,domain,boundary,tol));
            const unsigned boundaryPenalty=static_cast<unsigned>(boundaryEdges>=2);
            const Polygon2D ear{{a,b,c}};
            const double earArea=ear.area();
            const double earAngle=minimumInteriorAngle(ear);
            if (!selected || boundaryPenalty<selectedBoundaryPenalty ||
                (boundaryPenalty==selectedBoundaryPenalty && earArea>selectedArea) ||
                (boundaryPenalty==selectedBoundaryPenalty && earArea==selectedArea &&
                 earAngle>selectedAngle)) {
                selected=pos;
                selectedBoundaryPenalty=boundaryPenalty;
                selectedArea=earArea;
                selectedAngle=earAngle;
            }
        }
        if (!selected) return std::nullopt;
        const std::size_t pos=*selected;
        const std::size_t prev=remaining[(pos+remaining.size()-1)%remaining.size()];
        const std::size_t current=remaining[pos];
        const std::size_t next=remaining[(pos+1)%remaining.size()];
        triangles.push_back({{polygon.vertices[prev],polygon.vertices[current],
                               polygon.vertices[next]}});
        remaining.erase(remaining.begin()+static_cast<std::ptrdiff_t>(pos));
    }
    triangles.push_back({{polygon.vertices[remaining[0]],polygon.vertices[remaining[1]],
                           polygon.vertices[remaining[2]]}});
    // A triangular prism at a sharp physical corner can have only one
    // non-boundary side face. OpenFOAM correctly reports such a cell as
    // under-determined. Merge it with the first adjacent triangle whenever
    // their union is a strictly convex quadrilateral.
    for (std::size_t i=0;i<triangles.size();++i) {
        if (boundaryEdgeCount(triangles[i],domain,boundary,tol)<2) continue;
        for (std::size_t j=0;j<triangles.size();++j) {
            if (i==j) continue;
            unsigned shared=0;
            for (const auto& a:triangles[i].vertices) {
                for (const auto& b:triangles[j].vertices) {
                    if (a.x==b.x && a.y==b.y) ++shared;
                }
            }
            if (shared!=2) continue;
            const auto merged=mergeAdjacentPolygons(triangles[i],triangles[j]);
            if (!merged) continue;
            triangles[i]=*merged;
            triangles.erase(triangles.begin()+static_cast<std::ptrdiff_t>(j));
            if (j<i) --i;
            break;
        }
    }
    // Remove poor artificial diagonals introduced by ear clipping.  Among all
    // adjacent triangle pairs, deterministically merge the pair whose convex
    // quadrilateral gives the largest increase in minimum interior angle.
    // Physical boundary vertices and total area are unchanged.
    while (true) {
        std::optional<std::pair<std::size_t,std::size_t>> bestPair;
        std::optional<Polygon2D> bestPolygon;
        double bestGain=1.0e-12;
        for (std::size_t i=0;i<triangles.size();++i) {
            for (std::size_t j=i+1;j<triangles.size();++j) {
                unsigned shared=0;
                for (const auto& a:triangles[i].vertices) {
                    for (const auto& b:triangles[j].vertices) {
                        if (a.x==b.x && a.y==b.y) ++shared;
                    }
                }
                if (shared!=2) continue;
                const auto merged=mergeAdjacentPolygons(triangles[i],triangles[j]);
                if (!merged) continue;
                const double before=std::min(minimumInteriorAngle(triangles[i]),
                                             minimumInteriorAngle(triangles[j]));
                const double gain=minimumInteriorAngle(*merged)-before;
                if (gain>bestGain) {
                    bestGain=gain;
                    bestPair={{i,j}};
                    bestPolygon=merged;
                }
            }
        }
        if (!bestPair) break;
        triangles[bestPair->first]=*bestPolygon;
        triangles.erase(triangles.begin()+static_cast<std::ptrdiff_t>(bestPair->second));
    }
    return triangles;
}

[[nodiscard]] std::vector<std::size_t> mergedLineage(
    const std::vector<std::size_t>& first,
    const std::vector<std::size_t>& second) {
    std::vector<std::size_t> result=first;
    result.insert(result.end(),second.begin(),second.end());
    std::sort(result.begin(),result.end());
    result.erase(std::unique(result.begin(),result.end()),result.end());
    return result;
}

[[nodiscard]] Polygon2D removeArtificialCollinearVertices(
    Polygon2D polygon,const TolerancePolicy& tol) {
    bool changed=true;
    while (changed && polygon.vertices.size()>3U) {
        changed=false;
        for (std::size_t i=0;i<polygon.vertices.size();++i) {
            const auto n=polygon.vertices.size();
            const auto& previous=polygon.vertices[(i+n-1U)%n];
            const auto& current=polygon.vertices[i];
            const auto& next=polygon.vertices[(i+1U)%n];
            const Vector2D incoming=current-previous;
            const Vector2D outgoing=next-current;
            const double scale=std::sqrt(squaredNorm(incoming)*squaredNorm(outgoing));
            if (!(scale>0.0) || dot(incoming,outgoing)<0.0 ||
                std::abs(cross(incoming,outgoing))>
                    64.0*std::numeric_limits<double>::epsilon()*scale) continue;
            auto candidate=polygon;
            candidate.vertices.erase(candidate.vertices.begin()+
                static_cast<std::ptrdiff_t>(i));
            const double areaScale=std::max(polygon.area(),candidate.area());
            if (std::abs(candidate.area()-polygon.area())>
                tol.absolute*tol.absolute+tol.relative*areaScale) continue;
            polygon=std::move(candidate);
            changed=true;
            break;
        }
    }
    return polygon;
}

[[nodiscard]] CutCell2D makeCell(std::size_t id,const Polygon2D& polygon,
                                 const TolerancePolicy& tol,
                                 std::vector<std::size_t> sourceLineage={}) {
    CutCell2D cell;
    cell.sourceId=id;
    cell.sourceKey=id;
    cell.sourceLineage=std::move(sourceLineage);
    cell.backgroundBounds=polygon.bounds();
    cell.kind=CutCellKind::Cut;
    cell.fluidPolygon=polygon;
    cell.area=polygon.area();
    const double boxArea=(cell.backgroundBounds.max.x-cell.backgroundBounds.min.x)*
                         (cell.backgroundBounds.max.y-cell.backgroundBounds.min.y);
    cell.areaFraction=boxArea>0.0?cell.area/boxArea:0.0;
    cell.centroid=polygon.centroid(tol);
    return cell;
}

struct PartitionAttempt2D {
    TopologyMesh2D topology;
    std::vector<std::size_t> sourceForCell;
    std::vector<bool> immutableForCell;
    std::size_t partitionedSourceCellCount = 0;
    std::string error;
};

struct QualityScore2D {
    std::size_t issueCount = 0;
    double maximumSeverity = 0.0;
    double totalSeverity = 0.0;
};

struct LocalQualityRank2D {
    QualityScore2D issues;
    double maxNonOrthogonality = 0.0;
    double maxInternalSkewness = 0.0;
    double maxBoundarySkewness = 0.0;
    double maxCellAspect = 0.0;
    double negativeMinInteriorAngle = 0.0;
    double negativeMinFaceWeight = 0.0;
    double negativeMinVolumeRatio = 0.0;
};

[[nodiscard]] QualityScore2D qualityScore(const SolverQualityReport2D& quality) {
    QualityScore2D score;
    score.issueCount=quality.issues.size();
    for (const auto& issue:quality.issues) {
        const bool lowerBound=
            issue.code==SolverQualityIssueCode2D::ShortFace ||
            issue.code==SolverQualityIssueCode2D::SmallInteriorAngle ||
            issue.code==SolverQualityIssueCode2D::LowFaceWeight ||
            issue.code==SolverQualityIssueCode2D::LowVolumeRatio;
        const double severity=lowerBound
            ?issue.limit/(issue.measured+std::numeric_limits<double>::min())
            :issue.measured/(issue.limit+std::numeric_limits<double>::min());
        score.maximumSeverity=std::max(score.maximumSeverity,severity);
        score.totalSeverity+=severity;
    }
    return score;
}

[[nodiscard]] bool betterQualityScore(const QualityScore2D& candidate,
                                      const QualityScore2D& current) noexcept {
    return std::tie(candidate.issueCount,candidate.maximumSeverity,candidate.totalSeverity)<
           std::tie(current.issueCount,current.maximumSeverity,current.totalSeverity);
}

[[nodiscard]] bool betterLocalQualityRank(const LocalQualityRank2D& candidate,
                                          const LocalQualityRank2D& current) noexcept {
    return std::tie(candidate.issues.issueCount,
                    candidate.issues.maximumSeverity,
                    candidate.issues.totalSeverity,
                    candidate.maxNonOrthogonality,
                    candidate.maxInternalSkewness,
                    candidate.maxBoundarySkewness,
                    candidate.maxCellAspect,
                    candidate.negativeMinInteriorAngle,
                    candidate.negativeMinFaceWeight,
                    candidate.negativeMinVolumeRatio)<
           std::tie(current.issues.issueCount,
                    current.issues.maximumSeverity,
                    current.issues.totalSeverity,
                    current.maxNonOrthogonality,
                    current.maxInternalSkewness,
                    current.maxBoundarySkewness,
                    current.maxCellAspect,
                    current.negativeMinInteriorAngle,
                    current.negativeMinFaceWeight,
                    current.negativeMinVolumeRatio);
}

[[nodiscard]] std::optional<std::vector<Polygon2D>> partitionOnePolygon(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    if (strictlyConvex(polygon)) return std::vector<Polygon2D>{polygon};
    return triangulate(polygon,domain,boundary,tol);
}

[[nodiscard]] std::optional<LocalQualityRank2D> localPartitionQualityRank(
    const Polygon2D& outer,const std::vector<Polygon2D>& pieces,
    const TolerancePolicy& tol) {
    if (pieces.empty() || !(outer.signedArea()>0.0)) return std::nullopt;
    BoundaryLoop localLoop(outer.vertices);
    BoundaryRegion2D localBoundary(localLoop);
    const Domain2D localDomain{outer.bounds()};
    if (!localBoundary.diagnose(tol).valid() || !localDomain.valid(tol)) return std::nullopt;
    std::vector<CutCell2D> cells;
    cells.reserve(pieces.size());
    for (std::size_t i=0;i<pieces.size();++i) cells.push_back(makeCell(i,pieces[i],tol));
    const auto topology=buildGlobalTopology(cells,localDomain,localBoundary,tol);
    if (!topology.valid()) return std::nullopt;
    const auto quality=evaluateSolverQuality2D(topology,{},tol);
    SolverQualityReport2D filtered;
    filtered.policy=quality.policy;
    for (const auto& issue:quality.issues) {
        const bool cellIssue=issue.code==SolverQualityIssueCode2D::InvalidCell ||
            issue.code==SolverQualityIssueCode2D::ExcessiveConcavity ||
            issue.code==SolverQualityIssueCode2D::ExcessiveAspect ||
            issue.code==SolverQualityIssueCode2D::SmallInteriorAngle;
        const bool internalEdgeIssue=issue.edgeId<topology.edges.size() &&
            topology.edges[issue.edgeId].neighbour.has_value();
        if (cellIssue || internalEdgeIssue) filtered.issues.push_back(issue);
    }
    LocalQualityRank2D rank;
    rank.issues=qualityScore(filtered);
    rank.maxNonOrthogonality=quality.maxNonOrthogonalityDeg;
    rank.maxInternalSkewness=quality.maxInternalSkewness;
    rank.maxCellAspect=quality.maxCellAspect;
    rank.negativeMinInteriorAngle=-quality.minInteriorAngleDeg;
    rank.negativeMinFaceWeight=-quality.minFaceWeight;
    rank.negativeMinVolumeRatio=-quality.minVolumeRatio;
    return rank;
}

[[nodiscard]] bool internalInterfaceIssue(SolverQualityIssueCode2D code) noexcept {
    return code==SolverQualityIssueCode2D::ExcessiveNonOrthogonality ||
           code==SolverQualityIssueCode2D::ExcessiveSkewness ||
           code==SolverQualityIssueCode2D::LowFaceWeight ||
           code==SolverQualityIssueCode2D::LowVolumeRatio;
}

[[nodiscard]] PartitionAttempt2D partitionSourcePolygons(
    const std::vector<Polygon2D>& sourcePolygons,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile=nullptr,bool candidate=false,
    const std::vector<bool>& immutableSources={},
    const std::vector<bool>& preserveSources={},
    const std::vector<std::vector<std::size_t>>& sourceLineages={},
    std::shared_ptr<IntersectionRegistry2D> registry={}) {
    PartitionAttempt2D result;
    std::vector<CutCell2D> cells;
    cells.reserve(sourcePolygons.size());
    if (!sourceLineages.empty() && sourceLineages.size()!=sourcePolygons.size()) {
        result.error="solver source lineage count does not match source polygons";
        return result;
    }
    for (std::size_t source=0;source<sourcePolygons.size();++source) {
        const auto& polygon=sourcePolygons[source];
        const bool immutable=!immutableSources.empty() && immutableSources[source];
        const bool preserve=!preserveSources.empty() && preserveSources[source];
        if (immutable || preserve || strictlyConvex(polygon)) {
            cells.push_back(makeCell(
                cells.size(),polygon,tol,
                sourceLineages.empty()?std::vector<std::size_t>{source}:sourceLineages[source]));
            result.sourceForCell.push_back(source);
            result.immutableForCell.push_back(immutable);
            continue;
        }
        const auto pieces=triangulate(polygon,domain,boundary,tol);
        if (!pieces) {
            result.error="deterministic convex partition failed for source cell "+
                         std::to_string(source);
            return result;
        }
        ++result.partitionedSourceCellCount;
        for (const auto& piece:*pieces) {
            cells.push_back(makeCell(
                cells.size(),piece,tol,
                sourceLineages.empty()?std::vector<std::size_t>{source}:sourceLineages[source]));
            result.sourceForCell.push_back(source);
            result.immutableForCell.push_back(false);
        }
    }
    const auto topologyStart=ProfileClock::now();
    result.topology=buildGlobalTopology(cells,domain,boundary,tol,registry);
    if (profile) {
        const double elapsed=profileSeconds(topologyStart);
        profile->buildGlobalTopologySeconds+=elapsed;
        ++profile->globalTopologyRebuildCalls;
        if (candidate) {
            profile->candidateGlobalRebuildSeconds+=elapsed;
            ++profile->candidateTopologyCount;
        }
    }
    if (!result.topology.valid()) {
        result.error="partitioned solver topology audit failed";
        if (!result.topology.issues.empty()) {
            result.error+=": "+result.topology.issues.front().message;
        }
    }
    return result;
}

[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> sourceMergePairs(
    const PartitionAttempt2D& partition,const SolverQualityReport2D& quality,
    const std::vector<bool>& immutableSources) {
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    const auto addEdge=[&](const Edge2D& edge) {
        if (!edge.neighbour || edge.owner>=partition.sourceForCell.size() ||
            *edge.neighbour>=partition.sourceForCell.size()) return;
        const std::size_t first=partition.sourceForCell[edge.owner];
        const std::size_t second=partition.sourceForCell[*edge.neighbour];
        if ((!immutableSources.empty() && immutableSources[first]) ||
            (!immutableSources.empty() && immutableSources[second])) return;
        if (first!=second) pairs.emplace_back(std::min(first,second),std::max(first,second));
    };
    for (const auto& issue:quality.issues) {
        if (!internalInterfaceIssue(issue.code)) continue;
        if (issue.edgeId<partition.topology.edges.size()) {
            addEdge(partition.topology.edges[issue.edgeId]);
        }
        if (issue.cellId>=partition.sourceForCell.size() ||
            issue.cellId>=partition.topology.cells.size()) continue;
        for (const auto edgeId:partition.topology.cells[issue.cellId].edges) {
            if (edgeId<partition.topology.edges.size()) {
                addEdge(partition.topology.edges[edgeId]);
            }
        }
    }
    std::sort(pairs.begin(),pairs.end());
    pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    return pairs;
}

[[nodiscard]] Polygon2D topologyCellPolygon(const TopologyMesh2D& topology,
                                            std::size_t cellId) {
    Polygon2D polygon;
    if (cellId>=topology.cells.size()) return polygon;
    for (const auto vertex:topology.cells[cellId].vertices) {
        if (vertex>=topology.vertices.size()) return {};
        polygon.vertices.push_back(topology.vertices[vertex].point);
    }
    return polygon;
}

struct RepartitionBatch2D {
    TopologyMesh2D topology;
    std::vector<bool> immutableCells;
};

[[nodiscard]] std::vector<std::pair<Polygon2D,Polygon2D>> convexTwoPieceSplits(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    std::vector<std::pair<Polygon2D,Polygon2D>> splits;
    const std::size_t n=polygon.vertices.size();
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=i+1;j<n;++j) {
            if (j==i+1 || (i==0 && j+1==n)) continue;
            Polygon2D first,second;
            for (std::size_t k=i;;k=(k+1)%n) {
                first.vertices.push_back(polygon.vertices[k]);
                if (k==j) break;
            }
            for (std::size_t k=j;;k=(k+1)%n) {
                second.vertices.push_back(polygon.vertices[k]);
                if (k==i) break;
            }
            if (!strictlyConvex(first) || !strictlyConvex(second) ||
                underDeterminedBoundaryCell(first,domain,boundary,tol) ||
                underDeterminedBoundaryCell(second,domain,boundary,tol)) continue;
            const double areaError=std::abs(first.area()+second.area()-polygon.area());
            if (areaError>tol.absolute*tol.absolute+tol.relative*polygon.area()) continue;
            splits.emplace_back(std::move(first),std::move(second));
        }
    }
    return splits;
}

[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>>
solverRepartitionPairs(const TopologyMesh2D& topology,
                       const SolverQualityReport2D& quality,
                       const std::vector<bool>& immutableCells={}) {
    std::vector<std::size_t> affected;
    for (const auto& issue:quality.issues) {
        if (!internalInterfaceIssue(issue.code)) continue;
        if (issue.cellId<topology.cells.size()) affected.push_back(issue.cellId);
        if (issue.edgeId<topology.edges.size() && topology.edges[issue.edgeId].neighbour)
            affected.push_back(*topology.edges[issue.edgeId].neighbour);
    }
    std::sort(affected.begin(),affected.end());
    affected.erase(std::unique(affected.begin(),affected.end()),affected.end());
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    for (const auto cell:affected) {
        for (const auto edgeId:topology.cells[cell].edges) {
            if (edgeId>=topology.edges.size()) continue;
            const auto& edge=topology.edges[edgeId];
            if (!edge.neighbour) continue;
            const std::size_t other=edge.owner==cell?*edge.neighbour:edge.owner;
            if ((!immutableCells.empty() && immutableCells[cell]) ||
                (!immutableCells.empty() && immutableCells[other])) continue;
            pairs.emplace_back(std::min(cell,other),std::max(cell,other));
        }
    }
    std::sort(pairs.begin(),pairs.end());
    pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    return pairs;
}

[[nodiscard]] RepartitionBatch2D repartitionPair(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol,SolverTopologyProfile2D* profile,
    const std::vector<bool>& immutableCells) {
    std::vector<CutCell2D> cells;
    std::vector<bool> rebuiltImmutable;
    cells.reserve(topology.cells.size());
    const auto replacementLineage=mergedLineage(
        topology.cells[first].sourceLineage,topology.cells[second].sourceLineage);
    for (std::size_t cell=0;cell<topology.cells.size();++cell) {
        if (cell==first) {
            cells.push_back(makeCell(cells.size(),firstPiece,tol,replacementLineage));
            cells.push_back(makeCell(cells.size(),secondPiece,tol,replacementLineage));
            rebuiltImmutable.push_back(false);
            rebuiltImmutable.push_back(false);
        }
        if (cell==first || cell==second) continue;
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol,
                                 topology.cells[cell].sourceLineage));
        rebuiltImmutable.push_back(!immutableCells.empty() && immutableCells[cell]);
    }
    const auto topologyStart=ProfileClock::now();
    RepartitionBatch2D rebuilt;
    rebuilt.topology=buildGlobalTopology(cells,domain,boundary,tol,topology.constructionRegistry);
    rebuilt.immutableCells=std::move(rebuiltImmutable);
    if (profile) {
        const double elapsed=profileSeconds(topologyStart);
        profile->buildGlobalTopologySeconds+=elapsed;
        profile->candidateGlobalRebuildSeconds+=elapsed;
        ++profile->globalTopologyRebuildCalls;
        ++profile->candidateTopologyCount;
    }
    return rebuilt;
}

[[nodiscard]] RepartitionBatch2D agglomerateCellPair(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second,
    const Polygon2D& merged,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,const std::vector<bool>& immutableCells) {
    std::vector<CutCell2D> cells;
    std::vector<bool> rebuiltImmutable;
    cells.reserve(topology.cells.size()-1U);
    rebuiltImmutable.reserve(topology.cells.size()-1U);
    const auto replacementLineage=mergedLineage(
        topology.cells[first].sourceLineage,topology.cells[second].sourceLineage);
    for (std::size_t cell=0;cell<topology.cells.size();++cell) {
        if (cell==first) {
            cells.push_back(makeCell(cells.size(),merged,tol,replacementLineage));
            rebuiltImmutable.push_back(false);
        }
        if (cell==first || cell==second) continue;
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol,
                                 topology.cells[cell].sourceLineage));
        rebuiltImmutable.push_back(!immutableCells.empty() && immutableCells[cell]);
    }
    const auto topologyStart=ProfileClock::now();
    RepartitionBatch2D rebuilt;
    rebuilt.topology=buildGlobalTopology(cells,domain,boundary,tol,topology.constructionRegistry);
    rebuilt.immutableCells=std::move(rebuiltImmutable);
    if (profile) {
        const double elapsed=profileSeconds(topologyStart);
        profile->buildGlobalTopologySeconds+=elapsed;
        profile->candidateGlobalRebuildSeconds+=elapsed;
        ++profile->globalTopologyRebuildCalls;
        ++profile->candidateTopologyCount;
    }
    return rebuilt;
}

struct SourceMergeProposal2D {
    std::size_t first = 0;
    std::size_t second = 0;
    Polygon2D merged;
    LocalQualityRank2D rank;
    std::vector<std::size_t> halo;
};

struct RepartitionProposal2D {
    std::size_t first = 0;
    std::size_t second = 0;
    Polygon2D firstPiece;
    Polygon2D secondPiece;
    LocalQualityRank2D rank;
    std::vector<std::size_t> halo;
};

[[nodiscard]] bool sortedSetsIntersect(const std::vector<std::size_t>& first,
                                       const std::vector<std::size_t>& second) noexcept {
    std::size_t i=0,j=0;
    while (i<first.size() && j<second.size()) {
        if (first[i]==second[j]) return true;
        if (first[i]<second[j]) ++i;
        else ++j;
    }
    return false;
}

[[nodiscard]] std::vector<std::vector<std::size_t>> sourceAdjacency(
    const PartitionAttempt2D& partition,std::size_t sourceCount) {
    std::vector<std::vector<std::size_t>> adjacency(sourceCount);
    for (const auto& edge:partition.topology.edges) {
        if (!edge.neighbour || edge.owner>=partition.sourceForCell.size() ||
            *edge.neighbour>=partition.sourceForCell.size()) continue;
        const auto first=partition.sourceForCell[edge.owner];
        const auto second=partition.sourceForCell[*edge.neighbour];
        if (first==second || first>=sourceCount || second>=sourceCount) continue;
        adjacency[first].push_back(second);
        adjacency[second].push_back(first);
    }
    for (auto& neighbours:adjacency) {
        std::sort(neighbours.begin(),neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
    }
    return adjacency;
}

[[nodiscard]] std::vector<std::size_t> sourcePairHalo(
    std::size_t first,std::size_t second,
    const std::vector<std::vector<std::size_t>>& adjacency) {
    std::vector<std::size_t> halo{first,second};
    if (first<adjacency.size()) halo.insert(halo.end(),adjacency[first].begin(),adjacency[first].end());
    if (second<adjacency.size()) halo.insert(halo.end(),adjacency[second].begin(),adjacency[second].end());
    std::sort(halo.begin(),halo.end());
    halo.erase(std::unique(halo.begin(),halo.end()),halo.end());
    return halo;
}

[[nodiscard]] std::vector<std::size_t> cellPairHalo(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second) {
    std::vector<std::size_t> halo{first,second};
    const auto addNeighbours=[&](std::size_t cell) {
        if (cell>=topology.cells.size()) return;
        for (const auto edgeId:topology.cells[cell].edges) {
            if (edgeId>=topology.edges.size()) continue;
            const auto& edge=topology.edges[edgeId];
            if (edge.owner!=cell) halo.push_back(edge.owner);
            if (edge.neighbour && *edge.neighbour!=cell) halo.push_back(*edge.neighbour);
        }
    };
    addNeighbours(first);
    addNeighbours(second);
    std::sort(halo.begin(),halo.end());
    halo.erase(std::unique(halo.begin(),halo.end()),halo.end());
    return halo;
}

[[nodiscard]] std::optional<BoundaryRegion2D> patchBoundaryRegion(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& cells,
    const TolerancePolicy& tol) {
    const auto inside=[&](std::size_t cell) {
        return std::binary_search(cells.begin(),cells.end(),cell);
    };
    std::map<std::size_t,std::vector<std::size_t>> adjacency;
    for (const auto& edge:topology.edges) {
        const bool ownerInside=inside(edge.owner);
        const bool neighbourInside=edge.neighbour && inside(*edge.neighbour);
        if (ownerInside==neighbourInside) continue;
        adjacency[edge.v0].push_back(edge.v1);
        adjacency[edge.v1].push_back(edge.v0);
    }
    for (auto& [vertex,neighbours]:adjacency) {
        (void)vertex;
        std::sort(neighbours.begin(),neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
        if (neighbours.size()!=2) return std::nullopt;
    }
    std::map<std::pair<std::size_t,std::size_t>,bool> used;
    for (const auto& [vertex,neighbours]:adjacency) {
        for (const auto neighbour:neighbours) {
            used[{std::min(vertex,neighbour),std::max(vertex,neighbour)}]=false;
        }
    }
    std::vector<BoundaryLoop> loops;
    for (auto& [initialEdge,initialUsed]:used) {
        if (initialUsed) continue;
        const std::size_t start=initialEdge.first;
        std::size_t previous=start;
        std::size_t current=initialEdge.second;
        initialUsed=true;
        std::vector<Point2D> points{topology.vertices[start].point};
        for (std::size_t guard=0;guard<=used.size();++guard) {
            if (current==start) break;
            if (current>=topology.vertices.size()) return std::nullopt;
            points.push_back(topology.vertices[current].point);
            const auto found=adjacency.find(current);
            if (found==adjacency.end()) return std::nullopt;
            const std::size_t next=found->second[0]==previous
                ?found->second[1]:found->second[0];
            const auto key=std::make_pair(std::min(current,next),std::max(current,next));
            auto edgeUsed=used.find(key);
            if (edgeUsed==used.end() || edgeUsed->second) return std::nullopt;
            edgeUsed->second=true;
            previous=current;
            current=next;
        }
        if (current!=start || points.size()<3) return std::nullopt;
        BoundaryLoop loop(std::move(points));
        if (!loop.diagnose(tol).valid()) return std::nullopt;
        loops.push_back(std::move(loop));
    }
    if (loops.empty()) return std::nullopt;
    BoundaryRegion2D boundary(std::move(loops));
    if (!boundary.normalizeAlternating(tol) || !boundary.diagnose(tol).valid()) {
        return std::nullopt;
    }
    return boundary;
}

[[nodiscard]] std::optional<LocalQualityRank2D> localRepartitionQualityRank(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& halo,
    std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const TolerancePolicy& tol) {
    const auto boundary=patchBoundaryRegion(topology,halo,tol);
    if (!boundary) return std::nullopt;
    const Domain2D domain{boundary->bounds()};
    if (!domain.valid(tol)) return std::nullopt;
    std::vector<CutCell2D> cells;
    cells.reserve(halo.size());
    for (const auto cell:halo) {
        if (cell==first) {
            cells.push_back(makeCell(cells.size(),firstPiece,tol));
            cells.push_back(makeCell(cells.size(),secondPiece,tol));
        }
        if (cell==first || cell==second) continue;
        if (cell>=topology.cells.size()) return std::nullopt;
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol));
    }
    const auto patch=buildGlobalTopology(cells,domain,*boundary,tol);
    if (!patch.valid()) return std::nullopt;
    const auto quality=evaluateSolverQuality2D(patch,{},tol);
    SolverQualityReport2D filtered;
    filtered.policy=quality.policy;
    double physicalBoundarySkewness=0.0;
    const auto isPhysicalBoundaryEdge=[&](const Edge2D& patchEdge) {
        const auto& patchA=patch.vertices[patchEdge.v0].point;
        const auto& patchB=patch.vertices[patchEdge.v1].point;
        for (const auto cell:halo) {
            for (const auto edgeId:topology.cells[cell].edges) {
                const auto& edge=topology.edges[edgeId];
                if (edge.neighbour) continue;
                const auto& a=topology.vertices[edge.v0].point;
                const auto& b=topology.vertices[edge.v1].point;
                if ((a.x==patchA.x && a.y==patchA.y &&
                     b.x==patchB.x && b.y==patchB.y) ||
                    (a.x==patchB.x && a.y==patchB.y &&
                     b.x==patchA.x && b.y==patchA.y)) return true;
            }
        }
        return false;
    };
    for (const auto& issue:quality.issues) {
        if (issue.code!=SolverQualityIssueCode2D::ExcessiveBoundarySkewness ||
            (issue.edgeId<patch.edges.size() && isPhysicalBoundaryEdge(patch.edges[issue.edgeId]))) {
            filtered.issues.push_back(issue);
            if (issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness) {
                physicalBoundarySkewness=std::max(physicalBoundarySkewness,issue.measured);
            }
        }
    }
    LocalQualityRank2D rank;
    rank.issues=qualityScore(filtered);
    rank.maxNonOrthogonality=quality.maxNonOrthogonalityDeg;
    rank.maxInternalSkewness=quality.maxInternalSkewness;
    rank.maxBoundarySkewness=physicalBoundarySkewness;
    rank.maxCellAspect=quality.maxCellAspect;
    rank.negativeMinInteriorAngle=-quality.minInteriorAngleDeg;
    rank.negativeMinFaceWeight=-quality.minFaceWeight;
    rank.negativeMinVolumeRatio=-quality.minVolumeRatio;
    return rank;
}

template<class Proposal>
[[nodiscard]] std::vector<Proposal> selectIndependentProposals(
    std::vector<Proposal> proposals) {
    std::sort(proposals.begin(),proposals.end(),[](const Proposal& first,const Proposal& second) {
        if (betterLocalQualityRank(first.rank,second.rank)) return true;
        if (betterLocalQualityRank(second.rank,first.rank)) return false;
        return std::tie(first.first,first.second)<std::tie(second.first,second.second);
    });
    std::vector<std::vector<std::size_t>> halos;
    halos.reserve(proposals.size());
    for (const auto& proposal:proposals) halos.push_back(proposal.halo);
    const auto indices=selectIndependentSolverRepairPatches2D(halos);
    std::vector<Proposal> selected;
    selected.reserve(indices.size());
    for (const auto index:indices) selected.push_back(std::move(proposals[index]));
    return selected;
}

[[nodiscard]] std::vector<Polygon2D> applySourceMergeBatch(
    const std::vector<Polygon2D>& sourcePolygons,
    const std::vector<SourceMergeProposal2D>& selected) {
    std::vector<bool> removed(sourcePolygons.size(),false);
    std::vector<const Polygon2D*> replacements(sourcePolygons.size(),nullptr);
    for (const auto& proposal:selected) {
        if (proposal.first>=sourcePolygons.size() || proposal.second>=sourcePolygons.size()) continue;
        replacements[proposal.first]=&proposal.merged;
        removed[proposal.second]=true;
    }
    std::vector<Polygon2D> result;
    result.reserve(sourcePolygons.size()-selected.size());
    for (std::size_t source=0;source<sourcePolygons.size();++source) {
        if (removed[source]) continue;
        result.push_back(replacements[source]?*replacements[source]:sourcePolygons[source]);
    }
    return result;
}

[[nodiscard]] std::vector<bool> applySourceMergeProtection(
    const std::vector<bool>& immutableSources,
    const std::vector<SourceMergeProposal2D>& selected) {
    if (immutableSources.empty()) return {};
    std::vector<bool> removed(immutableSources.size(),false);
    for (const auto& proposal:selected) {
        if (proposal.second<removed.size()) removed[proposal.second]=true;
    }
    std::vector<bool> result;
    result.reserve(immutableSources.size()-selected.size());
    for (std::size_t source=0;source<immutableSources.size();++source) {
        if (!removed[source]) result.push_back(immutableSources[source]);
    }
    return result;
}

[[nodiscard]] std::vector<std::vector<std::size_t>> applySourceMergeLineage(
    const std::vector<std::vector<std::size_t>>& sourceLineages,
    const std::vector<SourceMergeProposal2D>& selected) {
    std::vector<bool> removed(sourceLineages.size(),false);
    auto merged=sourceLineages;
    for (const auto& proposal:selected) {
        if (proposal.first>=merged.size() || proposal.second>=merged.size()) continue;
        merged[proposal.first]=mergedLineage(
            merged[proposal.first],merged[proposal.second]);
        removed[proposal.second]=true;
    }
    std::vector<std::vector<std::size_t>> result;
    result.reserve(sourceLineages.size()-selected.size());
    for (std::size_t source=0;source<merged.size();++source) {
        if (!removed[source]) result.push_back(std::move(merged[source]));
    }
    return result;
}

[[nodiscard]] std::vector<bool> applySourceMergePreservation(
    const std::vector<bool>& preserveSources,
    const std::vector<SourceMergeProposal2D>& selected) {
    if (preserveSources.empty()) return {};
    std::vector<bool> removed(preserveSources.size(),false);
    std::vector<bool> resultFlags=preserveSources;
    for (const auto& proposal:selected) {
        if (proposal.first<resultFlags.size()) resultFlags[proposal.first]=false;
        if (proposal.second<removed.size()) removed[proposal.second]=true;
    }
    std::vector<bool> result;
    result.reserve(preserveSources.size()-selected.size());
    for (std::size_t source=0;source<preserveSources.size();++source) {
        if (!removed[source]) result.push_back(resultFlags[source]);
    }
    return result;
}

[[nodiscard]] RepartitionBatch2D applyRepartitionBatch(
    const TopologyMesh2D& topology,
    const std::vector<RepartitionProposal2D>& selected,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol,SolverTopologyProfile2D* profile,
    const std::vector<bool>& immutableCells) {
    std::vector<const RepartitionProposal2D*> replacements(topology.cells.size(),nullptr);
    std::vector<bool> removed(topology.cells.size(),false);
    for (const auto& proposal:selected) {
        if (proposal.first>=topology.cells.size() || proposal.second>=topology.cells.size()) continue;
        replacements[proposal.first]=&proposal;
        removed[proposal.second]=true;
    }
    std::vector<CutCell2D> cells;
    std::vector<bool> rebuiltImmutable;
    cells.reserve(topology.cells.size());
    rebuiltImmutable.reserve(topology.cells.size());
    for (std::size_t cell=0;cell<topology.cells.size();++cell) {
        if (removed[cell]) continue;
        if (replacements[cell]) {
            const auto lineage=mergedLineage(
                topology.cells[replacements[cell]->first].sourceLineage,
                topology.cells[replacements[cell]->second].sourceLineage);
            cells.push_back(makeCell(cells.size(),replacements[cell]->firstPiece,tol,lineage));
            cells.push_back(makeCell(cells.size(),replacements[cell]->secondPiece,tol,lineage));
            rebuiltImmutable.push_back(false);
            rebuiltImmutable.push_back(false);
        } else {
            cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol,
                                     topology.cells[cell].sourceLineage));
            rebuiltImmutable.push_back(!immutableCells.empty() && immutableCells[cell]);
        }
    }
    const auto topologyStart=ProfileClock::now();
    RepartitionBatch2D rebuilt;
    rebuilt.topology=buildGlobalTopology(cells,domain,boundary,tol,topology.constructionRegistry);
    rebuilt.immutableCells=std::move(rebuiltImmutable);
    if (profile) {
        const double elapsed=profileSeconds(topologyStart);
        profile->buildGlobalTopologySeconds+=elapsed;
        profile->candidateGlobalRebuildSeconds+=elapsed;
        ++profile->globalTopologyRebuildCalls;
        ++profile->candidateTopologyCount;
    }
    return rebuilt;
}

[[nodiscard]] SolverQualityReport2D timedFullQuality(
    const TopologyMesh2D& topology,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool candidate) {
    const auto start=ProfileClock::now();
    auto quality=evaluateSolverQuality2D(topology,{},tol);
    if (profile) {
        const double elapsed=profileSeconds(start);
        profile->fullQualitySeconds+=elapsed;
        ++profile->fullQualityCalls;
        if (candidate) {
            profile->candidateQualitySeconds+=elapsed;
            ++profile->candidateQualityEvaluationCount;
        } else {
            profile->maximumQualityIssueCount=std::max(
                profile->maximumQualityIssueCount,quality.issues.size());
        }
    }
    return quality;
}

SolverLocalRepartitionResult2D repartitionSolverTopologyByQualityImpl(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool useBatch,
    const std::vector<bool>& initialImmutableCells={}) {
    SolverLocalRepartitionResult2D result;
    result.topology=topology;
    result.immutableCells=initialImmutableCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back("local solver repartition requires valid inputs");
        return result;
    }
    for (std::size_t iteration=0;iteration<128;++iteration) {
        const auto quality=timedFullQuality(result.topology,tol,profile,false);
        if (quality.valid()) break;
        if (profile) ++profile->repartitionIterations;
        const auto generationStart=ProfileClock::now();
        const auto pairs=solverRepartitionPairs(
            result.topology,quality,result.immutableCells);
        if (profile) {
            profile->candidateGenerationSeconds+=profileSeconds(generationStart);
            profile->repartitionCandidatePairs+=pairs.size();
        }
        bool acceptedBatch=false;
        if (useBatch) {
            std::vector<RepartitionProposal2D> proposals;
            for (const auto& [first,second]:pairs) {
                const auto polygonStart=ProfileClock::now();
                const auto merged=mergeAdjacentPolygonsSimple(
                    topologyCellPolygon(result.topology,first),
                    topologyCellPolygon(result.topology,second),tol);
                if (!merged) {
                    if (profile) {
                        profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                    }
                    continue;
                }
                const auto splits=convexTwoPieceSplits(*merged,domain,boundary,tol);
                if (profile) profile->candidateSplits+=splits.size();
                std::optional<RepartitionProposal2D> best;
                const auto halo=cellPairHalo(result.topology,first,second);
                for (const auto& [firstPiece,secondPiece]:splits) {
                    const auto rank=localRepartitionQualityRank(
                        result.topology,halo,first,second,firstPiece,secondPiece,tol);
                    if (!rank) continue;
                    RepartitionProposal2D proposal;
                    proposal.first=first;
                    proposal.second=second;
                    proposal.firstPiece=firstPiece;
                    proposal.secondPiece=secondPiece;
                    proposal.rank=*rank;
                    proposal.halo=halo;
                    if (!best || betterLocalQualityRank(proposal.rank,best->rank)) {
                        best=std::move(proposal);
                    }
                }
                if (best) proposals.push_back(std::move(*best));
                if (profile) {
                    profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                }
            }
            auto selected=selectIndependentProposals(std::move(proposals));
            if (profile) profile->repairPatchCount+=selected.size();
            for (std::size_t batchSize=selected.size();batchSize>0;) {
                const std::vector<RepartitionProposal2D> batch(
                    selected.begin(),selected.begin()+static_cast<std::ptrdiff_t>(batchSize));
                auto candidate=applyRepartitionBatch(
                    result.topology,batch,domain,boundary,tol,profile,
                    result.immutableCells);
                if (candidate.topology.valid()) {
                    const auto candidateQuality=timedFullQuality(
                        candidate.topology,tol,profile,true);
                    if (betterQualityScore(qualityScore(candidateQuality),qualityScore(quality))) {
                        result.topology=std::move(candidate.topology);
                        result.immutableCells=std::move(candidate.immutableCells);
                        result.repartitionCount+=batchSize;
                        if (profile) {
                            profile->acceptedRepartitions+=batchSize;
                            ++profile->acceptedTopologyCommitCount;
                        }
                        acceptedBatch=true;
                        break;
                    }
                }
                if (batchSize==1) break;
                batchSize=(batchSize+1)/2;
            }
        }
        if (acceptedBatch) continue;

        // Fail closed to the H2 exhaustive global check if the locally ranked
        // independent batch does not strictly improve the authoritative score.
        std::optional<RepartitionBatch2D> bestTopology;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:pairs) {
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                topologyCellPolygon(result.topology,first),
                topologyCellPolygon(result.topology,second),tol);
            if (!merged) {
                if (profile) profile->candidatePolygonWorkSeconds+=
                    profileSeconds(polygonStart);
                continue;
            }
            // A deterministic convex partition may create a tiny triangle
            // beside a regular cell. If their exact union is already convex,
            // retaining two cells cannot improve interpolation weight or
            // volume ratio. Allow a true solver-cell agglomeration, subject to
            // the same full-topology validity and strict quality-score gate.
            if (strictlyConvex(*merged) &&
                !underDeterminedBoundaryCell(*merged,domain,boundary,tol)) {
                auto candidate=agglomerateCellPair(
                    result.topology,first,second,*merged,domain,boundary,tol,
                    profile,result.immutableCells);
                if (candidate.topology.valid()) {
                    const auto candidateQuality=timedFullQuality(
                        candidate.topology,tol,profile,true);
                    const auto candidateScore=qualityScore(candidateQuality);
                    if (betterQualityScore(candidateScore,bestScore)) {
                        bestScore=candidateScore;
                        bestTopology=std::move(candidate);
                    }
                }
            }
            const auto splits=convexTwoPieceSplits(*merged,domain,boundary,tol);
            if (profile) {
                profile->candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
            }
            for (const auto& [firstPiece,secondPiece]:splits) {
                auto candidate=repartitionPair(result.topology,first,second,
                                               firstPiece,secondPiece,
                                               domain,boundary,tol,profile,
                                               result.immutableCells);
                if (!candidate.topology.valid()) continue;
                const auto candidateQuality=timedFullQuality(
                    candidate.topology,tol,profile,true);
                const auto candidateScore=qualityScore(candidateQuality);
                if (betterQualityScore(candidateScore,bestScore)) {
                    bestScore=candidateScore;
                    bestTopology=std::move(candidate);
                }
            }
        }
        if (!bestTopology) break;
        result.topology=std::move(bestTopology->topology);
        result.immutableCells=std::move(bestTopology->immutableCells);
        ++result.repartitionCount;
        if (profile) ++profile->acceptedRepartitions;
        if (profile) ++profile->acceptedTopologyCommitCount;
    }
    return result;
}

} // namespace

std::vector<std::size_t> selectIndependentSolverRepairPatches2D(
    const std::vector<std::vector<std::size_t>>& halos) {
    std::vector<std::size_t> selected;
    for (std::size_t candidate=0;candidate<halos.size();++candidate) {
        bool conflict=false;
        for (const auto accepted:selected) {
            if (sortedSetsIntersect(halos[candidate],halos[accepted])) {
                conflict=true;
                break;
            }
        }
        if (!conflict) selected.push_back(candidate);
    }
    return selected;
}

SolverQualityReport2D evaluateSolverQualityPatch2D(
    const TopologyMesh2D& topology,
    const std::vector<std::size_t>& sortedCellIds,
    const TolerancePolicy& tol) {
    if (!topology.valid() || sortedCellIds.empty() ||
        !std::is_sorted(sortedCellIds.begin(),sortedCellIds.end()) ||
        std::adjacent_find(sortedCellIds.begin(),sortedCellIds.end())!=sortedCellIds.end() ||
        sortedCellIds.back()>=topology.cells.size()) {
        return evaluateSolverQuality2D(TopologyMesh2D{}, {}, tol);
    }
    const auto boundary=patchBoundaryRegion(topology,sortedCellIds,tol);
    if (!boundary) return evaluateSolverQuality2D(TopologyMesh2D{}, {}, tol);
    const Domain2D domain{boundary->bounds()};
    std::vector<CutCell2D> cells;
    cells.reserve(sortedCellIds.size());
    for (const auto cell:sortedCellIds) {
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol));
    }
    const auto patch=buildGlobalTopology(cells,domain,*boundary,tol);
    return evaluateSolverQuality2D(patch,{},tol);
}

SolverLocalRepartitionResult2D repartitionSolverTopologyByQuality2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return repartitionSolverTopologyByQualityImpl(
        topology,domain,boundary,tol,nullptr,true);
}

SolverLocalRepartitionResult2D
repartitionSolverTopologyByQualitySequentialReference2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return repartitionSolverTopologyByQualityImpl(
        topology,domain,boundary,tol,nullptr,false);
}

SolverShortFaceRepairResult2D repairSolverShortFaces2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const std::vector<bool>& immutableCells,
    const std::vector<double>& localBackgroundH,
    const std::vector<bool>& ratedCells,double minimumFaceOverLocalH,
    const TolerancePolicy& tol) {
    const auto repairStart=ProfileClock::now();
    SolverShortFaceRepairResult2D result;
    result.topology=topology;
    result.immutableCells=immutableCells;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        immutableCells.size()!=topology.cells.size() ||
        localBackgroundH.size()!=topology.cells.size() ||
        ratedCells.size()!=topology.cells.size() || !(minimumFaceOverLocalH>0.0)) {
        result.issues.push_back("revisioned short-face repair requires valid aligned inputs");
        return result;
    }
    struct FaceScore {
        std::size_t hardCount=0;
        double minimum=std::numeric_limits<double>::infinity();
        double maximumSeverity=0.0;
        double totalSeverity=0.0;
    };
    const auto scoreFaces=[&](const TopologyMesh2D& mesh,
                              const std::vector<double>& localH,
                              const std::vector<bool>& rated) {
        FaceScore score;
        for (const auto& edge:mesh.edges) {
            double h=0.0;
            if (rated.at(edge.owner)) h=std::max(h,localH.at(edge.owner));
            if (edge.neighbour && rated.at(*edge.neighbour))
                h=std::max(h,localH.at(*edge.neighbour));
            if (!(h>0.0)) continue;
            const double ratio=std::sqrt(squaredNorm(
                mesh.vertices[edge.v1].point-mesh.vertices[edge.v0].point))/h;
            score.minimum=std::min(score.minimum,ratio);
            if (ratio<minimumFaceOverLocalH) {
                ++score.hardCount;
                const double severity=minimumFaceOverLocalH/
                    (ratio+std::numeric_limits<double>::min());
                score.maximumSeverity=std::max(score.maximumSeverity,severity);
                score.totalSeverity+=severity;
            }
        }
        if (!std::isfinite(score.minimum)) score.minimum=0.0;
        return score;
    };
    const auto better=[](const FaceScore& lhs,const FaceScore& rhs) {
        return std::tie(lhs.hardCount,lhs.maximumSeverity,lhs.totalSeverity) <
               std::tie(rhs.hardCount,rhs.maximumSeverity,rhs.totalSeverity);
    };
    const auto legacyNoWorse=[](const SolverQualityReport2D& candidate,
                                const SolverQualityReport2D& current) {
        if (!candidate.valid()) return false;
        const auto upper=[](double after,double before) {
            return after-before<=1.0e-8*std::max(1.0,std::abs(before));
        };
        const auto lower=[](double after,double before) {
            return before-after<=1.0e-8*std::max(1.0,std::abs(before));
        };
        return upper(candidate.maxNonOrthogonalityDeg,current.maxNonOrthogonalityDeg) &&
               upper(candidate.maxInternalSkewness,current.maxInternalSkewness) &&
               upper(candidate.maxBoundarySkewness,current.maxBoundarySkewness) &&
               upper(candidate.maxConcavityDeg,current.maxConcavityDeg) &&
               upper(candidate.maxCellAspect,current.maxCellAspect) &&
               lower(candidate.minInteriorAngleDeg,current.minInteriorAngleDeg) &&
               lower(candidate.minFaceLength,current.minFaceLength) &&
               lower(candidate.minFaceWeight,current.minFaceWeight) &&
               lower(candidate.minVolumeRatio,current.minVolumeRatio) &&
               lower(candidate.minCompactness,current.minCompactness);
    };
    const auto currentQuality=evaluateSolverQuality2D(topology,{},tol);
    ++result.authoritativeFullQualityEvaluationCount;
    if (!currentQuality.valid()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const auto currentScore=scoreFaces(topology,localBackgroundH,ratedCells);
    result.hardFaceCountBefore=currentScore.hardCount;
    result.hardFaceCountAfter=currentScore.hardCount;
    result.minimumFaceOverLocalHBefore=currentScore.minimum;
    result.minimumFaceOverLocalHAfter=currentScore.minimum;
    if (currentScore.hardCount==0U) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }

    std::vector<std::tuple<double,std::size_t,std::size_t>> shortEdges;
    for (const auto& edge:topology.edges) {
        if (!edge.neighbour) continue;
        double h=0.0;
        if (ratedCells[edge.owner]) h=std::max(h,localBackgroundH[edge.owner]);
        if (ratedCells[*edge.neighbour]) h=std::max(h,localBackgroundH[*edge.neighbour]);
        if (!(h>0.0)) continue;
        const double ratio=std::sqrt(squaredNorm(
            topology.vertices[edge.v1].point-topology.vertices[edge.v0].point))/h;
        if (ratio<minimumFaceOverLocalH)
            shortEdges.emplace_back(ratio,std::min(edge.owner,*edge.neighbour),
                                    std::max(edge.owner,*edge.neighbour));
    }
    std::sort(shortEdges.begin(),shortEdges.end());
    if (shortEdges.empty()) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const std::vector<std::size_t> affected{std::get<1>(shortEdges.front()),
                                            std::get<2>(shortEdges.front())};
    // This R1 checkpoint migrates the narrow-gap Q2-B mode only: the failing
    // atomic face is a common-partition fragment on an immutable layer
    // support. Mutable/mutable sharp-tail patches still belong to the later
    // migration stage and must retain their pre-R1 behavior.
    if (!immutableCells[affected[0]] && !immutableCells[affected[1]]) {
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    result.applicable=true;
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    if (!immutableCells[affected[0]] && !immutableCells[affected[1]])
        pairs.emplace_back(affected[0],affected[1]);
    for (const auto cell:affected) {
        if (immutableCells[cell]) continue;
        for (const auto edgeId:topology.cells[cell].edges) {
            const auto& edge=topology.edges[edgeId];
            if (!edge.neighbour) continue;
            const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
            if (!immutableCells[other])
                pairs.emplace_back(std::min(cell,other),std::max(cell,other));
        }
    }
    std::sort(pairs.begin(),pairs.end());
    pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());

    const auto incidence=buildEdgeIncidenceStore2D(topology,0U);
    if (!incidence.valid()) {
        result.issues.push_back("short-face base incidence is invalid");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    std::map<std::pair<double,double>,TopologyDeltaVertex2D> identityByPoint;
    std::map<StableVertexId2D,TopologyDeltaVertex2D> identityById;
    for (std::size_t dense=0;dense<topology.vertices.size();++dense) {
        const auto id=incidence.stableVertexIds[dense];
        const TopologyDeltaVertex2D identity{id,
            {StableVertexKeyKind2D::LegacyCanonical,id,0U,0U,0U},
            topology.vertices[dense].point,true};
        identityByPoint.emplace(
            std::pair{topology.vertices[dense].point.x,topology.vertices[dense].point.y},
            identity);
        identityById.emplace(id,identity);
    }
    const auto sourceCell=[&](std::size_t cellId,const Polygon2D& polygon) {
        CutCell2D cell=makeCell(topology.cells[cellId].sourceId,polygon,tol,
                                topology.cells[cellId].sourceLineage);
        cell.sourceKey=topology.cells[cellId].sourceKey;
        cell.refinementLineageKeys=topology.cells[cellId].refinementLineageKeys;
        if (topology.constructionRegistry) {
            for (const auto& point:polygon.vertices)
                cell.canonicalVertexIds.push_back(
                    identityByPoint.at({point.x,point.y}).stableId);
            cell.canonicalRegistry=topology.constructionRegistry.get();
        }
        return cell;
    };
    std::vector<CutCell2D> baseSources;
    baseSources.reserve(topology.cells.size());
    for (const auto& cell:topology.cells)
        baseSources.push_back(sourceCell(cell.id,topologyCellPolygon(topology,cell.id)));

    // Phase 1: every candidate is generated, validated and ranked using only
    // the affected patch and its one-ring halo. No global topology build and no
    // full global quality evaluation may happen in this loop.
    struct ShortFaceCandidate2D {
        std::size_t first = 0;
        std::size_t second = 0;
        std::vector<std::size_t> selected;
        TopologyPatchTransaction2D transaction;
        std::vector<TopologyReplacementCell2D> replacements;
        PatchLocalQuality2D baseLocal;
        PatchLocalQuality2D candidateLocal;
    };
    std::optional<ShortFaceCandidate2D> winner;
    std::string localRejectReason;
    // Measured, not asserted: any accidental global rebuild or full-quality
    // evaluation inside the selection loop shows up in these deltas.
    const auto selectionStartGlobalBuilds=globalTopologyBuildCount2D();
    const auto selectionStartFullQuality=solverQualityEvaluationCount2D();
    for (const auto& [first,second]:pairs) {
        ++result.candidateCount;
        const auto merged=mergeAdjacentPolygonsSimple(
            topologyCellPolygon(topology,first),topologyCellPolygon(topology,second),tol);
        if (!merged) continue;
        const auto repairUnion=removeArtificialCollinearVertices(*merged,tol);
        if (!strictlyConvex(repairUnion) ||
            underDeterminedBoundaryCell(repairUnion,domain,boundary,tol)) continue;
        std::vector<std::size_t> selected{first,second};
        std::map<std::size_t,Polygon2D> simplifiedImmutable;
        std::set<std::size_t> localImmutable;
        for (const auto cell:{first,second,affected[0],affected[1]}) {
            if (immutableCells[cell]) localImmutable.insert(cell);
            for (const auto edgeId:topology.cells[cell].edges) {
                const auto& edge=topology.edges[edgeId];
                if (!edge.neighbour) continue;
                const auto other=edge.owner==cell?*edge.neighbour:edge.owner;
                if (immutableCells[other]) localImmutable.insert(other);
            }
        }
        for (const auto cell:localImmutable) {
            const auto original=topologyCellPolygon(topology,cell);
            const auto simplified=removeArtificialCollinearVertices(original,tol);
            if (simplified.vertices.size()!=original.vertices.size()) {
                selected.push_back(cell);
                simplifiedImmutable.emplace(cell,simplified);
            }
        }
        std::sort(selected.begin(),selected.end());
        selected.erase(std::unique(selected.begin(),selected.end()),selected.end());
        const auto transaction=prepareTopologyPatchTransaction2D(topology,incidence,selected);
        if (!transaction.valid()) continue;
        const auto identified=[&](std::size_t sourceId,Polygon2D polygon) {
            // Insert locked stable endpoints on any simplified outer patch
            // edge. Identity comes from the lock key; geometry is only used to
            // order already-identified endpoints along that local edge.
            Polygon2D locked;
            for (std::size_t i=0;i<polygon.vertices.size();++i) {
                const auto a=polygon.vertices[i];
                const auto b=polygon.vertices[(i+1U)%polygon.vertices.size()];
                std::vector<std::pair<double,TopologyDeltaVertex2D>> points;
                const auto add=[&](StableVertexId2D id) {
                    const auto found=identityById.find(id);
                    if (found==identityById.end() ||
                        !pointOnSegment(found->second.point,{a,b},tol)) return;
                    const auto direction=b-a;
                    const double denom=squaredNorm(direction);
                    const double t=denom>0.0?dot(found->second.point-a,direction)/denom:0.0;
                    points.emplace_back(t,found->second);
                };
                add(identityByPoint.at({a.x,a.y}).stableId);
                for (const auto& lock:transaction.boundaryLocks) {
                    add(lock.stableEdge.v0);
                    add(lock.stableEdge.v1);
                }
                std::sort(points.begin(),points.end(),[](const auto& lhs,const auto& rhs) {
                    return std::tie(lhs.first,lhs.second.stableId)<
                           std::tie(rhs.first,rhs.second.stableId);
                });
                points.erase(std::unique(points.begin(),points.end(),[](const auto& lhs,
                                                                        const auto& rhs) {
                    return lhs.second.stableId==rhs.second.stableId;
                }),points.end());
                for (const auto& [t,vertex]:points) {
                    (void)t;
                    if (locked.vertices.empty() || locked.vertices.back().x!=vertex.point.x ||
                        locked.vertices.back().y!=vertex.point.y)
                        locked.vertices.push_back(vertex.point);
                }
            }
            if (locked.vertices.size()>1U &&
                locked.vertices.front().x==locked.vertices.back().x &&
                locked.vertices.front().y==locked.vertices.back().y)
                locked.vertices.pop_back();
            TopologyReplacementCell2D replacement;
            replacement.cell=sourceCell(sourceId,locked);
            for (const auto& point:locked.vertices)
                replacement.vertices.push_back(identityByPoint.at({point.x,point.y}));
            return replacement;
        };
        std::vector<TopologyReplacementCell2D> replacements;
        std::vector<std::pair<double,bool>> replacementMetadata;
        replacements.push_back(identified(first,repairUnion));
        replacements.back().cell.sourceLineage=mergedLineage(
            topology.cells[first].sourceLineage,topology.cells[second].sourceLineage);
        replacementMetadata.emplace_back(
            std::max(localBackgroundH[first],localBackgroundH[second]),
            ratedCells[first] || ratedCells[second]);
        for (const auto& [cell,polygon]:simplifiedImmutable) {
            replacements.push_back(identified(cell,polygon));
            replacementMetadata.emplace_back(localBackgroundH[cell],ratedCells[cell]);
        }
        const auto localDelta=buildTopologyDelta2D(
            topology,incidence,transaction,replacements,tol);
        if (!localDelta.valid()) {
            if (!localDelta.issues.empty()) localRejectReason=localDelta.issues.front();
            continue;
        }
        ++result.localCandidateCount;

        const auto baseScope=buildPatchLocalScope2D(
            topology,incidence,selected,localBackgroundH,ratedCells);
        if (!baseScope.valid()) {
            localRejectReason=baseScope.issues.front();
            continue;
        }
        std::set<StableEdgeKey2D> physicalLocks;
        for (const auto& lock:transaction.boundaryLocks)
            if (lock.physicalBoundary) physicalLocks.insert(lock.stableEdge);
        std::vector<PatchLocalCell2D> candidateScope;
        for (std::size_t i=0;i<baseScope.cells.size();++i)
            if (!baseScope.cells[i].inPatch) candidateScope.push_back(baseScope.cells[i]);
        for (std::size_t i=0;i<replacements.size();++i) {
            const auto& replacement=replacements[i];
            PatchLocalCell2D entry;
            entry.inPatch=true;
            entry.localBackgroundH=replacementMetadata[i].first;
            entry.ratedForFaceLength=replacementMetadata[i].second;
            entry.polygon=replacement.cell.fluidPolygon;
            for (std::size_t local=0;local<replacement.vertices.size();++local) {
                const auto from=replacement.vertices[local].stableId;
                const auto to=replacement.vertices[
                    (local+1U)%replacement.vertices.size()].stableId;
                const auto endpoints=std::minmax(from,to);
                entry.loop.push_back(from);
                entry.physicalBoundaryFace.push_back(
                    physicalLocks.contains(StableEdgeKey2D{endpoints.first,endpoints.second}));
            }
            candidateScope.push_back(std::move(entry));
        }
        const auto baseLocal=evaluatePatchLocalQuality2D(
            baseScope.cells,minimumFaceOverLocalH,{},tol);
        const auto candidateLocal=evaluatePatchLocalQuality2D(
            candidateScope,minimumFaceOverLocalH,{},tol);
        result.localQualityEvaluationCount+=2U;
        if (!baseLocal.valid() || !candidateLocal.valid()) {
            localRejectReason=baseLocal.valid()?candidateLocal.issues.front()
                                               :baseLocal.issues.front();
            continue;
        }
        // Every excluded cell and face is identical in both scopes, so a local
        // improvement on a monotone aggregate is also a global improvement.
        if (!patchLocalQualityNoWorse2D(candidateLocal,baseLocal)) {
            localRejectReason="patch-local candidate quality is worse than the patch it replaces";
            continue;
        }
        if (!(std::tie(candidateLocal.hardShortFaceCount,
                       candidateLocal.maximumShortFaceSeverity,
                       candidateLocal.totalShortFaceSeverity)<
              std::tie(baseLocal.hardShortFaceCount,
                       baseLocal.maximumShortFaceSeverity,
                       baseLocal.totalShortFaceSeverity))) {
            localRejectReason="patch-local candidate does not reduce the local short-face score";
            continue;
        }
        ShortFaceCandidate2D accepted{first,second,selected,transaction,
                                      replacements,baseLocal,candidateLocal};
        if (!winner || patchLocalQualityRanksBetter2D(candidateLocal,
                                                     winner->candidateLocal))
            winner=std::move(accepted);
    }
    if (!winner) {
        result.candidateGlobalTopologyBuildCount=
            globalTopologyBuildCount2D()-selectionStartGlobalBuilds;
        result.candidateFullGlobalQualityEvaluationCount=
            solverQualityEvaluationCount2D()-selectionStartFullQuality;
        result.issues.push_back(localRejectReason.empty()
            ?"no patch-local short-face candidate passed the local quality gate"
            :localRejectReason);
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    result.candidateGlobalTopologyBuildCount=
        globalTopologyBuildCount2D()-selectionStartGlobalBuilds;
    result.candidateFullGlobalQualityEvaluationCount=
        solverQualityEvaluationCount2D()-selectionStartFullQuality;

    // Phase 2: exactly one global oracle/materialization and one authoritative
    // full-quality evaluation, for the selected winner only. A disagreement
    // fails the transaction closed rather than retrying other candidates.
    const auto& first=winner->first;
    const auto& second=winner->second;
    const auto& selected=winner->selected;
    const auto committed=evaluateTopologyPatchTransactionOracle2D(
        topology,incidence,winner->transaction,baseSources,winner->replacements,
        domain,boundary,{true,0.0},tol);
    result.globalOracleBuildCount+=committed.globalOracleBuildCount;
    if (!committed.accepted || !committed.deltaMatchesGlobalOracle) {
        result.issues.push_back(committed.issues.empty()
            ?"winner patch transaction failed the global oracle"
            :committed.issues.front());
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    std::map<TopologySourceIdentity2D,std::pair<double,bool>> metadata;
    std::set<TopologySourceIdentity2D> immutableIdentities;
    for (const auto& cell:topology.cells) {
        const TopologySourceIdentity2D identity{cell.sourceKey,cell.sourceId};
        metadata[identity]={localBackgroundH[cell.id],ratedCells[cell.id]};
        if (immutableCells[cell.id]) immutableIdentities.insert(identity);
    }
    metadata[{topology.cells[first].sourceKey,topology.cells[first].sourceId}]=
        {std::max(localBackgroundH[first],localBackgroundH[second]),
         ratedCells[first] || ratedCells[second]};
    std::vector<double> candidateH;
    std::vector<bool> candidateRated,candidateImmutable;
    for (const auto& cell:committed.topology.cells) {
        const auto identity=TopologySourceIdentity2D{cell.sourceKey,cell.sourceId};
        const auto item=metadata.at(identity);
        candidateH.push_back(item.first);
        candidateRated.push_back(item.second);
        candidateImmutable.push_back(immutableIdentities.contains(identity));
    }
    const auto candidateQuality=evaluateSolverQuality2D(committed.topology,{},tol);
    ++result.authoritativeFullQualityEvaluationCount;
    const auto candidateScore=scoreFaces(committed.topology,candidateH,candidateRated);
    result.localWinnerMatchesGlobalAuthority=
        legacyNoWorse(candidateQuality,currentQuality) &&
        better(candidateScore,currentScore);
    if (!result.localWinnerMatchesGlobalAuthority) {
        result.issues.push_back(
            "patch-local winner disagrees with the authoritative global quality result");
        result.repairSeconds=profileSeconds(repairStart);
        return result;
    }
    const std::set<std::size_t> selectedSet(selected.begin(),selected.end());
    bool outsideStable=true;
    const auto revisionedBase=buildRevisionedTopology2D(topology,incidence);
    for (const auto& cell:topology.cells) if (!selectedSet.contains(cell.id)) {
        const TopologySourceIdentity2D identity{cell.sourceKey,cell.sourceId};
        outsideStable=outsideStable &&
            committed.revisionedTopology.cells.at(identity).vertices==
            revisionedBase.cells.at(identity).vertices;
    }
    result.topology=committed.topology;
    result.immutableCells=std::move(candidateImmutable);
    result.accepted=true;
    result.hardFaceCountAfter=candidateScore.hardCount;
    result.minimumFaceOverLocalHAfter=candidateScore.minimum;
    result.patchOutsideStableIdsUnchanged=outsideStable;
    result.localDeltaMatchesGlobalOracle=committed.deltaMatchesGlobalOracle;
    result.repairSeconds=profileSeconds(repairStart);
    return result;
}

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    return buildSolverTopology2D(
        topology,domain,boundary,SolverTopologyConstraints2D{},tol);
}

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const SolverTopologyConstraints2D& constraints,
    const TolerancePolicy& tol) {
    SolverTopologyResult2D result;
    result.inputCellCount=topology.cells.size();
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        (!constraints.immutableInputCells.empty() &&
         constraints.immutableInputCells.size()!=topology.cells.size()) ||
        (!constraints.preserveInputCells.empty() &&
         constraints.preserveInputCells.size()!=topology.cells.size()) ||
        (!constraints.inputPolygonOverrides.empty() &&
         constraints.inputPolygonOverrides.size()!=topology.cells.size())) {
        result.issues.push_back("solver topology repair requires valid inputs");
        return result;
    }
    std::vector<Polygon2D> sourcePolygons;
    std::vector<std::vector<std::size_t>> sourceLineages;
    sourcePolygons.reserve(topology.cells.size());
    sourceLineages.reserve(topology.cells.size());
    for (std::size_t cellId=0;cellId<topology.cells.size();++cellId) {
        const auto& cell=topology.cells[cellId];
        sourceLineages.push_back(cell.sourceLineage.empty()
            ?std::vector<std::size_t>{cell.sourceId}:cell.sourceLineage);
        if (!constraints.inputPolygonOverrides.empty() &&
            constraints.inputPolygonOverrides[cellId].has_value()) {
            const auto& polygon=*constraints.inputPolygonOverrides[cellId];
            const BoundaryLoop loop(polygon.vertices);
            const double areaScale=std::max(std::abs(polygon.area()),
                                            std::abs(cell.geometryArea));
            if (!loop.diagnose(tol).valid() || !(polygon.signedArea()>0.0) ||
                std::abs(polygon.area()-cell.geometryArea)>
                    tol.absolute*tol.absolute+tol.relative*areaScale) {
                result.issues.push_back(
                    "solver topology polygon override does not match input cell");
                return result;
            }
            sourcePolygons.push_back(polygon);
            continue;
        }
        Polygon2D polygon;
        polygon.vertices.reserve(cell.vertices.size());
        for (const auto vertex:cell.vertices) {
            if (vertex>=topology.vertices.size()) {
                result.issues.push_back("solver topology cell has invalid vertex index");
                return result;
            }
            polygon.vertices.push_back(topology.vertices[vertex].point);
        }
        sourcePolygons.push_back(std::move(polygon));
    }

    const auto initialPartitionStart=ProfileClock::now();
    PartitionAttempt2D partition=partitionSourcePolygons(
        sourcePolygons,domain,boundary,tol,&result.profile,false,
        constraints.immutableInputCells,constraints.preserveInputCells,
        sourceLineages,topology.constructionRegistry);
    std::vector<bool> immutableSources=constraints.immutableInputCells;
    std::vector<bool> preserveSources=constraints.preserveInputCells;
    result.profile.initialPartitionSeconds=profileSeconds(initialPartitionStart);
    if (!partition.error.empty()) {
        result.issues.push_back(partition.error);
        return result;
    }
    const auto sourceRepairStart=ProfileClock::now();
    for (std::size_t iteration=0;iteration<32;++iteration) {
        const auto quality=timedFullQuality(partition.topology,tol,&result.profile,false);
        if (quality.valid()) break;
        ++result.profile.sourceRepairIterations;
        const auto generationStart=ProfileClock::now();
        const auto pairs=sourceMergePairs(partition,quality,immutableSources);
        result.profile.candidateGenerationSeconds+=profileSeconds(generationStart);
        result.profile.sourceCandidatePairs+=pairs.size();
        const auto adjacency=sourceAdjacency(partition,sourcePolygons.size());
        std::vector<SourceMergeProposal2D> proposals;
        for (const auto& [first,second]:pairs) {
            if (second>=sourcePolygons.size()) continue;
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                sourcePolygons[first],sourcePolygons[second],tol);
            if (merged) {
                const auto pieces=partitionOnePolygon(*merged,domain,boundary,tol);
                if (pieces) {
                    const auto rank=localPartitionQualityRank(*merged,*pieces,tol);
                    if (rank) {
                        SourceMergeProposal2D proposal;
                        proposal.first=first;
                        proposal.second=second;
                        proposal.merged=*merged;
                        proposal.rank=*rank;
                        proposal.halo=sourcePairHalo(first,second,adjacency);
                        proposals.push_back(std::move(proposal));
                    }
                }
            }
            result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
        }
        auto selected=selectIndependentProposals(std::move(proposals));
        result.profile.repairPatchCount+=selected.size();
        bool acceptedBatch=false;
        for (std::size_t batchSize=selected.size();batchSize>0;) {
            const std::vector<SourceMergeProposal2D> batch(
                selected.begin(),selected.begin()+static_cast<std::ptrdiff_t>(batchSize));
            auto candidatePolygons=applySourceMergeBatch(sourcePolygons,batch);
            auto candidateImmutable=applySourceMergeProtection(immutableSources,batch);
            auto candidatePreserve=applySourceMergePreservation(preserveSources,batch);
            auto candidateLineages=applySourceMergeLineage(sourceLineages,batch);
            auto candidate=partitionSourcePolygons(
                candidatePolygons,domain,boundary,tol,&result.profile,true,
                candidateImmutable,candidatePreserve,candidateLineages,
                topology.constructionRegistry);
            if (candidate.error.empty()) {
                const auto candidateQuality=timedFullQuality(
                    candidate.topology,tol,&result.profile,true);
                if (betterQualityScore(qualityScore(candidateQuality),qualityScore(quality))) {
                    partition=std::move(candidate);
                    sourcePolygons=std::move(candidatePolygons);
                    immutableSources=std::move(candidateImmutable);
                    preserveSources=std::move(candidatePreserve);
                    sourceLineages=std::move(candidateLineages);
                    result.qualityAgglomeratedSourceCellCount+=batchSize;
                    result.profile.acceptedSourceRepairs+=batchSize;
                    ++result.profile.acceptedTopologyCommitCount;
                    acceptedBatch=true;
                    break;
                }
            }
            if (batchSize==1) break;
            batchSize=(batchSize+1)/2;
        }
        if (acceptedBatch) continue;

        // Preserve the H2 exhaustive path as a correctness fallback when the
        // independent batch is invalid or does not improve the global score.
        std::optional<PartitionAttempt2D> bestPartition;
        std::optional<std::vector<Polygon2D>> bestPolygons;
        std::optional<std::vector<bool>> bestImmutable;
        std::optional<std::vector<bool>> bestPreserve;
        std::optional<std::vector<std::vector<std::size_t>>> bestLineages;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:pairs) {
            if (second>=sourcePolygons.size()) continue;
            const auto polygonStart=ProfileClock::now();
            const auto merged=mergeAdjacentPolygonsSimple(
                sourcePolygons[first],sourcePolygons[second],tol);
            if (!merged) {
                result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
                continue;
            }
            auto candidatePolygons=sourcePolygons;
            auto candidateImmutable=immutableSources;
            auto candidatePreserve=preserveSources;
            auto candidateLineages=sourceLineages;
            candidatePolygons[first]=*merged;
            candidatePolygons.erase(
                candidatePolygons.begin()+static_cast<std::ptrdiff_t>(second));
            if (!candidateImmutable.empty()) {
                candidateImmutable.erase(
                    candidateImmutable.begin()+static_cast<std::ptrdiff_t>(second));
            }
            if (!candidatePreserve.empty()) {
                candidatePreserve[first]=false;
                candidatePreserve.erase(
                    candidatePreserve.begin()+static_cast<std::ptrdiff_t>(second));
            }
            candidateLineages[first]=mergedLineage(
                candidateLineages[first],candidateLineages[second]);
            candidateLineages.erase(
                candidateLineages.begin()+static_cast<std::ptrdiff_t>(second));
            result.profile.candidatePolygonWorkSeconds+=profileSeconds(polygonStart);
            auto candidate=partitionSourcePolygons(
                candidatePolygons,domain,boundary,tol,&result.profile,true,
                candidateImmutable,candidatePreserve,candidateLineages,
                topology.constructionRegistry);
            if (!candidate.error.empty()) continue;
            const auto candidateQuality=timedFullQuality(
                candidate.topology,tol,&result.profile,true);
            const auto candidateScore=qualityScore(candidateQuality);
            if (betterQualityScore(candidateScore,bestScore)) {
                bestScore=candidateScore;
                bestPartition=std::move(candidate);
                bestPolygons=std::move(candidatePolygons);
                bestImmutable=std::move(candidateImmutable);
                bestPreserve=std::move(candidatePreserve);
                bestLineages=std::move(candidateLineages);
            }
        }
        if (!bestPartition) break;
        partition=std::move(*bestPartition);
        sourcePolygons=std::move(*bestPolygons);
        immutableSources=std::move(*bestImmutable);
        preserveSources=std::move(*bestPreserve);
        sourceLineages=std::move(*bestLineages);
        ++result.qualityAgglomeratedSourceCellCount;
        ++result.profile.acceptedSourceRepairs;
        ++result.profile.acceptedTopologyCommitCount;
    }
    result.profile.sourceRepairSeconds=profileSeconds(sourceRepairStart);
    const auto repartitionStart=ProfileClock::now();
    auto repartitioned=repartitionSolverTopologyByQualityImpl(
        partition.topology,domain,boundary,tol,&result.profile,true,
        partition.immutableForCell);
    result.profile.finalRepartitionSeconds=profileSeconds(repartitionStart);
    if (!repartitioned.valid()) {
        result.issues.insert(result.issues.end(),repartitioned.issues.begin(),
                             repartitioned.issues.end());
        return result;
    }
    partition.topology=std::move(repartitioned.topology);
    result.qualityRepartitionCount=repartitioned.repartitionCount;
    result.topology=std::move(partition.topology);
    result.immutableOutputCells=std::move(repartitioned.immutableCells);
    result.partitionedCellCount=partition.partitionedSourceCellCount;
    result.outputCellCount=result.topology.cells.size();
    return result;
}

} // namespace cartmesh2d
