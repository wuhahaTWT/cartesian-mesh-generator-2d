#include "SolverTopologyDetail2D.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>

namespace cartmesh2d::solver_topology_detail {

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

// Cartesian transition cells may have extra collinear vertices where a finer
// neighbour attaches. Keep those vertices and their shared faces, without adding
// diagonals inside a geometrically rectangular cell. Exact axis alignment is
// deliberate: this exception must not admit nearly-flat physical cut corners.
[[nodiscard]] bool cartesianRectangle(const Polygon2D& polygon) noexcept {
    if (polygon.vertices.size()<4 || !(polygon.signedArea()>0.0)) return false;
    unsigned corners=0;
    for (std::size_t i=0;i<polygon.vertices.size();++i) {
        const auto& a=polygon.vertices[i];
        const auto& b=polygon.vertices[(i+1)%polygon.vertices.size()];
        const auto& c=polygon.vertices[(i+2)%polygon.vertices.size()];
        if ((a.x!=b.x && a.y!=b.y) || (a.x==b.x && a.y==b.y)) return false;
        const int turn=orientationSign(a,b,c);
        if (turn<0) return false;
        if (turn>0) ++corners;
    }
    return corners==4;
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


namespace {

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygonBoundary(
    const Polygon2D& first,const Polygon2D& second,bool requireSingleSharedEdge) {
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
    if (requireSingleSharedEdge &&
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
            [&](const DirectedEdge& edge) { return same(edge.a,next); });
        if (following==remaining.end()) return std::nullopt;
        next=following->b;
        remaining.erase(following);
    }
    if (!same(next,start)) return std::nullopt;
    return Polygon2D{std::move(ordered)};
}

} // namespace

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygons(
    const Polygon2D& first,const Polygon2D& second) {
    auto merged=mergeAdjacentPolygonBoundary(first,second,true);
    if (!merged) return std::nullopt;
    const double sourceArea=first.area()+second.area();
    const double areaTolerance=1.0e-12+1.0e-10*sourceArea;
    return strictlyConvex(*merged) && std::abs(merged->area()-sourceArea)<=areaTolerance
        ?std::optional<Polygon2D>(std::move(*merged)):std::nullopt;
}

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentPolygonsSimple(
    const Polygon2D& first,const Polygon2D& second,const TolerancePolicy& tol) {
    auto merged=mergeAdjacentPolygonBoundary(first,second,false);
    if (!merged) return std::nullopt;
    const double sourceArea=first.area()+second.area();
    const double areaTolerance=tol.absolute*tol.absolute+tol.relative*sourceArea;
    BoundaryLoop loop(merged->vertices);
    if (!loop.diagnose(tol).valid() || !(merged->signedArea()>0.0) ||
        std::abs(merged->area()-sourceArea)>areaTolerance) return std::nullopt;
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
    Polygon2D polygon,const TolerancePolicy& tol,bool constructionRoundoff) {
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
            double crossBudget=64.0*std::numeric_limits<double>::epsilon()*scale;
            if (constructionRoundoff) {
                // Atomic intersection points are rounded in world coordinates.
                // An edge-angle-only epsilon loses those points' construction
                // error when a short fragment lies far from the origin.
                const double magnitude=std::max({std::abs(previous.x),std::abs(previous.y),
                    std::abs(current.x),std::abs(current.y),std::abs(next.x),std::abs(next.y)});
                const double roundoff=tol.constructionRoundoffScale(magnitude)*
                    (std::sqrt(squaredNorm(incoming))+std::sqrt(squaredNorm(outgoing)));
                crossBudget=std::max(crossBudget,std::min(tol.relative*scale,roundoff));
            }
            if (!(scale>0.0) || dot(incoming,outgoing)<0.0 ||
                std::abs(cross(incoming,outgoing))>crossBudget) continue;
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
                                 std::vector<std::size_t> sourceLineage) {
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

[[nodiscard]] QualityScore2D qualityScore(const SolverQualityReport2D& quality) {
    QualityScore2D score;
    score.issueCount=quality.issues.size();
    for (const auto& issue:quality.issues) {
        // Admissibility is categorical, with a zero dimensional limit. A
        // reversed signed distance must not become a negative severity or an
        // infinity that corrupts the candidate ordering.
        if (issue.code==SolverQualityIssueCode2D::InvalidFaceGeometry) {
            score.maximumSeverity=std::max(score.maximumSeverity,1.0);
            score.totalSeverity+=1.0;
            continue;
        }
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
    // Rank equal-quality patches independently of platform roundoff. Binning
    // each key preserves strict weak ordering, unlike pairwise epsilon tests.
    // This affects proposal order only; full quality acceptance uses raw values.
    const auto key=[](const LocalQualityRank2D& rank) {
        const auto bin=[](double value) {
            return std::round(value/TolerancePolicy{}.relative);
        };
        return std::tuple(rank.issueDelta,bin(rank.maximumSeverityDelta),
            bin(rank.totalSeverityDelta),rank.issues.issueCount,
            bin(rank.issues.maximumSeverity),bin(rank.issues.totalSeverity),
            bin(rank.maxNonOrthogonality),bin(rank.maxInternalSkewness),
            bin(rank.maxBoundarySkewness),bin(rank.maxCellAspect),
            bin(rank.negativeMinInteriorAngle),bin(rank.negativeMinFaceWeight),
            bin(rank.negativeMinVolumeRatio));
    };
    return key(candidate)<key(current);
}

void rankRelativeToBase(LocalQualityRank2D& rank,const LocalQualityRank2D& base) {
    rank.issueDelta=static_cast<std::ptrdiff_t>(rank.issues.issueCount)-
        static_cast<std::ptrdiff_t>(base.issues.issueCount);
    const auto finiteDelta=[](double candidate,double before) {
        const double delta=candidate-before;
        return std::isfinite(delta)?delta:0.0;
    };
    rank.maximumSeverityDelta=finiteDelta(rank.issues.maximumSeverity,base.issues.maximumSeverity);
    rank.totalSeverityDelta=finiteDelta(rank.issues.totalSeverity,base.issues.totalSeverity);
}

[[nodiscard]] std::optional<std::vector<Polygon2D>> partitionOnePolygon(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    // Candidate scoring and committed topology must use the same partition.
    if (cartesianRectangle(polygon) || strictlyConvex(polygon))
        return std::vector<Polygon2D>{polygon};
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
            issue.code==SolverQualityIssueCode2D::InvalidFaceGeometry ||
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
           code==SolverQualityIssueCode2D::InvalidFaceGeometry ||
           code==SolverQualityIssueCode2D::ExcessiveSkewness ||
           code==SolverQualityIssueCode2D::LowFaceWeight ||
           code==SolverQualityIssueCode2D::LowVolumeRatio;
}

[[nodiscard]] PartitionAttempt2D partitionSourcePolygons(
    const std::vector<Polygon2D>& sourcePolygons,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool candidate,
    const std::vector<bool>& immutableSources,
    const std::vector<bool>& preserveSources,
    const std::vector<std::vector<std::size_t>>& sourceLineages,
    std::shared_ptr<IntersectionRegistry2D> registry,
    const std::vector<Segment2D>& knownEmbeddedBoundary) {
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
        if (immutable || preserve || cartesianRectangle(polygon) || strictlyConvex(polygon)) {
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
    result.topology=buildGlobalTopology(
        cells,domain,boundary,tol,registry,knownEmbeddedBoundary);
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
        if (!internalInterfaceIssue(issue.code) &&
            issue.code!=SolverQualityIssueCode2D::ExcessiveBoundarySkewness) continue;
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

[[nodiscard]] std::vector<std::pair<Polygon2D,Polygon2D>> convexTwoPieceSplits(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    bool wallNormalCuts,bool allowCollinear) {
    std::vector<std::pair<Polygon2D,Polygon2D>> splits;
    const auto convex=[&](const Polygon2D& piece) {
        return strictlyConvex(allowCollinear?removeArtificialCollinearVertices(piece,tol):piece);
    };
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
            if (!convex(first) || !convex(second) ||
                underDeterminedBoundaryCell(first,domain,boundary,tol) ||
                underDeterminedBoundaryCell(second,domain,boundary,tol)) continue;
            const double areaError=std::abs(first.area()+second.area()-polygon.area());
            if (areaError>tol.absolute*tol.absolute+tol.relative*polygon.area()) continue;
            splits.emplace_back(std::move(first),std::move(second));
        }
    }
    if (wallNormalCuts) {
        // Short physical wall facets need a nearby cell centre. Existing-vertex
        // diagonals cannot always provide one: also cut inward from a wall
        // endpoint to the opposite edge, without moving the physical boundary.
        for (std::size_t wall=0;wall<n;++wall) {
            const auto a=polygon.vertices[wall];
            const auto b=polygon.vertices[(wall+1)%n];
            if (!boundaryEdge(a,b,domain,boundary,tol)) continue;
            const Vector2D normal{-(b.y-a.y),b.x-a.x};
            for (const auto start:{wall,(wall+1)%n}) {
                const auto p=polygon.vertices[start];
                for (std::size_t edge=0;edge<n;++edge) {
                    if (edge==start || (edge+1)%n==start) continue;
                    const auto c=polygon.vertices[edge];
                    const auto d=polygon.vertices[(edge+1)%n];
                    const auto tangent=d-c;
                    const double det=cross(normal,tangent);
                    if (std::abs(det)<=tol.relative*
                        std::sqrt(squaredNorm(normal)*squaredNorm(tangent))) continue;
                    const double distance=cross(c-p,tangent)/det;
                    const double fraction=cross(c-p,normal)/det;
                    if (distance<=tol.relative || fraction<=tol.relative ||
                        fraction>=1.0-tol.relative) continue;
                    const Point2D q{c.x+fraction*tangent.x,c.y+fraction*tangent.y};
                    Polygon2D first,second;
                    for (std::size_t k=start;;k=(k+1)%n) {
                        first.vertices.push_back(polygon.vertices[k]);
                        if (k==edge) break;
                    }
                    first.vertices.push_back(q);
                    second.vertices.push_back(q);
                    for (std::size_t k=(edge+1)%n;;k=(k+1)%n) {
                        second.vertices.push_back(polygon.vertices[k]);
                        if (k==start) break;
                    }
                    if (!convex(first) || !convex(second) ||
                        underDeterminedBoundaryCell(first,domain,boundary,tol) ||
                        underDeterminedBoundaryCell(second,domain,boundary,tol)) continue;
                    if (std::abs(first.area()+second.area()-polygon.area())>
                        tol.absolute*tol.absolute+tol.relative*polygon.area()) continue;
                    splits.emplace_back(std::move(first),std::move(second));
                }
            }
        }
    }
    return splits;
}

[[nodiscard]] std::vector<std::array<Polygon2D,3>> convexThreePieceFanSplits(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    std::vector<std::array<Polygon2D,3>> splits;
    const std::size_t n=polygon.vertices.size();
    if (n<5U) return splits;
    for (std::size_t anchor=0;anchor<n;++anchor) {
        std::vector<Point2D> rotated;
        rotated.reserve(n);
        for (std::size_t i=0;i<n;++i)
            rotated.push_back(polygon.vertices[(anchor+i)%n]);
        for (std::size_t firstCut=2U;firstCut+1U<n;++firstCut) {
            for (std::size_t secondCut=firstCut+1U;secondCut+1U<n;++secondCut) {
                std::array<Polygon2D,3> pieces;
                for (std::size_t i=0;i<=firstCut;++i)
                    pieces[0].vertices.push_back(rotated[i]);
                pieces[1].vertices.push_back(rotated[0]);
                for (std::size_t i=firstCut;i<=secondCut;++i)
                    pieces[1].vertices.push_back(rotated[i]);
                pieces[2].vertices.push_back(rotated[0]);
                for (std::size_t i=secondCut;i<n;++i)
                    pieces[2].vertices.push_back(rotated[i]);
                bool valid=true;
                double area=0.0;
                for (const auto& piece:pieces) {
                    valid=valid && strictlyConvex(piece) &&
                        !underDeterminedBoundaryCell(piece,domain,boundary,tol);
                    area+=piece.area();
                }
                const double areaError=std::abs(area-polygon.area());
                if (valid && areaError<=tol.absolute*tol.absolute+
                    tol.relative*polygon.area())
                    splits.push_back(std::move(pieces));
            }
        }
    }
    return splits;
}

[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>>
solverRepartitionPairs(const TopologyMesh2D& topology,
                       const SolverQualityReport2D& quality,
                       const std::vector<bool>& immutableCells) {
    std::vector<std::size_t> affected;
    for (const auto& issue:quality.issues) {
        if (!internalInterfaceIssue(issue.code) &&
            issue.code!=SolverQualityIssueCode2D::ExcessiveBoundarySkewness) continue;
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
    rebuilt.topology=buildGlobalTopology(
        cells,domain,boundary,tol,topology.constructionRegistry,
        topology.embeddedBoundaryFragments);
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
    rebuilt.topology=buildGlobalTopology(
        cells,domain,boundary,tol,topology.constructionRegistry,
        topology.embeddedBoundaryFragments);
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
    // Every boundary of this patch belongs to at least one selected cell.
    // Visit only that incidence, in the same edge-ID order as the full scan;
    // unrelated far-field edges cannot affect the patch boundary or ranking.
    std::vector<std::size_t> patchEdges;
    for (const auto cell:cells) {
        if (cell>=topology.cells.size()) return std::nullopt;
        const auto& edges=topology.cells[cell].edges;
        patchEdges.insert(patchEdges.end(),edges.begin(),edges.end());
    }
    std::sort(patchEdges.begin(),patchEdges.end());
    patchEdges.erase(std::unique(patchEdges.begin(),patchEdges.end()),patchEdges.end());
    std::map<std::size_t,std::vector<std::size_t>> adjacency;
    for (const auto edgeId:patchEdges) {
        if (edgeId>=topology.edges.size()) return std::nullopt;
        const auto& edge=topology.edges[edgeId];
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

[[nodiscard]] std::optional<LocalQualityRank2D> localReplacementQualityRank(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& halo,
    const std::vector<std::size_t>& removedCells,
    const std::vector<Polygon2D>& replacementPieces,
    const TolerancePolicy& tol,const SolverQualityPolicy2D& policy) {
    if (removedCells.empty() || replacementPieces.empty() ||
        !std::is_sorted(removedCells.begin(),removedCells.end()) ||
        !std::is_sorted(halo.begin(),halo.end())) return std::nullopt;
    for (const auto cell:removedCells) {
        if (!std::binary_search(halo.begin(),halo.end(),cell)) return std::nullopt;
    }
    const auto boundary=patchBoundaryRegion(topology,halo,tol);
    if (!boundary) return std::nullopt;
    const Domain2D domain{boundary->bounds()};
    if (!domain.valid(tol)) return std::nullopt;
    std::vector<CutCell2D> cells;
    cells.reserve(halo.size()+replacementPieces.size());
    for (const auto cell:halo) {
        if (cell==removedCells.front()) {
            for (const auto& piece:replacementPieces)
                cells.push_back(makeCell(cells.size(),piece,tol));
        }
        if (std::binary_search(removedCells.begin(),removedCells.end(),cell)) continue;
        if (cell>=topology.cells.size()) return std::nullopt;
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol));
    }
    const auto patch=buildGlobalTopology(cells,domain,*boundary,tol);
    if (!patch.valid()) return std::nullopt;
    const auto quality=evaluateSolverQuality2D(patch,policy,tol);
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
                const Segment2D sourceEdge{a,b};
                if (pointOnSegment(patchA,sourceEdge,tol) &&
                    pointOnSegment(patchB,sourceEdge,tol)) return true;
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
        // A physical-wall split replaces one cell by two pieces.
        if (proposal.second!=proposal.first) removed[proposal.second]=true;
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
            for (const auto& piece:replacements[cell]->pieces) {
                cells.push_back(makeCell(cells.size(),piece,tol,lineage));
                rebuiltImmutable.push_back(false);
            }
        } else {
            cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol,
                                     topology.cells[cell].sourceLineage));
            rebuiltImmutable.push_back(!immutableCells.empty() && immutableCells[cell]);
        }
    }
    const auto topologyStart=ProfileClock::now();
    RepartitionBatch2D rebuilt;
    rebuilt.topology=buildGlobalTopology(
        cells,domain,boundary,tol,topology.constructionRegistry,
        topology.embeddedBoundaryFragments);
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
    SolverTopologyProfile2D* profile,bool candidate,
    const SolverQualityPolicy2D& policy) {
    const auto start=ProfileClock::now();
    auto quality=evaluateSolverQuality2D(topology,policy,tol);
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

} // namespace cartmesh2d::solver_topology_detail
