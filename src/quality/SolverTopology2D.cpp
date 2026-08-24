#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <numbers>
#include <optional>
#include <tuple>

namespace cartmesh2d {
namespace {

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
        if (orientationSign(polygon.vertices[(i+n-1)%n],polygon.vertices[i],
                            polygon.vertices[(i+1)%n])<=0) return false;
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

[[nodiscard]] CutCell2D makeCell(std::size_t id,const Polygon2D& polygon,
                                 const TolerancePolicy& tol) {
    CutCell2D cell;
    cell.sourceId=id;
    cell.sourceKey=id;
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

} // namespace

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    SolverTopologyResult2D result;
    result.inputCellCount=topology.cells.size();
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back("solver topology repair requires valid inputs");
        return result;
    }
    std::vector<CutCell2D> cells;
    cells.reserve(topology.cells.size());
    for (const auto& cell:topology.cells) {
        Polygon2D polygon;
        polygon.vertices.reserve(cell.vertices.size());
        for (const auto vertex:cell.vertices) {
            if (vertex>=topology.vertices.size()) {
                result.issues.push_back("solver topology cell has invalid vertex index");
                return result;
            }
            polygon.vertices.push_back(topology.vertices[vertex].point);
        }
        if (strictlyConvex(polygon)) {
            cells.push_back(makeCell(cells.size(),polygon,tol));
            continue;
        }
        const auto triangles=triangulate(polygon,domain,boundary,tol);
        if (!triangles) {
            result.issues.push_back("deterministic convex partition failed for cell "+
                                    std::to_string(cell.id));
            return result;
        }
        ++result.partitionedCellCount;
        for (const auto& triangle:*triangles) cells.push_back(makeCell(cells.size(),triangle,tol));
    }
    result.topology=buildGlobalTopology(cells,domain,boundary,tol);
    result.outputCellCount=result.topology.cells.size();
    if (!result.topology.valid()) result.issues.push_back("partitioned solver topology audit failed");
    return result;
}

} // namespace cartmesh2d
