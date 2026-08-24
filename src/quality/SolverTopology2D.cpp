#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
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

[[nodiscard]] std::optional<Polygon2D> mergeAdjacentTriangles(
    const Polygon2D& first,const Polygon2D& second) {
    std::vector<Point2D> unique=first.vertices;
    for (const auto& point:second.vertices) {
        const bool exists=std::any_of(unique.begin(),unique.end(),[&](const Point2D& other) {
            return point.x==other.x && point.y==other.y;
        });
        if (!exists) unique.push_back(point);
    }
    if (unique.size()!=4) return std::nullopt;
    Point2D center{};
    for (const auto& point:unique) { center.x+=point.x; center.y+=point.y; }
    center.x/=4.0;
    center.y/=4.0;
    std::sort(unique.begin(),unique.end(),[&](const Point2D& lhs,const Point2D& rhs) {
        const double a=std::atan2(lhs.y-center.y,lhs.x-center.x);
        const double b=std::atan2(rhs.y-center.y,rhs.x-center.x);
        return a==b?std::tie(lhs.x,lhs.y)<std::tie(rhs.x,rhs.y):a<b;
    });
    Polygon2D merged{unique};
    if (merged.signedArea()<0.0) std::reverse(merged.vertices.begin(),merged.vertices.end());
    return strictlyConvex(merged)?std::optional<Polygon2D>(merged):std::nullopt;
}

[[nodiscard]] std::optional<std::vector<Polygon2D>> triangulate(
    const Polygon2D& polygon,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    if (polygon.vertices.size()<3 || !(polygon.signedArea()>0.0)) return std::nullopt;
    std::vector<std::size_t> remaining(polygon.vertices.size());
    std::iota(remaining.begin(),remaining.end(),0);
    std::vector<Polygon2D> triangles;
    triangles.reserve(polygon.vertices.size()-2);
    while (remaining.size()>3) {
        std::optional<std::size_t> selected;
        unsigned selectedBoundaryEdges=3;
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
            if (!selected || boundaryEdges<selectedBoundaryEdges) {
                selected=pos;
                selectedBoundaryEdges=boundaryEdges;
            }
            if (boundaryEdges==0) break;
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
            const auto merged=mergeAdjacentTriangles(triangles[i],triangles[j]);
            if (!merged) continue;
            triangles[i]=*merged;
            triangles.erase(triangles.begin()+static_cast<std::ptrdiff_t>(j));
            if (j<i) --i;
            break;
        }
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
