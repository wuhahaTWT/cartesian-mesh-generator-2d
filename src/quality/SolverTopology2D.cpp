#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"

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

struct PartitionAttempt2D {
    TopologyMesh2D topology;
    std::vector<std::size_t> sourceForCell;
    std::size_t partitionedSourceCellCount = 0;
    std::string error;
};

struct QualityScore2D {
    std::size_t issueCount = 0;
    double maximumSeverity = 0.0;
    double totalSeverity = 0.0;
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

[[nodiscard]] bool internalInterfaceIssue(SolverQualityIssueCode2D code) noexcept {
    return code==SolverQualityIssueCode2D::ExcessiveNonOrthogonality ||
           code==SolverQualityIssueCode2D::ExcessiveSkewness ||
           code==SolverQualityIssueCode2D::LowFaceWeight ||
           code==SolverQualityIssueCode2D::LowVolumeRatio;
}

[[nodiscard]] PartitionAttempt2D partitionSourcePolygons(
    const std::vector<Polygon2D>& sourcePolygons,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    PartitionAttempt2D result;
    std::vector<CutCell2D> cells;
    cells.reserve(sourcePolygons.size());
    for (std::size_t source=0;source<sourcePolygons.size();++source) {
        const auto& polygon=sourcePolygons[source];
        if (strictlyConvex(polygon)) {
            cells.push_back(makeCell(cells.size(),polygon,tol));
            result.sourceForCell.push_back(source);
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
            cells.push_back(makeCell(cells.size(),piece,tol));
            result.sourceForCell.push_back(source);
        }
    }
    result.topology=buildGlobalTopology(cells,domain,boundary,tol);
    if (!result.topology.valid()) result.error="partitioned solver topology audit failed";
    return result;
}

[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> sourceMergePairs(
    const PartitionAttempt2D& partition,const SolverQualityReport2D& quality) {
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    const auto addEdge=[&](const Edge2D& edge) {
        if (!edge.neighbour || edge.owner>=partition.sourceForCell.size() ||
            *edge.neighbour>=partition.sourceForCell.size()) return;
        const std::size_t first=partition.sourceForCell[edge.owner];
        const std::size_t second=partition.sourceForCell[*edge.neighbour];
        if (first!=second) pairs.emplace_back(std::min(first,second),std::max(first,second));
    };
    for (const auto& issue:quality.issues) {
        if (!internalInterfaceIssue(issue.code)) continue;
        if (issue.edgeId<partition.topology.edges.size()) {
            addEdge(partition.topology.edges[issue.edgeId]);
        }
        if (issue.cellId>=partition.sourceForCell.size()) continue;
        for (const auto& edge:partition.topology.edges) {
            if (edge.owner==issue.cellId ||
                (edge.neighbour && *edge.neighbour==issue.cellId)) {
                addEdge(edge);
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
                       const SolverQualityReport2D& quality) {
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
        for (const auto& edge:topology.edges) {
            if (!edge.neighbour) continue;
            if (edge.owner==cell) pairs.emplace_back(std::min(cell,*edge.neighbour),
                                                     std::max(cell,*edge.neighbour));
            else if (*edge.neighbour==cell) pairs.emplace_back(std::min(cell,edge.owner),
                                                                std::max(cell,edge.owner));
        }
    }
    std::sort(pairs.begin(),pairs.end());
    pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    return pairs;
}

[[nodiscard]] TopologyMesh2D repartitionPair(
    const TopologyMesh2D& topology,std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const Domain2D& domain,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol) {
    std::vector<CutCell2D> cells;
    cells.reserve(topology.cells.size());
    for (std::size_t cell=0;cell<topology.cells.size();++cell) {
        if (cell==first) {
            cells.push_back(makeCell(cells.size(),firstPiece,tol));
            cells.push_back(makeCell(cells.size(),secondPiece,tol));
        }
        if (cell==first || cell==second) continue;
        cells.push_back(makeCell(cells.size(),topologyCellPolygon(topology,cell),tol));
    }
    return buildGlobalTopology(cells,domain,boundary,tol);
}

} // namespace

SolverLocalRepartitionResult2D repartitionSolverTopologyByQuality2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    SolverLocalRepartitionResult2D result;
    result.topology=topology;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back("local solver repartition requires valid inputs");
        return result;
    }
    for (std::size_t iteration=0;iteration<32;++iteration) {
        const auto quality=evaluateSolverQuality2D(result.topology,{},tol);
        if (quality.valid()) break;
        const auto pairs=solverRepartitionPairs(result.topology,quality);
        std::optional<TopologyMesh2D> bestTopology;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:pairs) {
            const auto merged=mergeAdjacentPolygonsSimple(
                topologyCellPolygon(result.topology,first),
                topologyCellPolygon(result.topology,second),tol);
            if (!merged) continue;
            for (const auto& [firstPiece,secondPiece]:
                 convexTwoPieceSplits(*merged,domain,boundary,tol)) {
                auto candidate=repartitionPair(result.topology,first,second,
                                               firstPiece,secondPiece,
                                               domain,boundary,tol);
                if (!candidate.valid()) continue;
                const auto candidateQuality=evaluateSolverQuality2D(candidate,{},tol);
                const auto candidateScore=qualityScore(candidateQuality);
                if (betterQualityScore(candidateScore,bestScore)) {
                    bestScore=candidateScore;
                    bestTopology=std::move(candidate);
                }
            }
        }
        if (!bestTopology) break;
        result.topology=std::move(*bestTopology);
        ++result.repartitionCount;
    }
    return result;
}

SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    SolverTopologyResult2D result;
    result.inputCellCount=topology.cells.size();
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back("solver topology repair requires valid inputs");
        return result;
    }
    std::vector<Polygon2D> sourcePolygons;
    sourcePolygons.reserve(topology.cells.size());
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
        sourcePolygons.push_back(std::move(polygon));
    }

    PartitionAttempt2D partition=partitionSourcePolygons(
        sourcePolygons,domain,boundary,tol);
    if (!partition.error.empty()) {
        result.issues.push_back(partition.error);
        return result;
    }
    for (std::size_t iteration=0;iteration<32;++iteration) {
        const auto quality=evaluateSolverQuality2D(partition.topology,{},tol);
        if (quality.valid()) break;
        const auto pairs=sourceMergePairs(partition,quality);
        std::optional<PartitionAttempt2D> bestPartition;
        std::optional<std::vector<Polygon2D>> bestPolygons;
        QualityScore2D bestScore=qualityScore(quality);
        for (const auto& [first,second]:pairs) {
            if (second>=sourcePolygons.size()) continue;
            const auto merged=mergeAdjacentPolygonsSimple(
                sourcePolygons[first],sourcePolygons[second],tol);
            if (!merged) continue;
            auto candidatePolygons=sourcePolygons;
            candidatePolygons[first]=*merged;
            candidatePolygons.erase(
                candidatePolygons.begin()+static_cast<std::ptrdiff_t>(second));
            auto candidate=partitionSourcePolygons(
                candidatePolygons,domain,boundary,tol);
            if (!candidate.error.empty()) continue;
            const auto candidateQuality=evaluateSolverQuality2D(candidate.topology,{},tol);
            const auto candidateScore=qualityScore(candidateQuality);
            if (betterQualityScore(candidateScore,bestScore)) {
                bestScore=candidateScore;
                bestPartition=std::move(candidate);
                bestPolygons=std::move(candidatePolygons);
            }
        }
        if (!bestPartition) break;
        partition=std::move(*bestPartition);
        sourcePolygons=std::move(*bestPolygons);
        ++result.qualityAgglomeratedSourceCellCount;
    }
    auto repartitioned=repartitionSolverTopologyByQuality2D(
        partition.topology,domain,boundary,tol);
    if (!repartitioned.valid()) {
        result.issues.insert(result.issues.end(),repartitioned.issues.begin(),
                             repartitioned.issues.end());
        return result;
    }
    partition.topology=std::move(repartitioned.topology);
    result.qualityRepartitionCount=repartitioned.repartitionCount;
    result.topology=std::move(partition.topology);
    result.partitionedCellCount=partition.partitionedSourceCellCount;
    result.outputCellCount=result.topology.cells.size();
    return result;
}

} // namespace cartmesh2d
