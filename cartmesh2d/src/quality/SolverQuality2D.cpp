#include "cartmesh2d/quality/SolverQuality2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>

namespace cartmesh2d {
namespace {

[[nodiscard]] double degrees(double radians) noexcept {
    return radians * 180.0 / std::numbers::pi;
}

[[nodiscard]] double length(const Point2D& a, const Point2D& b) noexcept {
    return std::hypot(b.x-a.x,b.y-a.y);
}

[[nodiscard]] Polygon2D cellPolygon(const TopologyCell2D& cell,
                                    const TopologyMesh2D& topology,
                                    bool& valid) {
    Polygon2D polygon;
    polygon.vertices.reserve(cell.vertices.size());
    for (const auto vertexId:cell.vertices) {
        if (vertexId>=topology.vertices.size()) {
            valid=false;
            return polygon;
        }
        polygon.vertices.push_back(topology.vertices[vertexId].point);
    }
    valid=polygon.vertices.size()>=3;
    return polygon;
}

[[nodiscard]] std::string indent(int count) {
    return std::string(static_cast<std::size_t>(std::max(0,count)),' ');
}

} // namespace

SolverQualityReport2D evaluateSolverQuality2D(
    const TopologyMesh2D& topology, const SolverQualityPolicy2D& policy,
    const TolerancePolicy& tol) {
    SolverQualityReport2D report;
    report.policy=policy;
    if (!topology.valid()) {
        report.issues.push_back({SolverQualityIssueCode2D::InvalidTopology,0,0,0.0,0.0,
                                 "solver quality requires valid global topology"});
        return report;
    }

    std::vector<Point2D> centroids(topology.cells.size());
    std::vector<bool> centroidValid(topology.cells.size(),false);
    double minFace=std::numeric_limits<double>::infinity();

    for (const auto& cell:topology.cells) {
        bool valid=false;
        const Polygon2D polygon=cellPolygon(cell,topology,valid);
        const auto centroid=valid?polygon.centroid(tol):std::nullopt;
        if (!valid || !centroid || !(polygon.area()>0.0)) {
            report.issues.push_back({SolverQualityIssueCode2D::InvalidCell,cell.id,0,0.0,0.0,
                                     "cell polygon or centroid is invalid"});
            continue;
        }
        centroids[cell.id]=*centroid;
        centroidValid[cell.id]=true;

        double perimeter=0.0;
        double maxEdge=0.0;
        double minAngle=360.0;
        double maxConcavity=0.0;
        const std::size_t n=polygon.vertices.size();
        for (std::size_t i=0;i<n;++i) {
            const Point2D& previous=polygon.vertices[(i+n-1)%n];
            const Point2D& current=polygon.vertices[i];
            const Point2D& next=polygon.vertices[(i+1)%n];
            const double edgeLength=length(current,next);
            perimeter+=edgeLength;
            maxEdge=std::max(maxEdge,edgeLength);

            const Vector2D a=previous-current;
            const Vector2D b=next-current;
            const double denom=std::sqrt(squaredNorm(a)*squaredNorm(b));
            if (!(denom>tol.absolute*tol.absolute)) continue;
            const double base=std::acos(std::clamp(dot(a,b)/denom,-1.0,1.0));
            const double angle=orientationSign(previous,current,next)>=0
                ? degrees(base) : 360.0-degrees(base);
            minAngle=std::min(minAngle,angle);
            maxConcavity=std::max(maxConcavity,std::max(0.0,angle-180.0));
        }

        const double area=polygon.area();
        const double hydraulicRadius=2.0*area/perimeter;
        const double aspect=maxEdge/hydraulicRadius;
        const double compactness=4.0*std::numbers::pi*area/(perimeter*perimeter);
        report.maxCellAspect=std::max(report.maxCellAspect,aspect);
        report.maxConcavityDeg=std::max(report.maxConcavityDeg,maxConcavity);
        report.minInteriorAngleDeg=std::min(report.minInteriorAngleDeg,minAngle);
        report.minCompactness=std::min(report.minCompactness,compactness);

        if (aspect>policy.maxCellAspect) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveAspect,cell.id,0,
                                     aspect,policy.maxCellAspect,"cell hydraulic aspect exceeds limit"});
        }
        if (maxConcavity>policy.maxConcavityDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveConcavity,cell.id,0,
                                     maxConcavity,policy.maxConcavityDeg,"cell concavity exceeds limit"});
        }
        if (minAngle<policy.minInteriorAngleDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::SmallInteriorAngle,cell.id,0,
                                     minAngle,policy.minInteriorAngleDeg,"cell interior angle is below limit"});
        }
    }

    for (const auto& edge:topology.edges) {
        if (edge.v0>=topology.vertices.size() || edge.v1>=topology.vertices.size()) continue;
        const Point2D& a=topology.vertices[edge.v0].point;
        const Point2D& b=topology.vertices[edge.v1].point;
        const double edgeLength=length(a,b);
        minFace=std::min(minFace,edgeLength);
        if (edgeLength<=policy.minFaceLength) {
            report.issues.push_back({SolverQualityIssueCode2D::ShortFace,edge.owner,edge.id,
                                     edgeLength,policy.minFaceLength,"face length is below solver limit"});
        }
        if (!edge.neighbour || edge.owner>=centroids.size() ||
            *edge.neighbour>=centroids.size() || !centroidValid[edge.owner] ||
            !centroidValid[*edge.neighbour]) continue;

        const Point2D owner=centroids[edge.owner];
        const Point2D neighbour=centroids[*edge.neighbour];
        const Vector2D d=neighbour-owner;
        const Vector2D e=b-a;
        const Vector2D normal{e.y,-e.x};
        const double denom=std::sqrt(squaredNorm(d)*squaredNorm(normal));
        if (!(denom>0.0)) continue;
        const double nonOrth=degrees(std::acos(std::clamp(std::abs(dot(d,normal))/denom,0.0,1.0)));
        report.maxNonOrthogonalityDeg=std::max(report.maxNonOrthogonalityDeg,nonOrth);
        if (nonOrth>policy.maxNonOrthogonalityDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveNonOrthogonality,
                                     edge.owner,edge.id,nonOrth,policy.maxNonOrthogonalityDeg,
                                     "internal face non-orthogonality exceeds limit"});
        }

        const double lineDenominator=cross(e,d);
        double skewness=std::numeric_limits<double>::infinity();
        const double scale=std::sqrt(squaredNorm(e)*squaredNorm(d));
        const double parallelEps=tol.absolute*tol.absolute+tol.relative*scale;
        if (std::abs(lineDenominator)>parallelEps) {
            const double t=cross(owner-a,d)/lineDenominator;
            const Point2D intersection=a+e*t;
            const Point2D midpoint{0.5*(a.x+b.x),0.5*(a.y+b.y)};
            skewness=length(intersection,midpoint)/edgeLength;
        }
        report.maxInternalSkewness=std::max(report.maxInternalSkewness,skewness);
        if (skewness>policy.maxInternalSkewness) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveSkewness,
                                     edge.owner,edge.id,skewness,policy.maxInternalSkewness,
                                     "internal face skewness exceeds limit"});
        }
    }
    report.minFaceLength=std::isfinite(minFace)?minFace:0.0;
    return report;
}

std::string solverQualityReportToJson(const SolverQualityReport2D& report,
                                      int indentSpaces) {
    const int step=std::max(0,indentSpaces);
    const std::string i1=indent(step);
    const std::string i2=indent(step*2);
    std::ostringstream out;
    out<<std::setprecision(17);
    out<<"{\n";
    out<<i1<<"\"valid\": "<<(report.valid()?"true":"false")<<",\n";
    out<<i1<<"\"policy\": {\n";
    out<<i2<<"\"max_non_orthogonality_deg\": "<<report.policy.maxNonOrthogonalityDeg<<",\n";
    out<<i2<<"\"max_internal_skewness\": "<<report.policy.maxInternalSkewness<<",\n";
    out<<i2<<"\"max_concavity_deg\": "<<report.policy.maxConcavityDeg<<",\n";
    out<<i2<<"\"max_cell_aspect\": "<<report.policy.maxCellAspect<<",\n";
    out<<i2<<"\"min_interior_angle_deg\": "<<report.policy.minInteriorAngleDeg<<",\n";
    out<<i2<<"\"min_face_length\": "<<report.policy.minFaceLength<<"\n";
    out<<i1<<"},\n";
    out<<i1<<"\"metrics\": {\n";
    out<<i2<<"\"max_non_orthogonality_deg\": "<<report.maxNonOrthogonalityDeg<<",\n";
    out<<i2<<"\"max_internal_skewness\": "<<report.maxInternalSkewness<<",\n";
    out<<i2<<"\"max_concavity_deg\": "<<report.maxConcavityDeg<<",\n";
    out<<i2<<"\"max_cell_aspect\": "<<report.maxCellAspect<<",\n";
    out<<i2<<"\"min_interior_angle_deg\": "<<report.minInteriorAngleDeg<<",\n";
    out<<i2<<"\"min_face_length\": "<<report.minFaceLength<<",\n";
    out<<i2<<"\"min_compactness\": "<<report.minCompactness<<"\n";
    out<<i1<<"},\n";
    out<<i1<<"\"issue_count\": "<<report.issues.size()<<",\n";
    out<<i1<<"\"issues\": [";
    for (std::size_t i=0;i<report.issues.size();++i) {
        const auto& issue=report.issues[i];
        out<<(i==0?"\n":",\n")<<i2<<"{\"code\": "<<static_cast<int>(issue.code)
           <<", \"cell_id\": "<<issue.cellId<<", \"edge_id\": "<<issue.edgeId
           <<", \"measured\": "<<issue.measured<<", \"limit\": "<<issue.limit
           <<", \"message\": \""<<issue.message<<"\"}";
    }
    if (!report.issues.empty()) out<<'\n'<<i1;
    out<<"]\n}\n";
    return out.str();
}

} // namespace cartmesh2d
