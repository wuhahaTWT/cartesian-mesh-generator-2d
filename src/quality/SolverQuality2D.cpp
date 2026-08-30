#include "cartmesh2d/quality/SolverQuality2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>

namespace cartmesh2d {
namespace {

std::size_t solverQualityEvaluations=0U;

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

SolverCellMetrics2D evaluateSolverCellMetrics2D(const Polygon2D& polygon,
                                                const TolerancePolicy& tol) {
    SolverCellMetrics2D metrics;
    const auto centroid=polygon.vertices.size()>=3?polygon.centroid(tol):std::nullopt;
    if (!centroid || !(polygon.area()>0.0)) return metrics;
    metrics.centroid=*centroid;
    metrics.area=polygon.area();

    double maxEdge=0.0;
    const std::size_t n=polygon.vertices.size();
    for (std::size_t i=0;i<n;++i) {
        const Point2D& previous=polygon.vertices[(i+n-1)%n];
        const Point2D& current=polygon.vertices[i];
        const Point2D& next=polygon.vertices[(i+1)%n];
        const double edgeLength=length(current,next);
        metrics.perimeter+=edgeLength;
        maxEdge=std::max(maxEdge,edgeLength);

        const Vector2D a=previous-current;
        const Vector2D b=next-current;
        const double denom=std::sqrt(squaredNorm(a)*squaredNorm(b));
        if (!(denom>tol.absolute*tol.absolute)) continue;
        const double base=std::acos(std::clamp(dot(a,b)/denom,-1.0,1.0));
        const double angle=orientationSign(previous,current,next)>=0
            ? degrees(base) : 360.0-degrees(base);
        metrics.minInteriorAngleDeg=std::min(metrics.minInteriorAngleDeg,angle);
        metrics.maxConcavityDeg=std::max(metrics.maxConcavityDeg,
                                         std::max(0.0,angle-180.0));
    }
    const double hydraulicRadius=2.0*metrics.area/metrics.perimeter;
    metrics.hydraulicAspect=maxEdge/hydraulicRadius;
    metrics.compactness=4.0*std::numbers::pi*metrics.area/
                        (metrics.perimeter*metrics.perimeter);
    metrics.valid=true;
    return metrics;
}

SolverInternalFaceMetrics2D evaluateSolverInternalFaceMetrics2D(
    const Point2D& a,const Point2D& b,const Point2D& ownerCentroid,
    const Point2D& neighbourCentroid,double ownerArea,double neighbourArea,
    const TolerancePolicy& tol) {
    SolverInternalFaceMetrics2D metrics;
    metrics.length=length(a,b);
    const Vector2D d=neighbourCentroid-ownerCentroid;
    const Vector2D e=b-a;
    const Vector2D normal{e.y,-e.x};
    const double denom=std::sqrt(squaredNorm(d)*squaredNorm(normal));
    if (!(denom>0.0)) return metrics;
    metrics.orientationValid=true;
    metrics.nonOrthogonalityDeg=degrees(
        std::acos(std::clamp(std::abs(dot(d,normal))/denom,0.0,1.0)));

    const double lineDenominator=cross(e,d);
    metrics.skewness=std::numeric_limits<double>::infinity();
    const double scale=std::sqrt(squaredNorm(e)*squaredNorm(d));
    const double parallelEps=tol.absolute*tol.absolute+tol.relative*scale;
    if (std::abs(lineDenominator)>parallelEps) {
        const double t=cross(ownerCentroid-a,d)/lineDenominator;
        const Point2D intersection=a+e*t;
        const Point2D midpoint{0.5*(a.x+b.x),0.5*(a.y+b.y)};
        // Match OpenFOAM's internal-face normalization: the skew vector is
        // scaled by at least 0.2*|cell-centre connector| and otherwise by
        // the approximate face-centre-to-edge distance in the skew
        // direction. For an extruded 2D edge the latter is half its length.
        const double normalization=std::max(0.2*std::sqrt(squaredNorm(d)),
                                            0.5*metrics.length)+
                                   std::numeric_limits<double>::min();
        metrics.skewness=length(intersection,midpoint)/normalization;
    }

    // OpenFOAM face interpolation weight: project each cell-centre to
    // face-centre vector onto the face area vector, then divide the
    // smaller absolute distance by their sum.  The face-area magnitude
    // cancels in two dimensions, so the unnormalised edge normal is exact.
    const Point2D faceCentre{0.5*(a.x+b.x),0.5*(a.y+b.y)};
    const double dOwn=std::abs(dot(normal,faceCentre-ownerCentroid));
    const double dNei=std::abs(dot(normal,neighbourCentroid-faceCentre));
    metrics.faceWeight=std::min(dOwn,dNei)/(dOwn+dNei+
        std::numeric_limits<double>::min());

    // OpenFOAM volume-ratio check is min(Vowner,Vneighbour) / max(...).
    // Extrusion thickness is common and cancels, leaving the 2D areas.
    metrics.volumeRatio=std::min(ownerArea,neighbourArea)/
        (std::max(ownerArea,neighbourArea)+std::numeric_limits<double>::min());
    return metrics;
}

double evaluateSolverBoundaryFaceSkewness2D(const Point2D& a,const Point2D& b,
                                            const Point2D& ownerCentroid) {
    const Point2D faceCentre{0.5*(a.x+b.x),0.5*(a.y+b.y)};
    const Vector2D ownerToFace=faceCentre-ownerCentroid;
    const Vector2D edgeVector=b-a;
    const Vector2D normal{edgeVector.y,-edgeVector.x};
    const double normalMagnitude=std::sqrt(squaredNorm(normal));
    if (!(normalMagnitude>0.0)) return std::numeric_limits<double>::infinity();
    const Vector2D unitNormal=normal*(1.0/normalMagnitude);
    const Vector2D normalProjection=unitNormal*dot(ownerToFace,unitNormal);
    const Vector2D tangential{ownerToFace.x-normalProjection.x,
                              ownerToFace.y-normalProjection.y};
    const double skewMagnitude=std::sqrt(squaredNorm(tangential));
    double faceDirectionDistance=0.0;
    if (skewMagnitude>0.0) {
        const Vector2D skewDirection=tangential*(1.0/skewMagnitude);
        faceDirectionDistance=std::max(
            std::abs(dot(skewDirection,a-faceCentre)),
            std::abs(dot(skewDirection,b-faceCentre)));
    }
    // OpenFOAM's boundary-face check mirrors the owner cell and
    // normalises the tangential correction by the larger of
    // 0.4 times the normal owner-to-face distance and the
    // face-centre-to-edge distance in the skew direction.
    const double normalization=std::max(
        0.4*std::sqrt(squaredNorm(normalProjection)),
        faceDirectionDistance)+std::numeric_limits<double>::min();
    return skewMagnitude/normalization;
}

SolverQualityReport2D evaluateSolverQuality2D(
    const TopologyMesh2D& topology, const SolverQualityPolicy2D& policy,
    const TolerancePolicy& tol) {
    ++solverQualityEvaluations;
    SolverQualityReport2D report;
    report.policy=policy;
    if (!topology.valid()) {
        report.issues.push_back({SolverQualityIssueCode2D::InvalidTopology,0,0,0.0,0.0,
                                 "solver quality requires valid global topology"});
        return report;
    }

    std::vector<Point2D> centroids(topology.cells.size());
    std::vector<bool> centroidValid(topology.cells.size(),false);
    std::vector<double> areas(topology.cells.size(),0.0);
    double minFace=std::numeric_limits<double>::infinity();

    for (const auto& cell:topology.cells) {
        bool valid=false;
        const Polygon2D polygon=cellPolygon(cell,topology,valid);
        const auto metrics=valid?evaluateSolverCellMetrics2D(polygon,tol)
                               :SolverCellMetrics2D{};
        if (!metrics.valid) {
            report.issues.push_back({SolverQualityIssueCode2D::InvalidCell,cell.id,0,0.0,0.0,
                                     "cell polygon or centroid is invalid"});
            continue;
        }
        centroids[cell.id]=metrics.centroid;
        centroidValid[cell.id]=true;
        areas[cell.id]=metrics.area;

        report.maxCellAspect=std::max(report.maxCellAspect,metrics.hydraulicAspect);
        report.maxConcavityDeg=std::max(report.maxConcavityDeg,metrics.maxConcavityDeg);
        report.minInteriorAngleDeg=std::min(report.minInteriorAngleDeg,
                                            metrics.minInteriorAngleDeg);
        report.minCompactness=std::min(report.minCompactness,metrics.compactness);

        if (metrics.hydraulicAspect>policy.maxCellAspect) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveAspect,cell.id,0,
                                     metrics.hydraulicAspect,policy.maxCellAspect,
                                     "cell hydraulic aspect exceeds limit"});
        }
        if (metrics.maxConcavityDeg>policy.maxConcavityDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveConcavity,cell.id,0,
                                     metrics.maxConcavityDeg,policy.maxConcavityDeg,
                                     "cell concavity exceeds limit"});
        }
        if (metrics.minInteriorAngleDeg<policy.minInteriorAngleDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::SmallInteriorAngle,cell.id,0,
                                     metrics.minInteriorAngleDeg,policy.minInteriorAngleDeg,
                                     "cell interior angle is below limit"});
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
        if (edge.owner>=centroids.size() || !centroidValid[edge.owner]) continue;

        if (!edge.neighbour) {
            const double skewness=evaluateSolverBoundaryFaceSkewness2D(
                a,b,centroids[edge.owner]);
            report.maxBoundarySkewness=std::max(report.maxBoundarySkewness,skewness);
            if (skewness>policy.maxBoundarySkewness) {
                report.issues.push_back({SolverQualityIssueCode2D::ExcessiveBoundarySkewness,
                                         edge.owner,edge.id,skewness,
                                         policy.maxBoundarySkewness,
                                         "boundary face skewness exceeds limit"});
            }
            continue;
        }
        if (*edge.neighbour>=centroids.size() || !centroidValid[*edge.neighbour]) continue;

        const auto metrics=evaluateSolverInternalFaceMetrics2D(
            a,b,centroids[edge.owner],centroids[*edge.neighbour],
            areas[edge.owner],areas[*edge.neighbour],tol);
        if (!metrics.orientationValid) continue;
        report.maxNonOrthogonalityDeg=std::max(report.maxNonOrthogonalityDeg,
                                               metrics.nonOrthogonalityDeg);
        if (metrics.nonOrthogonalityDeg>policy.maxNonOrthogonalityDeg) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveNonOrthogonality,
                                     edge.owner,edge.id,metrics.nonOrthogonalityDeg,
                                     policy.maxNonOrthogonalityDeg,
                                     "internal face non-orthogonality exceeds limit"});
        }

        report.maxInternalSkewness=std::max(report.maxInternalSkewness,metrics.skewness);
        if (metrics.skewness>policy.maxInternalSkewness) {
            report.issues.push_back({SolverQualityIssueCode2D::ExcessiveSkewness,
                                     edge.owner,edge.id,metrics.skewness,
                                     policy.maxInternalSkewness,
                                     "internal face skewness exceeds limit"});
        }

        report.minFaceWeight=std::min(report.minFaceWeight,metrics.faceWeight);
        if (metrics.faceWeight<policy.minFaceWeight) {
            report.issues.push_back({SolverQualityIssueCode2D::LowFaceWeight,
                                     edge.owner,edge.id,metrics.faceWeight,
                                     policy.minFaceWeight,
                                     "internal face interpolation weight is below limit"});
        }

        report.minVolumeRatio=std::min(report.minVolumeRatio,metrics.volumeRatio);
        if (metrics.volumeRatio<policy.minVolumeRatio) {
            report.issues.push_back({SolverQualityIssueCode2D::LowVolumeRatio,
                                     edge.owner,edge.id,metrics.volumeRatio,
                                     policy.minVolumeRatio,
                                     "neighbouring cell volume ratio is below limit"});
        }
    }
    report.minFaceLength=std::isfinite(minFace)?minFace:0.0;
    return report;
}

std::size_t solverQualityEvaluationCount2D() noexcept {
    return solverQualityEvaluations;
}

std::string solverQualityReportToJson(const SolverQualityReport2D& report,
                                      int indentSpaces) {
    const int step=std::max(0,indentSpaces);
    const std::string i1=indent(step);
    const std::string i2=indent(step*2);
    std::ostringstream out;
    out<<std::setprecision(17);
    out<<"{\n";
    out<<i1<<"\"quality_class\": \"solver_quality\",\n";
    out<<i1<<"\"metric_semantics\": \"native 2D solver-topology metrics; OpenFOAM-like formulas are not external checkMesh evidence\",\n";
    out<<i1<<"\"valid\": "<<(report.valid()?"true":"false")<<",\n";
    out<<i1<<"\"policy\": {\n";
    out<<i2<<"\"max_non_orthogonality_deg\": "<<report.policy.maxNonOrthogonalityDeg<<",\n";
    out<<i2<<"\"max_internal_skewness\": "<<report.policy.maxInternalSkewness<<",\n";
    out<<i2<<"\"max_boundary_skewness\": "<<report.policy.maxBoundarySkewness<<",\n";
    out<<i2<<"\"max_concavity_deg\": "<<report.policy.maxConcavityDeg<<",\n";
    out<<i2<<"\"max_cell_aspect\": "<<report.policy.maxCellAspect<<",\n";
    out<<i2<<"\"min_interior_angle_deg\": "<<report.policy.minInteriorAngleDeg<<",\n";
    out<<i2<<"\"min_face_length\": "<<report.policy.minFaceLength<<",\n";
    out<<i2<<"\"min_face_weight\": "<<report.policy.minFaceWeight<<",\n";
    out<<i2<<"\"min_volume_ratio\": "<<report.policy.minVolumeRatio<<"\n";
    out<<i1<<"},\n";
    out<<i1<<"\"metrics\": {\n";
    out<<i2<<"\"max_non_orthogonality_deg\": "<<report.maxNonOrthogonalityDeg<<",\n";
    out<<i2<<"\"max_internal_skewness\": "<<report.maxInternalSkewness<<",\n";
    out<<i2<<"\"max_boundary_skewness\": "<<report.maxBoundarySkewness<<",\n";
    out<<i2<<"\"max_concavity_deg\": "<<report.maxConcavityDeg<<",\n";
    out<<i2<<"\"max_cell_aspect\": "<<report.maxCellAspect<<",\n";
    out<<i2<<"\"min_interior_angle_deg\": "<<report.minInteriorAngleDeg<<",\n";
    out<<i2<<"\"min_face_length\": "<<report.minFaceLength<<",\n";
    out<<i2<<"\"min_face_weight\": "<<report.minFaceWeight<<",\n";
    out<<i2<<"\"min_volume_ratio\": "<<report.minVolumeRatio<<",\n";
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
