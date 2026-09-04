#include "cartmesh2d/quality/QualityContract2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>

namespace cartmesh2d {
namespace {

[[nodiscard]] bool ordinary(QualityCellType2D type) noexcept {
    return type==QualityCellType2D::Cartesian ||
           type==QualityCellType2D::RemainderCut ||
           type==QualityCellType2D::Transition ||
           type==QualityCellType2D::Termination;
}

[[nodiscard]] const OrdinaryCellQualityLimits2D& limitsFor(
    const QualityContract2D& contract,QualityCellType2D type) noexcept {
    switch (type) {
    case QualityCellType2D::Cartesian: return contract.cartesian;
    case QualityCellType2D::RemainderCut: return contract.remainderCut;
    case QualityCellType2D::Transition: return contract.transition;
    case QualityCellType2D::Termination: return contract.termination;
    case QualityCellType2D::BoundaryLayer:
    case QualityCellType2D::Unknown: return contract.cartesian;
    }
    return contract.cartesian;
}

[[nodiscard]] Polygon2D cellPolygon(const TopologyMesh2D& topology,
                                    const TopologyCell2D& cell,bool& valid) {
    Polygon2D polygon;
    polygon.vertices.reserve(cell.vertices.size());
    for (const auto vertex:cell.vertices) {
        if (vertex>=topology.vertices.size()) {
            valid=false;
            return {};
        }
        polygon.vertices.push_back(topology.vertices[vertex].point);
    }
    valid=polygon.vertices.size()>=3U;
    return polygon;
}

[[nodiscard]] double pointDistance(const Point2D& a,const Point2D& b) noexcept {
    return std::hypot(b.x-a.x,b.y-a.y);
}

[[nodiscard]] double degrees(double radians) noexcept {
    return radians*180.0/std::numbers::pi;
}

[[nodiscard]] double percentile(const std::vector<double>& sorted,double fraction) {
    if (sorted.empty()) return 0.0;
    const double position=fraction*static_cast<double>(sorted.size()-1U);
    const auto lower=static_cast<std::size_t>(std::floor(position));
    const auto upper=static_cast<std::size_t>(std::ceil(position));
    const double weight=position-static_cast<double>(lower);
    return sorted[lower]*(1.0-weight)+sorted[upper]*weight;
}

[[nodiscard]] QualityMetricSummary2D summarize(
    const std::vector<QualityMetricSample2D>& samples,bool lowerIsWorse) {
    QualityMetricSummary2D result;
    result.count=samples.size();
    result.lowerIsWorse=lowerIsWorse;
    if (samples.empty()) return result;
    std::vector<double> values;
    values.reserve(samples.size());
    for (const auto& sample:samples) values.push_back(sample.value);
    std::sort(values.begin(),values.end());
    result.p50=percentile(values,0.50);
    result.p95=percentile(values,0.95);
    result.p99=percentile(values,0.99);
    const auto selected=lowerIsWorse
        ?std::min_element(samples.begin(),samples.end(),
            [](const auto& a,const auto& b){return a.value<b.value;})
        :std::max_element(samples.begin(),samples.end(),
            [](const auto& a,const auto& b){return a.value<b.value;});
    result.worst=selected->value;
    result.worstEntity=selected->entity;
    return result;
}

[[nodiscard]] bool violates(double value,double limit,bool lowerBound) noexcept {
    return lowerBound?value<limit:value>limit;
}

void checkSample(QualityContractReport2D& report,const std::string& metric,
                 const QualityMetricSample2D& sample,const QualityLimit2D& limit) {
    if (!ordinary(sample.entity.cellType)) return;
    auto& typed=report.byCellType[sample.entity.cellType];
    if (violates(sample.value,limit.hard,limit.lowerBound)) {
        ++typed.hardIssueCount;
        report.issues.push_back({QualityContractLevel2D::Hard,metric,sample.value,
                                 limit.hard,sample.entity});
    } else if (violates(sample.value,limit.preferred,limit.lowerBound)) {
        ++typed.preferredIssueCount;
        report.issues.push_back({QualityContractLevel2D::Preferred,metric,sample.value,
                                 limit.preferred,sample.entity});
    }
}

[[nodiscard]] std::string indent(int count) {
    return std::string(static_cast<std::size_t>(std::max(0,count)),' ');
}

void writeEntity(std::ostream& out,const QualityEntity2D& entity,
                 const std::string& indentation) {
    out<<"{\"cell_id\": "<<entity.cellId<<", \"edge_id\": ";
    if (entity.edgeId) out<<*entity.edgeId; else out<<"null";
    out<<", \"coordinates\": ["<<entity.coordinates.x<<", "
       <<entity.coordinates.y<<"], \"cell_type\": \""
       <<qualityCellTypeName(entity.cellType)<<"\", \"source_id\": "
       <<entity.sourceId<<", \"owner\": "<<entity.owner
       <<", \"neighbour\": ";
    if (entity.neighbour) out<<*entity.neighbour; else out<<"null";
    out<<", \"local_h\": "<<entity.localH<<'}';
    (void)indentation;
}

void writeLimit(std::ostream& out,const QualityLimit2D& limit) {
    out<<"{\"preferred\": "<<limit.preferred<<", \"hard\": "<<limit.hard
       <<", \"direction\": \""<<(limit.lowerBound?"min":"max")<<"\"}";
}

void writeOrdinaryLimits(std::ostream& out,const OrdinaryCellQualityLimits2D& limits,
                         const std::string& i) {
    const auto item=[&](const char* name,const QualityLimit2D& value,bool comma) {
        out<<i<<'"'<<name<<"\": ";
        writeLimit(out,value);
        out<<(comma?",\n":"\n");
    };
    out<<"{\n";
    item("non_orthogonality_deg",limits.nonOrthogonalityDeg,true);
    item("skewness",limits.skewness,true);
    item("face_weight",limits.faceWeight,true);
    item("volume_ratio",limits.volumeRatio,true);
    item("background_volume_ratio",limits.backgroundVolumeRatio,true);
    item("minimum_interior_angle_deg",limits.minimumInteriorAngleDeg,true);
    item("hydraulic_aspect",limits.hydraulicAspect,true);
    item("face_length_over_local_background_h",limits.faceOverLocalBackgroundH,true);
    item("face_length_over_sqrt_owner_area",limits.faceOverSqrtOwnerArea,true);
    item("face_length_over_sqrt_neighbour_area",limits.faceOverSqrtNeighbourArea,false);
    out<<i.substr(0,i.size()>=2U?i.size()-2U:0U)<<'}';
}

} // namespace

QualityContractStatus2D QualityContractReport2D::status() const noexcept {
    if (!validInput()) return QualityContractStatus2D::Fail;
    bool warn=false;
    for (const auto& [type,result]:byCellType) {
        (void)type;
        if (result.status()==QualityContractStatus2D::Fail) {
            return QualityContractStatus2D::Fail;
        }
        warn=warn || result.status()==QualityContractStatus2D::Warn;
    }
    return warn?QualityContractStatus2D::Warn:QualityContractStatus2D::Pass;
}

QualityContractReport2D evaluateQualityContract2D(
    const TopologyMesh2D& topology,
    const std::vector<QualityCellMetadata2D>& metadata,
    const BoundaryLayerQualitySamples2D& boundaryLayer,
    const QualityContract2D& contract,
    const SolverQualityReport2D* legacyHardSafety,
    const TolerancePolicy& tol) {
    QualityContractReport2D report;
    report.contract=contract;
    report.legacyHardSafety=legacyHardSafety?*legacyHardSafety:
        evaluateSolverQuality2D(topology,{},tol);
    if (!topology.valid()) report.inputIssues.push_back("quality contract requires valid topology");
    if (metadata.size()!=topology.cells.size()) {
        report.inputIssues.push_back("quality contract metadata count must equal cell count");
        return report;
    }

    std::vector<Polygon2D> polygons(topology.cells.size());
    std::vector<Point2D> centroids(topology.cells.size());
    std::vector<double> areas(topology.cells.size(),0.0);
    std::vector<bool> valid(topology.cells.size(),false);
    std::map<std::string,std::vector<QualityMetricSample2D>> samples;
    const auto entityForCell=[&](std::size_t cellId,Point2D point)->QualityEntity2D {
        return {cellId,std::nullopt,point,metadata[cellId].type,
                metadata[cellId].sourceId,cellId,std::nullopt,
                metadata[cellId].localBackgroundH};
    };
    for (const auto& cell:topology.cells) {
        bool polygonValid=false;
        polygons[cell.id]=cellPolygon(topology,cell,polygonValid);
        const auto centroid=polygonValid?polygons[cell.id].centroid(tol):std::nullopt;
        const double area=polygonValid?polygons[cell.id].area():0.0;
        if (!polygonValid || !centroid || !(area>0.0)) {
            report.inputIssues.push_back("quality contract encountered invalid cell geometry");
            continue;
        }
        if (metadata[cell.id].type==QualityCellType2D::Unknown) {
            report.inputIssues.push_back("quality contract encountered unknown cell type");
        }
        if (ordinary(metadata[cell.id].type) && !(metadata[cell.id].localBackgroundH>0.0)) {
            report.inputIssues.push_back("ordinary quality cell lacks positive local background h");
        }
        centroids[cell.id]=*centroid;
        areas[cell.id]=area;
        valid[cell.id]=true;
        auto& typedResult=report.byCellType[metadata[cell.id].type];
        ++typedResult.cellCount;
        if (metadata[cell.id].type==QualityCellType2D::BoundaryLayer) {
            // Q1 records the six layer-specific metrics but deliberately does
            // not invent acceptance thresholds that the requested contract
            // did not specify.
            typedResult.rated=false;
        }
        if (!ordinary(metadata[cell.id].type)) continue;

        double perimeter=0.0;
        double maxEdge=0.0;
        double minAngle=360.0;
        const auto& points=polygons[cell.id].vertices;
        for (std::size_t i=0;i<points.size();++i) {
            const Point2D& previous=points[(i+points.size()-1U)%points.size()];
            const Point2D& current=points[i];
            const Point2D& next=points[(i+1U)%points.size()];
            const double edgeLength=pointDistance(current,next);
            perimeter+=edgeLength;
            maxEdge=std::max(maxEdge,edgeLength);
            const Vector2D a=previous-current;
            const Vector2D b=next-current;
            const double denominator=std::sqrt(squaredNorm(a)*squaredNorm(b));
            if (!(denominator>0.0)) continue;
            const double base=std::acos(std::clamp(dot(a,b)/denominator,-1.0,1.0));
            const double angle=orientationSign(previous,current,next)>=0
                ?degrees(base):360.0-degrees(base);
            minAngle=std::min(minAngle,angle);
        }
        const auto entity=entityForCell(cell.id,*centroid);
        samples["hydraulic_aspect"].push_back(
            {maxEdge/(2.0*area/perimeter),entity});
        samples["minimum_interior_angle_deg"].push_back({minAngle,entity});
    }

    for (const auto& edge:topology.edges) {
        if (edge.v0>=topology.vertices.size() || edge.v1>=topology.vertices.size() ||
            edge.owner>=metadata.size() || !valid[edge.owner]) continue;
        const Point2D a=topology.vertices[edge.v0].point;
        const Point2D b=topology.vertices[edge.v1].point;
        const Point2D midpoint{0.5*(a.x+b.x),0.5*(a.y+b.y)};
        const double faceLength=pointDistance(a,b);
        std::optional<std::size_t> evaluationCell;
        if (ordinary(metadata[edge.owner].type)) evaluationCell=edge.owner;
        if (!evaluationCell && edge.neighbour && *edge.neighbour<metadata.size() &&
            ordinary(metadata[*edge.neighbour].type)) evaluationCell=*edge.neighbour;
        if (!evaluationCell) continue;
        QualityEntity2D edgeEntity{*evaluationCell,edge.id,midpoint,
            metadata[*evaluationCell].type,metadata[*evaluationCell].sourceId,
            edge.owner,edge.neighbour,metadata[*evaluationCell].localBackgroundH};
        double localH=metadata[*evaluationCell].localBackgroundH;
        if (ordinary(metadata[edge.owner].type)) {
            localH=std::max(localH,metadata[edge.owner].localBackgroundH);
        }
        if (edge.neighbour && *edge.neighbour<metadata.size() &&
            ordinary(metadata[*edge.neighbour].type)) {
            localH=std::max(localH,metadata[*edge.neighbour].localBackgroundH);
        }
        edgeEntity.localH=localH;
        samples["face_length_over_local_background_h"].push_back(
            {faceLength/localH,edgeEntity});

        if (ordinary(metadata[edge.owner].type)) {
            auto ownerEntity=edgeEntity;
            ownerEntity.cellId=edge.owner;
            ownerEntity.cellType=metadata[edge.owner].type;
            ownerEntity.sourceId=metadata[edge.owner].sourceId;
            ownerEntity.localH=metadata[edge.owner].localBackgroundH;
            samples["face_length_over_sqrt_owner_area"].push_back(
                {faceLength/std::sqrt(areas[edge.owner]),ownerEntity});
        }
        if (edge.neighbour && *edge.neighbour<metadata.size() && valid[*edge.neighbour] &&
            ordinary(metadata[*edge.neighbour].type)) {
            auto neighbourEntity=edgeEntity;
            neighbourEntity.cellId=*edge.neighbour;
            neighbourEntity.cellType=metadata[*edge.neighbour].type;
            neighbourEntity.sourceId=metadata[*edge.neighbour].sourceId;
            neighbourEntity.localH=metadata[*edge.neighbour].localBackgroundH;
            samples["face_length_over_sqrt_neighbour_area"].push_back(
                {faceLength/std::sqrt(areas[*edge.neighbour]),neighbourEntity});
        }

        const Vector2D edgeVector=b-a;
        const Vector2D normal{edgeVector.y,-edgeVector.x};
        if (!edge.neighbour) {
            const Vector2D ownerToFace=midpoint-centroids[edge.owner];
            const double normalMagnitude=std::sqrt(squaredNorm(normal));
            const Vector2D unitNormal=normal*(1.0/normalMagnitude);
            const Vector2D projection=unitNormal*dot(ownerToFace,unitNormal);
            const Vector2D tangent{ownerToFace.x-projection.x,
                                   ownerToFace.y-projection.y};
            const double skewMagnitude=std::sqrt(squaredNorm(tangent));
            double faceDirectionDistance=0.0;
            if (skewMagnitude>0.0) {
                const Vector2D direction=tangent*(1.0/skewMagnitude);
                faceDirectionDistance=std::max(
                    std::abs(dot(direction,a-midpoint)),
                    std::abs(dot(direction,b-midpoint)));
            }
            const double normalization=std::max(
                0.4*std::sqrt(squaredNorm(projection)),faceDirectionDistance)+
                std::numeric_limits<double>::min();
            samples["skewness"].push_back({skewMagnitude/normalization,edgeEntity});
            continue;
        }
        if (*edge.neighbour>=metadata.size() || !valid[*edge.neighbour]) continue;
        const Point2D owner=centroids[edge.owner];
        const Point2D neighbour=centroids[*edge.neighbour];
        const Vector2D d=neighbour-owner;
        const double denominator=std::sqrt(squaredNorm(d)*squaredNorm(normal));
        if (!(denominator>0.0)) continue;
        const double nonOrth=degrees(std::acos(std::clamp(
            std::abs(dot(d,normal))/denominator,0.0,1.0)));
        samples["non_orthogonality_deg"].push_back({nonOrth,edgeEntity});

        const double lineDenominator=cross(edgeVector,d);
        double skewness=std::numeric_limits<double>::infinity();
        const double scale=std::sqrt(squaredNorm(edgeVector)*squaredNorm(d));
        const double parallelEps=tol.absolute*tol.absolute+tol.relative*scale;
        if (std::abs(lineDenominator)>parallelEps) {
            const double t=cross(owner-a,d)/lineDenominator;
            const Point2D intersection=a+edgeVector*t;
            const double normalization=std::max(0.2*std::sqrt(squaredNorm(d)),
                                                0.5*faceLength)+
                                       std::numeric_limits<double>::min();
            skewness=pointDistance(intersection,midpoint)/normalization;
        }
        samples["skewness"].push_back({skewness,edgeEntity});
        const double ownerDistance=std::abs(dot(normal,midpoint-owner));
        const double neighbourDistance=std::abs(dot(normal,neighbour-midpoint));
        samples["face_weight"].push_back(
            {std::min(ownerDistance,neighbourDistance)/
             (ownerDistance+neighbourDistance+std::numeric_limits<double>::min()),
             edgeEntity});
        samples["volume_ratio"].push_back(
            {std::min(areas[edge.owner],areas[*edge.neighbour])/
             (std::max(areas[edge.owner],areas[*edge.neighbour])+
              std::numeric_limits<double>::min()),edgeEntity});
        // Divide each area by its own total background box area before comparing, so
        // the 2:1 level difference cancels instead of being scored as a defect.
        const double ownerBackground=metadata[edge.owner].backgroundArea;
        const double neighbourBackground=metadata[*edge.neighbour].backgroundArea;
        if (ownerBackground>0.0 && neighbourBackground>0.0) {
            const double ownerFraction=areas[edge.owner]/ownerBackground;
            const double neighbourFraction=areas[*edge.neighbour]/neighbourBackground;
            samples["background_volume_ratio"].push_back(
                {std::min(ownerFraction,neighbourFraction)/
                 (std::max(ownerFraction,neighbourFraction)+
                  std::numeric_limits<double>::min()),edgeEntity});
        }
    }

    const std::map<std::string,bool> lowerIsWorse{
        {"hydraulic_aspect",false},{"minimum_interior_angle_deg",true},
        {"face_length_over_local_background_h",true},
        {"face_length_over_sqrt_owner_area",true},
        {"face_length_over_sqrt_neighbour_area",true},
        {"non_orthogonality_deg",false},{"skewness",false},
        {"face_weight",true},{"volume_ratio",true},
        {"background_volume_ratio",true}};
    for (const auto& [metric,metricSamples]:samples) {
        report.ordinaryMetrics[metric]=summarize(metricSamples,lowerIsWorse.at(metric));
        for (const auto& sample:metricSamples) {
            const auto& limits=limitsFor(contract,sample.entity.cellType);
            if (metric=="hydraulic_aspect") checkSample(report,metric,sample,limits.hydraulicAspect);
            else if (metric=="minimum_interior_angle_deg") checkSample(report,metric,sample,limits.minimumInteriorAngleDeg);
            else if (metric=="face_length_over_local_background_h") checkSample(report,metric,sample,limits.faceOverLocalBackgroundH);
            else if (metric=="face_length_over_sqrt_owner_area") checkSample(report,metric,sample,limits.faceOverSqrtOwnerArea);
            else if (metric=="face_length_over_sqrt_neighbour_area") checkSample(report,metric,sample,limits.faceOverSqrtNeighbourArea);
            else if (metric=="non_orthogonality_deg") checkSample(report,metric,sample,limits.nonOrthogonalityDeg);
            else if (metric=="skewness") checkSample(report,metric,sample,limits.skewness);
            else if (metric=="face_weight") checkSample(report,metric,sample,limits.faceWeight);
            else if (metric=="volume_ratio") checkSample(report,metric,sample,limits.volumeRatio);
            // "background_volume_ratio" is computed and reported but not gated yet.
            // It is the grading-aware form and its limit is already certified
            // reachable by tests/quality_reachability_test.cpp, but on this topology
            // it cannot replace the raw gate: the contract is evaluated on the solver
            // convex-partitioned mesh, where one background box becomes several
            // triangles.  A piece inherits the whole parent lineage and therefore the
            // whole background area, so its area fraction is a partition artefact
            // rather than a sizing property.  Measured on circle at level 8 the gated
            // form reports 689 hard issues with a p50 of 1/3, against 0 hard and a p50
            // of 1/2 for the raw form.  Switching the gate requires evaluating the
            // contract on the stabilized topology instead, which is a separate change.
        }
    }

    report.boundaryLayerMetrics["wall_normal_orthogonality_error_deg"]=
        summarize(boundaryLayer.wallNormalOrthogonalityDeg,false);
    report.boundaryLayerMetrics["growth_ratio"]=summarize(boundaryLayer.growthRatio,false);
    report.boundaryLayerMetrics["tangential_normal_spacing_ratio"]=
        summarize(boundaryLayer.tangentialNormalSpacingRatio,false);
    report.boundaryLayerMetrics["adjacent_column_thickness_variation"]=
        summarize(boundaryLayer.adjacentColumnThicknessVariation,false);
    report.boundaryLayerMetrics["scaled_jacobian"]=summarize(boundaryLayer.scaledJacobian,true);
    report.boundaryLayerMetrics["first_layer_continuity"]=
        summarize(boundaryLayer.firstLayerContinuity,true);
    return report;
}

const char* qualityCellTypeName(QualityCellType2D type) noexcept {
    switch (type) {
    case QualityCellType2D::Cartesian: return "cartesian";
    case QualityCellType2D::RemainderCut: return "remainder_cut";
    case QualityCellType2D::Transition: return "transition";
    case QualityCellType2D::Termination: return "termination";
    case QualityCellType2D::BoundaryLayer: return "boundary_layer";
    case QualityCellType2D::Unknown: return "unknown";
    }
    return "unknown";
}

const char* qualityContractStatusName(QualityContractStatus2D status) noexcept {
    switch (status) {
    case QualityContractStatus2D::Pass: return "PASS";
    case QualityContractStatus2D::Warn: return "WARN";
    case QualityContractStatus2D::Fail: return "FAIL";
    }
    return "FAIL";
}

std::string qualityContractReportToJson(const QualityContractReport2D& report,
                                        int indentSpaces) {
    const int step=std::max(0,indentSpaces);
    const std::string i1=indent(step),i2=indent(2*step),i3=indent(3*step);
    std::ostringstream out;
    out<<std::setprecision(17);
    out<<"{\n"<<i1<<"\"quality_class\": \"solver_quality_contract\",\n"
       <<i1<<"\"dimensionless\": true,\n"
       <<i1<<"\"status\": \""<<qualityContractStatusName(report.status())<<"\",\n"
       <<i1<<"\"contract\": {\n";
    const auto writeType=[&](QualityCellType2D type,
                             const OrdinaryCellQualityLimits2D& limits,bool comma) {
        out<<i2<<'"'<<qualityCellTypeName(type)<<"\": ";
        writeOrdinaryLimits(out,limits,i3);
        out<<(comma?",\n":"\n");
    };
    writeType(QualityCellType2D::Cartesian,report.contract.cartesian,true);
    writeType(QualityCellType2D::RemainderCut,report.contract.remainderCut,true);
    writeType(QualityCellType2D::Transition,report.contract.transition,true);
    writeType(QualityCellType2D::Termination,report.contract.termination,false);
    out<<i1<<"},\n"<<i1<<"\"by_cell_type\": {";
    bool first=true;
    for (const auto& [type,result]:report.byCellType) {
        if (!first) out<<',';
        out<<"\n"<<i2<<'"'<<qualityCellTypeName(type)<<"\": {\"status\": \"";
        if (result.rated) out<<qualityContractStatusName(result.status());
        else out<<"OBSERVED";
        out<<"\", \"rated\": "<<(result.rated?"true":"false")
           <<", \"cell_count\": "
           <<result.cellCount<<", \"preferred_issue_count\": "
           <<result.preferredIssueCount<<", \"hard_issue_count\": "
           <<result.hardIssueCount<<'}';
        first=false;
    }
    if (!report.byCellType.empty()) out<<'\n'<<i1;
    out<<"},\n";
    const auto writeSummaries=[&](const char* name,
                                  const std::map<std::string,QualityMetricSummary2D>& metrics,
                                  bool comma) {
        out<<i1<<'"'<<name<<"\": {";
        bool firstMetric=true;
        for (const auto& [metric,summary]:metrics) {
            if (!firstMetric) out<<',';
            out<<"\n"<<i2<<'"'<<metric<<"\": {\"count\": "<<summary.count
               <<", \"p50\": "<<summary.p50<<", \"p95\": "<<summary.p95
               <<", \"p99\": "<<summary.p99<<", \"worst\": "<<summary.worst
               <<", \"worst_direction\": \""
               <<(summary.lowerIsWorse?"min":"max")<<"\", \"worst_entity\": ";
            if (summary.worstEntity) writeEntity(out,*summary.worstEntity,i3);
            else out<<"null";
            out<<'}';
            firstMetric=false;
        }
        if (!metrics.empty()) out<<'\n'<<i1;
        out<<'}'<<(comma?",\n":"\n");
    };
    writeSummaries("ordinary_metrics",report.ordinaryMetrics,true);
    writeSummaries("boundary_layer_metrics",report.boundaryLayerMetrics,true);
    out<<i1<<"\"legacy_hard_safety\": {\"valid\": "
       <<(report.legacyHardSafety.valid()?"true":"false")
       <<", \"issue_count\": "<<report.legacyHardSafety.issues.size()
       <<", \"min_face_length_absolute\": "<<report.legacyHardSafety.minFaceLength
       <<", \"max_non_orthogonality_deg\": "
       <<report.legacyHardSafety.maxNonOrthogonalityDeg
       <<", \"max_internal_skewness\": "<<report.legacyHardSafety.maxInternalSkewness
       <<", \"max_boundary_skewness\": "<<report.legacyHardSafety.maxBoundarySkewness
       <<", \"max_cell_aspect\": "<<report.legacyHardSafety.maxCellAspect
       <<", \"min_face_weight\": "<<report.legacyHardSafety.minFaceWeight
       <<", \"min_volume_ratio\": "<<report.legacyHardSafety.minVolumeRatio<<"},\n";
    out<<i1<<"\"issues\": [";
    for (std::size_t index=0;index<report.issues.size();++index) {
        const auto& issue=report.issues[index];
        out<<(index==0?"\n":",\n")<<i2<<"{\"level\": \""
           <<(issue.level==QualityContractLevel2D::Hard?"hard":"preferred")
           <<"\", \"metric\": \""<<issue.metric<<"\", \"measured\": "
           <<issue.measured<<", \"limit\": "<<issue.limit<<", \"entity\": ";
        writeEntity(out,issue.entity,i3);
        out<<'}';
    }
    if (!report.issues.empty()) out<<'\n'<<i1;
    out<<"],\n"<<i1<<"\"input_issues\": [";
    for (std::size_t index=0;index<report.inputIssues.size();++index) {
        if (index!=0U) out<<", ";
        out<<'"'<<report.inputIssues[index]<<'"';
    }
    out<<"]\n}\n";
    return out.str();
}

} // namespace cartmesh2d
