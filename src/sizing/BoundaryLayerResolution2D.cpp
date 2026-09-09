#include "cartmesh2d/sizing/BoundaryLayerResolution2D.hpp"
#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace cartmesh2d {
namespace {
bool samePolygon(const Polygon2D& source,const TopologyCell2D& cell,
                 const TopologyMesh2D& mesh,double reference,const TolerancePolicy& tol) {
    if (source.vertices.size()<3 || cell.vertices.size()<3) return false;
    const auto origin=source.vertices.front();
    const auto normalized=[&](Point2D p) { return Point2D{(p.x-origin.x)/reference,(p.y-origin.y)/reference}; };
    Polygon2D a,b;
    for (const auto p:source.vertices) a.vertices.push_back(normalized(p));
    for (const auto id:cell.vertices) b.vertices.push_back(normalized(mesh.vertices.at(id).point));
    if (!tol.nearlyEqual(a.area(),b.area())) return false;
    for (const auto p:b.vertices) {
        bool on=false;
        for (std::size_t i=0;i<a.vertices.size();++i)
            on=on || pointOnSegment(p,{a.vertices[i],a.vertices[(i+1)%a.vertices.size()]},tol);
        if (!on) return false;
    }
    for (const auto p:a.vertices) {
        bool present=false;
        for (const auto q:b.vertices) present=present || (tol.nearlyEqual(p.x,q.x) && tol.nearlyEqual(p.y,q.y));
        if (!present) return false;
    }
    return true;
}
}

BoundaryLayerResolution2D measureBoundaryLayerResolution2D(
    const TopologyMesh2D& mesh,double reference,std::optional<std::size_t> requested,
    const HybridMeshBuildResult2D* hybrid,const TolerancePolicy& tol,
    std::optional<double> requestedHeight) {
    if (!(reference>0) || !std::isfinite(reference) || (requested && *requested==0) ||
        (requestedHeight && (!std::isfinite(*requestedHeight) || !(*requestedHeight>0))))
        throw std::invalid_argument("invalid boundary layer resolution request");
    BoundaryLayerResolution2D result;
    if (requestedHeight) result.requestedFirstLayerHeight=*requestedHeight/reference;
    for (const auto& edge:mesh.edges) if (!edge.neighbour && edge.patch==BoundaryPatch2D::EmbeddedBoundary) {
        const auto p=mesh.vertices.at(edge.v0).point,q=mesh.vertices.at(edge.v1).point;
        result.wallLength+=std::hypot(p.x-q.x,p.y-q.y)/reference;
    }
    if (!hybrid) { result.status=requested?"no_layers_retained":"not_requested";return result; }
    result.status="evaluated";
    using Key=std::tuple<std::size_t,std::size_t,std::size_t>;
    std::map<Key,std::size_t> sourceAt;
    std::vector<std::vector<std::size_t>> finalForSource(hybrid->sourceCells.size());
    for (const auto& cell:mesh.cells) {
        auto lineage=cell.sourceLineage;
        if (lineage.empty()) lineage.push_back(cell.sourceId);
        std::sort(lineage.begin(),lineage.end());
        lineage.erase(std::unique(lineage.begin(),lineage.end()),lineage.end());
        if (lineage.size()==1 && lineage.front()<finalForSource.size())
            finalForSource[lineage.front()].push_back(cell.id);
    }
    std::map<std::size_t,Key> certifiedCells;
    for (const auto& source:hybrid->sourceCells) {
        if (source.kind!=HybridCellKind2D::BoundaryLayer) continue;
        ++result.constructedCells;
        if (!source.stripId || !source.wallSegment || !source.layerIndex ||
            *source.stripId>=hybrid->strips.size() ||
            *source.wallSegment>=hybrid->strips[*source.stripId].wallChain.segmentCount() ||
            source.id>=finalForSource.size()) {++result.sourceMismatches;continue;}
        const auto& finals=finalForSource[source.id];
        const Key key{*source.stripId,*source.wallSegment,*source.layerIndex};
        if (sourceAt.contains(key) || finals.size()!=1 ||
            !samePolygon(source.polygon,mesh.cells.at(finals.front()),mesh,reference,tol)) {
            ++result.sourceMismatches;continue;
        }
        sourceAt[key]=source.id;
        certifiedCells[finals.front()]=key;
        ++result.retainedCells;
    }
    std::set<std::pair<std::size_t,std::size_t>> neighbours;
    std::set<std::size_t> wallOwners;
    for (const auto& edge:mesh.edges) {
        if (edge.neighbour) neighbours.emplace(std::min(edge.owner,*edge.neighbour),std::max(edge.owner,*edge.neighbour));
        else if (edge.patch==BoundaryPatch2D::EmbeddedBoundary) wallOwners.insert(edge.owner);
    }
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> columnAt;
    for (std::size_t stripId=0;stripId<hybrid->strips.size();++stripId) {
        const auto& strip=hybrid->strips[stripId];
        const auto layers=requested.value_or(strip.parameters.nLayers);
        for (std::size_t segment=0;segment<strip.wallChain.segmentCount();++segment) {
            LayerColumnResolution2D column;
            column.stripId=stripId;column.wallSegment=segment;column.requestedLayers=layers;
            const auto wall=strip.wallChain.segments[segment];
            column.wall=wall;
            const double length=std::hypot(wall.b.x-wall.a.x,wall.b.y-wall.a.y);
            column.wallLengthOverReference=length/reference;
            result.requestedCells+=layers;
            for (std::size_t layer=0;layer<layers;++layer) {
                const auto found=sourceAt.find({stripId,segment,layer});
                if (found==sourceAt.end()) break;
                const auto cell=finalForSource[found->second].front();
                if (layer==0 && !wallOwners.contains(cell)) {++result.continuityMismatches;break;}
                if (layer>0) {
                    const auto previous=column.solverCellIds.back();
                    if (!neighbours.contains({std::min(previous,cell),std::max(previous,cell)})) {
                        ++result.continuityMismatches;break;
                    }
                }
                column.solverCellIds.push_back(cell);
                ++column.retainedLayers;
                if (layer==0 && length>0) {
                    const double sign=strip.wallChain.fluidSide==FluidSide2D::Left?1.:-1.;
                    for (const auto sample:hybrid->sourceCells[found->second].polygon.vertices) {
                        // Source corner presence was checked above; measure at
                        // its actual final vertex, excluding added collinear
                        // side nodes from the outer-front height statistic.
                        Point2D p=sample;double distance=std::numeric_limits<double>::infinity();
                        for (const auto vertex:mesh.cells.at(cell).vertices) {
                            const auto q=mesh.vertices.at(vertex).point;
                            const double candidate=std::hypot(q.x-sample.x,q.y-sample.y);
                            if (candidate<distance) {distance=candidate;p=q;}
                        }
                        const double normal=sign*((wall.b.x-wall.a.x)*(p.y-wall.a.y)-
                            (wall.b.y-wall.a.y)*(p.x-wall.a.x))/(length*reference);
                        if (normal>tol.scale(1.)) {
                            column.firstLayerNormalHeightMin=std::min(column.firstLayerNormalHeightMin.value_or(normal),normal);
                            column.firstLayerNormalHeightMax=std::max(column.firstLayerNormalHeightMax.value_or(normal),normal);
                        }
                    }
                }
            }
            columnAt[{stripId,segment}]=result.columns.size();
            result.columns.push_back(std::move(column));
        }
    }
    for (const auto& edge:mesh.edges) {
        if (edge.neighbour || edge.patch!=BoundaryPatch2D::EmbeddedBoundary) continue;
        const auto owner=certifiedCells.find(edge.owner);
        if (owner==certifiedCells.end()) continue;
        const auto [strip,segment,layer]=owner->second;
        if (layer!=0) continue;
        const auto& column=result.columns.at(columnAt.at({strip,segment}));
        if (!column.retainedLayers) continue;
        const auto p=mesh.vertices.at(edge.v0).point,q=mesh.vertices.at(edge.v1).point;
        const double length=std::hypot(q.x-p.x,q.y-p.y)/reference;
        result.firstLayerWallLength+=length;
        if (column.retainedLayers==column.requestedLayers) result.fullLayerWallLength+=length;
        if (result.requestedFirstLayerHeight && column.firstLayerNormalHeightMax &&
            *column.firstLayerNormalHeightMax>*result.requestedFirstLayerHeight &&
            !tol.nearlyEqual(*column.firstLayerNormalHeightMax,*result.requestedFirstLayerHeight))
            result.firstLayerHeightExceededWallLength+=length;
    }
    if (result.sourceMismatches) {
        result.status="incomplete";
        result.issues.push_back("Constructed layer cells could not be matched one-to-one to unchanged final solver polygons.");
    }
    if (result.continuityMismatches) {
        result.status="incomplete";
        result.issues.push_back("Retained layer cells lack final wall ownership or consecutive owner/neighbour connections.");
    }
    return result;
}

std::string boundaryLayerResolutionToJson2D(const BoundaryLayerResolution2D& r) {
    std::ostringstream out;out<<std::setprecision(17);
    const auto number=[&](std::optional<double> v) {if (v) out<<*v;else out<<"null";};
    out<<"{\"status\":\""<<r.status<<"\",\"requested_cells\":";
    if (r.status=="no_layers_retained") out<<"null";else out<<r.requestedCells;
    out<<",\"constructed_cells\":"<<r.constructedCells<<",\"retained_cells\":"<<r.retainedCells
       <<",\"source_mismatches\":"<<r.sourceMismatches<<",\"first_layer_wall_length_fraction\":";
    number(r.wallLength>0?std::optional<double>(r.firstLayerWallLength/r.wallLength):std::nullopt);
    out<<",\"full_requested_layers_wall_length_fraction\":";
    number(r.wallLength>0?std::optional<double>(r.fullLayerWallLength/r.wallLength):std::nullopt);
    out<<",\"first_layer_height_exceedance_wall_length_fraction\":";
    number(r.requestedFirstLayerHeight && r.wallLength>0
        ?std::optional<double>(r.firstLayerHeightExceededWallLength/r.wallLength):std::nullopt);
    out<<",\"continuity_mismatches\":"<<r.continuityMismatches<<",\"columns\":[";
    for (std::size_t i=0;i<r.columns.size();++i) {
        const auto& c=r.columns[i];out<<(i?",":"")<<"{\"strip_id\":"<<c.stripId
            <<",\"wall_segment\":"<<c.wallSegment<<",\"requested_layers\":"<<c.requestedLayers
            <<",\"retained_layers\":"<<c.retainedLayers<<",\"wall_length_over_reference\":"<<c.wallLengthOverReference
            <<",\"wall_segment_endpoints\":[["<<c.wall.a.x<<','<<c.wall.a.y<<"],["<<c.wall.b.x<<','<<c.wall.b.y
            <<"]],\"solver_cell_ids\":[";
        for (std::size_t j=0;j<c.solverCellIds.size();++j) out<<(j?",":"")<<c.solverCellIds[j];
        out<<"],\"first_layer_normal_height_min_over_reference\":";number(c.firstLayerNormalHeightMin);
        out<<",\"first_layer_normal_height_max_over_reference\":";number(c.firstLayerNormalHeightMax);out<<'}';
    }
    out<<"],\"issues\":[";
    for (std::size_t i=0;i<r.issues.size();++i) out<<(i?",":"")<<'"'<<r.issues[i]<<'"';
    out<<"]}";return out.str();
}
}
