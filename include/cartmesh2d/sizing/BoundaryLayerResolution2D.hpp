#pragma once
#include "cartmesh2d/topology/Topology2D.hpp"
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {
struct HybridMeshBuildResult2D;
struct LayerColumnResolution2D {
    std::size_t stripId=0, wallSegment=0, requestedLayers=0, retainedLayers=0;
    double wallLengthOverReference=0;
    Segment2D wall;
    std::vector<std::size_t> solverCellIds;
    std::optional<double> firstLayerNormalHeightMin,firstLayerNormalHeightMax;
};
struct BoundaryLayerResolution2D {
    std::string status="not_requested";
    std::size_t requestedCells=0,constructedCells=0,retainedCells=0,sourceMismatches=0;
    std::size_t continuityMismatches=0;
    double wallLength=0,firstLayerWallLength=0,fullLayerWallLength=0;
    double firstLayerHeightExceededWallLength=0;
    std::optional<double> requestedFirstLayerHeight;
    std::vector<LayerColumnResolution2D> columns;
    std::vector<std::string> issues;
};
// Layer identity alone is insufficient: remeasure final polygons, embedded wall
// ownership and the connected sequence of final cells in every column.
[[nodiscard]] BoundaryLayerResolution2D measureBoundaryLayerResolution2D(
    const TopologyMesh2D& mesh,double referenceLength,std::optional<std::size_t> requestedLayers,
    const HybridMeshBuildResult2D* hybrid=nullptr,const TolerancePolicy& tol={},
    std::optional<double> requestedFirstLayerHeight=std::nullopt);
[[nodiscard]] std::string boundaryLayerResolutionToJson2D(const BoundaryLayerResolution2D& report);
}
