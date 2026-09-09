#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"
#include <optional>
#include <string>

namespace cartmesh2d {

// Measurements are made on the final solver partition. Area-equivalent size is
// a descriptive statistic, never a substitute for boundary-layer directions.
struct MeshResolutionTargets2D {
    double referenceLength = 1.0;
    bool explicitReferenceLength = false;
    std::optional<double> wallSize;
    std::optional<double> backgroundSize;
    std::optional<double> firstLayerHeight;
};

[[nodiscard]] std::string meshResolutionReportToJson2D(
    const TopologyMesh2D& mesh, const MeshResolutionTargets2D& targets,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
