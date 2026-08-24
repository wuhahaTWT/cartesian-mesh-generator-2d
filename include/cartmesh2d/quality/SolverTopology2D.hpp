#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cartmesh2d {

struct SolverTopologyResult2D {
    TopologyMesh2D topology;
    std::size_t inputCellCount = 0;
    std::size_t outputCellCount = 0;
    std::size_t partitionedCellCount = 0;
    std::size_t qualityAgglomeratedSourceCellCount = 0;
    std::size_t qualityRepartitionCount = 0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid() && outputCellCount>0;
    }
};

struct SolverLocalRepartitionResult2D {
    TopologyMesh2D topology;
    std::size_t repartitionCount = 0;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid();
    }
};

// Replaces two adjacent solver cells by an area-identical alternative convex
// two-piece partition only when the complete solver-quality score improves.
[[nodiscard]] SolverLocalRepartitionResult2D repartitionSolverTopologyByQuality2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

// Retains strictly convex cells and deterministically ear-clips cells with
// reflex or collinear vertices.  The latter occur at conformal 2:1 Cartesian
// transitions and are reported as concave polyhedra by OpenFOAM when extruded.
[[nodiscard]] SolverTopologyResult2D buildSolverTopology2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
