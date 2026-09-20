#pragma once
#include "cartmesh2d/topology/Topology2D.hpp"

namespace cartmesh2d {
// Native 2D, filled rectangular fluid region. Explicit strictly increasing
// coordinate arrays allow independent streamwise and wall-normal grading.
// This is not a cut-cell replacement for solid geometry. All outer faces are
// DomainBoundary; callers apply physical boundary conditions separately.
// Throws on invalid input or failure of the unchanged default Solver quality gate.
[[nodiscard]] TopologyMesh2D makeRectilinearMesh2D(
    const std::vector<double>& x, const std::vector<double>& y);
}
