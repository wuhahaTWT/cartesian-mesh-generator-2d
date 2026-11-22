#pragma once
#include "FlowSolverDetail2D.hpp"

namespace cartmesh2d::fv::solver_detail {
// Called after the common controls, geometry, connectivity and guess checks.
FlowResult2D solveCoupledFlow2D(const FvMesh2D&,const FlowControls2D&,const Boundary&,
    const std::function<void(const FlowIteration2D&)>&,const FlowInitialGuess2D*);
}
