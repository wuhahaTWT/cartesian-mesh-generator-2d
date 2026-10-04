#pragma once
#include "cartmesh2d/fv/Incompressible2D.hpp"

namespace cartmesh2d::fv::detail {
FlowResult2D solveNewtonKrylovFlow2D(const FvMesh2D&,const FlowControls2D&,
    const FlowInitialGuess2D*,const std::function<void(const FlowIteration2D&)>&);
}
