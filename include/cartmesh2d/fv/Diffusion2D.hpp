#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <functional>

namespace cartmesh2d::fv {
struct DiffusionProblem2D {
    double diffusivity = 1;
    // -div(diffusivity grad(value)) = source. All boundary faces Dirichlet.
    std::function<double(Point2D)> source;
    std::function<double(Point2D, BoundaryPatch2D)> boundaryValue;
};
struct DiffusionControls2D {
    std::size_t maxCorrections = 400;
    std::size_t maxLinearIterations = 5000;
    double relativeTolerance = 1e-10;
    double absoluteTolerance = 1e-12;
    double linearRelativeTolerance = 1e-13;
    double linearAbsoluteTolerance = 1e-14;
    double relaxation = 0.7;
};
struct DiffusionIteration2D {
    std::size_t iteration = 0, linearIterations = 0;
    double relativeResidual = 0, residualNorm = 0, maxCellImbalance = 0;
};
struct DiffusionResult2D {
    bool converged = false;
    std::vector<double> values, fluxes, sourceIntegrals;
    std::vector<DiffusionIteration2D> history;
    double boundaryFlux = 0, sourceIntegral = 0, globalBalance = 0;
};
// Linear SPD two-point part solved by Jacobi-PCG, with deferred non-orthogonal
// correction from weighted least squares gradients. Failure never means PASS.
[[nodiscard]] DiffusionResult2D solveDiffusion2D(
    const FvMesh2D&, const DiffusionProblem2D&, const DiffusionControls2D& = {});
}
