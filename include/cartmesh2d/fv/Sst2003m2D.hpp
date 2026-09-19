#pragma once

#include "cartmesh2d/fv/ScalarTransport2D.hpp"

namespace cartmesh2d::fv {
// Incompressible, density-normalized SST-2003m. Definition:
// https://tmbwg.github.io/turbmodels/sst.html#sst-2003
// Fixed SI units with rho=1 kg/m^3 for the published CD floor 1e-10.
// This is a local closure/frozen transport building block, not a RANS solver.
struct Sst2003mPoint2D {
    double k = 0, omega = 1, nu = 1e-5, wallDistance = 1;
    // sqrt(2*Sij*Sij), including both xy and yx in the tensor contraction.
    double strainMagnitude = 0;
    Vector2D gradientK{}, gradientOmega{};
};
struct Sst2003mCoefficients2D {
    double f1 = 0, f2 = 0, turbulentViscosity = 0;
    double sigmaK = 0, sigmaOmega = 0, beta = 0, gamma = 0;
    double diffusivityK = 0, diffusivityOmega = 0;
    // source - lossRate * newScalar. At the supplied k/omega, these reproduce
    // the original nonlinear source exactly. Negative cross diffusion is a loss.
    double productionK = 0, productionOmega = 0, crossDiffusion = 0;
    double sourceK = 0, sourceOmega = 0, lossRateK = 0, lossRateOmega = 0;
};
[[nodiscard]] Sst2003mCoefficients2D evaluateSst2003m2D(const Sst2003mPoint2D&);

struct FrozenSst2003mProblem2D {
    double nu = 1e-5;
    std::vector<double> k, omega, wallDistance, strainMagnitude, volumeFlux;
    std::vector<Vector2D> gradientK, gradientOmega;
    // Caller supplies actual wall/inlet/outlet data, not a guessed wall function.
    std::vector<ScalarBoundary2D> boundaryK, boundaryOmega;
};
struct FrozenSst2003mResult2D {
    ScalarTransportResult2D k, omega;
    std::vector<Sst2003mCoefficients2D> coefficients;
};
// One frozen nonlinear iteration: evaluate closure at p.k/omega, linearly
// interpolate cell diffusivities to internal face centres with mesh weights,
// use owner diffusivity at boundary faces, solve k and omega with implicit loss.
// Caller supplies consistent gradients, nearest-wall distances and carrier flux.
// Empty previous fields mean steady; otherwise BOTH previous states and dt>0.
// Positive results/inner convergence do NOT establish nonlinear RANS convergence.
// No clipping: failed linear solves or negative k/nonpositive omega throw.
[[nodiscard]] FrozenSst2003mResult2D solveFrozenSst2003mTransport2D(
    const FvMesh2D&, const FrozenSst2003mProblem2D&,
    const ScalarTransportControls2D& = {},
    const std::vector<double>& previousK = {},
    const std::vector<double>& previousOmega = {}, double timeStep = 0);
}
