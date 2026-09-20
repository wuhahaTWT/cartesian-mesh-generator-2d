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
    // Optional resolved-wall mask. Such faces require k=0, finite positive
    // omega, exactly zero carrier flux; face turbulent diffusivity is zero.
    std::vector<bool> resolvedWalls;
};
struct FrozenSst2003mResult2D {
    ScalarTransportResult2D k, omega;
    std::vector<Sst2003mCoefficients2D> coefficients;
};
struct SstVelocityBoundary2D {
    Vector2D value{};
    // Each unset component has zero normal derivative, not zero value.
    bool fixedX = false, fixedY = false;
};
struct Sst2003mGradients2D {
    std::vector<Vector2D> k, omega, u, v;
    std::vector<double> strainMagnitude;
};
// Least-squares reconstruction from actual cell and face data. Turbulence BCs
// may be fixed values or ZERO diffusive flux only (nonzero flux requires known
// effective diffusivity and is explicitly rejected here). Velocity components
// independently accept fixed values or zero-normal derivatives.
[[nodiscard]] Sst2003mGradients2D reconstructSst2003mGradients2D(
    const FvMesh2D&, const FrozenSst2003mProblem2D&,
    const std::vector<Vector2D>& velocity,
    const std::vector<SstVelocityBoundary2D>& velocityBoundary);
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

// Set k_wall=0 and omega_wall=60*nu/(beta1*d_normal^2), beta1=.075,
// with each wall face's own owner-to-face normal spacing. Low-Re resolved-wall
// prescription from the TMR SST reference, NOT a y+ wall function or a mesh
// qualification. Does not overwrite inlet/outlet data or cell distance field.
void setSst2003mResolvedWalls2D(const FvMesh2D&, FrozenSst2003mProblem2D&,
    const std::vector<bool>& walls);
struct SstTransportControls2D {
    std::size_t maxIterations = 500;
    double relaxation = .5;
    ScalarTransportControls2D transport;
    // Experimental steady Picard/deferred-correction iteration. Zero preserves
    // complete frozen transport solves. Positive values start from current k/w
    // and take at most this many scalar corrections before updating the closure.
    // A correction is NOT a converged transport result; only the recomputed full
    // nonlinear equations with transport's original gates can accept the result.
    std::size_t scalarCorrectionsPerUpdate = 0;
};
struct SstTransportIteration2D {
    std::size_t iteration = 0;
    double kResidualNorm = 0, omegaResidualNorm = 0;
    double kCellResidual = 0, omegaCellResidual = 0;
};
struct SstTransportResult2D {
    ScalarTransportPerformance2D scalarSolves, scalarEvaluations;
    bool converged = false;
    // Coefficients/fluxes/residuals evaluated at the RETURNED fields, never
    // reused from the preceding frozen linear solve.
    FrozenSst2003mResult2D fields;
    std::vector<SstTransportIteration2D> history;
};
// Nonlinear SST transport on a FIXED, conservative velocity/flux field. Each
// iteration reconstructs k/omega gradients and updates the closure; acceptance
// evaluates the original nonlinear equations with current coefficients.
// No momentum/pressure feedback, wall functions, or RANS convergence claim.
[[nodiscard]] SstTransportResult2D solveSst2003mTransport2D(
    const FvMesh2D&, const FrozenSst2003mProblem2D& initial,
    const std::vector<Vector2D>& velocity,
    const std::vector<SstVelocityBoundary2D>& velocityBoundary,
    const SstTransportControls2D& = {},
    const std::vector<double>& previousK = {},
    const std::vector<double>& previousOmega = {}, double timeStep = 0);
}
