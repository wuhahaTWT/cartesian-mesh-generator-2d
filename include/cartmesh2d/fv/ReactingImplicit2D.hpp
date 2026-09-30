#pragma once

#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include <functional>
#include <array>
#include <memory>

namespace cartmesh2d::fv {

struct ReactingImplicitControls2D {
    // CVODES weighted local errors. Conserved variables use fixed initial
    // rho, rho*a, rho*cv*T scales; species use initial rho times mass fraction.
    // These are integration controls, not flame/experimental accuracy gates.
    double relativeTolerance = 1e-7;
    double absoluteConservedTolerance = 1e-12, absoluteSpeciesFraction = 1e-18;
    double maximumStep = 1e-5, initialStep = 0;
    unsigned spatialOrder = 2, maximumBdfOrder = 2;
    unsigned maximumNonlinearIterations = 3;
    // Experimental globalization option; the default rejects an unconverged
    // damped update and lets CVODES retry the time step.
    bool continueDampedNewton = false;
    // Try a species-bound reflection of an infeasible Newton candidate,
    // requiring a decrease of the true nonlinear residual before use.
    bool reflectSpeciesNewton = false;
    // Approximate Newton matrix only. A first-order advection linearization
    // avoids differencing discontinuous MUSCL limiter switches. The actual
    // residual and its convergence checks still use spatialOrder throughout.
    unsigned jacobianAdvectionOrder = 1;
    std::size_t maximumAcceptedSteps = 20000;
    std::size_t maximumResidualEvaluations = 100000;
    std::size_t maximumBandBytes = 256 * 1024 * 1024;
};

struct ReactingImplicitProgress2D {
    bool reachedEnd = false, canceled = false;
    std::string failure;
    ReactingState2D lastAccepted; // initial when no accepted step exists
    // Since session initialization, including all open boundaries and source
    // terms. These are separately integrated CVODES quadratures.
    ReactingConservative2D boundaryImpulse, chemistryChange, constraintChange;
    std::size_t rhsCalls = 0, rejectedRhsCalls = 0, bandHalfWidth = 0, bandBytes = 0;
    long internalSteps = 0, errorTestFailures = 0, linearSetups = 0;
    long jacobianEvaluations = 0, nonlinearIterations = 0, nonlinearConvergenceFailures = 0;
    std::array<long, 6> acceptedByBdfOrder{};
    int lastBdfOrder = 0;
    std::size_t dampedNewtonUpdates = 0, continuedDampedNewtonUpdates = 0;
    std::size_t reflectionAttempts = 0, reflectedNewtonUpdates = 0;
    double maximumReflectedSpeciesScaledChange = 0;
    double minimumNewtonFraction = 1;
    double lastNewtonCorrectionNorm = 0, lastNewtonResidualNorm = 0, lastNewtonTolerance = 0;
    double lastPredictorCorrectionNorm = 0, lastLocalErrorNorm = 0, maximumLocalErrorNorm = 0;
    std::string lastDampedTrialFailure;
    double minimumAcceptedStep = 0, maximumAcceptedStep = 0;
};

// Fully coupled BDF integration of the same native compressible finite-volume
// residual, including detailed reactions and molecular transport. No splitting
// or physical-model reduction. N-1 independent species and total density define
// the redundant species algebraically. RCM ordering changes only linear algebra.
// A banded direct solve is currently used; large bandwidth is explicitly rejected.
// The gas and residual contexts must outlive this serial session.
class ReactingImplicitIntegrator2D {
public:
    ReactingImplicitIntegrator2D(chemistry::DetailedGas&, ReactingFlowStepper2D&,
                                 const ReactingState2D&, const ReactingImplicitControls2D& = {});
    ~ReactingImplicitIntegrator2D();
    ReactingImplicitIntegrator2D(const ReactingImplicitIntegrator2D&) = delete;
    ReactingImplicitIntegrator2D& operator=(const ReactingImplicitIntegrator2D&) = delete;
    // Cancel between steps preserves this session. Cancel during an internal
    // trial, or a failure, retains the last accepted state and invalidates the
    // BDF history; restart from that state with a fresh session (not bitwise).
    [[nodiscard]] ReactingImplicitProgress2D advance(double endTime,
        const std::function<bool()>& cancel = {},
        const std::function<void(const ReactingImplicitProgress2D&)>& onAccepted = {},
        // Trial-evaluation observer, including startup before the first
        // accepted step. The reported state remains the last accepted state.
        const std::function<void(const ReactingImplicitProgress2D&, double)>& onEvaluation = {});
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cartmesh2d::fv
