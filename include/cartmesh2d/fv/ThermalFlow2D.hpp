#pragma once
#include "cartmesh2d/fv/ScalarTransport2D.hpp"
#include "cartmesh2d/fv/FlowTimeStep2D.hpp"

namespace cartmesh2d::fv {
// Fixed-in-time data are explicit so restart checks compare the actual fields,
// not a user-supplied label or unverifiable callback identity.
struct ThermalEvent2D {
    double time = 0; // right-continuous piecewise-constant data
    std::vector<double> sourceDensity;
    std::vector<ScalarBoundary2D> boundary;
};
struct ThermalSetup2D {
    double diffusivity = .01;
    std::vector<double> sourceDensity;
    std::vector<ScalarBoundary2D> boundary;
    std::vector<ThermalEvent2D> events; // complete snapshots; bound to restart
};
struct ThermalControllerHistory2D {
    std::vector<double> controls; // exact numerical signature; reset on change
    double nextStep=0, velocityRelaxation=0;
};
struct ThermalFlowState2D {
    FlowState2D flow; // the single physical clock for both fields
    std::vector<double> scalar;
    std::optional<ThermalControllerHistory2D> controller;
};
struct ThermalFlowResult2D {
    FlowResult2D flow;
    ScalarTransportResult2D scalar;
    // Present only after BOTH equations converge at the same new physical time.
    std::optional<ThermalFlowState2D> accepted;
    // False only when a converged flow still fails the unchanged scalar carrier
    // gate after a tighter pressure solve. The controller must reject this trial.
    bool carrierCompatible = true;
};
struct ThermalTimeControls2D {
    FlowTimeStepControls2D limits;
    bool estimateError = false;
    double temperatureScale = 1; // K: characteristic temperature CHANGE, not absolute T
    double velocityScale = 1; // m/s
    double relativeTolerance = .01; // error / characteristic scale
    double temperatureAbsoluteTolerance = .001; // K
    double velocityAbsoluteTolerance = .0001; // m/s
};
struct ThermalAttempt2D {
    double startTime=0, timeStep=0, error=0, courant=0, velocityRelaxation=0;
    std::string reason;
};
struct ThermalControlledResult2D {
    ThermalFlowResult2D step;
    std::vector<ThermalAttempt2D> attempts;
};
// Accepted recommendation/history is explicit and checkpointed. A change to
// numerical controls resets it; physical compatibility remains strict. BE full step is
// accepted; two half steps estimate its error, without Richardson extrapolation.
[[nodiscard]] ThermalControlledResult2D advanceControlledThermalFlow2D(
    const FvMesh2D&, const FlowControls2D&, const ThermalSetup2D&,
    const ScalarTransportControls2D&, const ThermalFlowState2D&, const ThermalTimeControls2D&,
    const std::function<bool()>& cancelled = {},
    const std::function<void(const ThermalAttempt2D&)>& observeAttempt = {});
void validateThermalSetup2D(const FvMesh2D&, const ThermalSetup2D&);
// One-way constant-property coupling: flow at t+dt, then conservative scalar
// using precisely that flow's face flux. No buoyancy or scalar feedback.
// Input is immutable on nonconvergence, exceptions, or cancellation callback.
[[nodiscard]] ThermalFlowResult2D advanceThermalFlow2D(
    const FvMesh2D&, const FlowControls2D&, const ThermalSetup2D&,
    const ScalarTransportControls2D&, const ThermalFlowState2D&, double timeStep,
    const std::function<void(const FlowIteration2D&)>& progress = {});
}
