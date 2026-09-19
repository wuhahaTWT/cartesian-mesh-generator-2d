#pragma once
#include "cartmesh2d/fv/ScalarTransport2D.hpp"

namespace cartmesh2d::fv {
// Fixed-in-time data are explicit so restart checks compare the actual fields,
// not a user-supplied label or unverifiable callback identity.
struct ThermalSetup2D {
    double diffusivity = .01;
    std::vector<double> sourceDensity;
    std::vector<ScalarBoundary2D> boundary;
};
struct ThermalFlowState2D {
    FlowState2D flow; // the single physical clock for both fields
    std::vector<double> scalar;
};
struct ThermalFlowResult2D {
    FlowResult2D flow;
    ScalarTransportResult2D scalar;
    // Present only after BOTH equations converge at the same new physical time.
    std::optional<ThermalFlowState2D> accepted;
};
void validateThermalSetup2D(const FvMesh2D&, const ThermalSetup2D&);
// One-way constant-property coupling: flow at t+dt, then conservative scalar
// using precisely that flow's face flux. No buoyancy or scalar feedback.
// Input is immutable on nonconvergence, exceptions, or cancellation callback.
[[nodiscard]] ThermalFlowResult2D advanceThermalFlow2D(
    const FvMesh2D&, const FlowControls2D&, const ThermalSetup2D&,
    const ScalarTransportControls2D&, const ThermalFlowState2D&, double timeStep,
    const std::function<void(const FlowIteration2D&)>& progress = {});
}
