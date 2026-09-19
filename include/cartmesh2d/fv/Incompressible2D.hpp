#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"

#include <functional>
#include <string>

namespace cartmesh2d::fv {

enum class ConvectionScheme2D { Upwind, LimitedLinearUpwind };

enum class PressurePreconditioner2D { Jacobi, IncompleteCholesky0, Aggregation };
enum class ViscousStress2D { Laplacian, Symmetric };
enum class OutletBackflow2D { Reject, NormalInlet };

struct FlowControls2D {
    std::string scenario = "external";
    double nu = .01;
    double speed = 1;
    double tolerance = 1e-6;
    std::size_t maxIterations = 1500;
    double velocityRelaxation = .6;
    double pressureRelaxation = .25;
    bool profile = false;
    double manufacturedPressureSlope = 0; // verification-only linear pressure addition
    ViscousStress2D viscousStress = ViscousStress2D::Symmetric;
    ConvectionScheme2D convection = ConvectionScheme2D::Upwind;
    PressurePreconditioner2D pressurePreconditioner = PressurePreconditioner2D::IncompleteCholesky0;
    // Opt-in pressure outlet: reverse flow has zero tangential velocity and
    // zero normal velocity gradient. The conservative flux is never clipped.
    OutletBackflow2D outletBackflow = OutletBackflow2D::Reject;
    // Prescribed kinematic viscosity on every shared face, including boundaries.
    // Empty uses nu. Positive finite values; no implicit material averaging.
    std::vector<double> faceViscosity;
    double manufacturedViscositySlope = 0; // verification only: nu(x)=nu*(1+slope*x)
};

// Accepted state at a physical time, including the conservative face flux.
struct FlowState2D {
    double time = 0;
    std::vector<double> u, v, p, flux;
};

struct FlowPerformance2D {
    std::size_t momentumSolves = 0;
    std::size_t momentumIterations = 0;
    std::size_t maxMomentumIterations = 0;
    std::size_t pressureSolves = 0;
    std::size_t pressureCorrectionPassesSkipped = 0;
    std::size_t pressureFactorizations = 0;
    std::size_t pressureFactorReuses = 0;
    std::size_t pressureHierarchyBuilds = 0;
    std::size_t pressureHierarchyReuses = 0;
    std::size_t pressureHierarchyRefreshes = 0;
    std::size_t maxPressureHierarchyLevels = 0;
    std::size_t maxPressureCoarseCells = 0;
    std::size_t pressureIterations = 0;
    std::size_t maxPressureIterations = 0;
    double momentumLinearSolveSeconds = 0;
    double pressureLinearSolveSeconds = 0;
    double solveSeconds = 0;
};

struct FlowIteration2D {
    std::size_t iteration = 0;
    double momentumResidual = 0;
    double continuity = 0;
    double velocityChange = 0;
    double pressureChange = 0;
};

struct FaceMomentum2D {
    bool wall = false; // impermeable no-slip wall or moving lid, excludes slip/inlet/outlet
    double pressure = 0;
    Vector2D advection{};
    Vector2D diffusion{}; // same viscous flux as momentum, selected by viscousStress
};

struct FlowResult2D {
    double time = 0;
    double timeStep = 0; // zero for the steady solver
    double maxCourant = 0;
    std::vector<double> previousU, previousV;
    std::vector<Vector2D> temporalIntegrals;
    bool converged = false;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> p;
    std::vector<double> flux;
    std::vector<FaceMomentum2D> faceMomentum;
    std::vector<Vector2D> sourceIntegrals; // populated only for manufactured verification
    std::vector<FlowIteration2D> history;
    double globalImbalance = 0;
    double globalRelativeImbalance = 0;
    std::size_t outletBackflowFaces = 0;
    double outletInflow = 0; // positive inward volume flux per unit depth
    double forceX = 0;
    double forceY = 0;
    double pressureForceX = 0;
    double pressureForceY = 0;
    double discreteForceX = 0;
    double discreteForceY = 0;
    double reconstructedForceX = 0; // legacy cell-gradient Newtonian diagnostic
    double reconstructedForceY = 0;
    double wallForceX = 0; // all no-slip walls and lid, fluid on boundary
    double wallForceY = 0;
    double wallViscousForceX = 0;
    double wallViscousForceY = 0;
    double domainHeight = 0;
    FlowPerformance2D performance;
};

// Fixed-grid, constant-density laminar SIMPLE, kinematic pressure.
// Viscosity is uniform or explicitly prescribed on faces; no material feedback.
// No turbulence, heat transport, moving mesh or compressibility.
[[nodiscard]] FlowResult2D solveIncompressible2D(
    const FvMesh2D&,
    const FlowControls2D&,
    const std::function<void(const FlowIteration2D&)>& progress = {});

// Backward Euler on a fixed mesh; SIMPLE iterations converge each time step.
// Temporal and inner-iteration face-flux defects retain momentum interpolation.
// The caller must not accept a result unless converged is true.
[[nodiscard]] FlowResult2D advanceIncompressible2D(
    const FvMesh2D&, const FlowControls2D&, const FlowState2D&, double timeStep,
    const std::function<void(const FlowIteration2D&)>& progress = {});
// Physical cases start at rest, with prescribed boundary velocities switched on
// for t>0. taylor-green is a verification-only exact initial vortex on [0,1]^2.
[[nodiscard]] FlowState2D initialIncompressibleState2D(const FvMesh2D&, const FlowControls2D&);

}
