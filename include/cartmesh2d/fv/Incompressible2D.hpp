#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"

#include <functional>
#include <string>

namespace cartmesh2d::fv {

enum class ConvectionScheme2D { Upwind, LimitedLinearUpwind, FaceLimitedLinearUpwind };

enum class PressurePreconditioner2D { Jacobi, IncompleteCholesky0, Aggregation, SystemCholesky };
enum class ViscousStress2D { Laplacian, Symmetric };
enum class OutletBackflow2D { Reject, NormalInlet };
enum class FlatPlateTop2D { PressureFarfield, Symmetry };
enum class SteadyAcceleration2D { None, Anderson };
enum class FlowConvergence2D { Strict, Engineering };

// PressureOpening prescribes static kinematic pressure on axis-aligned faces.
// Normal velocity is free; incoming tangential velocity is zero.
enum class FlowBoundaryKind2D { VelocityInlet, PressureOutlet, Wall, MovingWall, SmoothMovingWall, PressureOpening, Symmetry };

// Explicit conditions refer to boundary face IDs in the final FvMesh2D only.
// Velocity and kinematic pressure are physical values, not multiples of speed.
// MovingWall has a constant trace per face. SmoothMovingWall treats supplied
// face-centre velocities as samples of a smooth wall field and reconstructs
// its tangential derivative from the adjacent cell gradient. Both are
// impermeable on the actual polygon faces; neither moves the mesh.
struct FlowBoundaryCondition2D {
    std::size_t face = 0;
    FlowBoundaryKind2D kind = FlowBoundaryKind2D::Wall;
    Vector2D velocity{};
    double pressure = 0;
    std::string name;
};

struct FlowControls2D {
    // duct: planar vertical x-extrema are the inlet/outlet; all remaining
    // boundaries are stationary no-slip walls, including curved walls/holes.
    // It must be selected explicitly and is not a general patch-BC interface.
    std::string scenario = "external";
    // scenario="custom": exactly one entry per boundary face, no internal
    // faces. Pressure openings and symmetry must be axis aligned; other kinds support arbitrary
    // orientations. Moving walls must be tangential.
    // Named groups may have spatially varying values but one physical kind.
    std::vector<FlowBoundaryCondition2D> boundaryConditions;
    double nu = .01;
    double speed = 1;
    double tolerance = 1e-6;
    std::size_t maxIterations = 1500;
    double velocityRelaxation = .6;
    std::size_t pressureCorrectionPasses = 4; // Non-orthogonal pressure corrections per SIMPLE iteration (1..4).
    double pressureRelaxation = .25;
    // Optional safeguarded fixed-point extrapolation. Steady laminar only.
    SteadyAcceleration2D steadyAcceleration = SteadyAcceleration2D::None;
    // Explicit opt-ins preserve existing API/checkpoint and verification cases.
    // Engineering stopping is steady laminar only; it also requires a 50-step
    // window of field/physical-monitor stability and a strict final linear step.
    FlowConvergence2D convergence = FlowConvergence2D::Strict;
    // Laminar steady or transient. Every accepted state requires a strict
    // final linear step; each physical time step resets the forcing sequence.
    // Constitutive material coupling is deliberately unsupported.
    bool adaptiveLinear = false;
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
    // Experimental steady flat-plate case: left uniform inlet, right p=0,
    // bottom symmetry upstream of this x and no-slip downstream. Top defaults
    // to p=0 with free normal velocity, freestream tangential velocity on inflow
    // and zero-gradient velocity on outflow; symmetry is an explicit alternative.
    // Leading edge must coincide with a bottom face endpoint, never bisect a face.
    // Other scenarios require zero. Not yet supported by transient/checkpoint APIs.
    double flatPlateLeadingEdge = 0;
    FlatPlateTop2D flatPlateTop = FlatPlateTop2D::PressureFarfield;
    // Optional cooperative stop, checked after a COMPLETE SIMPLE iteration.
    // A result already meeting every convergence gate wins over this request.
    // Otherwise return the current diagnostic field with stopped=true and
    // converged=false. Does not interrupt a running linear solve; callers still
    // need a hard timeout. Callback exceptions propagate. Empty is unchanged.
    std::function<bool()> stopRequested;
};

// Validates coverage, ownership, names, physical types and impermeability.
// Does not infer or repair missing conditions. The mesh must be validated first.
void validateFlowBoundaryConditions2D(const FvMesh2D&, const FlowControls2D&);

// Accepted state at a physical time, including the conservative face flux.
struct FlowState2D {
    double time = 0;
    std::vector<double> u, v, p, flux;
};

struct FlowPerformance2D {
    std::size_t accelerationCandidates = 0, accelerationAccepted = 0, accelerationRejected = 0;
    std::size_t momentumSolves = 0;
    std::size_t momentumIterations = 0;
    std::size_t maxMomentumIterations = 0;
    std::size_t pressureSolves = 0;
    std::size_t pressureCorrectionPassesSkipped = 0;
    std::size_t pressureFactorizations = 0;
    std::size_t pressureFactorReuses = 0;
    std::size_t pressureCholeskyBuilds = 0, pressureCholeskyRefactors = 0, pressureCholeskyReuses = 0;
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
    // Original unrelaxed momentum residual at the cell setting momentumResidual.
    // Signed components share the same (diagU+diagV)*speed normalization.
    std::size_t momentumWorstCell = 0;
    double momentumResidualX = 0, momentumResidualY = 0;
    // profile only: max true relaxed predictor |b-Au| / (diag*speed),
    // before pressure correction; distinct from nonlinear convergence.
    double momentumPredictorResidual = 0;
    std::size_t momentumPredictorWorstCell = 0;
    // profile only: true PCG residual norm / (speed * shortest face length).
    double pressureLinearResidual = 0;
    double linearRelativeTolerance = 1e-11;
    bool strictLinearStep = true;
    double globalRelativeImbalance = 0;
    std::vector<double> monitors;
};

struct FaceMomentum2D {
    bool wall = false; // impermeable no-slip wall or moving lid, excludes slip/inlet/outlet
    double pressure = 0;
    Vector2D advection{};
    Vector2D diffusion{}; // same viscous flux as momentum, selected by viscousStress
};

// Custom no-slip patches only. Fluid-on-boundary loads divided by density
// and unit depth; torque is about the fixed Cartesian origin (0,0).
struct FlowWallLoad2D {
    std::string name;
    std::size_t faces = 0;
    double length = 0;
    Vector2D pressure{}, viscous{};
    double pressureTorque = 0, viscousTorque = 0;
};

struct FlowResult2D {
    double convergenceReference = 0;
    std::vector<std::string> monitorNames;
    double time = 0;
    double timeStep = 0; // zero for the steady solver
    double maxCourant = 0;
    std::vector<double> previousU, previousV;
    std::vector<Vector2D> temporalIntegrals;
    bool converged = false;
    bool stopped = false;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> p;
    std::vector<double> flux;
    std::vector<FaceMomentum2D> faceMomentum;
    std::vector<FlowWallLoad2D> namedWallLoads;
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

// Optional steady initial iterate, not a physical initial condition or restart.
// An empty flux vector requests target-mesh reconstruction. Optional face
// iterates preserve same-mesh work; prescribed boundary fluxes are validated.
// Every original acceptance gate applies, including to unconverged iterates.
struct FlowInitialGuess2D { std::vector<double> u, v, p, flux; };
[[nodiscard]] FlowResult2D solveIncompressibleFromGuess2D(
    const FvMesh2D&, const FlowControls2D&, const FlowInitialGuess2D&,
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
