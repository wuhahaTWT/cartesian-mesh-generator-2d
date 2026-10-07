#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>
#include <functional>
#include <optional>
#include <string>
namespace cartmesh2d::fv {
// Explicit compatible steady solver. It does not change the existing solver's
// defaults or certify spatial accuracy. Mesh must pass makeFvMesh2D first.
// One connected fluid component; all original cells/faces are retained.
enum class CompatibleBoundaryKind2D { Velocity, Traction, PseudoTraction };
struct CompatibleBoundary2D {
    std::size_t face=0;
    CompatibleBoundaryKind2D kind=CompatibleBoundaryKind2D::Velocity;
    // Physical coordinates. Velocity in m/s; prescribed traction/density in
    // m^2/s^2. Empty is zero. Traction=(nu*(G+G^T)-pI)n;
    // PseudoTraction=(nu*G-pI)n, with the implicit nu*G^T*n retained.
    // Exactly one entry per boundary face; no inferred patch conditions.
    std::function<Vector2D(Point2D)> value;
};
enum class CompatibleEquation2D { Stokes, NavierStokes };
enum class CompatibleGlobalization2D { Backtracking, PseudoTime };
enum class CompatiblePressureInverse2D { ViscousMass, DiagonalSchur };
struct CompatibleFlowState2D {
    // Per cell: u0,uX,uY,v0,vX,vY,p0,pX,pY in basis (1,dx/h,dy/h).
    // Physical velocity and kinematic pressure. h is the actual diameter.
    std::vector<std::array<double,9>> cells;
    // Per original face: u0,us,v0,vs in s in [-1/2,1/2], oriented by
    // point=face.centre+s*(-Sy,Sx). These are independent P1 trace DOFs.
    std::vector<std::array<double,4>> faces;
};
struct CompatibleFlowMetrics2D {
    // All dimensionless after x/Lref, u/Uref, p/Uref^2 conversion.
    // Cell weak momentum / area; free-face weak momentum / face length;
    // maximum P1 weak divergence at the original fan vertices.
    double cellMomentum=0,faceMomentum=0,divergence=0;
    double stateChange=0; // max normalized cell/face polynomial coefficient change
    double residualNorm=0; // raw free weak-equation L2 norm, for globalization only
};
struct CompatibleFlowIteration2D {
    std::size_t iteration=0,linearRestarts=0,matrixProducts=0,trials=0;
    double linearRelativeResidual=0,pseudoStep=0,alpha=1;
    bool accepted=false;
    std::optional<CompatibleFlowMetrics2D> metrics; // absent until original equations were evaluated
};
struct CompatibleFlowControls2D {
    double viscosity=.01,referenceLength=1,referenceVelocity=1;
    CompatibleEquation2D equation=CompatibleEquation2D::NavierStokes;
    CompatibleGlobalization2D globalization=CompatibleGlobalization2D::PseudoTime;
    CompatiblePressureInverse2D pressureInverse=CompatiblePressureInverse2D::ViscousMass;
    // Iteration targets inherited from the compatible research solver,
    // applied separately to the normalized quantities above. Not accuracy bounds.
    double equationTolerance=1e-9,stateTolerance=1e-9,linearTolerance=1e-13;
    std::size_t maximumIterations=40,maximumLinearRestarts=50,krylovDirections=60,maximumBacktracks=15;
    // Pseudo-time is in Lref/Uref units; it is not a physical time advance.
    double initialPseudoStep=.1,maximumPseudoStep=1e6,armijo=1e-4;
    int quadratureOrder=6;
    std::vector<CompatibleBoundary2D> boundaries;
    std::function<Vector2D(Point2D)> acceleration; // physical m/s^2, empty is zero
    // Checked during preparation, assembly, Krylov and trial evaluation.
    // Callback exceptions propagate; cancellation does not accept an unfinished trial.
    std::function<bool()> stopRequested;
    std::function<void(const CompatibleFlowIteration2D&)> iterationAccepted;
};
enum class CompatibleFlowStop2D { Converged, Cancelled, NonlinearBudget, LinearBudget, BacktrackingBudget, NumericalFailure };
struct CompatibleFlowResult2D {
    CompatibleFlowStop2D stop=CompatibleFlowStop2D::NonlinearBudget;
    std::string reason;
    // A seed is never called an accepted flow or a physical checkpoint.
    std::optional<CompatibleFlowState2D> seed,lastAccepted,lastRejected;
    std::vector<CompatibleFlowIteration2D> iterations;
    bool converged()const{return stop==CompatibleFlowStop2D::Converged;}
};
// initial is an explicit algebraic seed, not a physical restart/checkpoint.
// Invalid inputs throw. Numerical failure/cancellation/budget results preserve
// the last fully accepted iterate; only Converged satisfies all numerical gates.
[[nodiscard]] CompatibleFlowResult2D solveCompatibleIncompressible2D(
    const FvMesh2D&,const CompatibleFlowControls2D&,
    const std::optional<CompatibleFlowState2D>& initial=std::nullopt);
}
