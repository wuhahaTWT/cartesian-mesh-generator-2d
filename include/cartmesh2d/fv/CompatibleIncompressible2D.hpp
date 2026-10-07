#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>
#include <functional>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
namespace cartmesh2d::fv {
// Explicit compatible steady solver. It does not change the existing solver's
// defaults or certify spatial accuracy. Mesh must pass makeFvMesh2D first.
// One connected fluid component; all original cells/faces are retained.
enum class CompatibleBoundaryKind2D { Velocity, Traction, PseudoTraction, NormalVelocity };
struct CompatibleBoundary2D {
    std::size_t face=0;
    CompatibleBoundaryKind2D kind=CompatibleBoundaryKind2D::Velocity;
    // Physical coordinates. Velocity in m/s; prescribed traction/density in
    // m^2/s^2. Empty is zero. Traction=(nu*(G+G^T)-pI)n;
    // PseudoTraction=(nu*G-pI)n, with the implicit nu*G^T*n retained.
    // Exactly one entry per boundary face; no inferred patch conditions.
    std::function<Vector2D(Point2D)> value;
    // NormalVelocity prescribes u.n and tangential physical traction t.sigma.n,
    // with n the outward unit normal, t=(-n.y,n.x). Empty values are zero:
    // the ordinary impermeable free-slip symmetry condition. value must be
    // empty for this kind; scalar callbacks must be empty for other kinds.
    std::function<double(Point2D)> normalVelocity; // m/s
    std::function<double(Point2D)> tangentialTraction; // m^2/s^2
    // Explicit ordinary-outlet policy. Valid only for Traction/PseudoTraction.
    // Check P1 endpoint normal velocities against -1e-12*Uref, inherited
    // from the existing outlet numerical allowance. Never clips a flux.
    bool rejectBackflow=false;

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
namespace detail {struct CompatibleCheckpointData2D;struct CompatibleCheckpointAccess2D;}
// Immutable accepted steady iterate. This is a numerical continuation, not a
// physical-time checkpoint or an accuracy certificate. Only the native solver
// and the checked reader create it; seeds/rejected candidates cannot be saved.
class CompatibleFlowCheckpoint2D {
public:
    std::size_t acceptedIterations()const;
    double nextPseudoStep()const;
    CompatibleFlowMetrics2D lastAcceptedMetrics()const;
private:
    std::shared_ptr<const detail::CompatibleCheckpointData2D> data_;
    explicit CompatibleFlowCheckpoint2D(std::shared_ptr<const detail::CompatibleCheckpointData2D> data):data_(std::move(data)){}
    friend struct detail::CompatibleCheckpointAccess2D;
};
// Bit-exact double storage and FNV-1a integrity check (not authentication).
// Read bounds allocation by the supplied mesh; resume validates full context.
void writeCompatibleFlowCheckpoint2D(std::ostream&,const CompatibleFlowCheckpoint2D&);
[[nodiscard]] CompatibleFlowCheckpoint2D readCompatibleFlowCheckpoint2D(std::istream&,const FvMesh2D&);
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
    // Runs after complete acceptance and before progress notification.
    // It may publish a checkpoint atomically; exceptions propagate.
    std::function<void(const CompatibleFlowCheckpoint2D&)> checkpointAccepted;
};
// Integrals against 1 and face coordinate s in [-1/2,1/2], in the same
// orientation as the face state. Positive traction acts ON THE FLUID with
// its outward normal. Force on the adjacent body is its negative.
struct CompatibleBoundaryLoad2D {
    std::size_t face=0;
    std::array<Vector2D,2> tractionMoments{},pressureMoments{}; // m^3/s^2
    // tractionMoments - pressureMoments includes viscosity, stabilization,
    // transport reconstruction and lifted-load effects; not pure wall shear.
    Vector2D momentumFlux{}; // m^3/s^2; zero for Stokes
    double volumeFlux=0; // m^2/s
    double torqueOnFluid=0; // m^4/s^2 about momentOrigin
};
struct CompatibleFlowLoads2D {
    std::size_t acceptedIterations=0;
    Point2D momentOrigin{};
    bool absolutePressureReference=false;
    std::vector<CompatibleBoundaryLoad2D> boundaries;
    Vector2D bodyForce{},boundaryTraction{},boundaryMomentumFlux{},momentumImbalance{};
    double boundaryVolumeFlux=0,boundaryTorqueOnFluid=0;
};
// Re-evaluate the ORIGINAL steady equations of a context-matching accepted
// checkpoint, including incomplete accepted iterates. Never accepts a seed or
// rejected trial and does not grant convergence/accuracy qualification. Uses
// kinematic loads per unit depth: multiply forces by density for N/m and
// torques by density for N. Closed-domain pressure uses the saved gauge.
// Cancellation throws; no partial report is returned. User exceptions propagate.
[[nodiscard]] CompatibleFlowLoads2D evaluateCompatibleFlowLoads2D(
    const FvMesh2D&,const CompatibleFlowControls2D&,const CompatibleFlowCheckpoint2D&,
    Point2D momentOrigin={});
enum class CompatibleFlowStop2D { Converged, Cancelled, NonlinearBudget, LinearBudget, BacktrackingBudget, NumericalFailure, BoundaryFailure };
struct CompatibleFlowResult2D {
    CompatibleFlowStop2D stop=CompatibleFlowStop2D::NonlinearBudget;
    std::string reason;
    // A seed is never called an accepted flow or a physical checkpoint.
    std::optional<CompatibleFlowState2D> seed,lastAccepted,lastRejected;
    std::vector<CompatibleFlowIteration2D> iterations;
    std::optional<CompatibleFlowCheckpoint2D> checkpoint;
    std::size_t acceptedIterationsBefore=0; // validated restart history, not this call's attempts
    bool resumed=false;
    bool converged()const{return stop==CompatibleFlowStop2D::Converged;}
};
// initial is an explicit algebraic seed, not a physical restart/checkpoint.
// Invalid inputs throw. Numerical failure/cancellation/budget results preserve
// the last fully accepted iterate; only Converged satisfies all numerical gates.
[[nodiscard]] CompatibleFlowResult2D solveCompatibleIncompressible2D(
    const FvMesh2D&,const CompatibleFlowControls2D&,
    const std::optional<CompatibleFlowState2D>& initial=std::nullopt);
// Restores normalized coefficients and the next pseudo-step without seed
// projection/gauge reset. Binds mesh/order, scales, equation, quadrature,
// globalization, boundary kinds/policies, projected traces and native loads.
// Iteration/linear targets and budgets may change; convergence is re-evaluated.
[[nodiscard]] CompatibleFlowResult2D resumeCompatibleIncompressible2D(
    const FvMesh2D&,const CompatibleFlowControls2D&,const CompatibleFlowCheckpoint2D&);
}
