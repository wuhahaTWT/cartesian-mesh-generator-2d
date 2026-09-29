#pragma once

#include "cartmesh2d/fv/ReactingDiffusion2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <iosfwd>

namespace cartmesh2d::fv {

// rho, rho*u, rho*v, rho*E, rho*Y_0 ... rho*Y_(N-1), all per volume.
// E includes formation + sensible internal + kinetic energy and may be negative.
using ReactingConservative2D = std::vector<double>;
struct ReactingPrimitive2D {
    chemistry::GasState gas;
    chemistry::GasProperties properties;
    Vector2D velocity;
    double soundSpeed = 0; // frozen-composition sound speed, m/s
};
[[nodiscard]] ReactingConservative2D reactingConservative2D(const chemistry::GasState&, Vector2D velocity = {});
[[nodiscard]] ReactingPrimitive2D reactingPrimitive2D(chemistry::DetailedGas&, const ReactingConservative2D&);
struct ReactingFaceFlux2D {
    ReactingConservative2D integratedFlux;
    double waveSpeed = 0;
    bool hlleFallback = false;
};
// HLLC with endpoint acoustic bounds and explicit thermodynamic admissibility
// checks on star states. HLLE fallback is reported. No constant-gamma EOS.
[[nodiscard]] ReactingFaceFlux2D reactingFaceFlux2D(chemistry::DetailedGas&,
    const ReactingConservative2D&, const ReactingConservative2D&, Vector2D areaVector);

enum class ReactingBoundaryKind2D { SlipWall, NoSlipWall, Reservoir, ExtrapolatedOutflow };
struct ReactingBoundary2D {
    std::size_t face = 0;
    ReactingBoundaryKind2D kind = ReactingBoundaryKind2D::SlipWall;
    Vector2D velocity; // reservoir velocity or tangential no-slip wall motion
    std::optional<chemistry::GasState> reservoir;
    double wallTemperature = 0; // zero: adiabatic; positive K: isothermal
};
struct ReactingPhysics2D {
    bool chemistry = true, molecularTransport = true; // false only for explicitly frozen/inviscid cases
};
struct ReactingState2D {
    double time = 0;
    std::size_t steps = 0;
    std::string binding; // full resolved mechanism, geometry, boundaries and physical model
    std::vector<ReactingConservative2D> cells;
};
struct ReactingStepControls2D {
    double maximumStep = 1e-6, minimumStep = 1e-15, courant = .35;
    double predictionSafety = .9; // reserve room for wave-speed changes after the first reaction half-step
    std::size_t maximumRetries = 16;
    unsigned order = 2; // MUSCL / SSPRK2; 1 is explicit first-order verification mode
    std::optional<double> endTime;
    chemistry::ChemistryControls chemistry;
};
struct ReactingStepResult2D {
    std::optional<ReactingState2D> accepted;
    std::string failure;
    std::vector<std::string> rejectedReasons;
    double step = 0, combinedCourant = 0;
    std::size_t sourceCalls = 0, hlleFallbacks = 0;
    long successfulSourceInternalSteps = 0;
    // Stage-averaged shared numerical flux, including molecular transport.
    std::vector<ReactingConservative2D> faceFlux;
    ReactingConservative2D beforeIntegral, afterIntegral, boundaryFlux, chemistryChange, balanceError;
    // Integrated, stage-weighted N-1 composition closure. Chemistry closure
    // is already included in chemistryChange; transport closure is separate.
    // Physical mass/element budgets must still include these roundoff effects.
    ReactingConservative2D transportMassClosureChange, chemistryMassClosureChange;
    double maximumMassClosureFraction = 0;
    double absoluteMassClosureIntegral = 0; // sum of weighted |change| * cell area
    // after - before + dt*boundary, divided by integrated initial mass, per element.
    std::vector<double> elementalBalanceError;
};

struct ReactingResidual2D {
    // Instantaneous semi-discrete operators, not an accepted time step.
    // Derivatives are per-volume per-second; face fluxes are edge-integrated.
    std::vector<ReactingConservative2D> faceFlux, transportDerivative, chemistryDerivative, derivative;
    ReactingConservative2D boundaryFlux, chemistryIntegral;
    std::vector<double> transportRate;
    std::size_t hlleFallbacks = 0;
};

// Serial native planar reacting Navier-Stokes development kernel. Holds a
// reference to one DetailedGas context, which must outlive this stepper.
// Chemistry (dt/2) -> conservative SSPRK transport (dt) -> chemistry (dt/2).
// All trial work is private; failure never advances the caller's accepted clock.
class ReactingFlowStepper2D {
public:
    ReactingFlowStepper2D(chemistry::DetailedGas&, FvMesh2D,
        std::vector<ReactingBoundary2D>, ReactingPhysics2D = {});
    [[nodiscard]] ReactingState2D initialState(std::vector<ReactingConservative2D>) const;
    [[nodiscard]] ReactingResidual2D evaluateResidual(const ReactingState2D&, unsigned order = 2);
    [[nodiscard]] ReactingStepResult2D advance(const ReactingState2D&, const ReactingStepControls2D& = {});
    void writeCheckpoint(std::ostream&, const ReactingState2D&) const;
    [[nodiscard]] ReactingState2D readCheckpoint(std::istream&) const;
private:
    struct Stage {
        std::vector<ReactingConservative2D> faceFlux, residual;
        std::vector<double> rate;
        std::size_t fallbacks = 0;
    };
    void validate(const ReactingState2D&) const;
    [[nodiscard]] Stage spatial(const std::vector<ReactingConservative2D>&, unsigned order);
    void react(std::vector<ReactingConservative2D>&, double,
        const chemistry::ChemistryControls&, ReactingStepResult2D&, ReactingConservative2D& change,
        ReactingConservative2D& closureChange, double& maximumClosure, double& absoluteClosure);
    [[nodiscard]] ReactingConservative2D integral(const std::vector<ReactingConservative2D>&) const;
    chemistry::DetailedGas& gas_;
    FvMesh2D mesh_;
    std::vector<std::optional<ReactingBoundary2D>> boundaries_;
    ReactingPhysics2D physics_;
    std::string binding_, mechanismDefinition_, backendVersion_;
    std::optional<ReactingDiffusionOperator2D> diffusion_;
    std::optional<ViscousStressOperator2D> viscous_;
};
} // namespace cartmesh2d::fv
