#pragma once

#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <array>
#include <string>
#include <functional>
#include <memory>

namespace cartmesh2d::fv {

struct EulerTransport2D { double thermalConductivity=0, dynamicViscosity=0; }; // W/(m K), Pa s; constant, nonnegative.
struct IdealGas2D { double gamma=1.4, gasConstant=287.05; };
struct EulerPrimitive2D { double density=1, u=0, v=0, pressure=1; };
// Per-volume densities: rho, rho*u, rho*v, rho*E (total internal + kinetic).
using EulerConservative2D = std::array<double,4>;
// Append kinds: checkpoint bindings persist the existing integer values.
enum class EulerBoundaryKind2D { SlipWall, Transmissive, Farfield, Periodic, NoSlipWall, PressureOutlet, TotalInlet };
struct EulerBoundary2D {
    std::size_t face=0;
    EulerBoundaryKind2D kind=EulerBoundaryKind2D::SlipWall;
    EulerPrimitive2D reference{};
    std::optional<std::size_t> partner;
    std::string name;
    HeatBoundaryKind2D thermalKind=HeatBoundaryKind2D::Insulated;
    double thermalValue=0;
    Vector2D wallVelocity{}; // Tangential only on a static no-slip wall.
};
struct EulerState2D {
    double time=0;
    std::size_t steps=0;
    std::vector<EulerConservative2D> cells;
};
enum class EulerFluxScheme2D { Rusanov, Hllc };
enum class EulerTimeStepControl2D { Legacy, StageGuarded };
enum class EulerTimeIntegrator2D { Explicit, Sdirk2 };
// FrozenFluxIlu0 approximates only the Newton right preconditioner. The true
// residual, accepted RK flux and nonlinear tolerance are unchanged.
enum class EulerImplicitPreconditioner2D { Diagonal, FrozenFluxIlu0 };
// Optional research transport forms. HybridHeat keeps corrected viscous stress;
// Hybrid also solves shared velocity traces. Neither changes the default.
enum class EulerDiffusionScheme2D { Corrected, HybridHeat, Hybrid };
struct EulerFaceFlux2D {
    EulerConservative2D integratedFlux{};
    double waveSpeed=0;
    bool hllcFallback=false;
};
// A single outward-owner flux per actual edge; reversing states and the area
// vector reverses the flux. HLLC falls back visibly if its star states are invalid.
[[nodiscard]] EulerFaceFlux2D eulerFaceFlux2D(const EulerConservative2D& left,
    const EulerConservative2D& right, Vector2D areaVector, const IdealGas2D&, EulerFluxScheme2D, double contactRestoration=1);
struct EulerStepControls2D {
    // Legacy name: with k>0 this caps the combined acoustic + thermal + viscous rate.
    double maximumStep=1, minimumStep=1e-14, acousticCourant=.4;
    std::size_t maximumRetries=12;
    EulerFluxScheme2D fluxScheme=EulerFluxScheme2D::Rusanov;
    unsigned order=1; // 1: constant/forward Euler; 2: limited linear/SSPRK(2,2).
    WallGradient2D wallGradient=WallGradient2D::Linear;
    EulerDiffusionScheme2D diffusionScheme=EulerDiffusionScheme2D::Corrected;
    // Optional exact integration horizon; splits the penultimate step if the
    // remaining tail would fall below minimumStep. No accepted clock snapping.
    std::optional<double> endTime;
    // Optional 5% headroom for changing stage rates and rate-based CFL retries.
    // Stage positivity/CFL limits are unchanged; no history is needed on restart.
    EulerTimeStepControl2D timeStepControl=EulerTimeStepControl2D::Legacy;
    EulerTimeIntegrator2D integrator=EulerTimeIntegrator2D::Explicit;
    EulerImplicitPreconditioner2D implicitPreconditioner=EulerImplicitPreconditioner2D::Diagonal;
    // SDIRK uses maximumStep as a uniform physical step, never a local clock.
    // Dimensionless stage defect relative to fixed rho/rho*c/rho*E scales.
    double nonlinearTolerance=2e-14;
    std::size_t maximumNewtonIterations=16, maximumKrylovIterations=240;
    // Optional cancellation/budget callback, checked inside implicit residuals.
    std::function<bool()> interrupted;
};
struct EulerStepResult2D {
    EulerState2D state;
    std::vector<EulerConservative2D> previousCells;
    // Unique outward-owner numerical flux integrated over each actual edge.
    // Periodic partners hold exact opposite fluxes, evaluated only once.
    std::vector<EulerConservative2D> faceFlux;
    std::vector<double> faceWaveSpeed,faceHeatFlux,cellHeatRate,cellViscousRate;
    std::vector<std::array<double,3>> faceViscousFlux;
    double thermalCourant=0,viscousCourant=0,combinedCourant=0,boundaryHeat=0,boundaryViscousWork=0;
    std::size_t heatNonMonotoneRows=0,quadraticHeatWalls=0,quadraticViscousWalls=0;
    // The legacy row count is meaningful only when this is true. A hybrid
    // Schur operator is nonlocal and has not been enumerated for this diagnostic.
    bool heatMonotonicityAssessed=true;
    // Bits 0/1 identify HLLC fallback at the first/second RK stage. Periodic
    // partners share the mask, but evaluation counts include each pair once.
    std::vector<unsigned char> faceHllcFallbackStages;
    std::size_t hllcFallbackEvaluations=0, reconstructionFallbackCells=0;
    double minimumContactRestoration=1;
    EulerConservative2D beforeIntegral{},afterIntegral{},boundaryFlux{},balanceError{};
    double step=0, acousticCourant=0, minimumDensity=0, minimumPressure=0;
    double maximumCellBalanceError=0;
    std::size_t rejectedCandidates=0;
    std::string lastRejectedReason;
    std::size_t cflRejectedCandidates=0,spatialEvaluations=0;
    std::size_t nonlinearIterations=0,linearIterations=0;
};
struct EulerResidualDiagnostics2D {
    double rate=0; // s^-1, maximum component residual with local physical scales.
    EulerConservative2D boundaryFlux{},integralScale{};
};

void validateIdealGas2D(const IdealGas2D&);
void validateEulerTransport2D(const EulerTransport2D&,const std::vector<EulerBoundary2D>&);
[[nodiscard]] EulerConservative2D eulerConservative2D(const EulerPrimitive2D&,const IdealGas2D& = {});
[[nodiscard]] EulerPrimitive2D eulerPrimitive2D(const EulerConservative2D&,const IdealGas2D& = {});
[[nodiscard]] double eulerSoundSpeed2D(const EulerPrimitive2D&,const IdealGas2D& = {});
// Subsonic outflow: prescribed static pressure, outgoing acoustic invariant,
// interior entropy and tangential velocity. Supersonic outflow extrapolates.
// Backflow/choking of the subsonic boundary is rejected, never clipped.
[[nodiscard]] EulerPrimitive2D eulerPressureOutletState2D(const EulerPrimitive2D&,
    double pressure, Vector2D outwardArea, const IdealGas2D& = {});
// Three incoming conditions: reference entropy, total enthalpy and tangential
// velocity. Outgoing J+ comes from the interior. Reference is a reservoir state,
// not a prescription of all four boundary primitive variables.
[[nodiscard]] EulerPrimitive2D eulerTotalInletState2D(const EulerPrimitive2D&,
    const EulerPrimitive2D& reference,Vector2D outwardArea,const IdealGas2D& = {});
void validateEulerBoundaries2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,const IdealGas2D& = {});
// Selectable Rusanov/HLLC and first/second-order spatial and temporal methods.
// Explicit RK stages obey combined acoustic/heat/viscous rates and cell areas.
// SDIRK2 solves both coupled stages with feasible Newton-Krylov and the same residual.
// Reconstruction
// limits slopes; no accepted conserved density, pressure or energy is clipped.
// A failed trial is retried with smaller dt, leaving the input untouched.
[[nodiscard]] EulerStepResult2D advanceEuler2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,
    const IdealGas2D&,const EulerState2D&,const EulerStepControls2D&,const EulerTransport2D& = {});
// A solver owns an immutable snapshot of mesh, physics and boundaries. Local
// transport geometry is reused; hybrid heat currently assembles and solves its
// global face system on every residual evaluation, including Newton products.
struct EulerDiffusionOperators2D;
class EulerStepper2D {
public:
    EulerStepper2D(FvMesh2D,std::vector<EulerBoundary2D>,IdealGas2D = {},EulerTransport2D = {},WallGradient2D = WallGradient2D::Linear,
        EulerDiffusionScheme2D = EulerDiffusionScheme2D::Corrected);
    [[nodiscard]] EulerStepResult2D advance(const EulerState2D&,const EulerStepControls2D&) const;
    // Maximum |R/(V scale)|, in s^-1, scales rho, rho*c, rho*c, rho*E.
    [[nodiscard]] double residualRate(const EulerState2D&,const EulerStepControls2D&) const;
    [[nodiscard]] EulerResidualDiagnostics2D diagnostics(const EulerState2D&,const EulerStepControls2D&) const;
private:
    WallGradient2D wallGradient_=WallGradient2D::Linear;
    EulerDiffusionScheme2D diffusionScheme_=EulerDiffusionScheme2D::Corrected;
    FvMesh2D mesh_;
    std::vector<EulerBoundary2D> boundaries_;
    IdealGas2D gas_;
    std::shared_ptr<const EulerDiffusionOperators2D> diffusion_;
};

} // namespace cartmesh2d::fv
