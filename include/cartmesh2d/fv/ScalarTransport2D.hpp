#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <optional>
#include <memory>

namespace cartmesh2d::fv {
// Explicit serial workspace; owns its sparse pattern, matrix and Krylov arrays.
// Each call still validates the current mesh/inputs and reassembles all numeric
// coefficients. Exact connectivity changes rebuild the pattern. Never share
// concurrently, or move/destroy during a call; recursive use is rejected.
class ScalarTransportWorkspace2D {
public:
    ScalarTransportWorkspace2D();
    ~ScalarTransportWorkspace2D();
    ScalarTransportWorkspace2D(ScalarTransportWorkspace2D&&) noexcept;
    ScalarTransportWorkspace2D& operator=(ScalarTransportWorkspace2D&&) noexcept;
    ScalarTransportWorkspace2D(const ScalarTransportWorkspace2D&)=delete;
    ScalarTransportWorkspace2D& operator=(const ScalarTransportWorkspace2D&)=delete;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    friend struct ScalarTransportWorkspaceAccess2D;
};
enum class ScalarBoundaryKind2D { Value, DiffusiveFlux };
enum class ScalarPreconditioner2D { Jacobi, ILU0 };
struct ScalarBoundary2D {
    ScalarBoundaryKind2D kind = ScalarBoundaryKind2D::Value;
    // Value: scalar at face. DiffusiveFlux: outward -D grad(s).n per unit length.
    // For temperature, D=k/(rho*cp), so prescribed physical heat flux is q/(rho*cp).
    double value = 0;
    // Required on an inflowing DiffusiveFlux boundary; never invent a backflow value.
    std::optional<double> inflowValue;
};
struct ScalarTransportProblem2D {
    // d(s)/dt + div(U*s - D grad(s)) + sinkRate*s = source, fixed density.
    double diffusivity = .01;
    std::vector<double> volumeFlux; // one owner-outward integrated U.S per mesh face
    // Supply callbacks OR indexed data, never both. Indexed data supports
    // exact restart validation of spatially varying fixed source/BC data.
    std::vector<double> sourceDensity;
    std::vector<ScalarBoundary2D> boundaryData;
    std::function<double(Point2D)> source;
    std::function<ScalarBoundary2D(std::size_t, const Face&)> boundary;
    // Optional prescribed face coefficients, including boundary faces. Empty
    // uses diffusivity; otherwise exactly one finite positive D per face.
    // The caller defines material-interface interpolation (no hidden averaging).
    // diffusivity remains a finite positive reference value in either mode.
    std::vector<double> faceDiffusivity;
    // Optional prescribed nonnegative loss rate [1/time] per cell. Empty is
    // zero. Integrated implicitly, and included in the reported true balance.
    // A caller linearizing nonlinear losses must also check nonlinear convergence.
    std::vector<double> sinkRate;
};
struct ScalarTransportControls2D {
    ConvectionScheme2D convection = ConvectionScheme2D::Upwind;
    std::size_t maxCorrections = 1000;
    double relaxation = .8;
    double relativeTolerance = 1e-9, absoluteTolerance = 1e-12;
    double cellTolerance = 1e-9; // imbalance / unrelaxed diagonal, in scalar units
    double carrierRelativeTolerance = 1e-8, carrierAbsoluteTolerance = 1e-12;
    bool profile = false; // diagnostic only; no numerical effect
    ScalarPreconditioner2D preconditioner = ScalarPreconditioner2D::Jacobi;
};
struct ScalarTransportIteration2D {
    std::size_t iteration = 0, linearIterations = 0;
    double residualNorm = 0, relativeResidual = 0, maxCellImbalance = 0;
    double maxDiagonalScaledImbalance = 0;
    // Additional original-operator check when an apparently converged double
    // field requests accuracy near its representation floor. Face-flux balance
    // remains independently required; a precise linear candidate is insufficient.
    bool matrixAudited = false;
    double matrixResidualNorm = 0, matrixMaxDiagonalScaledImbalance = 0;
};
// totalSeconds includes all subphases. setup includes validation and assembly;
// faceFluxSeconds includes all reconstructions, before and after corrections.
struct ScalarTransportPerformance2D {
    std::size_t calls=0, patternBuilds=0, patternReuses=0, linearIterations=0, ilu0Builds=0, ilu0Reuses=0;
    double totalSeconds=0, setupSeconds=0, linearSeconds=0, faceFluxSeconds=0;
    void add(const ScalarTransportPerformance2D& p) {
        calls+=p.calls;patternBuilds+=p.patternBuilds;patternReuses+=p.patternReuses;linearIterations+=p.linearIterations;
        ilu0Builds+=p.ilu0Builds;ilu0Reuses+=p.ilu0Reuses;
        totalSeconds+=p.totalSeconds;setupSeconds+=p.setupSeconds;
        linearSeconds+=p.linearSeconds;faceFluxSeconds+=p.faceFluxSeconds;
    }
};
struct ScalarTransportResult2D {
    ScalarTransportPerformance2D performance;
    bool converged = false;
    std::vector<double> values, advectiveFlux, diffusiveFlux;
    std::vector<double> sourceIntegrals, temporalIntegrals;
    std::vector<ScalarTransportIteration2D> history;
    double boundaryFlux = 0, sourceIntegral = 0, temporalIntegral = 0, globalBalance = 0;
    double maxCarrierImbalance = 0, minValue = 0, maxValue = 0, maxCourant = 0;
    std::vector<double> sinkIntegrals;
    double sinkIntegral = 0;
};
// Conservative transport on a validated final mesh. No clipping or hidden sinks.
// Empty previous means steady; otherwise backward Euler with dt>0. Callbacks are
// sampled at the new physical time by the caller. Do not accept unconverged output.
// A steady component needs a value boundary, prescribed inflow, or positive loss.
// Outflow advection uses owner reconstruction; boundary values constrain diffusion.
[[nodiscard]] ScalarTransportResult2D solveScalarTransport2D(
    const FvMesh2D&, const ScalarTransportProblem2D&,
    const ScalarTransportControls2D& = {},
    const std::vector<double>& previous = {}, double timeStep = 0, ScalarTransportWorkspace2D* workspace = nullptr);
// Iterative initial guess for a STEADY equation only; not a previous physical
// time layer. Original source, operator and convergence requirements apply.
// Explicit initial values must contain exactly one finite value per cell.
[[nodiscard]] ScalarTransportResult2D solveSteadyScalarTransportFromInitial2D(
    const FvMesh2D&, const ScalarTransportProblem2D&,
    const std::vector<double>& initial, const ScalarTransportControls2D& = {},
    ScalarTransportWorkspace2D* workspace = nullptr);
// Rebuild constitutive face fluxes and the original equation residual at supplied
// values without solving or modifying them. Uses the same validation/stopping
// definitions; history contains one entry with iteration=linearIterations=0.
[[nodiscard]] ScalarTransportResult2D evaluateScalarTransport2D(
    const FvMesh2D&, const ScalarTransportProblem2D&, const std::vector<double>& values,
    const ScalarTransportControls2D& = {},
    const std::vector<double>& previous = {}, double timeStep = 0, ScalarTransportWorkspace2D* workspace = nullptr);
}
