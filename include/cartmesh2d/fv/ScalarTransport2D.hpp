#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <optional>

namespace cartmesh2d::fv {
enum class ScalarBoundaryKind2D { Value, DiffusiveFlux };
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
};
struct ScalarTransportIteration2D {
    std::size_t iteration = 0, linearIterations = 0;
    double residualNorm = 0, relativeResidual = 0, maxCellImbalance = 0;
    double maxDiagonalScaledImbalance = 0;
};
struct ScalarTransportResult2D {
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
    const std::vector<double>& previous = {}, double timeStep = 0);
// Rebuild constitutive face fluxes and the original equation residual at supplied
// values without solving or modifying them. Uses the same validation/stopping
// definitions; history contains one entry with iteration=linearIterations=0.
[[nodiscard]] ScalarTransportResult2D evaluateScalarTransport2D(
    const FvMesh2D&, const ScalarTransportProblem2D&, const std::vector<double>& values,
    const ScalarTransportControls2D& = {},
    const std::vector<double>& previous = {}, double timeStep = 0);
}
