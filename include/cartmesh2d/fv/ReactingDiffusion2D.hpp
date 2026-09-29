#pragma once

#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"

namespace cartmesh2d::fv {

enum class ReactingDiffusionBoundaryKind2D { AdiabaticWall, IsothermalWall, Reservoir, ZeroFluxOutflow };
struct ReactingDiffusionBoundary2D {
    std::size_t face = 0;
    ReactingDiffusionBoundaryKind2D kind = ReactingDiffusionBoundaryKind2D::AdiabaticWall;
    double temperature = 0; // K, isothermal wall only
    std::optional<chemistry::GasState> reservoir; // fixed thermodynamic face state only
};

// Integrated over a real 2D edge, per unit out-of-plane depth, owner-outward.
struct ReactingDiffusionFace2D {
    std::vector<double> species; // kg/(m s)
    double conduction = 0, speciesEnthalpy = 0, energy = 0; // W/m
    double rawMassResidual = 0, closureCorrection = 0; // kg/(m s)
    double zeroLimitCorrection = 0; // kg/(m s), exact absent/pure-species limits
};
struct ReactingDiffusionResult2D {
    std::vector<ReactingDiffusionFace2D> faces;
    // Outward sums; conservative time derivative is -residual / true cell area.
    std::vector<std::vector<double>> speciesResidual; // kg/(m s)
    std::vector<double> energyResidual; // W/m
    std::vector<double> boundarySpeciesFlux; // kg/(m s)
    double boundaryEnergyFlux = 0; // W/m
    // Optional absolute row-sum bound [1/s] of the frozen-coefficient
    // thermochemical Jacobian at fixed density, in delta T/T / independent Y.
    // Includes nonorthogonal reconstruction; not a nonlinear stability proof.
    std::vector<double> rate;
};

// Spatial constitutive operator for the reacting-flow solver, not a stand-alone
// flow/time integrator. One shared flux per internal face. LS gradients and the
// existing nonorthogonal face correction use real FvMesh2D geometry. Convex face
// T, log(p), and Y interpolation never clips species. Wall fluxes are imposed
// as boundary constraints; they are not obtained by zeroing the Soret term.
// All states and reservoirs must use the supplied DetailedGas species order.
class ReactingDiffusionOperator2D {
public:
    ReactingDiffusionOperator2D(FvMesh2D, std::vector<ReactingDiffusionBoundary2D>);
    [[nodiscard]] ReactingDiffusionResult2D evaluate(chemistry::DetailedGas&,
        const std::vector<chemistry::GasState>&, bool estimateRate = false) const;
private:
    struct Sample { std::size_t index; bool boundary; Vector2D weight; };
    using Stencil = std::vector<std::vector<Sample>>;
    [[nodiscard]] Stencil makeStencil(bool temperature) const;
    [[nodiscard]] std::vector<Vector2D> gradients(const Stencil&,
        const std::vector<double>&, const std::vector<double>&) const;
    FvMesh2D mesh_;
    std::vector<std::optional<ReactingDiffusionBoundary2D>> boundaries_;
    Stencil temperatureStencil_, compositionStencil_;
};

} // namespace cartmesh2d::fv
