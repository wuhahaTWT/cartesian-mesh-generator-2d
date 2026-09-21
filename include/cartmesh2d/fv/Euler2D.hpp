#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>
#include <string>

namespace cartmesh2d::fv {

struct IdealGas2D { double gamma=1.4, gasConstant=287.05; };
struct EulerPrimitive2D { double density=1, u=0, v=0, pressure=1; };
// Per-volume densities: rho, rho*u, rho*v, rho*E (total internal + kinetic).
using EulerConservative2D = std::array<double,4>;
enum class EulerBoundaryKind2D { SlipWall, Transmissive, Farfield, Periodic };
struct EulerBoundary2D {
    std::size_t face=0;
    EulerBoundaryKind2D kind=EulerBoundaryKind2D::SlipWall;
    EulerPrimitive2D reference{};
    std::optional<std::size_t> partner;
    std::string name;
};
struct EulerState2D {
    double time=0;
    std::size_t steps=0;
    std::vector<EulerConservative2D> cells;
};
struct EulerStepControls2D {
    double maximumStep=1, minimumStep=1e-14, acousticCourant=.4;
    std::size_t maximumRetries=12;
};
struct EulerStepResult2D {
    EulerState2D state;
    std::vector<EulerConservative2D> previousCells;
    // Unique outward-owner numerical flux integrated over each actual edge.
    // Periodic partners hold exact opposite fluxes, evaluated only once.
    std::vector<EulerConservative2D> faceFlux;
    std::vector<double> faceWaveSpeed;
    EulerConservative2D beforeIntegral{},afterIntegral{},boundaryFlux{},balanceError{};
    double step=0, acousticCourant=0, minimumDensity=0, minimumPressure=0;
    double maximumCellBalanceError=0;
    std::size_t rejectedCandidates=0;
};

void validateIdealGas2D(const IdealGas2D&);
[[nodiscard]] EulerConservative2D eulerConservative2D(const EulerPrimitive2D&,const IdealGas2D& = {});
[[nodiscard]] EulerPrimitive2D eulerPrimitive2D(const EulerConservative2D&,const IdealGas2D& = {});
[[nodiscard]] double eulerSoundSpeed2D(const EulerPrimitive2D&,const IdealGas2D& = {});
void validateEulerBoundaries2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,const IdealGas2D& = {});
// First-order Rusanov flux and forward Euler. Chooses dt from actual face
// acoustic speeds and cell areas; never clips density, pressure or energy.
// A failed trial is retried with smaller dt, leaving the input untouched.
[[nodiscard]] EulerStepResult2D advanceEuler2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,
    const IdealGas2D&,const EulerState2D&,const EulerStepControls2D&);

} // namespace cartmesh2d::fv
