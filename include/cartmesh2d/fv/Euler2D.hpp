#pragma once

#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <array>
#include <string>

namespace cartmesh2d::fv {

struct EulerTransport2D { double thermalConductivity=0, dynamicViscosity=0; }; // W/(m K), Pa s; constant, nonnegative.
struct IdealGas2D { double gamma=1.4, gasConstant=287.05; };
struct EulerPrimitive2D { double density=1, u=0, v=0, pressure=1; };
// Per-volume densities: rho, rho*u, rho*v, rho*E (total internal + kinetic).
using EulerConservative2D = std::array<double,4>;
enum class EulerBoundaryKind2D { SlipWall, Transmissive, Farfield, Periodic, NoSlipWall };
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
    // Optional exact integration horizon; splits the penultimate step if the
    // remaining tail would fall below minimumStep. No accepted clock snapping.
    std::optional<double> endTime;
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
    // Bits 0/1 identify HLLC fallback at the first/second RK stage. Periodic
    // partners share the mask, but evaluation counts include each pair once.
    std::vector<unsigned char> faceHllcFallbackStages;
    std::size_t hllcFallbackEvaluations=0, reconstructionFallbackCells=0;
    double minimumContactRestoration=1;
    EulerConservative2D beforeIntegral{},afterIntegral{},boundaryFlux{},balanceError{};
    double step=0, acousticCourant=0, minimumDensity=0, minimumPressure=0;
    double maximumCellBalanceError=0;
    std::size_t rejectedCandidates=0;
};

void validateIdealGas2D(const IdealGas2D&);
void validateEulerTransport2D(const EulerTransport2D&,const std::vector<EulerBoundary2D>&);
[[nodiscard]] EulerConservative2D eulerConservative2D(const EulerPrimitive2D&,const IdealGas2D& = {});
[[nodiscard]] EulerPrimitive2D eulerPrimitive2D(const EulerConservative2D&,const IdealGas2D& = {});
[[nodiscard]] double eulerSoundSpeed2D(const EulerPrimitive2D&,const IdealGas2D& = {});
void validateEulerBoundaries2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,const IdealGas2D& = {});
// Selectable Rusanov/HLLC and first/second-order spatial and temporal methods.
// Both RK stages obey combined acoustic/heat/viscous rates and cell areas. Reconstruction
// limits slopes; no accepted conserved density, pressure or energy is clipped.
// A failed trial is retried with smaller dt, leaving the input untouched.
[[nodiscard]] EulerStepResult2D advanceEuler2D(const FvMesh2D&,const std::vector<EulerBoundary2D>&,
    const IdealGas2D&,const EulerState2D&,const EulerStepControls2D&,const EulerTransport2D& = {});
// A solver owns an immutable snapshot of mesh, physics and boundaries. Geometry
// and the full thermal Jacobian are prepared once, then reused for every stage.
class EulerStepper2D {
public:
    EulerStepper2D(FvMesh2D,std::vector<EulerBoundary2D>,IdealGas2D = {},EulerTransport2D = {},WallGradient2D = WallGradient2D::Linear);
    [[nodiscard]] EulerStepResult2D advance(const EulerState2D&,const EulerStepControls2D&) const;
private:
    WallGradient2D wallGradient_=WallGradient2D::Linear;
    FvMesh2D mesh_;
    std::vector<EulerBoundary2D> boundaries_;
    IdealGas2D gas_;
    std::optional<HeatConductionOperator2D> heat_;
    std::optional<ViscousStressOperator2D> viscous_;
};

} // namespace cartmesh2d::fv
