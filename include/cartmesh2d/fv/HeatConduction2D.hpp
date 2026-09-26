#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"

namespace cartmesh2d::fv {
enum class HeatBoundaryKind2D { Insulated, Temperature, OutwardFlux, Periodic };
struct HeatBoundary2D {
    std::size_t face=0;
    HeatBoundaryKind2D kind=HeatBoundaryKind2D::Insulated;
    double value=0; // K for Temperature; W/m^2 OUT of fluid for OutwardFlux.
    std::optional<std::size_t> partner;
};
struct HeatConductionResult2D {
    std::vector<double> faceHeatFlux; // W/m, outward owner, per unit depth.
    std::vector<double> cellResidual; // Sum of outward shared face heat fluxes.
    std::vector<double> rate; // Half absolute row sum / (area * volumetric cv), 1/s.
};
// Immutable geometry/physical-boundary cache. The evaluated Fourier flux uses
// least-squares temperature gradients with an over-relaxed nonorthogonal split.
// This corrected operator is linear-exact, but NOT generally an M-matrix. The
// rate bounds its full row norm; callers must still check stage admissibility.
class HeatConductionOperator2D {
public:
    HeatConductionOperator2D(const FvMesh2D&,const std::vector<HeatBoundary2D>&,double conductivity);
    [[nodiscard]] HeatConductionResult2D evaluate(const std::vector<double>& temperature,
                                                 const std::vector<double>& volumetricHeatCapacity) const;
    [[nodiscard]] std::size_t nonMonotoneRows() const { return nonMonotoneRows_; }
private:
    struct Sample { std::optional<std::size_t> cell; double value=0; Vector2D weight{}; };
    struct Gradient { Vector2D constant{}; std::vector<Sample> samples; };
    struct FaceData {
        std::size_t owner=0;std::optional<std::size_t> neighbour,partner;
        HeatBoundaryKind2D kind=HeatBoundaryKind2D::Insulated;
        double value=0,length=0,transmissibility=0,weight=0;
        Vector2D correction{};
    };
    double conductivity_=0;
    std::vector<double> areas_,rowNorm_;
    std::vector<Gradient> gradients_;
    std::vector<FaceData> faces_;
    std::size_t nonMonotoneRows_=0;
};
} // namespace cartmesh2d::fv
