#pragma once
#include "cartmesh2d/fv/HybridViscousCell2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <functional>
#include <memory>

namespace cartmesh2d::fv {
struct HybridViscousResult2D {
    std::vector<Vector2D> trace;
    // Unique actual-face flux, outward from owner. Periodic partners are
    // exactly opposite. Residuals and mechanical heating use THESE fluxes.
    std::vector<std::array<double,3>> faceFlux,cellResidual;
    std::vector<double> cellDissipation,cellMechanicalHeating,rate;
    double dissipation=0,mechanicalHeating=0,maximumCellWorkDefect=0;
    double maximumTractionJump=0,maximumBoundaryConstraintResidual=0;
    double traceResidualNorm=0,traceResidualTarget=0; // N/m, per unit depth
    std::size_t iterations=0,traceUnknowns=0;
};

// Research global trace operator for constant-viscosity planar Newtonian gas.
// Every connected fluid component MUST have a prescribed-velocity face. The
// local rigid-rotation kernels can otherwise leave unconstrained trace modes
// (e.g. checkerboard rotations on a periodic Cartesian mesh). No diagonal
// shift, artificial rotation damping, or silent boundary substitution is used.
// Slip means zero normal trace velocity and zero tangential traction; open
// boundaries mean zero traction. Translational periodic pairs are supported.
// This class does not yet connect to the full compressible time integrator.
class HybridViscousOperator2D {
public:
    HybridViscousOperator2D(const FvMesh2D&,const std::vector<ViscousBoundary2D>&,
        double dynamicViscosity);
    [[nodiscard]] HybridViscousResult2D evaluate(const std::vector<Vector2D>& velocity,
        const std::vector<double>& density,double relativeTolerance=1e-12) const;
    [[nodiscard]] std::size_t traceUnknowns() const;
    // Actual symmetric trace matrix entries, in Pa s. Visits diagonal and
    // both off-diagonal triangles; no independently rebuilt PDE matrix.
    void visitTraceMatrix(const std::function<void(std::size_t,std::size_t,double)>&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
} // namespace cartmesh2d::fv
