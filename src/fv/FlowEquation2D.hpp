#pragma once
#include "FlowSolverDetail2D.hpp"

namespace cartmesh2d::fv::solver_detail {

// A reconstruction belongs to one field, one boundary trace and one advective
// flux. Iteration methods may freeze those inputs; the equation does not.
struct FlowReconstruction2D {
    std::vector<Vector2D> gu, gv, gp, pressureForce, stress, faceVelocity;
    Vec pressureFaces;
};

// Geometry caches are owned by a solve, never by a mesh-address/global cache.
// SIMPLE, coupled linear actions, nonlinear checks and postprocessing use the
// same reconstruction and assembly. This class owns no accepted solution.
class FlowEquation2D {
    const FvMesh2D& mesh_;
    detail::ChangingFlowGradientStencil2D velocityU_, velocityV_;
    detail::FlowGradientStencil2D pressure_;
    Vec noFluxDefect_, productU_, productV_;
public:
    FlowEquation2D(const FvMesh2D&, const std::vector<bool>& fixedPressure);
    const detail::FlowGradientStencil2D& pressureStencil() const { return pressure_; }
    std::vector<Vector2D> pressureForce(const Vec&, const Vec& boundaryPressure,
        const std::vector<bool>& fixedPressure) const;
    std::vector<Vector2D> velocityGradient(const Vec&, const Boundary&, bool y);
    FlowReconstruction2D reconstruct(const FlowResult2D&, const Boundary&, const Vec& boundaryPressure);
    FlowReconstruction2D reconstruct(const FlowResult2D&, const FlowControls2D&,
        const Boundary&, const Vec& boundaryPressure, const Vec& advectiveFlux);
    void assembleMomentum(System& u, System& v, const FlowResult2D&,
        const FlowControls2D&, const Boundary&, const FlowReconstruction2D&,
        const Vec& advectiveFlux, const std::vector<Vector2D>& sources,
        const FlowState2D* previous = nullptr, double timeStep = 0) const;
    Vec steadyFlux(const FlowResult2D&, const FlowControls2D&, const Boundary&,
        const FlowReconstruction2D&, const Vec& response, Vec& coefficients);
    // Integrated momentum and divergence, before any row/reference scaling.
    // The gauge row is replaced only for a linear solve; physical checks retain
    // the actual divergence of every cell, including the reference cell.
    void frozenResidual(System& u, System& v, const FlowResult2D&,
        const FlowControls2D&, const Boundary&, const Vec& boundaryPressure,
        const Vec& advectiveFlux, const Vec& response,
        const std::vector<Vector2D>& sources, Vec& residual, Vec& coefficients,
        bool pinPressure);
};

// Shared orthogonal pressure block: SIMPLE correction and coupled block
// preconditioning use the same face coefficient, boundary and gauge rules.
void assemblePressureBlock2D(System&, const FvMesh2D&, const Boundary&, const Vec& coefficients);

} // namespace cartmesh2d::fv::solver_detail
