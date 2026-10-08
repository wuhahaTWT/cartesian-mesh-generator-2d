#include "FlowEquation2D.hpp"

namespace cartmesh2d::fv::solver_detail {

FlowEquation2D::FlowEquation2D(const FvMesh2D& mesh, const std::vector<bool>& fixedPressure)
    : mesh_(mesh), velocityU_(mesh), velocityV_(mesh),
      pressure_(detail::buildFlowGradientStencil2D(mesh, fixedPressure, true)) {}

std::vector<Vector2D> FlowEquation2D::pressureForce(const Vec& pressure,
    const Vec& boundaryPressure, const std::vector<bool>& fixedPressure) const {
    const auto gradient = pressure_.apply(pressure, boundaryPressure);
    return detail::conservativePressureGradient(mesh_,
        detail::pressureFaceValues(mesh_, pressure, gradient, boundaryPressure, fixedPressure));
}

std::vector<Vector2D> FlowEquation2D::velocityGradient(const Vec& field, const Boundary& boundary, bool y) {
    return y ? velocityV_.apply(field, boundary.v, boundary.fixedV, boundary.gradientGroups)
             : velocityU_.apply(field, boundary.u, boundary.fixedU, boundary.gradientGroups);
}

FlowReconstruction2D FlowEquation2D::reconstruct(const FlowResult2D& field,
    const Boundary& boundary, const Vec& boundaryPressure) {
    FlowReconstruction2D result;
    result.gp = pressure_.apply(field.p, boundaryPressure);
    result.gu = velocityGradient(field.u, boundary, false);
    result.gv = velocityGradient(field.v, boundary, true);
    result.pressureFaces = detail::pressureFaceValues(mesh_, field.p, result.gp,
        boundaryPressure, boundary.fixedP);
    result.pressureForce = detail::conservativePressureGradient(mesh_, result.pressureFaces);
    return result;
}

FlowReconstruction2D FlowEquation2D::reconstruct(const FlowResult2D& field,
    const FlowControls2D& controls, const Boundary& boundary,
    const Vec& boundaryPressure, const Vec& advectiveFlux) {
    auto result = reconstruct(field, boundary, boundaryPressure);
    if (controls.viscousStress == ViscousStress2D::Symmetric)
        result.stress = detail::symmetricViscousCorrection(mesh_, field.u, field.v,
            result.gu, result.gv, boundary.u, boundary.v, boundary.fixedU, boundary.fixedV,
            boundary.constantU, boundary.constantV, controls.nu, controls.faceViscosity);
    if (controls.convection == ConvectionScheme2D::FaceLimitedLinearUpwind)
        result.faceVelocity = detail::faceFrameVelocityValues(mesh_, advectiveFlux,
            field.u, field.v, result.gu, result.gv, boundary.u, boundary.v,
            boundary.fixedU, boundary.fixedV);
    return result;
}

void FlowEquation2D::assembleMomentum(System& u, System& v, const FlowResult2D& field,
    const FlowControls2D& controls, const Boundary& boundary, const FlowReconstruction2D& reconstruction,
    const Vec& advectiveFlux, const std::vector<Vector2D>& sources,
    const FlowState2D* previous, double timeStep) const {
    momentum(u, mesh_, controls, boundary, field.u, advectiveFlux, reconstruction.gu,
        reconstruction.pressureForce, sources, reconstruction.stress, reconstruction.faceVelocity,
        false, previous ? &previous->u : nullptr, timeStep);
    momentum(v, mesh_, controls, boundary, field.v, advectiveFlux, reconstruction.gv,
        reconstruction.pressureForce, sources, reconstruction.stress, reconstruction.faceVelocity,
        true, previous ? &previous->v : nullptr, timeStep);
}

Vec FlowEquation2D::steadyFlux(const FlowResult2D& field, const FlowControls2D& controls,
    const Boundary& boundary, const FlowReconstruction2D& reconstruction,
    const Vec& response, Vec& coefficients) {
    // response is V/aP of the unrelaxed equation. The steady fixed-point flux
    // does not contain a lagged under-relaxation defect.
    noFluxDefect_.resize(mesh_.faces.size(), 0);
    return steadyRhieChowFlux(mesh_, controls, boundary, field, response,
        reconstruction.gu, reconstruction.gv, reconstruction.gp,
        reconstruction.pressureForce, noFluxDefect_, coefficients);
}

void FlowEquation2D::frozenResidual(System& u, System& v, const FlowResult2D& field,
    const FlowControls2D& controls, const Boundary& boundary, const Vec& boundaryPressure,
    const Vec& advectiveFlux, const Vec& response, const std::vector<Vector2D>& sources,
    Vec& residual, Vec& coefficients, bool pinPressure) {
    const auto reconstruction = reconstruct(field, controls, boundary, boundaryPressure, advectiveFlux);
    assembleMomentum(u, v, field, controls, boundary, reconstruction, advectiveFlux, sources);
    u.apply(field.u, productU_);
    v.apply(field.v, productV_);
    const auto n = mesh_.cells.size();
    residual.assign(3*n, 0);
    for (std::size_t i=0; i<n; ++i) {
        residual[i] = productU_[i] - u.rhs[i];
        residual[n+i] = productV_[i] - v.rhs[i];
    }
    const auto flux = steadyFlux(field, controls, boundary, reconstruction, response, coefficients);
    for (std::size_t id=0; id<mesh_.faces.size(); ++id) {
        const auto& face = mesh_.faces[id];
        residual[2*n+face.owner] += flux[id];
        if (face.neighbour) residual[2*n+*face.neighbour] -= flux[id];
    }
    if (pinPressure && boundary.closed) residual[2*n] = field.p[0];
}

void assemblePressureBlock2D(System& system, const FvMesh2D& mesh,
    const Boundary& boundary, const Vec& coefficients) {
    ensure(coefficients.size() == mesh.faces.size(), "Pressure coefficient dimensions mismatch");
    system.reset();
    for (std::size_t id=0; id<mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        const auto i = face.owner;
        const double coefficient = coefficients[id];
        if (face.neighbour) {
            const auto j = *face.neighbour;
            system.diag[i] += coefficient;
            system.diag[j] += coefficient;
            system.add(i, j, -coefficient);
            system.add(j, i, -coefficient);
        } else if (boundary.fixedP[id]) system.diag[i] += coefficient;
    }
    if (boundary.closed) system.pin(0);
}

} // namespace cartmesh2d::fv::solver_detail
