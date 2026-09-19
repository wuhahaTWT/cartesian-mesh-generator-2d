#include "cartmesh2d/fv/ThermalFlow2D.hpp"
#include <cmath>
#include <stdexcept>

namespace cartmesh2d::fv {
void validateThermalSetup2D(const FvMesh2D& mesh,const ThermalSetup2D& setup) {
    validateFvMesh2D(mesh);
    if (!std::isfinite(setup.diffusivity) || setup.diffusivity<=0 ||
        setup.sourceDensity.size()!=mesh.cells.size() || setup.boundary.size()!=mesh.faces.size())
        throw std::invalid_argument("Thermal setup invalid diffusivity or field dimensions");
    for (double s:setup.sourceDensity)
        if (!std::isfinite(s)) throw std::invalid_argument("Thermal setup nonfinite source");
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        if (mesh.faces[id].neighbour) continue;
        const auto& b=setup.boundary[id];
        if ((b.kind!=ScalarBoundaryKind2D::Value && b.kind!=ScalarBoundaryKind2D::DiffusiveFlux) ||
            !std::isfinite(b.value) || (b.inflowValue && !std::isfinite(*b.inflowValue)))
            throw std::invalid_argument("Thermal setup invalid boundary");
    }
}
ThermalFlowResult2D advanceThermalFlow2D(const FvMesh2D& mesh,
    const FlowControls2D& flowControls,const ThermalSetup2D& setup,
    const ScalarTransportControls2D& scalarControls,const ThermalFlowState2D& previous,
    double timeStep,const std::function<void(const FlowIteration2D&)>& progress) {
    validateThermalSetup2D(mesh,setup);
    if (previous.scalar.size()!=mesh.cells.size())
        throw std::invalid_argument("Thermal state invalid scalar dimensions");
    for (double s:previous.scalar)
        if (!std::isfinite(s)) throw std::invalid_argument("Thermal state nonfinite scalar");
    ThermalFlowResult2D result;
    result.flow=advanceIncompressible2D(mesh,flowControls,previous.flow,timeStep,progress);
    if (!result.flow.converged) return result;
    ScalarTransportProblem2D scalar;
    scalar.diffusivity=setup.diffusivity;
    scalar.sourceDensity=setup.sourceDensity;
    scalar.boundaryData=setup.boundary;
    scalar.volumeFlux=result.flow.flux;
    result.scalar=solveScalarTransport2D(mesh,scalar,scalarControls,previous.scalar,timeStep);
    if (result.scalar.converged)
        result.accepted=ThermalFlowState2D{
            {result.flow.time,result.flow.u,result.flow.v,result.flow.p,result.flow.flux},result.scalar.values};
    return result;
}
}
