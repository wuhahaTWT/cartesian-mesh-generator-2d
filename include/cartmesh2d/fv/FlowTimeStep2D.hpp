#pragma once
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace cartmesh2d::fv {

// Controller limits, not a temporal accuracy estimate. Every proposed step
// must still pass the original nonlinear/continuity gates AND its measured CFL.
struct FlowTimeStepControls2D {
    double maximumStep = .01;
    double minimumStep = .01/1024;
    double targetTime = 1;
    double maximumCourant = 1;
    std::size_t maximumRetries = 10;
    std::size_t maximumAcceptedSteps = 100000;
};

inline void validateFlowTimeStepControls2D(const FlowTimeStepControls2D& c) {
    if (!std::isfinite(c.maximumStep) || !std::isfinite(c.minimumStep) ||
        !std::isfinite(c.targetTime) || !std::isfinite(c.maximumCourant) ||
        c.minimumStep<=0 || c.maximumStep<c.minimumStep || c.targetTime<=0 || c.maximumCourant<=0 ||
        c.maximumRetries>30 || c.maximumAcceptedSteps<1 || c.maximumAcceptedSteps>1000000)
        throw std::runtime_error("Invalid adaptive time-step controls");
}

inline double nextAdaptiveFlowTimeStep2D(const FvMesh2D& mesh, const FlowState2D& state,
                                        const FlowTimeStepControls2D& c) {
    validateFlowTimeStepControls2D(c);
    validateFvMesh2D(mesh);
    if (!std::isfinite(state.time) || state.time<0 || state.flux.size()!=mesh.faces.size())
        throw std::runtime_error("Invalid accepted state for adaptive time-step prediction");
    const double remaining=c.targetTime-state.time;
    if (!(remaining>0)) return 0;
    std::vector<double> absoluteFlux(mesh.cells.size(),0);
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        if (!std::isfinite(state.flux[id])) throw std::runtime_error("Nonfinite accepted flux in time-step prediction");
        const auto& face=mesh.faces[id];const double flux=std::abs(state.flux[id]);
        absoluteFlux[face.owner]+=flux;
        if (face.neighbour) absoluteFlux[*face.neighbour]+=flux;
    }
    double rate=0;
    for (std::size_t i=0;i<mesh.cells.size();++i) {
        const double local=.5*absoluteFlux[i]/mesh.cells[i].area;
        if (!std::isfinite(local)) throw std::runtime_error("Adaptive Courant rate overflow");
        rate=std::max(rate,local);
    }
    // The 0.8 factor leaves room for velocity changes. It cannot replace the
    // acceptance check on the trial's actual conservative face fluxes.
    double step=rate>0?std::min(c.maximumStep,.8*c.maximumCourant/rate):c.maximumStep;
    step=std::min(remaining,std::max(c.minimumStep,step));
    if (!(state.time+step>state.time) || !std::isfinite(state.time+step))
        throw std::runtime_error("Adaptive time step cannot advance floating-point time");
    return step;
}

inline std::optional<double> reducedFlowTimeStep2D(double step, double candidateCourant,
    double remaining, const FlowTimeStepControls2D& c) {
    validateFlowTimeStepControls2D(c);
    if (!std::isfinite(step) || step<=0 || !std::isfinite(candidateCourant) || candidateCourant<0 ||
        !std::isfinite(remaining) || remaining<=0)
        throw std::runtime_error("Invalid rejected adaptive time step");
    const double floor=std::min(c.minimumStep,remaining); // exact final remainder may be smaller
    const double factor=candidateCourant>0?std::min(.5,.8*c.maximumCourant/candidateCourant):.5;
    const double reduced=std::max(floor,step*factor);
    if (!(reduced<step*(1-64*std::numeric_limits<double>::epsilon()))) return {};
    return reduced;
}
}
