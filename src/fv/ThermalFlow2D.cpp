#include "cartmesh2d/fv/ThermalFlow2D.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace cartmesh2d::fv {
void validateThermalSetup2D(const FvMesh2D& mesh,const ThermalSetup2D& setup) {
    validateFvMesh2D(mesh);
    if (!std::isfinite(setup.diffusivity) || setup.diffusivity<=0 ||
        setup.sourceDensity.size()!=mesh.cells.size() || setup.boundary.size()!=mesh.faces.size())
        throw std::invalid_argument("Thermal setup invalid diffusivity or field dimensions");
    double last=-1;
    for (const auto& event:setup.events) {
        if (!std::isfinite(event.time) || event.time<0 || event.time<=last)
            throw std::invalid_argument("Thermal events must have finite strictly increasing nonnegative times");
        last=event.time;
        ThermalSetup2D data{setup.diffusivity,event.sourceDensity,event.boundary,{}};
        validateThermalSetup2D(mesh,data);
    }
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
    const std::vector<double>* sources=&setup.sourceDensity;
    const std::vector<ScalarBoundary2D>* boundary=&setup.boundary;
    for (const auto& event:setup.events) {
        if (event.time<=previous.flow.time) {sources=&event.sourceDensity;boundary=&event.boundary;}
        else if (timeStep>event.time-previous.flow.time &&
                 timeStep-(event.time-previous.flow.time)>8*std::numeric_limits<double>::epsilon()*
                     std::max(event.time,timeStep))
            throw std::invalid_argument("Thermal step crosses event; split at event time");
    }
    ThermalFlowResult2D result;
    result.flow=advanceIncompressible2D(mesh,flowControls,previous.flow,timeStep,progress);
    if (!result.flow.converged) return result;
    ScalarTransportProblem2D scalar;
    scalar.diffusivity=setup.diffusivity;
    scalar.sourceDensity=*sources;
    scalar.boundaryData=*boundary;
    scalar.volumeFlux=result.flow.flux;
    result.scalar=solveScalarTransport2D(mesh,scalar,scalarControls,previous.scalar,timeStep);
    if (result.scalar.converged)
        result.accepted=ThermalFlowState2D{
            {result.flow.time,result.flow.u,result.flow.v,result.flow.p,result.flow.flux},result.scalar.values,{}};
    return result;
}
}

namespace cartmesh2d::fv {
ThermalControlledResult2D advanceControlledThermalFlow2D(const FvMesh2D& mesh,
    const FlowControls2D& fc,const ThermalSetup2D& setup,const ScalarTransportControls2D& sc,
    const ThermalFlowState2D& previous,const ThermalTimeControls2D& c,const std::function<bool()>& cancelled,
    const std::function<void(const ThermalAttempt2D&)>& observeAttempt) {
    validateThermalSetup2D(mesh,setup);
    validateFlowTimeStepControls2D(c.limits);
    for (double v:{c.temperatureScale,c.velocityScale,c.relativeTolerance,
                   c.temperatureAbsoluteTolerance,c.velocityAbsoluteTolerance})
        if (!std::isfinite(v) || v<=0) throw std::invalid_argument("Invalid thermal time-error scale/tolerance");
    auto limits=c.limits;
    for (const auto& event:setup.events)
        if (event.time>previous.flow.time) {limits.targetTime=std::min(limits.targetTime,event.time);break;}
    const std::vector<double> signature={c.limits.maximumStep,c.limits.minimumStep,c.limits.maximumCourant,
        c.estimateError?1.:0.,c.temperatureScale,c.velocityScale,c.relativeTolerance,
        c.temperatureAbsoluteTolerance,c.velocityAbsoluteTolerance,fc.velocityRelaxation};
    const bool reuse=previous.controller && previous.controller->controls==signature;
    double dt=nextAdaptiveFlowTimeStep2D(mesh,previous.flow,limits);
    if(reuse) {
        const auto& h=*previous.controller;
        if(!std::isfinite(h.nextStep)||h.nextStep<=0 || !std::isfinite(h.velocityRelaxation)||h.velocityRelaxation<=0||h.velocityRelaxation>1)
            throw std::invalid_argument("Invalid accepted thermal controller history");
        dt=std::min(dt,std::max(std::min(limits.minimumStep,limits.targetTime-previous.flow.time),h.nextStep));
    }
    if (dt==0) throw std::invalid_argument("Thermal target must exceed accepted time");
    ThermalControlledResult2D out;
    const auto checkCancel=[&]() {
        if (cancelled && cancelled()) throw std::runtime_error("Thermal calculation cancelled; accepted state retained");
    };
    auto trialFlow=fc;
    if(reuse)trialFlow.velocityRelaxation=std::min(fc.velocityRelaxation,previous.controller->velocityRelaxation);
    auto trialScalar=sc;
    trialScalar.stopRequested=[&] {return (cancelled && cancelled()) || (sc.stopRequested && sc.stopRequested());};
    const auto run=[&](const ThermalFlowState2D& start,double h) {
        checkCancel();
        auto r=advanceThermalFlow2D(mesh,trialFlow,setup,trialScalar,start,h,[&](const FlowIteration2D&){checkCancel();});
        checkCancel(); return r;
    };
    for (std::size_t retry=0;;++retry) {
        auto full=run(previous,dt);
        double error=0, courant=std::max(full.flow.maxCourant,full.scalar.maxCourant);
        std::string reason=!full.flow.converged?"flow":!full.scalar.converged?"scalar":"accepted";
        if (full.accepted && (!std::isfinite(courant) || courant>limits.maximumCourant)) reason="courant";
        if (reason=="accepted" && c.estimateError) {
            const auto half=run(previous,.5*dt);
            if (!half.accepted) reason=half.flow.converged?"scalar-half":"flow-half";
            else {
                const auto second=run(*half.accepted,.5*dt);
                if (!second.accepted) reason=second.flow.converged?"scalar-half":"flow-half";
                else {
                    courant=std::max({courant,half.flow.maxCourant,second.flow.maxCourant});
                    const double ts=c.temperatureAbsoluteTolerance+c.relativeTolerance*c.temperatureScale;
                    const double vs=c.velocityAbsoluteTolerance+c.relativeTolerance*c.velocityScale;
                    for (std::size_t i=0;i<mesh.cells.size();++i) {
                        // BE full-step local defect estimate: 2*(full - two halves).
                        error=std::max({error,2*std::abs(full.accepted->scalar[i]-second.accepted->scalar[i])/ts,
                            2*std::abs(full.accepted->flow.u[i]-second.accepted->flow.u[i])/vs,
                            2*std::abs(full.accepted->flow.v[i]-second.accepted->flow.v[i])/vs});
                    }
                    if (!std::isfinite(error)) reason="nonfinite-error";
                    else if (error>1) reason="time-error";
                    else if (courant>limits.maximumCourant) reason="courant";
                }
            }
        }
        out.attempts.push_back({previous.flow.time,dt,error,courant,trialFlow.velocityRelaxation,reason});
        if(observeAttempt)observeAttempt(out.attempts.back());
        if (reason=="accepted") {
            if (dt==limits.targetTime-previous.flow.time) {
                full.flow.time=limits.targetTime;full.accepted->flow.time=limits.targetTime;
            }
            // BE local defect scales as dt^2. Grow at most twice, and leave
            // a 10% error-budget margin. CFL predictor still caps the next trial.
            const double factor=c.estimateError && error>0?std::clamp(.9/std::sqrt(error),.5,2.):2.;
            const double suggested=std::clamp(dt*factor,c.limits.minimumStep,c.limits.maximumStep);
            double relaxation=trialFlow.velocityRelaxation;
            if(full.flow.history.size()<20)relaxation=std::min(fc.velocityRelaxation,1.2*relaxation);
            full.accepted->controller=ThermalControllerHistory2D{signature,suggested,relaxation};
            out.step=std::move(full);return out;
        }
        out.step=std::move(full);out.step.accepted.reset();
        if (out.step.flow.stopped) throw std::runtime_error("Thermal flow stopped; accepted state retained");
        if (retry>=limits.maximumRetries) return out;
        // Dampen failed nonlinear solves before shrinking physical time. This
        // changes the fixed-point path only, retaining every original residual
        // gate and the same BE equation, flux and physical setup.
        if (reason.starts_with("flow") && trialFlow.velocityRelaxation>.15) {
            trialFlow.velocityRelaxation=std::max(.15,.5*trialFlow.velocityRelaxation);
            continue;
        }
        auto smaller=reducedFlowTimeStep2D(dt,std::isfinite(courant)?courant:0,
            limits.targetTime-previous.flow.time,limits);
        if (!smaller) return out;
        dt=*smaller;
        if (!(previous.flow.time+dt>previous.flow.time) || (c.estimateError && !(previous.flow.time+.5*dt>previous.flow.time)))
            return out;
    }
}
}
