#include "cartmesh2d/fv/detail/FlowNewtonKrylov2D.hpp"
#include "cartmesh2d/fv/detail/NewtonKrylov2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include <chrono>

namespace cartmesh2d::fv::detail {
namespace {
using Vec=LinearVector2D;
struct BudgetReached {};
struct StopRequested {};

void addWork(FlowPerformance2D& total,const FlowPerformance2D& part) {
    total.momentumSolves+=part.momentumSolves;
    total.momentumIterations+=part.momentumIterations;
    total.maxMomentumIterations=std::max(total.maxMomentumIterations,part.maxMomentumIterations);
    total.pressureSolves+=part.pressureSolves;
    total.pressureCorrectionPassesSkipped+=part.pressureCorrectionPassesSkipped;
    total.pressureFactorizations+=part.pressureFactorizations;
    total.pressureFactorReuses+=part.pressureFactorReuses;
    total.pressureCholeskyBuilds+=part.pressureCholeskyBuilds;
    total.pressureCholeskyRefactors+=part.pressureCholeskyRefactors;
    total.pressureCholeskyReuses+=part.pressureCholeskyReuses;
    total.pressureHierarchyBuilds+=part.pressureHierarchyBuilds;
    total.pressureHierarchyReuses+=part.pressureHierarchyReuses;
    total.pressureHierarchyRefreshes+=part.pressureHierarchyRefreshes;
    total.maxPressureHierarchyLevels=std::max(total.maxPressureHierarchyLevels,part.maxPressureHierarchyLevels);
    total.maxPressureCoarseCells=std::max(total.maxPressureCoarseCells,part.maxPressureCoarseCells);
    total.pressureIterations+=part.pressureIterations;
    total.maxPressureIterations=std::max(total.maxPressureIterations,part.maxPressureIterations);
    total.momentumLinearSolveSeconds+=part.momentumLinearSolveSeconds;
    total.pressureLinearSolveSeconds+=part.pressureLinearSolveSeconds;
}
}

FlowResult2D solveNewtonKrylovFlow2D(const FvMesh2D& mesh,const FlowControls2D& control,
    const FlowInitialGuess2D* initialGuess,const std::function<void(const FlowIteration2D&)>& progress) {
    linearEnsure(control.maxIterations>0,"Invalid flow iteration budget");
    linearEnsure(std::isfinite(control.tolerance) && control.tolerance>0,"Invalid flow tolerance");
    linearEnsure(!control.adaptiveLinear && control.convergence==FlowConvergence2D::Strict,
        "Newton-Krylov requires strict linear solving and strict convergence");
    using Clock=std::chrono::steady_clock;const auto start=Clock::now();
    auto one=control;one.steadyAcceleration=SteadyAcceleration2D::None;
    one.maxIterations=1;one.stopRequested={};
    // Existing strict RHS-relative accuracy is 1e-11. The original momentum
    // row target is .01*tolerance*Uref*alphaU; ensure that it is also no larger
    // than 1e-11*Uref*alphaU when probing a Jacobian. This tightens only the
    // inner map evaluation, never the user's final convergence gates.
    constexpr double mapError=1e-11;
    one.tolerance=std::min(control.tolerance,100*mapError);
    FlowPerformance2D work;
    auto result=initialGuess ? solveIncompressibleFromGuess2D(mesh,one,*initialGuess)
                             : solveIncompressible2D(mesh,one);
    work.coupledEvaluations=1;addWork(work,result.performance);
    std::vector<FlowIteration2D> history=result.history;
    const auto finish=[&](FlowResult2D value) {
        value.history=history;
        value.history.back().iteration=work.coupledEvaluations;
        work.solveSeconds=std::chrono::duration<double>(Clock::now()-start).count();
        value.performance=work;
        return value;
    };
    if(progress)progress(history.back());
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    double area=0,length=0;
    for(const auto& cell:mesh.cells)area+=cell.area;
    for(const auto& face:mesh.faces)length+=std::hypot(face.areaVector.x,face.areaVector.y);
    const double pressureScale=control.speed*control.speed+control.nu*control.speed/result.domainHeight;
    Vec scale(3*n+nf);
    for(std::size_t i=0;i<n;++i) {
        const double weight=std::sqrt(mesh.cells[i].area/area);
        scale[i]=scale[n+i]=weight/control.speed;scale[2*n+i]=weight/pressureScale;
    }
    for(std::size_t i=0;i<nf;++i) {
        const auto& face=mesh.faces[i];const double size=std::hypot(face.areaVector.x,face.areaVector.y);
        scale[3*n+i]=std::sqrt(size/length)/(control.speed*size);
    }
    const auto pack=[&](const FlowResult2D& value) {
        Vec state;state.reserve(scale.size());
        for(const auto* field:{&value.u,&value.v,&value.p,&value.flux})state.insert(state.end(),field->begin(),field->end());
        for(std::size_t i=0;i<state.size();++i)state[i]=linearFinite(state[i]*scale[i]);
        return state;
    };
    const auto unpack=[&](const Vec& state) {
        FlowInitialGuess2D value;
        value.u.resize(n);value.v.resize(n);value.p.resize(n);value.flux.resize(nf);
        std::size_t i=0;
        for(auto* field:{&value.u,&value.v,&value.p,&value.flux})for(auto& x:*field){x=linearFinite(state[i]/scale[i]);++i;}
        return value;
    };
    const auto checkBudget=[&] {
        if(control.stopRequested && control.stopRequested())throw StopRequested{};
        if(work.coupledEvaluations>=control.maxIterations)throw BudgetReached{};
    };
    struct Evaluation {Vec residual;FlowResult2D mapped;};
    const auto evaluate=[&](const Vec& state,bool alreadyChecked=false,bool certificate=false) {
        if(!alreadyChecked)checkBudget();
        ++work.coupledEvaluations;
        try {
            auto settings=one;if(certificate)settings.tolerance=control.tolerance;
            auto mapped=solveIncompressibleFromGuess2D(mesh,settings,unpack(state));addWork(work,mapped.performance);
            auto residual=pack(mapped);
            for(std::size_t i=0;i<residual.size();++i)residual[i]=linearFinite(state[i]-residual[i]);
            return Evaluation{std::move(residual),std::move(mapped)};
        } catch(const std::runtime_error& error) {
            ++work.coupledFailedEvaluations;work.coupledLastFailure=error.what();throw;
        }
    };
    auto state=pack(result);
    try {
        auto base=evaluate(state);result=base.mapped;
        auto record=[&] {
            auto step=result.history.back();step.iteration=work.coupledEvaluations;
            history.push_back(step);if(progress)progress(step);
        };
        record();
        for(;;) {
            checkBudget();const auto& step=base.mapped.history.back();
            // Certify the actual returned field with a fresh, unaccelerated
            // map at the user's original strict linear settings. The shared
            // momentum, full update and conservation gates all remain required.
            // Ten chained SIMPLE updates instead test contraction of G, which
            // is unnecessary for Newton's root of x-G(x) and can amplify noise.
            if(strictFlowResidualsAccepted2D(step,control.tolerance)) {
                std::optional<Evaluation> certificate;
                try {certificate=evaluate(pack(result),false,true);}
                catch(const std::runtime_error&) {
                    // Already counted and diagnosed by evaluate. A failed
                    // certificate must not replace the accepted coupled pair.
                }
                if(certificate && certificate->mapped.history.back().strictLinearStep &&
                   strictFlowResidualsAccepted2D(certificate->mapped.history.back(),control.tolerance)) {
                    result=std::move(certificate->mapped);result.converged=true;
                    record();return finish(std::move(result));
                }
            }
            const double merit=linearNorm(base.residual);
            const double epsilon=std::sqrt(mapError)*std::sqrt(1+linearNorm(state));
            const auto jacobian=[&](const Vec& direction) {
                auto trial=state;
                for(std::size_t i=0;i<trial.size();++i)trial[i]+=epsilon*direction[i];
                auto residual=evaluate(trial).residual;
                for(std::size_t i=0;i<residual.size();++i)residual[i]=linearFinite((residual[i]-base.residual[i])/epsilon);
                return residual;
            };
            auto rhs=base.residual;for(auto& value:rhs)value=-value;
            const auto direction=newtonKrylovDirection2D(jacobian,rhs);
            if(!direction){result.nonlinearStagnated=true;break;}
            bool accepted=false;double weight=1;
            // The dimensionless Armijo decrease compares the actual scaled
            // u,v,p,face-flux fixed-point residual. Its 1e-4 slope and twelve
            // halvings bound globalization work; neither replaces a CFD gate.
            for(int backtrack=0;backtrack<12;++backtrack,weight*=.5) {
                checkBudget();
                auto trial=state;
                for(std::size_t i=0;i<trial.size();++i)trial[i]+=weight*(*direction)[i];
                std::optional<Evaluation> value;
                try {value=evaluate(trial,true);}
                catch(const std::runtime_error&) {
                    // A bad Newton candidate is not the accepted result.
                }
                ++work.accelerationCandidates;
                if(value && linearNorm(value->residual)<(1-1e-4*weight)*merit) {
                    state=std::move(trial);base=std::move(*value);result=base.mapped;
                    ++work.accelerationAccepted;++work.coupledOuterIterations;
                    accepted=true;break;
                }
                ++work.accelerationRejected;
            }
            if(!accepted){result.nonlinearStagnated=true;break;}
            record();
        }
    } catch(const BudgetReached&) {
        // Return only the last fully evaluated, accepted field/flux pair.
    } catch(const StopRequested&) {result.stopped=true;}
    return finish(std::move(result));
}
} // namespace cartmesh2d::fv::detail
