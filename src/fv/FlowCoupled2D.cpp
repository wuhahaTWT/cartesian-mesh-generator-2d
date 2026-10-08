#include "FlowCoupled2D.hpp"
#include "FlowEquation2D.hpp"
#include "FlowKrylov2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include <algorithm>
#include <chrono>
#include <limits>

namespace cartmesh2d::fv::solver_detail {
namespace {
using Clock=std::chrono::steady_clock;

detail::LinearPressureMethod2D pressureMethod(PressurePreconditioner2D method) {
    switch(method) {
    case PressurePreconditioner2D::Jacobi:return detail::LinearPressureMethod2D::Jacobi;
    case PressurePreconditioner2D::IncompleteCholesky0:return detail::LinearPressureMethod2D::IC0;
    case PressurePreconditioner2D::Aggregation:return detail::LinearPressureMethod2D::Aggregation;
    case PressurePreconditioner2D::SystemCholesky:return detail::LinearPressureMethod2D::SystemCholesky;
    }
    throw std::runtime_error("Invalid coupled pressure preconditioner");
}

struct CoupledMetrics {
    FlowIteration2D step;
    double merit=0,netFlux=0;
    std::size_t reverseOutletFaces=0;
};
} // namespace

FlowResult2D solveCoupledFlow2D(const FvMesh2D& mesh,const FlowControls2D& input,const Boundary& boundary,
    const std::function<void(const FlowIteration2D&)>& progress,const FlowInitialGuess2D* guess) {
    const auto start=Clock::now();
    ensure(input.convergence==FlowConvergence2D::Strict,"Coupled solving currently requires strict convergence");
    ensure(input.steadyAcceleration==SteadyAcceleration2D::None,"Coupled solving does not support Anderson acceleration");
    ensure(input.outletBackflow==OutletBackflow2D::Reject,"Coupled solving requires outlet backflow rejection");
    ensure(input.faceViscosity.empty() && input.manufacturedViscositySlope==0,
           "Coupled solving currently requires constant viscosity");
    for(const auto role:boundary.role)
        ensure(role!=Role::Opening && role!=Role::Farfield,"Coupled solving does not support dynamic pressure-opening/farfield boundaries");

    // The steady relaxation-defect Rhie--Chow equation has the unrelaxed
    // coefficient V/max(aU,aV) at its fixed point. Solve that equation directly;
    // the line search supplies nonlinear damping without changing the flux law.
    const bool consistent=input.coupling==FlowCoupling2D::SimpleConsistent;
    const double alphaU=consistent?input.velocityRelaxation:1.;
    const double alphaP=consistent?input.pressureRelaxation:1.;
    auto controls=input;controls.velocityRelaxation=1.;controls.pressureRelaxation=1.;
    auto lowOrder=controls;lowOrder.convection=ConvectionScheme2D::Upwind;
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    const double height=boundary.ymax-boundary.ymin;
    const double pressureScale=finite(controls.speed*controls.speed+controls.nu*controls.speed/height);
    ensure(pressureScale>0,"Coupled pressure reference underflow");
    const double consistencyTolerance=std::min(controls.tolerance,1e-8);
    Vec zeros(nf),zeroState(3*n);
    const std::vector<Vector2D> noSources;
    const auto& pressureBoundary=controls.scenario=="custom"?boundary.p:zeros;
    auto homogeneous=boundary;homogeneous.u=zeros;homogeneous.v=zeros;homogeneous.p=zeros;
    FlowEquation2D equation(mesh,boundary.fixedP);
    std::vector<std::pair<std::size_t,std::size_t>> connections;
    connections.reserve(nf);for(const auto& face:mesh.faces)if(face.neighbour)connections.emplace_back(face.owner,*face.neighbour);
    const detail::SparsePattern2D pattern(n,connections);
    System momentumU(pattern),momentumV(pattern),preconditionU(pattern),preconditionV(pattern),pressure(pattern);
    detail::LinearWorkspace2D pressureWorkspace(n);
    CoupledKrylov krylov(3*n);
    Vec state(3*n),ra(n,1.),candidateRa(n,1.),df(nf),rowScale(3*n),rhs(3*n),base(3*n),highResidual(3*n),lowResidual(3*n);
    Vec divergence(n),linearProduct(3*n),ru(n),rv(n),preU(n),preV(n),preP(n);
    Vec relaxationU(n),relaxationV(n),faceDefect(nf),noResponse(n);
    FlowResult2D result,scratch;
    for(auto* target:{&result,&scratch}){target->u.resize(n);target->v.resize(n);target->p.resize(n);target->flux.resize(nf);}
    result.domainHeight=height;initializeCaseSources(mesh,controls,result);
    initializeCaseVelocity(mesh,controls,boundary,result);
    if(guess) {
        result.u=guess->u;result.v=guess->v;result.p=guess->p;
        if(boundary.closed){const double gauge=result.p[0];for(auto& value:result.p)value=finite(value-gauge);}
    }
    for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];
        result.flux[id]=face.neighbour?interpolate(face,result.u)*face.areaVector.x+interpolate(face,result.v)*face.areaVector.y
            :boundary.role[id]==Role::Inlet?boundary.u[id]*face.areaVector.x+boundary.v[id]*face.areaVector.y
            :boundary.role[id]==Role::Outlet?result.u[face.owner]*face.areaVector.x+result.v[face.owner]*face.areaVector.y:0.;
    }
    if(guess && !guess->flux.empty()) {
        result.flux=guess->flux;
        for(std::size_t id=0;id<nf;++id) {
            const auto& face=mesh.faces[id];if(face.neighbour || boundary.role[id]==Role::Outlet)continue;
            const double expected=boundary.role[id]==Role::Inlet?boundary.u[id]*face.areaVector.x+boundary.v[id]*face.areaVector.y:0.;
            const double scale=controls.speed*std::hypot(face.areaVector.x,face.areaVector.y)+std::abs(expected);
            ensure(std::abs(result.flux[id]-expected)<=64*std::numeric_limits<double>::epsilon()*scale,
                   "Initial face flux violates prescribed boundary flux");
            result.flux[id]=expected;
        }
    }
    for(std::size_t i=0;i<n;++i){state[i]=result.u[i]/controls.speed;state[n+i]=result.v[i]/controls.speed;state[2*n+i]=result.p[i]/pressureScale;}
    double domainArea=0,totalLength=0;
    for(const auto& cell:mesh.cells)domainArea+=cell.area;
    for(const auto& face:mesh.faces)totalLength+=std::hypot(face.areaVector.x,face.areaVector.y);
    std::vector<std::size_t> monitorGroup(nf);Vec monitorLengths;
    prepareMonitors(mesh,controls,boundary,result,monitorGroup,monitorLengths);

    const auto unpack=[&](const Vec& x) {
        for(std::size_t i=0;i<n;++i){scratch.u[i]=finite(controls.speed*x[i]);scratch.v[i]=finite(controls.speed*x[n+i]);scratch.p[i]=finite(pressureScale*x[2*n+i]);}
    };
    const auto forceGradient=[&](const Vec& p,const Vec& pb) {
        return equation.pressureForce(p,pb,boundary.fixedP);
    };
    const auto flux=[&](const Boundary& b,const Vec& response,const Vec& pb) {
        const auto reconstruction=equation.reconstruct(scratch,b,pb);
        return equation.steadyFlux(scratch,controls,b,reconstruction,response,df);
    };
    // With frozen advective flux, the upwind operator is affine in U,V,p.
    // Physical sources and prescribed traces belong only to its affine part.
    const auto residual=[&](const Vec& x,Vec& out,const Boundary& b,const FlowControls2D& c,
                            const Vec& advectiveFlux,const Vec& response,bool physical) {
        unpack(x);
        equation.frozenResidual(momentumU,momentumV,scratch,c,b,physical?pressureBoundary:zeros,
            advectiveFlux,response,physical?result.sourceIntegrals:noSources,out,df,false);
        if(boundary.closed)out[2*n]=x[2*n];
    };
    const auto feasible=[&](const Vec& q) {
        for(std::size_t id=0;id<nf;++id)if(!mesh.faces[id].neighbour && boundary.role[id]==Role::Outlet && q[id]<0)return false;
        return true;
    };
    ensure(feasible(result.flux),"Coupled initial face flux violates outlet backflow rejection");
    const auto evaluate=[&](const Vec& x,const Vec& q) {
        ensure(feasible(q),"Coupled trial face flux violates outlet backflow rejection");
        residual(x,linearProduct,boundary,controls,q,ra,true);
        for(std::size_t i=0;i<n;++i) {
            const double diagonal=std::max(momentumU.diag[i],momentumV.diag[i]);
            ensure(diagonal>0 && std::isfinite(diagonal),"Coupled momentum diagonal is not positive");
            candidateRa[i]=mesh.cells[i].area/diagonal;
        }
        const auto reconstructed=flux(boundary,candidateRa,pressureBoundary);
        std::fill(divergence.begin(),divergence.end(),0.);CoupledMetrics metrics;double incoming=0;long double energy=0;
        for(std::size_t id=0;id<nf;++id) {
            const auto& face=mesh.faces[id];const double length=std::hypot(face.areaVector.x,face.areaVector.y);
            divergence[face.owner]+=q[id];if(face.neighbour)divergence[*face.neighbour]-=q[id];
            else{metrics.netFlux+=q[id];incoming+=std::max(0.,-q[id]);if(boundary.role[id]==Role::Outlet && q[id]<0)++metrics.reverseOutletFaces;}
            const double defect=finite((q[id]-reconstructed[id])/(controls.speed*length));
            metrics.step.fluxConsistency=std::max(metrics.step.fluxConsistency,std::abs(defect));
            energy+=length/totalLength*defect*defect;
        }
        for(std::size_t i=0;i<n;++i) {
            const double scale=(momentumU.diag[i]+momentumV.diag[i])*controls.speed;
            const double xResidual=linearProduct[i]/scale,yResidual=linearProduct[n+i]/scale;
            const double magnitude=std::hypot(xResidual,yResidual);
            if(magnitude>metrics.step.momentumResidual){metrics.step.momentumResidual=magnitude;metrics.step.momentumWorstCell=i;metrics.step.momentumResidualX=xResidual;metrics.step.momentumResidualY=yResidual;}
            const double continuity=divergence[i]/(controls.speed*std::sqrt(mesh.cells[i].area));
            metrics.step.continuity=std::max(metrics.step.continuity,std::abs(continuity));
            energy+=mesh.cells[i].area/domainArea*(xResidual*xResidual+yResidual*yResidual+continuity*continuity);
        }
        const double throughput=boundary.closed?controls.speed*height:incoming;
        ensure(throughput>0,"Coupled flow has no positive reference throughput");
        metrics.step.globalRelativeImbalance=finite(std::abs(metrics.netFlux)/throughput);
        metrics.merit=finite(std::sqrt(static_cast<double>(energy)));
        return metrics;
    };
    const auto converged=[&](const FlowIteration2D& step) {
        return detail::strictFlowConverged2D(step.iteration,step.momentumResidual,step.velocityChange,
            step.pressureChange,step.continuity,step.globalRelativeImbalance,controls.tolerance,true)
            && step.fluxConsistency<consistencyTolerance;
    };
    auto accepted=evaluate(state,result.flux);
    accepted.step.strictLinearStep=false;
    bool certifyNext=false;
    for(std::size_t iteration=1;iteration<=controls.maxIterations;++iteration) {
        const Vec old=state,oldFlux=result.flux;
        try {
            const auto before=evaluate(old,oldFlux);ra=candidateRa;
            residual(old,lowResidual,boundary,lowOrder,oldFlux,ra,true);
            // A diagonal pressure correction discards precisely the off-row
            // viscous response present in the predictor. Solve the full frozen
            // response instead. The added diagonal and its old-field source
            // cancel at the fixed point, including unequal slip diagonals.
            for(std::size_t i=0;i<n;++i) {
                const double common=std::max(momentumU.diag[i],momentumV.diag[i]);
                relaxationU[i]=consistent?common/alphaU-momentumU.diag[i]:0.;
                relaxationV[i]=consistent?common/alphaU-momentumV.diag[i]:0.;
                ra[i]*=alphaU;
            }
            if(consistent) {
                unpack(old);
                const auto interpolated=flux(boundary,noResponse,pressureBoundary);
                for(std::size_t id=0;id<nf;++id) {
                    const auto& face=mesh.faces[id];
                    faceDefect[id]=(face.neighbour || boundary.role[id]==Role::Outlet)
                        ?(1-alphaU)*(oldFlux[id]-interpolated[id]):0.;
                }
            }
            preconditionU.diag=momentumU.diag;preconditionU.off=momentumU.off;
            preconditionV.diag=momentumV.diag;preconditionV.off=momentumV.off;
            for(std::size_t i=0;i<n;++i){preconditionU.diag[i]+=relaxationU[i];preconditionV.diag[i]+=relaxationV[i];}
            preconditionU.factorILU0();preconditionV.factorILU0();
            for(std::size_t id=0;id<nf;++id)df[id]=interpolate(mesh.faces[id],ra)*mesh.faces[id].transmissibility;
            assemblePressureBlock2D(pressure,mesh,boundary,df);
            for(std::size_t i=0;i<n;++i){rowScale[i]=preconditionU.diag[i]*controls.speed;rowScale[n+i]=preconditionV.diag[i]*controls.speed;rowScale[2*n+i]=pressure.diag[i]*pressureScale;}
            if(boundary.closed)rowScale[2*n]=1.;
            residual(old,highResidual,boundary,controls,oldFlux,ra,true);
            if(consistent)residual(old,lowResidual,boundary,lowOrder,oldFlux,ra,true);
            residual(zeroState,base,boundary,lowOrder,oldFlux,ra,true);
            for(std::size_t i=0;i<3*n;++i)rhs[i]=finite((-base[i]-highResidual[i]+lowResidual[i])/rowScale[i]);
            for(std::size_t i=0;i<n;++i) {
                rhs[i]+=relaxationU[i]*controls.speed*old[i]/rowScale[i];
                rhs[n+i]+=relaxationV[i]*controls.speed*old[n+i]/rowScale[n+i];
            }
            for(std::size_t id=0;id<nf;++id) {
                const auto& face=mesh.faces[id];
                rhs[2*n+face.owner]-=faceDefect[id]/rowScale[2*n+face.owner];
                if(face.neighbour)rhs[2*n+*face.neighbour]+=faceDefect[id]/rowScale[2*n+*face.neighbour];
            }
            if(boundary.closed)rhs[2*n]=0.;
            const auto apply=[&](const Vec& x,Vec& y){
                residual(x,y,homogeneous,lowOrder,oldFlux,ra,false);
                for(std::size_t i=0;i<n;++i){y[i]+=relaxationU[i]*controls.speed*x[i];y[n+i]+=relaxationV[i]*controls.speed*x[n+i];}
                for(std::size_t i=0;i<3*n;++i)y[i]/=rowScale[i];
            };
            const auto precondition=[&](const Vec& x,Vec& y) {
                for(std::size_t i=0;i<n;++i){ru[i]=x[i]*rowScale[i];rv[i]=x[n+i]*rowScale[n+i];}
                preconditionU.preconditionILU0(ru,preU);preconditionV.preconditionILU0(rv,preV);
                scratch.u=preU;scratch.v=preV;std::fill(scratch.p.begin(),scratch.p.end(),0.);
                const auto q=flux(homogeneous,ra,zeros);
                for(std::size_t i=0;i<n;++i)pressure.rhs[i]=x[2*n+i]*rowScale[2*n+i];
                for(std::size_t id=0;id<nf;++id){const auto& face=mesh.faces[id];pressure.rhs[face.owner]-=q[id];if(face.neighbour)pressure.rhs[*face.neighbour]+=q[id];}
                if(boundary.closed)pressure.rhs[0]=0.;
                std::fill(preP.begin(),preP.end(),0.);
                const auto pressureStart=controls.profile?Clock::now():Clock::time_point{};
                const auto builds=pressure.ic0Builds(),reuses=pressure.ic0Reuses();
                const auto cb=pressure.choleskyBuilds(),cr=pressure.choleskyRefactors(),cu=pressure.choleskyReuses();
                const auto hb=pressure.hierarchyBuilds(),hr=pressure.hierarchyReuses(),hf=pressure.hierarchyRefreshes();
                const auto iterations=pressure.solvePressure(preP,pressureWorkspace,pressureMethod(controls.pressurePreconditioner),1e-3);
                auto& performance=result.performance;++performance.pressureSolves;performance.pressureIterations+=iterations;
                performance.maxPressureIterations=std::max(performance.maxPressureIterations,iterations);
                performance.pressureFactorizations+=pressure.ic0Builds()-builds;performance.pressureFactorReuses+=pressure.ic0Reuses()-reuses;
                performance.pressureCholeskyBuilds+=pressure.choleskyBuilds()-cb;performance.pressureCholeskyRefactors+=pressure.choleskyRefactors()-cr;performance.pressureCholeskyReuses+=pressure.choleskyReuses()-cu;
                performance.pressureHierarchyBuilds+=pressure.hierarchyBuilds()-hb;performance.pressureHierarchyReuses+=pressure.hierarchyReuses()-hr;performance.pressureHierarchyRefreshes+=pressure.hierarchyRefreshes()-hf;
                performance.maxPressureHierarchyLevels=std::max(performance.maxPressureHierarchyLevels,pressure.hierarchyLevels());performance.maxPressureCoarseCells=std::max(performance.maxPressureCoarseCells,pressure.hierarchyCoarseCells());
                if(controls.profile)performance.pressureLinearSolveSeconds+=std::chrono::duration<double>(Clock::now()-pressureStart).count();
                const auto g=forceGradient(preP,zeros);
                for(std::size_t i=0;i<n;++i){y[i]=(preU[i]-ra[i]*g[i].x)/controls.speed;y[n+i]=(preV[i]-ra[i]*g[i].y)/controls.speed;y[2*n+i]=preP[i]/pressureScale;}
            };
            bool strict=!controls.adaptiveLinear || certifyNext ||
                (before.step.momentumResidual<100*controls.tolerance && before.step.fluxConsistency<100*consistencyTolerance && before.step.continuity<1e-8);
            certifyNext=false;
            bool found=false;Vec proposed=old,trial(3*n),trialFlux(nf);CoupledMetrics candidate;
            double weight=1.,relativeTolerance=1e-13;
            for(unsigned attempt=0;attempt<2 && !found;++attempt) {
                relativeTolerance=strict?1e-13:std::clamp(.01*before.merit,1e-13,1e-3);
                const auto linearStart=controls.profile?Clock::now():Clock::time_point{};
                const auto count=result.performance.coupledIterations;++result.performance.coupledSolves;
                const auto finishLinear=[&]{result.performance.maxCoupledIterations=std::max(result.performance.maxCoupledIterations,result.performance.coupledIterations-count);if(controls.profile)result.performance.coupledLinearSolveSeconds+=std::chrono::duration<double>(Clock::now()-linearStart).count();};
                proposed=old;
                try{krylov.solve(apply,precondition,proposed,rhs,relativeTolerance,strict?.01*controls.tolerance:std::numeric_limits<double>::infinity(),result.performance);}
                catch(...){finishLinear();throw;}finishLinear();
                unpack(proposed);auto proposedFlux=flux(boundary,ra,pressureBoundary);
                for(std::size_t id=0;id<nf;++id)proposedFlux[id]+=faceDefect[id];
                // The pressure in the block solve is p_old + pc. Only the
                // stored pressure is relaxed; U and the conservative face flux
                // share the complete, unrelaxed correction response.
                if(consistent && alphaP!=1.)for(std::size_t i=0;i<n;++i)proposed[2*n+i]=old[2*n+i]+alphaP*(proposed[2*n+i]-old[2*n+i]);
                weight=1.;
                for(std::size_t id=0;id<nf;++id)if(!mesh.faces[id].neighbour && boundary.role[id]==Role::Outlet && proposedFlux[id]<0)
                    weight=std::min(weight,.95*oldFlux[id]/(oldFlux[id]-proposedFlux[id]));
                ensure(weight>0 && std::isfinite(weight),"Coupled update has no feasible direction at the outlet");
                for(unsigned backtrack=0;backtrack<24;++backtrack) {
                    ++result.performance.lineSearchTrials;
                    for(std::size_t i=0;i<3*n;++i)trial[i]=finite(old[i]+weight*(proposed[i]-old[i]));
                    for(std::size_t id=0;id<nf;++id)trialFlux[id]=finite(oldFlux[id]+weight*(proposedFlux[id]-oldFlux[id]));
                    candidate=evaluate(trial,trialFlux);
                    candidate.step.iteration=iteration;candidate.step.velocityChange=0;candidate.step.pressureChange=0;
                    for(std::size_t i=0;i<n;++i){candidate.step.velocityChange=std::max(candidate.step.velocityChange,std::hypot(trial[i]-old[i],trial[n+i]-old[n+i]));candidate.step.pressureChange=std::max(candidate.step.pressureChange,std::abs(trial[2*n+i]-old[2*n+i]));}
                    // The merit chooses a step only. It never replaces the
                    // maximum-cell, field-change, continuity or strict gates.
                    auto roundoffCandidate=candidate.step;roundoffCandidate.iteration=std::max<std::size_t>(10,iteration);
                    if(candidate.merit<=before.merit*(1.-1e-4*weight) || converged(roundoffCandidate)){found=true;break;}
                    ++result.performance.lineSearchRejected;weight*=.5;
                }
                if(!found && !strict){strict=true;continue;}
                if(!found)throw std::runtime_error("Coupled line search found no residual-decreasing feasible step");
            }
            ensure(found,"Coupled line search failed");
            unpack(trial);
            if(!result.monitorNames.empty()) {
                scratch.flux=trialFlux;scratch.monitorNames=result.monitorNames;
                const auto reconstructed=equation.reconstruct(scratch,controls,boundary,pressureBoundary,trialFlux);
                candidate.step.monitors=physicalMonitors(mesh,controls,boundary,scratch,monitorGroup,monitorLengths,
                    pressureScale,reconstructed.pressureFaces,reconstructed.gu,reconstructed.gv,reconstructed.stress);
            }
            // Commit only a fully evaluated candidate. Scratch fields used by
            // Krylov, preconditioning and rejected trials never own the result.
            state=trial;result.flux=trialFlux;
            result.u=scratch.u;result.v=scratch.v;result.p=scratch.p;
            accepted=candidate;accepted.step.strictLinearStep=strict;
            accepted.step.linearRelativeTolerance=relativeTolerance;
            accepted.step.coupledLinearResidual=krylov.residualNorm;accepted.step.lineSearchWeight=weight;
            result.globalImbalance=accepted.netFlux;result.globalRelativeImbalance=accepted.step.globalRelativeImbalance;
            if(iteration<=10)result.convergenceReference=std::max(result.convergenceReference,accepted.step.momentumResidual);
            result.history.push_back(accepted.step);
        } catch(const std::runtime_error& error) {
            result.failureReason=error.what();
            if(result.history.empty())result.history.push_back(accepted.step);
            break;
        }
        const auto& step=result.history.back();
        if(progress && (iteration==1 || iteration%10==0))progress(step);
        if(converged(step)) {
            if(step.strictLinearStep){result.converged=true;break;}
            certifyNext=true;
        }
        if(controls.stopRequested && controls.stopRequested()){result.stopped=true;break;}
    }
    // Independently refresh rAU from the retained flux before any success is
    // returned. The divergence includes the gauge cell, not its pinned row.
    const auto final=evaluate(state,result.flux);
    result.globalImbalance=final.netFlux;result.globalRelativeImbalance=final.step.globalRelativeImbalance;
    if(result.converged) {
        auto finalStep=final.step;const auto& last=result.history.back();
        finalStep.iteration=last.iteration;finalStep.velocityChange=last.velocityChange;finalStep.pressureChange=last.pressureChange;
        if(!converged(finalStep) || final.reverseOutletFaces!=0){result.converged=false;result.failureReason="Coupled final residual/flux consistency certification failed";}
    }
    const auto reconstructed=equation.reconstruct(result,controls,boundary,pressureBoundary,result.flux);
    postprocessForces(mesh,controls,boundary,result,pressureBoundary,
        reconstructed.gp,reconstructed.gu,reconstructed.gv,reconstructed.stress);
    if(controls.profile)result.performance.solveSeconds=std::chrono::duration<double>(Clock::now()-start).count();
    return result;
}
} // namespace cartmesh2d::fv::solver_detail
