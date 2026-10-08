#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/TaylorGreen2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/fv/detail/FlowMaterial2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include "cartmesh2d/fv/detail/Anderson2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

#include "FlowSolverDetail2D.hpp"
#include "FlowPressure2D.hpp"
#include "FlowCoupled2D.hpp"
namespace cartmesh2d::fv {
using namespace solver_detail;
static FlowResult2D solveFlow(
    const FvMesh2D& m, const FlowControls2D& input,
    const std::function<void(const FlowIteration2D&)>& progress,
    const FlowState2D* previous, double timeStep,
    const detail::MaterialUpdate2D& material = {}, const FlowInitialGuess2D* guess = nullptr) {
    auto c=input;
    using Clock = std::chrono::steady_clock;
    const auto solveStart = c.profile ? Clock::now() : Clock::time_point{};
    validateFvMesh2D(m);
    ensure((c.scenario == "external" || c.scenario == "channel" || c.scenario == "duct" || c.scenario == "custom" || c.scenario == "cavity" || c.scenario == "manufactured" || c.scenario == "taylor-green" || c.scenario == "counterflow" || c.scenario == "flatplate") &&
               std::isfinite(c.nu) && c.nu > 0 && std::isfinite(c.speed) && c.speed > 0 &&
               std::isfinite(c.tolerance) && c.tolerance > 0 && c.maxIterations > 0,
           "Invalid flow controls");
    ensure(std::isfinite(c.flatPlateLeadingEdge) &&
           (c.scenario=="flatplate" || c.flatPlateLeadingEdge==0),
           "Leading edge is only supported by the flat plate case");
    ensure(c.flatPlateTop==FlatPlateTop2D::PressureFarfield ||
           (c.scenario=="flatplate" && c.flatPlateTop==FlatPlateTop2D::Symmetry),
           "Invalid flat plate upper boundary");
    ensure(!previous || c.scenario!="flatplate","Transient flat plate is not implemented");
    ensure(std::isfinite(c.manufacturedPressureSlope) &&
               (c.scenario == "manufactured" || c.manufacturedPressureSlope == 0),
           "Manufactured pressure slope is only valid for the verification case");
    ensure(c.velocityRelaxation > 0 && c.velocityRelaxation <= 1 &&
               c.pressureRelaxation > 0 && c.pressureRelaxation <= 1,
           "Invalid SIMPLE relaxation");
    ensure(c.coupling==FlowCoupling2D::Simple ||
           ((c.coupling==FlowCoupling2D::SimpleC || c.coupling==FlowCoupling2D::Coupled) && !previous && !material),
           "SIMPLEC and coupled solving require steady constant-property laminar flow");
    ensure(c.pressureCorrectionPasses>=1 && c.pressureCorrectionPasses<=4,
           "Pressure corrections must be in [1,4]");
    ensure(!material || c.pressureCorrectionPasses==4,
           "Variable pressure corrections support constant-property laminar flow only");
    ensure(c.convergence==FlowConvergence2D::Strict ||
           (c.convergence==FlowConvergence2D::Engineering && !previous && !material),
           "Engineering convergence requires steady laminar flow");
    ensure(!c.adaptiveLinear || !material,"Adaptive linear solving does not support material coupling");
    ensure(c.steadyAcceleration==SteadyAcceleration2D::None ||
           (c.steadyAcceleration==SteadyAcceleration2D::Anderson && !previous && !material),
           "Anderson acceleration requires steady laminar flow");
    ensure(c.andersonHistory>=1 && c.andersonHistory<=32 && c.andersonStart>=1,
           "Invalid Anderson history/start controls");
    ensure(c.pressurePreconditioner == PressurePreconditioner2D::Jacobi ||
               c.pressurePreconditioner == PressurePreconditioner2D::IncompleteCholesky0 ||
               c.pressurePreconditioner == PressurePreconditioner2D::Aggregation ||
               c.pressurePreconditioner == PressurePreconditioner2D::SystemCholesky,
           "Invalid pressure preconditioner");
    ensure(c.pressurePreconditioner!=PressurePreconditioner2D::SystemCholesky || detail::systemCholeskyAvailable2D(),
           "System sparse Cholesky is available only on macOS");
    ensure(c.pressurePreconditioner!=PressurePreconditioner2D::SystemCholesky || !material,
           "System sparse Cholesky currently supports constant-property laminar flow only");
    ensure(c.convection == ConvectionScheme2D::Upwind ||
               c.convection == ConvectionScheme2D::LimitedLinearUpwind ||
               c.convection == ConvectionScheme2D::FaceLimitedLinearUpwind,
           "Invalid convection scheme");
    ensure(c.momentumInertia == 1 || (c.momentumInertia == 0 && c.scenario == "custom" && !previous && !material),
           "Momentum inertia must be 1, or 0 for steady custom laminar Stokes only");
    ensure(c.viscousStress == ViscousStress2D::Symmetric ||
               c.viscousStress == ViscousStress2D::Laplacian,
           "Invalid viscous stress form");
    ensure(c.outletBackflow == OutletBackflow2D::Reject || c.outletBackflow == OutletBackflow2D::NormalInlet,
           "Invalid outlet backflow model");
    validateViscosity(m,c);
    auto b = boundaries(m, c);
    const auto n = m.cells.size();
    const auto nf = m.faces.size();
    if(guess) {
        ensure(!previous && !material,"Steady initial guess cannot initialize time or material coupling");
        ensure(guess->u.size()==n && guess->v.size()==n && guess->p.size()==n,"Initial guess size differs from mesh");
        ensure(guess->flux.empty() || guess->flux.size()==nf,"Initial face flux size differs from mesh");
        for(const auto* field:{&guess->u,&guess->v,&guess->p,&guess->flux})for(double value:*field)finite(value);
    }
    if (previous) {
        ensure(std::isfinite(timeStep) && timeStep>0 && std::isfinite(previous->time) && previous->time>=0 &&
                   std::isfinite(previous->time+timeStep) && previous->time+timeStep>previous->time,
               "Invalid physical time step");
        ensure(c.scenario!="manufactured" && c.scenario!="counterflow", "Transient forced manufactured case not implemented");
        ensure(previous->u.size()==n && previous->v.size()==n && previous->p.size()==n && previous->flux.size()==nf,
               "Transient state size differs from mesh");
        for (const auto* field : {&previous->u,&previous->v,&previous->p,&previous->flux})
            for (double value : *field) finite(value);
    } else ensure(c.scenario!="taylor-green", "Taylor-Green requires physical time stepping");

    // A single pressure gauge is only valid for one connected fluid region.
    // Reject disconnected cavities instead of silently selecting arbitrary gauges.
    std::vector<bool> visited(n, false);
    std::vector<std::size_t> pending{0};
    visited[0] = true;
    for (std::size_t k = 0; k < pending.size(); ++k) {
        const auto cell = pending[k];
        for (auto id : m.cells[cell].faces) {
            const auto& f = m.faces[id];
            if (!f.neighbour) {
                continue;
            }
            const auto other = f.owner == cell ? *f.neighbour : f.owner;
            if (!visited[other]) {
                visited[other] = true;
                pending.push_back(other);
            }
        }
    }
    ensure(pending.size() == n,
           "Flow requires one connected fluid region with an unambiguous pressure reference");

    if(c.coupling==FlowCoupling2D::Coupled) {
        auto coupled=solveCoupledFlow2D(m,c,b,progress,guess);
        if(c.profile)coupled.performance.solveSeconds=std::chrono::duration<double>(Clock::now()-solveStart).count();
        return coupled;
    }

    std::vector<std::pair<std::size_t, std::size_t>> connections;
    connections.reserve(nf);
    for (const auto& face : m.faces)
        if (face.neighbour) connections.emplace_back(face.owner, *face.neighbour);
    const detail::SparsePattern2D pattern(n, connections);
    System au(pattern), av(pattern), ap(pattern), checkU(pattern), checkV(pattern);
    detail::LinearWorkspace2D workspace(n);
    Vec mu(n), mv(n);
    FlowResult2D r;
    // Global RHS-relative stopping can mask a tiny cut-cell residual when
    // large far-field cells carry the time term. Also require each row's
    // residual/diagonal in velocity units to be <=1% of the nonlinear target.
    // The existing global linear and nonlinear gates both remain in force.
    // Apply this to steady solves too: an RHS-relative linear tolerance can
    // otherwise return zero updates above a stricter nonlinear stop forever.
    const double momentumScaledStop=finite(.01*c.tolerance*c.speed*c.velocityRelaxation);
    ensure(momentumScaledStop>0,"Flow linear residual scale underflow");
    // Diagnostic normalization only; does not change pressure stopping.
    double shortestFace=std::numeric_limits<double>::infinity();
    for(const auto& face:m.faces)shortestFace=std::min(shortestFace,std::hypot(face.areaVector.x,face.areaVector.y));
    // Profiling observes the same solves and stopping rules, including zero-step
    // solves. Timing includes each linear solver's setup, but not assembly.
    detail::AdaptiveLinear2D linearPolicy;
    bool certifyNext=false, strictLinearIteration=true;
    double linearRelativeTolerance=1e-11, momentumRowStop=momentumScaledStop;
    double predictorResidual=0,pressureLinearResidual=0;
    std::size_t predictorWorstCell=0;
    auto linearSolve = [&](const System& system, Vec& field, bool pressure) {
        const auto start = c.profile ? Clock::now() : Clock::time_point{};
        const auto oldBuilds = system.ic0Builds(), oldReuses = system.ic0Reuses();
        const auto oldHierarchies=system.hierarchyBuilds(), oldHierarchyReuses=system.hierarchyReuses();
        const auto oldHierarchyRefreshes=system.hierarchyRefreshes();
        const auto oldCholeskyBuilds=system.choleskyBuilds(),oldCholeskyRefactors=system.choleskyRefactors(),oldCholeskyReuses=system.choleskyReuses();
        const auto iterations = pressure ? system.solvePressure(field, workspace,
            c.pressurePreconditioner == PressurePreconditioner2D::SystemCholesky ? detail::LinearPressureMethod2D::SystemCholesky :
            c.pressurePreconditioner == PressurePreconditioner2D::Aggregation ? detail::LinearPressureMethod2D::Aggregation :
            (c.pressurePreconditioner == PressurePreconditioner2D::IncompleteCholesky0 ? detail::LinearPressureMethod2D::IC0 : detail::LinearPressureMethod2D::Jacobi),linearRelativeTolerance)
            : system.solve(field, workspace, momentumRowStop,std::numeric_limits<double>::infinity(),detail::LinearSolveMethod2D::Jacobi,linearRelativeTolerance);
        if(c.profile && pressure) {
            for(std::size_t i=0;i<n;++i)workspace.ax[i]=system.compensatedResidualRow(i,field);
            pressureLinearResidual=std::max(pressureLinearResidual,
                detail::linearNorm(workspace.ax)/(c.speed*shortestFace));
        }
        if(c.profile && !pressure) {
            for(std::size_t i=0;i<n;++i) {
                const double scaled=std::abs(system.compensatedResidualRow(i,field))/(system.diag[i]*c.speed);
                if(scaled>predictorResidual) {predictorResidual=scaled;predictorWorstCell=i;}
            }
        }
        if (c.profile) {
            auto& p = r.performance;
            if (pressure) {
                p.pressureFactorizations += system.ic0Builds() - oldBuilds;
                p.pressureFactorReuses += system.ic0Reuses() - oldReuses;
                p.pressureCholeskyBuilds+=system.choleskyBuilds()-oldCholeskyBuilds;
                p.pressureCholeskyRefactors+=system.choleskyRefactors()-oldCholeskyRefactors;
                p.pressureCholeskyReuses+=system.choleskyReuses()-oldCholeskyReuses;
                p.pressureHierarchyBuilds += system.hierarchyBuilds() - oldHierarchies;
                p.pressureHierarchyReuses += system.hierarchyReuses() - oldHierarchyReuses;
                p.pressureHierarchyRefreshes += system.hierarchyRefreshes() - oldHierarchyRefreshes;
                p.maxPressureHierarchyLevels=std::max(p.maxPressureHierarchyLevels,system.hierarchyLevels());
                p.maxPressureCoarseCells=std::max(p.maxPressureCoarseCells,system.hierarchyCoarseCells());
            }
            (pressure ? p.pressureSolves : p.momentumSolves) += 1;
            (pressure ? p.pressureIterations : p.momentumIterations) += iterations;
            auto& maximum = pressure ? p.maxPressureIterations : p.maxMomentumIterations;
            maximum = std::max(maximum, iterations);
            (pressure ? p.pressureLinearSolveSeconds : p.momentumLinearSolveSeconds) +=
                std::chrono::duration<double>(Clock::now() - start).count();
        }
        return iterations;
    };
    r.u.resize(n);
    r.v.resize(n);
    r.p.resize(n);
    r.flux.resize(nf);
    r.domainHeight = b.ymax - b.ymin;
    initializeCaseSources(m,c,r);
    Vec zeros(nf);
    const auto& pressureBoundary = c.scenario == "custom" ? b.p : zeros;
    // fixedP is constant. Velocity caches detect every mask change caused by
    // backflow, including candidate rejection and physical-time initialization.
    FlowEquation2D equation(m,b.fixedP);
    const auto& pressureGradientStencil=equation.pressureStencil();
    Vec ra(n), correctionResponse(c.coupling==FlowCoupling2D::SimpleC?n:0);
    Vec pc(n);
    Vec df(nf);
    const double h = b.ymax - b.ymin;
    const double pressureScale = finite(c.speed * c.speed + c.nu * c.speed / h);
    ensure(pressureScale > 0, "Flow reference pressure scale underflow");
    initializeCaseVelocity(m,c,b,r);
    if(guess) {
        r.u=guess->u;r.v=guess->v;r.p=guess->p;
        if(b.closed) {const double gauge=r.p.front();for(double& p:r.p)p=finite(p-gauge);}
    }
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = m.faces[id];
        r.flux[id] =
            f.neighbour
                ? interpolate(f, r.u) * f.areaVector.x
                : (b.role[id] == Role::Inlet
                       ? b.u[id] * f.areaVector.x
                       : ((b.role[id] == Role::Outlet || b.role[id] == Role::Farfield) ? r.u[f.owner] * f.areaVector.x : 0.));
        if (c.scenario == "custom")
            r.flux[id] = f.neighbour ? interpolate(f, r.u)*f.areaVector.x + interpolate(f, r.v)*f.areaVector.y
                : b.role[id] == Role::Inlet ? b.u[id]*f.areaVector.x + b.v[id]*f.areaVector.y
                : (b.role[id] == Role::Outlet || b.role[id] == Role::Opening) ? r.u[f.owner]*f.areaVector.x + r.v[f.owner]*f.areaVector.y : 0;
    }

    if(guess)for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id];
        r.flux[id]=f.neighbour ? interpolate(f,r.u)*f.areaVector.x+interpolate(f,r.v)*f.areaVector.y
            : b.role[id]==Role::Inlet ? b.u[id]*f.areaVector.x+b.v[id]*f.areaVector.y
            : (b.role[id]==Role::Outlet || b.role[id]==Role::Opening || b.role[id]==Role::Farfield)
                ? r.u[f.owner]*f.areaVector.x+r.v[f.owner]*f.areaVector.y : 0.;
        finite(r.flux[id]);
    }

    if(guess && !guess->flux.empty()) {
        r.flux=guess->flux;
        for(std::size_t id=0;id<nf;++id) {
            const auto& face=m.faces[id];
            if(face.neighbour || b.role[id]==Role::Outlet ||
               b.role[id]==Role::Opening || b.role[id]==Role::Farfield) continue;
            const double expected=b.role[id]==Role::Inlet
                ? b.u[id]*face.areaVector.x+b.v[id]*face.areaVector.y : 0.;
            const double scale=c.speed*std::hypot(face.areaVector.x,face.areaVector.y)+std::abs(expected);
            ensure(std::abs(r.flux[id]-expected)<=64*std::numeric_limits<double>::epsilon()*scale,
                   "Initial face flux violates prescribed boundary flux");
            r.flux[id]=expected;
        }
    }

    Vec oldFluxDefect(nf);
    if (previous) {
        r.u=previous->u; r.v=previous->v; r.p=previous->p; r.flux=previous->flux;
        updateOutletBoundary(b,m,c,r.flux);
        r.time=previous->time+timeStep; r.timeStep=timeStep;
        r.previousU=previous->u; r.previousV=previous->v;
        const auto oldGu=equation.velocityGradient(r.u,b,false),oldGv=equation.velocityGradient(r.v,b,true);
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=m.faces[id];
            if (!f.neighbour) {
                if (b.role[id]==Role::Outlet || b.role[id]==Role::Opening)
                    oldFluxDefect[id]=r.flux[id]-r.u[f.owner]*f.areaVector.x-r.v[f.owner]*f.areaVector.y;
                continue; // fixed-velocity and impermeable boundaries impose their new-time flux
            }
            const auto i=f.owner,j=*f.neighbour;
            const double w=f.neighbourWeight;
            const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,
                                m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
            const auto skew=f.centre-point;
            const double uf=interpolate(f,r.u)+dot(interpolateGradient(f,oldGu),skew);
            const double vf=interpolate(f,r.v)+dot(interpolateGradient(f,oldGv),skew);
            oldFluxDefect[id]=r.flux[id]-uf*f.areaVector.x-vf*f.areaVector.y;
        }
    }

    // The unrelaxed systems used to measure the accepted current iterate's
    // residual are exactly the systems needed at the start of the next SIMPLE
    // iteration. Refresh after every field/flux/boundary update, then transfer
    // numeric storage and apply relaxation without rebuilding the same rows.
    FlowReconstruction2D reconstruction;
    auto& gp=reconstruction.gp;auto& gu=reconstruction.gu;auto& gv=reconstruction.gv;
    auto& forceGradient=reconstruction.pressureForce;auto& stressCorrection=reconstruction.stress;
    auto& pressureFaces=reconstruction.pressureFaces;
    bool materialConverged=!material;
    const auto refreshMomentum = [&](bool updateMaterial=false) {
        updateOutletBoundary(b,m,c,r.flux);
        if (material && updateMaterial) {
            std::vector<detail::MaterialBoundary2D> snapshot(nf);
            for(std::size_t id=0;id<nf;++id)if(!m.faces[id].neighbour)
                snapshot[id]={{b.u[id],b.v[id]},b.fixedU[id],b.fixedV[id],
                    b.role[id]==Role::Wall||b.role[id]==Role::Lid,
                    b.role[id]==Role::Inlet || ((b.role[id]==Role::Farfield || b.role[id]==Role::Opening) && r.flux[id]<0),
                    b.role[id]==Role::Outlet || ((b.role[id]==Role::Farfield || b.role[id]==Role::Opening) && r.flux[id]>=0)};
            auto materialState=material(r,snapshot);
            materialConverged=materialState.converged;
            c.faceViscosity=std::move(materialState.faceViscosity);
            ensure(c.faceViscosity.size()==nf,"Material update must supply every face viscosity");
            validateViscosity(m,c);
        }
        reconstruction=equation.reconstruct(r,c,b,pressureBoundary,r.flux);
        equation.assembleMomentum(checkU,checkV,r,c,b,reconstruction,r.flux,r.sourceIntegrals,previous,timeStep);
    };
    const auto takeRelaxed = [&](System& destination,System& source,const Vec& field) {
        destination.diag.swap(source.diag);
        destination.off.swap(source.off);
        destination.rhs.swap(source.rhs);
        for (std::size_t i=0;i<n;++i) {
            const double old=destination.diag[i];
            destination.diag[i]/=c.velocityRelaxation;
            destination.rhs[i]+=(destination.diag[i]-old)*field[i];
        }
    };
    // Fixed boundary groups keep monitoring deterministic and avoid cancelling
    // distinct named wall loads. Values are dimensionless using fixed reference
    // scales, so zero lift/flow never causes a relative division by zero.
    std::vector<std::size_t> monitorGroup(nf);
    Vec monitorLengths;
    prepareMonitors(m,c,b,r,monitorGroup,monitorLengths);
    const bool accelerated=c.steadyAcceleration==SteadyAcceleration2D::Anderson;
    detail::Anderson2D accelerator(c.andersonHistory);
    Vec stateScale;
    if(accelerated) {
        double area=0,length=0;
        for(const auto& cell:m.cells)area+=cell.area;
        for(const auto& face:m.faces)length+=std::hypot(face.areaVector.x,face.areaVector.y);
        stateScale.resize(3*n+nf);
        for(std::size_t i=0;i<n;++i) {
            const double weight=std::sqrt(m.cells[i].area/area);
            stateScale[i]=stateScale[n+i]=weight/c.speed;stateScale[2*n+i]=weight/pressureScale;
        }
        for(std::size_t i=0;i<nf;++i) {
            const double size=std::hypot(m.faces[i].areaVector.x,m.faces[i].areaVector.y);
            stateScale[3*n+i]=std::sqrt(size/length)/(c.speed*size);
        }
    }
    const auto packState=[&] {
        Vec state;state.reserve(3*n+nf);
        for(const auto* field:{&r.u,&r.v,&r.p,&r.flux})state.insert(state.end(),field->begin(),field->end());
        for(std::size_t i=0;i<state.size();++i)state[i]*=stateScale[i];
        return state;
    };
    const auto unpackState=[&](const Vec& state) {
        std::size_t i=0;
        for(auto* field:{&r.u,&r.v,&r.p,&r.flux})for(auto& value:*field){value=state[i]/stateScale[i];++i;}
    };
    refreshMomentum();
    for (std::size_t it = 1; it <= c.maxIterations; ++it) {
        const bool finalCertification=certifyNext;
        strictLinearIteration=!c.adaptiveLinear || certifyNext;
        certifyNext=false;
        linearRelativeTolerance=linearPolicy.relative(strictLinearIteration);
        momentumRowStop=linearPolicy.row(momentumScaledStop,strictLinearIteration,c.speed*c.velocityRelaxation);
        const Vec oldU = r.u;
        const Vec oldV = r.v;
        const Vec oldP = r.p;
        const Vec previousScaled=accelerated ? packState() : Vec{};
        takeRelaxed(au,checkU,r.u);
        takeRelaxed(av,checkV,r.v);
        // Use one pressure response for both components. Slip constraints can
        // give different diagonals; extra implicit relaxation preserves each
        // original fixed-point equation while making rAU scalar and consistent.
        for (std::size_t i = 0; i < n; ++i) {
            const double common = std::max(au.diag[i], av.diag[i]);
            au.rhs[i] += (common - au.diag[i]) * r.u[i];
            av.rhs[i] += (common - av.diag[i]) * r.v[i];
            au.diag[i] = av.diag[i] = common;
        }
        predictorResidual=0;predictorWorstCell=0;pressureLinearResidual=0;
        linearSolve(au, r.u, false);linearSolve(av, r.v, false);
        // Both components share the scalar pressure response away from slip walls.
        for(std::size_t i=0;i<n;++i)ra[i]=m.cells[i].area/au.diag[i];
        const auto gup=equation.velocityGradient(r.u,b,false),gvp=equation.velocityGradient(r.v,b,true);
        const auto predicted=rhieChowFlux(m,c,b,r,ra,oldU,oldV,gu,gv,gup,gvp,gp,forceGradient,
            oldFluxDefect,previous!=nullptr,timeStep,df);
        if(c.coupling==FlowCoupling2D::SimpleC) {
            // Approximate neighbouring velocity corrections by the local
            // correction: aP' = aP + sum(aPN), with signed off-diagonals.
            // This changes only the correction preconditioner, not the
            // predicted face equation or the residual used for acceptance.
            for(std::size_t i=0;i<n;++i) {
                double du=au.diag[i],dv=av.diag[i];
                for(auto k=pattern.rows[i];k<pattern.rows[i+1];++k){du+=au.off[k];dv+=av.off[k];}
                const double diagonal=std::max(du,dv);
                ensure(std::isfinite(diagonal) && diagonal>0,
                       "SIMPLEC correction diagonal is not positive; reduce velocity relaxation");
                correctionResponse[i]=m.cells[i].area/diagonal;
            }
            for(std::size_t id=0;id<nf;++id)df[id]=interpolate(m.faces[id],correctionResponse)*m.faces[id].transmissibility;
        }
        const auto& pressureResponse=c.coupling==FlowCoupling2D::SimpleC?correctionResponse:ra;
        // rAU and the orthogonal pressure coefficients stay fixed across the
        // four non-orthogonal corrections. Assemble and pin once; only the
        // explicit correction/RHS changes. The next SIMPLE iteration resets
        // the matrix and invalidates its IC(0) factorization.
        const auto correction=solvePressureCorrection(m,c,b,r,ap,pc,pressureResponse,df,predicted,pressureGradientStencil,zeros,linearSolve);
        const auto gc=equation.pressureForce(pc,zeros,b.fixedP);
        double du=0,dp=0;
        for(std::size_t i=0;i<n;++i){r.u[i]-=pressureResponse[i]*gc[i].x;r.v[i]-=pressureResponse[i]*gc[i].y;r.p[i]+=c.pressureRelaxation*pc[i];
            finite(r.u[i]);finite(r.v[i]);finite(r.p[i]);du=std::max(du,std::hypot(r.u[i]-oldU[i],r.v[i]-oldV[i])/c.speed);dp=std::max(dp,std::abs(r.p[i]-oldP[i])/pressureScale);}
        Vec div(n);r.globalImbalance=0;
        for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];double flux=predicted[id];
            if(f.neighbour||b.fixedP[id])flux+=df[id]*(pc[f.owner]-(f.neighbour?pc[*f.neighbour]:0))+correction[id];
            r.flux[id]=finite(flux);div[f.owner]+=flux;if(f.neighbour)div[*f.neighbour]-=flux;else r.globalImbalance+=flux;
        }
        double continuity=0;for(std::size_t i=0;i<n;++i)continuity=std::max(continuity,std::abs(div[i])/(c.speed*std::sqrt(m.cells[i].area)));
        double inflow=0;for(std::size_t id=0;id<nf;++id)if(!m.faces[id].neighbour)inflow+=std::max(0.,-r.flux[id]);
        const double flowScale=b.closed?finite(c.speed*h):finite(inflow);
        // An exactly quiescent pressure-driven domain has no throughput to
        // normalize. Only an exactly zero imbalance is admissible in that case.
        ensure(flowScale>0 || (b.pressureOpenings && r.globalImbalance==0),"Flow has no positive reference throughput");
        r.globalRelativeImbalance=flowScale>0 ? finite(std::abs(r.globalImbalance)/flowScale) : 0;
        refreshMomentum(true);
        checkU.apply(r.u,mu);checkV.apply(r.v,mv);double mr=0;
        std::size_t worstCell=0;double worstX=0,worstY=0;
        for(std::size_t i=0;i<n;++i){const double scale=finite((checkU.diag[i]+checkV.diag[i])*c.speed);
            ensure(scale>0,"Flow momentum scale underflow");
            const double residual=std::hypot(mu[i]-checkU.rhs[i],mv[i]-checkV.rhs[i])/scale;
            if(residual>mr) {mr=residual;worstCell=i;worstX=(mu[i]-checkU.rhs[i])/scale;worstY=(mv[i]-checkV.rhs[i])/scale;}}
        bool acceleratedCandidate=false;
        const bool baseConverged=detail::strictFlowConverged2D(it,mr,du,dp,continuity,r.globalRelativeImbalance,c.tolerance,materialConverged);
        if(accelerated && !finalCertification && it>=c.andersonStart && !baseConverged && !(c.adaptiveLinear && strictLinearIteration)) {
            const auto candidate=accelerator.propose(previousScaled,packState());
            if(candidate) {
                ++r.performance.accelerationCandidates;
                const Vec baseU=r.u,baseV=r.v,baseP=r.p,baseFlux=r.flux;
                try {
                    unpackState(*candidate);
                    Vec candidateDivergence(n);double imbalance=0,incoming=0,candidateContinuity=0;
                    bool allowed=true;
                    for(std::size_t id=0;id<nf;++id) {
                        const auto& face=m.faces[id];const double flux=r.flux[id];
                        allowed=allowed && std::isfinite(flux);
                        candidateDivergence[face.owner]+=flux;
                        if(face.neighbour)candidateDivergence[*face.neighbour]-=flux;
                        else {
                            imbalance+=flux;incoming+=std::max(0.,-flux);
                            if(b.role[id]==Role::Outlet && c.outletBackflow==OutletBackflow2D::Reject && flux<0)allowed=false;
                        }
                    }
                    for(std::size_t i=0;i<n;++i)candidateContinuity=std::max(candidateContinuity,
                        std::abs(candidateDivergence[i])/(c.speed*std::sqrt(m.cells[i].area)));
                    const double reference=b.closed?c.speed*h:incoming;
                    const double relative=reference>0 ? std::abs(imbalance)/reference : (imbalance==0?0:1);
                    allowed=allowed && std::isfinite(candidateContinuity) && candidateContinuity<1e-8 && relative<1e-8;
                    if(allowed) {
                        refreshMomentum();checkU.apply(r.u,mu);checkV.apply(r.v,mv);
                        double candidateResidual=0,candidateX=0,candidateY=0;std::size_t candidateWorst=0;
                        for(std::size_t i=0;i<n;++i) {
                            const double scale=(checkU.diag[i]+checkV.diag[i])*c.speed;
                            const double x=(mu[i]-checkU.rhs[i])/scale,y=(mv[i]-checkV.rhs[i])/scale;
                            const double residual=std::hypot(x,y);
                            allowed=allowed && std::isfinite(residual);
                            if(residual>candidateResidual){candidateResidual=residual;candidateWorst=i;candidateX=x;candidateY=y;}
                        }
                        if(allowed && candidateResidual<mr) {
                            acceleratedCandidate=true;
                            mr=candidateResidual;worstCell=candidateWorst;worstX=candidateX;worstY=candidateY;
                            continuity=candidateContinuity;r.globalImbalance=imbalance;r.globalRelativeImbalance=relative;
                            du=dp=0;
                            for(std::size_t i=0;i<n;++i) {
                                du=std::max(du,std::hypot(r.u[i]-oldU[i],r.v[i]-oldV[i])/c.speed);
                                dp=std::max(dp,std::abs(r.p[i]-oldP[i])/pressureScale);
                            }
                        }
                    }
                } catch(const std::runtime_error&) {
                    // A failed extrapolation is not a failed SIMPLE step.
                    // Restore its exact fields, then rebuild and recheck the
                    // original equations; failures of that restore propagate.
                }
                if(acceleratedCandidate)++r.performance.accelerationAccepted;
                else {
                    ++r.performance.accelerationRejected;
                    r.u=baseU;r.v=baseV;r.p=baseP;r.flux=baseFlux;refreshMomentum();accelerator.clear();
                }
            }
        }
        FlowIteration2D step{it,finite(mr),finite(continuity),finite(du),finite(dp)};
        step.momentumWorstCell=worstCell;step.momentumResidualX=worstX;step.momentumResidualY=worstY;
        step.momentumPredictorResidual=finite(predictorResidual);step.momentumPredictorWorstCell=predictorWorstCell;
        step.pressureLinearResidual=finite(pressureLinearResidual);
        step.linearRelativeTolerance=linearRelativeTolerance;
        step.strictLinearStep=strictLinearIteration && !acceleratedCandidate;
        step.globalRelativeImbalance=r.globalRelativeImbalance;
        step.monitors=physicalMonitors(m,c,b,r,monitorGroup,monitorLengths,pressureScale,pressureFaces,gu,gv,stressCorrection);
        if(it<=10)r.convergenceReference=std::max(r.convergenceReference,step.momentumResidual);
        r.history.push_back(step);
        if(progress&&(it==1||it%10==0))progress(step);
        // Only an ordinary step can certify the selected convergence mode.
        // An extrapolated field is never the final
        // convergence proof, even when its momentum residual is small.
        const bool stoppingCandidate=c.convergence==FlowConvergence2D::Engineering
            ? detail::engineeringWindow2D(r.history,r.convergenceReference,c.tolerance).accepted : baseConverged;
        if(stoppingCandidate) {
            if(strictLinearIteration && !acceleratedCandidate){r.converged=true;break;}
            certifyNext=true;
        }
        if(c.adaptiveLinear)linearPolicy.observe(std::max({mr,du,dp,continuity}));
        if(c.stopRequested && c.stopRequested()){r.stopped=true;break;}
    }
    if (previous) {
        r.temporalIntegrals.resize(n);
        Vec absoluteFlux(n);
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=m.faces[id]; const double q=std::abs(r.flux[id]);
            absoluteFlux[f.owner]+=q;
            if (f.neighbour) absoluteFlux[*f.neighbour]+=q;
        }
        for (std::size_t i=0;i<n;++i) {
            r.temporalIntegrals[i]={finite(m.cells[i].area*(r.u[i]-previous->u[i])/timeStep),
                                    finite(m.cells[i].area*(r.v[i]-previous->v[i])/timeStep)};
            r.maxCourant=std::max(r.maxCourant,finite(.5*timeStep*absoluteFlux[i]/m.cells[i].area));
        }
    }
    postprocessForces(m,c,b,r,pressureBoundary,gp,gu,gv,stressCorrection);
    if (c.profile) {
        r.performance.solveSeconds = std::chrono::duration<double>(Clock::now() - solveStart).count();
    }
    return r;
}

FlowResult2D solveIncompressible2D(const FvMesh2D& m, const FlowControls2D& c,
    const std::function<void(const FlowIteration2D&)>& progress) {
    return solveFlow(m,c,progress,nullptr,0);
}
FlowResult2D solveIncompressibleFromGuess2D(const FvMesh2D& m,const FlowControls2D& c,
    const FlowInitialGuess2D& guess,const std::function<void(const FlowIteration2D&)>& progress) {
    return solveFlow(m,c,progress,nullptr,0,{},&guess);
}
FlowResult2D detail::solveMaterialFlow2D(const FvMesh2D& m,const FlowControls2D& c,
    const MaterialUpdate2D& material,const std::function<void(const FlowIteration2D&)>& progress) {
    ensure(bool(material),"Material flow requires a constitutive update");
    ensure(c.scenario=="channel"||c.scenario=="cavity"||c.scenario=="external"||c.scenario=="flatplate",
        "Material flow supports steady physical cases only");
    return solveFlow(m,c,progress,nullptr,0,material);
}
FlowResult2D advanceIncompressible2D(const FvMesh2D& m, const FlowControls2D& c,
    const FlowState2D& previous, double timeStep,
    const std::function<void(const FlowIteration2D&)>& progress) {
    return solveFlow(m,c,progress,&previous,timeStep);
}

}
