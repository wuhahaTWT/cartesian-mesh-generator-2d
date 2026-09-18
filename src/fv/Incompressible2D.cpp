#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/TaylorGreen2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {

using Vec = std::vector<double>;

void ensure(bool valid, const char* message) {
    if (!valid) throw std::runtime_error(message);
}

double finite(double x) {
    ensure(std::isfinite(x), "Flow numerical range exceeded");
    return x;
}

using System = detail::SparseSystem2D;

enum class Role { Wall, Inlet, Outlet, Slip, Lid };

struct Boundary {
    std::vector<Role> role;
    Vec u;
    Vec v;
    std::vector<bool> fixedU;
    std::vector<bool> fixedV;
    std::vector<bool> fixedP;
    std::vector<bool> constantU;
    std::vector<bool> constantV;
    double xmin = 0;
    double xmax = 0;
    double ymin = 0;
    double ymax = 0;
    bool closed = false;
};
Boundary boundaries(const FvMesh2D& m, const FlowControls2D& c) {
    const auto nf = m.faces.size();
    Boundary b;
    b.role.assign(nf, Role::Wall);
    b.u.resize(nf);
    b.v.resize(nf);
    b.fixedU.assign(nf, false);
    b.fixedV.assign(nf, false);
    b.fixedP.assign(nf, false);
    b.constantU.assign(nf, false);
    b.constantV.assign(nf, false);
    b.xmin = b.ymin = std::numeric_limits<double>::infinity();
    b.xmax = b.ymax = -b.xmin;
    for (const auto& f : m.faces) {
        if (!f.neighbour) {
            b.xmin = std::min(b.xmin, f.centre.x);
            b.xmax = std::max(b.xmax, f.centre.x);
            b.ymin = std::min(b.ymin, f.centre.y);
            b.ymax = std::max(b.ymax, f.centre.y);
        }
    }
    const double width = b.xmax - b.xmin;
    const double height = b.ymax - b.ymin;
    ensure(width > 0 && height > 0, "Flow empty domain");
    const double eps = TolerancePolicy{}.scale(std::max(width, height));
    const auto equal = [&](double a, double d) { return std::abs(a - d) <= eps; };

    std::size_t inlets = 0;
    std::size_t outlets = 0;
    std::size_t walls = 0;
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = m.faces[id];
        if (f.neighbour) {
            continue;
        }
        const bool left = equal(f.centre.x, b.xmin);
        const bool right = equal(f.centre.x, b.xmax);
        const bool top = equal(f.centre.y, b.ymax);
        const bool bottom = equal(f.centre.y, b.ymin);
        const bool embedded =
            c.scenario == "external" && f.patch == BoundaryPatch2D::EmbeddedBoundary;
        if (embedded) {
            b.role[id] = Role::Wall;
            ++walls;
        } else {
            ensure(left || right || top || bottom,
                   "Selected flow case requires rectangular outer boundary");
            // An axis-position match alone is insufficient for a sloping edge.
            ensure((left || right) ? std::abs(f.areaVector.y) <= eps
                                   : std::abs(f.areaVector.x) <= eps,
                   "Flow outer edge is not axis aligned");
            if (c.scenario == "taylor-green") {
                b.role[id] = Role::Slip;
            } else if (c.scenario == "manufactured") {
                b.role[id] = Role::Wall;
            } else if (c.scenario == "cavity") {
                b.role[id] = top ? Role::Lid : Role::Wall;
            } else if (left) {
                b.role[id] = Role::Inlet;
                ++inlets;
            } else if (right) {
                b.role[id] = Role::Outlet;
                ++outlets;
            } else {
                b.role[id] = c.scenario == "external" ? Role::Slip : Role::Wall;
            }
        }

        b.constantU[id]=b.role[id]==Role::Wall || b.role[id]==Role::Lid;
        b.constantV[id]=b.constantU[id];
        switch (b.role[id]) {
        case Role::Wall:
            b.fixedU[id] = b.fixedV[id] = true;
            break;
        case Role::Lid:
            b.fixedU[id] = b.fixedV[id] = true;
            b.u[id] = c.speed;
            break;
        case Role::Inlet:
            b.fixedU[id] = b.fixedV[id] = true;
            b.u[id] = c.scenario == "channel"
                          ? 4 * c.speed * (f.centre.y - b.ymin) * (b.ymax - f.centre.y) /
                                (height * height)
                          : c.speed;
            break;
        case Role::Outlet:
            b.fixedP[id] = true;
            break;
        case Role::Slip:
            // Outer boundaries were checked to be axis aligned above.
            if (left || right) b.fixedU[id] = b.constantU[id] = true;
            else b.fixedV[id] = b.constantV[id] = true;
            break;
        }
    }
    b.closed = c.scenario == "cavity" || c.scenario == "manufactured" || c.scenario == "taylor-green";
    if (c.scenario == "manufactured" || c.scenario == "taylor-green")
        ensure(equal(b.xmin,0) && equal(b.ymin,0) && equal(b.xmax,1) && equal(b.ymax,1),
               "Verification flow requires the unit square [0,1]^2");
    if (!b.closed) {
        ensure(inlets && outlets, "Flow needs left inlet and right pressure outlet");
    }
    if (c.scenario == "external") {
        ensure(walls > 0, "External case requires embedded solid wall");
    } else {
        double area = 0;
        for (const auto& cell : m.cells) {
            area += cell.area;
        }
        ensure(std::abs(area - width * height) <
                   TolerancePolicy{}.areaScale(std::max(width, height)) * 10,
               "Channel/cavity requires a filled rectangle without holes");
    }
    return b;
}

using detail::flowGradient;

Vector2D interpolateGradient(const Face& f, const std::vector<Vector2D>& g) {
    auto result = g[f.owner];
    if (f.neighbour) {
        const auto neighbour = g[*f.neighbour];
        const double weight = f.neighbourWeight;
        result = {result.x * (1 - weight) + neighbour.x * weight,
                  result.y * (1 - weight) + neighbour.y * weight};
    }
    return result;
}

double interpolate(const Face& f, const Vec& x) {
    return f.neighbour ? x[f.owner] * (1 - f.neighbourWeight) +
                             x[*f.neighbour] * f.neighbourWeight
                       : x[f.owner];
}

void momentum(System& a,
                const FvMesh2D& m,
                const FlowControls2D& c,
                const Boundary& b,
                const Vec& field,
                const Vec& flux,
                const std::vector<Vector2D>& gradField,
                const std::vector<Vector2D>& gp,
                const std::vector<Vector2D>& source,
                const std::vector<Vector2D>& stressCorrection,
                bool y,
                bool relaxed,
                const Vec* previous = nullptr,
                double timeStep = 0) {
    a.reset();
    const auto& bc = y ? b.v : b.u;
    const auto& fixed = y ? b.fixedV : b.fixedU;
    const auto limiter = c.convection == ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m, field, gradField, bc, fixed) : Vec{};
    for (std::size_t i = 0; i < m.cells.size(); ++i) {
        a.rhs[i] = -m.cells[i].area * (y ? gp[i].y : gp[i].x);
        if (!source.empty()) a.rhs[i] += y ? source[i].y : source[i].x;
        if (previous) {
            const double mass=finite(m.cells[i].area/timeStep);
            a.diag[i]+=mass;
            a.rhs[i]+=finite(mass*(*previous)[i]);
        }
    }
    for (std::size_t id = 0; id < m.faces.size(); ++id) {
        const auto& f = m.faces[id];
        const auto i = f.owner;
        const double q = flux[id];
        const double d = c.nu * f.transmissibility;
        if (!stressCorrection.empty()) {
            const double extra=y?stressCorrection[id].y:stressCorrection[id].x;
            a.rhs[i]-=extra;
            if (f.neighbour) a.rhs[*f.neighbour]+=extra;
        }
        if (f.neighbour) {
            const auto j = *f.neighbour;
            a.diag[i] += d + std::max(q, 0.);
            a.add(i, j, -d + std::min(q, 0.));
            a.diag[j] += d + std::max(-q, 0.);
            a.add(j, i, -d - std::max(q, 0.));
            const double correction =
                c.nu * dot(interpolateGradient(f, gradField), f.correction);
            const double upwind = q >= 0 ? field[i] : field[j];
            const double deferred = q * (detail::upwindFaceValue(m, id, q, field, gradField, limiter)-upwind);
            a.rhs[i] += correction-deferred;
            a.rhs[j] -= correction-deferred;
        } else {
            if (fixed[id]) {
                a.diag[i] += d;
                a.rhs[i] += d * bc[id] + c.nu * dot(gradField[i], f.correction);
                a.rhs[i] -= q * bc[id];
            } else {
                ensure(q >= -1e-12 * c.speed * std::hypot(f.areaVector.x, f.areaVector.y),
                       "Flow outlet backflow unsupported in this laminar prototype");
                a.diag[i] += q;
                a.rhs[i] -= q * (detail::upwindFaceValue(m, id, q, field, gradField, limiter)-field[i]);
            }
        }
    }
    if (relaxed) {
        for (std::size_t i = 0; i < m.cells.size(); ++i) {
            const double old = a.diag[i];
            a.diag[i] /= c.velocityRelaxation;
            a.rhs[i] += (a.diag[i] - old) * field[i];
        }
    }
}

} // namespace

static FlowResult2D solveFlow(
    const FvMesh2D& m, const FlowControls2D& c,
    const std::function<void(const FlowIteration2D&)>& progress,
    const FlowState2D* previous, double timeStep) {
    using Clock = std::chrono::steady_clock;
    const auto solveStart = c.profile ? Clock::now() : Clock::time_point{};
    validateFvMesh2D(m);
    ensure((c.scenario == "external" || c.scenario == "channel" || c.scenario == "cavity" || c.scenario == "manufactured" || c.scenario == "taylor-green") &&
               std::isfinite(c.nu) && c.nu > 0 && std::isfinite(c.speed) && c.speed > 0 &&
               std::isfinite(c.tolerance) && c.tolerance > 0 && c.maxIterations > 0,
           "Invalid flow controls");
    ensure(std::isfinite(c.manufacturedPressureSlope) &&
               (c.scenario == "manufactured" || c.manufacturedPressureSlope == 0),
           "Manufactured pressure slope is only valid for the verification case");
    ensure(c.velocityRelaxation > 0 && c.velocityRelaxation <= 1 &&
               c.pressureRelaxation > 0 && c.pressureRelaxation <= 1,
           "Invalid SIMPLE relaxation");
    ensure(c.pressurePreconditioner == PressurePreconditioner2D::Jacobi ||
               c.pressurePreconditioner == PressurePreconditioner2D::IncompleteCholesky0,
           "Invalid pressure preconditioner");
    ensure(c.convection == ConvectionScheme2D::Upwind ||
               c.convection == ConvectionScheme2D::LimitedLinearUpwind,
           "Invalid convection scheme");
    ensure(c.viscousStress == ViscousStress2D::Symmetric ||
               c.viscousStress == ViscousStress2D::Laplacian,
           "Invalid viscous stress form");
    const auto b = boundaries(m, c);
    const auto n = m.cells.size();
    const auto nf = m.faces.size();
    if (previous) {
        ensure(std::isfinite(timeStep) && timeStep>0 && std::isfinite(previous->time) && previous->time>=0 &&
                   std::isfinite(previous->time+timeStep) && previous->time+timeStep>previous->time,
               "Invalid physical time step");
        ensure(c.scenario!="manufactured", "Transient forced manufactured case not implemented");
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
    const double momentumScaledStop=previous?finite(.01*c.tolerance*c.speed*c.velocityRelaxation)
        :std::numeric_limits<double>::infinity();
    ensure(momentumScaledStop>0,"Transient linear residual scale underflow");
    // Profiling observes the same solves and stopping rules, including zero-step
    // solves. Timing includes each linear solver's setup, but not assembly.
    auto linearSolve = [&](const System& system, Vec& field, bool pressure) {
        const auto start = c.profile ? Clock::now() : Clock::time_point{};
        const auto iterations = pressure ? system.solvePressure(field, workspace,
            c.pressurePreconditioner == PressurePreconditioner2D::IncompleteCholesky0)
            : system.solve(field, workspace, momentumScaledStop);
        if (c.profile) {
            auto& p = r.performance;
            (pressure ? p.pressureSolves : p.momentumSolves) += 1;
            (pressure ? p.pressureIterations : p.momentumIterations) += iterations;
            auto& maximum = pressure ? p.maxPressureIterations : p.maxMomentumIterations;
            maximum = std::max(maximum, iterations);
            (pressure ? p.pressureLinearSolveSeconds : p.momentumLinearSolveSeconds) +=
                std::chrono::duration<double>(Clock::now() - start).count();
        }
    };
    r.u.resize(n);
    r.v.resize(n);
    r.p.resize(n);
    r.flux.resize(nf);
    r.domainHeight = b.ymax - b.ymin;
    if (c.scenario == "manufactured") {
        r.sourceIntegrals.reserve(n);
        for (const auto& cell : m.cells) {
            const auto acceleration=manufacturedFlow2D(cell.centre,c.speed,c.nu,c.manufacturedPressureSlope).acceleration;
            r.sourceIntegrals.push_back({finite(cell.area*acceleration.x),finite(cell.area*acceleration.y)});
        }
    }
    Vec zeros(nf);
    Vec ra(n);
    Vec pc(n);
    Vec df(nf);
    const double h = b.ymax - b.ymin;
    const double pressureScale = finite(c.speed * c.speed + c.nu * c.speed / h);
    ensure(pressureScale > 0, "Flow reference pressure scale underflow");
    for (std::size_t i = 0; i < n; ++i) {
        const double y = m.cells[i].centre.y;
        r.u[i] = b.closed ? 0
                          : (c.scenario == "channel"
                                 ? 4 * c.speed * (y - b.ymin) * (b.ymax - y) / (h * h)
                                 : c.speed);
    }
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = m.faces[id];
        r.flux[id] =
            f.neighbour
                ? interpolate(f, r.u) * f.areaVector.x
                : (b.role[id] == Role::Inlet
                       ? b.u[id] * f.areaVector.x
                       : (b.role[id] == Role::Outlet ? r.u[f.owner] * f.areaVector.x : 0.));
    }

    Vec oldFluxDefect(nf);
    if (previous) {
        r.u=previous->u; r.v=previous->v; r.p=previous->p; r.flux=previous->flux;
        r.time=previous->time+timeStep; r.timeStep=timeStep;
        r.previousU=previous->u; r.previousV=previous->v;
        const auto oldGu=flowGradient(m,r.u,b.u,b.fixedU),oldGv=flowGradient(m,r.v,b.v,b.fixedV);
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=m.faces[id];
            if (!f.neighbour) {
                if (b.role[id]==Role::Outlet)
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

    for (std::size_t it = 1; it <= c.maxIterations; ++it) {
        const Vec oldU = r.u;
        const Vec oldV = r.v;
        const Vec oldP = r.p;
        const auto gp = flowGradient(m, r.p, zeros, b.fixedP, true);
        const auto forceGradient = detail::conservativePressureGradient(m,
            detail::pressureFaceValues(m, r.p, gp, zeros, b.fixedP));
        const auto gu = flowGradient(m, r.u, b.u, b.fixedU);
        const auto gv = flowGradient(m, r.v, b.v, b.fixedV);
        const auto stressCorrection=c.viscousStress==ViscousStress2D::Symmetric
            ? detail::symmetricViscousCorrection(m,r.u,r.v,gu,gv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu)
            : std::vector<Vector2D>{};
        momentum(au, m, c, b, r.u, r.flux, gu, forceGradient, r.sourceIntegrals, stressCorrection, false, true, previous?&previous->u:nullptr, timeStep);
        momentum(av, m, c, b, r.v, r.flux, gv, forceGradient, r.sourceIntegrals, stressCorrection, true, true, previous?&previous->v:nullptr, timeStep);
        // Use one pressure response for both components. Slip constraints can
        // give different diagonals; extra implicit relaxation preserves each
        // original fixed-point equation while making rAU scalar and consistent.
        for (std::size_t i = 0; i < n; ++i) {
            const double common = std::max(au.diag[i], av.diag[i]);
            au.rhs[i] += (common - au.diag[i]) * r.u[i];
            av.rhs[i] += (common - av.diag[i]) * r.v[i];
            au.diag[i] = av.diag[i] = common;
        }
        linearSolve(au, r.u, false);linearSolve(av, r.v, false);
        // Both components share the scalar pressure response away from slip walls.
        for(std::size_t i=0;i<n;++i)ra[i]=m.cells[i].area/au.diag[i];
        const auto gup=flowGradient(m,r.u,b.u,b.fixedU),gvp=flowGradient(m,r.v,b.v,b.fixedV);
        Vec predicted(nf);
        for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
            const double rf=interpolate(f,ra);df[id]=rf*f.transmissibility;
            if(f.neighbour){const auto j=*f.neighbour;const double w=f.neighbourWeight;
                const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
                const auto skew=f.centre-point;
                const double uf=interpolate(f,r.u)+dot(interpolateGradient(f,gup),skew),vf=interpolate(f,r.v)+dot(interpolateGradient(f,gvp),skew);
                const Vector2D rag{(1-w)*ra[i]*forceGradient[i].x+w*ra[j]*forceGradient[j].x,(1-w)*ra[i]*forceGradient[i].y+w*ra[j]*forceGradient[j].y};
                predicted[id]=uf*f.areaVector.x+vf*f.areaVector.y+dot(rag,f.areaVector)-rf*(f.transmissibility*(r.p[j]-r.p[i])+dot(interpolateGradient(f,gp),f.correction));
                if (previous) {
                    const double oldUf=interpolate(f,oldU)+dot(interpolateGradient(f,gu),skew);
                    const double oldVf=interpolate(f,oldV)+dot(interpolateGradient(f,gv),skew);
                    // Backward Euler transports the accepted old face flux.
                    // Correct both its cell-interpolation defect and the inner
                    // under-relaxation defect; otherwise the dt -> 0 response
                    // depends on the arbitrary inner relaxation factor.
                    predicted[id]+=rf/timeStep*oldFluxDefect[id]
                        +(1-c.velocityRelaxation)*(r.flux[id]-oldUf*f.areaVector.x-oldVf*f.areaVector.y);
                }
            }else if(b.role[id]==Role::Outlet){predicted[id]=r.u[i]*f.areaVector.x+r.v[i]*f.areaVector.y+ra[i]*dot(forceGradient[i],f.areaVector)-rf*(-f.transmissibility*r.p[i]+dot(gp[i],f.correction));
                if (previous) predicted[id]+=rf/timeStep*oldFluxDefect[id]
                    +(1-c.velocityRelaxation)*(r.flux[id]-oldU[i]*f.areaVector.x-oldV[i]*f.areaVector.y);}
            else if(b.role[id]==Role::Inlet)predicted[id]=b.u[id]*f.areaVector.x+b.v[id]*f.areaVector.y;
            else predicted[id]=0;
        }
        std::fill(pc.begin(),pc.end(),0);Vec correction(nf);
        for(int pass=0;pass<4;++pass){
            ap.reset();const auto gc=flowGradient(m,pc,zeros,b.fixedP,true);
            for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
                correction[id]=(f.neighbour||b.fixedP[id])?-interpolate(f,ra)*dot(interpolateGradient(f,gc),f.correction):0;
                ap.rhs[i]-=predicted[id]+correction[id];
                if(f.neighbour){const auto j=*f.neighbour;ap.rhs[j]+=predicted[id]+correction[id];ap.diag[i]+=df[id];ap.diag[j]+=df[id];ap.add(i,j,-df[id]);ap.add(j,i,-df[id]);}
                else if(b.fixedP[id])ap.diag[i]+=df[id];
            }
            if(b.closed){ap.pin(0);pc[0]=0;}
            linearSolve(ap, pc, true);
        }
        const auto correctionGradient=flowGradient(m,pc,zeros,b.fixedP,true);
        const auto gc=detail::conservativePressureGradient(m,
            detail::pressureFaceValues(m,pc,correctionGradient,zeros,b.fixedP));
        double du=0,dp=0;
        for(std::size_t i=0;i<n;++i){r.u[i]-=ra[i]*gc[i].x;r.v[i]-=ra[i]*gc[i].y;r.p[i]+=c.pressureRelaxation*pc[i];
            finite(r.u[i]);finite(r.v[i]);finite(r.p[i]);du=std::max(du,std::hypot(r.u[i]-oldU[i],r.v[i]-oldV[i])/c.speed);dp=std::max(dp,std::abs(r.p[i]-oldP[i])/pressureScale);}
        Vec div(n);r.globalImbalance=0;
        for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];double flux=predicted[id];
            if(f.neighbour||b.fixedP[id])flux+=df[id]*(pc[f.owner]-(f.neighbour?pc[*f.neighbour]:0))+correction[id];
            r.flux[id]=finite(flux);div[f.owner]+=flux;if(f.neighbour)div[*f.neighbour]-=flux;else r.globalImbalance+=flux;
        }
        double continuity=0;for(std::size_t i=0;i<n;++i)continuity=std::max(continuity,std::abs(div[i])/(c.speed*std::sqrt(m.cells[i].area)));
        double inflow=0;for(std::size_t id=0;id<nf;++id)if(!m.faces[id].neighbour)inflow+=std::max(0.,-r.flux[id]);
        const double flowScale=b.closed?finite(c.speed*h):finite(inflow);
        ensure(flowScale>0,"Flow has no positive reference throughput");
        r.globalRelativeImbalance=finite(std::abs(r.globalImbalance)/flowScale);
        const auto newGp=flowGradient(m,r.p,zeros,b.fixedP,true),newGu=flowGradient(m,r.u,b.u,b.fixedU),newGv=flowGradient(m,r.v,b.v,b.fixedV);
        const auto newForceGradient=detail::conservativePressureGradient(m,
            detail::pressureFaceValues(m,r.p,newGp,zeros,b.fixedP));
        const auto newStressCorrection=c.viscousStress==ViscousStress2D::Symmetric
            ? detail::symmetricViscousCorrection(m,r.u,r.v,newGu,newGv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu)
            : std::vector<Vector2D>{};
        momentum(checkU,m,c,b,r.u,r.flux,newGu,newForceGradient,r.sourceIntegrals,newStressCorrection,false,false,previous?&previous->u:nullptr,timeStep);
        momentum(checkV,m,c,b,r.v,r.flux,newGv,newForceGradient,r.sourceIntegrals,newStressCorrection,true,false,previous?&previous->v:nullptr,timeStep);
        checkU.apply(r.u,mu);checkV.apply(r.v,mv);double mr=0;
        for(std::size_t i=0;i<n;++i){const double scale=finite((checkU.diag[i]+checkV.diag[i])*c.speed);
            ensure(scale>0,"Flow momentum scale underflow");
            mr=std::max(mr,std::hypot(mu[i]-checkU.rhs[i],mv[i]-checkV.rhs[i])/scale);}
        FlowIteration2D step{it,finite(mr),finite(continuity),finite(du),finite(dp)};r.history.push_back(step);
        if(progress&&(it==1||it%10==0))progress(step);
        if(it>=10&&mr<c.tolerance&&du<c.tolerance&&dp<c.tolerance&&continuity<1e-8&&r.globalRelativeImbalance<1e-8){r.converged=true;break;}
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
    const auto gu=flowGradient(m,r.u,b.u,b.fixedU),gv=flowGradient(m,r.v,b.v,b.fixedV),gp=flowGradient(m,r.p,zeros,b.fixedP,true);
    const auto pf=detail::pressureFaceValues(m,r.p,gp,zeros,b.fixedP);
    const auto lu=c.convection==ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m,r.u,gu,b.u,b.fixedU) : Vec{};
    const auto lv=c.convection==ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m,r.v,gv,b.v,b.fixedV) : Vec{};
    const auto stressCorrection=c.viscousStress==ViscousStress2D::Symmetric
        ? detail::symmetricViscousCorrection(m,r.u,r.v,gu,gv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu)
        : std::vector<Vector2D>{};
    r.faceMomentum.resize(nf);
    for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id]; const auto i=f.owner;
        auto& fm=r.faceMomentum[id]; fm.pressure=finite(pf[id]);
        auto component=[&](const Vec& value,const std::vector<Vector2D>& g,
                           const Vec& bc,const std::vector<bool>& fixed,const Vec& limiter) {
            const double faceValue=(!f.neighbour && fixed[id]) ? bc[id]
                : detail::upwindFaceValue(m,id,r.flux[id],value,g,limiter);
            double diffusion=0;
            if(f.neighbour || fixed[id]) {
                const double other=f.neighbour ? value[*f.neighbour] : bc[id];
                diffusion=-c.nu*(f.transmissibility*(other-value[i])+dot(interpolateGradient(f,g),f.correction));
            }
            return std::pair{finite(r.flux[id]*faceValue),finite(diffusion)};
        };
        auto [ax,dx]=component(r.u,gu,b.u,b.fixedU,lu);
        auto [ay,dy]=component(r.v,gv,b.v,b.fixedV,lv);
        if (!stressCorrection.empty()) {dx+=stressCorrection[id].x;dy+=stressCorrection[id].y;}
        fm.advection={ax,ay}; fm.diffusion={dx,dy};
        fm.wall=!f.neighbour && (b.role[id]==Role::Wall || b.role[id]==Role::Lid);
        if (fm.wall) {
            r.wallForceX+=pf[id]*f.areaVector.x+dx;
            r.wallForceY+=pf[id]*f.areaVector.y+dy;
            r.wallViscousForceX+=dx; r.wallViscousForceY+=dy;
        }
        if(f.neighbour || f.patch!=BoundaryPatch2D::EmbeddedBoundary || b.role[id]!=Role::Wall) continue;
        const double px=pf[id]*f.areaVector.x,py=pf[id]*f.areaVector.y;
        r.pressureForceX+=px; r.pressureForceY+=py;
        r.discreteForceX+=px+dx; r.discreteForceY+=py+dy;
        // Preserve the old cell-gradient diagnostic for explicit comparison;
        // symmetric mode reports the actual shared-face stress above as force.
        r.reconstructedForceX+=px-c.nu*(2*gu[i].x*f.areaVector.x+(gu[i].y+gv[i].x)*f.areaVector.y);
        r.reconstructedForceY+=py-c.nu*((gu[i].y+gv[i].x)*f.areaVector.x+2*gv[i].y*f.areaVector.y);
    }
    r.forceX=c.viscousStress==ViscousStress2D::Symmetric?r.discreteForceX:r.reconstructedForceX;
    r.forceY=c.viscousStress==ViscousStress2D::Symmetric?r.discreteForceY:r.reconstructedForceY;
    finite(r.reconstructedForceX);finite(r.reconstructedForceY);
    finite(r.wallForceX);finite(r.wallForceY);
    finite(r.wallViscousForceX);finite(r.wallViscousForceY);
    finite(r.globalImbalance);finite(r.forceX);finite(r.forceY);
    finite(r.pressureForceX);finite(r.pressureForceY);finite(r.discreteForceX);finite(r.discreteForceY);
    if (c.profile) {
        r.performance.solveSeconds = std::chrono::duration<double>(Clock::now() - solveStart).count();
    }
    return r;
}

FlowResult2D solveIncompressible2D(const FvMesh2D& m, const FlowControls2D& c,
    const std::function<void(const FlowIteration2D&)>& progress) {
    return solveFlow(m,c,progress,nullptr,0);
}
FlowResult2D advanceIncompressible2D(const FvMesh2D& m, const FlowControls2D& c,
    const FlowState2D& previous, double timeStep,
    const std::function<void(const FlowIteration2D&)>& progress) {
    return solveFlow(m,c,progress,&previous,timeStep);
}
FlowState2D initialIncompressibleState2D(const FvMesh2D& m, const FlowControls2D& c) {
    validateFvMesh2D(m);
    ensure((c.scenario=="external" || c.scenario=="channel" || c.scenario=="cavity" || c.scenario=="taylor-green") &&
           std::isfinite(c.nu) && c.nu>0 && std::isfinite(c.speed) && c.speed>0,
           "Invalid transient initial-state controls");
    const auto b=boundaries(m,c); (void)b;
    FlowState2D s;
    s.u.resize(m.cells.size());s.v.resize(m.cells.size());s.p.resize(m.cells.size());s.flux.resize(m.faces.size());
    if (c.scenario=="taylor-green") {
        ensure(std::isfinite(c.nu) && c.nu>0 && std::isfinite(c.speed) && c.speed>0,"Invalid vortex controls");
        const double gauge=taylorGreen2D(m.cells.front().centre,0,c.speed,c.nu).pressure;
        for (std::size_t i=0;i<m.cells.size();++i) {
            const auto q=taylorGreen2D(m.cells[i].centre,0,c.speed,c.nu);
            s.u[i]=q.velocity.x;s.v[i]=q.velocity.y;s.p[i]=q.pressure-gauge;
        }
        // Streamfunction differences give an exactly conservative integral flux.
        for (std::size_t id=0;id<m.faces.size();++id) {
            const auto& f=m.faces[id];
            if (!f.neighbour) continue;
            const Point2D a{f.centre.x+.5*f.areaVector.y,f.centre.y-.5*f.areaVector.x};
            const Point2D z{f.centre.x-.5*f.areaVector.y,f.centre.y+.5*f.areaVector.x};
            s.flux[id]=taylorGreen2D(z,0,c.speed,c.nu).streamfunction-taylorGreen2D(a,0,c.speed,c.nu).streamfunction;
        }
    }
    return s;
}
}
