#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
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

double product(const Vec& a, const Vec& b) {
    long double sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        sum += static_cast<long double>(a[i]) * b[i];
    return finite(static_cast<double>(sum));
}

double norm(const Vec& a) {
    long double result = 0;
    for (double value : a)
        result = std::hypot(result, static_cast<long double>(finite(value)));
    return finite(static_cast<double>(result));
}

struct System {
    Vec diag;
    Vec rhs;
    std::vector<std::map<std::size_t, double>> off;

    explicit System(std::size_t n) : diag(n), rhs(n), off(n) {}

    Vec apply(const Vec& x) const {
        Vec y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = diag[i] * x[i];
            for (auto [j, a] : off[i]) y[i] += a * x[j];
            finite(y[i]);
        }
        return y;
    }

    void pin(std::size_t id) {
        for (auto& row : off) row.erase(id);
        off[id].clear();
        rhs[id] = 0;
    }

    std::size_t solvePressure(Vec& x) const {
        // Pressure correction is a symmetric positive graph Laplacian after
        // outlet Dirichlet terms or a symmetric gauge elimination. Preserve
        // conjugacy instead of restarting the nonsymmetric momentum method.
        for (double d : diag) {
            ensure(d > 0 && std::isfinite(d), "Flow pressure matrix diagonal invalid");
        }
        Vec residual(x.size()), z(x.size()), direction(x.size());
        auto ax = apply(x);
        for (std::size_t i = 0; i < x.size(); ++i) {
            residual[i] = rhs[i] - ax[i];
            direction[i] = z[i] = residual[i] / diag[i];
        }
        const double stop = finite(1e-13 + 1e-11 * norm(rhs));
        if (norm(residual) <= stop) {
            return 0;
        }
        double rz = product(residual, z);
        for (std::size_t iteration = 0; iteration < 3000; ++iteration) {
            const auto ad = apply(direction);
            const double denominator = product(direction, ad);
            ensure(denominator > 0 && rz > 0, "Flow pressure PCG lost positive definiteness");
            const double alpha = rz / denominator;
            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] += alpha * direction[i];
                residual[i] -= alpha * ad[i];
            }
            // Only an explicitly recomputed residual can terminate the solve.
            if (norm(residual) <= stop) {
                ax = apply(x);
                for (std::size_t i = 0; i < x.size(); ++i) residual[i] = rhs[i] - ax[i];
                if (norm(residual) <= stop) {
                    return iteration + 1;
                }
                for (std::size_t i = 0; i < x.size(); ++i) {
                    direction[i] = z[i] = residual[i] / diag[i];
                }
                rz = product(residual, z);
                continue;
            }
            for (std::size_t i = 0; i < x.size(); ++i) {
                z[i] = residual[i] / diag[i];
            }
            const double next = product(residual, z);
            const double beta = next / rz;
            for (std::size_t i = 0; i < x.size(); ++i) {
                direction[i] = z[i] + beta * direction[i];
            }
            rz = next;
        }
        ax = apply(x);
        for (std::size_t i = 0; i < x.size(); ++i) residual[i] = rhs[i] - ax[i];
        std::ostringstream message;
        message << "Flow pressure PCG iteration limit reached: true residual=" << norm(residual)
                << ", target=" << stop << ", rhs=" << norm(rhs);
        throw std::runtime_error(message.str());
    }

    std::size_t solve(Vec& x) const {
        for (double d : diag) {
            ensure(d > 0 && std::isfinite(d), "Flow singular/nonpositive matrix diagonal");
        }
        auto ax = apply(x);
        Vec r(x.size()), r0, p(x.size()), v(x.size()), s(x.size()), t, z(x.size()), zs(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) r[i] = rhs[i] - ax[i];
        r0 = r;
        const double stop = finite(1e-13 + 1e-11 * norm(rhs));
        if (norm(r) <= stop) return 0;
        double rhoOld = 1, alpha = 1, omega = 1;
        for (std::size_t step = 0; step < 3000; ++step) {
            const double rho = product(r0, r);
            ensure(rho != 0 && omega != 0, "Flow BiCGStab breakdown");
            const double beta = (rho / rhoOld) * (alpha / omega);
            for (std::size_t i = 0; i < x.size(); ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
                z[i] = p[i] / diag[i];
            }
            v = apply(z);
            const double rv = product(r0, v);
            ensure(rv != 0, "Flow BiCGStab singular projection");
            alpha = rho / rv;
            for (std::size_t i = 0; i < x.size(); ++i) s[i] = r[i] - alpha * v[i];
            if (norm(s) <= stop) {
                for (std::size_t i = 0; i < x.size(); ++i) x[i] += alpha * z[i];
            } else {
                for (std::size_t i = 0; i < x.size(); ++i) zs[i] = s[i] / diag[i];
                t = apply(zs);
                const double tt = product(t, t);
                ensure(tt > 0, "Flow BiCGStab null update");
                omega = product(t, s) / tt;
                for (std::size_t i = 0; i < x.size(); ++i)
                    x[i] += alpha * z[i] + omega * zs[i];
            }
            ax = apply(x);
            for (std::size_t i = 0; i < x.size(); ++i) r[i] = rhs[i] - ax[i];
            if (norm(r) <= stop) return step + 1;
            // Restart with the true residual periodically; never accept recurrence alone.
            if (step % 40 == 39) {
                r0 = r;
                std::fill(p.begin(), p.end(), 0);
                std::fill(v.begin(), v.end(), 0);
                rhoOld = alpha = omega = 1;
            } else {
                rhoOld = rho;
            }
        }
        throw std::runtime_error("Flow linear solver iteration limit reached");
    }
};

enum class Role { Wall, Inlet, Outlet, Slip, Lid };

struct Boundary {
    std::vector<Role> role;
    Vec u;
    Vec v;
    std::vector<bool> fixedU;
    std::vector<bool> fixedV;
    std::vector<bool> fixedP;
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
            if (c.scenario == "cavity") {
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
            b.fixedV[id] = true;
            break;
        }
    }
    b.closed = c.scenario == "cavity";
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

std::vector<Vector2D> gradient(const FvMesh2D& m,
                               const Vec& u,
                               const Vec& bc,
                               const std::vector<bool>& fixed) {
    std::vector<Vector2D> g(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) {
        double xx = 0;
        double xy = 0;
        double yy = 0;
        double bx = 0;
        double by = 0;
        for (auto id : m.cells[i].faces) {
            const auto& f = m.faces[id];
            const auto j =
                f.owner == i ? f.neighbour : std::optional<std::size_t>(f.owner);
            Vector2D d;
            double value = 0;
            if (j || fixed[id]) {
                d = (j ? m.cells[*j].centre : f.centre) - m.cells[i].centre;
                const double length = std::hypot(d.x, d.y);
                ensure(length > 0, "Flow gradient degenerate stencil");
                value = ((j ? u[*j] : bc[id]) - u[i]) / length;
                d = d * (1 / length);
            } else {
                d = f.areaVector;
                const double length = std::hypot(d.x, d.y);
                d = d * (1 / length);
            }
            xx += d.x * d.x;
            xy += d.x * d.y;
            yy += d.y * d.y;
            bx += d.x * value;
            by += d.y * value;
        }
        const double det = xx * yy - xy * xy;
        ensure(det > 64 * std::numeric_limits<double>::epsilon() * (xx + yy) * (xx + yy),
               "Flow gradient rank deficient");
        g[i] = {finite((yy * bx - xy * by) / det),
                finite((xx * by - xy * bx) / det)};
    }
    return g;
}

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

System momentum(const FvMesh2D& m,
                const FlowControls2D& c,
                const Boundary& b,
                const Vec& field,
                const Vec& flux,
                const std::vector<Vector2D>& gradField,
                const std::vector<Vector2D>& gp,
                bool y,
                bool relaxed) {
    System a(m.cells.size());
    const auto& bc = y ? b.v : b.u;
    const auto& fixed = y ? b.fixedV : b.fixedU;
    for (std::size_t i = 0; i < m.cells.size(); ++i) {
        a.rhs[i] = -m.cells[i].area * (y ? gp[i].y : gp[i].x);
    }
    for (std::size_t id = 0; id < m.faces.size(); ++id) {
        const auto& f = m.faces[id];
        const auto i = f.owner;
        const double q = flux[id];
        const double d = c.nu * f.transmissibility;
        if (f.neighbour) {
            const auto j = *f.neighbour;
            a.diag[i] += d + std::max(q, 0.);
            a.off[i][j] += -d + std::min(q, 0.);
            a.diag[j] += d + std::max(-q, 0.);
            a.off[j][i] += -d - std::max(q, 0.);
            const double correction =
                c.nu * dot(interpolateGradient(f, gradField), f.correction);
            a.rhs[i] += correction;
            a.rhs[j] -= correction;
        } else {
            if (fixed[id]) {
                a.diag[i] += d;
                a.rhs[i] += d * bc[id] + c.nu * dot(gradField[i], f.correction);
                a.rhs[i] -= q * bc[id];
            } else {
                ensure(q >= -1e-12 * c.speed * std::hypot(f.areaVector.x, f.areaVector.y),
                       "Flow outlet backflow unsupported in this laminar prototype");
                a.diag[i] += q;
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
    return a;
}

} // namespace

FlowResult2D solveIncompressible2D(
    const FvMesh2D& m,
    const FlowControls2D& c,
    const std::function<void(const FlowIteration2D&)>& progress) {
    using Clock = std::chrono::steady_clock;
    const auto solveStart = c.profile ? Clock::now() : Clock::time_point{};
    validateFvMesh2D(m);
    ensure((c.scenario == "external" || c.scenario == "channel" || c.scenario == "cavity") &&
               std::isfinite(c.nu) && c.nu > 0 && std::isfinite(c.speed) && c.speed > 0 &&
               std::isfinite(c.tolerance) && c.tolerance > 0 && c.maxIterations > 0,
           "Invalid flow controls");
    ensure(c.velocityRelaxation > 0 && c.velocityRelaxation <= 1 &&
               c.pressureRelaxation > 0 && c.pressureRelaxation <= 1,
           "Invalid SIMPLE relaxation");
    const auto b = boundaries(m, c);
    const auto n = m.cells.size();
    const auto nf = m.faces.size();

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

    FlowResult2D r;
    // Profiling observes the same solves and stopping rules, including zero-step
    // solves. Timing includes each linear solver's setup, but not assembly.
    auto linearSolve = [&](const System& system, Vec& field, bool pressure) {
        const auto start = c.profile ? Clock::now() : Clock::time_point{};
        const auto iterations = pressure ? system.solvePressure(field) : system.solve(field);
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

    for (std::size_t it = 1; it <= c.maxIterations; ++it) {
        const Vec oldU = r.u;
        const Vec oldV = r.v;
        const Vec oldP = r.p;
        const auto gp = gradient(m, r.p, zeros, b.fixedP);
        const auto gu = gradient(m, r.u, b.u, b.fixedU);
        const auto gv = gradient(m, r.v, b.v, b.fixedV);
        auto au = momentum(m, c, b, r.u, r.flux, gu, gp, false, true);
        auto av = momentum(m, c, b, r.v, r.flux, gv, gp, true, true);
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
        const auto gup=gradient(m,r.u,b.u,b.fixedU),gvp=gradient(m,r.v,b.v,b.fixedV);
        Vec predicted(nf);
        for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
            const double rf=interpolate(f,ra);df[id]=rf*f.transmissibility;
            if(f.neighbour){const auto j=*f.neighbour;const double w=f.neighbourWeight;
                const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
                const auto skew=f.centre-point;
                const double uf=interpolate(f,r.u)+dot(interpolateGradient(f,gup),skew),vf=interpolate(f,r.v)+dot(interpolateGradient(f,gvp),skew);
                const Vector2D rag{(1-w)*ra[i]*gp[i].x+w*ra[j]*gp[j].x,(1-w)*ra[i]*gp[i].y+w*ra[j]*gp[j].y};
                predicted[id]=uf*f.areaVector.x+vf*f.areaVector.y+dot(rag,f.areaVector)-rf*(f.transmissibility*(r.p[j]-r.p[i])+dot(interpolateGradient(f,gp),f.correction));
            }else if(b.role[id]==Role::Outlet){predicted[id]=r.u[i]*f.areaVector.x+r.v[i]*f.areaVector.y+ra[i]*dot(gp[i],f.areaVector)-rf*(-f.transmissibility*r.p[i]+dot(gp[i],f.correction));}
            else if(b.role[id]==Role::Inlet)predicted[id]=b.u[id]*f.areaVector.x+b.v[id]*f.areaVector.y;
            else predicted[id]=0;
        }
        std::fill(pc.begin(),pc.end(),0);Vec correction(nf);
        for(int pass=0;pass<4;++pass){
            System ap(n);const auto gc=gradient(m,pc,zeros,b.fixedP);
            for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
                correction[id]=(f.neighbour||b.fixedP[id])?-interpolate(f,ra)*dot(interpolateGradient(f,gc),f.correction):0;
                ap.rhs[i]-=predicted[id]+correction[id];
                if(f.neighbour){const auto j=*f.neighbour;ap.rhs[j]+=predicted[id]+correction[id];ap.diag[i]+=df[id];ap.diag[j]+=df[id];ap.off[i][j]-=df[id];ap.off[j][i]-=df[id];}
                else if(b.fixedP[id])ap.diag[i]+=df[id];
            }
            if(b.closed){ap.pin(0);pc[0]=0;}
            linearSolve(ap, pc, true);
        }
        const auto gc=gradient(m,pc,zeros,b.fixedP);
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
        const auto newGp=gradient(m,r.p,zeros,b.fixedP),newGu=gradient(m,r.u,b.u,b.fixedU),newGv=gradient(m,r.v,b.v,b.fixedV);
        const auto checkU=momentum(m,c,b,r.u,r.flux,newGu,newGp,false,false),checkV=momentum(m,c,b,r.v,r.flux,newGv,newGp,true,false);
        const auto mu=checkU.apply(r.u),mv=checkV.apply(r.v);double mr=0;
        for(std::size_t i=0;i<n;++i){const double scale=finite((checkU.diag[i]+checkV.diag[i])*c.speed);
            ensure(scale>0,"Flow momentum scale underflow");
            mr=std::max(mr,std::hypot(mu[i]-checkU.rhs[i],mv[i]-checkV.rhs[i])/scale);}
        FlowIteration2D step{it,finite(mr),finite(continuity),finite(du),finite(dp)};r.history.push_back(step);
        if(progress&&(it==1||it%10==0))progress(step);
        if(it>=10&&mr<c.tolerance&&du<c.tolerance&&dp<c.tolerance&&continuity<1e-8&&r.globalRelativeImbalance<1e-8){r.converged=true;break;}
    }
    const auto gu=gradient(m,r.u,b.u,b.fixedU),gv=gradient(m,r.v,b.v,b.fixedV),gp=gradient(m,r.p,zeros,b.fixedP);
    for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];if(f.neighbour||f.patch!=BoundaryPatch2D::EmbeddedBoundary||b.role[id]!=Role::Wall)continue;
        const auto i=f.owner;const auto d=f.centre-m.cells[i].centre;
        const double wallP=r.p[i]+dot(gp[i],d);
        r.forceX+=wallP*f.areaVector.x-c.nu*(2*gu[i].x*f.areaVector.x+(gu[i].y+gv[i].x)*f.areaVector.y);
        r.forceY+=wallP*f.areaVector.y-c.nu*((gu[i].y+gv[i].x)*f.areaVector.x+2*gv[i].y*f.areaVector.y);
    }
    finite(r.globalImbalance);finite(r.forceX);finite(r.forceY);
    if (c.profile) {
        r.performance.solveSeconds = std::chrono::duration<double>(Clock::now() - solveStart).count();
    }
    return r;
}
}
