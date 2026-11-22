#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <map>
#include <limits>
namespace cartmesh2d::fv::solver_detail {
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

double faceNu(const FlowControls2D& c, std::size_t id) {
    return c.faceViscosity.empty()?c.nu:c.faceViscosity[id];
}
void validateViscosity(const FvMesh2D& m, const FlowControls2D& c) {
    ensure(c.faceViscosity.empty() || c.faceViscosity.size()==m.faces.size(),
           "Face viscosity dimensions mismatch");
    for (double value:c.faceViscosity)
        ensure(std::isfinite(value)&&value>0,"Face viscosity must be finite positive");
    ensure(std::isfinite(c.manufacturedViscositySlope) &&
           (c.manufacturedViscositySlope==0 || (c.scenario=="manufactured" && !c.faceViscosity.empty() && c.manufacturedViscositySlope>-1)),
           "Manufactured viscosity slope requires manufactured case and positive face field");
    if (!c.faceViscosity.empty()) {
        ensure(c.scenario!="taylor-green" && c.scenario!="counterflow",
               "Verification case does not support prescribed face viscosity");
        if (c.scenario=="manufactured")
            for (std::size_t id=0;id<m.faces.size();++id)
                ensure(c.faceViscosity[id]==c.nu*(1+c.manufacturedViscositySlope*m.faces[id].centre.x),
                       "Manufactured viscosity field differs from analytic definition");
    }
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
                const std::vector<Vector2D>& faceVelocity,
                bool y,
                const Vec* previous,
                double timeStep) {
    a.reset();
    const auto& bc = y ? b.v : b.u;
    const auto& fixed = y ? b.fixedV : b.fixedU;
    const auto limiter = c.convection == ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m, field, gradField, bc, fixed) : Vec{};
    const auto advectiveValue=[&](std::size_t id,double q) {
        return faceVelocity.empty() ? detail::upwindFaceValue(m,id,q,field,gradField,limiter)
                                    : y ? faceVelocity[id].y : faceVelocity[id].x;
    };
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
        const double q = c.momentumInertia * flux[id];
        const double viscosity=faceNu(c,id);
        const double d = viscosity * f.transmissibility;
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
                viscosity * dot(interpolateGradient(f, gradField), f.correction);
            const double upwind = q >= 0 ? field[i] : field[j];
            const double deferred = q * (advectiveValue(id,q)-upwind);
            a.rhs[i] += correction-deferred;
            a.rhs[j] -= correction-deferred;
        } else {
            if (fixed[id]) {
                a.diag[i] += d;
                a.rhs[i] += d * bc[id] + viscosity * dot(gradField[i], f.correction);
                a.rhs[i] -= q * bc[id];
            } else {
                if (q < 0 && (b.role[id] == Role::Farfield ||
                    b.role[id] == Role::Opening ||
                    (b.role[id] == Role::Outlet && c.outletBackflow == OutletBackflow2D::NormalInlet))) {
                    // Same zero-gradient normal advective flux q*Uowner,
                    // lagged in the linear solve to retain a positive diagonal.
                    // The unrelaxed residual uses the current field exactly.
                    a.rhs[i] -= q * field[i];
                    continue;
                }
                ensure(q >= -1e-12 * c.speed * std::hypot(f.areaVector.x, f.areaVector.y),
                       "Flow outlet backflow unsupported in this laminar prototype");
                a.diag[i] += q;
                a.rhs[i] -= q * (advectiveValue(id,q)-field[i]);
            }
        }
    }
}

}
