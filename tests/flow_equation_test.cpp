#include "FlowEquation2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <iostream>
#include <limits>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::solver_detail;

namespace {
void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
}

int main() {
    try {
        const auto mesh = fv_test::rectangle(6, 4, 2, true);
        const auto n = mesh.cells.size(), nf = mesh.faces.size();
        Boundary boundary;
        boundary.role.assign(nf, Role::Wall);
        boundary.u.resize(nf); boundary.v.resize(nf); boundary.p.resize(nf);
        boundary.fixedU.assign(nf, true); boundary.fixedV.assign(nf, true);
        boundary.fixedP.assign(nf, false);
        boundary.constantU.assign(nf, false); boundary.constantV.assign(nf, false);
        for (std::size_t id=0; id<nf; ++id) {
            const auto& f = mesh.faces[id];
            if (f.neighbour) continue;
            boundary.u[id] = .5 + .1*f.centre.y;
            boundary.v[id] = .03*f.centre.x;
            if (f.areaVector.x < -.1) boundary.role[id] = Role::Inlet;
            if (f.areaVector.x > .1) {
                boundary.role[id] = Role::Outlet;
                boundary.fixedU[id] = boundary.fixedV[id] = false;
                boundary.fixedP[id] = true;
                boundary.p[id] = .4 + .05*f.centre.y;
            }
        }
        FlowControls2D controls;
        controls.scenario = "custom"; controls.nu = .1;
        controls.convection = ConvectionScheme2D::Upwind;
        std::vector<std::pair<std::size_t, std::size_t>> edges;
        for (const auto& f:mesh.faces) if (f.neighbour) edges.emplace_back(f.owner, *f.neighbour);
        const detail::SparsePattern2D pattern(n, edges);
        System u(pattern), v(pattern), pressure(pattern);
        FlowEquation2D equation(mesh, boundary.fixedP);
        FlowResult2D x, y, combination, zero;
        for (auto* state:{&x, &y, &combination, &zero}) {
            state->u.resize(n); state->v.resize(n); state->p.resize(n); state->flux.resize(nf);
        }
        const double a=.37, b=-.21;
        for (std::size_t i=0; i<n; ++i) {
            const auto& centre = mesh.cells[i].centre;
            x.u[i]=.6+.2*std::sin(centre.x); x.v[i]=.1*centre.y; x.p[i]=.2*centre.x*centre.y;
            y.u[i]=.2*centre.y; y.v[i]=-.1*centre.x; y.p[i]=.3*std::cos(centre.x);
            combination.u[i]=a*x.u[i]+b*y.u[i];
            combination.v[i]=a*x.v[i]+b*y.v[i];
            combination.p[i]=a*x.p[i]+b*y.p[i];
        }
        Vec advective(nf), response(n, .07), coefficients(nf);
        std::vector<Vector2D> sources(n, {.012, -.009});
        for (std::size_t id=0; id<nf; ++id) {
            const auto& f=mesh.faces[id];
            if (f.neighbour || boundary.role[id]==Role::Inlet || boundary.role[id]==Role::Outlet)
                advective[id]=.5*f.areaVector.x;
        }
        // The frozen upwind equation must be affine even with nonzero boundary
        // pressure, velocity and volume sources. This is the algebra assumed by
        // every Krylov action, not a physical-accuracy threshold.
        double maximumAffineError=0;
        for (const auto stress:{ViscousStress2D::Laplacian, ViscousStress2D::Symmetric}) {
            controls.viscousStress=stress;
            Vec rx, ry, rz, rc;
            const auto evaluate=[&](const auto& field, auto& residual) {
                equation.frozenResidual(u, v, field, controls, boundary, boundary.p,
                    advective, response, sources, residual, coefficients, false);
            };
            evaluate(x, rx); evaluate(y, ry); evaluate(zero, rz); evaluate(combination, rc);
            for (std::size_t i=0; i<rc.size(); ++i) {
                const double expected=a*rx[i]+b*ry[i]+(1-a-b)*rz[i];
                maximumAffineError=std::max(maximumAffineError,
                    std::abs(rc[i]-expected)/(1+std::abs(rc[i])+std::abs(expected)));
            }
            const auto reconstructed=equation.reconstruct(x, controls, boundary, boundary.p, advective);
            const auto flux=equation.steadyFlux(x, controls, boundary, reconstructed, response, coefficients);
            double divergenceSum=0, boundaryFlux=0;
            for (std::size_t i=0; i<n; ++i) divergenceSum+=rx[2*n+i];
            for (std::size_t id=0; id<nf; ++id) if (!mesh.faces[id].neighbour) boundaryFlux+=flux[id];
            require(std::abs(divergenceSum-boundaryFlux)<4096*std::numeric_limits<double>::epsilon(),
                "physical divergence omitted a cell or double-counted an internal face");
            assemblePressureBlock2D(pressure, mesh, boundary, coefficients);
            Vec ones(n, 1), product;
            pressure.apply(ones, product);
            Vec expected(n);
            for (std::size_t id=0; id<nf; ++id) if (!mesh.faces[id].neighbour && boundary.fixedP[id])
                expected[mesh.faces[id].owner]+=coefficients[id];
            for (std::size_t i=0; i<n; ++i) require(std::abs(product[i]-expected[i])<4096*std::numeric_limits<double>::epsilon(),
                "pressure block changed zero-gradient or fixed-pressure boundary");
        }
        require(maximumAffineError<4096*std::numeric_limits<double>::epsilon(), "frozen equation is not affine");
        // Changing boundary roles rebuilds only geometry masks; restoring a
        // role and changing a live trace cannot reuse a rejected trial's values.
        const auto original=equation.velocityGradient(x.u, boundary, false);
        auto changed=boundary;
        for (std::size_t id=0; id<nf; ++id) if (!mesh.faces[id].neighbour && boundary.fixedU[id]) {
            changed.fixedU[id]=false; break;
        }
        (void)equation.velocityGradient(y.u, changed, false);
        const auto restored=equation.velocityGradient(x.u, boundary, false);
        for (std::size_t i=0; i<n; ++i) require(original[i].x==restored[i].x && original[i].y==restored[i].y,
            "restored boundary used stale geometry or rejected field values");
        boundary.u[0]+=.7;
        const auto live=equation.velocityGradient(x.u, boundary, false);
        const auto direct=detail::flowGradient(mesh, x.u, boundary.u, boundary.fixedU);
        for (std::size_t i=0; i<n; ++i) require(live[i].x==direct[i].x && live[i].y==direct[i].y,
            "cached geometry froze a boundary value");
        std::cout << "frozen affine residual error=" << maximumAffineError << '\n';
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
