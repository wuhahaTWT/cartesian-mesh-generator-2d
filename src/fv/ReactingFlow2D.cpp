#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include "cartmesh2d/chemistry/SpeciesMassClosure.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::fv {
namespace {
void require(bool test, const std::string& message) {
    if (!test) throw std::runtime_error("reacting flow: " + message);
}
double finite(double value) { require(std::isfinite(value), "nonfinite arithmetic"); return value; }
double kinetic(const ReactingConservative2D& u) { return .5 * (u[1] * u[1] + u[2] * u[2]) / u[0]; }
ReactingConservative2D flux(const ReactingConservative2D& u, const ReactingPrimitive2D& p, Vector2D n) {
    const double un = dot(p.velocity, n);
    ReactingConservative2D f(u.size());
    for (std::size_t k = 0; k < u.size(); ++k) f[k] = u[k] * un;
    f[1] += p.properties.pressure * n.x; f[2] += p.properties.pressure * n.y;
    f[3] += p.properties.pressure * un;
    return f;
}
std::vector<double> primitiveValues(const ReactingPrimitive2D& p) {
    std::vector<double> v{p.properties.temperature, std::log(p.properties.pressure), p.velocity.x, p.velocity.y};
    v.insert(v.end(), p.properties.massFractions.begin(), p.properties.massFractions.end());
    return v;
}
ReactingConservative2D fromValues(chemistry::DetailedGas& gas, const std::vector<double>& v) {
    std::vector<double> y(v.begin() + 4, v.end());
    (void)chemistry::closeSpeciesMassRoundoff(1, y);
    return reactingConservative2D(gas.fromMassFractions(v[0], std::exp(v[1]), y), {v[2], v[3]});
}
bool wall(ReactingBoundaryKind2D k) { return k == ReactingBoundaryKind2D::SlipWall || k == ReactingBoundaryKind2D::NoSlipWall; }
}

ReactingConservative2D reactingConservative2D(const chemistry::GasState& q, Vector2D u) {
    require(q.density > 0 && std::isfinite(q.density), "invalid density");
    finite(u.x); finite(u.y);
    ReactingConservative2D result{q.density, q.density * u.x, q.density * u.y,
        q.internalEnergyDensity + .5 * q.density * (u.x * u.x + u.y * u.y)};
    result.insert(result.end(), q.speciesDensities.begin(), q.speciesDensities.end());
    for (double v : result) finite(v);
    return result;
}
ReactingPrimitive2D reactingPrimitive2D(chemistry::DetailedGas& gas, const ReactingConservative2D& u) {
    require(u.size() == gas.mechanism().species.size() + 4, "conservative species count mismatch");
    for (double v : u) finite(v);
    require(u[0] > 0, "nonpositive density");
    ReactingPrimitive2D p;
    p.gas = {u[0], u[3] - kinetic(u), {u.begin() + 4, u.end()}};
    p.properties = gas.properties(p.gas);
    p.velocity = {u[1] / u[0], u[2] / u[0]};
    p.soundSpeed = finite(std::sqrt(p.properties.cp / p.properties.cv * p.properties.pressure / u[0]));
    require(p.soundSpeed > 0, "nonpositive acoustic speed");
    return p;
}

ReactingFaceFlux2D reactingFaceFlux2D(chemistry::DetailedGas& gas,
    const ReactingConservative2D& ul, const ReactingConservative2D& ur, Vector2D s) {
    const auto l = reactingPrimitive2D(gas, ul), r = reactingPrimitive2D(gas, ur);
    const double length = finite(std::hypot(s.x, s.y)); require(length > 0, "zero face area vector");
    const Vector2D n{s.x / length, s.y / length};
    const double vl = dot(l.velocity, n), vr = dot(r.velocity, n);
    const double sl = std::min(vl - l.soundSpeed, vr - r.soundSpeed);
    const double sr = std::max(vl + l.soundSpeed, vr + r.soundSpeed);
    const auto fl = flux(ul, l, n), fr = flux(ur, r, n);
    ReactingFaceFlux2D out; out.integratedFlux.resize(ul.size());
    out.waveSpeed = std::max(std::abs(sl), std::abs(sr));
    if (sl >= 0) out.integratedFlux = fl;
    else if (sr <= 0) out.integratedFlux = fr;
    else {
        const double dl = ul[0] * (sl - vl), dr = ur[0] * (sr - vr);
        const double sm = (r.properties.pressure - l.properties.pressure + dl * vl - dr * vr) / (dl - dr);
        bool usable = std::isfinite(sm) && sl < sm && sm < sr;
        ReactingConservative2D star[2];
        double starPressure[2]{};
        for (unsigned side = 0; side < 2 && usable; ++side) {
            const auto& u = side ? ur : ul; const auto& p = side ? r : l;
            const double speed = side ? sr : sl, vn = side ? vr : vl;
            const double ratio = (speed - vn) / (speed - sm);
            const double ps = p.properties.pressure + u[0] * (speed - vn) * (sm - vn);
            starPressure[side] = ps;
            auto& us = star[side]; us.resize(u.size());
            us[0] = u[0] * ratio;
            us[1] = us[0] * (p.velocity.x + (sm - vn) * n.x);
            us[2] = us[0] * (p.velocity.y + (sm - vn) * n.y);
            us[3] = ((speed - vn) * u[3] - p.properties.pressure * vn + ps * sm) / (speed - sm);
            for (std::size_t k = 4; k < u.size(); ++k) us[k] = u[k] * ratio;
            usable = std::isfinite(ps) && ps > 0 && ratio > 0;
            try { if (usable) (void)reactingPrimitive2D(gas, us); } catch (const std::exception&) { usable = false; }
        }
        if (usable) {
            const auto side = sm >= 0 ? 0U : 1U;
            const auto& u = star[side]; const auto& p = side ? r : l;
            // Algebraically the HLLC jump flux, evaluated in the star region.
            // Avoid subtracting nearly equal densities separately for every
            // species when the normal mass flow is close to zero.
            out.integratedFlux[0] = u[0] * sm;
            out.integratedFlux[1] = u[1] * sm + starPressure[side] * n.x;
            out.integratedFlux[2] = u[2] * sm + starPressure[side] * n.y;
            out.integratedFlux[3] = (u[3] + starPressure[side]) * sm;
            for (std::size_t k = 4; k < ul.size(); ++k)
                out.integratedFlux[k] = out.integratedFlux[0] * p.properties.massFractions[k - 4];
        } else for (std::size_t k = 0; k < ul.size(); ++k)
            out.integratedFlux[k] = (sr * fl[k] - sl * fr[k] + sl * sr * (ur[k] - ul[k])) / (sr - sl);
        out.hlleFallback = !usable;
    }
    for (double& v : out.integratedFlux) v = finite(v * length);
    finite(out.waveSpeed);
    return out;
}

ReactingFlowStepper2D::ReactingFlowStepper2D(chemistry::DetailedGas& gas, FvMesh2D mesh,
    std::vector<ReactingBoundary2D> boundaries, ReactingPhysics2D physics)
    : gas_(gas), mesh_(std::move(mesh)), physics_(physics) {
    mechanismDefinition_ = gas.mechanism().resolvedDefinition;
    backendVersion_ = gas.mechanism().backendVersion;
    validateFvMesh2D(mesh_);
    boundaries_.resize(mesh_.faces.size());
    std::vector<ReactingDiffusionBoundary2D> diff;
    std::vector<ViscousBoundary2D> visc;
    for (auto& b : boundaries) {
        require(b.face < boundaries_.size() && !mesh_.faces[b.face].neighbour && !boundaries_[b.face], "invalid/duplicate boundary");
        require(wall(b.kind) || b.kind == ReactingBoundaryKind2D::Reservoir || b.kind == ReactingBoundaryKind2D::ExtrapolatedOutflow, "unknown boundary kind");
        finite(b.velocity.x); finite(b.velocity.y); finite(b.wallTemperature);
        require((b.kind == ReactingBoundaryKind2D::Reservoir) == b.reservoir.has_value(), "reservoir state missing/inactive");
        require(wall(b.kind) ? b.wallTemperature >= 0 : b.wallTemperature == 0, "invalid/inactive wall temperature");
        require(b.wallTemperature == 0 || (b.wallTemperature >= gas.mechanism().minimumTemperature
                && b.wallTemperature <= gas.mechanism().maximumTemperature), "wall temperature outside mechanism range");
        if (b.reservoir) (void)gas.properties(*b.reservoir);
        const auto& f = mesh_.faces[b.face]; const double length = std::hypot(f.areaVector.x, f.areaVector.y);
        if (wall(b.kind)) require(std::abs(dot(b.velocity, f.areaVector)) <= 64 * std::numeric_limits<double>::epsilon() * length * std::hypot(b.velocity.x, b.velocity.y), "moving normal wall on a static mesh");
        if (b.kind == ReactingBoundaryKind2D::SlipWall || b.kind == ReactingBoundaryKind2D::ExtrapolatedOutflow)
            require(b.velocity.x == 0 && b.velocity.y == 0, "inactive boundary velocity");
        ReactingDiffusionBoundary2D db; db.face = b.face;
        if (b.reservoir) { db.kind = ReactingDiffusionBoundaryKind2D::Reservoir; db.reservoir = b.reservoir; }
        else if (b.kind == ReactingBoundaryKind2D::ExtrapolatedOutflow) db.kind = ReactingDiffusionBoundaryKind2D::ZeroFluxOutflow;
        else if (b.wallTemperature > 0) { db.kind = ReactingDiffusionBoundaryKind2D::IsothermalWall; db.temperature = b.wallTemperature; }
        diff.push_back(db);
        ViscousBoundary2D vb; vb.face = b.face;
        if (b.kind == ReactingBoundaryKind2D::SlipWall) vb.kind = ViscousBoundaryKind2D::Slip;
        else if (b.kind == ReactingBoundaryKind2D::ExtrapolatedOutflow) vb.kind = ViscousBoundaryKind2D::ZeroTraction;
        else { vb.kind = ViscousBoundaryKind2D::Velocity; vb.velocity = b.velocity; }
        visc.push_back(vb);
        boundaries_[b.face] = std::move(b);
    }
    for (std::size_t f = 0; f < mesh_.faces.size(); ++f) require(mesh_.faces[f].neighbour || boundaries_[f], "missing flow boundary");
    if (physics_.molecularTransport) {
        diffusion_.emplace(mesh_, std::move(diff));
        viscous_.emplace(mesh_, std::move(visc), 1.0);
    }
    std::ostringstream binding; binding.imbue(std::locale::classic()); binding << std::setprecision(17);
    binding << "CM2D_REACTING_BINDING 1\n" << gas.mechanism().backendVersion << '\n' << gas.mechanism().resolvedDefinition
            << "\nPHYSICS " << physics.chemistry << ' ' << physics.molecularTransport << "\nCELLS " << mesh_.cells.size() << '\n';
    for (const auto& c : mesh_.cells) { binding << c.centre.x << ' ' << c.centre.y << ' ' << c.area << ' ' << c.faces.size(); for (auto f : c.faces) binding << ' ' << f; binding << '\n'; }
    binding << "FACES " << mesh_.faces.size() << '\n';
    for (std::size_t id = 0; id < mesh_.faces.size(); ++id) {
        const auto& f = mesh_.faces[id];
        binding << f.owner << ' '; if (f.neighbour) binding << *f.neighbour; else binding << '-';
        binding << ' ' << static_cast<int>(f.patch) << ' ' << f.centre.x << ' ' << f.centre.y << ' ' << f.areaVector.x << ' ' << f.areaVector.y
                << ' ' << f.transmissibility << ' ' << f.neighbourWeight << ' ' << f.correction.x << ' ' << f.correction.y << '\n';
        if (boundaries_[id]) {
            const auto& b = *boundaries_[id]; binding << "BC " << static_cast<int>(b.kind) << ' ' << b.velocity.x << ' ' << b.velocity.y << ' ' << b.wallTemperature;
            if (b.reservoir) { binding << ' ' << b.reservoir->density << ' ' << b.reservoir->internalEnergyDensity; for (double v : b.reservoir->speciesDensities) binding << ' ' << v; }
            binding << '\n';
        }
    }
    binding_ = binding.str();
}

void ReactingFlowStepper2D::validate(const ReactingState2D& state) const {
    require(gas_.mechanism().resolvedDefinition == mechanismDefinition_ && gas_.mechanism().backendVersion == backendVersion_,
            "chemistry context was replaced after solver construction");
    require(state.binding == binding_, "incompatible resolved mechanism, mesh, boundary or physical model");
    require(std::isfinite(state.time) && state.time >= 0 && state.steps < std::numeric_limits<std::size_t>::max(), "invalid clock/step counter");
    require(state.cells.size() == mesh_.cells.size(), "cell count mismatch");
    for (const auto& u : state.cells) (void)reactingPrimitive2D(gas_, u);
}
ReactingState2D ReactingFlowStepper2D::initialState(std::vector<ReactingConservative2D> cells) const {
    ReactingState2D s; s.binding = binding_; s.cells = std::move(cells); validate(s); return s;
}

ReactingResidual2D ReactingFlowStepper2D::evaluateResidual(const ReactingState2D& state, unsigned order) {
    validate(state);
    require(order == 1 || order == 2, "invalid residual reconstruction order");
    auto stage = spatial(state.cells, order);
    const auto nv = gas_.mechanism().species.size() + 4;
    ReactingResidual2D out;
    out.faceFlux = std::move(stage.faceFlux); out.transportRate = std::move(stage.rate);
    out.hlleFallbacks = stage.fallbacks;
    out.transportDerivative = std::move(stage.residual);
    out.chemistryDerivative.assign(state.cells.size(), ReactingConservative2D(nv));
    for (std::size_t i = 0; i < state.cells.size(); ++i) {
        for (double& value : out.transportDerivative[i]) value = finite(-value / mesh_.cells[i].area);
        if (physics_.chemistry) {
            const auto p = reactingPrimitive2D(gas_, state.cells[i]);
            for (std::size_t k = 4; k < nv; ++k) out.chemistryDerivative[i][k] = p.properties.massProductionRates[k - 4];
        }
    }
    out.derivative = out.transportDerivative;
    for (std::size_t i = 0; i < state.cells.size(); ++i)
        for (std::size_t k = 0; k < nv; ++k) out.derivative[i][k] = finite(out.derivative[i][k] + out.chemistryDerivative[i][k]);
    out.chemistryIntegral = integral(out.chemistryDerivative);
    out.boundaryFlux.assign(nv, 0);
    for (std::size_t f = 0; f < mesh_.faces.size(); ++f) if (!mesh_.faces[f].neighbour)
        for (std::size_t k = 0; k < nv; ++k) out.boundaryFlux[k] += out.faceFlux[f][k];
    return out;
}

ReactingFlowStepper2D::Stage ReactingFlowStepper2D::spatial(const std::vector<ReactingConservative2D>& cells, unsigned order) {
    const auto nc = cells.size(), nf = mesh_.faces.size(), nv = gas_.mechanism().species.size() + 4;
    Stage out; out.faceFlux.assign(nf, ReactingConservative2D(nv)); out.residual.assign(nc, ReactingConservative2D(nv)); out.rate.resize(nc);
    std::vector<ReactingPrimitive2D> primitive;
    std::vector<std::vector<double>> phi, boundaryPhi(nf);
    for (const auto& u : cells) { primitive.push_back(reactingPrimitive2D(gas_, u)); phi.push_back(primitiveValues(primitive.back())); }
    for (std::size_t f = 0; f < nf; ++f) if (boundaries_[f]) {
        const auto& b = *boundaries_[f]; const auto owner = mesh_.faces[f].owner;
        boundaryPhi[f] = b.reservoir ? primitiveValues(reactingPrimitive2D(gas_, reactingConservative2D(*b.reservoir, b.velocity))) : phi[owner];
        if (wall(b.kind)) {
            const auto& s = mesh_.faces[f].areaVector; const double length = std::hypot(s.x, s.y); const Vector2D n{s.x / length, s.y / length};
            const double un = dot(primitive[owner].velocity, n);
            boundaryPhi[f][2] = b.kind == ReactingBoundaryKind2D::NoSlipWall ? b.velocity.x : phi[owner][2] - un * n.x;
            boundaryPhi[f][3] = b.kind == ReactingBoundaryKind2D::NoSlipWall ? b.velocity.y : phi[owner][3] - un * n.y;
            if (b.wallTemperature > 0) boundaryPhi[f][0] = b.wallTemperature;
        }
    }
    std::vector<std::vector<Vector2D>> gradient(nv, std::vector<Vector2D>(nc));
    if (order == 2) {
        for (std::size_t k = 0; k < nv; ++k) {
            std::vector<double> v(nc), bc(nf);
            for (std::size_t i = 0; i < nc; ++i) v[i] = phi[i][k];
            for (std::size_t f = 0; f < nf; ++f) if (boundaries_[f]) bc[f] = boundaryPhi[f][k];
            gradient[k] = reconstructGradient(mesh_, v, bc);
        }
        for (std::size_t i = 0; i < nc; ++i) {
            const auto dependent = static_cast<std::size_t>(std::max_element(phi[i].begin() + 4, phi[i].end()) - phi[i].begin());
            gradient[dependent][i] = {};
            for (std::size_t k = 4; k < nv; ++k) if (k != dependent) {
                gradient[dependent][i].x -= gradient[k][i].x; gradient[dependent][i].y -= gradient[k][i].y;
            }
            std::vector<double> theta(nv, 1);
            for (std::size_t k = 0; k < nv; ++k) {
                double lo = phi[i][k], hi = lo;
                for (auto id : mesh_.cells[i].faces) {
                    const auto& face = mesh_.faces[id]; const auto other = face.owner == i ? face.neighbour : std::optional<std::size_t>(face.owner);
                    const double value = other ? phi[*other][k] : boundaryPhi[id][k]; lo = std::min(lo, value); hi = std::max(hi, value);
                }
                for (auto id : mesh_.cells[i].faces) {
                    const double delta = dot(gradient[k][i], mesh_.faces[id].centre - mesh_.cells[i].centre);
                    if (delta > 0) theta[k] = std::min(theta[k], (hi - phi[i][k]) / delta);
                    if (delta < 0) theta[k] = std::min(theta[k], (lo - phi[i][k]) / delta);
                }
            }
            const double speciesTheta = *std::min_element(theta.begin() + 4, theta.end());
            // Keep reconstructed face values inside their convex bounds after
            // multiply/add roundoff. Without this arithmetic reserve a trace
            // value exactly limited to zero can become, e.g., -2.5e-60.
            // One common species factor preserves sum(grad Y)=0. No face or
            // conserved cell value is clipped and no positive trace is erased.
            const double inward = 1 - 64 * std::numeric_limits<double>::epsilon();
            for (std::size_t k = 0; k < nv; ++k)
                gradient[k][i] = gradient[k][i] * (inward * std::clamp(k < 4 ? theta[k] : speciesTheta, 0., 1.));
        }
    }
    const auto reconstructed = [&](std::size_t i, Point2D point) {
        if (order == 1) return cells[i];
        auto values = phi[i];
        for (std::size_t k = 0; k < nv; ++k) values[k] += dot(gradient[k][i], point - mesh_.cells[i].centre);
        try { return fromValues(gas_, values); }
        catch (const std::exception& error) {
            std::ostringstream message; message << std::setprecision(17) << "MUSCL face from cell " << i << ": " << error.what();
            for (std::size_t k = 4; k < nv; ++k) if (values[k] < 0)
                message << "; " << gas_.mechanism().species[k - 4] << "=" << values[k];
            throw std::runtime_error(message.str());
        }
    };
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = mesh_.faces[id]; const auto* b = boundaries_[id] ? &*boundaries_[id] : nullptr;
        const double length = std::hypot(f.areaVector.x, f.areaVector.y); const Vector2D n{f.areaVector.x / length, f.areaVector.y / length};
        const auto left = reconstructed(f.owner, f.centre);
        auto right = f.neighbour ? reconstructed(*f.neighbour, f.centre) : left;
        if (b && b->reservoir) right = reactingConservative2D(*b->reservoir, b->velocity);
        if (b && wall(b->kind)) {
            const auto p = reactingPrimitive2D(gas_, left); const double un = dot(p.velocity, n);
            right = reactingConservative2D(p.gas, {p.velocity.x - 2 * un * n.x, p.velocity.y - 2 * un * n.y});
        }
        if (b && b->kind == ReactingBoundaryKind2D::ExtrapolatedOutflow)
            require(left[1] * n.x + left[2] * n.y >= 0, "outflow backflow needs a specified reservoir state");
        auto numerical = reactingFaceFlux2D(gas_, left, right, f.areaVector);
        if (b && wall(b->kind)) {
            const double pressureForce = numerical.integratedFlux[1] * n.x + numerical.integratedFlux[2] * n.y;
            numerical.integratedFlux.assign(nv, 0); numerical.integratedFlux[1] = pressureForce * n.x; numerical.integratedFlux[2] = pressureForce * n.y;
        }
        out.faceFlux[id] = std::move(numerical.integratedFlux); out.fallbacks += numerical.hlleFallback ? 1 : 0;
        out.rate[f.owner] += numerical.waveSpeed * length / mesh_.cells[f.owner].area;
        if (f.neighbour) out.rate[*f.neighbour] += numerical.waveSpeed * length / mesh_.cells[*f.neighbour].area;
    }
    if (diffusion_) {
        std::vector<chemistry::GasState> state; std::vector<Vector2D> velocity; std::vector<double> density;
        for (const auto& p : primitive) { state.push_back(p.gas); velocity.push_back(p.velocity); density.push_back(p.gas.density); }
        const auto d = diffusion_->evaluate(gas_, state, true);
        std::vector<double> mu(nf);
        for (std::size_t id = 0; id < nf; ++id) {
            const auto& f = mesh_.faces[id]; const auto* b = boundaries_[id] ? &*boundaries_[id] : nullptr;
            const auto& lp = primitive[f.owner].properties;
            auto rp = f.neighbour ? primitive[*f.neighbour].properties : b->reservoir ? gas_.properties(*b->reservoir) : lp;
            if (b && b->wallTemperature > 0) rp.temperature = b->wallTemperature;
            const double w = f.neighbour ? f.neighbourWeight : 1;
            auto y = lp.massFractions; for (std::size_t k = 0; k < y.size(); ++k) y[k] = (1 - w) * y[k] + w * rp.massFractions[k];
            (void)chemistry::closeSpeciesMassRoundoff(1, y);
            const auto q = gas_.fromMassFractions((1 - w) * lp.temperature + w * rp.temperature,
                std::exp((1 - w) * std::log(lp.pressure) + w * std::log(rp.pressure)), y);
            mu[id] = gas_.viscosity(q);
            out.faceFlux[id][3] += d.faces[id].energy;
            for (std::size_t k = 4; k < nv; ++k) out.faceFlux[id][k] += d.faces[id].species[k - 4];
        }
        const auto v = viscous_->evaluate(velocity, density, mu);
        for (std::size_t f = 0; f < nf; ++f) for (std::size_t k = 0; k < 3; ++k) out.faceFlux[f][k + 1] += v.faceFlux[f][k];
        for (std::size_t i = 0; i < nc; ++i) out.rate[i] += d.rate[i] + v.rate[i];
    }
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = mesh_.faces[id];
        for (std::size_t k = 0; k < nv; ++k) { const double value = finite(out.faceFlux[id][k]); out.residual[f.owner][k] += value; if (f.neighbour) out.residual[*f.neighbour][k] -= value; }
    }
    for (double v : out.rate) require(std::isfinite(v) && v > 0, "invalid transport rate");
    return out;
}

ReactingConservative2D ReactingFlowStepper2D::integral(const std::vector<ReactingConservative2D>& cells) const {
    std::vector<long double> sums(gas_.mechanism().species.size() + 4);
    for (std::size_t i = 0; i < cells.size(); ++i) for (std::size_t k = 0; k < sums.size(); ++k) sums[k] += static_cast<long double>(cells[i][k]) * mesh_.cells[i].area;
    return {sums.begin(), sums.end()};
}
void ReactingFlowStepper2D::react(std::vector<ReactingConservative2D>& cells, double dt,
    const chemistry::ChemistryControls& controls, ReactingStepResult2D& result, ReactingConservative2D& change,
    ReactingConservative2D& closureChange, double& maximumClosure, double& absoluteClosure) {
    if (!physics_.chemistry) return;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto before = cells[i]; const auto p = reactingPrimitive2D(gas_, before);
        ++result.sourceCalls;
        const auto source = gas_.advanceConstantVolume(p.gas, dt, controls);
        require(source.accepted.has_value(), "cell " + std::to_string(i) + " chemistry: " + source.failure);
        result.successfulSourceInternalSteps += source.internalSteps;
        auto& u = cells[i]; u[0] = source.accepted->density;
        u[3] = source.accepted->internalEnergyDensity + kinetic(u);
        for (std::size_t k = 4; k < u.size(); ++k) u[k] = source.accepted->speciesDensities[k - 4];
        for (std::size_t k = 4; k < u.size(); ++k) {
            closureChange[k] += source.massClosureChange[k - 4] * mesh_.cells[i].area;
            maximumClosure = std::max(maximumClosure, std::abs(source.massClosureChange[k - 4] / u[0]));
            absoluteClosure += std::abs(source.massClosureChange[k - 4]) * mesh_.cells[i].area;
        }
        for (std::size_t k = 0; k < u.size(); ++k) change[k] += (u[k] - before[k]) * mesh_.cells[i].area;
        (void)reactingPrimitive2D(gas_, u);
    }
}

ReactingStepResult2D ReactingFlowStepper2D::advance(const ReactingState2D& initial, const ReactingStepControls2D& c) {
    ReactingStepResult2D result;
    try {
        validate(initial);
        require(std::isfinite(c.maximumStep) && std::isfinite(c.minimumStep) && c.minimumStep > 0 && c.maximumStep >= c.minimumStep,
                "invalid time-step bounds");
        require(std::isfinite(c.courant) && c.courant > 0 && c.courant <= .45 && c.maximumRetries <= 30 && (c.order == 1 || c.order == 2), "invalid time/space controls");
        require(std::isfinite(c.predictionSafety) && c.predictionSafety > 0 && c.predictionSafety <= 1, "invalid step prediction safety");
        require(!c.endTime || (std::isfinite(*c.endTime) && *c.endTime > initial.time), "end time is not ahead of accepted clock");
        const auto nc = initial.cells.size(), nv = gas_.mechanism().species.size() + 4;
        const auto estimate = spatial(initial.cells, c.order);
        double dt = std::min(c.maximumStep, c.predictionSafety * c.courant / *std::max_element(estimate.rate.begin(), estimate.rate.end()));
        if (c.endTime) dt = std::min(dt, *c.endTime - initial.time);
        result.beforeIntegral = integral(initial.cells);
        for (std::size_t attempt = 0; attempt <= c.maximumRetries; ++attempt) {
            require(dt >= c.minimumStep && std::isfinite(dt) && initial.time + dt > initial.time, "required step below minimum or clock resolution");
            auto cells = initial.cells; ReactingConservative2D chemicalChange(nv), chemicalClosure(nv), transportClosure(nv);
            double maximumClosure = 0, absoluteClosure = 0;
            const auto closeStage = [&](std::size_t i, double weight) {
                const auto closure = chemistry::closeSpeciesMassRoundoff(cells[i][0], cells[i], 4);
                transportClosure[closure.species + 4] += weight * closure.change * mesh_.cells[i].area;
                maximumClosure = std::max(maximumClosure, std::abs(closure.change / cells[i][0]));
                absoluteClosure += weight * std::abs(closure.change) * mesh_.cells[i].area;
            };
            try {
                react(cells, dt / 2, c.chemistry, result, chemicalChange, chemicalClosure, maximumClosure, absoluteClosure);
                const auto base = cells;
                auto first = spatial(base, c.order);
                double courant = dt * *std::max_element(first.rate.begin(), first.rate.end());
                require(courant <= c.courant * (1 + 16 * std::numeric_limits<double>::epsilon()), "post-chemistry transport CFL increased");
                for (std::size_t i = 0; i < nc; ++i) {
                    for (std::size_t k = 0; k < nv; ++k) cells[i][k] -= dt / mesh_.cells[i].area * first.residual[i][k];
                    try { closeStage(i, c.order == 2 ? .5 : 1); (void)reactingPrimitive2D(gas_, cells[i]); }
                    catch (const std::exception& e) {
                        std::ostringstream message; message << std::setprecision(17) << "first transport stage, cell " << i << ": " << e.what();
                        for (std::size_t k = 4; k < nv; ++k) if (cells[i][k] < 0)
                            message << "; rhoY(" << gas_.mechanism().species[k - 4] << ")=" << cells[i][k]
                                    << ", base=" << base[i][k] << ", derivative=" << -first.residual[i][k] / mesh_.cells[i].area;
                        throw std::runtime_error(message.str());
                    }
                }
                if (c.order == 2) {
                    const auto second = spatial(cells, c.order);
                    courant = std::max(courant, dt * *std::max_element(second.rate.begin(), second.rate.end()));
                    require(courant <= c.courant * (1 + 16 * std::numeric_limits<double>::epsilon()), "second-stage transport CFL increased");
                    for (std::size_t i = 0; i < nc; ++i) {
                        for (std::size_t k = 0; k < nv; ++k) cells[i][k] = .5 * (base[i][k] + cells[i][k] - dt / mesh_.cells[i].area * second.residual[i][k]);
                        try { closeStage(i, 1); (void)reactingPrimitive2D(gas_, cells[i]); }
                        catch (const std::exception& e) {
                            std::ostringstream message; message << std::setprecision(17) << "second transport stage, cell " << i << ": " << e.what();
                            for (std::size_t k = 4; k < nv; ++k) if (cells[i][k] < 0)
                                message << "; rhoY(" << gas_.mechanism().species[k - 4] << ")=" << cells[i][k];
                            throw std::runtime_error(message.str());
                        }
                    }
                    for (std::size_t f = 0; f < first.faceFlux.size(); ++f) for (std::size_t k = 0; k < nv; ++k) first.faceFlux[f][k] = .5 * (first.faceFlux[f][k] + second.faceFlux[f][k]);
                    first.fallbacks += second.fallbacks;
                }
                react(cells, dt / 2, c.chemistry, result, chemicalChange, chemicalClosure, maximumClosure, absoluteClosure);
                ReactingState2D accepted{initial.time + dt, initial.steps + 1, binding_, std::move(cells)};
                validate(accepted);
                require(!c.endTime || accepted.time <= *c.endTime, "candidate exceeds integration horizon");
                result.step = dt; result.combinedCourant = courant; result.hlleFallbacks = first.fallbacks;
                result.faceFlux = std::move(first.faceFlux); result.chemistryChange = std::move(chemicalChange);
                result.transportMassClosureChange = std::move(transportClosure);
                result.chemistryMassClosureChange = std::move(chemicalClosure);
                result.maximumMassClosureFraction = maximumClosure;
                result.absoluteMassClosureIntegral = absoluteClosure;
                result.afterIntegral = integral(accepted.cells); result.boundaryFlux.assign(nv, 0); result.balanceError.resize(nv);
                for (std::size_t f = 0; f < mesh_.faces.size(); ++f) if (!mesh_.faces[f].neighbour)
                    for (std::size_t k = 0; k < nv; ++k) result.boundaryFlux[k] += result.faceFlux[f][k];
                for (std::size_t k = 0; k < nv; ++k) result.balanceError[k] = finite(result.afterIntegral[k] - result.beforeIntegral[k]
                    + dt * result.boundaryFlux[k] - result.chemistryChange[k] - result.transportMassClosureChange[k]);
                const auto& m = gas_.mechanism(); result.elementalBalanceError.resize(m.elements.size());
                for (std::size_t e = 0; e < m.elements.size(); ++e) {
                    long double sum = 0;
                    for (std::size_t k = 0; k < m.species.size(); ++k) sum += (result.afterIntegral[k + 4] - result.beforeIntegral[k + 4] + dt * result.boundaryFlux[k + 4])
                        * static_cast<long double>(m.atomicWeights[e] * m.atomCounts[k * m.elements.size() + e] / m.molecularWeights[k]);
                    result.elementalBalanceError[e] = finite(static_cast<double>(sum / result.beforeIntegral[0]));
                }
                result.accepted = std::move(accepted); return result;
            } catch (const std::exception& e) {
                result.rejectedReasons.push_back(e.what());
                if (attempt == c.maximumRetries) throw std::runtime_error("reacting flow retry budget exhausted: " + std::string(e.what()));
                dt *= .5;
            }
        }
    } catch (const std::exception& e) { result.failure = e.what(); }
    return result;
}

void ReactingFlowStepper2D::writeCheckpoint(std::ostream& out, const ReactingState2D& state) const {
    validate(state); require(state.steps > 0 && state.time > 0, "initial condition is not an accepted physical checkpoint");
    std::ostringstream data; data.imbue(std::locale::classic()); data << std::setprecision(17);
    data << "CM2D_REACTING_CHECKPOINT 1\n" << binding_.size() << '\n' << binding_ << "STATE " << state.time << ' ' << state.steps << '\n';
    for (const auto& u : state.cells) { for (double v : u) data << v << ' '; data << '\n'; }
    data << "END\n"; out << data.str(); require(static_cast<bool>(out), "checkpoint write failed");
}
ReactingState2D ReactingFlowStepper2D::readCheckpoint(std::istream& in) const {
    in.imbue(std::locale::classic());
    std::string line; require(static_cast<bool>(std::getline(in, line)) && line == "CM2D_REACTING_CHECKPOINT 1", "invalid checkpoint header");
    require(static_cast<bool>(std::getline(in, line)) && line == std::to_string(binding_.size()), "incompatible checkpoint binding size");
    std::string binding(binding_.size(), '\0'); in.read(binding.data(), static_cast<std::streamsize>(binding.size()));
    require(static_cast<bool>(in) && binding == binding_, "checkpoint mechanism/geometry/boundary/model mismatch");
    ReactingState2D s; s.binding = binding_; std::string token, steps;
    require(static_cast<bool>(in >> token >> s.time >> steps) && token == "STATE" && !steps.empty()
        && steps.find_first_not_of("0123456789") == std::string::npos, "invalid checkpoint clock");
    const auto count = std::stoull(steps); require(count > 0 && count < std::numeric_limits<std::size_t>::max() && s.time > 0, "invalid accepted step count/time");
    s.steps = static_cast<std::size_t>(count); s.cells.assign(mesh_.cells.size(), ReactingConservative2D(gas_.mechanism().species.size() + 4));
    for (auto& u : s.cells) for (double& v : u) require(static_cast<bool>(in >> v), "truncated checkpoint state");
    require(static_cast<bool>(in >> token) && token == "END" && !(in >> token), "invalid checkpoint terminator/trailing content");
    validate(s); return s;
}
} // namespace cartmesh2d::fv
