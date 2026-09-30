#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::chemistry;
namespace {
void check(bool condition, const char* text) { if (!condition) throw std::runtime_error(text); }
template<class Action> void rejects(Action action, const char* text) {
    bool failed = false; try { action(); } catch (const std::exception&) { failed = true; } check(failed, text);
}
std::vector<double> composition(const DetailedGas& gas, std::initializer_list<std::pair<std::string,double>> input) {
    std::vector<double> result(gas.mechanism().species.size());
    for (const auto& [name, amount] : input) {
        const auto it = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), name);
        check(it != gas.mechanism().species.end(), "test species missing");
        result[static_cast<std::size_t>(it - gas.mechanism().species.begin())] = amount;
    }
    return result;
}
FvMesh2D scaledMesh(int nx, int ny, double size) {
    auto m = fv_test::rectangle(nx, ny, 1, true);
    for (auto& c : m.cells) { c.centre.x *= size; c.centre.y *= size; c.area *= size * size; }
    for (auto& f : m.faces) { f.centre.x *= size; f.centre.y *= size; f.areaVector = f.areaVector * size; f.correction = f.correction * size; }
    return m;
}
std::vector<ReactingBoundary2D> sealed(const FvMesh2D& m) {
    std::vector<ReactingBoundary2D> b;
    for (std::size_t f = 0; f < m.faces.size(); ++f) if (!m.faces[f].neighbour) b.push_back({f, ReactingBoundaryKind2D::SlipWall, {}, {}, 0});
    return b;
}
void accepted(const ReactingStepResult2D& step) { if (!step.accepted) throw std::runtime_error(step.failure); }
double conservation(const ReactingStepResult2D& s) {
    double error = 0;
    for (std::size_t k = 0; k < s.balanceError.size(); ++k) {
        // Keep each equation's own units. Gross face activity supplies a
        // roundoff scale even when opposite momentum fluxes cancel globally.
        double grossTransport = 0;
        for (const auto& flux : s.faceFlux) grossTransport += s.step * std::abs(flux[k]);
        const double scale = std::max({std::abs(s.beforeIntegral[k]), std::abs(s.afterIntegral[k]),
            grossTransport, std::abs(s.chemistryChange[k])});
        if (scale == 0) check(s.balanceError[k] == 0, "inactive conservation equation has nonzero error");
        else error = std::max(error, std::abs(s.balanceError[k]) / scale);
    }
    check(error < 2e-11, "coupled finite-volume balance failed");
    for (double v : s.elementalBalanceError) check(std::abs(v) < 2e-8, "coupling violated elemental conservation");
    return error;
}
double homogeneous(DetailedGas& gas, const std::vector<double>& x, double temperature, double end,
    bool testRestart) {
    const auto m = scaledMesh(2, 2, 100);
    auto b = sealed(m);
    const auto q = gas.fromMoleAmounts(temperature, 101325, x);
    ReactingFlowStepper2D solver(gas, m, b);
    const auto initial = solver.initialState(std::vector<ReactingConservative2D>(m.cells.size(), reactingConservative2D(q)));
    ReactingStepControls2D controls; controls.maximumStep = end; controls.endTime = end;
    const auto step = solver.advance(initial, controls); accepted(step);
    check(step.accepted->time == end && step.accepted->steps == 1, "uniform source case did not reach requested time");
    const auto reference = gas.advanceConstantVolume(q, end); check(reference.accepted.has_value(), "reference source integration failed");
    const auto rp = gas.properties(*reference.accepted);
    double error = 0;
    for (const auto& u : step.accepted->cells) {
        const auto p = reactingPrimitive2D(gas, u);
        error = std::max(error, std::abs(p.properties.temperature / rp.temperature - 1));
        check(std::hypot(p.velocity.x, p.velocity.y) < 1e-8, "uniform reaction generated spurious flow");
        for (std::size_t k = 0; k < x.size(); ++k) error = std::max(error, std::abs(p.properties.massFractions[k] - rp.massFractions[k]));
    }
    // CVODES restart/local-error difference in a homogeneous exact splitting
    // case. This does not validate flame speed or a combustion mechanism.
    check(error < 2e-6, "split flow/source integration disagrees with homogeneous reactor");
    check(std::abs(step.afterIntegral[3] - step.beforeIntegral[3]) < 2e-11 * std::abs(step.beforeIntegral[3]), "reaction counted formation heat twice");
    (void)conservation(step);
    check(initial.time == 0 && initial.steps == 0 && initial.cells[0] == reactingConservative2D(q), "advance mutated accepted input");
    if (testRestart) {
        rejects([&] { std::ostringstream out; solver.writeCheckpoint(out, initial); }, "initial condition exported as accepted checkpoint");
        std::stringstream stream; solver.writeCheckpoint(stream, *step.accepted);
        const auto bytes = stream.str(); const auto restored = solver.readCheckpoint(stream);
        check(restored.cells == step.accepted->cells && restored.time == step.accepted->time && restored.steps == step.accepted->steps, "checkpoint round trip lost state/clock");
        // The full flame produced positive AR densities down to 4.3e-321.
        // libc++ stream extraction flags these as underflow despite returning
        // a representable value. Preserve them, including the smallest double.
        auto trace = *step.accepted;
        const auto argon = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), "AR");
        check(argon != gas.mechanism().species.end(), "checkpoint trace regression needs argon");
        const auto slot = 4 + static_cast<std::size_t>(argon - gas.mechanism().species.begin());
        const double tiny[] = {4.3529079602675429e-310, 4.3230744011109073e-321,
                               std::numeric_limits<double>::denorm_min(), 1e-308};
        for (std::size_t i = 0; i < trace.cells.size(); ++i) trace.cells[i][slot] = tiny[i % 4];
        std::stringstream traceStream; solver.writeCheckpoint(traceStream, trace);
        const auto traceBytes = traceStream.str();
        const auto traceRestored = solver.readCheckpoint(traceStream);
        check(traceRestored.cells == trace.cells && traceRestored.time == trace.time && traceRestored.steps == trace.steps,
              "checkpoint changed a representable subnormal species density");
        const auto stateStart = traceBytes.rfind("\nSTATE ");
        check(stateStart != std::string::npos, "checkpoint scalar fixture missing state");
        const auto firstScalar = traceBytes.find('\n', stateStart + 1) + 1;
        const auto scalarEnd = traceBytes.find(' ', firstScalar);
        for (const auto* invalid : {"nan", "inf", "1e9999", "1e-9999", "1.0tail"}) {
            auto corrupt = traceBytes; corrupt.replace(firstScalar, scalarEnd - firstScalar, invalid);
            std::stringstream reader(corrupt);
            rejects([&] { (void)solver.readCheckpoint(reader); }, "invalid or unrepresentable checkpoint scalar accepted");
        }
        controls.maximumStep = 1e-6; controls.endTime.reset();
        const auto continued = solver.advance(*step.accepted, controls), resumed = solver.advance(restored, controls); accepted(continued); accepted(resumed);
        check(continued.accepted->cells == resumed.accepted->cells && continued.accepted->time == resumed.accepted->time, "restart continuation not reproducible");
        auto damaged = initial; damaged.cells[0][4] = -1e-30;
        const auto rejected = solver.advance(damaged, controls);
        check(!rejected.accepted && !rejected.failure.empty() && damaged.time == 0, "negative species advanced clock");
        auto wrong = restored; wrong.binding += "changed";
        check(!solver.advance(wrong, controls).accepted, "changed mechanism/geometry binding accepted");
        std::stringstream truncated(bytes.substr(0, bytes.size() - 5));
        rejects([&] { (void)solver.readCheckpoint(truncated); }, "truncated checkpoint accepted");
        std::stringstream trailing(bytes + "extra"); rejects([&] { (void)solver.readCheckpoint(trailing); }, "trailing checkpoint data accepted");
        auto changedBoundary = b; changedBoundary.front().wallTemperature = 800;
        ReactingFlowStepper2D other(gas, m, changedBoundary);
        std::stringstream incompatible(bytes); rejects([&] { (void)other.readCheckpoint(incompatible); }, "different thermal boundary accepted on restart");
        ReactingFlowStepper2D frozen(gas, m, b, {false, true});
        std::stringstream physics(bytes); rejects([&] { (void)frozen.readCheckpoint(physics); }, "different reaction model accepted on restart");
        auto failedControls = controls; failedControls.maximumRetries = 0; failedControls.chemistry.maximumInternalSteps = 1;
        const auto failed = solver.advance(initial, failedControls);
        check(!failed.accepted && !failed.rejectedReasons.empty(), "chemistry budget failure was accepted");
        const auto recovered = solver.advance(restored, controls); accepted(recovered);
        check(recovered.accepted->cells == resumed.accepted->cells, "failed candidate contaminated later restart");
    }
    return error;
}
double acoustic(DetailedGas& gas, const std::vector<double>& x, int nx) {
    const auto m = fv_test::rectangle(nx, 2, 1, false);
    ReactingFlowStepper2D solver(gas, m, sealed(m), {false, false});
    const auto base = gas.fromMoleAmounts(1200, 101325, x); const auto bp = gas.properties(base);
    const double a = std::sqrt(bp.cp / bp.cv * bp.pressure / base.density), pi = std::acos(-1.), amplitude = 1e-4;
    std::vector<ReactingConservative2D> cells;
    for (const auto& c : m.cells) {
        const double dp = amplitude * bp.pressure * std::cos(pi * c.centre.x);
        const double p = bp.pressure + dp, rho = base.density + dp / (a * a);
        const double t = p * bp.meanMolecularWeight / (8314.46261815324 * rho);
        cells.push_back(reactingConservative2D(gas.fromMoleAmounts(t, p, x)));
    }
    auto state = solver.initialState(cells);
    const double end = .25 / a;
    ReactingStepControls2D controls; controls.maximumStep = end; controls.endTime = end;
    while (state.time < end) {
        const auto step = solver.advance(state, controls); accepted(step); (void)conservation(step);
        state = *step.accepted; check(state.steps < 2000, "acoustic run exceeded explicit budget");
    }
    double error = 0, norm = 0;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto p = reactingPrimitive2D(gas, state.cells[i]);
        const double expected = amplitude * bp.pressure * std::cos(pi * m.cells[i].centre.x) * std::cos(pi * a * end);
        error += m.cells[i].area * std::pow((p.properties.pressure - bp.pressure - expected) / (amplitude * bp.pressure), 2);
        norm += m.cells[i].area;
    }
    return std::sqrt(error / norm);
}
void variableViscosity() {
    const auto m = fv_test::rectangle(4, 3, 1, true);
    const auto velocity = [](Point2D p) { return Vector2D{.2 + .3 * p.x - .1 * p.y, -.1 + .4 * p.x + .6 * p.y}; };
    std::vector<Vector2D> u; for (const auto& c : m.cells) u.push_back(velocity(c.centre));
    std::vector<double> rho(m.cells.size(), 1.2), mu;
    std::vector<ViscousBoundary2D> bc;
    for (std::size_t f = 0; f < m.faces.size(); ++f) {
        mu.push_back(.2 + .1 * m.faces[f].centre.x);
        if (!m.faces[f].neighbour) bc.push_back({f, ViscousBoundaryKind2D::Velocity, velocity(m.faces[f].centre), {}});
    }
    ViscousStressOperator2D op(m, bc, 1);
    const auto r = op.evaluate(u, rho, mu);
    for (std::size_t f = 0; f < m.faces.size(); ++f) {
        const auto& s = m.faces[f].areaVector; const double div = .9;
        const Vector2D traction{mu[f] * ((.6 - 2 * div / 3) * s.x + .3 * s.y),
            mu[f] * (.3 * s.x + (1.2 - 2 * div / 3) * s.y)};
        const auto v = velocity(m.faces[f].centre);
        check(std::abs(r.faceFlux[f][0] + traction.x) < 2e-12 && std::abs(r.faceFlux[f][1] + traction.y) < 2e-12
            && std::abs(r.faceFlux[f][2] + dot(v, traction)) < 2e-12, "variable viscosity traction/work inconsistent");
    }
    std::vector<std::array<double,2>> row(m.cells.size());
    for (std::size_t column = 0; column < 2 * u.size(); ++column) {
        auto plus = u, minus = u;
        (column % 2 ? plus[column / 2].y : plus[column / 2].x) += 1;
        (column % 2 ? minus[column / 2].y : minus[column / 2].x) -= 1;
        const auto p = op.evaluate(plus, rho, mu), n = op.evaluate(minus, rho, mu);
        for (std::size_t i = 0; i < u.size(); ++i) for (std::size_t k = 0; k < 2; ++k) row[i][k] += .5 * std::abs(p.cellResidual[i][k] - n.cellResidual[i][k]);
    }
    for (std::size_t i = 0; i < u.size(); ++i)
        check(r.rate[i] >= .5 * std::max(row[i][0], row[i][1]) / (m.cells[i].area * rho[i]) * (1 - 2e-12), "variable viscosity rate underestimates Jacobian");
}
double diffusionRate(DetailedGas& gas) {
    const auto m = fv_test::rectangle(3, 2, 1, true);
    const auto ns = gas.mechanism().species.size();
    const auto q = gas.fromMassFractions(1100, 101325, std::vector<double>(ns, 1. / static_cast<double>(ns)));
    const auto p = gas.properties(q);
    std::vector<ReactingDiffusionBoundary2D> b;
    for (std::size_t f = 0; f < m.faces.size(); ++f) if (!m.faces[f].neighbour)
        b.push_back({f, ReactingDiffusionBoundaryKind2D::Reservoir, 0, q});
    ReactingDiffusionOperator2D op(m, b);
    const std::vector<GasState> states(m.cells.size(), q);
    const auto base = op.evaluate(gas, states, true);
    std::vector<std::vector<double>> norm(m.cells.size(), std::vector<double>(ns + 1));
    const double epsilon = 1e-6, energyScale = q.density * p.cv * p.temperature;
    auto speciesEnergy = p.speciesEnthalpies;
    for (std::size_t k = 0; k < ns; ++k)
        speciesEnergy[k] -= 8314.46261815324 * p.temperature / gas.mechanism().molecularWeights[k];
    const auto ref = static_cast<std::size_t>(std::max_element(p.massFractions.begin(), p.massFractions.end()) - p.massFractions.begin());
    for (std::size_t c = 0; c < states.size(); ++c) for (std::size_t column = 0; column <= ns; ++column) if (column != ref) {
        auto plus = states, minus = states;
        if (column == ns) { plus[c].internalEnergyDensity += energyScale * epsilon; minus[c].internalEnergyDensity -= energyScale * epsilon; }
        else {
            plus[c].speciesDensities[column] += q.density * epsilon; plus[c].speciesDensities[ref] -= q.density * epsilon;
            minus[c].speciesDensities[column] -= q.density * epsilon; minus[c].speciesDensities[ref] += q.density * epsilon;
            // Independent composition perturbations keep temperature fixed.
            const double de = q.density * epsilon * (speciesEnergy[column] - speciesEnergy[ref]);
            plus[c].internalEnergyDensity += de; minus[c].internalEnergyDensity -= de;
        }
        const auto a = op.evaluate(gas, plus), z = op.evaluate(gas, minus);
        for (std::size_t i = 0; i < states.size(); ++i) {
            for (std::size_t k = 0; k < ns; ++k) norm[i][k] += std::abs(a.speciesResidual[i][k] - z.speciesResidual[i][k]) / (2 * epsilon * m.cells[i].area * q.density);
            long double thermal = a.energyResidual[i] - z.energyResidual[i];
            for (std::size_t k = 0; k < ns; ++k)
                thermal -= speciesEnergy[k] * (a.speciesResidual[i][k] - z.speciesResidual[i][k]);
            norm[i][ns] += static_cast<double>(std::abs(thermal)) / (2 * epsilon * m.cells[i].area * energyScale);
        }
    }
    double maximumRatio = 0;
    for (std::size_t i = 0; i < states.size(); ++i) {
        const double measured = *std::max_element(norm[i].begin(), norm[i].end());
        maximumRatio = std::max(maximumRatio, measured / base.rate[i]);
        // At uniform fields coefficient derivatives multiply zero gradients;
        // this finite difference is the full frozen thermochemical Jacobian.
        check(measured <= base.rate[i] * (1 + 5e-5), "diffusion rate underestimates measured Jacobian");
    }
    return maximumRatio;
}
void mechanismBinding(const char* file) {
    DetailedGas gas(file);
    const auto m = fv_test::rectangle(1, 1, 1);
    ReactingFlowStepper2D solver(gas, m, sealed(m));
    const auto x = composition(gas, {{"H2", 2}, {"O2", 1}, {"N2", 3.76}});
    const auto initial = solver.initialState({reactingConservative2D(gas.fromMoleAmounts(1100, 101325, x))});
    auto definition = gas.mechanism().resolvedDefinition;
    const auto a = definition.find("A: "); check(a != std::string::npos, "serialized mechanism missing rates");
    std::size_t used = 0; const double old = std::stod(definition.substr(a + 3), &used);
    std::ostringstream value; value << std::setprecision(17) << old * 1.01;
    definition.replace(a + 3, used, value.str());
    const auto path = std::filesystem::temp_directory_path() / ("cartmesh2d-mechanism-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".yaml");
    { std::ofstream out(path); out << definition; check(static_cast<bool>(out), "cannot write mechanism regression fixture"); }
    DetailedGas changed(path.string());
    check(changed.mechanism().species == gas.mechanism().species && changed.mechanism().reactions == gas.mechanism().reactions,
          "mechanism-change fixture changed species or reaction count");
    check(changed.mechanism().resolvedDefinition != gas.mechanism().resolvedDefinition, "changed kinetic rate missing from binding");
    gas = std::move(changed);
    check(!solver.advance(initial).accepted, "replaced same-size chemistry context accepted");
    ReactingFlowStepper2D other(gas, m, sealed(m));
    check(!other.advance(initial).accepted, "changed reaction rates accepted with old state binding");
    std::filesystem::remove(path);
}
}
int main(int argc, char** argv) {
    try {
        check(argc == 3, "expected hydrogen and methane mechanism paths");
        DetailedGas gas(argv[1]);
        const auto x = composition(gas, {{"H2", 2}, {"O2", 1}, {"N2", 3.76}});
        const auto left = reactingConservative2D(gas.fromMoleAmounts(1100, 151325, x), {130, -20});
        const auto right = reactingConservative2D(gas.fromMoleAmounts(700, 101325, x), {-40, 10});
        const auto f = reactingFaceFlux2D(gas, left, right, {.2, .1});
        const auto reverse = reactingFaceFlux2D(gas, right, left, {-.2, -.1});
        for (std::size_t k = 0; k < left.size(); ++k)
            check(std::abs(f.integratedFlux[k] + reverse.integratedFlux[k]) < 2e-12 * (1 + std::abs(f.integratedFlux[k])), "HLLC face reversal failed");
        const auto mixture = gas.properties(gas.fromMoleAmounts(1100, 101325, x)).massFractions;
        for (double delta : {-1e-12, -1e-15, 1e-15, 1e-14, 1e-12}) {
            const auto a = reactingConservative2D(gas.fromMoleAmounts(1100, 101325 * (1 + delta), x));
            const auto b = reactingConservative2D(gas.fromMoleAmounts(1100, 101325, x));
            const auto tiny = reactingFaceFlux2D(gas, a, b, {1, 0});
            check(!tiny.hlleFallback, "near-zero-flow regression unexpectedly uses HLLE");
            for (std::size_t k = 0; k < mixture.size(); ++k) {
                const double expected = tiny.integratedFlux[0] * mixture[k];
                check(std::abs(tiny.integratedFlux[k + 4] - expected) <= 64 * std::numeric_limits<double>::epsilon() * std::abs(expected),
                      "HLLC species and density use inconsistent near-zero mass flux");
            }
        }
        const auto water = reactingConservative2D(gas.fromMoleAmounts(800, 101325, composition(gas, {{"H2O", 1}})), {20, 30});
        check(water[3] < 0 && std::abs(reactingPrimitive2D(gas, water).properties.temperature - 800) < 2e-9, "formation-inclusive total energy rejected");
        std::cerr << "hydrogen coupled source/restart\n";
        const double hydrogenError = homogeneous(gas, x, 1100, .001, true);
        std::cerr << "acoustics and variable viscosity\n";
        const double coarse = acoustic(gas, x, 12), fine = acoustic(gas, x, 24);
        check(fine < .06 && fine < .7 * coarse, "acoustic wave does not improve under refinement");
        variableViscosity();
        const double rateRatio = diffusionRate(gas);
        mechanismBinding(argv[1]);
        std::cerr << "methane coupled source\n";
        DetailedGas methane(argv[2]);
        const double methaneError = homogeneous(methane, composition(methane, {{"CH4", 1}, {"O2", 2}, {"N2", 7.52}}), 1400, .002, false);
        DetailedGas same(argv[1]);
        check(gas.mechanism().resolvedDefinition == same.mechanism().resolvedDefinition, "mechanism binding depends on mutable state/time");
        std::cout << std::setprecision(17) << "{\"passed\":true,\"hydrogenHomogeneousError\":" << hydrogenError
                  << ",\"methaneHomogeneousError\":" << methaneError << ",\"acousticCoarseL2\":" << coarse
                  << ",\"acousticFineL2\":" << fine << ",\"diffusionJacobianToBound\":" << rateRatio << "}\n";
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
