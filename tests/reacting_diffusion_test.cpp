#include "cartmesh2d/fv/ReactingDiffusion2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::chemistry;
namespace {
void check(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
template<class Action> void rejects(Action action, const char* message) {
    bool failed = false;
    try { action(); } catch (const std::exception&) { failed = true; }
    check(failed, message);
}
std::size_t species(const DetailedGas& gas, const std::string& name) {
    const auto& names = gas.mechanism().species;
    const auto it = std::find(names.begin(), names.end(), name);
    check(it != names.end(), "missing test species");
    return static_cast<std::size_t>(it - names.begin());
}
GasState field(DetailedGas& gas, Point2D point) {
    std::vector<double> x(gas.mechanism().species.size(), .0001);
    x[species(gas, "H2")] = .15 + .05 * point.x;
    x[species(gas, "O2")] = .2 + .02 * point.y;
    x[species(gas, "N2")] = .65 - .05 * point.x - .02 * point.y;
    return gas.fromMoleAmounts(800 + 100 * point.x + 40 * point.y,
        101325 * std::exp(.2 * point.x - .1 * point.y), x);
}
std::vector<ReactingDiffusionBoundary2D> walls(const FvMesh2D& mesh) {
    std::vector<ReactingDiffusionBoundary2D> b;
    for (std::size_t f = 0; f < mesh.faces.size(); ++f) if (!mesh.faces[f].neighbour)
        b.push_back({f, ReactingDiffusionBoundaryKind2D::AdiabaticWall, 0, {}});
    return b;
}
// Exact algebraic conservation is checked against roundoff, scaled by gross
// edge transport. This numerical assembly check is not physical flame accuracy.
double budget(const ReactingDiffusionResult2D& r) {
    long double maxError = 0, energySum = 0, energyScale = 0;
    for (double q : r.energyResidual) energySum += q;
    for (const auto& f : r.faces) energyScale += std::abs(f.energy);
    const long double energyError = std::abs(energySum - r.boundaryEnergyFlux);
    check(energyError <= 4096 * std::numeric_limits<double>::epsilon() * energyScale, "energy face assembly is not conservative");
    if (energyScale > 0) maxError = energyError / energyScale;
    for (std::size_t k = 0; k < r.boundarySpeciesFlux.size(); ++k) {
        long double sum = 0, scale = 0;
        for (const auto& row : r.speciesResidual) sum += row[k];
        for (const auto& f : r.faces) scale += std::abs(f.species[k]);
        const long double error = std::abs(sum - r.boundarySpeciesFlux[k]);
        check(error <= 4096 * std::numeric_limits<double>::epsilon() * scale, "species face assembly is not conservative");
        if (scale > 0) maxError = std::max(maxError, error / scale);
    }
    return static_cast<double>(maxError);
}
}

int main(int argc, char** argv) {
    try {
        check(argc == 2, "expected hydrogen mechanism");
        DetailedGas gas(argv[1]);
        const auto ns = gas.mechanism().species.size();
        const auto q = field(gas, {.2, .4});
        GasDiffusionGradient d{2e4, 50, std::vector<double>(ns)};
        d.moleFractions[species(gas, "H2")] = 50;
        d.moleFractions[species(gas, "N2")] = -50;
        const auto f = gas.diffusiveFlux(q, d);
        d.temperature *= -1; d.logPressure *= -1;
        for (double& x : d.moleFractions) x *= -1;
        const auto reverse = gas.diffusiveFlux(q, d);
        for (std::size_t k = 0; k < ns; ++k) check(f.species[k] == -reverse.species[k], "normal reversal changed diffusion law");
        check(f.energy == -reverse.energy, "energy flux is not antisymmetric");
        check(f.speciesEnthalpy != 0 && f.conduction != 0, "energy transport term missing");
        d.moleFractions[0] += 1;
        rejects([&] { (void)gas.diffusiveFlux(q, d); }, "inconsistent gradients were normalized");
        d.moleFractions.resize(ns - 1);
        rejects([&] { (void)gas.diffusiveFlux(q, d); }, "wrong gradient size accepted");
        d.moleFractions.assign(ns, 0); d.temperature = std::numeric_limits<double>::quiet_NaN();
        rejects([&] { (void)gas.diffusiveFlux(q, d); }, "nonfinite gradient accepted");
        d.temperature = 0; d.logPressure = 0;
        const auto zero = gas.diffusiveFlux(q, d);
        check(zero.energy == 0 && std::all_of(zero.species.begin(), zero.species.end(), [](double v) { return v == 0; }),
              "constant thermodynamic state diffuses");
        // Minimal failure: pure water, dT/dn=100000 K/m. Backend concentration
        // regularization previously produced false Soret flux / closure failure.
        std::vector<double> waterAmounts(ns);
        waterAmounts[species(gas, "H2O")] = 1;
        const auto water = gas.fromMoleAmounts(800, 101325, waterAmounts);
        d.temperature = 100000;
        const auto pure = gas.diffusiveFlux(water, d);
        check(pure.pureSpeciesLimit && pure.conduction < 0 && pure.speciesEnthalpy == 0,
              "pure gas diffusion limit lost Fourier conduction");
        for (double value : pure.species) check(value == 0, "pure gas separated into absent species");
        std::vector<double> binaryAmounts(ns);
        binaryAmounts[species(gas, "H2")] = 1; binaryAmounts[species(gas, "N2")] = 1;
        const auto binary = gas.fromMoleAmounts(900, 101325, binaryAmounts);
        d.moleFractions[species(gas, "H2")] = 20; d.moleFractions[species(gas, "N2")] = -20;
        const auto binaryFlux = gas.diffusiveFlux(binary, d);
        for (std::size_t k = 0; k < ns; ++k) if (binaryAmounts[k] == 0)
            check(binaryFlux.species[k] == 0, "zero concentration with zero gradient acquired spurious trace flux");
        // A zero concentration WITH a gradient is not blocked (e.g. a reservoir
        // face into which a previously absent constituent diffuses).
        d.temperature = 0; d.moleFractions.assign(ns, 0);
        d.moleFractions[species(gas, "H2")] = 1;
        d.moleFractions[species(gas, "H2O")] = -1;
        const auto entering = gas.diffusiveFlux(water, d);
        check(entering.species[species(gas, "H2")] < 0, "zero-species limit incorrectly blocks a nonzero gradient");

        const auto mesh = fv_test::rectangle(8, 6, 1, true);
        auto boundary = walls(mesh);
        const ReactingDiffusionOperator2D sealed(mesh, boundary);
        std::vector<GasState> cells;
        for (const auto& c : mesh.cells) cells.push_back(field(gas, c.centre));
        const auto original = cells;
        const auto closed = sealed.evaluate(gas, cells);
        const double closedBudget = budget(closed);
        check(closed.boundaryEnergyFlux == 0, "adiabatic wall leaked energy");
        check(std::all_of(closed.boundarySpeciesFlux.begin(), closed.boundarySpeciesFlux.end(), [](double v) { return v == 0; }),
              "impermeable wall leaked species");
        // Exercise the actual conservative residual with a small explicit test
        // substep. This is a harness, not a stability-certified time integrator.
        const double dt = 1e-4;
        long double beforeEnergy = 0, afterEnergy = 0;
        std::vector<long double> beforeSpecies(ns), afterSpecies(ns);
        for (std::size_t i = 0; i < cells.size(); ++i) {
            const double area = mesh.cells[i].area;
            cells[i].internalEnergyDensity -= dt * closed.energyResidual[i] / area;
            for (std::size_t k = 0; k < ns; ++k) {
                cells[i].speciesDensities[k] -= dt * closed.speciesResidual[i][k] / area;
                beforeSpecies[k] += original[i].speciesDensities[k] * area;
                afterSpecies[k] += cells[i].speciesDensities[k] * area;
            }
            (void)gas.properties(cells[i]); // positivity, closure and energy inversion
            beforeEnergy += original[i].internalEnergyDensity * area;
            afterEnergy += cells[i].internalEnergyDensity * area;
        }
        check(std::abs(afterEnergy - beforeEnergy) < 1e-12L * std::abs(beforeEnergy), "diffusion substep changed closed energy");
        for (std::size_t k = 0; k < ns; ++k)
            check(std::abs(afterSpecies[k] - beforeSpecies[k]) < 1e-12L * beforeSpecies[k], "diffusion substep changed global species");

        for (auto& b : boundary) {
            b.kind = ReactingDiffusionBoundaryKind2D::Reservoir;
            b.reservoir = field(gas, mesh.faces[b.face].centre);
        }
        const ReactingDiffusionOperator2D open(mesh, boundary);
        const auto transported = open.evaluate(gas, original);
        const double openBudget = budget(transported);
        // On this affine T/log(p)/X field, normal reconstruction must reproduce
        // the specified gradient, even on nonorthogonal/skew faces. Coefficients
        // deliberately use the documented convex thermodynamic face state.
        double affineError = 0;
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            const double length = std::hypot(face.areaVector.x, face.areaVector.y), nx = face.areaVector.x / length, ny = face.areaVector.y / length;
            const auto left = gas.properties(original[face.owner]);
            const auto right = face.neighbour ? gas.properties(original[*face.neighbour]) : gas.properties(field(gas, face.centre));
            const double w = face.neighbour ? face.neighbourWeight : 1;
            std::vector<double> y(ns);
            for (std::size_t k = 0; k < ns; ++k) y[k] = left.massFractions[k] * (1 - w) + right.massFractions[k] * w;
            const auto faceState = gas.fromMassFractions(left.temperature * (1 - w) + right.temperature * w,
                std::exp(std::log(left.pressure) * (1 - w) + std::log(right.pressure) * w), y);
            const double moleTotal = 1 + static_cast<double>(ns - 3) * .0001;
            GasDiffusionGradient exact{100 * nx + 40 * ny, .2 * nx - .1 * ny, std::vector<double>(ns)};
            exact.moleFractions[species(gas, "H2")] = .05 * nx / moleTotal;
            exact.moleFractions[species(gas, "O2")] = .02 * ny / moleTotal;
            exact.moleFractions[species(gas, "N2")] = -(.05 * nx + .02 * ny) / moleTotal;
            const auto expected = gas.diffusiveFlux(faceState, exact);
            double magnitude = 0, error = 0;
            for (std::size_t k = 0; k < ns; ++k) {
                magnitude += std::abs(expected.species[k] * length);
                error = std::max(error, std::abs(expected.species[k] * length - transported.faces[id].species[k]));
            }
            affineError = std::max(affineError, error / magnitude);
            // Interface arithmetic/reconstruction threshold; no claim about
            // truncation error on nonlinear fields or physical transport data.
            check(error / magnitude < 2e-10, "affine multicomponent normal gradient failed");
            check(std::abs(expected.energy * length - transported.faces[id].energy)
                    < 2e-10 * length * (std::abs(expected.conduction) + std::abs(expected.speciesEnthalpy)), "affine enthalpy flux failed");
        }
        for (auto& b : boundary) {
            b.kind = ReactingDiffusionBoundaryKind2D::IsothermalWall;
            b.reservoir.reset(); b.temperature = 1000;
        }
        const ReactingDiffusionOperator2D heated(mesh, boundary);
        const auto heating = heated.evaluate(gas, original);
        (void)budget(heating);
        check(heating.boundaryEnergyFlux < 0, "hot wall failed to heat cooler gas");
        for (const auto& b : boundary) for (double j : heating.faces[b.face].species) check(j == 0, "isothermal wall leaked species");
        auto invalid = original; invalid[0].speciesDensities[0] = -1e-20;
        rejects([&] { (void)sealed.evaluate(gas, invalid); }, "negative species accepted in spatial operator");
        auto incomplete = boundary; incomplete.pop_back();
        rejects([&] { ReactingDiffusionOperator2D bad(mesh, incomplete); }, "missing diffusion boundary accepted");
        auto duplicate = boundary; duplicate.push_back(duplicate.front());
        rejects([&] { ReactingDiffusionOperator2D bad(mesh, duplicate); }, "duplicate boundary accepted");
        const auto recovered = sealed.evaluate(gas, original);
        check(recovered.energyResidual == closed.energyResidual, "failed call changed later accepted residual");
        for (std::size_t i = 0; i < original.size(); ++i)
            check(original[i].speciesDensities == field(gas, mesh.cells[i].centre).speciesDensities,
                  "spatial operator mutated input");
        std::cout << std::setprecision(17) << "{\"passed\":true,\"cells\":" << mesh.cells.size()
                  << ",\"faces\":" << mesh.faces.size() << ",\"closedBudgetRelative\":" << closedBudget
                  << ",\"openBudgetRelative\":" << openBudget << ",\"affineFluxRelative\":" << affineError << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
