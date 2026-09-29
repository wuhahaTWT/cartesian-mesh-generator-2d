#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d::chemistry;
namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Action> void rejects(Action action, const char* message) {
    bool failed = false;
    try { action(); } catch (const std::exception&) { failed = true; }
    check(failed, message);
}
std::vector<double> mixture(const DetailedGas& gas,
                           std::initializer_list<std::pair<std::string, double>> values) {
    std::vector<double> x(gas.mechanism().species.size());
    for (const auto& [name, value] : values) {
        const auto it = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), name);
        check(it != gas.mechanism().species.end(), "test species missing from mechanism");
        x[static_cast<std::size_t>(it - gas.mechanism().species.begin())] = value;
    }
    return x;
}
void ratesConserveElements(const DetailedGas& gas, const GasProperties& p) {
    const auto& m = gas.mechanism();
    for (std::size_t e = 0; e < m.elements.size(); ++e) {
        long double sum = 0, magnitude = 0;
        for (std::size_t k = 0; k < m.species.size(); ++k) {
            const long double rate = p.massProductionRates[k] * m.atomCounts[k * m.elements.size() + e]
                * m.atomicWeights[e] / m.molecularWeights[k];
            sum += rate;
            magnitude += p.massRateActivities[k] * m.atomCounts[k * m.elements.size() + e]
                * m.atomicWeights[e] / m.molecularWeights[k];
        }
        // Floating-point evaluation of an exact mechanism invariant, scaled by
        // gross creation/destruction. The near-cancelled net source is not a
        // rounding-error scale (the methane nitrogen regression catches this).
        // This is not a flame accuracy threshold.
        check(magnitude == 0 || std::abs(sum) / magnitude < 1024 * std::numeric_limits<double>::epsilon(),
              "reaction source violates elemental balance");
    }
}
}

int main(int argc, char** argv) {
    try {
        check(argc == 3, "expected explicit hydrogen and methane mechanism files");
        DetailedGas hydrogen(argv[1]);
        check(hydrogen.mechanism().species.size() == 10 && hydrogen.mechanism().reactions == 29,
              "hydrogen detailed mechanism was changed/reduced");
        const auto x = mixture(hydrogen, {{"H2", 2}, {"O2", 1}, {"N2", 3.76}});
        for (double temperature : {300., 800., 999., 1001., 2500.}) {
            for (double pressure : {101325., 10 * 101325.}) {
                const auto state = hydrogen.fromMoleAmounts(temperature, pressure, x);
                const auto p = hydrogen.properties(state);
                // Energy inversion should preserve the requested state to an
                // arithmetic allowance, separately from mechanism accuracy.
                check(std::abs(p.temperature - temperature) < 2e-9, "energy inversion changed temperature");
                check(std::abs(p.pressure / pressure - 1) < 2e-11, "energy inversion changed pressure");
                check(std::abs(p.cp - p.cv - 8314.46261815324 / p.meanMolecularWeight) < 1e-8,
                      "ideal mixture cp-cv is inconsistent");
                ratesConserveElements(hydrogen, p);
            }
        }
        const auto initial = hydrogen.fromMoleAmounts(1100, 101325, x);
        const auto saved = initial;
        const auto transport = hydrogen.transport(initial);
        check(transport.multicomponentDiffusion.size() == 100 && transport.thermalDiffusion.size() == 10,
              "multicomponent/Soret coefficients missing");
        check(transport.viscosity > 0 && transport.thermalConductivity > 0, "invalid transport properties");
        const auto step = hydrogen.advanceConstantVolume(initial, .001);
        if (!step.accepted) throw std::runtime_error(step.failure);
        const auto burned = hydrogen.properties(*step.accepted);
        check(burned.temperature > 2500 && burned.temperature < 3300, "hydrogen did not ignite");
        check(step.reachedTime == .001 && step.internalSteps > 0, "wrong chemistry clock or no integration");
        check(step.maximumElementDrift < 1e-10 && step.relativeEnergyDrift < 1e-10,
              "closed reactor conservation failed");
        ratesConserveElements(hydrogen, burned);
        check(initial.density == saved.density && initial.internalEnergyDensity == saved.internalEnergyDensity
              && initial.speciesDensities == saved.speciesDensities, "source step mutated accepted input");

        ChemistryControls insufficient;
        insufficient.maximumInternalSteps = 1;
        const auto failed = hydrogen.advanceConstantVolume(initial, .001, insufficient);
        check(!failed.accepted && !failed.failure.empty(), "step limit was accepted as successful chemistry");
        const auto again = hydrogen.advanceConstantVolume(initial, .001);
        check(again.accepted.has_value(), "failed source step poisoned next run");
        check(std::abs(hydrogen.properties(*again.accepted).temperature - burned.temperature) < 1e-7,
              "retry from accepted input is not reproducible");
        const auto badTime = hydrogen.advanceConstantVolume(initial, -1);
        check(!badTime.accepted && !badTime.failure.empty(), "negative chemistry time accepted");

        auto negative = initial;
        negative.speciesDensities[0] = -1e-30;
        rejects([&] { (void)hydrogen.properties(negative); }, "negative species was clipped/accepted");
        auto missingMass = initial;
        missingMass.speciesDensities[0] *= .5;
        rejects([&] { (void)hydrogen.properties(missingMass); }, "composition silently renormalized");
        auto impossible = initial;
        impossible.internalEnergyDensity = 1e20;
        rejects([&] { (void)hydrogen.properties(impossible); }, "thermodynamic extrapolation accepted");
        rejects([&] { (void)hydrogen.fromMoleAmounts(100, 101325, x); }, "out-of-range input temperature accepted");
        rejects([&] { (void)hydrogen.fromMoleAmounts(1000, -1, x); }, "negative pressure accepted");
        const auto water = hydrogen.fromMoleAmounts(800, 101325, mixture(hydrogen, {{"H2O", 1}}));
        check(water.internalEnergyDensity < 0, "species formation energy was discarded");
        check(std::abs(hydrogen.properties(water).temperature - 800) < 2e-9,
              "negative formation-inclusive internal energy was rejected");

        DetailedGas methane(argv[2]);
        check(methane.mechanism().species.size() == 53 && methane.mechanism().reactions == 325,
              "methane mechanism was reduced");
        const auto ch4 = methane.fromMoleAmounts(1400, 101325,
            mixture(methane, {{"CH4", 1}, {"O2", 2}, {"N2", 7.52}}));
        const auto methaneStep = methane.advanceConstantVolume(ch4, .002);
        if (!methaneStep.accepted) throw std::runtime_error(methaneStep.failure);
        ratesConserveElements(methane, methane.properties(*methaneStep.accepted));
        std::cout << "Detailed-gas invariants, full mechanisms, multicomponent transport, "
                     "stiff integration and rejection/retry passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
