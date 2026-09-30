#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include "cartmesh2d/chemistry/SpeciesMassClosure.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
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
void sameThermodynamics(const GasThermodynamics& actual, const GasThermodynamics& expected) {
    check(actual.temperature == expected.temperature && actual.pressure == expected.pressure
          && actual.cp == expected.cp && actual.cv == expected.cv
          && actual.meanMolecularWeight == expected.meanMolecularWeight
          && actual.massFractions == expected.massFractions
          && actual.elementalMassFractions == expected.elementalMassFractions
          && actual.speciesEnthalpies == expected.speciesEnthalpies,
          "query order or exact-state reuse changed thermodynamics");
}
void thermodynamicQueries(DetailedGas& gas, const GasState& input) {
    const auto expected = gas.properties(input);
    sameThermodynamics(gas.thermodynamics(input), expected);
    // A direct T/p/composition setter invalidates an earlier conservative key.
    const auto other = gas.fromMassFractions(1500, 2 * expected.pressure, expected.massFractions);
    sameThermodynamics(gas.thermodynamics(input), expected);
    const auto transport = gas.transport(input);
    (void)gas.thermodynamics(other);
    const auto revisited = gas.transport(input);
    check(transport.viscosity == revisited.viscosity
          && transport.thermalConductivity == revisited.thermalConductivity
          && transport.multicomponentDiffusion == revisited.multicomponentDiffusion
          && transport.binaryDiffusion == revisited.binaryDiffusion
          && transport.thermalDiffusion == revisited.thermalDiffusion,
          "thermal-only query contaminated full multicomponent transport");
    const auto rates = gas.properties(input);
    check(rates.massProductionRates == expected.massProductionRates
          && rates.massRateActivities == expected.massRateActivities
          && rates.enthalpyReleaseRate == expected.enthalpyReleaseRate,
          "thermal-only query changed full reaction source");
    // A finite out-of-range target visits the bracket endpoints before failing;
    // an overflowing e = (rho*e)/rho must fail before an infinite error allowance
    // can make the temperature solve appear converged. Both invalidate the key.
    for (double energy : {input.density * 1e12, std::numeric_limits<double>::max() / 4}) {
        auto invalidEnergy = input;
        invalidEnergy.internalEnergyDensity = energy;
        rejects([&] { (void)gas.thermodynamics(invalidEnergy); }, "invalid thermal energy accepted");
        sameThermodynamics(gas.thermodynamics(input), expected);
    }
    // CVODES uses the same phase. A failed source solve can leave it at a trial
    // state even though the caller still holds the original accepted state.
    ChemistryControls limited; limited.maximumInternalSteps = 1;
    const auto failed = gas.advanceConstantVolume(input, .01, limited);
    check(!failed.accepted && !failed.failure.empty(), "source-budget invalidation regression did not fail");
    sameThermodynamics(gas.thermodynamics(input), expected);
    const auto argon = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), "AR");
    check(argon != gas.mechanism().species.end(), "exact-state key regression needs argon");
    const auto index = static_cast<std::size_t>(argon - gas.mechanism().species.begin());
    auto trace = input;
    trace.speciesDensities[index] = 1e-250;
    const auto before = gas.thermodynamics(trace);
    trace.speciesDensities[index] *= 2;
    const auto after = gas.thermodynamics(trace);
    check(after.massFractions[index] == trace.speciesDensities[index] / trace.density
          && after.massFractions[index] == 2 * before.massFractions[index],
          "exact-state key ignored a positive trace species change");
    trace.speciesDensities[index] = -1e-250;
    rejects([&] { (void)gas.thermodynamics(trace); }, "thermal-only query repaired negative trace species");
    sameThermodynamics(gas.thermodynamics(input), expected);
}
}

int main(int argc, char** argv) {
    try {
        check(argc == 3, "expected explicit hydrogen and methane mechanism files");
        DetailedGas hydrogen(argv[1]);
        check(hydrogen.mechanism().species.size() == 10 && hydrogen.mechanism().reactions == 29,
              "hydrogen detailed mechanism was changed/reduced");
        const auto x = mixture(hydrogen, {{"H2", 2}, {"O2", 1}, {"N2", 3.76}});
        // The repeated-source failure originated in a cold upstream cell.
        // Require bounded composition error after many genuine CVODES
        // restarts, not merely a one-shot algebraic closure test.
        auto cold = hydrogen.fromMoleAmounts(400, 101325, x);
        const auto coldInitial = cold;
        double maximumClosure = 0;
        for (unsigned i = 0; i < 2048; ++i) {
            const auto source = hydrogen.advanceConstantVolume(cold, 2.2e-10);
            if (!source.accepted) throw std::runtime_error(source.failure);
            cold = *source.accepted;
            const auto y = hydrogen.properties(cold).massFractions;
            check(std::abs(std::accumulate(y.begin(), y.end(), 0.) - 1) <= 4 * static_cast<double>(y.size()) * std::numeric_limits<double>::epsilon(),
                  "composition closure accumulated during repeated chemistry");
            for (double change : source.massClosureChange) maximumClosure = std::max(maximumClosure, std::abs(change / cold.density));
        }
        check(maximumClosure < 64 * 11 * std::numeric_limits<double>::epsilon(), "source closure exceeded roundoff scale");
        check(std::abs(cold.density / coldInitial.density - 1) < 2e-11, "repeated constant-volume source changed mass");
        std::vector<double> traceMass{.2, 1e-310, .8 + 1e-15, 0};
        const auto traceBefore = traceMass;
        const auto closed = closeSpeciesMassRoundoff(1, traceMass);
        check(closed.species == 2 && traceMass[0] == traceBefore[0] && traceMass[1] == traceBefore[1]
              && traceMass[3] == 0, "mass closure changed independent or trace species");
        for (const auto& invalid : {std::vector<double>{.2, .8 + 1e-8}, std::vector<double>{-1e-310, 1.}}) {
            auto candidate = invalid;
            rejects([&] { (void)closeSpeciesMassRoundoff(1, candidate); }, "invalid composition was repaired");
            check(candidate == invalid, "failed mass closure modified its input");
        }
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
        thermodynamicQueries(hydrogen, initial);
        const auto saved = initial;
        const auto transport = hydrogen.transport(initial);
        check(transport.multicomponentDiffusion.size() == 100 && transport.thermalDiffusion.size() == 10,
              "multicomponent/Soret coefficients missing");
        check(transport.viscosity > 0 && transport.thermalConductivity > 0, "invalid transport properties");
        // A -> viscosity-only B -> A used to leave the thermal cache's
        // temperature key valid while binary self-diffusion was overwritten.
        const auto anotherTemperature = hydrogen.fromMoleAmounts(1200, 101325, x);
        (void)hydrogen.viscosity(anotherTemperature);
        const auto revisited = hydrogen.transport(initial);
        check(transport.thermalConductivity == revisited.thermalConductivity
            && transport.thermalDiffusion == revisited.thermalDiffusion
            && transport.multicomponentDiffusion == revisited.multicomponentDiffusion,
            "transport depends on intervening viscosity-only temperature query");
        // A trace concentration gradient must survive inverse-matrix
        // cancellation and the backend's internal mole-fraction floor.
        const auto trace = hydrogen.fromMoleAmounts(1000, 101325, mixture(hydrogen, {{"H2", 1e-24}, {"N2", 1}}));
        const auto tp = hydrogen.properties(trace); const auto td = hydrogen.transport(trace);
        const auto& mechanism = hydrogen.mechanism();
        const auto h2 = static_cast<std::size_t>(std::find(mechanism.species.begin(), mechanism.species.end(), "H2") - mechanism.species.begin());
        const auto n2 = static_cast<std::size_t>(std::find(mechanism.species.begin(), mechanism.species.end(), "N2") - mechanism.species.begin());
        GasDiffusionGradient gradient; gradient.moleFractions.resize(mechanism.species.size());
        gradient.moleFractions[h2] = 1e-20; gradient.moleFractions[n2] = -1e-20;
        const auto traceFlux = hydrogen.diffusiveFlux(trace, gradient);
        const double binary = -trace.density * mechanism.molecularWeights[h2] * mechanism.molecularWeights[n2]
            / (tp.meanMolecularWeight * tp.meanMolecularWeight) * td.binaryDiffusion[h2 + mechanism.species.size() * n2] * 1e-20;
        check(std::abs((traceFlux.species[h2] - binary) / binary) < 5e-12, "trace binary analytic diffusion limit failed");
        // IEEE gradual-underflow regression: keep positive species and solve
        // the physical flux even when its storage has an absolute half-ulp
        // error instead of the normal relative-epsilon error model.
        for (double dx : {1e-310, 1e-318}) {
            gradient.moleFractions[h2] = dx; gradient.moleFractions[n2] = -dx;
            const auto tinyFlux = hydrogen.diffusiveFlux(trace, gradient);
            const double expected = (binary / 1e-20) * dx;
            const double rounding = 64 * std::numeric_limits<double>::epsilon() * std::abs(expected)
                + 2 * std::numeric_limits<double>::denorm_min();
            check(std::abs(tinyFlux.species[h2] - expected) <= rounding, "subnormal binary flux lost its representable value");
        }
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
        thermodynamicQueries(methane, ch4);
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
