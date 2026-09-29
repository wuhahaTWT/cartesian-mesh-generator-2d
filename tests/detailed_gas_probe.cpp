#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d::chemistry;
namespace {
template<class Values> void numbers(const Values& values) {
    std::cout << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << value;
    }
    std::cout << ']';
}
void output(DetailedGas& gas, const GasState& state) {
    const auto p = gas.properties(state);
    const auto t = gas.transport(state);
    std::cout << "{\"temperature\":" << p.temperature << ",\"pressure\":" << p.pressure
              << ",\"density\":" << state.density << ",\"internalEnergyDensity\":" << state.internalEnergyDensity
              << ",\"cp\":" << p.cp << ",\"cv\":" << p.cv << ",\"molecularWeight\":" << p.meanMolecularWeight
              << ",\"viscosity\":" << t.viscosity << ",\"conductivity\":" << t.thermalConductivity
              << ",\"enthalpyReleaseRate\":" << p.enthalpyReleaseRate << ",\"Y\":";
    numbers(p.massFractions);
    std::cout << ",\"speciesEnthalpies\":"; numbers(p.speciesEnthalpies);
    std::cout << ",\"massProductionRates\":"; numbers(p.massProductionRates);
    std::cout << ",\"massRateActivities\":"; numbers(p.massRateActivities);
    std::cout << ",\"multicomponentDiffusion\":"; numbers(t.multicomponentDiffusion);
    std::cout << ",\"thermalDiffusion\":"; numbers(t.thermalDiffusion);
    std::cout << '}';
}
}
int main(int argc, char** argv) {
    try {
        if (argc < 6) throw std::invalid_argument("mechanism T[K] p[Pa] dt[s] SPECIES=MOLE_AMOUNT ...");
        DetailedGas gas(argv[1]);
        const double temperature = std::stod(argv[2]), pressure = std::stod(argv[3]), dt = std::stod(argv[4]);
        std::vector<double> x(gas.mechanism().species.size());
        for (int i = 5; i < argc; ++i) {
            const std::string arg(argv[i]);
            const auto separator = arg.find('=');
            if (separator == std::string::npos) throw std::invalid_argument("composition must be NAME=AMOUNT");
            const auto name = arg.substr(0, separator);
            const auto it = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), name);
            if (it == gas.mechanism().species.end()) throw std::invalid_argument("unknown species: " + name);
            x[static_cast<std::size_t>(it - gas.mechanism().species.begin())] += std::stod(arg.substr(separator + 1));
        }
        const auto initial = gas.fromMoleAmounts(temperature, pressure, x);
        std::cout << std::setprecision(17) << "{\"initial\":";
        output(gas, initial);
        if (dt != 0) {
            const auto step = gas.advanceConstantVolume(initial, dt);
            if (!step.accepted) throw std::runtime_error(step.failure);
            std::cout << ",\"final\":"; output(gas, *step.accepted);
            std::cout << ",\"time\":" << step.reachedTime << ",\"internalSteps\":" << step.internalSteps
                      << ",\"elementDrift\":" << step.maximumElementDrift
                      << ",\"energyDrift\":" << step.relativeEnergyDrift;
        }
        std::cout << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
