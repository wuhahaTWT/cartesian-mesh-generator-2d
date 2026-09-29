#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d::chemistry;
int main(int argc, char** argv) {
    try {
        const bool batch = argc == 3 && std::string(argv[2]) == "--batch";
        if (argc != 2 && !batch) throw std::invalid_argument("expected mechanism path [--batch]; stdin: T p dT/dn dlnp/dn then N mole amounts and N mole-fraction gradients; one complete record per line in batch mode");
        DetailedGas gas(argv[1]);
        do {
            std::string line;
            if (batch && !std::getline(std::cin, line)) break;
            std::istringstream lineInput(line);
            std::istream& input = batch ? static_cast<std::istream&>(lineInput) : std::cin;
            double temperature = 0, pressure = 0;
            GasDiffusionGradient d;
            input >> temperature >> pressure >> d.temperature >> d.logPressure;
            std::vector<double> amounts(gas.mechanism().species.size());
            d.moleFractions.resize(amounts.size());
            for (double& x : amounts) input >> x;
            for (double& x : d.moleFractions) input >> x;
            if (!input) throw std::invalid_argument("incomplete numeric input");
            std::string extra;
            if (input >> extra) throw std::invalid_argument("unexpected trailing input");
            const auto state = gas.fromMoleAmounts(temperature, pressure, amounts);
            const auto flux = gas.diffusiveFlux(state, d);
            std::cout << std::setprecision(17) << "{\"species\":[";
            for (std::size_t k = 0; k < flux.species.size(); ++k) {
                if (k) std::cout << ',';
                std::cout << flux.species[k];
            }
            std::cout << "],\"conduction\":" << flux.conduction << ",\"speciesEnthalpy\":" << flux.speciesEnthalpy
                      << ",\"energy\":" << flux.energy << ",\"rawMassResidual\":" << flux.rawMassResidual
                      << ",\"closureCorrection\":" << flux.closureCorrection
                      << ",\"zeroLimitCorrection\":" << flux.zeroLimitCorrection
                      << ",\"concentrationSolveResidual\":" << flux.concentrationSolveResidual
                      << ",\"matrixFormDifferenceL1\":" << flux.matrixFormDifferenceL1
                      << ",\"pureSpeciesLimit\":" << (flux.pureSpeciesLimit ? "true" : "false") << "}\n" << std::flush;
        } while (batch);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
