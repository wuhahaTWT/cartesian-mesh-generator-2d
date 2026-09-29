#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d::chemistry;
int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("expected mechanism path; stdin: T p dT/dn dlnp/dn then N mole amounts and N mole-fraction gradients");
        DetailedGas gas(argv[1]);
        double temperature = 0, pressure = 0;
        GasDiffusionGradient d;
        std::cin >> temperature >> pressure >> d.temperature >> d.logPressure;
        std::vector<double> amounts(gas.mechanism().species.size());
        d.moleFractions.resize(amounts.size());
        for (double& x : amounts) std::cin >> x;
        for (double& x : d.moleFractions) std::cin >> x;
        if (!std::cin) throw std::invalid_argument("incomplete numeric input");
        std::string extra;
        if (std::cin >> extra) throw std::invalid_argument("unexpected trailing input");
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
                  << ",\"pureSpeciesLimit\":" << (flux.pureSpeciesLimit ? "true" : "false") << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
