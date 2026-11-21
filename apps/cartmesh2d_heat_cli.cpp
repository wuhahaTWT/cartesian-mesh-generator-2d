#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc, char** argv) {
    using namespace cartmesh2d::structured;
    const bool surfaceOnly = argc == 5 && std::string_view(argv[1]) == "--surface";
    if (argc != 3 && !surfaceOnly) {
        std::cerr << "Usage: cartmesh2d_heat_cli INPUT.cfg OUTPUT_DIRECTORY\n"
                     "       cartmesh2d_heat_cli --surface INPUT.cfg FIELD.csv SURFACE.csv\n";
        return 1;
    }
    try {
        const auto parameters = readParameters(argv[surfaceOnly ? 2 : 1]);
        const auto grid = makeGrid(parameters);
        if (surfaceOnly) {
            if (parameters.problem != Case::Square) throw std::invalid_argument("Surface profiles require a square cylinder case");
            writeSquareSurface(argv[4], grid, readFieldCsv(argv[3], grid, false));
            return 0;
        }
        const auto result = solve(grid, parameters, [](bool thermal, const Residual& r) {
            if (r.iteration == 1 || r.iteration % 100 == 0)
                std::cout << (thermal ? "thermal " : "flow ") << r.iteration << ' ' << r.continuity
                    << ' ' << r.u << ' ' << r.v << ' ' << r.temperature << std::endl;
        });
        writeResult(argv[2], grid, parameters, result);
        std::cout << result.stop << ": " << result.flowHistory.size() << " flow iterations ("
            << result.flowSeconds << " s), " << result.thermalHistory.size() << " thermal iterations ("
            << result.thermalSeconds << " s)\n";
        return result.flowConverged && (parameters.problem != Case::Square || result.thermalConverged) ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
