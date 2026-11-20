#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <iostream>

int main(int argc, char** argv) {
    using namespace cartmesh2d::structured;
    if (argc != 3) {
        std::cerr << "Usage: cartmesh2d_heat_cli INPUT.cfg OUTPUT_DIRECTORY\n";
        return 1;
    }
    try {
        const auto parameters = readParameters(argv[1]);
        const auto grid = makeGrid(parameters);
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
