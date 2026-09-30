#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::chemistry;
using namespace cartmesh2d::fv;
namespace {
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
}
int main(int argc, char** argv) {
    try {
        check(argc == 2, "expected H2 mechanism");
        DetailedGas gas(argv[1]); const auto& names = gas.mechanism().species;
        const auto index = [&](const char* name) {
            const auto it = std::find(names.begin(), names.end(), name);
            check(it != names.end(), "test species missing"); return static_cast<std::size_t>(it - names.begin());
        };
        const auto h2 = index("H2"), o2 = index("O2"), n2 = index("N2"), ar = index("AR");
        auto mesh = fv_test::rectangle(8, 4, 1, true);
        std::vector<ReactingBoundary2D> boundaries;
        for (std::size_t i = 0; i < mesh.faces.size(); ++i) if (!mesh.faces[i].neighbour)
            boundaries.push_back({i, ReactingBoundaryKind2D::SlipWall, {}, {}, 0});
        // Isolate the reconstruction/advective operator, without reaction or
        // molecular transport obscuring a trace-dependent limiter jump.
        ReactingFlowStepper2D solver(gas, mesh, boundaries, {false, false});
        std::vector<ReactingConservative2D> cells;
        for (const auto& cell : mesh.cells) {
            std::vector<double> x(names.size()); x[h2] = 1.6 * (1 + .25 * cell.centre.x); x[o2] = 1; x[n2] = 3.76;
            // Keep cell and reconstructed face temperatures away from the
            // mechanism's 1000 K NASA polynomial switch. The supplied fits
            // have a small enthalpy mismatch there, which is a separate data
            // continuity issue and must not masquerade as a limiter defect.
            cells.push_back(reactingConservative2D(gas.fromMoleAmounts(1150 + 200 * cell.centre.x, 101325, x), {2, .2}));
        }
        const auto base = solver.initialState(cells);
        const auto original = solver.evaluateResidual(base, 2, false);
        const auto firstOrder = solver.evaluateResidual(base, 1, false);
        bool discontinuous = false;
        for (double amplitude : {1e-12, 1e-16}) {
            auto perturbed = base;
            for (std::size_t i = 0; i < cells.size(); ++i) {
                const double fraction = amplitude * (i % 3 == 0 ? 1. : .1);
                const double added = cells[i][0] * fraction;
                perturbed.cells[i][ar + 4] += added; perturbed.cells[i][n2 + 4] -= added;
            }
            const auto changed = solver.evaluateResidual(perturbed, 2, false);
            const auto firstChanged = solver.evaluateResidual(perturbed, 1, false);
            double dp = 0, dt = 0, firstDifference = 0, firstScale = 0;
            for (std::size_t i = 0; i < cells.size(); ++i) {
                const auto a = reactingPrimitive2D(gas, base.cells[i]), b = reactingPrimitive2D(gas, perturbed.cells[i]);
                dp = std::max(dp, std::abs(a.properties.pressure - b.properties.pressure));
                dt = std::max(dt, std::abs(a.properties.temperature - b.properties.temperature));
            }
            for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
                firstDifference += std::abs(firstChanged.faceFlux[f][0] - firstOrder.faceFlux[f][0]);
                firstScale += std::abs(firstOrder.faceFlux[f][0]);
            }
            std::cout << "cell changes: p=" << dp << " T=" << dt << " first-order mass response=" << firstDifference / firstScale
                << " HLLE=" << original.hlleFallbacks << '/' << changed.hlleFallbacks << '\n';
            for (std::size_t component : {std::size_t{0}, std::size_t{3}, h2 + 4}) {
                double difference = 0, scale = 0;
                for (std::size_t face = 0; face < mesh.faces.size(); ++face) {
                    difference += std::abs(changed.faceFlux[face][component] - original.faceFlux[face][component]);
                    scale += std::abs(original.faceFlux[face][component]);
                }
                const double relative = difference / scale;
                std::cout << "trace amplitude=" << amplitude << " component=" << component << " relative flux response=" << relative << '\n';
                // The perturbation has the stated mass-fraction amplitude.
                // A broad 1000x Lipschitz/roundoff allowance detects a finite
                // bulk-flux jump as a physically irrelevant trace tends to zero;
                // it is not a flame-accuracy or universal residual threshold.
                discontinuous = discontinuous || relative > 1000 * (amplitude + std::numeric_limits<double>::epsilon());
            }
        }
        check(!discontinuous, "vanishing trace caused a finite bulk-flux jump");
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
