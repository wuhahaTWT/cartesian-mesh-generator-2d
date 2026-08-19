#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace cartmesh2d;

int main() {
    constexpr double pi = 3.14159265358979323846;
    std::vector<Point2D> points;
    points.reserve(64);
    for (int i = 0; i < 64; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 64.0;
        points.push_back({std::cos(angle), std::sin(angle)});
    }

    BoundaryLoop circle(points);
    if (!circle.normalizeCounterClockwise()) return EXIT_FAILURE;

    const Domain2D domain{{{-2.0, -2.0}, {2.0, 2.0}}};
    Quadtree2D tree(domain, 6, circle);
    QuadtreeRefinementPolicy2D policy;
    policy.boundaryLevel = 5;
    policy.distanceBands = {{0.20, 4}};
    tree.refine(circle, policy);
    const auto balance = tree.enforceTwoToOneBalance(circle);
    if (balance.violationsAfter != 0) return EXIT_FAILURE;

    double totalFluidArea = 0.0;
    std::size_t unsupported = 0;
    std::size_t invalidPolygons = 0;
    std::size_t badAreaFractions = 0;

    for (const auto& leaf : tree.leaves()) {
        const auto cut = buildCutCell(leaf, circle);
        if (cut.kind == CutCellKind::Unsupported) {
            ++unsupported;
            continue;
        }
        if (cut.kind == CutCellKind::Empty) continue;

        totalFluidArea += cut.area;
        if (!(cut.areaFraction > 0.0 && cut.areaFraction <= 1.0 + 1.0e-12)) {
            ++badAreaFractions;
        }
        BoundaryLoop polygonLoop(cut.fluidPolygon.vertices);
        if (!polygonLoop.diagnose().valid()) ++invalidPolygons;
    }

    const double expectedArea = circle.polygon().area();
    const double error = std::abs(totalFluidArea - expectedArea);

    if (unsupported != 0 || invalidPolygons != 0 || badAreaFractions != 0 ||
        error > 1.0e-9) {
        std::cerr << "unsupported=" << unsupported
                  << " invalidPolygons=" << invalidPolygons
                  << " badAreaFractions=" << badAreaFractions
                  << " areaError=" << error << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "cartmesh2d 2D-3 quadtree cut-cell conservation audit passed\n";
    return EXIT_SUCCESS;
}
