#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {
int failures = 0;
void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}
}

int main() {
    const AABB2D unit{{0.0, 0.0}, {1.0, 1.0}};

    // The closed loop intersects the unit cell in two disconnected interior
    // fluid islands of equal area. The legacy single-cell API must reject the
    // ambiguity, while the solver API must preserve both components.
    BoundaryLoop disconnected({
        {0.82, -0.30}, {0.78, -0.30}, {0.22, -0.30}, {0.18, -0.30},
        {0.18, -0.20}, {0.18, 0.20}, {0.10, 0.20}, {0.10, 0.80},
        {0.30, 0.80}, {0.30, 0.20}, {0.22, 0.20}, {0.22, -0.20},
        {0.78, -0.20}, {0.78, 0.20}, {0.70, 0.20}, {0.70, 0.80},
        {0.90, 0.80}, {0.90, 0.20}, {0.82, 0.20}, {0.82, -0.20}});

    const auto single = buildCutCell(unit, CellClass::Intersected, disconnected,
                                     FluidRegion2D::Interior);
    check(!single.valid() && single.kind == CutCellKind::Unsupported,
          "single-component API rejects disconnected fluid islands");

    auto components = buildCutCells(unit, CellClass::Intersected, disconnected,
                                    FluidRegion2D::Interior);
    check(components.size() == 2,
          "solver Cut-cell API preserves both disconnected fluid components");
    std::vector<double> areas;
    double totalArea = 0.0;
    for (const auto& component : components) {
        check(component.valid() && component.kind == CutCellKind::Cut,
              "every disconnected component is a valid Cut-cell");
        check(component.centroid.has_value(),
              "every disconnected component has a centroid");
        check(component.fluidPolygon.signedArea() > 0.0,
              "every disconnected component has CCW fluid orientation");
        areas.push_back(component.area);
        totalArea += component.area;
    }
    std::sort(areas.begin(), areas.end());
    if (areas.size() == 2) {
        check(std::abs(areas[0] - 0.128) <= 1.0e-12 &&
              std::abs(areas[1] - 0.128) <= 1.0e-12,
              "both interior islands preserve their analytic area");
    }
    check(std::abs(totalArea - 0.256) <= 1.0e-12,
          "multi-component solver emission conserves total fluid area");

    // On the exterior side the same solid geometry leaves one connected fluid
    // component, so the default product API remains a single valid Cut-cell.
    const auto exterior = buildCutCells(unit, CellClass::Intersected, disconnected);
    check(exterior.size() == 1 && exterior.front().valid(),
          "default exterior fluid stays connected for this fixture");
    if (exterior.size() == 1) {
        check(std::abs(exterior.front().area - 0.744) <= 1.0e-12,
              "exterior complement area is exact");
    }

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d multi-component Cut-cell tests passed\n";
    return EXIT_SUCCESS;
}
