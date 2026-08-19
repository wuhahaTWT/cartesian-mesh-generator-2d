#include "cartmesh2d/grid/CartesianGrid2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <stdexcept>
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

void near(double a, double b, double eps, const std::string& message) {
    check(std::abs(a - b) <= eps, message);
}
} // namespace

int main() {
    const Domain2D domain{{{0.0, 0.0}, {4.0, 3.0}}};
    UniformCartesianGrid2D grid(domain, 4, 3);
    check(grid.cells().size() == 12, "4x3 grid has 12 cells");
    near(grid.dx(), 1.0, 1e-12, "dx=1");
    near(grid.dy(), 1.0, 1e-12, "dy=1");
    check(grid.cell(0, 0).id == 0 && grid.cell(3, 0).id == 3 && grid.cell(0, 1).id == 4,
          "row-major deterministic IDs");
    near(grid.cell(2, 1).center().x, 2.5, 1e-12, "cell center x");
    near(grid.cell(2, 1).center().y, 1.5, 1e-12, "cell center y");

    double coveredArea = 0.0;
    for (const auto& cell : grid.cells()) {
        coveredArea += (cell.bounds.max.x - cell.bounds.min.x) *
                       (cell.bounds.max.y - cell.bounds.min.y);
    }
    near(coveredArea, 12.0, 1e-12, "cells exactly cover domain area");
    near(grid.cell(0, 0).bounds.max.x, grid.cell(1, 0).bounds.min.x, 1e-12,
         "horizontal neighbors meet without overlap/gap");
    near(grid.cell(0, 0).bounds.max.y, grid.cell(0, 1).bounds.min.y, 1e-12,
         "vertical neighbors meet without overlap/gap");

    auto spacingGrid = UniformCartesianGrid2D::fromTargetSpacing(domain, 1.5, 1.1);
    check(spacingGrid.nx() == 3 && spacingGrid.ny() == 3,
          "spacing-based dimensions use ceil");
    check(spacingGrid.dx() <= 1.5 + 1e-12 && spacingGrid.dy() <= 1.1 + 1e-12,
          "actual spacing does not exceed target");

    const AABB2D unit{{0.0, 0.0}, {1.0, 1.0}};
    check(segmentIntersectsClosedAABB({{-1.0, 0.5}, {2.0, 0.5}}, unit),
          "crossing segment hits box");
    check(segmentIntersectsClosedAABB({{-1.0, 1.0}, {2.0, 1.0}}, unit),
          "segment on box edge counts as hit");
    check(!segmentIntersectsClosedAABB({{-1.0, 2.0}, {2.0, 2.0}}, unit),
          "separated segment misses box");
    check(segmentIntersectsClosedAABB({{-1.0, -1.0}, {0.0, 0.0}}, unit),
          "corner tangent counts as hit");
    check(segmentIntersectsClosedAABB({{0.2, 0.2}, {0.8, 0.8}}, unit),
          "segment fully inside box hits box");

    BoundaryLoop rect({{1.0, 1.0}, {3.0, 1.0}, {3.0, 3.0}, {1.0, 3.0}});
    UniformCartesianGrid2D classificationGrid({{{0.0, 0.0}, {4.0, 4.0}}}, 8, 8);
    const auto summary = classifyGrid(classificationGrid, rect);
    check(summary.inside == 4, "aligned rectangle has 4 strictly interior cells");
    check(summary.intersected == 32,
          "closed-AABB rule marks both sides of grid-line boundary");
    check(summary.outside == 28, "remaining cells are outside");

    UniformCartesianGrid2D classificationGridRepeat({{{0.0, 0.0}, {4.0, 4.0}}}, 8, 8);
    const auto repeatSummary = classifyGrid(classificationGridRepeat, rect);
    check(repeatSummary.inside == summary.inside && repeatSummary.outside == summary.outside &&
              repeatSummary.intersected == summary.intersected,
          "repeat classification summary stable");
    bool repeatStable = true;
    for (std::size_t k = 0; k < classificationGrid.cells().size(); ++k) {
        if (classificationGrid.cells()[k].id != classificationGridRepeat.cells()[k].id ||
            classificationGrid.cells()[k].classification !=
                classificationGridRepeat.cells()[k].classification) {
            repeatStable = false;
            break;
        }
    }
    check(repeatStable, "repeat IDs and classifications stable");

    check(classificationGrid.cell(1, 1).classification == CellClass::Intersected,
          "outside-adjacent grid-line cell is intersected");
    check(classificationGrid.cell(3, 3).classification == CellClass::Inside,
          "strictly interior cell remains inside");
    check(classificationGrid.cell(0, 0).classification == CellClass::Outside,
          "far cell remains outside");

    BoundaryLoop diamond({{1.0, 1.0}, {2.0, 0.5}, {3.0, 1.0}, {2.0, 2.0}});
    UniformCartesianGrid2D tangentGrid({{{0.0, 0.0}, {4.0, 4.0}}}, 4, 4);
    const auto tangentSummary = classifyGrid(tangentGrid, diamond);
    check(tangentSummary.intersected > 0,
          "diagonal/tangent boundary produces intersected cells");
    check(tangentGrid.cell(0, 0).classification == CellClass::Intersected,
          "boundary touching corner (1,1) marks cell as intersected");

    bool threw = false;
    try {
        UniformCartesianGrid2D bad(domain, 0, 2);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "zero grid dimension rejected");

    std::vector<Point2D> circleVertices;
    for (int k = 0; k < 64; ++k) {
        const double theta = 2.0 * std::numbers::pi * static_cast<double>(k) / 64.0;
        circleVertices.push_back({std::cos(theta), std::sin(theta)});
    }
    BoundaryLoop circle(std::move(circleVertices));
    UniformCartesianGrid2D circleGrid({{{-2.0, -2.0}, {2.0, 2.0}}}, 32, 32);
    const auto circleSummary = classifyGrid(circleGrid, circle);
    check(circleSummary.inside == 164 && circleSummary.intersected == 68 &&
              circleSummary.outside == 792,
          "64-segment circle classification counts are reproducible");

    threw = false;
    try {
        BoundaryLoop bowTie({{0.5, 0.5}, {3.5, 3.5}, {0.5, 3.5}, {3.5, 0.5}});
        UniformCartesianGrid2D invalidBoundaryGrid({{{0.0, 0.0}, {4.0, 4.0}}}, 4, 4);
        (void)classifyGrid(invalidBoundaryGrid, bowTie);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "self-intersecting boundary rejected before classification");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-1 grid/classification tests passed\n";
    return EXIT_SUCCESS;
}
