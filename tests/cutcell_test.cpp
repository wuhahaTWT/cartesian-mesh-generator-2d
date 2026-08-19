#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

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

bool pointOnInputBoundary(const Point2D& point, const BoundaryLoop& boundary,
                          const TolerancePolicy& tol = {}) {
    const auto& vertices = boundary.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (pointOnSegment(point, {vertices[i], vertices[(i + 1) % vertices.size()]}, tol)) {
            return true;
        }
    }
    return false;
}
} // namespace

int main() {
    const AABB2D unit{{0.0, 0.0}, {1.0, 1.0}};

    BoundaryLoop diagonalTriangle({{-1.0, -1.0}, {2.0, -1.0}, {-1.0, 2.0}});
    const auto diagonal = buildCutCell(unit, CellClass::Intersected, diagonalTriangle);
    check(diagonal.valid(), "diagonal analytic cut is valid");
    check(diagonal.kind == CutCellKind::Cut, "diagonal analytic case is a real cut cell");
    check(diagonal.fluidPolygon.vertices.size() == 3, "diagonal cut produces a triangle");
    near(diagonal.area, 0.5, 1.0e-12, "diagonal triangle area = 0.5");
    near(diagonal.areaFraction, 0.5, 1.0e-12, "diagonal area fraction = 0.5");
    check(diagonal.centroid.has_value(), "diagonal cut has centroid");
    if (diagonal.centroid) {
        near(diagonal.centroid->x, 1.0 / 3.0, 1.0e-12, "diagonal centroid x");
        near(diagonal.centroid->y, 1.0 / 3.0, 1.0e-12, "diagonal centroid y");
    }
    check(diagonal.embeddedBoundary.size() == 1,
          "diagonal cut has one embedded boundary fragment");
    if (!diagonal.embeddedBoundary.empty()) {
        const auto& fragment = diagonal.embeddedBoundary.front();
        check(pointOnInputBoundary(fragment.a, diagonalTriangle),
              "embedded fragment start lies on input boundary");
        check(pointOnInputBoundary(fragment.b, diagonalTriangle),
              "embedded fragment end lies on input boundary");
        near(std::hypot(fragment.b.x - fragment.a.x, fragment.b.y - fragment.a.y),
             std::sqrt(2.0), 1.0e-12, "diagonal embedded fragment length");
    }

    BoundaryLoop verticalCut({{-1.0, -1.0}, {0.25, -1.0}, {0.25, 2.0}, {-1.0, 2.0}});
    const auto quarter = buildCutCell(unit, CellClass::Intersected, verticalCut);
    check(quarter.valid() && quarter.kind == CutCellKind::Cut,
          "vertical analytic case is a valid cut cell");
    near(quarter.area, 0.25, 1.0e-12, "vertical cut area = 0.25");
    near(quarter.areaFraction, 0.25, 1.0e-12, "vertical area fraction = 0.25");
    check(quarter.centroid.has_value(), "vertical cut has centroid");
    if (quarter.centroid) {
        near(quarter.centroid->x, 0.125, 1.0e-12, "vertical centroid x");
        near(quarter.centroid->y, 0.5, 1.0e-12, "vertical centroid y");
    }
    check(quarter.fluidPolygon.signedArea() > 0.0,
          "cut polygon orientation is CCW");

    const auto full = buildCutCell(unit, CellClass::Inside, verticalCut);
    check(full.valid() && full.kind == CutCellKind::Full, "inside leaf becomes full fluid cell");
    near(full.area, 1.0, 1.0e-12, "full cell area");
    near(full.areaFraction, 1.0, 1.0e-12, "full cell alpha");
    check(full.fluidPolygon.vertices.size() == 4, "full cell stores real rectangle polygon");

    const auto empty = buildCutCell(unit, CellClass::Outside, verticalCut);
    check(empty.valid() && empty.kind == CutCellKind::Empty, "outside leaf becomes empty");
    near(empty.areaFraction, 0.0, 1.0e-12, "empty alpha = 0");
    check(empty.fluidPolygon.vertices.empty(), "empty cell has no fake polygon");

    BoundaryLoop tangent({{-1.0, 1.0}, {0.0, 1.0}, {-1.0, 2.0}});
    const auto tangentResult = buildCutCell(unit, CellClass::Intersected, tangent);
    check(tangentResult.valid() && tangentResult.kind == CutCellKind::Empty,
          "zero-area tangent contact does not create a fake cut polygon");

    BoundaryLoop clockwise({{-1.0, 2.0}, {0.25, 2.0}, {0.25, -1.0}, {-1.0, -1.0}});
    const auto normalized = buildCutCell(unit, CellClass::Intersected, clockwise);
    check(normalized.valid() && normalized.kind == CutCellKind::Cut,
          "clockwise valid input is normalized safely");
    near(normalized.area, 0.25, 1.0e-12, "clockwise normalized area");
    check(normalized.fluidPolygon.signedArea() > 0.0,
          "normalized cut polygon is CCW");

    BoundaryLoop bowTie({{-1.0, -1.0}, {2.0, 2.0}, {-1.0, 2.0}, {2.0, -1.0}});
    const auto invalidBoundary = buildCutCell(unit, CellClass::Intersected, bowTie);
    check(!invalidBoundary.valid() && invalidBoundary.kind == CutCellKind::Unsupported,
          "self-intersecting boundary is rejected");

    BoundaryLoop disconnectedIntersection({
        {0.82, -0.30}, {0.78, -0.30}, {0.22, -0.30}, {0.18, -0.30},
        {0.18, -0.20}, {0.18, 0.20}, {0.10, 0.20}, {0.10, 0.80},
        {0.30, 0.80}, {0.30, 0.20}, {0.22, 0.20}, {0.22, -0.20},
        {0.78, -0.20}, {0.78, 0.20}, {0.70, 0.20}, {0.70, 0.80},
        {0.90, 0.80}, {0.90, 0.20}, {0.82, 0.20}, {0.82, -0.20}});
    const auto multi = buildCutCell(unit, CellClass::Intersected, disconnectedIntersection);
    check(!multi.valid() && multi.kind == CutCellKind::Unsupported,
          "multi-component cell intersection fails explicitly");
    bool sawMultiComponentIssue = false;
    for (const auto& issue : multi.issues) {
        if (issue.code == CutCellIssueCode::MultipleEmbeddedComponents ||
            issue.code == CutCellIssueCode::DegeneratePolygon) {
            sawMultiComponentIssue = true;
        }
    }
    check(sawMultiComponentIssue,
          "multi-component failure reports a topology/geometry reason");

    QuadtreeLeaf2D leaf;
    leaf.id = 17;
    leaf.key = 12345;
    leaf.bounds = unit;
    leaf.classification = CellClass::Intersected;
    const auto fromLeaf = buildCutCell(leaf, diagonalTriangle);
    check(fromLeaf.sourceId == 17 && fromLeaf.sourceKey == 12345,
          "quadtree source identity is preserved");
    near(fromLeaf.areaFraction, 0.5, 1.0e-12,
         "quadtree wrapper preserves cut geometry");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-3 cut-cell tests passed\n";
    return EXIT_SUCCESS;
}
