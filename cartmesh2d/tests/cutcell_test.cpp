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

    // The default CFD convention is EXTERIOR fluid. This triangle occupies
    // the lower-left half of the unit cell, so the retained fluid is the
    // upper-right half.
    BoundaryLoop diagonalTriangle({{-1.0, -1.0}, {2.0, -1.0}, {-1.0, 2.0}});
    const auto diagonal = buildCutCell(unit, CellClass::Intersected, diagonalTriangle);
    check(diagonal.valid(), "default exterior diagonal cut is valid");
    check(diagonal.kind == CutCellKind::Cut, "default exterior case is a real cut cell");
    check(diagonal.fluidPolygon.vertices.size() == 3, "exterior diagonal cut produces a triangle");
    near(diagonal.area, 0.5, 1.0e-12, "exterior diagonal area = 0.5");
    near(diagonal.areaFraction, 0.5, 1.0e-12, "exterior diagonal area fraction = 0.5");
    check(diagonal.centroid.has_value(), "exterior diagonal cut has centroid");
    if (diagonal.centroid) {
        near(diagonal.centroid->x, 2.0 / 3.0, 1.0e-12, "exterior diagonal centroid x");
        near(diagonal.centroid->y, 2.0 / 3.0, 1.0e-12, "exterior diagonal centroid y");
    }
    check(diagonal.embeddedBoundary.size() == 1,
          "exterior diagonal cut has one embedded boundary fragment");
    if (!diagonal.embeddedBoundary.empty()) {
        const auto& fragment = diagonal.embeddedBoundary.front();
        check(pointOnInputBoundary(fragment.a, diagonalTriangle),
              "embedded fragment start lies on input boundary");
        check(pointOnInputBoundary(fragment.b, diagonalTriangle),
              "embedded fragment end lies on input boundary");
        near(std::hypot(fragment.b.x - fragment.a.x, fragment.b.y - fragment.a.y),
             std::sqrt(2.0), 1.0e-12, "diagonal embedded fragment length");
    }

    BoundaryLoop verticalSolid({{-1.0, -1.0}, {0.25, -1.0}, {0.25, 2.0}, {-1.0, 2.0}});
    const auto exteriorQuarter = buildCutCell(unit, CellClass::Intersected, verticalSolid);
    check(exteriorQuarter.valid() && exteriorQuarter.kind == CutCellKind::Cut,
          "default exterior vertical case is a valid cut cell");
    near(exteriorQuarter.area, 0.75, 1.0e-12, "exterior vertical fluid area = 0.75");
    near(exteriorQuarter.areaFraction, 0.75, 1.0e-12,
         "exterior vertical fluid fraction = 0.75");
    check(exteriorQuarter.centroid.has_value(), "exterior vertical cut has centroid");
    if (exteriorQuarter.centroid) {
        near(exteriorQuarter.centroid->x, 0.625, 1.0e-12, "exterior vertical centroid x");
        near(exteriorQuarter.centroid->y, 0.5, 1.0e-12, "exterior vertical centroid y");
    }
    check(exteriorQuarter.fluidPolygon.signedArea() > 0.0,
          "exterior cut polygon orientation is CCW");

    // Explicit internal-flow mode remains available, but it must never be the
    // default product interpretation.
    const auto interiorQuarter = buildCutCell(unit, CellClass::Intersected, verticalSolid,
                                              FluidRegion2D::Interior);
    check(interiorQuarter.valid() && interiorQuarter.kind == CutCellKind::Cut,
          "explicit interior vertical case is a valid cut cell");
    near(interiorQuarter.area, 0.25, 1.0e-12, "interior vertical fluid area = 0.25");
    near(interiorQuarter.areaFraction, 0.25, 1.0e-12,
         "interior vertical fluid fraction = 0.25");
    if (interiorQuarter.centroid) {
        near(interiorQuarter.centroid->x, 0.125, 1.0e-12, "interior vertical centroid x");
        near(interiorQuarter.centroid->y, 0.5, 1.0e-12, "interior vertical centroid y");
    } else {
        check(false, "interior vertical cut has centroid");
    }

    // Geometric INSIDE means solid for the default exterior CFD product.
    const auto solidInterior = buildCutCell(unit, CellClass::Inside, verticalSolid);
    check(solidInterior.valid() && solidInterior.kind == CutCellKind::Empty,
          "geometric inside leaf is empty solid in default exterior mode");
    near(solidInterior.areaFraction, 0.0, 1.0e-12,
         "solid interior has zero fluid fraction");

    // Geometric OUTSIDE means full fluid for the default exterior CFD product.
    const auto exteriorFull = buildCutCell(unit, CellClass::Outside, verticalSolid);
    check(exteriorFull.valid() && exteriorFull.kind == CutCellKind::Full,
          "geometric outside leaf is full fluid in default exterior mode");
    near(exteriorFull.area, 1.0, 1.0e-12, "exterior full cell area");
    near(exteriorFull.areaFraction, 1.0, 1.0e-12, "exterior full cell alpha");
    check(exteriorFull.fluidPolygon.vertices.size() == 4,
          "exterior full cell stores real rectangle polygon");

    const auto explicitInteriorFull = buildCutCell(unit, CellClass::Inside, verticalSolid,
                                                   FluidRegion2D::Interior);
    check(explicitInteriorFull.valid() && explicitInteriorFull.kind == CutCellKind::Full,
          "explicit interior mode still keeps geometric inside as fluid");

    BoundaryLoop tangent({{-1.0, 1.0}, {0.0, 1.0}, {-1.0, 2.0}});
    const auto tangentResult = buildCutCell(unit, CellClass::Intersected, tangent);
    check(tangentResult.valid() && tangentResult.kind == CutCellKind::Full,
          "zero-area solid tangent contact leaves the exterior fluid cell full");

    BoundaryLoop clockwise({{-1.0, 2.0}, {0.25, 2.0}, {0.25, -1.0}, {-1.0, -1.0}});
    const auto normalized = buildCutCell(unit, CellClass::Intersected, clockwise);
    check(normalized.valid() && normalized.kind == CutCellKind::Cut,
          "clockwise valid input is normalized safely in exterior mode");
    near(normalized.area, 0.75, 1.0e-12, "clockwise normalized exterior area");
    check(normalized.fluidPolygon.signedArea() > 0.0,
          "normalized exterior cut polygon is CCW");

    BoundaryLoop bowTie({{-1.0, -1.0}, {2.0, 2.0}, {-1.0, 2.0}, {2.0, -1.0}});
    const auto invalidBoundary = buildCutCell(unit, CellClass::Intersected, bowTie);
    check(!invalidBoundary.valid() && invalidBoundary.kind == CutCellKind::Unsupported,
          "self-intersecting boundary is rejected");

    // Keep the old multi-component geometry regression in explicit Interior
    // mode so it tests topology robustness rather than product-side semantics.
    BoundaryLoop disconnectedIntersection({
        {0.82, -0.30}, {0.78, -0.30}, {0.22, -0.30}, {0.18, -0.30},
        {0.18, -0.20}, {0.18, 0.20}, {0.10, 0.20}, {0.10, 0.80},
        {0.30, 0.80}, {0.30, 0.20}, {0.22, 0.20}, {0.22, -0.20},
        {0.78, -0.20}, {0.78, 0.20}, {0.70, 0.20}, {0.70, 0.80},
        {0.90, 0.80}, {0.90, 0.20}, {0.82, 0.20}, {0.82, -0.20}});
    const auto multi = buildCutCell(unit, CellClass::Intersected, disconnectedIntersection,
                                    FluidRegion2D::Interior);
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
         "quadtree wrapper preserves exterior cut geometry");
    if (fromLeaf.centroid) {
        near(fromLeaf.centroid->x, 2.0 / 3.0, 1.0e-12,
             "quadtree wrapper uses exterior fluid side");
    } else {
        check(false, "quadtree exterior cut has centroid");
    }

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-3 exterior-fluid Cut-cell tests passed\n";
    return EXIT_SUCCESS;
}
