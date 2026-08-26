#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
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

void near(double actual, double expected, double tolerance,
          const std::string& message) {
    check(std::abs(actual - expected) <= tolerance,
          message + " actual=" + std::to_string(actual) +
          " expected=" + std::to_string(expected));
}

std::vector<Point2D> ellipse(std::size_t count, double a, double b) {
    std::vector<Point2D> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double angle = 2.0 * std::numbers::pi *
                             static_cast<double>(i) / static_cast<double>(count);
        points.push_back({a * std::cos(angle), b * std::sin(angle)});
    }
    return points;
}

LayerParameters2D firstThickness(std::size_t layers, double first,
                                 double ratio) {
    return {layers, LayerThicknessMode2D::FirstLayerThickness, first, ratio};
}

bool sameStrip(const BoundaryLayerStrip2D& lhs,
               const BoundaryLayerStrip2D& rhs) {
    if (lhs.vertices.size() != rhs.vertices.size() ||
        lhs.cells.size() != rhs.cells.size()) return false;
    for (std::size_t i = 0; i < lhs.vertices.size(); ++i) {
        if (lhs.vertices[i].id != rhs.vertices[i].id ||
            lhs.vertices[i].ring != rhs.vertices[i].ring ||
            lhs.vertices[i].chainVertex != rhs.vertices[i].chainVertex ||
            lhs.vertices[i].point.x != rhs.vertices[i].point.x ||
            lhs.vertices[i].point.y != rhs.vertices[i].point.y) return false;
    }
    for (std::size_t i = 0; i < lhs.cells.size(); ++i) {
        if (lhs.cells[i].id != rhs.cells[i].id ||
            lhs.cells[i].layer != rhs.cells[i].layer ||
            lhs.cells[i].wallSegment != rhs.cells[i].wallSegment ||
            lhs.cells[i].vertices != rhs.cells[i].vertices ||
            lhs.cells[i].area != rhs.cells[i].area) return false;
    }
    return true;
}

bool samePoints(const std::vector<Point2D>& lhs,
                const std::vector<Point2D>& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i].x != rhs[i].x || lhs[i].y != rhs[i].y) return false;
    }
    return true;
}

} // namespace

int main() {
    const TolerancePolicy tolerance{};

    const auto resolvedOne = resolveLayerParameters2D(
        {4U, LayerThicknessMode2D::FirstLayerThickness, 0.1, 1.0});
    check(resolvedOne.success(), "growth ratio one resolves without cancellation");
    if (resolvedOne.success()) {
        near(resolvedOne.parameters->totalThickness, 0.4, 1.0e-14,
             "ratio-one total thickness");
    }
    const auto resolvedTotal = resolveLayerParameters2D(
        {3U, LayerThicknessMode2D::TotalThickness, 0.7, 2.0});
    check(resolvedTotal.success(), "total-thickness parameterization resolves");
    if (resolvedTotal.success()) {
        near(resolvedTotal.parameters->firstLayerThickness, 0.1, 1.0e-14,
             "total thickness converts to first thickness");
    }
    check(!resolveLayerParameters2D(
        {0U, LayerThicknessMode2D::FirstLayerThickness, 0.1, 1.0}).success(),
        "zero layer count fails closed");
    check(!resolveLayerParameters2D(
        {2U, LayerThicknessMode2D::FirstLayerThickness,
         std::numeric_limits<double>::quiet_NaN(), 1.0}).success(),
        "NaN thickness fails closed");

    const auto open = makeOpenWallChain2D(
        {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}, 7U, "open_wall",
        FluidSide2D::Left, tolerance);
    check(open.success(), "ordered straight open chain is valid");
    if (open.success()) {
        const auto result = buildBoundaryLayerStrip2D(
            *open.chain, firstThickness(3U, 0.1, 2.0));
        check(result.success(), "straight open chain produces a strip");
        if (result.success()) {
            const auto& strip = result.strips.front();
            check(strip.metrics.cellCount == 6U && strip.metrics.vertexCount == 12U,
                  "open strip has deterministic cell and vertex counts");
            near(strip.vertices[strip.ringVertexIds[1][1]].point.y, 0.1, 1.0e-14,
                 "first open layer advances to fluid side");
            near(strip.vertices[strip.ringVertexIds[3][1]].point.y, 0.7, 1.0e-14,
                 "open hair-edge geometric spacing reaches total thickness");
        }
    }

    const auto circlePoints = ellipse(64U, 1.0, 1.0);
    const auto circleChain = makeClosedWallChain2D(
        BoundaryLoop(circlePoints), 0U, "circle_wall");
    check(circleChain.success(), "circle wall chain is valid");
    const auto circleResult = circleChain.success()
        ? buildBoundaryLayerStrip2D(*circleChain.chain,
                                    firstThickness(4U, 0.02, 1.2))
        : BoundaryLayerBuildResult2D{};
    check(circleResult.success(), "convex circle produces four complete layers");
    if (circleResult.success()) {
        const auto& strip = circleResult.strips.front();
        check(strip.wallChain.orientation == WallChainOrientation2D::CounterClockwise &&
              strip.wallChain.fluidSide == FluidSide2D::Right,
              "exterior circle keeps solid left and fluid right");
        check(strip.metrics.cellCount == 256U && strip.metrics.vertexCount == 320U,
              "circle layer topology counts are exact");
        check(strip.metrics.minCellArea > 0.0 &&
              strip.metrics.maxCellArea >= strip.metrics.minCellArea,
              "circle layer quads have positive audited areas");
        check(strip.outerEnvelope().size() == circlePoints.size(),
              "circle outer envelope is a complete closed ring");
    }

    auto rotatedCircle = circlePoints;
    std::rotate(rotatedCircle.begin(), rotatedCircle.begin() + 17,
                rotatedCircle.end());
    std::reverse(rotatedCircle.begin(), rotatedCircle.end());
    const auto canonicalCircle = makeClosedWallChain2D(
        BoundaryLoop(rotatedCircle), 0U, "circle_wall");
    check(canonicalCircle.success() && circleChain.success() &&
          samePoints(canonicalCircle.chain->vertices, circleChain.chain->vertices),
          "closed wall ordering is invariant to cyclic start and input direction");
    if (canonicalCircle.success() && circleResult.success()) {
        const auto repeated = buildBoundaryLayerStrip2D(
            *canonicalCircle.chain, firstThickness(4U, 0.02, 1.2));
        check(repeated.success() && sameStrip(repeated.strips.front(),
                                             circleResult.strips.front()),
              "canonical input produces byte-value deterministic strip IDs and geometry");
    }

    const auto ellipseChain = makeClosedWallChain2D(
        BoundaryLoop(ellipse(96U, 2.0, 1.0)), 1U, "ellipse_wall");
    const auto ellipseResult = ellipseChain.success()
        ? buildBoundaryLayerStrip2D(*ellipseChain.chain,
                                    firstThickness(5U, 0.01, 1.15))
        : BoundaryLayerBuildResult2D{};
    check(ellipseResult.success(),
          "non-uniform-curvature convex ellipse produces a legal strip");

    const auto rectangleChain = makeClosedWallChain2D(
        BoundaryLoop({{-2.0, -1.0}, {2.0, -1.0},
                      {2.0, 1.0}, {-2.0, 1.0}}),
        2U, "rectangle_wall");
    const auto rectangleResult = rectangleChain.success()
        ? buildBoundaryLayerStrip2D(*rectangleChain.chain,
                                    firstThickness(3U, 0.04, 1.25))
        : BoundaryLayerBuildResult2D{};
    check(rectangleResult.success(), "mild 90-degree convex corners are supported");
    if (rectangleResult.success()) {
        check(std::all_of(rectangleResult.strips.front().wallVertexKinds.begin(),
                          rectangleResult.strips.front().wallVertexKinds.end(),
                          [](WallVertexKind2D kind) {
                              return kind == WallVertexKind2D::MildConvex;
                          }), "rectangle corners are explicitly classified mild-convex");
    }

    const auto concaveChain = makeClosedWallChain2D(
        BoundaryLoop({{0.0, 0.0}, {3.0, 0.0}, {3.0, 3.0},
                      {2.0, 3.0}, {2.0, 1.0}, {1.0, 1.0},
                      {1.0, 3.0}, {0.0, 3.0}}),
        3U, "concave_wall");
    const auto concaveResult = concaveChain.success()
        ? buildBoundaryLayerStrip2D(*concaveChain.chain,
                                    firstThickness(3U, 0.02, 1.1))
        : BoundaryLayerBuildResult2D{};
    check(!concaveResult.success() &&
          concaveResult.failure.reason == BoundaryLayerFailureReason2D::ConcaveCorner,
          "severe concave geometry fails closed with a structured reason");

    const auto sharpChain = makeClosedWallChain2D(
        BoundaryLoop({{0.0, 0.0}, {3.0, 0.0}, {0.05, 0.1}}),
        4U, "sharp_wall");
    const auto sharpResult = sharpChain.success()
        ? buildBoundaryLayerStrip2D(*sharpChain.chain,
                                    firstThickness(2U, 0.01, 1.0))
        : BoundaryLayerBuildResult2D{};
    check(!sharpResult.success() &&
          sharpResult.failure.reason == BoundaryLayerFailureReason2D::SharpCorner,
          "trailing-edge-like sharp corner fails closed");

    const auto leftGap = makeClosedWallChain2D(
        BoundaryLoop({{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}),
        10U, "left_gap_wall");
    const auto rightGap = makeClosedWallChain2D(
        BoundaryLoop({{1.2, 0.0}, {2.2, 0.0}, {2.2, 1.0}, {1.2, 1.0}}),
        11U, "right_gap_wall");
    check(leftGap.success() && rightGap.success(), "narrow-gap wall chains are valid inputs");
    if (leftGap.success() && rightGap.success()) {
        const auto savedLeft = leftGap.chain->vertices;
        const auto gapResult = buildBoundaryLayerStrips2D(
            {*leftGap.chain, *rightGap.chain}, firstThickness(1U, 0.11, 1.0));
        check(!gapResult.success() &&
              gapResult.failure.reason ==
                  BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit &&
              gapResult.failure.safeThickness.has_value(),
              "requested narrow-gap thickness fails against an explicit safe limit");
        check(samePoints(leftGap.chain->vertices, savedLeft),
              "failed candidate leaves the existing wall chain unchanged");
    }

    const auto degenerate = makeOpenWallChain2D(
        {{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}}, 20U,
        "degenerate_wall", FluidSide2D::Right);
    check(!degenerate.success() &&
          (degenerate.failureReason == WallChainFailureReason2D::DuplicateVertex ||
           degenerate.failureReason == WallChainFailureReason2D::DegenerateSegment),
          "duplicate consecutive point is rejected before layer construction");

    if (failures == 0) {
        std::cout << "cartmesh2d H4-1 boundary-layer tests: PASS\n";
        return EXIT_SUCCESS;
    }
    std::cerr << failures << " H4-1 test(s) failed\n";
    return EXIT_FAILURE;
}
