// R2/W0: top-level H4 attribution invariants.
//
// R1F section 7 withdrew every performance claim because the existing solver
// sub-phase timings could not explain the end-to-end wall time. These tests fix
// the *structure* of the new attribution so it cannot silently regress:
// the deterministic call counters must match the process counter, and the stage
// timers must never over-account for the total.
//
// Wall-clock magnitudes are deliberately not asserted; they are not
// reproducible and a loaded CI machine must not be able to fail this test.

#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

[[nodiscard]] BoundaryLoop circleLoop(std::size_t segments, double radius) {
    std::vector<Point2D> points;
    points.reserve(segments);
    for (std::size_t i = 0; i < segments; ++i) {
        const double angle = 2.0 * std::numbers::pi *
                             static_cast<double>(i) / static_cast<double>(segments);
        points.push_back({radius * std::cos(angle), radius * std::sin(angle)});
    }
    return BoundaryLoop(std::move(points));
}

// The counter must advance by exactly one per conformal hybrid construction,
// so a caller can prove build counts by measurement rather than by assertion.
void conformalBuildCounterAdvancesPerBuild() {
    const auto loop = circleLoop(24U, 1.0);
    BoundaryRegion2D walls({loop});
    auto chain = makeClosedWallChain2D(loop, 0U, "wall_0");
    check(chain.success(), "wall chain for the counter case must build");
    if (!chain.success()) return;

    LayerParameters2D layers;
    layers.nLayers = 2U;
    layers.thicknessMode = LayerThicknessMode2D::FirstLayerThickness;
    layers.thickness = 0.04;
    layers.growthRatio = 1.2;
    const auto strips = buildBoundaryLayerStrips2D({*chain.chain}, layers);
    check(strips.success(), "boundary-layer strips for the counter case must build");
    if (!strips.success()) return;

    const auto bounds = walls.bounds();
    const Domain2D domain{{{bounds.min.x - 1.0, bounds.min.y - 1.0},
                           {bounds.max.x + 1.0, bounds.max.y + 1.0}}};
    QuadtreeRefinementPolicy2D refinement;
    refinement.minimumLevel = 3U;
    refinement.boundaryLevel = 5U;

    const auto before = conformalHybridBuildCount2D();
    const auto result = buildAutomaticHybridWithConstruction2D(
        strips, domain, walls, 5U, refinement, HybridMeshPolicy2D{});
    const auto after = conformalHybridBuildCount2D();
    check(after > before,
          "one automatic hybrid attempt must consume at least one conformal build");
    check(result.metrics.transitionRingCount >= 3U,
          "the automatic transition plan keeps at least three rings");
}

// buildRobustH4Mesh2D must report exactly the conformal builds it consumed, and
// must never claim more attributed time than it measured in total.
void robustProfileMatchesMeasuredCounters() {
    const auto loop = circleLoop(24U, 1.0);
    BoundaryRegion2D walls({loop});
    auto chain = makeClosedWallChain2D(loop, 0U, "wall_0");
    check(chain.success(), "wall chain for the robust case must build");
    if (!chain.success()) return;

    LayerParameters2D layers;
    layers.nLayers = 2U;
    layers.thicknessMode = LayerThicknessMode2D::FirstLayerThickness;
    layers.thickness = 0.04;
    layers.growthRatio = 1.2;

    const auto bounds = walls.bounds();
    const Domain2D domain{{{bounds.min.x - 1.0, bounds.min.y - 1.0},
                           {bounds.max.x + 1.0, bounds.max.y + 1.0}}};
    QuadtreeRefinementPolicy2D refinement;
    refinement.minimumLevel = 3U;
    refinement.boundaryLevel = 5U;

    const auto before = conformalHybridBuildCount2D();
    const auto robust = buildRobustH4Mesh2D(
        {*chain.chain}, layers, domain, walls, 5U, refinement, {},
        HybridMeshPolicy2D{});
    const auto measured = conformalHybridBuildCount2D() - before;

    check(robust.profile.conformalHybridBuildCalls == measured,
          "the reported conformal build count must equal the measured counter delta");
    check(robust.profile.conformalHybridBuildCalls >= 1U,
          "a robust build always performs at least one conformal build");

    const double attributed = robust.profile.requestedLayerSeconds +
                              robust.profile.requestedHybridSeconds +
                              robust.profile.localLayerSeconds +
                              robust.profile.localHybridSeconds +
                              robust.profile.pureCutCellFallbackSeconds;
    check(attributed <= robust.profile.totalSeconds + 1.0e-9,
          "stage timers must not attribute more than the measured total");
    check(robust.profile.unattributedSeconds() >= 0.0,
          "unattributed seconds are never negative");
    check(std::isfinite(robust.profile.totalSeconds),
          "the total is a finite wall-time measurement");

    // Stage attempt counters must agree with the committed mode: a hybrid mesh
    // never runs the pure Cut-cell fallback.
    if (robust.mode == H4MeshMode2D::Hybrid) {
        check(robust.profile.pureCutCellFallbackAttempts == 0U,
              "a committed hybrid mesh must not have entered the fallback stage");
        check(robust.profile.requestedHybridAttempts +
                  robust.profile.localHybridAttempts >= 1U,
              "a committed hybrid mesh entered at least one hybrid attempt");
    }
    if (robust.profile.localHybridAttempts > 0U) {
        check(robust.profile.requestedHybridAttempts +
                  robust.profile.localHybridAttempts >= 2U ||
                  !robust.requestedLayerCandidate.success(),
              "the local attempt only runs after the requested candidate failed");
    }
}

} // namespace

int main() {
    conformalBuildCounterAdvancesPerBuild();
    robustProfileMatchesMeasuredCounters();
    if (failures != 0) {
        std::cerr << failures << " refinement-attribution checks failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "refinement attribution tests passed\n";
    return EXIT_SUCCESS;
}
