#pragma once
#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"
#include "cartmesh2d/grid/CartesianGrid2D.hpp"

namespace cartmesh2d {

struct TransitionCanonicalizationPolicy2D {
    // Kept equal to the Q1 hard face/local_h limit by a regression test.
    double minimumFaceFraction = 0.01;
    double featureTurnRadians = 0.523598775598298873;
    double clearanceFraction = 0.2;
};

// Only the mutable outer transition front is resampled. The wall and H4 layer
// interface are immutable. All changed common edges are replaced on both sides
// before Cut-cell clipping; this is not a topology repair or cell deletion.
[[nodiscard]] bool canonicalizeTransitionEnvelope2D(
    std::vector<BoundaryLoop>& remainderBoundaryLoops,
    std::vector<Polygon2D>& transitionPolygons,
    const BoundaryRegion2D& outerRegion, const Domain2D& domain,
    std::size_t boundaryLevel, double transitionRingThickness,
    bool localTermination, IntersectionRegistry2D& envelopeRegistry,
    const TolerancePolicy& tolerance, std::string& error);

} // namespace cartmesh2d
