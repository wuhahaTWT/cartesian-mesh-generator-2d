#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"

namespace cartmesh2d::fv {
struct WallDistanceResult2D {
    std::vector<double> distance;
    std::vector<std::size_t> nearestFace;
    std::size_t segmentTests = 0; // actual leaf distance evaluations, diagnostic
};
// Nearest finite straight wall segment, not face-centre or supporting-line
// distance. Caller explicitly identifies physical walls; inlet/outlet/slip
// boundaries are not inferred from geometry patches. At least one wall required.
// Uses the validated final FV segment endpoints recovered as C +/- tangent(S)/2.
// Deterministic AABB hierarchy; no approximate search or distance floor.
[[nodiscard]] WallDistanceResult2D computeWallDistance2D(
    const FvMesh2D&, const std::vector<bool>& wallFaces);
}
