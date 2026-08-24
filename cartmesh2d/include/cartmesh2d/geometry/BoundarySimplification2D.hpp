#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <optional>
#include <string>

namespace cartmesh2d {

struct BoundarySimplificationReport2D {
    double requestedDeviation = 0.0;
    double appliedDeviation = 0.0;
    double measuredMaxDeviation = 0.0;
    double originalArea = 0.0;
    double simplifiedArea = 0.0;
    std::size_t originalVertexCount = 0;
    std::size_t simplifiedVertexCount = 0;
    std::size_t attempts = 0;
};

struct BoundarySimplificationResult2D {
    std::optional<BoundaryRegion2D> boundary;
    BoundarySimplificationReport2D report;
    std::string error;

    [[nodiscard]] bool valid() const noexcept {
        return boundary.has_value() && error.empty();
    }
};

// Deterministic closed-loop Douglas-Peucker simplification. The two anchors
// are selected geometrically (lexicographic minimum and its farthest vertex),
// so cyclic input rotation does not change the result. If the requested
// tolerance changes loop/region validity, it is halved until the original
// topology is retained or the operation fails closed.
[[nodiscard]] BoundarySimplificationResult2D simplifyBoundaryRegion2D(
    const BoundaryRegion2D& source, double maximumDeviation,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
