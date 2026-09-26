#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"

namespace cartmesh2d::fv {
enum class WallGradient2D { Linear, Quadratic };
// Geometry-only recovery of a gradient at a prescribed-value wall face.
// Input values are CENTROID POINT estimates, as in the existing second-order
// primitive-variable transport discretization, not exact volume averages.
// The target wall value anchors the constant term. Every sample is evaluated
// as (value - targetValue), so constants cancel before any weighted sum.
struct WallGradientStencil2D {
    struct Sample {
        std::size_t index=0;
        bool boundary=false;
        Vector2D weight{};
    };
    std::vector<Sample> samples;
    unsigned rings=0;
    double pivotRatio=0; // min/max diagonal of scaled, column-pivoted QR.
};
// Connectivity expansion never crosses a physical boundary. Periodic images
// are translated through reciprocal partners. Rank failure is explicit: an
// under-resolved stencil must not silently be advertised as quadratic.
[[nodiscard]] WallGradientStencil2D quadraticWallGradient2D(
    const FvMesh2D&,std::size_t face,const std::vector<bool>& prescribed,
    const std::vector<std::optional<std::size_t>>& partners);
} // namespace cartmesh2d::fv
