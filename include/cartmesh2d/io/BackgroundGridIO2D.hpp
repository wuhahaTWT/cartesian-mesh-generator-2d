#pragma once

#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include <filesystem>
#include <string>

namespace cartmesh2d {
// Full, uncut leaves including the solid interior. This is deliberately not
// CM2D solver topology: adaptive hanging nodes and immersed walls need their own
// discretization before these cells can be used as a fluid domain.
[[nodiscard]] bool writeBackgroundGrid2D(
    const Quadtree2D& tree, const BoundaryRegion2D& boundary,
    const std::filesystem::path& prefix, bool uniform, std::string* error = nullptr);
}
