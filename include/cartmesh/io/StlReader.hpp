#pragma once

#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <filesystem>

namespace cartmesh {

[[nodiscard]] SurfaceMesh read_stl(const std::filesystem::path& path);

} // 命名空间 cartmesh
