#pragma once

#include "cartmesh/geometry/Triangle.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace cartmesh {

enum class SurfaceFormat : std::uint8_t {
    ascii_stl,
    binary_stl,
};

class SurfaceMesh {
  public:
    SurfaceMesh(std::vector<Triangle> triangles, SurfaceFormat format, std::string name = {});

    [[nodiscard]] const std::vector<Triangle>& triangles() const noexcept { return triangles_; }
    [[nodiscard]] const AABB& bounds() const noexcept { return bounds_; }
    [[nodiscard]] SurfaceFormat format() const noexcept { return format_; }
    [[nodiscard]] const std::string& name() const noexcept { return name_; }

  private:
    std::vector<Triangle> triangles_;
    AABB bounds_;
    SurfaceFormat format_;
    std::string name_;
};

[[nodiscard]] constexpr const char* surface_format_name(SurfaceFormat format) noexcept {
    return format == SurfaceFormat::ascii_stl ? "ascii_stl" : "binary_stl";
}

} // 命名空间 cartmesh
