#pragma once

#include <string>

namespace cartmesh2d {

enum class BoundaryConditionRole2D {
    Wall,
    Inlet,
    Outlet,
    Slip,
    Symmetry
};

struct EmbeddedBoundaryPatch2D {
    std::string name;
    std::string type;
    BoundaryConditionRole2D role = BoundaryConditionRole2D::Wall;
    std::string sourceLayer;

    [[nodiscard]] bool valid() const noexcept {
        return !name.empty() && !type.empty();
    }
};

[[nodiscard]] const char* boundaryConditionRoleName(
    BoundaryConditionRole2D role) noexcept;

[[nodiscard]] bool parseBoundaryConditionRole(
    const std::string& value, BoundaryConditionRole2D& role) noexcept;

} // namespace cartmesh2d
