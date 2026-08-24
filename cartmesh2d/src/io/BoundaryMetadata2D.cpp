#include "cartmesh2d/io/BoundaryMetadata2D.hpp"

namespace cartmesh2d {

const char* boundaryConditionRoleName(BoundaryConditionRole2D role) noexcept {
    switch (role) {
    case BoundaryConditionRole2D::Wall: return "wall";
    case BoundaryConditionRole2D::Inlet: return "inlet";
    case BoundaryConditionRole2D::Outlet: return "outlet";
    case BoundaryConditionRole2D::Slip: return "slip";
    case BoundaryConditionRole2D::Symmetry: return "symmetry";
    }
    return "wall";
}

bool parseBoundaryConditionRole(const std::string& value,
                                BoundaryConditionRole2D& role) noexcept {
    if (value=="wall") role=BoundaryConditionRole2D::Wall;
    else if (value=="inlet") role=BoundaryConditionRole2D::Inlet;
    else if (value=="outlet") role=BoundaryConditionRole2D::Outlet;
    else if (value=="slip") role=BoundaryConditionRole2D::Slip;
    else if (value=="symmetry") role=BoundaryConditionRole2D::Symmetry;
    else return false;
    return true;
}

} // namespace cartmesh2d
