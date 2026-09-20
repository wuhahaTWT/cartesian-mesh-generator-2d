#pragma once
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <iosfwd>

namespace cartmesh2d::fv {

// Boundary configurations bind face IDs to owner, face centre and area vector
// on a specific final mesh. A matching face count alone is insufficient.
// Every boundary must occur once; internal faces and trailing data are rejected.
[[nodiscard]] std::vector<FlowBoundaryCondition2D> readFlowBoundaryConditions2D(
    std::istream&, const FvMesh2D&, const FlowControls2D&);
void writeFlowBoundaryConditions2D(std::ostream&, const FvMesh2D&, const FlowControls2D&);

// Template for concentric regular polygon approximations to circular walls
// (even numbers of sides >=16); inner wall moves CCW at the supplied speed.
// Rejects squares, eccentric/irregular rings and extra boundary components.
[[nodiscard]] std::vector<FlowBoundaryCondition2D> rotatingAnnulusBoundaryPreset2D(
    const FvMesh2D&, double innerSurfaceSpeed);

// Converts supported physical presets (channel, duct, cavity, annulus) to explicit data.
// Verification forcing and unsupported slip/farfield conditions are rejected.
[[nodiscard]] std::vector<FlowBoundaryCondition2D> explicitFlowBoundaryPreset2D(
    const FvMesh2D&, const FlowControls2D&);

[[nodiscard]] const char* flowBoundaryKindName2D(FlowBoundaryKind2D);
[[nodiscard]] FlowBoundaryKind2D flowBoundaryKindFromName2D(const std::string&);

}
