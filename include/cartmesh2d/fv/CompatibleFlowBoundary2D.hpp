#pragma once
#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
namespace cartmesh2d::fv {
// Translate validated existing explicit conditions. Values remain facewise
// constant. PressureOutlet keeps pseudo-traction and backflow rejection;
// Symmetry is zero normal velocity plus zero tangential physical traction.
// SmoothMovingWall and PressureOpening require their own higher-order/dynamic
// policies and are explicitly rejected; no silent boundary-kind substitution.
[[nodiscard]] std::vector<CompatibleBoundary2D> compatibleFlowBoundaries2D(
    const FvMesh2D&,const std::vector<FlowBoundaryCondition2D>&,double referenceSpeed);
}
