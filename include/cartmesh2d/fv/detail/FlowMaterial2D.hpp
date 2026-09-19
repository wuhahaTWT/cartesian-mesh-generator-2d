#pragma once
#include "cartmesh2d/fv/Incompressible2D.hpp"

namespace cartmesh2d::fv::detail {
// Internal constitutive coupling point. These are the ACTUAL current velocity
// constraints, not another geometric boundary classifier in a turbulence model.
struct MaterialBoundary2D {
    Vector2D velocity;
    bool fixedX=false, fixedY=false, wall=false, inlet=false, outlet=false;
};
using MaterialUpdate2D = std::function<std::vector<double>(
    const FlowResult2D&, const std::vector<MaterialBoundary2D>&)>;
// Called after each pressure/flux correction, before current-state momentum
// assembly and acceptance. Returned coefficients must cover all faces and be
// finite positive. Exceptions abort; the old laminar entry points use no hook.
[[nodiscard]] FlowResult2D solveMaterialFlow2D(const FvMesh2D&,const FlowControls2D&,
    const MaterialUpdate2D&,const std::function<void(const FlowIteration2D&)>& = {});
}
