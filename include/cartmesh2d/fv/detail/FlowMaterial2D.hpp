#pragma once
#include "cartmesh2d/fv/Incompressible2D.hpp"

namespace cartmesh2d::fv::detail {
// Internal constitutive coupling point. These are the ACTUAL current velocity
// constraints, not another geometric boundary classifier in a turbulence model.
struct MaterialBoundary2D {
    Vector2D velocity;
    bool fixedX=false, fixedY=false, wall=false, inlet=false, outlet=false;
};
struct MaterialState2D {
    std::vector<double> faceViscosity;
    bool converged=false; // current constitutive equations, not a change heuristic
};
using MaterialUpdate2D = std::function<MaterialState2D(
    const FlowResult2D&, const std::vector<MaterialBoundary2D>&)>;
// Called after each pressure/flux correction, before current-state momentum
// assembly and acceptance. Returned coefficients must cover all faces and be
// finite positive. Unconverged constitutive updates keep SIMPLE running even
// when momentum has converged. Exceptions abort; laminar entry points use no hook.
[[nodiscard]] FlowResult2D solveMaterialFlow2D(const FvMesh2D&,const FlowControls2D&,
    const MaterialUpdate2D&,const std::function<void(const FlowIteration2D&)>& = {});
}
