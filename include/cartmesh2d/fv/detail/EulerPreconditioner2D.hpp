#pragma once
#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"

namespace cartmesh2d::fv::detail {
using EulerMatrix4 = std::array<std::array<double,4>,4>;
// Physical inviscid flux Jacobian integrated over the actual oriented edge.
[[nodiscard]] EulerMatrix4 eulerNormalFluxJacobian2D(const EulerConservative2D&,
    Vector2D area,const IdealGas2D&);

// Approximate RIGHT preconditioner only. The Newton product, line search and
// accepted time update still evaluate the original conservative residual.
// Frozen Rusanov face speed + normal two-point heat/Newtonian stress, scaled
// with the same fixed conservative variables as the Newton stage; ILU(0)
// preserves the actual cell-adjacency graph, with no diagonal shifts or pivot
// replacement. This is not a new transport discretization or physical model.
class EulerFrozenPreconditioner2D {
public:
    using BoundaryMap=std::function<EulerConservative2D(const EulerConservative2D&,
        const EulerBoundary2D&,Vector2D outwardUnitNormal)>;
    EulerFrozenPreconditioner2D(const FvMesh2D&,
        const std::vector<const EulerBoundary2D*>& lookup,const IdealGas2D&,
        const EulerTransport2D&,const std::vector<EulerConservative2D>&,
        double stageStep,const LinearVector2D& scales,const BoundaryMap&);
    EulerFrozenPreconditioner2D(const EulerFrozenPreconditioner2D&)=delete;
    EulerFrozenPreconditioner2D& operator=(const EulerFrozenPreconditioner2D&)=delete;
    [[nodiscard]] LinearVector2D apply(const LinearVector2D&) const;
    [[nodiscard]] LinearVector2D applyMatrix(const LinearVector2D&) const;
    [[nodiscard]] std::size_t nonzeros() const{return system_.diag.size()+system_.off.size();}
private:
    SparsePattern2D pattern_;
    SparseSystem2D system_;
};
} // namespace cartmesh2d::fv::detail
