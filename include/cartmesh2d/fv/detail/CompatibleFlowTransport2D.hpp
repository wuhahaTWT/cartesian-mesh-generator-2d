#pragma once
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
namespace cartmesh2d::fv::detail::compatible {
// Conservative RT1 transport and its analytic advector derivative. Every
// original atomic face is retained; split the P1 normal flux at its zero.
// open is indexed by cell-local face order and enables the natural advective trace.
Mat transportMatrix(const FvMesh2D&,std::size_t,const P1Local&,const Lift&,
                    const Vec& advector,const std::vector<bool>& open,int order);
Mat transportAdvectorDerivative(const FvMesh2D&,std::size_t,const P1Local&,const Lift&,
                               const Vec& advector,const std::vector<bool>& open,int order);
}
