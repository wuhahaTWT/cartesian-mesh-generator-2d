#pragma once

#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <functional>
#include <optional>

namespace cartmesh2d::fv::detail {

// A bounded GMRES subspace solves the Newton direction for a scaled coupled
// residual. No SPD assumption or normal equations. Two orthogonalization
// passes retain the small Krylov directions. The .1 forcing term is an inner
// work target, not a nonlinear/physical acceptance tolerance. The caller must
// test the actual nonlinear residual and retain its original convergence gates.
// Independent implementation of standard inexact Newton--Krylov principles:
// https://petsc.org/release/manual/snes/#inexact-newton-like-methods
inline std::optional<LinearVector2D> newtonKrylovDirection2D(
    const std::function<LinearVector2D(const LinearVector2D&)>& apply,
    const LinearVector2D& rhs, std::size_t maximumDirections=60) {
    const auto count=std::min(maximumDirections,rhs.size());
    const double beta=linearNorm(rhs);
    if(count==0 || beta==0)return {};
    std::vector<LinearVector2D> basis;
    auto first=rhs;for(auto& value:first)value/=beta;
    basis.push_back(std::move(first));
    std::vector<LinearVector2D> h(count+1,LinearVector2D(count));
    LinearVector2D g(count+1),cosine(count),sine(count);g[0]=beta;
    std::size_t used=0;
    for(std::size_t j=0;j<count;++j) {
        auto w=apply(basis[j]);
        linearEnsure(w.size()==rhs.size(),"Coupled Jacobian product size mismatch");
        for(int pass=0;pass<2;++pass)for(std::size_t i=0;i<=j;++i) {
            const double value=linearProduct(basis[i],w);h[i][j]+=value;
            for(std::size_t k=0;k<w.size();++k)w[k]-=value*basis[i][k];
        }
        const double remainder=linearNorm(w);h[j+1][j]=remainder;
        if(remainder>0) {
            for(auto& value:w)value/=remainder;
            basis.push_back(std::move(w));
        }
        for(std::size_t i=0;i<j;++i) {
            const double value=cosine[i]*h[i][j]+sine[i]*h[i+1][j];
            h[i+1][j]=-sine[i]*h[i][j]+cosine[i]*h[i+1][j];h[i][j]=value;
        }
        const double diagonal=std::hypot(h[j][j],h[j+1][j]);
        if(!(diagonal>0) || !std::isfinite(diagonal))return {};
        cosine[j]=h[j][j]/diagonal;sine[j]=h[j+1][j]/diagonal;
        h[j][j]=diagonal;h[j+1][j]=0;
        g[j+1]=-sine[j]*g[j];g[j]*=cosine[j];used=j+1;
        if(std::abs(g[j+1])<=.1*beta || remainder==0)break;
    }
    LinearVector2D y(used),solution(rhs.size());
    for(std::size_t i=used;i-- >0;) {
        double value=g[i];
        for(std::size_t j=i+1;j<used;++j)value-=h[i][j]*y[j];
        y[i]=linearFinite(value/h[i][i]);
    }
    for(std::size_t j=0;j<used;++j)for(std::size_t k=0;k<solution.size();++k)
        solution[k]=linearFinite(solution[k]+y[j]*basis[j][k]);
    return solution;
}
} // namespace cartmesh2d::fv::detail
