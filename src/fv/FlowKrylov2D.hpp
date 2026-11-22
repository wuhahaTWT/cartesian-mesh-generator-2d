#pragma once
#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <sstream>

namespace cartmesh2d::fv::solver_detail {

// Flexible right-preconditioned GMRES. Reorthogonalization and the final
// explicitly recomputed residual are mandatory: the inner pressure solve is
// deliberately inexact and therefore is not a fixed linear preconditioner.
class CoupledKrylov {
    const std::size_t size_,width_;
    std::vector<Vec> basis_,preconditioned_;
    Vec product_,residual_,work_,hessenberg_,cosine_,sine_,right_,solution_;
public:
    double residualNorm=0,maximumResidual=0;
    explicit CoupledKrylov(std::size_t size):size_(size),width_(std::min<std::size_t>(100,size)),
        basis_(width_+1,Vec(size)),preconditioned_(width_,Vec(size)),product_(size),
        residual_(size),work_(size),hessenberg_((width_+1)*width_),cosine_(width_),
        sine_(width_),right_(width_+1),solution_(width_) {}

    template<class Apply,class Precondition>
    void solve(Apply&& apply,Precondition&& precondition,Vec& x,const Vec& rhs,
               double relativeTolerance,double maximumRowResidual,FlowPerformance2D& performance) {
        constexpr std::size_t maximumIterations=2000;
        const double rhsNorm=detail::linearNorm(rhs),target=relativeTolerance*rhsNorm;
        const auto measure=[&] {
            apply(x,product_);double maximum=0;
            for(std::size_t i=0;i<size_;++i) {
                residual_[i]=rhs[i]-product_[i];maximum=std::max(maximum,std::abs(residual_[i]));
            }
            maximumResidual=maximum;residualNorm=detail::linearNorm(residual_);
            return residualNorm<=target && maximum<=maximumRowResidual;
        };
        // The exactly homogeneous affine system has the exact zero solution;
        // this also avoids demanding a relative residual against a zero RHS.
        if(rhsNorm==0)std::fill(x.begin(),x.end(),0.);
        if(measure())return;
        std::size_t iterations=0;
        while(iterations<maximumIterations) {
            for(std::size_t i=0;i<size_;++i)basis_[0][i]=residual_[i]/residualNorm;
            std::fill(hessenberg_.begin(),hessenberg_.end(),0.);
            std::fill(right_.begin(),right_.end(),0.);right_[0]=residualNorm;
            std::size_t used=0;
            for(std::size_t j=0;j<width_ && iterations<maximumIterations;++j) {
                precondition(basis_[j],preconditioned_[j]);apply(preconditioned_[j],work_);
                ++iterations;++performance.coupledIterations;
                for(unsigned pass=0;pass<2;++pass)for(std::size_t k=0;k<=j;++k) {
                    const double projection=detail::linearProduct(work_,basis_[k]);
                    hessenberg_[k*width_+j]+=projection;
                    for(std::size_t i=0;i<size_;++i)work_[i]-=projection*basis_[k][i];
                }
                const double next=detail::linearNorm(work_);hessenberg_[(j+1)*width_+j]=next;
                if(next>0)for(std::size_t i=0;i<size_;++i)basis_[j+1][i]=work_[i]/next;
                for(std::size_t k=0;k<j;++k) {
                    const double a=hessenberg_[k*width_+j],b=hessenberg_[(k+1)*width_+j];
                    hessenberg_[k*width_+j]=cosine_[k]*a+sine_[k]*b;
                    hessenberg_[(k+1)*width_+j]=-sine_[k]*a+cosine_[k]*b;
                }
                const double a=hessenberg_[j*width_+j],b=hessenberg_[(j+1)*width_+j];
                const double length=std::hypot(a,b);ensure(length>0,"Coupled FGMRES breakdown");
                cosine_[j]=a/length;sine_[j]=b/length;
                hessenberg_[j*width_+j]=length;hessenberg_[(j+1)*width_+j]=0;
                right_[j+1]=-sine_[j]*right_[j];right_[j]*=cosine_[j];used=j+1;
                if(next==0 || std::abs(right_[j+1])<target*.5)break;
            }
            for(std::size_t j=used;j-->0;) {
                double value=right_[j];
                for(std::size_t k=j+1;k<used;++k)value-=hessenberg_[j*width_+k]*solution_[k];
                solution_[j]=finite(value/hessenberg_[j*width_+j]);
            }
            for(std::size_t j=0;j<used;++j)for(std::size_t i=0;i<size_;++i)
                x[i]+=solution_[j]*preconditioned_[j][i];
            if(measure())return;
        }
        std::ostringstream message;
        message<<"Coupled FGMRES iteration limit: true scaled residual="<<residualNorm<<", target="<<target
            <<", maximum row residual="<<maximumResidual<<", row target="<<maximumRowResidual;
        throw std::runtime_error(message.str());
    }
};

} // namespace cartmesh2d::fv::solver_detail
