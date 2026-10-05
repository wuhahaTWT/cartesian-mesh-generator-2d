#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace cartmesh2d::fv::detail {
using NewtonVector = std::vector<double>;
inline double newtonDot(const NewtonVector& a,const NewtonVector& b) {
    long double sum=0;for(std::size_t i=0;i<a.size();++i)sum+=static_cast<long double>(a[i])*b[i];
    return static_cast<double>(sum);
}
inline double newtonNorm(const NewtonVector& a) {return std::sqrt(newtonDot(a,a));}
// Restarted right-preconditioned GMRES, two-pass modified Gram-Schmidt.
// Return only after recomputing the actual linear residual; no false convergence
// from the Hessenberg estimate. Independently implemented Arnoldi/Givens algebra.
inline NewtonVector eulerGmres(const std::function<NewtonVector(const NewtonVector&)>& apply,
    const NewtonVector& rhs,std::size_t maximum,std::size_t& iterations) {
    const auto n=rhs.size();NewtonVector x(n,0),r=rhs;
    const double target=.03*newtonNorm(rhs);
    const std::size_t restart=std::min<std::size_t>(40,n);
    for(std::size_t used=0;used<maximum;) {
        const double beta=newtonNorm(r);if(beta<=target)return x;
        std::vector<NewtonVector> v(restart+1,NewtonVector(n)),h(restart+1,NewtonVector(restart));
        NewtonVector cs(restart),sn(restart),g(restart+1);g[0]=beta;
        for(std::size_t i=0;i<n;++i)v[0][i]=r[i]/beta;
        std::size_t m=0;
        for(;m<restart&&used<maximum;++m) {
            auto w=apply(v[m]);++used;++iterations;
            for(unsigned pass=0;pass<2;++pass)for(std::size_t j=0;j<=m;++j) {
                const double d=newtonDot(w,v[j]);h[j][m]+=d;
                for(std::size_t i=0;i<n;++i)w[i]-=d*v[j][i];
            }
            h[m+1][m]=newtonNorm(w);
            if(h[m+1][m]>0)for(std::size_t i=0;i<n;++i)v[m+1][i]=w[i]/h[m+1][m];
            for(std::size_t j=0;j<m;++j) {
                const double a=cs[j]*h[j][m]+sn[j]*h[j+1][m];
                h[j+1][m]=-sn[j]*h[j][m]+cs[j]*h[j+1][m];h[j][m]=a;
            }
            const double d=std::hypot(h[m][m],h[m+1][m]);
            if(!(d>0&&std::isfinite(d)))throw std::runtime_error("Euler implicit: singular Arnoldi system");
            cs[m]=h[m][m]/d;sn[m]=h[m+1][m]/d;h[m][m]=d;
            g[m+1]=-sn[m]*g[m];g[m]*=cs[m];
            if(std::abs(g[m+1])<=target){++m;break;}
        }
        NewtonVector y(m);
        for(std::size_t j=m;j-->0;) {
            double value=g[j];for(std::size_t k=j+1;k<m;++k)value-=h[j][k]*y[k];y[j]=value/h[j][j];
        }
        for(std::size_t j=0;j<m;++j)for(std::size_t i=0;i<n;++i)x[i]+=v[j][i]*y[j];
        const auto ax=apply(x);for(std::size_t i=0;i<n;++i)r[i]=rhs[i]-ax[i];
    }
    if(newtonNorm(r)<=target)return x;
    throw std::runtime_error("Euler implicit: GMRES iteration budget exhausted");
}
} // namespace cartmesh2d::fv::detail
