#include "cartmesh2d/fv/Diffusion2D.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
using Values=std::vector<double>;
double checked(double value) {
    if(!std::isfinite(value)) throw std::runtime_error("FVM diffusion: numerical range exceeded");
    return value;
}
double inner(const Values& a,const Values& b) {
    long double sum=0; for(std::size_t i=0;i<a.size();++i) sum+=static_cast<long double>(a[i])*b[i];
    return checked(static_cast<double>(sum));
}
double norm(const Values& a) {
    long double result=0;
    for(double v:a) result=std::hypot(result,static_cast<long double>(checked(v)));
    return checked(static_cast<double>(result));
}
struct Matrix {
    const FvMesh2D& mesh;
    double k;
    Values diagonal;
    Values apply(const Values& x) const {
        Values y(x.size());
        // Difference form preserves constants on internal faces to roundoff.
        for (const auto& f:mesh.faces) {
            const double a=k*f.transmissibility;
            if(f.neighbour) { const auto j=*f.neighbour; const double q=a*(x[f.owner]-x[j]); y[f.owner]+=q; y[j]-=q; }
            else y[f.owner]+=a*x[f.owner];
        }
        return y;
    }
};
std::size_t pcg(const Matrix& a,const Values& b,Values& x,const DiffusionControls2D& c) {
    auto ax=a.apply(x); Values r(b.size()),z(b.size()),p(b.size());
    const double stop=checked(c.linearAbsoluteTolerance+c.linearRelativeTolerance*norm(b));
    for(std::size_t i=0;i<b.size();++i) { r[i]=b[i]-ax[i]; z[i]=r[i]/a.diagonal[i]; }
    if(norm(r)<=stop) return 0;
    p=z; double rz=inner(r,z);
    for(std::size_t it=1;it<=c.maxLinearIterations;++it) {
        const auto ap=a.apply(p); const double pap=inner(p,ap);
        if (!(pap>0) || !std::isfinite(pap) || !(rz>0)) throw std::runtime_error("FVM PCG: matrix breakdown");
        const double alpha=rz/pap;
        for(std::size_t i=0;i<x.size();++i) { x[i]+=alpha*p[i]; r[i]-=alpha*ap[i]; }
        if(norm(r)<=stop || it%50==0) {
            // Check the true algebraic residual, not just the recurrence.
            ax=a.apply(x);
            for(std::size_t i=0;i<x.size();++i) r[i]=b[i]-ax[i];
            if(norm(r)<=stop) return it;
            for(std::size_t i=0;i<x.size();++i) z[i]=r[i]/a.diagonal[i];
            p=z; rz=inner(r,z); continue;
        }
        for(std::size_t i=0;i<x.size();++i) z[i]=r[i]/a.diagonal[i];
        const double next=inner(r,z),beta=next/rz;
        for(std::size_t i=0;i<x.size();++i) p[i]=z[i]+beta*p[i];
        rz=next;
    }
    throw std::runtime_error("FVM PCG: iteration limit reached");
}
Values correctionFlux(const FvMesh2D& mesh,const Values& u,const Values& bc,double k) {
    const auto gradients=reconstructGradient(mesh,u,bc); Values correction(mesh.faces.size());
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id]; auto g=gradients[f.owner];
        if(f.neighbour) { const auto gn=gradients[*f.neighbour]; const double w=f.neighbourWeight;
            g={g.x*(1-w)+gn.x*w,g.y*(1-w)+gn.y*w}; }
        correction[id]=-k*dot(g,f.correction);
    }
    return correction;
}
}
DiffusionResult2D solveDiffusion2D(const FvMesh2D& mesh,const DiffusionProblem2D& problem,const DiffusionControls2D& c) {
    validateFvMesh2D(mesh);
    if(mesh.cells.empty() || !problem.source || !problem.boundaryValue ||
       !std::isfinite(problem.diffusivity) || problem.diffusivity<=0 ||
       !std::isfinite(c.relaxation) || c.relaxation<=0 || c.relaxation>1 ||
       c.maxCorrections==0 || c.maxLinearIterations==0)
        throw std::invalid_argument("FVM diffusion: invalid problem or controls");
    for(double t:{c.relativeTolerance,c.absoluteTolerance,c.linearRelativeTolerance,c.linearAbsoluteTolerance})
        if(!std::isfinite(t) || t<=0) throw std::invalid_argument("FVM diffusion: tolerances must be finite positive");
    DiffusionResult2D result;
    const auto n=mesh.cells.size(); const double k=problem.diffusivity;
    result.values.assign(n,0); result.sourceIntegrals.resize(n);
    Matrix a{mesh,k,Values(n)}; Values base(n),bc(mesh.faces.size());
    std::vector<bool> anchored(n,false); std::vector<std::size_t> queue;
    for(std::size_t i=0;i<n;++i) {
        base[i]=problem.source(mesh.cells[i].centre)*mesh.cells[i].area;
        if(!std::isfinite(base[i])) throw std::runtime_error("FVM diffusion: nonfinite source");
        result.sourceIntegrals[i]=base[i]; result.sourceIntegral=checked(result.sourceIntegral+base[i]);
    }
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id]; const double t=k*f.transmissibility;
        a.diagonal[f.owner]+=t;
        if(f.neighbour) a.diagonal[*f.neighbour]+=t;
        else {
            bc[id]=problem.boundaryValue(f.centre,f.patch);
            if(!std::isfinite(bc[id])) throw std::runtime_error("FVM diffusion: nonfinite boundary value");
            base[f.owner]+=t*bc[id];
            if(!anchored[f.owner]) { anchored[f.owner]=true; queue.push_back(f.owner); }
        }
    }
    for(std::size_t head=0;head<queue.size();++head) {
        const auto i=queue[head];
        for(auto id:mesh.cells[i].faces) { const auto& f=mesh.faces[id]; if(!f.neighbour) continue;
            const auto j=f.owner==i?*f.neighbour:f.owner;
            if(!anchored[j]) { anchored[j]=true; queue.push_back(j); }
        }
    }
    if(queue.size()!=n) throw std::runtime_error("FVM diffusion: unanchored connected component");
    for(double d:a.diagonal) if(!(d>0) || !std::isfinite(d)) throw std::runtime_error("FVM diffusion: invalid diagonal");
    const double scale=norm(base),stop=checked(c.absoluteTolerance+c.relativeTolerance*scale);
    for(std::size_t it=1;it<=c.maxCorrections;++it) {
        const auto correction=correctionFlux(mesh,result.values,bc,k); auto rhs=base;
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id]; rhs[f.owner]-=correction[id];
            if(f.neighbour) rhs[*f.neighbour]+=correction[id];
        }
        auto candidate=result.values; const auto linear=pcg(a,rhs,candidate,c);
        for(std::size_t i=0;i<n;++i) {
            result.values[i]+=c.relaxation*(candidate[i]-result.values[i]);
            if(!std::isfinite(result.values[i])) throw std::runtime_error("FVM diffusion: nonfinite solution");
        }
        const auto finalCorrection=correctionFlux(mesh,result.values,bc,k);
        Values residual(n); for(std::size_t i=0;i<n;++i) residual[i]=-result.sourceIntegrals[i];
        result.fluxes.resize(mesh.faces.size()); result.boundaryFlux=0;
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];
            const double other=f.neighbour?result.values[*f.neighbour]:bc[id];
            const double flux=k*f.transmissibility*(result.values[f.owner]-other)+finalCorrection[id];
            if(!std::isfinite(flux)) throw std::runtime_error("FVM diffusion: nonfinite face flux");
            result.fluxes[id]=flux; residual[f.owner]+=flux;
            if(f.neighbour) residual[*f.neighbour]-=flux; else result.boundaryFlux+=flux;
        }
        double maximum=0; for(double r:residual) maximum=std::max(maximum,std::abs(r));
        const double residualNorm=norm(residual);
        result.history.push_back({it,linear,checked(residualNorm/std::max(scale,c.absoluteTolerance)),residualNorm,maximum});
        result.boundaryFlux=checked(result.boundaryFlux);
        result.globalBalance=checked(result.boundaryFlux-result.sourceIntegral);
        if(residualNorm<=stop) { result.converged=true; break; }
    }
    return result;
}
}
