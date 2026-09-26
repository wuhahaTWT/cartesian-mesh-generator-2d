#include "cartmesh2d/fv/Euler2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
EulerConservative2D physicalFlux(const EulerConservative2D& u,const EulerPrimitive2D& q,Vector2D n) {
    const double un=q.u*n.x+q.v*n.y;
    return {u[0]*un,u[1]*un+q.pressure*n.x,u[2]*un+q.pressure*n.y,(u[3]+q.pressure)*un};
}
bool validStar(const EulerConservative2D& u,const IdealGas2D& gas) {
    try {(void)eulerPrimitive2D(u,gas);return true;}catch(const std::runtime_error&){return false;}
}
}
EulerFaceFlux2D eulerFaceFlux2D(const EulerConservative2D& ul,const EulerConservative2D& ur,
    Vector2D areaVector,const IdealGas2D& gas,EulerFluxScheme2D scheme,double contactRestoration) {
    if(!std::isfinite(contactRestoration)||contactRestoration<0||contactRestoration>1)
        throw std::runtime_error("Euler flux: contact restoration must be in [0,1]");
    const auto l=eulerPrimitive2D(ul,gas),r=eulerPrimitive2D(ur,gas);
    const double length=std::hypot(areaVector.x,areaVector.y);
    if(!std::isfinite(length)||length<=0)throw std::runtime_error("Euler flux: invalid face area vector");
    if(scheme!=EulerFluxScheme2D::Rusanov&&scheme!=EulerFluxScheme2D::Hllc)
        throw std::runtime_error("Euler flux: unknown flux scheme");
    const Vector2D n{areaVector.x/length,areaVector.y/length};
    const double vl=l.u*n.x+l.v*n.y,vr=r.u*n.x+r.v*n.y;
    const double al=eulerSoundSpeed2D(l,gas),ar=eulerSoundSpeed2D(r,gas);
    const auto fl=physicalFlux(ul,l,n),fr=physicalFlux(ur,r,n);
    EulerFaceFlux2D result;result.waveSpeed=std::max(std::abs(vl)+al,std::abs(vr)+ar);
    const auto rusanov=[&]{
        for(std::size_t k=0;k<4;++k)
            result.integratedFlux[k]=.5*length*(fl[k]+fr[k]-result.waveSpeed*(ur[k]-ul[k]));
    };
    if(scheme==EulerFluxScheme2D::Rusanov)rusanov();
    else {
        // Einfeldt/Roe estimates enlarged to include BOTH endpoint acoustic
        // cones. No claim that these estimates bound every exact nonlinear wave.
        // Toro's HLLC Rankine-Hugoniot star construction, with an explicit
        // admissibility guard rather than pressure or conserved-energy floors.
        const double rl=std::sqrt(l.density),rr=std::sqrt(r.density),weight=rl+rr;
        const double ux=(rl*l.u+rr*r.u)/weight,uy=(rl*l.v+rr*r.v)/weight;
        const double enthalpy=(rl*(ul[3]+l.pressure)/l.density+rr*(ur[3]+r.pressure)/r.density)/weight;
        const double a2=(gas.gamma-1)*(enthalpy-.5*(ux*ux+uy*uy));
        bool usable=std::isfinite(a2)&&a2>0;
        const double roeSound=usable?std::sqrt(a2):0,roeNormal=ux*n.x+uy*n.y;
        const double sl=std::min({vl-al,vr-ar,roeNormal-roeSound});
        const double sr=std::max({vl+al,vr+ar,roeNormal+roeSound});
        if(usable)result.waveSpeed=std::max({result.waveSpeed,std::abs(sl),std::abs(sr)});
        if(usable&&sl>=0)for(std::size_t k=0;k<4;++k)result.integratedFlux[k]=length*fl[k];
        else if(usable&&sr<=0)for(std::size_t k=0;k<4;++k)result.integratedFlux[k]=length*fr[k];
        else {
            const double dl=l.density*(sl-vl),dr=r.density*(sr-vr);
            const double sm=(r.pressure-l.pressure+dl*vl-dr*vr)/(dl-dr);
            usable=usable&&std::isfinite(sm)&&sl<sm&&sm<sr;
            EulerConservative2D stars[2]{};
            for(unsigned side=0;side<2&&usable;++side) {
                const auto& q=side?r:l;const auto& u=side?ur:ul;
                const double wave=side?sr:sl,vn=side?vr:vl;
                const double ps=q.pressure+q.density*(wave-vn)*(sm-vn);
                const double rho=q.density*(wave-vn)/(wave-sm);
                auto& star=stars[side];
                star={rho,rho*(q.u+(sm-vn)*n.x),rho*(q.v+(sm-vn)*n.y),
                    ((wave-vn)*u[3]-q.pressure*vn+ps*sm)/(wave-sm)};
                usable=std::isfinite(ps)&&ps>0&&validStar(star,gas);
            }
            if(usable) {
                const bool left=sm>=0;const auto& f=left?fl:fr;const auto& u=left?ul:ur;
                const auto& star=stars[left?0:1];const double wave=left?sl:sr;
                for(std::size_t k=0;k<4;++k) {
                    const double hllc=f[k]+wave*(star[k]-u[k]);
                    const double hlle=(sr*fl[k]-sl*fr[k]+sl*sr*(ur[k]-ul[k]))/(sr-sl);
                    result.integratedFlux[k]=length*(contactRestoration==1?hllc:hlle+contactRestoration*(hllc-hlle));
                }
            }else {result.hllcFallback=true;rusanov();}
        }
    }
    if(!std::isfinite(result.waveSpeed)||result.waveSpeed<=0)
        throw std::runtime_error("Euler flux: non-finite acoustic speed");
    for(double value:result.integratedFlux)if(!std::isfinite(value))
        throw std::runtime_error("Euler flux: non-finite numerical flux");
    return result;
}
} // namespace cartmesh2d::fv
