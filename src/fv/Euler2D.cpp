#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/HybridHeat2D.hpp"
#include "cartmesh2d/fv/HybridViscous2D.hpp"
#include "cartmesh2d/fv/detail/EulerNewton2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace cartmesh2d::fv {
struct EulerDiffusionOperators2D {
    std::optional<HeatConductionOperator2D> heat;
    std::optional<ViscousStressOperator2D> viscous;
    std::optional<HybridHeatOperator2D> hybridHeat;
    std::optional<HybridViscousOperator2D> hybridViscous;
    std::vector<double> heatCellStiffness;
    [[nodiscard]] bool active() const{return heat||viscous||hybridHeat||hybridViscous;}
};
namespace {
void require(bool valid,const std::string& message) { if(!valid)throw std::runtime_error("Euler: "+message); }
double finite(double x) {require(std::isfinite(x),"non-finite arithmetic");return x;}
void primitiveValid(const EulerPrimitive2D& q) {
    require(std::isfinite(q.density)&&q.density>0&&std::isfinite(q.pressure)&&q.pressure>0&&
            std::isfinite(q.u)&&std::isfinite(q.v),"density and pressure must be strictly positive; velocity must be finite");
}
EulerPrimitive2D ghost(const EulerPrimitive2D& inside,const EulerBoundary2D& b,Vector2D n,const IdealGas2D& gas) {
    const double normal=inside.u*n.x+inside.v*n.y;
    if((b.kind==EulerBoundaryKind2D::SlipWall||b.kind==EulerBoundaryKind2D::NoSlipWall))
        return {inside.density,inside.u-2*normal*n.x,inside.v-2*normal*n.y,inside.pressure};
    if(b.kind==EulerBoundaryKind2D::Transmissive)return inside;
    if(b.kind==EulerBoundaryKind2D::PressureOutlet)return eulerPressureOutletState2D(inside,b.reference.pressure,n,gas);
    if(b.kind==EulerBoundaryKind2D::TotalInlet)return eulerTotalInletState2D(inside,b.reference,n,gas);
    require(b.kind==EulerBoundaryKind2D::Farfield,"unexpected boundary ghost type");
    const double sound=eulerSoundSpeed2D(inside,gas);
    if(normal>=sound)return inside;
    if(normal<=-sound)return b.reference;
    const auto& far=b.reference;
    const double farNormal=far.u*n.x+far.v*n.y,farSound=eulerSoundSpeed2D(far,gas);
    const double outgoing=normal+2*sound/(gas.gamma-1),incoming=farNormal-2*farSound/(gas.gamma-1);
    const double boundaryNormal=.5*(outgoing+incoming),boundarySound=.25*(gas.gamma-1)*(outgoing-incoming);
    require(std::isfinite(boundarySound)&&boundarySound>0,"farfield characteristics imply non-positive sound speed");
    const auto& entropyState=boundaryNormal<0 ? far : inside;
    const double entropy=entropyState.pressure/std::pow(entropyState.density,gas.gamma);
    const double density=finite(std::pow(boundarySound*boundarySound/(gas.gamma*entropy),1/(gas.gamma-1)));
    const double oldNormal=entropyState.u*n.x+entropyState.v*n.y;
    EulerPrimitive2D result{density,entropyState.u+(boundaryNormal-oldNormal)*n.x,
        entropyState.v+(boundaryNormal-oldNormal)*n.y,finite(entropy*std::pow(density,gas.gamma))};
    primitiveValid(result);return result;
}
}

void validateEulerTransport2D(const EulerTransport2D& transport,const std::vector<EulerBoundary2D>& boundaries) {
    require(std::isfinite(transport.dynamicViscosity)&&transport.dynamicViscosity>=0,"dynamic viscosity must be finite and nonnegative");
    require(std::isfinite(transport.thermalConductivity)&&transport.thermalConductivity>=0,"thermal conductivity must be finite and nonnegative");
    for(const auto& b:boundaries) {
        require(b.kind!=EulerBoundaryKind2D::NoSlipWall||transport.dynamicViscosity>0,"no-slip wall requires positive dynamic viscosity");
        require(b.thermalKind==HeatBoundaryKind2D::Insulated||b.thermalKind==HeatBoundaryKind2D::Temperature||b.thermalKind==HeatBoundaryKind2D::OutwardFlux,"unknown thermal condition");
        require(std::isfinite(b.thermalValue),"nonfinite thermal boundary value");
        require(b.thermalKind!=HeatBoundaryKind2D::Temperature||b.thermalValue>0,"temperature boundary requires positive Kelvin");
        require(b.thermalKind!=HeatBoundaryKind2D::Insulated||b.thermalValue==0,"insulated boundary must have zero thermal value");
        require(b.kind!=EulerBoundaryKind2D::Periodic||b.thermalKind==HeatBoundaryKind2D::Insulated,"periodic thermal coupling is automatic, not a wall condition");
        require(transport.thermalConductivity>0||b.thermalKind==HeatBoundaryKind2D::Insulated,"active thermal boundary requires positive conductivity");
    }
}
void validateIdealGas2D(const IdealGas2D& gas) {
    require(std::isfinite(gas.gamma)&&gas.gamma>1&&std::isfinite(gas.gasConstant)&&gas.gasConstant>0,
        "ideal gas requires gamma > 1 and R > 0");
}
EulerConservative2D eulerConservative2D(const EulerPrimitive2D& q,const IdealGas2D& gas) {
    validateIdealGas2D(gas);primitiveValid(q);
    const double energy=finite(q.pressure/(gas.gamma-1)+.5*q.density*(q.u*q.u+q.v*q.v));
    return {q.density,finite(q.density*q.u),finite(q.density*q.v),energy};
}
EulerPrimitive2D eulerPrimitive2D(const EulerConservative2D& q,const IdealGas2D& gas) {
    validateIdealGas2D(gas);for(double value:q)finite(value);
    require(q[0]>0,"non-positive conserved density");
    const double u=finite(q[1]/q[0]),v=finite(q[2]/q[0]);
    // Subtract kinetic energy at extended precision, but never manufacture
    // positive internal energy when it is not represented in the input.
    const long double kinetic=.5L*(static_cast<long double>(q[1])*q[1]+static_cast<long double>(q[2])*q[2])/q[0];
    const double pressure=finite(static_cast<double>((gas.gamma-1)*(static_cast<long double>(q[3])-kinetic)));
    EulerPrimitive2D result{q[0],u,v,pressure};primitiveValid(result);return result;
}
double eulerSoundSpeed2D(const EulerPrimitive2D& q,const IdealGas2D& gas) {
    validateIdealGas2D(gas);primitiveValid(q);return finite(std::sqrt(gas.gamma*q.pressure/q.density));
}
EulerPrimitive2D eulerPressureOutletState2D(const EulerPrimitive2D& inside,double pressure,Vector2D area,const IdealGas2D& gas) {
    const double sound=eulerSoundSpeed2D(inside,gas),length=std::hypot(area.x,area.y);
    require(std::isfinite(length)&&length>0,"pressure-outlet requires a finite nonzero outward area");
    require(std::isfinite(pressure)&&pressure>0,"pressure-outlet requires positive absolute pressure");
    const Vector2D n{area.x/length,area.y/length};
    const double normal=inside.u*n.x+inside.v*n.y;
    require(normal>=0,"pressure-outlet backflow is outside the supported unidirectional regime");
    if(normal>=sound)return inside; // No incoming characteristic can impose p.
    const double ratio=pressure/inside.pressure;
    const double density=finite(inside.density*std::pow(ratio,1/gas.gamma));
    const double boundarySound=finite(sound*std::pow(ratio,(gas.gamma-1)/(2*gas.gamma)));
    // J+ = un + 2c/(gamma-1), entropy p/rho^gamma, and ut leave the domain.
    const double boundaryNormal=finite(normal+2*(sound-boundarySound)/(gas.gamma-1));
    require(boundaryNormal>=0,"pressure-outlet back pressure induces unsupported backflow");
    require(boundaryNormal<boundarySound,"pressure-outlet back pressure induces choking; extend the domain or change the condition");
    EulerPrimitive2D result{density,inside.u+(boundaryNormal-normal)*n.x,
        inside.v+(boundaryNormal-normal)*n.y,pressure};
    primitiveValid(result);return result;
}
EulerPrimitive2D eulerTotalInletState2D(const EulerPrimitive2D& inside,const EulerPrimitive2D& ref,Vector2D area,const IdealGas2D& gas) {
    primitiveValid(ref);const double length=std::hypot(area.x,area.y);
    require(std::isfinite(length)&&length>0,"total inlet requires finite nonzero outward area");
    const Vector2D n{area.x/length,area.y/length};const double a=gas.gamma-1;
    const double un=dot(Vector2D{inside.u,inside.v},n),c=eulerSoundSpeed2D(inside,gas);
    require(un<c,"total inlet cannot constrain supersonic outflow");
    const double refn=dot(Vector2D{ref.u,ref.v},n),utx=ref.u-refn*n.x,uty=ref.v-refn*n.y;
    const double h0=gas.gamma/a*ref.pressure/ref.density+.5*(ref.u*ref.u+ref.v*ref.v);
    const double available=h0-.5*(utx*utx+uty*uty),outgoing=un+2*c/a;
    require(available>0,"total inlet has invalid total enthalpy");
    double lo=-std::sqrt(2*a*available/(gas.gamma+1)),hi=0;
    const auto sound=[&](double u){return std::sqrt(a*(available-.5*u*u));};
    const auto invariant=[&](double u){return u+2*sound(u)/a;};
    require(outgoing>invariant(lo)&&outgoing<invariant(hi),"total inlet implies choking or flow reversal; change reservoir conditions");
    for(unsigned i=0;i<80;++i){const double mid=.5*(lo+hi);if(invariant(mid)<outgoing)lo=mid;else hi=mid;}
    const double normal=.5*(lo+hi),cb=sound(normal),entropy=ref.pressure/std::pow(ref.density,gas.gamma);
    const double rho=std::pow(cb*cb/(gas.gamma*entropy),1/a);
    EulerPrimitive2D q{rho,utx+normal*n.x,uty+normal*n.y,entropy*std::pow(rho,gas.gamma)};primitiveValid(q);return q;
}
void validateEulerBoundaries2D(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,const IdealGas2D& gas) {
    validateIdealGas2D(gas);std::vector<const EulerBoundary2D*> lookup(mesh.faces.size(),nullptr);
    std::map<std::string,EulerBoundaryKind2D> groups;
    for(const auto& b:boundaries) {
        require(b.face<mesh.faces.size()&&!mesh.faces[b.face].neighbour&&!lookup[b.face],"boundary face missing, duplicated or internal");
        require(!b.name.empty(),"boundary name is empty");
        const auto [it,inserted]=groups.emplace(b.name,b.kind);
        require(inserted||it->second==b.kind,"one boundary name has multiple physical kinds");
        require(b.kind==EulerBoundaryKind2D::SlipWall||b.kind==EulerBoundaryKind2D::Transmissive||
                b.kind==EulerBoundaryKind2D::Farfield||b.kind==EulerBoundaryKind2D::Periodic||b.kind==EulerBoundaryKind2D::NoSlipWall||
                b.kind==EulerBoundaryKind2D::PressureOutlet||b.kind==EulerBoundaryKind2D::TotalInlet,"unknown boundary kind");
        require(std::isfinite(b.wallVelocity.x)&&std::isfinite(b.wallVelocity.y),"nonfinite wall velocity");
        require(b.kind==EulerBoundaryKind2D::NoSlipWall||(b.wallVelocity.x==0&&b.wallVelocity.y==0),"only no-slip walls accept wall velocity");
        const auto normal=mesh.faces[b.face].areaVector;
        require(std::abs(dot(b.wallVelocity,normal))<=64*std::numeric_limits<double>::epsilon()*std::hypot(normal.x,normal.y)*std::hypot(b.wallVelocity.x,b.wallVelocity.y),"static mesh requires tangential wall velocity");
        primitiveValid(b.reference);
        require((b.kind==EulerBoundaryKind2D::Periodic)==b.partner.has_value(),"only periodic faces require a partner");
        lookup[b.face]=&b;
    }
    std::map<std::string,Vector2D> translations;
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& face=mesh.faces[id];if(face.neighbour)continue;
        require(lookup[id]!=nullptr,"boundary coverage is incomplete");
        const auto& b=*lookup[id];if(!b.partner)continue;
        const auto partner=*b.partner;
        require(partner<mesh.faces.size()&&partner!=id&&lookup[partner]&&lookup[partner]->partner==id&&
                lookup[partner]->kind==EulerBoundaryKind2D::Periodic&&lookup[partner]->name==b.name,"periodic pairing is not reciprocal");
        const auto& other=mesh.faces[partner];const auto normal=face.areaVector,opposite=other.areaVector;
        const double length=std::hypot(normal.x,normal.y);
        require(std::hypot(normal.x+opposite.x,normal.y+opposite.y)<=1e-10*length,"periodic faces have unequal or non-opposite normals");
        if(id>partner)continue;
        Vector2D translation{other.centre.x-face.centre.x,other.centre.y-face.centre.y};
        const double distance=std::hypot(translation.x,translation.y);
        require(distance>0,"periodic translation is zero");
        if(translation.x<0 || (translation.x==0&&translation.y<0)){translation.x=-translation.x;translation.y=-translation.y;}
        const auto [it,inserted]=translations.emplace(b.name,translation);
        require(inserted||std::hypot(translation.x-it->second.x,translation.y-it->second.y)<=1e-10*distance,
                "one periodic group has inconsistent translations");
    }
}

namespace {
using PrimitiveValues=std::array<double,4>;
using Gradient=std::array<Vector2D,4>;
PrimitiveValues values(const EulerPrimitive2D& q){return {q.density,q.u,q.v,q.pressure};}
EulerPrimitive2D primitiveValues(const PrimitiveValues& q){return {q[0],q[1],q[2],q[3]};}
struct SpatialOperator {
    std::vector<EulerConservative2D> faceFlux,residual;
    std::vector<double> speed,spectral,heatFlux,heatRate,viscousRate;
    std::vector<std::array<double,3>> viscousFlux;
    std::vector<unsigned char> fallback;
    std::size_t fallbackEvaluations=0,reconstructionFallbackCells=0;
    double minimumContactRestoration=1;
};
SpatialOperator spatialOperator(const FvMesh2D& mesh,const std::vector<const EulerBoundary2D*>& lookup,
    const IdealGas2D& gas,const std::vector<EulerConservative2D>& cells,const EulerStepControls2D& control,const EulerDiffusionOperators2D& diffusion) {
    const auto* heat=diffusion.heat?&*diffusion.heat:nullptr;const auto* viscous=diffusion.viscous?&*diffusion.viscous:nullptr;
    const auto nc=mesh.cells.size(),nf=mesh.faces.size();
    std::vector<EulerPrimitive2D> primitive;primitive.reserve(nc);
    for(const auto& u:cells)primitive.push_back(eulerPrimitive2D(u,gas));
    SpatialOperator out;out.faceFlux.resize(nf);out.residual.resize(nc);out.speed.resize(nf);
    out.spectral.resize(nc);out.fallback.resize(nf);out.heatFlux.resize(nf);out.heatRate.resize(nc);out.viscousFlux.resize(nf);out.viscousRate.resize(nc);
    // Multidimensional pressure sensor inspired by Simon & Mandal (2018),
    // eqs. 49-50, alpha=3. Our polygon stencil includes every face incident on
    // either adjacent cell and blends ALL conserved flux components, rather
    // than claiming to reproduce the paper's structured, selective ADC scheme.
    std::vector<double> contactWeight(nc,1.);
    if(control.fluxScheme==EulerFluxScheme2D::Hllc)for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];const auto* bc=lookup[id];auto neighbour=face.neighbour;
        if(bc&&bc->partner)neighbour=mesh.faces[*bc->partner].owner;
        const double length=std::hypot(face.areaVector.x,face.areaVector.y);
        const double pl=primitive[face.owner].pressure,pr=neighbour?primitive[*neighbour].pressure:
            ghost(primitive[face.owner],*bc,{face.areaVector.x/length,face.areaVector.y/length},gas).pressure;
        const double ratio=std::min(pl,pr)/std::max(pl,pr),weight=ratio*ratio*ratio;
        contactWeight[face.owner]=std::min(contactWeight[face.owner],weight);
        if(neighbour)contactWeight[*neighbour]=std::min(contactWeight[*neighbour],weight);
    }
    std::vector<Gradient> gradient(control.order==2?nc:0);
    if(control.order==2)for(std::size_t i=0;i<nc;++i) {
        const auto& cell=mesh.cells[i];
        // Limit velocity components in a frame attached to this cell's first
        // edge, not to the global x/y axes. A rigid mesh rotation rotates the
        // limiter frame as well, avoiding coordinate-dependent velocity limiting.
        const auto axis=mesh.faces[cell.faces.front()].areaVector;
        const double axisLength=std::hypot(axis.x,axis.y),ex=axis.x/axisLength,ey=axis.y/axisLength;
        const auto localValues=[&](const EulerPrimitive2D& q){return PrimitiveValues{q.density,q.u*ex+q.v*ey,-q.u*ey+q.v*ex,q.pressure};};
        const auto centre=localValues(primitive[i]);
        auto minimum=centre,maximum=centre;double xx=0,xy=0,yy=0;
        std::array<double,4> bx{},by{};
        for(const auto id:cell.faces) {
            const auto& f=mesh.faces[id];const auto* bc=lookup[id];
            Vector2D d{};EulerPrimitive2D other{};
            if(f.neighbour) {
                const auto j=i==f.owner?*f.neighbour:f.owner;
                d={mesh.cells[j].centre.x-cell.centre.x,mesh.cells[j].centre.y-cell.centre.y};other=primitive[j];
            }else if(bc->partner) {
                const auto& partner=mesh.faces[*bc->partner];const auto& c=mesh.cells[partner.owner];
                d={c.centre.x+f.centre.x-partner.centre.x-cell.centre.x,
                    c.centre.y+f.centre.y-partner.centre.y-cell.centre.y};other=primitive[partner.owner];
            }else {
                const double length=std::hypot(f.areaVector.x,f.areaVector.y);
                const Vector2D normal{f.areaVector.x/length,f.areaVector.y/length};
                const double distance=(f.centre.x-cell.centre.x)*normal.x+(f.centre.y-cell.centre.y)*normal.y;
                d={2*distance*normal.x,2*distance*normal.y};other=ghost(primitive[i],*bc,normal,gas);
                // The pressure outlet supplies a trace at the actual face,
                // unlike the mirrored state used by solid walls.
                if(bc->kind==EulerBoundaryKind2D::PressureOutlet||bc->kind==EulerBoundaryKind2D::TotalInlet)d={f.centre.x-cell.centre.x,f.centre.y-cell.centre.y};
            }
            const double distance2=finite(d.x*d.x+d.y*d.y);require(distance2>0,"zero reconstruction stencil distance");
            const double w=1/distance2;xx+=w*d.x*d.x;xy+=w*d.x*d.y;yy+=w*d.y*d.y;
            const auto neighbour=localValues(other);
            for(std::size_t k=0;k<4;++k) {
                minimum[k]=std::min(minimum[k],neighbour[k]);maximum[k]=std::max(maximum[k],neighbour[k]);
                bx[k]+=w*d.x*(neighbour[k]-centre[k]);by[k]+=w*d.y*(neighbour[k]-centre[k]);
            }
        }
        const double determinant=finite(xx*yy-xy*xy),trace=finite(xx+yy);
        // Relative conditioning diagnostic in dimensionless least-squares geometry.
        // An unresolved gradient becomes an explicitly counted constant reconstruction.
        if(determinant<=64*std::numeric_limits<double>::epsilon()*trace*trace) {
            ++out.reconstructionFallbackCells;continue;
        }
        auto& grad=gradient[i];PrimitiveValues theta{1,1,1,1};
        for(std::size_t k=0;k<4;++k)grad[k]={finite((yy*bx[k]-xy*by[k])/determinant),finite((xx*by[k]-xy*bx[k])/determinant)};
        // Barth-Jespersen face limiter: reconstruct primitive variables inside
        // the cell/neighbour extrema, so positive density/pressure remain positive.
        for(const auto id:cell.faces) {
            const auto c=mesh.faces[id].centre;const Vector2D d{c.x-cell.centre.x,c.y-cell.centre.y};
            for(std::size_t k=0;k<4;++k) {
                const double increment=finite(grad[k].x*d.x+grad[k].y*d.y);
                if(increment>0)theta[k]=std::min(theta[k],(maximum[k]-centre[k])/increment);
                else if(increment<0)theta[k]=std::min(theta[k],(minimum[k]-centre[k])/increment);
            }
        }
        for(std::size_t k=0;k<4;++k){theta[k]=std::clamp(theta[k],0.,1.);grad[k].x*=theta[k];grad[k].y*=theta[k];}
        bool positive=true;
        for(const auto id:cell.faces) {
            auto q=centre;const auto c=mesh.faces[id].centre;
            for(std::size_t k=0;k<4;++k)q[k]+=grad[k].x*(c.x-cell.centre.x)+grad[k].y*(c.y-cell.centre.y);
            try {primitiveValid(primitiveValues(q));}catch(const std::runtime_error&){positive=false;break;}
        }
        // Floating-point cancellation at extreme contrasts must not become a floor.
        if(!positive){grad={};++out.reconstructionFallbackCells;}
        const auto normalGradient=grad[1],tangentGradient=grad[2];
        grad[1]={ex*normalGradient.x-ey*tangentGradient.x,ex*normalGradient.y-ey*tangentGradient.y};
        grad[2]={ey*normalGradient.x+ex*tangentGradient.x,ey*normalGradient.y+ex*tangentGradient.y};
    }
    const auto reconstructed=[&](std::size_t i,Point2D faceCentre) {
        if(control.order==1)return cells[i];
        auto q=values(primitive[i]);const auto c=mesh.cells[i].centre;
        for(std::size_t k=0;k<4;++k)q[k]+=gradient[i][k].x*(faceCentre.x-c.x)+gradient[i][k].y*(faceCentre.y-c.y);
        return eulerConservative2D(primitiveValues(q),gas);
    };
    for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];const auto* boundary=lookup[id];
        if(boundary&&boundary->partner&&id>*boundary->partner)continue;
        const auto owner=face.owner;auto neighbour=face.neighbour;auto rightCentre=face.centre;
        if(boundary&&boundary->partner) {
            const auto& other=mesh.faces[*boundary->partner];neighbour=other.owner;rightCentre=other.centre;
        }
        const auto left=reconstructed(owner,face.centre);
        const double length=std::hypot(face.areaVector.x,face.areaVector.y);
        const Vector2D normal{face.areaVector.x/length,face.areaVector.y/length};
        const auto right=neighbour?reconstructed(*neighbour,rightCentre):
            eulerConservative2D(ghost(eulerPrimitive2D(left,gas),*boundary,normal,gas),gas);
        const double restoration=neighbour?std::min(contactWeight[owner],contactWeight[*neighbour]):contactWeight[owner];
        auto flux=eulerFaceFlux2D(left,right,face.areaVector,gas,control.fluxScheme,restoration);
        if(boundary&&(boundary->kind==EulerBoundaryKind2D::PressureOutlet||boundary->kind==EulerBoundaryKind2D::TotalInlet)) {
            // Apply the characteristic trace directly. A second Riemann solve
            // against the interior would weaken the specified static pressure.
            const double speed=flux.waveSpeed;
            flux=eulerFaceFlux2D(right,right,face.areaVector,gas,control.fluxScheme);
            flux.waveSpeed=std::max(speed,flux.waveSpeed);
        }
        if(boundary&&(boundary->kind==EulerBoundaryKind2D::SlipWall||boundary->kind==EulerBoundaryKind2D::NoSlipWall)) {
            // The mirror Riemann problem has exactly zero mass/energy flux and
            // purely normal pressure traction. Enforce that analytical symmetry
            // before assembling the residual, rather than leaving cancellation
            // of large SI energy terms to floating-point star-state arithmetic.
            const double traction=flux.integratedFlux[1]*normal.x+flux.integratedFlux[2]*normal.y;
            flux.integratedFlux={0,traction*normal.x,traction*normal.y,0};
        }
        out.minimumContactRestoration=std::min(out.minimumContactRestoration,restoration);
        out.faceFlux[id]=flux.integratedFlux;out.speed[id]=flux.waveSpeed;
        out.fallback[id]=static_cast<unsigned char>(flux.hllcFallback);
        if(flux.hllcFallback)++out.fallbackEvaluations;
        out.spectral[owner]=finite(out.spectral[owner]+flux.waveSpeed*length);
        if(neighbour)out.spectral[*neighbour]=finite(out.spectral[*neighbour]+flux.waveSpeed*length);
        for(std::size_t k=0;k<4;++k) {
            const double value=flux.integratedFlux[k];out.residual[owner][k]=finite(out.residual[owner][k]+value);
            if(neighbour)out.residual[*neighbour][k]=finite(out.residual[*neighbour][k]-value);
            if(boundary&&boundary->partner)out.faceFlux[*boundary->partner][k]=-value;
        }
        if(boundary&&boundary->partner) {
            out.speed[*boundary->partner]=flux.waveSpeed;out.fallback[*boundary->partner]=out.fallback[id];
        }
    }
    if(heat) {
        std::vector<double> temperatures(nc),capacity(nc);
        for(std::size_t i=0;i<nc;++i) {
            temperatures[i]=primitive[i].pressure/(primitive[i].density*gas.gasConstant);
            capacity[i]=primitive[i].density*gas.gasConstant/(gas.gamma-1);
        }
        auto conduction=heat->evaluate(temperatures,capacity);
        out.heatFlux=std::move(conduction.faceHeatFlux);out.heatRate=std::move(conduction.rate);
        for(std::size_t id=0;id<nf;++id)out.faceFlux[id][3]=finite(out.faceFlux[id][3]+out.heatFlux[id]);
        for(std::size_t i=0;i<nc;++i)out.residual[i][3]=finite(out.residual[i][3]+conduction.cellResidual[i]);
    }
    if(viscous) {
        std::vector<Vector2D> velocity(nc);std::vector<double> density(nc);
        for(std::size_t i=0;i<nc;++i){velocity[i]={primitive[i].u,primitive[i].v};density[i]=primitive[i].density;}
        auto stress=viscous->evaluate(velocity,density);
        out.viscousFlux=std::move(stress.faceFlux);out.viscousRate=std::move(stress.rate);
        for(std::size_t id=0;id<nf;++id)for(std::size_t k=0;k<3;++k)out.faceFlux[id][k+1]=finite(out.faceFlux[id][k+1]+out.viscousFlux[id][k]);
        for(std::size_t i=0;i<nc;++i)for(std::size_t k=0;k<3;++k)out.residual[i][k+1]=finite(out.residual[i][k+1]+stress.cellResidual[i][k]);
    }
    if(diffusion.hybridHeat) {
        std::vector<double> temperatures(nc);
        for(std::size_t i=0;i<nc;++i)temperatures[i]=primitive[i].pressure/(primitive[i].density*gas.gasConstant);
        const auto conduction=diffusion.hybridHeat->evaluateAtCells(temperatures,2e-14);
        out.heatFlux=conduction.faceHeatFlux;
        for(std::size_t id=0;id<nf;++id)out.faceFlux[id][3]=finite(out.faceFlux[id][3]+out.heatFlux[id]);
        for(std::size_t i=0;i<nc;++i) {
            out.residual[i][3]=finite(out.residual[i][3]+conduction.cellResidual[i]);
            const double capacity=primitive[i].density*gas.gasConstant/(gas.gamma-1);
            out.heatRate[i]=finite(diffusion.heatCellStiffness[i]/(mesh.cells[i].area*capacity));
        }
    }
    if(diffusion.hybridViscous) {
        std::vector<Vector2D> velocity(nc);std::vector<double> density(nc);
        for(std::size_t i=0;i<nc;++i){velocity[i]={primitive[i].u,primitive[i].v};density[i]=primitive[i].density;}
        auto stress=diffusion.hybridViscous->evaluate(velocity,density,2e-14);
        out.viscousFlux=std::move(stress.faceFlux);out.viscousRate=std::move(stress.rate);
        for(std::size_t id=0;id<nf;++id)for(std::size_t k=0;k<3;++k)out.faceFlux[id][k+1]=finite(out.faceFlux[id][k+1]+out.viscousFlux[id][k]);
        for(std::size_t i=0;i<nc;++i)for(std::size_t k=0;k<3;++k)out.residual[i][k+1]=finite(out.residual[i][k+1]+stress.cellResidual[i][k]);
    }
    return out;
}
bool positiveUpdate(const FvMesh2D& mesh,const std::vector<EulerConservative2D>& old,
    const std::vector<EulerConservative2D>& residual,double dt,const IdealGas2D& gas,
    std::vector<EulerConservative2D>& next) {
    next=old;
    for(std::size_t i=0;i<old.size();++i) {
        for(std::size_t k=0;k<4;++k)next[i][k]=old[i][k]-dt/mesh.cells[i].area*residual[i][k];
        try {(void)eulerPrimitive2D(next[i],gas);}catch(const std::runtime_error&){return false;}
    }
    return true;
}
}
namespace {
struct EulerInterrupted : std::exception {
    const char* what() const noexcept override {return "Euler calculation interrupted by cancellation or wall-time budget; previous accepted state retained";}
};
struct ImplicitStage {std::vector<EulerConservative2D> cells;SpatialOperator op;};
ImplicitStage implicitStage(const FvMesh2D& mesh,const std::vector<const EulerBoundary2D*>& lookup,
    const IdealGas2D& gas,const std::vector<EulerConservative2D>& base,
    const std::vector<EulerConservative2D>& guess,double h,const EulerStepControls2D& control,
    const EulerDiffusionOperators2D& diffusion,EulerStepResult2D& work) {
    using namespace detail;const auto n=base.size();NewtonVector scale(4*n),x(4*n);
    for(std::size_t i=0;i<n;++i) {
        const auto p=eulerPrimitive2D(guess[i],gas);const double c=eulerSoundSpeed2D(p,gas);
        scale[4*i]=p.density;scale[4*i+1]=scale[4*i+2]=p.density*c;scale[4*i+3]=guess[i][3];
        for(std::size_t k=0;k<4;++k)x[4*i+k]=guess[i][k]/scale[4*i+k];
    }
    auto evaluate=[&](const NewtonVector& y) {
        if(control.interrupted&&control.interrupted())throw EulerInterrupted{};
        ImplicitStage s;s.cells.resize(n);
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<4;++k)s.cells[i][k]=y[4*i+k]*scale[4*i+k];
        ++work.spatialEvaluations;s.op=spatialOperator(mesh,lookup,gas,s.cells,control,diffusion);return s;
    };
    auto defect=[&](const ImplicitStage& s) {
        NewtonVector f(4*n);
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<4;++k)
            f[4*i+k]=(s.cells[i][k]-base[i][k]+h/mesh.cells[i].area*s.op.residual[i][k])/scale[4*i+k];
        return f;
    };
    auto current=evaluate(x);auto f=defect(current);
    for(std::size_t it=0;it<=control.maximumNewtonIterations;++it) {
        double largest=0;for(double a:f)largest=std::max(largest,std::abs(a));
        if(largest<=control.nonlinearTolerance)return current;
        require(it<control.maximumNewtonIterations,"implicit Newton budget exhausted; previous accepted state retained");
        ++work.nonlinearIterations;
        NewtonVector diagonal(4*n),rhs=f;
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<4;++k) {
            diagonal[4*i+k]=1+h*(current.op.spectral[i]/mesh.cells[i].area+current.op.heatRate[i]+current.op.viscousRate[i]);
            rhs[4*i+k]=-rhs[4*i+k];
        }
        const auto apply=[&](const NewtonVector& v) {
            NewtonVector direction(v.size());for(std::size_t j=0;j<v.size();++j)direction[j]=v[j]/diagonal[j];
            const double length=newtonNorm(direction);if(length==0)return direction;
            double epsilon=std::sqrt(std::numeric_limits<double>::epsilon())*(1+newtonNorm(x))/length;
            ImplicitStage shifted;bool feasible=false;
            for(unsigned trial=0;trial<12;++trial) {
                auto y=x;for(std::size_t j=0;j<y.size();++j)y[j]+=epsilon*direction[j];
                try {shifted=evaluate(y);feasible=true;break;}catch(const std::runtime_error&){epsilon*=-.5;}
            }
            require(feasible,"implicit Jacobian perturbation is inadmissible");
            NewtonVector product(v.size());
            for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<4;++k) {
                const auto j=4*i+k;
                product[j]=direction[j]+h/mesh.cells[i].area*(shifted.op.residual[i][k]-current.op.residual[i][k])/(epsilon*scale[j]);
            }
            return product;
        };
        auto delta=eulerGmres(apply,rhs,control.maximumKrylovIterations,work.linearIterations);
        for(std::size_t j=0;j<delta.size();++j)delta[j]/=diagonal[j];
        const double oldNorm=newtonNorm(f);bool accepted=false;
        for(double fraction=1;fraction>=1./4096;fraction*=.5) {
            auto y=x;for(std::size_t j=0;j<y.size();++j)y[j]+=fraction*delta[j];
            try {
                auto candidate=evaluate(y);auto next=defect(candidate);
                if(newtonNorm(next)<oldNorm*(1-1e-4*fraction)) {
                    x=std::move(y);current=std::move(candidate);f=std::move(next);accepted=true;break;
                }
            }catch(const std::runtime_error&){} // Reject physical/domain-invalid candidates, never repair them.
        }
        require(accepted,"implicit feasible line search failed; previous accepted state retained");
    }
    throw std::runtime_error("Euler implicit: unreachable Newton exit");
}
}
static EulerStepResult2D advanceEulerImpl(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,
    const IdealGas2D& gas,const EulerState2D& initial,const EulerStepControls2D& control,const EulerDiffusionOperators2D& diffusion) {
    const auto* heat=diffusion.heat?&*diffusion.heat:nullptr;const auto* viscous=diffusion.viscous?&*diffusion.viscous:nullptr;
    validateFvMesh2D(mesh);validateEulerBoundaries2D(mesh,boundaries,gas);
    require(initial.cells.size()==mesh.cells.size()&&std::isfinite(initial.time)&&initial.time>=0,
        "initial state does not match mesh or physical time");
    require(initial.steps<std::numeric_limits<std::size_t>::max(),"step counter overflow");
    require(std::isfinite(control.maximumStep)&&control.maximumStep>0&&std::isfinite(control.minimumStep)&&control.minimumStep>0&&
            control.minimumStep<=control.maximumStep&&std::isfinite(control.acousticCourant)&&control.acousticCourant>0&&
            control.acousticCourant<=.45&&control.maximumRetries<=30,"invalid explicit acoustic time controls");
    require(control.wallGradient==WallGradient2D::Linear||control.wallGradient==WallGradient2D::Quadratic||control.wallGradient==WallGradient2D::FaceQuadratic,"invalid wall gradient scheme");
    require(control.order==1||control.order==2,"spatial/time order must be 1 or 2");
    require(control.timeStepControl==EulerTimeStepControl2D::Legacy||control.timeStepControl==EulerTimeStepControl2D::StageGuarded,"unknown time-step control");
    require(control.integrator==EulerTimeIntegrator2D::Explicit||control.integrator==EulerTimeIntegrator2D::Sdirk2,"unknown time integrator");
    const bool implicit=control.integrator==EulerTimeIntegrator2D::Sdirk2;
    if(implicit)require(std::isfinite(control.nonlinearTolerance)&&control.nonlinearTolerance>0&&control.nonlinearTolerance<=2e-14&&control.maximumNewtonIterations>0&&control.maximumKrylovIterations>0,"invalid implicit solver controls");
    const auto nc=mesh.cells.size(),nf=mesh.faces.size();
    std::vector<const EulerBoundary2D*> lookup(nf,nullptr);for(const auto& b:boundaries)lookup[b.face]=&b;
    require(!control.endTime||(std::isfinite(*control.endTime)&&*control.endTime>initial.time),"invalid integration end time");
    const auto first=spatialOperator(mesh,lookup,gas,initial.cells,control,diffusion);
    const bool guarded=control.timeStepControl==EulerTimeStepControl2D::StageGuarded;
    double cflStep=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<nc;++i)cflStep=std::min(cflStep,diffusion.active()?control.acousticCourant/(first.spectral[i]/mesh.cells[i].area+first.heatRate[i]+first.viscousRate[i]):control.acousticCourant*mesh.cells[i].area/first.spectral[i]);
    double dt=implicit?control.maximumStep:std::min(control.maximumStep,(guarded&&control.order==2?.95:1)*cflStep);
    const auto limitHorizon=[&](double step) {
      if(control.endTime) {
        const double remaining=*control.endTime-initial.time;step=std::min(step,remaining);
        // Keep every step within BOTH user limits. A roundoff-sized final tail
        // is avoided by taking two ordinary smaller steps, never by advancing
        // the clock without flux or silently reducing the declared minimum.
        if(remaining>step&&remaining-step<control.minimumStep)step=.5*remaining;
      }
      return step;
    };
    dt=limitHorizon(dt);
    require(std::isfinite(dt)&&dt>=control.minimumStep,"combined acoustic/heat/viscous CFL requires a step below the declared minimum");
    EulerStepResult2D result;result.spatialEvaluations=1;result.quadraticHeatWalls=heat?heat->quadraticWalls():0;result.quadraticViscousWalls=viscous?viscous->quadraticWalls():0;result.previousCells=initial.cells;
    std::vector<EulerConservative2D> accepted,stage,secondEuler,residual,absolute(nc);
    std::vector<double> spectral;
    std::string failure="non-positive stage state";
    for(std::size_t attempt=0;;++attempt) {
        require(std::isfinite(initial.time+dt)&&initial.time+dt>initial.time,"physical time cannot advance at this step");
        bool valid=implicit?false:positiveUpdate(mesh,initial.cells,first.residual,dt,gas,stage);
        bool cflExceeded=false,positiveStages=valid;
        double stageBound=dt;
        SpatialOperator combined=first;
        if(implicit) {
            try {
                const double gamma=1-1/std::sqrt(2.);
                const auto a=implicitStage(mesh,lookup,gas,initial.cells,initial.cells,gamma*dt,control,diffusion,result);
                auto base=initial.cells;
                for(std::size_t i=0;i<nc;++i)for(std::size_t k=0;k<4;++k)base[i][k]-=(1-gamma)*dt/mesh.cells[i].area*a.op.residual[i][k];
                auto b=implicitStage(mesh,lookup,gas,base,a.cells,gamma*dt,control,diffusion,result);
                combined=b.op;
                std::fill(combined.spectral.begin(),combined.spectral.end(),0.);
                for(std::size_t id=0;id<nf;++id) {
                    for(std::size_t k=0;k<4;++k)combined.faceFlux[id][k]=(1-gamma)*a.op.faceFlux[id][k]+gamma*b.op.faceFlux[id][k];
                    combined.heatFlux[id]=(1-gamma)*a.op.heatFlux[id]+gamma*b.op.heatFlux[id];
                    for(std::size_t k=0;k<3;++k)combined.viscousFlux[id][k]=(1-gamma)*a.op.viscousFlux[id][k]+gamma*b.op.viscousFlux[id][k];
                    combined.speed[id]=std::max(a.op.speed[id],b.op.speed[id]);
                    // Match the exported per-face stage envelope, as for SSPRK2.
                    const auto& face=mesh.faces[id];
                    const double contribution=combined.speed[id]*std::hypot(face.areaVector.x,face.areaVector.y);
                    combined.spectral[face.owner]+=contribution;
                    if(face.neighbour)combined.spectral[*face.neighbour]+=contribution;
                    combined.fallback[id]=static_cast<unsigned char>(a.op.fallback[id]|(b.op.fallback[id]<<1));
                }
                for(std::size_t i=0;i<nc;++i) {
                    for(std::size_t k=0;k<4;++k)combined.residual[i][k]=(1-gamma)*a.op.residual[i][k]+gamma*b.op.residual[i][k];
                    combined.heatRate[i]=std::max(a.op.heatRate[i],b.op.heatRate[i]);
                    combined.viscousRate[i]=std::max(a.op.viscousRate[i],b.op.viscousRate[i]);
                }
                combined.fallbackEvaluations+=a.op.fallbackEvaluations;combined.reconstructionFallbackCells+=a.op.reconstructionFallbackCells;
                combined.minimumContactRestoration=std::min(a.op.minimumContactRestoration,b.op.minimumContactRestoration);
                // RK output is ALWAYS the conservative quadrature of converged
                // stage fluxes. Copying a finite-tolerance Newton stage would
                // accumulate algebraic mass/energy defects, especially in tiny
                // final steps. This is the RK update, not a positivity repair.
                valid=positiveUpdate(mesh,initial.cells,combined.residual,dt,gas,accepted);
                require(valid,"implicit conservative RK output is non-positive");
                for(std::size_t i=0;i<nc;++i) {
                    const auto p=eulerPrimitive2D(initial.cells[i],gas);const double c=eulerSoundSpeed2D(p,gas);
                    const EulerConservative2D scale{p.density,p.density*c,p.density*c,initial.cells[i][3]};
                    // Sum of two stage defects: bound is (1+(1-g)/g)*tol.
                    for(std::size_t k=0;k<4;++k)require(std::abs(accepted[i][k]-b.cells[i][k])<=8*control.nonlinearTolerance*scale[k],"implicit stage/output defect exceeds nonlinear solve tolerance");
                }
            }catch(const std::runtime_error& e){valid=false;failure=e.what();}
        }else if(valid&&control.order==2) {
            try {
                ++result.spatialEvaluations;
                const auto second=spatialOperator(mesh,lookup,gas,stage,control,diffusion);
                // SSPRK(2,2): U1=Un+dt L(Un); Un+1=1/2 Un+1/2 [U1+dt L(U1)].
                // Require admissibility of each FE stage, not just of the mixture.
                valid=positiveUpdate(mesh,stage,second.residual,dt,gas,secondEuler);
                positiveStages=valid;
                combined.fallbackEvaluations+=second.fallbackEvaluations;
                combined.reconstructionFallbackCells+=second.reconstructionFallbackCells;
                combined.minimumContactRestoration=std::min(combined.minimumContactRestoration,second.minimumContactRestoration);
                std::fill(combined.spectral.begin(),combined.spectral.end(),0.);
                std::fill(combined.residual.begin(),combined.residual.end(),EulerConservative2D{});
                for(std::size_t id=0;id<nf;++id) {
                    const auto& face=mesh.faces[id];const double length=std::hypot(face.areaVector.x,face.areaVector.y);
                    combined.speed[id]=std::max(first.speed[id],second.speed[id]);
                    combined.heatFlux[id]=.5*(first.heatFlux[id]+second.heatFlux[id]);
                    for(std::size_t k=0;k<3;++k)combined.viscousFlux[id][k]=.5*(first.viscousFlux[id][k]+second.viscousFlux[id][k]);
                    combined.fallback[id]=static_cast<unsigned char>(first.fallback[id]|(second.fallback[id]<<1));
                    combined.spectral[face.owner]+=combined.speed[id]*length;
                    if(face.neighbour)combined.spectral[*face.neighbour]+=combined.speed[id]*length;
                    for(std::size_t k=0;k<4;++k) {
                        const double value=.5*(first.faceFlux[id][k]+second.faceFlux[id][k]);combined.faceFlux[id][k]=value;
                        combined.residual[face.owner][k]+=value;
                        if(face.neighbour)combined.residual[*face.neighbour][k]-=value;
                    }
                }
                // Bound the combined acoustic and full corrected thermal operator
                // in BOTH stages. Density-dependent volumetric cv changes at U1.
                for(std::size_t i=0;i<nc;++i) {
                    combined.heatRate[i]=std::max(first.heatRate[i],second.heatRate[i]);
                    combined.viscousRate[i]=std::max(first.viscousRate[i],second.viscousRate[i]);
                    const double rate=combined.spectral[i]/mesh.cells[i].area+combined.heatRate[i]+combined.viscousRate[i];
                    stageBound=std::min(stageBound,control.acousticCourant/rate);
                    if(dt*rate>control.acousticCourant*(1+8*std::numeric_limits<double>::epsilon())) {
                        valid=false;cflExceeded=true;failure="second-stage combined acoustic/heat/viscous CFL exceeded";
                    }
                }
                if(valid)valid=positiveUpdate(mesh,initial.cells,combined.residual,dt,gas,accepted);
            }catch(const std::runtime_error& e){valid=false;positiveStages=false;failure=e.what();}
        }else if(valid)accepted=stage;
        if(valid) {
            result.faceFlux=std::move(combined.faceFlux);result.faceWaveSpeed=std::move(combined.speed);
            result.faceHeatFlux=std::move(combined.heatFlux);result.cellHeatRate=std::move(combined.heatRate);
            result.faceViscousFlux=std::move(combined.viscousFlux);result.cellViscousRate=std::move(combined.viscousRate);
            result.heatNonMonotoneRows=heat?heat->nonMonotoneRows():0;result.heatMonotonicityAssessed=!diffusion.hybridHeat;
            result.faceHllcFallbackStages=std::move(combined.fallback);
            result.hllcFallbackEvaluations=combined.fallbackEvaluations;
            result.reconstructionFallbackCells=combined.reconstructionFallbackCells;
            result.minimumContactRestoration=combined.minimumContactRestoration;
            residual=std::move(combined.residual);spectral=std::move(combined.spectral);break;
        }
        require(attempt<control.maximumRetries,"stage positivity/CFL retry budget exhausted ("+failure+"); previous accepted state retained");
        ++result.rejectedCandidates;if(cflExceeded)++result.cflRejectedCandidates;
        // Only a CFL-only rejection may use the measured rate. Any positivity
        // or boundary/operator failure retains the conservative half-step retry.
        dt=guarded&&cflExceeded&&positiveStages?limitHorizon(.95*stageBound):.5*dt;
        require(dt>=control.minimumStep,"stage positivity/CFL requires a step below the declared minimum; previous accepted state retained");
    }
    result.step=dt;result.state={initial.time+dt,initial.steps+1,std::move(accepted)};
    result.minimumDensity=std::numeric_limits<double>::infinity();result.minimumPressure=result.minimumDensity;
    for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];const auto* bc=lookup[id];
        if(!face.neighbour&&!bc->partner){result.boundaryHeat+=result.faceHeatFlux[id];result.boundaryViscousWork+=result.faceViscousFlux[id][2];}
        for(std::size_t k=0;k<4;++k) {
            const double value=result.faceFlux[id][k];absolute[face.owner][k]+=std::abs(value);
            if(face.neighbour)absolute[*face.neighbour][k]+=std::abs(value);
            else if(!bc->partner)result.boundaryFlux[k]=finite(result.boundaryFlux[k]+value);
        }
    }
    for(std::size_t i=0;i<nc;++i) {
        const double area=mesh.cells[i].area;result.acousticCourant=std::max(result.acousticCourant,dt*spectral[i]/area);
        result.thermalCourant=std::max(result.thermalCourant,dt*result.cellHeatRate[i]);
        result.viscousCourant=std::max(result.viscousCourant,dt*result.cellViscousRate[i]);
        result.combinedCourant=std::max(result.combinedCourant,dt*(spectral[i]/area+result.cellHeatRate[i]+result.cellViscousRate[i]));
        const auto q=eulerPrimitive2D(result.state.cells[i],gas);
        result.minimumDensity=std::min(result.minimumDensity,q.density);result.minimumPressure=std::min(result.minimumPressure,q.pressure);
        for(std::size_t k=0;k<4;++k) {
            result.beforeIntegral[k]=finite(result.beforeIntegral[k]+area*initial.cells[i][k]);
            result.afterIntegral[k]=finite(result.afterIntegral[k]+area*result.state.cells[i][k]);
            const double balance=area*(result.state.cells[i][k]-initial.cells[i][k])+dt*residual[i][k];
            const double scale=area*(std::abs(result.state.cells[i][k])+std::abs(initial.cells[i][k]))+dt*absolute[i][k];
            const double relative=scale>0?std::abs(balance)/scale:std::abs(balance);
            result.maximumCellBalanceError=std::max(result.maximumCellBalanceError,finite(relative));
        }
    }
    for(std::size_t k=0;k<4;++k)result.balanceError[k]=finite(result.afterIntegral[k]-result.beforeIntegral[k]+dt*result.boundaryFlux[k]);
    require(result.maximumCellBalanceError<1e-12,"conservative update failed its cell balance audit");
    return result;
}
namespace {
std::optional<HeatConductionOperator2D> prepareHeat(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,const EulerTransport2D& transport,WallGradient2D wallGradient) {
    validateEulerTransport2D(transport,boundaries);
    if(transport.thermalConductivity==0)return {};
    std::vector<HeatBoundary2D> thermal;thermal.reserve(boundaries.size());
    for(const auto& b:boundaries)thermal.push_back({b.face,b.partner?HeatBoundaryKind2D::Periodic:b.thermalKind,b.thermalValue,b.partner});
    return HeatConductionOperator2D(mesh,thermal,transport.thermalConductivity,wallGradient);
}
std::optional<ViscousStressOperator2D> prepareViscous(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,const EulerTransport2D& transport,WallGradient2D wallGradient) {
    validateEulerTransport2D(transport,boundaries);
    if(transport.dynamicViscosity==0)return {};
    std::vector<ViscousBoundary2D> bc;bc.reserve(boundaries.size());
    for(const auto& b:boundaries) {
        const auto kind=b.partner?ViscousBoundaryKind2D::Periodic:b.kind==EulerBoundaryKind2D::NoSlipWall?ViscousBoundaryKind2D::Velocity:
            b.kind==EulerBoundaryKind2D::SlipWall?ViscousBoundaryKind2D::Slip:ViscousBoundaryKind2D::ZeroTraction;
        bc.push_back({b.face,kind,b.wallVelocity,b.partner});
    }
    return ViscousStressOperator2D(mesh,bc,transport.dynamicViscosity,wallGradient);
}
std::shared_ptr<const EulerDiffusionOperators2D> prepareDiffusion(const FvMesh2D& mesh,
    const std::vector<EulerBoundary2D>& boundaries,const EulerTransport2D& transport,
    WallGradient2D wallGradient,EulerDiffusionScheme2D scheme) {
    validateEulerTransport2D(transport,boundaries);
    require(scheme==EulerDiffusionScheme2D::Corrected||scheme==EulerDiffusionScheme2D::HybridHeat||scheme==EulerDiffusionScheme2D::Hybrid,"invalid diffusion scheme");
    require(scheme==EulerDiffusionScheme2D::Corrected||wallGradient==WallGradient2D::Linear,"hybrid diffusion requires the linear wall option; it uses its own geometric form");
    auto out=std::make_shared<EulerDiffusionOperators2D>();
    if(scheme==EulerDiffusionScheme2D::Corrected)out->heat=prepareHeat(mesh,boundaries,transport,wallGradient);
    else if(transport.thermalConductivity>0) {
        std::vector<HeatBoundary2D> bc;for(const auto& b:boundaries)bc.push_back({b.face,b.partner?HeatBoundaryKind2D::Periodic:b.thermalKind,b.thermalValue,b.partner});
        out->hybridHeat.emplace(mesh,bc,transport.thermalConductivity);
        // The cell Schur complement minimizes the nonnegative hybrid energy,
        // hence its mass-scaled spectrum is bounded by the local cell pivots.
        // This bound is not an M-matrix or full compressible positivity claim.
        out->heatCellStiffness.resize(mesh.cells.size());
        out->hybridHeat->visitLocalMatrices([&](std::size_t c,std::size_t,std::size_t,double a){out->heatCellStiffness[c]+=a;});
        for(double c:out->heatCellStiffness)require(std::isfinite(c)&&c>0,"invalid hybrid heat cell bound");
    }
    if(scheme!=EulerDiffusionScheme2D::Hybrid)out->viscous=prepareViscous(mesh,boundaries,transport,wallGradient);
    else if(transport.dynamicViscosity>0) {
        std::vector<ViscousBoundary2D> bc;
        for(const auto& b:boundaries) {
            const auto kind=b.partner?ViscousBoundaryKind2D::Periodic:b.kind==EulerBoundaryKind2D::NoSlipWall?ViscousBoundaryKind2D::Velocity:b.kind==EulerBoundaryKind2D::SlipWall?ViscousBoundaryKind2D::Slip:ViscousBoundaryKind2D::ZeroTraction;
            bc.push_back({b.face,kind,b.wallVelocity,b.partner});
        }
        out->hybridViscous.emplace(mesh,bc,transport.dynamicViscosity);
    }
    return out;
}
}
EulerStepper2D::EulerStepper2D(FvMesh2D mesh,std::vector<EulerBoundary2D> boundaries,IdealGas2D gas,EulerTransport2D transport,WallGradient2D wallGradient,EulerDiffusionScheme2D diffusionScheme)
    :wallGradient_(wallGradient),diffusionScheme_(diffusionScheme),mesh_(std::move(mesh)),boundaries_(std::move(boundaries)),gas_(gas) {
    validateFvMesh2D(mesh_);validateEulerBoundaries2D(mesh_,boundaries_,gas_);
    require(wallGradient==WallGradient2D::Linear||wallGradient==WallGradient2D::Quadratic||wallGradient==WallGradient2D::FaceQuadratic,"invalid wall gradient scheme");
    diffusion_=prepareDiffusion(mesh_,boundaries_,transport,wallGradient,diffusionScheme);
}
EulerStepResult2D EulerStepper2D::advance(const EulerState2D& initial,const EulerStepControls2D& controls) const {
    require(controls.wallGradient==wallGradient_,"wall gradient setting differs from prepared solver");
    require(controls.diffusionScheme==diffusionScheme_,"diffusion scheme differs from prepared solver");
    return advanceEulerImpl(mesh_,boundaries_,gas_,initial,controls,*diffusion_);
}
EulerResidualDiagnostics2D EulerStepper2D::diagnostics(const EulerState2D& state,const EulerStepControls2D& control) const {
    require(state.cells.size()==mesh_.cells.size(),"residual state size differs");
    require(control.diffusionScheme==diffusionScheme_&&control.wallGradient==wallGradient_,"residual controls differ from prepared solver");
    std::vector<const EulerBoundary2D*> lookup(mesh_.faces.size(),nullptr);for(const auto& b:boundaries_)lookup[b.face]=&b;
    const auto op=spatialOperator(mesh_,lookup,gas_,state.cells,control,*diffusion_);
    EulerResidualDiagnostics2D result;
    for(std::size_t face=0;face<mesh_.faces.size();++face)if(!mesh_.faces[face].neighbour&&!lookup[face]->partner)for(std::size_t k=0;k<4;++k)result.boundaryFlux[k]+=op.faceFlux[face][k];
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto p=eulerPrimitive2D(state.cells[i],gas_);const double c=eulerSoundSpeed2D(p,gas_);
        const EulerConservative2D scale{p.density,p.density*c,p.density*c,state.cells[i][3]};
        for(std::size_t k=0;k<4;++k) {result.rate=std::max(result.rate,std::abs(op.residual[i][k])/(mesh_.cells[i].area*scale[k]));result.integralScale[k]+=mesh_.cells[i].area*scale[k];}
    }
    return result;
}
double EulerStepper2D::residualRate(const EulerState2D& state,const EulerStepControls2D& control) const {return diagnostics(state,control).rate;}
EulerStepResult2D advanceEuler2D(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,
    const IdealGas2D& gas,const EulerState2D& initial,const EulerStepControls2D& controls,const EulerTransport2D& transport) {
    const auto diffusion=prepareDiffusion(mesh,boundaries,transport,controls.wallGradient,controls.diffusionScheme);
    return advanceEulerImpl(mesh,boundaries,gas,initial,controls,*diffusion);
}
} // namespace cartmesh2d::fv
