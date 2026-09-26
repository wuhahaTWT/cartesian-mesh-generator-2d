#include "cartmesh2d/fv/Euler2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool valid,const std::string& message) { if(!valid)throw std::runtime_error("Euler: "+message); }
double finite(double x) {require(std::isfinite(x),"non-finite arithmetic");return x;}
void primitiveValid(const EulerPrimitive2D& q) {
    require(std::isfinite(q.density)&&q.density>0&&std::isfinite(q.pressure)&&q.pressure>0&&
            std::isfinite(q.u)&&std::isfinite(q.v),"density and pressure must be strictly positive; velocity must be finite");
}
EulerPrimitive2D ghost(const EulerPrimitive2D& inside,const EulerBoundary2D& b,Vector2D n,const IdealGas2D& gas) {
    const double normal=inside.u*n.x+inside.v*n.y;
    if(b.kind==EulerBoundaryKind2D::SlipWall)
        return {inside.density,inside.u-2*normal*n.x,inside.v-2*normal*n.y,inside.pressure};
    if(b.kind==EulerBoundaryKind2D::Transmissive)return inside;
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
void validateEulerBoundaries2D(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,const IdealGas2D& gas) {
    validateIdealGas2D(gas);std::vector<const EulerBoundary2D*> lookup(mesh.faces.size(),nullptr);
    std::map<std::string,EulerBoundaryKind2D> groups;
    for(const auto& b:boundaries) {
        require(b.face<mesh.faces.size()&&!mesh.faces[b.face].neighbour&&!lookup[b.face],"boundary face missing, duplicated or internal");
        require(!b.name.empty(),"boundary name is empty");
        const auto [it,inserted]=groups.emplace(b.name,b.kind);
        require(inserted||it->second==b.kind,"one boundary name has multiple physical kinds");
        require(b.kind==EulerBoundaryKind2D::SlipWall||b.kind==EulerBoundaryKind2D::Transmissive||
                b.kind==EulerBoundaryKind2D::Farfield||b.kind==EulerBoundaryKind2D::Periodic,"unknown boundary kind");
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
    std::vector<double> speed,spectral;
    std::vector<unsigned char> fallback;
    std::size_t fallbackEvaluations=0,reconstructionFallbackCells=0;
    double minimumContactRestoration=1;
};
SpatialOperator spatialOperator(const FvMesh2D& mesh,const std::vector<const EulerBoundary2D*>& lookup,
    const IdealGas2D& gas,const std::vector<EulerConservative2D>& cells,const EulerStepControls2D& control) {
    const auto nc=mesh.cells.size(),nf=mesh.faces.size();
    std::vector<EulerPrimitive2D> primitive;primitive.reserve(nc);
    for(const auto& u:cells)primitive.push_back(eulerPrimitive2D(u,gas));
    SpatialOperator out;out.faceFlux.resize(nf);out.residual.resize(nc);out.speed.resize(nf);
    out.spectral.resize(nc);out.fallback.resize(nf);
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
        if(boundary&&boundary->kind==EulerBoundaryKind2D::SlipWall) {
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
EulerStepResult2D advanceEuler2D(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,
    const IdealGas2D& gas,const EulerState2D& initial,const EulerStepControls2D& control) {
    validateFvMesh2D(mesh);validateEulerBoundaries2D(mesh,boundaries,gas);
    require(initial.cells.size()==mesh.cells.size()&&std::isfinite(initial.time)&&initial.time>=0,
        "initial state does not match mesh or physical time");
    require(initial.steps<std::numeric_limits<std::size_t>::max(),"step counter overflow");
    require(std::isfinite(control.maximumStep)&&control.maximumStep>0&&std::isfinite(control.minimumStep)&&control.minimumStep>0&&
            control.minimumStep<=control.maximumStep&&std::isfinite(control.acousticCourant)&&control.acousticCourant>0&&
            control.acousticCourant<=.45&&control.maximumRetries<=30,"invalid explicit acoustic time controls");
    require(control.order==1||control.order==2,"spatial/time order must be 1 or 2");
    const auto nc=mesh.cells.size(),nf=mesh.faces.size();
    std::vector<const EulerBoundary2D*> lookup(nf,nullptr);for(const auto& b:boundaries)lookup[b.face]=&b;
    const auto first=spatialOperator(mesh,lookup,gas,initial.cells,control);
    double dt=control.maximumStep;
    for(std::size_t i=0;i<nc;++i)dt=std::min(dt,control.acousticCourant*mesh.cells[i].area/first.spectral[i]);
    require(std::isfinite(dt)&&dt>=control.minimumStep,"acoustic CFL requires a step below the declared minimum");
    EulerStepResult2D result;result.previousCells=initial.cells;
    std::vector<EulerConservative2D> accepted,stage,secondEuler,residual,absolute(nc);
    std::vector<double> spectral;
    std::string failure="non-positive stage state";
    for(std::size_t attempt=0;;++attempt) {
        require(std::isfinite(initial.time+dt)&&initial.time+dt>initial.time,"physical time cannot advance at this step");
        bool valid=positiveUpdate(mesh,initial.cells,first.residual,dt,gas,stage);
        SpatialOperator combined=first;
        if(valid&&control.order==2) {
            try {
                const auto second=spatialOperator(mesh,lookup,gas,stage,control);
                // SSPRK(2,2): U1=Un+dt L(Un); Un+1=1/2 Un+1/2 [U1+dt L(U1)].
                // Require admissibility of each FE stage, not just of the mixture.
                valid=positiveUpdate(mesh,stage,second.residual,dt,gas,secondEuler);
                combined.fallbackEvaluations+=second.fallbackEvaluations;
                combined.reconstructionFallbackCells+=second.reconstructionFallbackCells;
                combined.minimumContactRestoration=std::min(combined.minimumContactRestoration,second.minimumContactRestoration);
                std::fill(combined.spectral.begin(),combined.spectral.end(),0.);
                std::fill(combined.residual.begin(),combined.residual.end(),EulerConservative2D{});
                for(std::size_t id=0;id<nf;++id) {
                    const auto& face=mesh.faces[id];const double length=std::hypot(face.areaVector.x,face.areaVector.y);
                    combined.speed[id]=std::max(first.speed[id],second.speed[id]);
                    combined.fallback[id]=static_cast<unsigned char>(first.fallback[id]|(second.fallback[id]<<1));
                    combined.spectral[face.owner]+=combined.speed[id]*length;
                    if(face.neighbour)combined.spectral[*face.neighbour]+=combined.speed[id]*length;
                    for(std::size_t k=0;k<4;++k) {
                        const double value=.5*(first.faceFlux[id][k]+second.faceFlux[id][k]);combined.faceFlux[id][k]=value;
                        combined.residual[face.owner][k]+=value;
                        if(face.neighbour)combined.residual[*face.neighbour][k]-=value;
                    }
                }
                // The exported face maximum is conservative for both stages.
                for(std::size_t i=0;i<nc;++i)if(dt*combined.spectral[i]/mesh.cells[i].area>control.acousticCourant*(1+8*std::numeric_limits<double>::epsilon())) {
                    valid=false;failure="second-stage acoustic CFL exceeded";break;
                }
                if(valid)valid=positiveUpdate(mesh,initial.cells,combined.residual,dt,gas,accepted);
            }catch(const std::runtime_error& e){valid=false;failure=e.what();}
        }else if(valid)accepted=stage;
        if(valid) {
            result.faceFlux=std::move(combined.faceFlux);result.faceWaveSpeed=std::move(combined.speed);
            result.faceHllcFallbackStages=std::move(combined.fallback);
            result.hllcFallbackEvaluations=combined.fallbackEvaluations;
            result.reconstructionFallbackCells=combined.reconstructionFallbackCells;
            result.minimumContactRestoration=combined.minimumContactRestoration;
            residual=std::move(combined.residual);spectral=std::move(combined.spectral);break;
        }
        require(attempt<control.maximumRetries,"stage positivity/CFL retry budget exhausted ("+failure+"); previous accepted state retained");
        ++result.rejectedCandidates;dt*=.5;
        require(dt>=control.minimumStep,"stage positivity/CFL requires a step below the declared minimum; previous accepted state retained");
    }
    result.step=dt;result.state={initial.time+dt,initial.steps+1,std::move(accepted)};
    result.minimumDensity=std::numeric_limits<double>::infinity();result.minimumPressure=result.minimumDensity;
    for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];const auto* bc=lookup[id];
        for(std::size_t k=0;k<4;++k) {
            const double value=result.faceFlux[id][k];absolute[face.owner][k]+=std::abs(value);
            if(face.neighbour)absolute[*face.neighbour][k]+=std::abs(value);
            else if(!bc->partner)result.boundaryFlux[k]=finite(result.boundaryFlux[k]+value);
        }
    }
    for(std::size_t i=0;i<nc;++i) {
        const double area=mesh.cells[i].area;result.acousticCourant=std::max(result.acousticCourant,dt*spectral[i]/area);
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
} // namespace cartmesh2d::fv
