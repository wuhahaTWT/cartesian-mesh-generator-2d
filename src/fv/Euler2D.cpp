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
EulerConservative2D flux(const EulerConservative2D& q,const EulerPrimitive2D& p,Vector2D n) {
    const double velocity=finite(p.u*n.x+p.v*n.y);
    return {finite(q[0]*velocity),finite(q[1]*velocity+p.pressure*n.x),
        finite(q[2]*velocity+p.pressure*n.y),finite((q[3]+p.pressure)*velocity)};
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

EulerStepResult2D advanceEuler2D(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,
    const IdealGas2D& gas,const EulerState2D& initial,const EulerStepControls2D& control) {
    validateFvMesh2D(mesh);validateEulerBoundaries2D(mesh,boundaries,gas);
    require(initial.cells.size()==mesh.cells.size()&&std::isfinite(initial.time)&&initial.time>=0,
        "initial state does not match mesh or physical time");
    require(initial.steps<std::numeric_limits<std::size_t>::max(),"step counter overflow");
    require(std::isfinite(control.maximumStep)&&control.maximumStep>0&&std::isfinite(control.minimumStep)&&control.minimumStep>0&&
            control.minimumStep<=control.maximumStep&&std::isfinite(control.acousticCourant)&&control.acousticCourant>0&&
            control.acousticCourant<=.45&&control.maximumRetries<=30,"invalid explicit acoustic time controls");
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    std::vector<const EulerBoundary2D*> lookup(nf,nullptr);
    for(const auto& b:boundaries)lookup[b.face]=&b;
    std::vector<EulerPrimitive2D> primitive;primitive.reserve(n);
    for(const auto& q:initial.cells)primitive.push_back(eulerPrimitive2D(q,gas));
    EulerStepResult2D result;result.faceFlux.resize(nf);result.faceWaveSpeed.resize(nf);result.previousCells=initial.cells;
    std::vector<EulerConservative2D> residual(n),absoluteFlux(n);
    std::vector<double> spectral(n);
    for(std::size_t id=0;id<nf;++id) {
        const auto& face=mesh.faces[id];const auto* boundary=lookup[id];
        if(boundary&&boundary->partner&&id>*boundary->partner)continue;
        const auto owner=face.owner;auto neighbour=face.neighbour;
        if(boundary&&boundary->partner)neighbour=mesh.faces[*boundary->partner].owner;
        const double length=std::hypot(face.areaVector.x,face.areaVector.y);
        const Vector2D normal{face.areaVector.x/length,face.areaVector.y/length};
        const auto left=primitive[owner];const auto right=neighbour?primitive[*neighbour]:ghost(left,*boundary,normal,gas);
        const auto qleft=initial.cells[owner],qright=neighbour?initial.cells[*neighbour]:eulerConservative2D(right,gas);
        const auto fleft=flux(qleft,left,normal),fright=flux(qright,right,normal);
        const double speed=finite(std::max(std::abs(left.u*normal.x+left.v*normal.y)+eulerSoundSpeed2D(left,gas),
            std::abs(right.u*normal.x+right.v*normal.y)+eulerSoundSpeed2D(right,gas)));
        require(speed>0,"non-positive acoustic speed");result.faceWaveSpeed[id]=speed;
        spectral[owner]=finite(spectral[owner]+speed*length);
        if(neighbour)spectral[*neighbour]=finite(spectral[*neighbour]+speed*length);
        for(std::size_t k=0;k<4;++k) {
            const double value=finite(.5*length*(fleft[k]+fright[k]-speed*(qright[k]-qleft[k])));
            result.faceFlux[id][k]=value;residual[owner][k]=finite(residual[owner][k]+value);absoluteFlux[owner][k]+=std::abs(value);
            if(neighbour){residual[*neighbour][k]=finite(residual[*neighbour][k]-value);absoluteFlux[*neighbour][k]+=std::abs(value);}
            else result.boundaryFlux[k]=finite(result.boundaryFlux[k]+value);
            if(boundary&&boundary->partner)result.faceFlux[*boundary->partner][k]=-value;
        }
        if(boundary&&boundary->partner)result.faceWaveSpeed[*boundary->partner]=speed;
    }
    double dt=control.maximumStep;
    for(std::size_t i=0;i<n;++i)dt=std::min(dt,control.acousticCourant*mesh.cells[i].area/spectral[i]);
    require(std::isfinite(dt)&&dt>=control.minimumStep,"acoustic CFL requires a step below the declared minimum");
    for(std::size_t attempt=0;;++attempt) {
        require(std::isfinite(initial.time+dt)&&initial.time+dt>initial.time,"physical time cannot advance at this step");
        EulerState2D candidate{initial.time+dt,initial.steps+1,initial.cells};bool positive=true;
        double density=std::numeric_limits<double>::infinity(),pressure=density;
        for(std::size_t i=0;i<n;++i) {
            for(std::size_t k=0;k<4;++k)candidate.cells[i][k]=initial.cells[i][k]-dt/mesh.cells[i].area*residual[i][k];
            try {const auto q=eulerPrimitive2D(candidate.cells[i],gas);density=std::min(density,q.density);pressure=std::min(pressure,q.pressure);}
            catch(const std::runtime_error&){positive=false;break;}
        }
        if(positive){result.state=std::move(candidate);result.minimumDensity=density;result.minimumPressure=pressure;break;}
        require(attempt<control.maximumRetries,"positivity retry budget exhausted; previous accepted state retained");
        ++result.rejectedCandidates;dt*=.5;
        require(dt>=control.minimumStep,"positivity requires a step below the declared minimum; previous accepted state retained");
    }
    result.step=dt;
    for(std::size_t i=0;i<n;++i) {
        const double area=mesh.cells[i].area;result.acousticCourant=std::max(result.acousticCourant,dt*spectral[i]/area);
        for(std::size_t k=0;k<4;++k) {
            result.beforeIntegral[k]=finite(result.beforeIntegral[k]+area*initial.cells[i][k]);
            result.afterIntegral[k]=finite(result.afterIntegral[k]+area*result.state.cells[i][k]);
            const double balance=area*(result.state.cells[i][k]-initial.cells[i][k])+dt*residual[i][k];
            const double scale=area*(std::abs(result.state.cells[i][k])+std::abs(initial.cells[i][k]))+dt*absoluteFlux[i][k];
            const double relative=scale>0?std::abs(balance)/scale:std::abs(balance);
            result.maximumCellBalanceError=std::max(result.maximumCellBalanceError,finite(relative));
        }
    }
    for(std::size_t k=0;k<4;++k)result.balanceError[k]=finite(result.afterIntegral[k]-result.beforeIntegral[k]+dt*result.boundaryFlux[k]);
    require(result.maximumCellBalanceError<1e-12,"conservative update failed its cell balance audit");
    return result;
}
} // namespace cartmesh2d::fv
