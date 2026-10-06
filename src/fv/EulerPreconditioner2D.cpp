#include "cartmesh2d/fv/detail/EulerPreconditioner2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::fv::detail {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("Euler frozen-flux preconditioner: ")+message);}
double checked(double x){require(std::isfinite(x),"nonfinite arithmetic");return x;}
std::vector<std::pair<std::size_t,std::size_t>> graph(const FvMesh2D& mesh,const std::vector<const EulerBoundary2D*>& lookup) {
    require(lookup.size()==mesh.faces.size(),"invalid boundary lookup");
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(std::size_t c=0;c<mesh.cells.size();++c)for(std::size_t i=0;i<4;++i)for(std::size_t j=i+1;j<4;++j)edges.emplace_back(4*c+i,4*c+j);
    for(std::size_t f=0;f<mesh.faces.size();++f) {
        const auto& face=mesh.faces[f];const auto* b=lookup[f];auto neighbour=face.neighbour;
        require(neighbour||b,"missing boundary");
        if(b&&b->partner){require(*b->partner<mesh.faces.size(),"invalid periodic partner");neighbour=mesh.faces[*b->partner].owner;}
        if(neighbour&&*neighbour!=face.owner)for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)edges.emplace_back(4*face.owner+i,4*(*neighbour)+j);
    }
    return edges;
}
EulerMatrix4 boundaryJacobian(const EulerConservative2D& state,const EulerBoundary2D& b,
    Vector2D normal,const IdealGas2D& gas,const LinearVector2D& scale,std::size_t cell,
    const EulerFrozenPreconditioner2D::BoundaryMap& map) {
    EulerMatrix4 result{};
    if(b.kind==EulerBoundaryKind2D::SlipWall||b.kind==EulerBoundaryKind2D::NoSlipWall) {
        result[0][0]=result[3][3]=1;
        result[1][1]=1-2*normal.x*normal.x;result[1][2]=-2*normal.x*normal.y;
        result[2][1]=result[1][2];result[2][2]=1-2*normal.y*normal.y;return result;
    }
    if(b.kind==EulerBoundaryKind2D::Transmissive){for(std::size_t i=0;i<4;++i)result[i][i]=1;return result;}
    const auto base=map(state,b,normal);
    for(std::size_t j=0;j<4;++j) {
        double step=std::cbrt(std::numeric_limits<double>::epsilon())*scale[4*cell+j];
        bool complete=false;
        for(unsigned trial=0;trial<12;++trial) {
            EulerConservative2D plus=state,minus=state,fp{},fm{};plus[j]+=step;minus[j]-=step;
            bool positive=false,negative=false;
            try{(void)eulerPrimitive2D(plus,gas);fp=map(plus,b,normal);positive=true;}catch(const std::runtime_error&){}
            try{(void)eulerPrimitive2D(minus,gas);fm=map(minus,b,normal);negative=true;}catch(const std::runtime_error&){}
            if(positive||negative) {
                for(std::size_t i=0;i<4;++i)result[i][j]=checked(positive&&negative?(fp[i]-fm[i])/(2*step):positive?(fp[i]-base[i])/step:(base[i]-fm[i])/step);
                complete=true;break;
            }
            step*=.5;
        }
        require(complete,"boundary derivative has no admissible perturbation");
    }
    return result;
}
std::array<double,4> temperatureDerivative(const EulerPrimitive2D& p,const IdealGas2D& gas) {
    const double factor=(gas.gamma-1)/(p.density*gas.gasConstant),temperature=p.pressure/(p.density*gas.gasConstant);
    return {factor*.5*(p.u*p.u+p.v*p.v)-temperature/p.density,-factor*p.u,-factor*p.v,factor};
}
std::array<std::array<double,4>,2> velocityDerivative(const EulerPrimitive2D& p) {
    return {{{-p.u/p.density,1/p.density,0,0},{-p.v/p.density,0,1/p.density,0}}};
}
}
EulerMatrix4 eulerNormalFluxJacobian2D(const EulerConservative2D& state,Vector2D s,const IdealGas2D& gas) {
    const auto p=eulerPrimitive2D(state,gas);const double normal=p.u*s.x+p.v*s.y;
    const std::array<double,4> pressure{(gas.gamma-1)*.5*(p.u*p.u+p.v*p.v),-(gas.gamma-1)*p.u,-(gas.gamma-1)*p.v,gas.gamma-1};
    const std::array<double,4> velocity{-normal/p.density,s.x/p.density,s.y/p.density,0};
    EulerMatrix4 a{};a[0]={0,s.x,s.y,0};
    for(std::size_t j=0;j<4;++j) {
        a[1][j]=state[1]*velocity[j]+pressure[j]*s.x+(j==1?normal:0);
        a[2][j]=state[2]*velocity[j]+pressure[j]*s.y+(j==2?normal:0);
        a[3][j]=(state[3]+p.pressure)*velocity[j]+(pressure[j]+(j==3?1:0))*normal;
        for(std::size_t i=0;i<4;++i)checked(a[i][j]);
    }
    return a;
}
EulerFrozenPreconditioner2D::EulerFrozenPreconditioner2D(const FvMesh2D& mesh,
    const std::vector<const EulerBoundary2D*>& lookup,const IdealGas2D& gas,
    const EulerTransport2D& transport,const std::vector<EulerConservative2D>& cells,
    double step,const LinearVector2D& scales,const BoundaryMap& boundaryMap)
    :pattern_(4*mesh.cells.size(),graph(mesh,lookup)),system_(pattern_) {
    require(cells.size()==mesh.cells.size()&&scales.size()==4*cells.size(),"invalid cell/scaling dimensions");
    require(std::isfinite(step)&&step>0,"invalid stage step");
    require(std::isfinite(transport.dynamicViscosity)&&transport.dynamicViscosity>=0&&std::isfinite(transport.thermalConductivity)&&transport.thermalConductivity>=0,"invalid transport coefficients");
    for(double x:scales)require(std::isfinite(x)&&x>0,"invalid variable scale");
    std::fill(system_.diag.begin(),system_.diag.end(),1.);
    const auto addBlock=[&](std::size_t rowCell,std::size_t columnCell,const EulerMatrix4& block,double sign) {
        for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j) {
            const auto row=4*rowCell+i,column=4*columnCell+j;
            const double value=checked(sign*step/mesh.cells[rowCell].area*block[i][j]*scales[column]/scales[row]);
            if(row==column)system_.diag[row]+=value;else system_.add(row,column,value);
        }
    };
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& face=mesh.faces[id];const auto* b=lookup[id];const auto owner=face.owner;
        if(b&&b->partner&&id>*b->partner)continue;
        auto neighbour=face.neighbour;if(b&&b->partner)neighbour=mesh.faces[*b->partner].owner;
        const auto s=face.areaVector;const double length=std::hypot(s.x,s.y);const Vector2D normal{s.x/length,s.y/length};
        const auto& left=cells[owner];const auto right=neighbour?cells[*neighbour]:boundaryMap(left,*b,normal);
        const auto pl=eulerPrimitive2D(left,gas),pr=eulerPrimitive2D(right,gas);
        const double wave=length*std::max(std::abs(pl.u*normal.x+pl.v*normal.y)+eulerSoundSpeed2D(pl,gas),std::abs(pr.u*normal.x+pr.v*normal.y)+eulerSoundSpeed2D(pr,gas));
        auto l=eulerNormalFluxJacobian2D(left,s,gas),r=eulerNormalFluxJacobian2D(right,s,gas);
        for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j){l[i][j]=.5*(l[i][j]+(i==j?wave:0));r[i][j]=.5*(r[i][j]-(i==j?wave:0));}
        if(!neighbour) {
            const auto derivative=boundaryJacobian(left,*b,normal,gas,scales,owner,boundaryMap);
            for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)for(std::size_t k=0;k<4;++k)l[i][j]+=r[i][k]*derivative[k][j];
            r={};
        }
        // Frozen normal-gradient transport approximation, sharing this face's
        // geometry and physical boundary type. The actual residual retains all
        // nonorthogonal/high-order or hybrid terms; this sparse form only helps
        // the Krylov solve. Face velocity in viscous work is frozen here.
        Vector2D distance{face.centre.x-mesh.cells[owner].centre.x,face.centre.y-mesh.cells[owner].centre.y};
        if(neighbour) {
            const auto centre=b&&b->partner?mesh.faces[*b->partner].centre:face.centre;
            distance.x-=centre.x-mesh.cells[*neighbour].centre.x;distance.y-=centre.y-mesh.cells[*neighbour].centre.y;
        }
        const double dn=dot(distance,normal);require(std::isfinite(dn)&&dn>0,"nonpositive normal distance");
        const auto duL=velocityDerivative(pl),duR=velocityDerivative(pr);
        const bool wall=b&&(b->kind==EulerBoundaryKind2D::NoSlipWall||b->kind==EulerBoundaryKind2D::SlipWall);
        if(transport.dynamicViscosity>0&&(neighbour||wall)) {
            const bool slip=!neighbour&&b->kind==EulerBoundaryKind2D::SlipWall;const double mu=transport.dynamicViscosity*length/dn;
            const double n[2]{normal.x,normal.y};double k[2][2]{};
            for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)k[i][j]=mu*((slip?4./3:1./3)*n[i]*n[j]+(!slip&&i==j?1:0));
            const Vector2D velocity=neighbour?Vector2D{.5*(pl.u+pr.u),.5*(pl.v+pr.v)}:slip?Vector2D{}:b->wallVelocity;
            for(std::size_t column=0;column<4;++column) {
                double ownerForce[2]{},neighbourForce[2]{};
                for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j){ownerForce[i]+=k[i][j]*duL[j][column];if(neighbour)neighbourForce[i]-=k[i][j]*duR[j][column];}
                for(std::size_t i=0;i<2;++i){l[i+1][column]+=ownerForce[i];r[i+1][column]+=neighbourForce[i];}
                l[3][column]+=velocity.x*ownerForce[0]+velocity.y*ownerForce[1];r[3][column]+=velocity.x*neighbourForce[0]+velocity.y*neighbourForce[1];
            }
        }
        if(transport.thermalConductivity>0&&(neighbour||b->thermalKind==HeatBoundaryKind2D::Temperature)) {
            const double k=transport.thermalConductivity*length/dn;const auto tl=temperatureDerivative(pl,gas),tr=temperatureDerivative(pr,gas);
            for(std::size_t column=0;column<4;++column){l[3][column]+=k*tl[column];if(neighbour)r[3][column]-=k*tr[column];}
        }
        addBlock(owner,owner,l,1);
        if(neighbour){addBlock(owner,*neighbour,r,1);addBlock(*neighbour,owner,l,-1);addBlock(*neighbour,*neighbour,r,-1);}
    }
    system_.factorILU0();
}
LinearVector2D EulerFrozenPreconditioner2D::apply(const LinearVector2D& rhs) const {
    LinearVector2D out(rhs.size());system_.preconditionILU0(rhs,out);return out;
}
LinearVector2D EulerFrozenPreconditioner2D::applyMatrix(const LinearVector2D& x) const {
    LinearVector2D out;system_.apply(x,out);return out;
}
} // namespace cartmesh2d::fv::detail
