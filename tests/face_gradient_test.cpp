#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include "cartmesh2d/fv/Euler2D.hpp"
#include <iomanip>
#include <iostream>
#include <numeric>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
constexpr double roundoff=4096*std::numeric_limits<double>::epsilon();
double fluxError=0,rateError=0,balanceError=0,orientationError=0;
Vector2D rotate(Vector2D p,double a){return {std::cos(a)*p.x-std::sin(a)*p.y,std::sin(a)*p.x+std::cos(a)*p.y};}
void polynomial(double scale,double angle) {
    auto mesh=fv_test::rectangle(7,5,2.3,true);
    for(auto& c:mesh.cells){auto p=rotate({c.centre.x,c.centre.y},angle);c.centre={p.x*scale,p.y*scale};c.area*=scale*scale;}
    for(auto& f:mesh.faces){auto p=rotate({f.centre.x,f.centre.y},angle),s=rotate(f.areaVector,angle),c=rotate(f.correction,angle);f.centre={scale*p.x,scale*p.y};f.areaVector={scale*s.x,scale*s.y};f.correction={scale*c.x,scale*c.y};}
    auto value=[&](Point2D p){double x=p.x/scale,y=p.y/scale;return 3+.1*x-.2*y+.3*x*x-.15*x*y+.4*y*y;};
    auto grad=[&](Point2D p){return Vector2D{(.1+.6*p.x/scale-.15*p.y/scale)/scale,(-.2-.15*p.x/scale+.8*p.y/scale)/scale};};
    std::vector<HeatBoundary2D> hb;std::vector<ViscousBoundary2D> vb;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){double t=value(mesh.faces[id].centre);hb.push_back({id,HeatBoundaryKind2D::Temperature,t,{}});vb.push_back({id,ViscousBoundaryKind2D::Velocity,{t,-.3*t},{}});}
    std::vector<double> t,capacity(mesh.cells.size(),2),rho(mesh.cells.size(),1.7);std::vector<Vector2D> u;
    for(const auto& c:mesh.cells){t.push_back(value(c.centre));u.push_back({t.back(),-.3*t.back()});}
    const HeatConductionOperator2D heat(mesh,hb,.37,WallGradient2D::FaceQuadratic);
    const ViscousStressOperator2D stress(mesh,vb,.23,WallGradient2D::FaceQuadratic);
    const auto h=heat.evaluate(t,capacity);const auto v=stress.evaluate(u,rho);
    require(heat.quadraticWalls()==hb.size()&&stress.quadraticWalls()==vb.size(),"wrong prescribed face count");
    double wallHeat=0;std::array<double,3> wallStress{};
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];auto g=grad(f.centre);
        double q=-.37*dot(g,f.areaVector);fluxError=std::max(fluxError,std::abs(h.faceHeatFlux[id]-q)/(1+std::abs(q)));
        double div=g.x-.3*g.y,xx=.23*(2*g.x-2./3*div),yy=.23*(-.6*g.y-2./3*div),xy=.23*(g.y-.3*g.x);
        const Vector2D force{-(xx*f.areaVector.x+xy*f.areaVector.y),-(xy*f.areaVector.x+yy*f.areaVector.y)};
        const std::array<double,3> exact{force.x,force.y,value(f.centre)*(force.x-.3*force.y)};
        for(std::size_t k=0;k<3;++k)fluxError=std::max(fluxError,std::abs(v.faceFlux[id][k]-exact[k])/(1+std::abs(exact[k])));
        if(!f.neighbour){wallHeat+=h.faceHeatFlux[id];for(std::size_t k=0;k<3;++k)wallStress[k]+=v.faceFlux[id][k];}
    }
    balanceError=std::max(balanceError,std::abs(std::accumulate(h.cellResidual.begin(),h.cellResidual.end(),0.)-wallHeat)/(1+std::abs(wallHeat)));
    for(std::size_t k=0;k<3;++k){double total=0;for(const auto& r:v.cellResidual)total+=r[k];balanceError=std::max(balanceError,std::abs(total-wallStress[k])/(1+std::abs(wallStress[k])));}
    // Perturb actual operator inputs; include every nonlocal coefficient and
    // boundary-sample subtraction from the owner reference in the row norm.
    std::vector<double> hn(t.size());std::vector<std::array<double,2>> vn(t.size());
    for(std::size_t j=0;j<t.size();++j) {
        auto tp=t,tm=t;tp[j]+=.5;tm[j]-=.5;const auto hp=heat.evaluate(tp,capacity),hm=heat.evaluate(tm,capacity);
        for(std::size_t i=0;i<t.size();++i)hn[i]+=std::abs(hp.cellResidual[i]-hm.cellResidual[i]);
        for(std::size_t a=0;a<2;++a){auto up=u,um=u;(a?up[j].y:up[j].x)+=.5;(a?um[j].y:um[j].x)-=.5;const auto vp=stress.evaluate(up,rho),vm=stress.evaluate(um,rho);for(std::size_t i=0;i<t.size();++i)for(std::size_t k=0;k<2;++k)vn[i][k]+=std::abs(vp.cellResidual[i][k]-vm.cellResidual[i][k]);}
    }
    for(std::size_t i=0;i<t.size();++i){double hr=.5*hn[i]/(mesh.cells[i].area*capacity[i]),vr=.5*std::max(vn[i][0],vn[i][1])/(mesh.cells[i].area*rho[i]);rateError=std::max({rateError,std::abs(h.rate[i]-hr)/hr,std::abs(v.rate[i]-vr)/vr});}
    // Shared reconstruction must not depend on which adjacent cell is owner.
    for(auto& f:mesh.faces)if(f.neighbour){auto old=f.owner;f.owner=*f.neighbour;f.neighbour=old;f.areaVector={-f.areaVector.x,-f.areaVector.y};f.correction={-f.correction.x,-f.correction.y};f.neighbourWeight=1-f.neighbourWeight;}
    for(std::size_t i=0;i<t.size();++i){t[i]=3+.1*std::sin(double(i));u[i]={std::cos(double(i)),std::sin(.7*double(i))};}
    const auto h0=heat.evaluate(t,capacity),h1=HeatConductionOperator2D(mesh,hb,.37,WallGradient2D::FaceQuadratic).evaluate(t,capacity);
    const auto v0=stress.evaluate(u,rho),v1=ViscousStressOperator2D(mesh,vb,.23,WallGradient2D::FaceQuadratic).evaluate(u,rho);
    for(std::size_t i=0;i<t.size();++i){orientationError=std::max(orientationError,std::abs(h0.cellResidual[i]-h1.cellResidual[i])/(1+std::abs(h0.cellResidual[i])));for(std::size_t k=0;k<3;++k)orientationError=std::max(orientationError,std::abs(v0.cellResidual[i][k]-v1.cellResidual[i][k])/(1+std::abs(v0.cellResidual[i][k])));}
}
void periodic() {
    const auto mesh=fv_test::rectangle(6,5,2,true);
    std::vector<HeatBoundary2D> hb;std::vector<ViscousBoundary2D> vb;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour) {
        const auto& f=mesh.faces[id];const bool x=std::abs(f.areaVector.x)>std::abs(f.areaVector.y);std::optional<std::size_t> partner;
        for(std::size_t j=0;j<mesh.faces.size();++j){const auto& p=mesh.faces[j];if(!p.neighbour&&dot(f.areaVector,p.areaVector)<0&&std::abs(x?f.centre.y-p.centre.y:f.centre.x-p.centre.x)<1e-14)partner=j;}
        require(partner.has_value(),"test periodic partner absent");hb.push_back({id,HeatBoundaryKind2D::Periodic,0,partner});vb.push_back({id,ViscousBoundaryKind2D::Periodic,{},partner});
    }
    HeatConductionOperator2D heat(mesh,hb,.37,WallGradient2D::FaceQuadratic);ViscousStressOperator2D stress(mesh,vb,.23,WallGradient2D::FaceQuadratic);
    std::vector<double> t(mesh.cells.size(),3),rho(t.size(),1);std::vector<Vector2D> u(t.size(),{.7,-.2});
    auto h=heat.evaluate(t,rho);auto v=stress.evaluate(u,rho);
    for(double q:h.faceHeatFlux)require(q==0,"constant heat field not exact");for(const auto& f:v.faceFlux)for(double q:f)require(q==0,"constant velocity not exact");
    for(std::size_t i=0;i<t.size();++i){t[i]+=.1*std::sin(double(i));u[i].x+=.1*std::cos(double(i));}
    h=heat.evaluate(t,rho);v=stress.evaluate(u,rho);
    for(const auto& b:hb){require(h.faceHeatFlux[b.face]==-h.faceHeatFlux[*b.partner],"periodic heat is not a shared flux");for(std::size_t k=0;k<3;++k)require(v.faceFlux[b.face][k]==-v.faceFlux[*b.partner][k],"periodic stress/work is not shared");}
    balanceError=std::max(balanceError,std::abs(std::accumulate(h.cellResidual.begin(),h.cellResidual.end(),0.)));
    for(std::size_t k=0;k<3;++k){double total=0;for(const auto& r:v.cellResidual)total+=r[k];balanceError=std::max(balanceError,std::abs(total));}
    std::vector<EulerBoundary2D> eb;
    for(const auto& b:hb){EulerBoundary2D e;e.face=b.face;e.kind=EulerBoundaryKind2D::Periodic;e.partner=b.partner;const auto s=mesh.faces[b.face].areaVector;e.name=std::abs(s.x)>std::abs(s.y)?"periodic-x":"periodic-y";eb.push_back(e);}
    const IdealGas2D gas{1.4,1};const EulerTransport2D transport{.02,.02};
    const EulerStepper2D solver(mesh,eb,gas,transport,WallGradient2D::FaceQuadratic);
    EulerState2D initial;
    for(const auto& c:mesh.cells){double wave=std::sin(std::acos(-1.)*c.centre.x);initial.cells.push_back(eulerConservative2D({1+.001*wave,.02*wave,.01*wave,1+.0014*wave},gas));}
    for(const auto method:{EulerTimeIntegrator2D::Explicit,EulerTimeIntegrator2D::Sdirk2}) {
        EulerStepControls2D controls;controls.order=2;controls.wallGradient=WallGradient2D::FaceQuadratic;
        controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.integrator=method;controls.maximumStep=.002;
        auto state=initial;
        for(int step=0;step<3;++step){auto next=solver.advance(state,controls);require(next.minimumDensity>0&&next.minimumPressure>0,"coupled candidate lost positivity");require(next.maximumCellBalanceError<roundoff,"coupled candidate conservation failure");for(const auto& b:eb)for(std::size_t k=0;k<4;++k)require(next.faceFlux[b.face][k]==-next.faceFlux[*b.partner][k],"coupled periodic flux not shared");state=std::move(next.state);}
        for(std::size_t k=0;k<4;++k){double difference=0,scale=0;for(std::size_t i=0;i<mesh.cells.size();++i){difference+=mesh.cells[i].area*(state.cells[i][k]-initial.cells[i][k]);scale+=mesh.cells[i].area*std::abs(initial.cells[i][k]);}balanceError=std::max(balanceError,std::abs(difference)/(1+scale));}
    }
}
}
int main(){try {
    for(double scale:{1e-3,1.,1e3})for(double angle:{0.,.731})polynomial(scale,angle);
    periodic();
    auto mesh=fv_test::rectangle(1,1,1);std::vector<HeatBoundary2D> bc;for(std::size_t id=0;id<mesh.faces.size();++id)bc.push_back({id,HeatBoundaryKind2D::Temperature,2,{}});
    bool rejected=false;try{HeatConductionOperator2D insufficient(mesh,bc,1,WallGradient2D::FaceQuadratic);}catch(const std::exception& e){rejected=std::string(e.what()).find("rank-deficient")!=std::string::npos;}
    // Existing 4096-epsilon algebra budget, dimensionless after normalization:
    // no new physical accuracy, spectral stability, or positivity claim.
    std::cout<<std::setprecision(17)<<"{\"allFacePolynomialFluxError\":"<<fluxError<<",\"fullRowNormRelativeError\":"<<rateError<<",\"conservationError\":"<<balanceError<<",\"ownerOrientationError\":"<<orientationError<<",\"rankDeficiencyRejected\":"<<(rejected?"true":"false")<<"}\n";
    require(fluxError<roundoff,"shared-face polynomial reproduction failed");require(rateError<roundoff,"shared-face row norm failed");require(balanceError<roundoff,"shared-face conservation failed");require(orientationError<roundoff,"shared-face owner orientation dependence");require(rejected,"rank failure silently accepted");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
