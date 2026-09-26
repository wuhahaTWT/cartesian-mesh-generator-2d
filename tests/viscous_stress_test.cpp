#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* s){if(!ok)throw std::runtime_error(s);}
template<class F> void rejects(F fn){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}require(rejected,"invalid viscous input accepted");}
Vector2D rotate(Vector2D p,double a){return {std::cos(a)*p.x-std::sin(a)*p.y,std::sin(a)*p.x+std::cos(a)*p.y};}
FvMesh2D meshAt(double a){auto m=fv_test::rectangle(8,4,2,true);for(auto& c:m.cells){auto p=rotate({c.centre.x,c.centre.y},a);c.centre={p.x,p.y};}for(auto& f:m.faces){auto p=rotate({f.centre.x,f.centre.y},a);f.centre={p.x,p.y};f.areaVector=rotate(f.areaVector,a);f.correction=rotate(f.correction,a);}return m;}
}
int main(){try {
    double maximumAffineError=0,maximumRateError=0;
    // Affine fields prescribe the exact Dirichlet trace on actual skew faces.
    // Rigid rotation must produce no stress; dilation tests the gas 2/3 term.
    for(double angle:{0.,.731})for(const auto a:{std::array<double,4>{0,0,0,0},{0,-2,2,0},{1,0,0,1},{0,3,0,0},{.4,-.7,1.2,-.9}}) {
        auto mesh=meshAt(angle);const double mu=.37;
        auto velocity=[&](Point2D p){return Vector2D{.3+a[0]*p.x+a[1]*p.y,-.6+a[2]*p.x+a[3]*p.y};};
        std::vector<ViscousBoundary2D> bc;std::vector<Vector2D> u;std::vector<double> rho(mesh.cells.size(),1.7);
        for(const auto& c:mesh.cells)u.push_back(velocity(c.centre));
        for(std::size_t i=0;i<mesh.faces.size();++i)if(!mesh.faces[i].neighbour)bc.push_back({i,ViscousBoundaryKind2D::Velocity,velocity(mesh.faces[i].centre),{}});
        ViscousStressOperator2D op(mesh,bc,mu);const auto result=op.evaluate(u,rho);
        const double div=a[0]+a[3],xx=mu*(2*a[0]-2./3*div),yy=mu*(2*a[3]-2./3*div),xy=mu*(a[1]+a[2]);
        for(std::size_t i=0;i<mesh.faces.size();++i){const auto& f=mesh.faces[i];const auto uv=velocity(f.centre);const Vector2D t{xx*f.areaVector.x+xy*f.areaVector.y,xy*f.areaVector.x+yy*f.areaVector.y};const std::array<double,3> exact{-t.x,-t.y,-dot(uv,t)};for(std::size_t k=0;k<3;++k)maximumAffineError=std::max(maximumAffineError,std::abs(result.faceFlux[i][k]-exact[k])/(1+std::abs(exact[k])));}
        if(a==std::array<double,4>{.4,-.7,1.2,-.9}) {
            std::vector<std::array<double,2>> rowSums(u.size());
            for(std::size_t j=0;j<2*u.size();++j) {
                auto plus=u,minus=u;double& pv=j%2?plus[j/2].y:plus[j/2].x;double& mv=j%2?minus[j/2].y:minus[j/2].x;pv+=1;mv-=1;
                const auto p=op.evaluate(plus,rho),m=op.evaluate(minus,rho);
                for(std::size_t i=0;i<u.size();++i)for(std::size_t k=0;k<2;++k)rowSums[i][k]+=.5*std::abs(p.cellResidual[i][k]-m.cellResidual[i][k]);
            }
            for(std::size_t i=0;i<u.size();++i){const double expected=.5*std::max(rowSums[i][0],rowSums[i][1])/(mesh.cells[i].area*rho[i]);maximumRateError=std::max(maximumRateError,std::abs(result.rate[i]-expected)/expected);}
        }
        auto bad=u;bad[0].x=INFINITY;rejects([&]{(void)op.evaluate(bad,rho);});
    }
    // 1024 eps is an accumulated floating-point operator check, not a flow
    // accuracy target. It is normalized by the physical flux / row rate.
    require(maximumAffineError<1024*std::numeric_limits<double>::epsilon(),"affine stress/work consistency failed");
    require(maximumRateError<1024*std::numeric_limits<double>::epsilon(),"full momentum Jacobian bound differs");
    const auto mesh=fv_test::rectangle(8,4,2,true);std::vector<EulerBoundary2D> bc;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour)bc.push_back({id,EulerBoundaryKind2D::NoSlipWall,{}, {},"wall"});
    const IdealGas2D gas{1.4,1};const EulerTransport2D transport{.5,.2};EulerStepper2D solver(mesh,bc,gas,transport);
    EulerState2D state{0,0,std::vector<EulerConservative2D>(mesh.cells.size(),eulerConservative2D({1,.1,-.2,1},gas))};
    EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=2;
    const auto initial=state;double total0=0;for(std::size_t i=0;i<state.cells.size();++i)total0+=mesh.cells[i].area*state.cells[i][3];
    for(int step=0;step<30;++step){const auto next=solver.advance(state,controls);for(const auto& b:bc)require(next.faceFlux[b.face][0]==0&&next.faceFlux[b.face][3]==0,"stationary insulated wall leaks mass/energy");require(next.combinedCourant<=controls.acousticCourant*(1+1e-13),"combined CFL invalid");state=next.state;}
    double total=0;for(std::size_t i=0;i<state.cells.size();++i)total+=mesh.cells[i].area*state.cells[i][3];require(std::abs(total-total0)<1e-13*total0,"sealed energy not conserved");
    std::stringstream checkpoint;writeEulerCheckpoint2D(checkpoint,mesh,bc,gas,state,"viscous",transport);const auto text=checkpoint.str();require(text.starts_with("CM2D_EULER_CHECKPOINT 3\n"),"viscous checkpoint version");
    const auto restored=readEulerCheckpoint2D(checkpoint,mesh,bc,gas,"viscous",transport);require(restored.cells==state.cells,"restart state differs");
    rejects([&]{std::istringstream in(text);readEulerCheckpoint2D(in,mesh,bc,gas,"viscous",{.5,.3});});
    rejects([&]{EulerStepper2D bad(mesh,bc,gas,{.5,0});});
    auto wrong=bc;wrong[0].wallVelocity=mesh.faces[wrong[0].face].areaVector;rejects([&]{EulerStepper2D bad(mesh,wrong,gas,transport);});
    rejects([&]{EulerStepper2D bad(mesh,bc,gas,{.5,-1});});
    std::cout<<std::setprecision(17)<<"{\"affineStressWorkError\":"<<maximumAffineError<<",\"momentumRowNormRelativeError\":"<<maximumRateError<<",\"sealedEnergyRelativeError\":"<<std::abs(total-total0)/total0<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
