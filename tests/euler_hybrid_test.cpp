#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <iomanip>
#include <iostream>
#include <numbers>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"invalid/failed hybrid Euler candidate accepted");}
std::vector<EulerBoundary2D> walls(const FvMesh2D& m){std::vector<EulerBoundary2D> b;for(std::size_t f=0;f<m.faces.size();++f)if(!m.faces[f].neighbour){EulerBoundary2D x;x.face=f;x.kind=EulerBoundaryKind2D::NoSlipWall;x.name="wall";b.push_back(x);}return b;}
double worstBalance=0,worstClosedEnergy=0;
void inspect(const FvMesh2D& m,const EulerStepResult2D& r) {
    require(r.minimumPressure>0&&r.minimumDensity>0,"hybrid Euler positivity failed");
    require(!r.heatMonotonicityAssessed,"hybrid heat incorrectly claims counted monotonicity rows");
    worstBalance=std::max(worstBalance,r.maximumCellBalanceError);
    worstClosedEnergy=std::max(worstClosedEnergy,std::abs(r.afterIntegral[3]-r.beforeIntegral[3])/(std::abs(r.beforeIntegral[3])+std::abs(r.afterIntegral[3])));
    // Reconstruct the ACTUAL exported conservative increment, not a second PDE.
    std::vector<EulerConservative2D> residual(m.cells.size());
    for(std::size_t f=0;f<m.faces.size();++f)for(std::size_t k=0;k<4;++k){residual[m.faces[f].owner][k]+=r.faceFlux[f][k];if(m.faces[f].neighbour)residual[*m.faces[f].neighbour][k]-=r.faceFlux[f][k];}
    for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k){const double delta=r.state.cells[c][k]-r.previousCells[c][k]+r.step/m.cells[c].area*residual[c][k];require(std::abs(delta)<1e-12*(std::abs(r.state.cells[c][k])+std::abs(r.previousCells[c][k])+1),"hybrid Euler output not exported RK flux update");}
}
void coupled() {
    const auto m=fv_test::rectangle(4,3,2,true);const auto bc=walls(m);const IdealGas2D gas{1.4,1};const EulerTransport2D physics{.02,.02};
    EulerState2D initial;
    for(const auto& c:m.cells){const double x=std::numbers::pi*c.centre.x/2,y=std::numbers::pi*c.centre.y;initial.cells.push_back(eulerConservative2D({1,.02*std::sin(x)*std::sin(y),.01*std::cos(x)*std::sin(y),1+.003*std::cos(x)*std::cos(y)},gas));}
    EulerStepControls2D control;control.diffusionScheme=EulerDiffusionScheme2D::Hybrid;control.integrator=EulerTimeIntegrator2D::Sdirk2;control.order=2;control.fluxScheme=EulerFluxScheme2D::Hllc;control.endTime=.02;
    const EulerStepper2D solver(m,bc,gas,physics,control.wallGradient,control.diffusionScheme);
    const auto evolve=[&](double dt) {
        auto s=initial;control.maximumStep=dt;
        while(s.time<*control.endTime){const auto next=solver.advance(s,control);inspect(m,next);
            if(s.steps==1){std::stringstream file;writeEulerCheckpoint2D(file,m,bc,gas,s,"hybrid-micro",physics);const auto restored=readEulerCheckpoint2D(file,m,bc,gas,"hybrid-micro",physics);const auto resumed=EulerStepper2D(m,bc,gas,physics,control.wallGradient,control.diffusionScheme).advance(restored,control);require(next.state.cells==resumed.state.cells&&next.state.time==resumed.state.time,"hybrid restart has hidden solver history");}
            s=next.state;
        }
        return s;
    };
    const auto coarse=evolve(.004),medium=evolve(.002),fine=evolve(.001),reference=evolve(.00025),checked=evolve(.0005);
    const auto error=[&](const EulerState2D& s){double total=0;for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k)total+=m.cells[c].area*std::pow(s.cells[c][k]-reference.cells[c][k],2);return std::sqrt(total);};
    const double ec=error(coarse),em=error(medium),ef=error(fine),order=std::log2(em/ef);
    std::cout<<std::setprecision(17)<<"{\"coarseError\":"<<ec<<",\"mediumError\":"<<em<<",\"fineError\":"<<ef<<",\"observedOrder\":"<<order<<",\"referenceCheckError\":"<<error(checked)<<"}\n";
    require(ec>em&&order>1.8&&error(checked)<.3*ef,"coupled hybrid SDIRK2 time refinement failed");
    auto explicitControl=control;explicitControl.integrator=EulerTimeIntegrator2D::Explicit;explicitControl.maximumStep=.0001;
    inspect(m,solver.advance(initial,explicitControl));
    explicitControl.diffusionScheme=EulerDiffusionScheme2D::HybridHeat;
    inspect(m,EulerStepper2D(m,bc,gas,physics,WallGradient2D::Linear,explicitControl.diffusionScheme).advance(initial,explicitControl));
    rejects([&]{(void)solver.advance(initial,explicitControl);});
    control.maximumStep=1;control.endTime.reset();control.maximumNewtonIterations=1;control.maximumRetries=0;const auto before=initial;
    rejects([&]{(void)solver.advance(initial,control);});require(initial.cells==before.cells&&initial.time==before.time,"failed hybrid stage overwrote accepted state");
    auto slip=bc;for(auto& b:slip)b.kind=EulerBoundaryKind2D::SlipWall;
    rejects([&]{EulerStepper2D unsupported(m,slip,gas,physics,WallGradient2D::Linear,EulerDiffusionScheme2D::Hybrid);});
}
void siScale() {
    auto m=fv_test::rectangle(3,2,2,true);for(auto& c:m.cells){c.centre.x*=1e-4;c.centre.y*=1e-4;c.area*=1e-8;}for(auto& f:m.faces){f.centre.x*=1e-4;f.centre.y*=1e-4;f.areaVector.x*=1e-4;f.areaVector.y*=1e-4;f.correction.x*=1e-4;f.correction.y*=1e-4;}
    const auto bc=walls(m);const IdealGas2D gas;const EulerTransport2D physics{.025759,1.846e-5};EulerState2D s;
    for(const auto& c:m.cells){const double t=300+std::sin(c.centre.x*10000);s.cells.push_back(eulerConservative2D({101325/(gas.gasConstant*t),2*std::sin(c.centre.y*10000),0,101325},gas));}
    EulerStepControls2D control;control.diffusionScheme=EulerDiffusionScheme2D::Hybrid;control.integrator=EulerTimeIntegrator2D::Sdirk2;control.order=2;control.maximumStep=1e-8;
    const EulerStepper2D solver(m,bc,gas,physics,control.wallGradient,control.diffusionScheme);
    inspect(m,solver.advance(s,control));
}
}
int main(){try{coupled();siScale();std::cout<<"{\"maximumCellBalance\":"<<worstBalance<<",\"maximumClosedEnergyDriftPerStep\":"<<worstClosedEnergy<<"}\n";require(worstBalance<1e-12&&worstClosedEnergy<1e-12,"coupled hybrid conservative balance failed");return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
