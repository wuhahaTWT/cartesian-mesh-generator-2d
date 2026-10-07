#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <limits>
#include <sstream>
#include <iostream>
using namespace cartmesh2d;using namespace cartmesh2d::fv;
namespace {
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F> void rejects(F action){bool failed=false;try{action();}catch(const std::exception&){failed=true;}require(failed,"invalid/failed implicit control accepted");}
}
int main(){try{
 const auto m=fv_test::rectangle(4,3,1,true);const IdealGas2D gas{1.4,1};const EulerTransport2D transport{.02,.02};
 std::vector<EulerBoundary2D> bc;for(std::size_t f=0;f<m.faces.size();++f)if(!m.faces[f].neighbour){EulerBoundary2D b;b.face=f;b.kind=EulerBoundaryKind2D::NoSlipWall;b.name="wall";b.thermalKind=HeatBoundaryKind2D::Temperature;b.thermalValue=1.1;bc.push_back(b);}
 EulerState2D initial;for(const auto& c:m.cells)initial.cells.push_back(eulerConservative2D({1+.02*c.centre.x,.03*c.centre.y,-.01*c.centre.x,1+.03*c.centre.y},gas));
 EulerStepper2D solver(m,bc,gas,transport);EulerStepControls2D controls;controls.order=2;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.integrator=EulerTimeIntegrator2D::Sdirk2;controls.maximumStep=.01;
 auto strict=initial,loose=initial;std::size_t converged=0;double maxGap=0;
 for(unsigned n=0;n<4;++n){controls.implicitIteration={};controls.nonlinearTolerance=2e-14;const auto reference=solver.advance(strict,controls);strict=reference.state;
  controls.nonlinearTolerance=1e-10;controls.implicitIteration=[&](const EulerImplicitIteration2D& e){if(e.converged){require(e.maximumScaledDefect<=controls.nonlinearTolerance,"unconverged stage claimed convergence");++converged;}};
  const auto before=loose;const auto step=solver.advance(loose,controls);require(loose.cells==before.cells&&loose.time==before.time,"advance mutated input");require(step.minimumDensity>0&&step.minimumPressure>0,"non-positive accepted output");require(step.maximumAcceptedStageDefect<=controls.nonlinearTolerance&&step.maximumStageOutputDefect<=8*controls.nonlinearTolerance,"stage budget/envelope violated");
  // Direct per-cell shared-flux RK identity, independent of scalar diagnostics.
  std::vector<EulerConservative2D> residual(m.cells.size());for(std::size_t f=0;f<m.faces.size();++f)for(std::size_t k=0;k<4;++k){residual[m.faces[f].owner][k]+=step.faceFlux[f][k];if(m.faces[f].neighbour)residual[*m.faces[f].neighbour][k]-=step.faceFlux[f][k];}
  for(std::size_t i=0;i<m.cells.size();++i)for(std::size_t k=0;k<4;++k){const double expected=before.cells[i][k]-step.step/m.cells[i].area*residual[i][k];require(std::abs(expected-step.state.cells[i][k])<=64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(expected)),"shared conservative RK update changed");maxGap=std::max(maxGap,std::abs(step.state.cells[i][k]-strict.cells[i][k]));}
  loose=step.state;
 }
 // Retry attempts can converge stage 1 before rejecting stage 2.
 require(converged>=8,"missing converged stages");
 // Unit-scale four-step empirical regression allowance, not a universal bound.
 require(maxGap<1e-8,"small-case algebraic error regression exceeded budget");
 controls.implicitIteration={};std::stringstream checkpoint;writeEulerCheckpoint2D(checkpoint,m,bc,gas,loose,"budget",transport);const auto restored=readEulerCheckpoint2D(checkpoint,m,bc,gas,"budget",transport);const auto a=solver.advance(loose,controls),b=solver.advance(restored,controls);require(a.state.cells==b.state.cells&&a.state.time==b.state.time,"restart changed same-control trajectory");
 for(double bad:{0.,-1.,1e-9,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){controls.nonlinearTolerance=bad;rejects([&]{(void)solver.advance(loose,controls);});}
 controls.nonlinearTolerance=1e-10;controls.maximumNewtonIterations=1;controls.maximumRetries=0;controls.maximumStep=.1;rejects([&]{(void)solver.advance(initial,controls);});require(initial.time==0&&initial.steps==0,"failed stage changed accepted input");
 controls.interrupted=[] {return true;};rejects([&]{(void)solver.advance(loose,controls);});require(loose.cells==restored.cells&&loose.time==restored.time,"cancelled stage changed accepted input");
 std::cout<<"max four-step algebraic field gap="<<maxGap<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
