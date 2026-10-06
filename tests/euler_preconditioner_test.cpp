#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include "cartmesh2d/fv/detail/EulerPreconditioner2D.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::detail;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::exception&){failed=true;}require(failed,"invalid implicit candidate accepted");}
LinearVector2D scales(const std::vector<EulerConservative2D>& cells,const IdealGas2D& gas) {
    LinearVector2D result;
    for(const auto& u:cells){const auto p=eulerPrimitive2D(u,gas);const double c=eulerSoundSpeed2D(p,gas);result.insert(result.end(),{p.density,p.density*c,p.density*c,u[3]});}
    return result;
}
std::vector<EulerBoundary2D> periodic(const FvMesh2D& m) {
    std::vector<EulerBoundary2D> result;
    for(std::size_t i=0;i<m.faces.size();++i)if(!m.faces[i].neighbour) {
        const auto& f=m.faces[i];EulerBoundary2D b;b.face=i;b.kind=EulerBoundaryKind2D::Periodic;
        b.name=f.areaVector.x!=0?"periodic-x":"periodic-y";
        for(std::size_t j=0;j<m.faces.size();++j)if(i!=j&&!m.faces[j].neighbour) {
            const auto& g=m.faces[j];
            if(dot(f.areaVector,g.areaVector)<0&&
               ((f.areaVector.x!=0&&f.centre.y==g.centre.y)||(f.areaVector.y!=0&&f.centre.x==g.centre.x)))b.partner=j;
        }
        require(b.partner.has_value(),"periodic fixture lacks partner");result.push_back(b);
    }
    return result;
}
std::vector<const EulerBoundary2D*> lookup(const FvMesh2D& m,const std::vector<EulerBoundary2D>& b) {
    std::vector<const EulerBoundary2D*> result(m.faces.size());for(const auto& x:b)result[x.face]=&x;return result;
}
EulerConservative2D noGhost(const EulerConservative2D&,const EulerBoundary2D&,Vector2D){throw std::runtime_error("periodic face incorrectly requested a ghost");}
double jacobianError=0,periodicMatrixError=0,sharedMatrixBalance=0,acceptedFieldGap=0,acceptedFluxBalance=0;
void fluxDerivative() {
    // Central differences of the native physical flux. Scaling has units of
    // flux, so density/momentum/energy and SI/nondimensional gas are comparable.
    for(const auto gas:{IdealGas2D{1.4,1},IdealGas2D{}})
        for(const auto primitive:{EulerPrimitive2D{1,.25,-.13,1},EulerPrimitive2D{1.17,65,-12,101325}})
            for(const auto area:{Vector2D{.2,.7},Vector2D{-.7,.2}}) {
                const auto u=eulerConservative2D(primitive,gas);const auto scale=scales({u},gas);
                const auto a=eulerNormalFluxJacobian2D(u,area,gas);
                const double speed=std::hypot(primitive.u,primitive.v)+eulerSoundSpeed2D(primitive,gas),length=std::hypot(area.x,area.y);
                for(std::size_t column=0;column<4;++column) {
                    const double delta=std::cbrt(std::numeric_limits<double>::epsilon())*scale[column];auto plus=u,minus=u;
                    plus[column]+=delta;minus[column]-=delta;
                    const auto fp=eulerFaceFlux2D(plus,plus,area,gas,EulerFluxScheme2D::Rusanov).integratedFlux;
                    const auto fm=eulerFaceFlux2D(minus,minus,area,gas,EulerFluxScheme2D::Rusanov).integratedFlux;
                    for(std::size_t row=0;row<4;++row)jacobianError=std::max(jacobianError,std::abs(a[row][column]-(fp[row]-fm[row])/(2*delta))*scale[column]/(scale[row]*speed*length));
                }
            }
    // O(eps^(2/3)) central-difference uncertainty, not a physical accuracy gate.
    require(jacobianError<512*std::pow(std::numeric_limits<double>::epsilon(),2./3),"physical Euler flux derivative mismatch");
}
void periodicMatrix() {
    const auto m=fv_test::rectangle(3,2,2);const auto bc=periodic(m);const auto map=lookup(m,bc);const IdealGas2D gas{1.4,1};
    validateEulerBoundaries2D(m,bc,gas);
    EulerState2D state;state.cells.assign(m.cells.size(),eulerConservative2D({1,.25,.13,1},gas));
    const auto scale=scales(state.cells,gas);const double step=.04;
    const EulerFrozenPreconditioner2D matrix(m,map,gas,{},state.cells,step,scale,noGhost);
    LinearVector2D direction(scale.size());for(std::size_t j=0;j<direction.size();++j)direction[j]=std::sin(double(j)+.3);
    const auto action=matrix.applyMatrix(direction);
    // At a uniform state the derivative of the Rusanov speed multiplies a zero
    // jump; this special case has an exact frozen-speed Jacobian. Read the true
    // native first-order face residual, without another PDE implementation.
    const double epsilon=std::cbrt(std::numeric_limits<double>::epsilon());auto plus=state,minus=state;
    for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k){plus.cells[c][k]+=epsilon*scale[4*c+k]*direction[4*c+k];minus.cells[c][k]-=epsilon*scale[4*c+k]*direction[4*c+k];}
    EulerStepControls2D control;control.maximumStep=1e-6;
    const auto fp=advanceEuler2D(m,bc,gas,plus,control).faceFlux,fm=advanceEuler2D(m,bc,gas,minus,control).faceFlux;
    std::vector<EulerConservative2D> derivative(m.cells.size());
    for(std::size_t f=0;f<m.faces.size();++f)for(std::size_t k=0;k<4;++k) {
        const double d=(fp[f][k]-fm[f][k])/(2*epsilon);derivative[m.faces[f].owner][k]+=d;
        if(m.faces[f].neighbour)derivative[*m.faces[f].neighbour][k]-=d;
    }
    for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k){const auto j=4*c+k;periodicMatrixError=std::max(periodicMatrixError,std::abs(action[j]-direction[j]-step*derivative[c][k]/(m.cells[c].area*scale[j])));}
    // The spectral-radius max/abs is nonsmooth at equal states; its central
    // difference contributes O(epsilon) through the perturbed state jump.
    require(periodicMatrixError<8*epsilon,"assembled periodic inviscid Jacobian mismatch");
    const EulerFrozenPreconditioner2D transport(m,map,gas,{.2,.1},state.cells,step,scale,noGhost);
    const auto coupled=transport.applyMatrix(direction);
    for(std::size_t k=0;k<4;++k) {
        long double sum=0,absolute=0;
        for(std::size_t c=0;c<m.cells.size();++c){const double value=m.cells[c].area*scale[4*c+k]*(coupled[4*c+k]-direction[4*c+k]);sum+=value;absolute+=std::abs(value);}
        sharedMatrixBalance=std::max(sharedMatrixBalance,static_cast<double>(std::abs(sum)/absolute));
    }
    require(sharedMatrixBalance<256*std::numeric_limits<double>::epsilon(),"preconditioner lost shared periodic conservation");
    const auto z=transport.apply(direction);for(double x:z)require(std::isfinite(x),"ILU application nonfinite");
    rejects([&]{(void)transport.apply({1});});
}
void inspect(const FvMesh2D& m,const EulerStepResult2D& result,const IdealGas2D& gas) {
    require(result.minimumDensity>0&&result.minimumPressure>0,"preconditioned step lost positivity");
    const auto scale=scales(result.previousCells,gas);std::vector<EulerConservative2D> residual(m.cells.size());
    for(std::size_t f=0;f<m.faces.size();++f)for(std::size_t k=0;k<4;++k){residual[m.faces[f].owner][k]+=result.faceFlux[f][k];if(m.faces[f].neighbour)residual[*m.faces[f].neighbour][k]-=result.faceFlux[f][k];}
    for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k)acceptedFluxBalance=std::max(acceptedFluxBalance,std::abs(result.state.cells[c][k]-result.previousCells[c][k]+result.step*residual[c][k]/m.cells[c].area)/scale[4*c+k]);
    require(acceptedFluxBalance<256*std::numeric_limits<double>::epsilon(),"preconditioned output differs from exported conservative flux");
}
void paired(const char* name,const FvMesh2D& m,const std::vector<EulerBoundary2D>& bc,
    const IdealGas2D& gas,const EulerTransport2D& physics,const EulerState2D& initial,EulerStepControls2D control) {
    const EulerStepper2D solver(m,bc,gas,physics,control.wallGradient,control.diffusionScheme);const auto scale=scales(initial.cells,gas);
    control.integrator=EulerTimeIntegrator2D::Sdirk2;control.maximumRetries=0;
    auto diagonal=initial,coupled=initial;std::size_t diagonalIterations=0,coupledIterations=0,diagonalResiduals=0,coupledResiduals=0;
    for(unsigned count=0;count<4;++count) {
        control.implicitPreconditioner=EulerImplicitPreconditioner2D::Diagonal;const auto a=solver.advance(diagonal,control);inspect(m,a,gas);
        control.implicitPreconditioner=EulerImplicitPreconditioner2D::FrozenFluxIlu0;const auto b=solver.advance(coupled,control);inspect(m,b,gas);
        require(a.step==b.step&&a.rejectedCandidates==0&&b.rejectedCandidates==0,"preconditioner changed accepted clock");
        diagonalIterations+=a.linearIterations;coupledIterations+=b.linearIterations;diagonalResiduals+=a.spatialEvaluations;coupledResiduals+=b.spatialEvaluations;
        if(count==1) {
            std::stringstream data;writeEulerCheckpoint2D(data,m,bc,gas,coupled,"preconditioner",physics);
            const auto restored=readEulerCheckpoint2D(data,m,bc,gas,"preconditioner",physics);
            const auto resumed=EulerStepper2D(m,bc,gas,physics,control.wallGradient,control.diffusionScheme).advance(restored,control);
            require(resumed.state.cells==b.state.cells&&resumed.state.time==b.state.time,"ILU restart has hidden history");
        }
        for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k)acceptedFieldGap=std::max(acceptedFieldGap,std::abs(a.state.cells[c][k]-b.state.cells[c][k])/scale[4*c+k]);
        diagonal=a.state;coupled=b.state;
    }
    // Compare the same converged nonlinear equations. This budget allows the
    // existing two-stage tolerance to accumulate through four physical steps.
    require(acceptedFieldGap<400*control.nonlinearTolerance,"preconditioner changed converged physical update");
    const auto before=coupled;control.maximumNewtonIterations=1;control.maximumStep*=100;
    rejects([&]{(void)solver.advance(coupled,control);});require(coupled.cells==before.cells&&coupled.time==before.time,"failed ILU stage overwrote last accepted state");
    control.interrupted=[] {return true;};rejects([&]{(void)solver.advance(coupled,control);});
    std::cout<<"{\"case\":\""<<name<<"\",\"diagonalLinearIterations\":"<<diagonalIterations<<",\"iluLinearIterations\":"<<coupledIterations<<",\"diagonalResiduals\":"<<diagonalResiduals<<",\"iluResiduals\":"<<coupledResiduals<<"}\n";
}
void coupledStages() {
    const auto m=fv_test::rectangle(6,3,2,true);const auto bc=periodic(m);const IdealGas2D gas{1.4,1};EulerState2D initial;
    for(const auto& c:m.cells){const double x=std::numbers::pi*c.centre.x,y=2*std::numbers::pi*c.centre.y;initial.cells.push_back(eulerConservative2D({1+.003*std::sin(x),.02*std::sin(x),.01*std::cos(y),1+.005*std::cos(x)},gas));}
    EulerStepControls2D control;control.maximumStep=.04;control.order=2;control.fluxScheme=EulerFluxScheme2D::Hllc;
    paired("periodic-coupled",m,bc,gas,{.02,.02},initial,control);
    auto walls=bc;for(auto& b:walls){b.kind=EulerBoundaryKind2D::NoSlipWall;b.partner.reset();b.name="wall";b.thermalKind=HeatBoundaryKind2D::Temperature;b.thermalValue=1.1;}
    control.diffusionScheme=EulerDiffusionScheme2D::Hybrid;control.maximumStep=.01;
    paired("hybrid-cavity",m,walls,gas,{.02,.02},initial,control);
}
void openChannel() {
    auto m=fv_test::rectangle(3,2,2,true);for(auto& c:m.cells){c.centre.x*=1e-4;c.centre.y*=1e-4;c.area*=1e-8;}for(auto& f:m.faces){f.centre.x*=1e-4;f.centre.y*=1e-4;f.areaVector.x*=1e-4;f.areaVector.y*=1e-4;f.correction.x*=1e-4;f.correction.y*=1e-4;}
    const IdealGas2D gas;const EulerPrimitive2D q{101325/(300*gas.gasConstant),60,0,101325};const EulerTransport2D physics{.025759,1.846e-5};std::vector<EulerBoundary2D> bc;
    for(std::size_t f=0;f<m.faces.size();++f)if(!m.faces[f].neighbour) {
        EulerBoundary2D b;b.face=f;b.name="wall";b.kind=EulerBoundaryKind2D::NoSlipWall;b.thermalKind=HeatBoundaryKind2D::Temperature;b.thermalValue=330;b.reference=q;
        if(m.faces[f].areaVector.x!=0){b.kind=m.faces[f].areaVector.x<0?EulerBoundaryKind2D::TotalInlet:EulerBoundaryKind2D::PressureOutlet;b.name=m.faces[f].areaVector.x<0?"inlet":"outlet";b.thermalKind=HeatBoundaryKind2D::Insulated;b.thermalValue=0;}bc.push_back(b);
    }
    EulerState2D state;state.cells.assign(m.cells.size(),eulerConservative2D(q,gas));EulerStepControls2D control;control.maximumStep=2e-8;control.order=2;control.fluxScheme=EulerFluxScheme2D::Hllc;
    paired("SI-total-inlet-heated-channel",m,bc,gas,physics,state,control);
}
}
int main() {
    try {
        std::cout<<std::setprecision(17);fluxDerivative();periodicMatrix();coupledStages();openChannel();
        std::cout<<"{\"physicalJacobianScaledError\":"<<jacobianError<<",\"periodicMatrixScaledError\":"<<periodicMatrixError<<",\"sharedMatrixBalance\":"<<sharedMatrixBalance<<",\"acceptedFieldGap\":"<<acceptedFieldGap<<",\"acceptedFluxBalance\":"<<acceptedFluxBalance<<"}\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
