#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include "cartmesh2d/fv/detail/FlowMaterial2D.hpp"
#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "../src/fv/FlowEquation2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}
FvMesh2D rectangle(std::size_t nx,std::size_t ny,double length=1) {
    std::vector<double> x(nx+1),y(ny+1);
    for(std::size_t i=0;i<=nx;++i)x[i]=length*static_cast<double>(i)/static_cast<double>(nx);
    for(std::size_t i=0;i<=ny;++i)y[i]=static_cast<double>(i)/static_cast<double>(ny);
    return makeFvMesh2D(makeRectilinearMesh2D(x,y));
}
void strict(const FlowResult2D& result,const FlowControls2D& c) {
    require(result.converged,result.failureReason.empty()?"coupled did not converge":result.failureReason);
    const auto& h=result.history.back();
    require(h.strictLinearStep && detail::strictFlowConverged2D(h.iteration,h.momentumResidual,h.velocityChange,
        h.pressureChange,h.continuity,result.globalRelativeImbalance,c.tolerance,true),"coupled bypassed strict gate");
    require(h.fluxConsistency<std::min(c.tolerance,1e-8),"retained flux does not satisfy refreshed Rhie-Chow equation");
    require(result.outletBackflowFaces==0 && result.outletInflow==0,"coupled accepted outlet backflow");
    require(result.performance.coupledSolves>0 && result.performance.coupledIterations>0,"missing total Krylov cost");
}
double difference(const FlowResult2D& a,const FlowResult2D& b) {
    double error=0;
    for(std::size_t i=0;i<a.u.size();++i)error=std::max({error,std::abs(a.u[i]-b.u[i]),std::abs(a.v[i]-b.v[i]),std::abs(a.p[i]-b.p[i])});
    for(std::size_t i=0;i<a.flux.size();++i)error=std::max(error,std::abs(a.flux[i]-b.flux[i]));
    return error;
}
void sameFixedPoint(const FlowResult2D& a,const FlowResult2D& b,const std::string& label) {
    // Reuse the established same-equation iteration-error tolerance, not a
    // spatial-accuracy gate. Include the conservative flux and final loads.
    require(difference(a,b)<2e-6,label+" changed fields or conservative face flux");
    for(const auto& pair:{std::pair{a.forceX,b.forceX},std::pair{a.forceY,b.forceY},
        std::pair{a.wallForceX,b.wallForceX},std::pair{a.wallForceY,b.wallForceY},
        std::pair{a.pressureForceX,b.pressureForceX},std::pair{a.pressureForceY,b.pressureForceY}})
        require(std::abs(pair.first-pair.second)<2e-6,label+" changed physical loads");
    require(a.namedWallLoads.size()==b.namedWallLoads.size(),label+" changed named wall groups");
    for(std::size_t i=0;i<a.namedWallLoads.size();++i) {
        const auto& x=a.namedWallLoads[i];const auto& y=b.namedWallLoads[i];
        require(x.name==y.name && x.faces==y.faces && x.length==y.length,label+" changed named wall geometry");
        require(std::abs(x.pressureTorque+x.viscousTorque-y.pressureTorque-y.viscousTorque)<2e-6,
            label+" changed named wall torque");
    }
}
void consistentRelaxations(const FvMesh2D& mesh,FlowControls2D c,const FlowResult2D& reference,
    bool identicalDiagonals) {
    c.coupling=FlowCoupling2D::SimpleConsistent;c.maxIterations=1000;c.adaptiveLinear=true;
    for(const auto& relaxations:{std::pair{.6,.25},std::pair{1.,1.},std::pair{.45,.6},
        std::pair{1.,.25},std::pair{.6,1.}}) {
        c.velocityRelaxation=relaxations.first;c.pressureRelaxation=relaxations.second;
        const auto result=solveIncompressible2D(mesh,c);strict(result,c);
        sameFixedPoint(result,reference,"SIMPLE full response");
        if(identicalDiagonals && relaxations.first==1 && relaxations.second==1) {
            require(difference(result,reference)<2e-12,"unit relaxation did not recover the coupled response");
            require(result.history.size()==reference.history.size(),"unit relaxation changed coupled outer iteration count");
        }
        std::cout<<"consistent alphaU="<<relaxations.first<<" alphaP="<<relaxations.second
                 <<" iterations="<<result.history.size()<<" difference="<<difference(result,reference)<<'\n';
    }
}
FlowControls2D custom(const FvMesh2D&,bool);
void asymmetricSlipResponse() {
    const auto mesh=rectangle(7,5,2);auto c=custom(mesh,false);
    c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
    for(auto& condition:c.boundaryConditions) {
        const auto& face=mesh.faces[condition.face];
        if(face.areaVector.y>0)condition={condition.face,FlowBoundaryKind2D::Symmetry,{},0,"top-slip"};
        if(condition.kind==FlowBoundaryKind2D::VelocityInlet)
            condition.velocity={.8+.3*face.centre.y,.08*std::sin(std::acos(-1.)*face.centre.y)};
    }
    c.coupling=FlowCoupling2D::Simple;c.maxIterations=1500;
    const auto segregated=solveIncompressible2D(mesh,c);
    require(segregated.converged,"asymmetric slip SIMPLE reference failed");
    c.coupling=FlowCoupling2D::Coupled;c.maxIterations=500;c.adaptiveLinear=true;
    const auto reference=solveIncompressible2D(mesh,c);strict(reference,c);
    sameFixedPoint(reference,segregated,"asymmetric slip reference");
    // Inspect the production compact assembly solely to certify that this
    // fixture actually exercises unequal U/V diagonals, rather than assuming
    // a rectangular mesh implies the required response asymmetry.
    const auto boundary=solver_detail::boundaries(mesh,c);
    solver_detail::FlowEquation2D equation(mesh,boundary.fixedP);
    std::vector<std::pair<std::size_t,std::size_t>> connections;
    for(const auto& face:mesh.faces)if(face.neighbour)connections.emplace_back(face.owner,*face.neighbour);
    const detail::SparsePattern2D pattern(mesh.cells.size(),connections);
    solver_detail::System u(pattern),v(pattern);
    const auto reconstruction=equation.reconstruct(reference,c,boundary,boundary.p,reference.flux);
    equation.assembleMomentum(u,v,reference,c,boundary,reconstruction,reference.flux,reference.sourceIntegrals);
    std::size_t unequal=0;
    for(std::size_t i=0;i<mesh.cells.size();++i)if(u.diag[i]!=v.diag[i])++unequal;
    require(unequal>0,"asymmetric slip fixture did not exercise unequal momentum diagonals");
    std::cout<<"asymmetric slip unequal momentum diagonal rows="<<unequal<<'\n';
    // At alpha=1 the original common-diagonal SIMPLE mechanism still adds
    // D-aU/V on slip rows. It preserves the fixed point, but does not promise
    // the identical coupled iteration trajectory in this fixture.
    consistentRelaxations(mesh,c,reference,false);
}
FlowResult2D compare(const FvMesh2D& mesh,FlowControls2D c) {
    c.coupling=FlowCoupling2D::Simple;c.maxIterations=1500;
    const auto reference=solveIncompressible2D(mesh,c);
    require(reference.converged,"SIMPLE comparison failed: "+c.scenario);
    c.coupling=FlowCoupling2D::Coupled;c.maxIterations=250;c.adaptiveLinear=true;
    const auto result=solveIncompressible2D(mesh,c);strict(result,c);
    const double error=difference(result,reference);
    // Same discrete equations on the same mesh: this is an iteration-error
    // comparison, not a new spatial-accuracy threshold.
    require(error<2e-6,"coupled changed the discrete fixed point: "+c.scenario);
    require(result.sourceIntegrals.size()==reference.sourceIntegrals.size(),"source integration size differs from SIMPLE");
    for(std::size_t i=0;i<result.sourceIntegrals.size();++i)
        require(result.sourceIntegrals[i].x==reference.sourceIntegrals[i].x && result.sourceIntegrals[i].y==reference.sourceIntegrals[i].y,"source integration differs from SIMPLE");
    require(std::abs(result.wallForceX-reference.wallForceX)<2e-6,"shared force postprocessing differs");
    std::cout<<c.scenario<<" SIMPLE="<<reference.history.size()<<" coupled="<<result.history.size()
        <<" Krylov="<<result.performance.coupledIterations<<" difference="<<error<<'\n';
    return result;
}
template<class Function> void rejects(Function function,const std::string& label) {
    bool rejected=false;try{function();}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"unsupported controls were accepted: "+label);
}
FlowControls2D custom(const FvMesh2D& mesh,bool fixedOutlet) {
    FlowControls2D c;c.scenario="custom";c.nu=.1;c.tolerance=1e-9;
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];if(f.neighbour)continue;
        FlowBoundaryCondition2D b{id,FlowBoundaryKind2D::Wall,{},0,"wall"};
        if(f.areaVector.x<0 || (fixedOutlet && f.areaVector.x>0))
            b={id,f.areaVector.x<0?FlowBoundaryKind2D::VelocityInlet:FlowBoundaryKind2D::VelocityOutlet,
                {4*f.centre.y*(1-f.centre.y),0},0,f.areaVector.x<0?"inlet":"outlet"};
        else if(f.areaVector.x>0)b={id,FlowBoundaryKind2D::PressureOutlet,{},.4,"outlet"};
        c.boundaryConditions.push_back(b);
    }
    return c;
}
}
int main() {
    const auto channel=rectangle(12,6,2),box=rectangle(10,10);
    FlowControls2D c;c.scenario="channel";c.nu=.1;c.tolerance=1e-9;
    for(const auto scheme:{ConvectionScheme2D::Upwind,ConvectionScheme2D::LimitedLinearUpwind,ConvectionScheme2D::FaceLimitedLinearUpwind}) {
        c.convection=scheme;c.viscousStress=scheme==ConvectionScheme2D::Upwind?ViscousStress2D::Laplacian:ViscousStress2D::Symmetric;
        (void)compare(channel,c);
    }
    c.scenario="manufactured";c.manufacturedPressureSlope=.2;
    const auto manufactured=compare(box,c);
    auto consistentManufactured=c;consistentManufactured.coupling=FlowCoupling2D::SimpleConsistent;
    consistentManufactured.maxIterations=1000;consistentManufactured.adaptiveLinear=true;
    const auto manufacturedResponse=solveIncompressible2D(box,consistentManufactured);
    strict(manufacturedResponse,consistentManufactured);sameFixedPoint(manufacturedResponse,manufactured,"manufactured full response");
    c.scenario="cavity";c.manufacturedPressureSlope=0;(void)compare(box,c);
    auto physical=custom(channel,false);const auto accepted=compare(channel,physical);
    consistentRelaxations(channel,physical,accepted,true);
    asymmetricSlipResponse();
    auto prescribed=custom(channel,true);(void)compare(channel,prescribed);
    prescribed.momentumInertia=0;(void)compare(channel,prescribed);

    auto tight=physical;tight.coupling=FlowCoupling2D::Coupled;tight.maxIterations=250;
    const auto allStrict=solveIncompressible2D(channel,tight);strict(allStrict,tight);
    require(difference(allStrict,accepted)<2e-8,"adaptive forcing changed the coupled fixed point");
    require(std::all_of(allStrict.history.begin(),allStrict.history.end(),[](const auto& h){return h.strictLinearStep;}),"strict policy loosened a coupled linear solve");
    for(const auto method:{PressurePreconditioner2D::Jacobi,PressurePreconditioner2D::Aggregation}) {
        tight.pressurePreconditioner=method;
        const auto preconditioned=solveIncompressible2D(channel,tight);strict(preconditioned,tight);
        require(difference(preconditioned,allStrict)<2e-8,"pressure preconditioner changed the coupled equations");
    }

    const FlowInitialGuess2D guess{accepted.u,accepted.v,accepted.p,accepted.flux};
    for(const auto coupling:{FlowCoupling2D::Coupled,FlowCoupling2D::SimpleConsistent}) {
        physical.coupling=coupling;physical.adaptiveLinear=true;
        physical.maxIterations=9;const auto tooShort=solveIncompressibleFromGuess2D(channel,physical,guess);
        require(!tooShort.converged && tooShort.history.back().iteration==9,"minimum iteration gate bypassed by exact initial guess");
        physical.maxIterations=30;physical.stopRequested=[] {return true;};
        const auto cancelled=solveIncompressibleFromGuess2D(channel,physical,guess);
        require(cancelled.stopped && !cancelled.converged && cancelled.history.back().iteration==1,"cancellation failed to preserve accepted iterate");
        physical.stopRequested={};
        auto invalidGuess=guess;
        for(std::size_t id=0;id<channel.faces.size();++id)if(!channel.faces[id].neighbour && channel.faces[id].areaVector.x>0)invalidGuess.flux[id]=-1;
        rejects([&]{(void)solveIncompressibleFromGuess2D(channel,physical,invalidGuess);},"outlet backflow");
        invalidGuess=guess;
        for(std::size_t id=0;id<channel.faces.size();++id)if(!channel.faces[id].neighbour && channel.faces[id].areaVector.x<0)invalidGuess.flux[id]=0;
        rejects([&]{(void)solveIncompressibleFromGuess2D(channel,physical,invalidGuess);},"prescribed inlet flux");
        for(int kind=0;kind<6;++kind) {
            auto bad=physical;
            if(kind==0)bad.convergence=FlowConvergence2D::Engineering;
            if(kind==1)bad.steadyAcceleration=SteadyAcceleration2D::Anderson;
            if(kind==2)bad.outletBackflow=OutletBackflow2D::NormalInlet;
            if(kind==3)bad.faceViscosity.assign(channel.faces.size(),bad.nu);
            if(kind==4)for(auto& b:bad.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOutlet)b.kind=FlowBoundaryKind2D::PressureOpening;
            if(kind==5){bad.scenario="manufactured";bad.boundaryConditions.clear();bad.manufacturedViscositySlope=.1;}
            rejects([&]{(void)solveIncompressible2D(kind==5?box:channel,bad);},std::to_string(kind));
        }
        const auto initial=initialIncompressibleState2D(channel,physical);
        rejects([&]{(void)advanceIncompressible2D(channel,physical,initial,.1);},"transient");
        rejects([&]{(void)detail::solveMaterialFlow2D(channel,physical,[](const auto&,const auto&){return detail::MaterialState2D{};});},"material/turbulence");
    }

    // An unattainable strict tolerance must fail without accepting a Krylov
    // candidate. Compare the returned state with an independently rerun prefix
    // ending at exactly its last accepted outer iteration.
    const auto small=rectangle(4,4);auto demanding=custom(small,false);
    for(const auto coupling:{FlowCoupling2D::Coupled,FlowCoupling2D::SimpleConsistent}) {
        demanding.coupling=coupling;demanding.adaptiveLinear=true;demanding.tolerance=1e-30;demanding.maxIterations=300;
        const auto failed=solveIncompressible2D(small,demanding);
        require(!failed.converged && !failed.failureReason.empty(),"unattainable tolerance was silently accepted");
        require(failed.history.back().iteration>0,"failure test did not reach an accepted prefix");
        demanding.maxIterations=failed.history.back().iteration;
        const auto prefix=solveIncompressible2D(small,demanding);
        require(prefix.failureReason.empty() && difference(failed,prefix)==0,"failed trial overwrote last accepted fields or flux");
        demanding.adaptiveLinear=false;
        const FlowInitialGuess2D failedGuess{prefix.u,prefix.v,prefix.p,prefix.flux};
        const auto firstFailure=solveIncompressibleFromGuess2D(small,demanding,failedGuess);
        require(!firstFailure.converged && !firstFailure.failureReason.empty() && firstFailure.history.back().iteration==0,
            "failed first solve is missing its iteration-zero diagnostic");
        require(difference(firstFailure,prefix)==0,"failed first solve changed initialization");
        std::cout<<"failure preserved iteration "<<failed.history.back().iteration<<": "<<failed.failureReason<<'\n';
    }
}
