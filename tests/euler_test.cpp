#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
EulerStepControls2D activeControls;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F> void rejects(F fn){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}require(rejected,"invalid Euler input accepted");}
using fv_test::rectangle;

std::vector<EulerBoundary2D> boundaries(const FvMesh2D& mesh,EulerBoundaryKind2D kind,EulerPrimitive2D reference={}) {
    std::vector<EulerBoundary2D> result;
    for(std::size_t i=0;i<mesh.faces.size();++i)if(!mesh.faces[i].neighbour)
        result.push_back({i,kind,reference,{},"boundary"});
    return result;
}
EulerState2D constant(const FvMesh2D& mesh,EulerPrimitive2D q) {
    return {0,0,std::vector<EulerConservative2D>(mesh.cells.size(),eulerConservative2D(q))};
}
void same(const EulerState2D& a,const EulerState2D& b,double tolerance=1e-12) {
    require(a.cells.size()==b.cells.size(),"state size differs");
    for(std::size_t i=0;i<a.cells.size();++i)for(std::size_t k=0;k<4;++k)
        require(std::abs(a.cells[i][k]-b.cells[i][k])<tolerance,"conserved state differs");
}
void pressureOutletChecks() {
    const IdealGas2D gas{1.4,1};
    const EulerPrimitive2D inside{1.3,.4,-.2,1.2};
    const auto outlet=eulerPressureOutletState2D(inside,1.1,{3,0},gas);
    const double ci=eulerSoundSpeed2D(inside,gas),cb=eulerSoundSpeed2D(outlet,gas);
    // 128 eps after nondimensional normalization checks characteristic algebra,
    // not a physical accuracy threshold for the computed flow.
    const double roundoff=128*std::numeric_limits<double>::epsilon();
    require(outlet.pressure==1.1&&outlet.v==inside.v,"outlet did not preserve pressure/tangential velocity");
    require(std::abs((outlet.u+2*cb/(gas.gamma-1))/(inside.u+2*ci/(gas.gamma-1))-1)<roundoff,"outgoing outlet invariant changed");
    require(std::abs(outlet.pressure/inside.pressure*std::pow(inside.density/outlet.density,gas.gamma)-1)<roundoff,"outlet entropy changed");
    const double angle=.73,cosine=std::cos(angle),sine=std::sin(angle);
    const auto turn=[&](const EulerPrimitive2D& q){return EulerPrimitive2D{q.density,q.u*cosine-q.v*sine,q.u*sine+q.v*cosine,q.pressure};};
    const auto turned=eulerPressureOutletState2D(turn(inside),1.1,{cosine,sine},gas),expected=turn(outlet);
    require(std::hypot(turned.u-expected.u,turned.v-expected.v)<roundoff,"outlet depends on global coordinate axes");
    const EulerPrimitive2D supersonic{1,2,0,1};
    require(eulerPressureOutletState2D(supersonic,.9,{1,0},gas).pressure==1,"supersonic outlet imposed an incoming pressure");
    rejects([&]{(void)eulerPressureOutletState2D(inside,0,{1,0},gas);});
    rejects([&]{(void)eulerPressureOutletState2D(inside,1.1,{0,0},gas);});
    rejects([&]{(void)eulerPressureOutletState2D(inside,2,{1,0},gas);});
    rejects([&]{(void)eulerPressureOutletState2D(inside,.01,{1,0},gas);});
    rejects([&]{(void)eulerPressureOutletState2D({1,-.1,0,1},1,{1,0},gas);});

    const auto mesh=rectangle(12,4,2,true);auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    const EulerPrimitive2D uniform{1,.3,0,1};
    for(auto& b:bc)if(mesh.faces[b.face].areaVector.x!=0) {
        const bool right=mesh.faces[b.face].areaVector.x>0;
        b.kind=right?EulerBoundaryKind2D::PressureOutlet:EulerBoundaryKind2D::Farfield;
        b.name=right?"outlet":"inlet";b.reference=uniform;
    }
    const EulerStepper2D solver(mesh,bc,gas);auto state=constant(mesh,uniform);const auto initial=state;
    for(unsigned i=0;i<40;++i)state=solver.advance(state,activeControls).state;
    same(initial,state);
    auto driven=bc;for(auto& b:driven)if(b.kind==EulerBoundaryKind2D::PressureOutlet)b.reference.pressure=.99;
    auto firstOrder=activeControls;firstOrder.order=1;
    const auto drivenStep=EulerStepper2D(mesh,driven,gas).advance(initial,firstOrder);
    const auto trace=eulerPressureOutletState2D(uniform,.99,{1,0},gas);
    for(const auto& b:driven)if(b.kind==EulerBoundaryKind2D::PressureOutlet) {
        const auto& face=mesh.faces[b.face];const auto& flux=drivenStep.faceFlux[b.face];
        const double length=std::hypot(face.areaVector.x,face.areaVector.y);
        require(std::abs(flux[0]/length-trace.density*trace.u)<roundoff,"pressure outlet mass flux does not use its characteristic trace");
        require(std::abs((flux[1]-flux[0]*trace.u)/length-.99)<roundoff,"Riemann flux weakened the requested outlet pressure");
    }
    std::stringstream checkpoint;writeEulerCheckpoint2D(checkpoint,mesh,bc,gas,state,"pressure-outlet");
    const auto saved=checkpoint.str();const auto restored=readEulerCheckpoint2D(checkpoint,mesh,bc,gas,"pressure-outlet");
    require(restored.cells==state.cells,"pressure outlet checkpoint differs");
    auto changed=bc;for(auto& b:changed)if(b.kind==EulerBoundaryKind2D::PressureOutlet)b.reference.pressure=.99;
    rejects([&]{std::istringstream in(saved);(void)readEulerCheckpoint2D(in,mesh,changed,gas,"pressure-outlet");});
    // A runtime regime failure must preserve the latest accepted state.
    auto backflow=state;for(auto& q:backflow.cells)q[1]=-q[1];const auto before=backflow;
    rejects([&]{(void)solver.advance(backflow,activeControls);});
    require(backflow.cells==before.cells&&backflow.time==before.time&&backflow.steps==before.steps,"outlet failure overwrote accepted state");
}
void restAndFreeStream() {
    const auto mesh=rectangle(12,6,2,true);
    for(auto kind:{EulerBoundaryKind2D::SlipWall,EulerBoundaryKind2D::Farfield,EulerBoundaryKind2D::Transmissive}) {
        const EulerPrimitive2D q{1.3,kind==EulerBoundaryKind2D::SlipWall?0:.7,kind==EulerBoundaryKind2D::SlipWall?0:-.2,2};
        const auto bc=boundaries(mesh,kind,q);const auto initial=constant(mesh,q);auto state=initial;
        for(int step=0;step<40;++step) {
            const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);
            require(next.acousticCourant<=.400000000000001&&next.maximumCellBalanceError<1e-14,"CFL / conservation failure");
            state=next.state;
        }
        same(initial,state);
    }
    for(double speed:{-5.,5.}) {
        const EulerPrimitive2D q{1,speed,.3,1};const auto state=constant(mesh,q);
        same(state,advanceEuler2D(mesh,boundaries(mesh,EulerBoundaryKind2D::Farfield,q),{},state,activeControls).state);
    }
}
std::vector<EulerBoundary2D> periodic(const FvMesh2D& mesh) {
    auto bc=boundaries(mesh,EulerBoundaryKind2D::Periodic);
    for(auto& b:bc) {
        const auto& f=mesh.faces[b.face];const bool x=std::abs(f.areaVector.x)>std::abs(f.areaVector.y);b.name=x?"periodic-x":"periodic-y";
        for(const auto& other:bc) {
            const auto& g=mesh.faces[other.face];
            if(other.face!=b.face&&std::hypot(f.areaVector.x+g.areaVector.x,f.areaVector.y+g.areaVector.y)<1e-12&&
               std::abs(x?f.centre.y-g.centre.y:f.centre.x-g.centre.x)<1e-12){b.partner=other.face;break;}
        }
    }
    return bc;
}
void periodicConservation() {
    const auto mesh=rectangle(20,10,2,true);const auto bc=periodic(mesh);auto state=constant(mesh,{1,.3,-.2,1});
    const double pi=std::acos(-1.);
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto c=mesh.cells[i].centre;
        state.cells[i]=eulerConservative2D({1+.2*std::sin(pi*c.x)*std::cos(2*pi*c.y),.3+.1*std::cos(pi*c.x),-.2,1});
    }
    EulerConservative2D before{},after{};
    for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)before[k]+=mesh.cells[i].area*state.cells[i][k];
    for(int i=0;i<80;++i) {
        const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);
        for(const auto& b:bc)for(std::size_t k=0;k<4;++k)require(next.faceFlux[b.face][k]==-next.faceFlux[*b.partner][k],"periodic action-reaction differs");
        state=next.state;
    }
    for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)after[k]+=mesh.cells[i].area*state.cells[i][k];
    for(std::size_t k=0;k<4;++k)require(std::abs(after[k]-before[k])<2e-14,"periodic total conservation failed");
    auto bad=bc;bad[0].partner=bad[0].face;rejects([&]{validateEulerBoundaries2D(mesh,bad);});
    bad=bc;bad[0].name="unpaired";rejects([&]{validateEulerBoundaries2D(mesh,bad);});
}
FvMesh2D rotated(FvMesh2D mesh,double angle) {
    const auto vector=[&](Vector2D p){return Vector2D{p.x*std::cos(angle)-p.y*std::sin(angle),p.x*std::sin(angle)+p.y*std::cos(angle)};};
    for(auto& c:mesh.cells){const auto r=vector({c.centre.x,c.centre.y});c.centre={r.x,r.y};}
    for(auto& f:mesh.faces){const auto r=vector({f.centre.x,f.centre.y});f.centre={r.x,r.y};f.areaVector=vector(f.areaVector);f.correction=vector(f.correction);}
    return mesh;
}
void wallAndUnitScaling() {
    const auto mesh=rotated(rectangle(10,6,2,true),.37);
    const auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    auto state=constant(mesh,{1,0,0,1}),scaled=state;
    const double rhoScale=1.225,pressureScale=101325,velocityScale=std::sqrt(pressureScale/rhoScale);
    const std::array<double,4> units{rhoScale,rhoScale*velocityScale,rhoScale*velocityScale,pressureScale};
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto c=mesh.cells[i].centre;
        state.cells[i]=eulerConservative2D({1+.1*std::sin(c.x),.4+.05*std::sin(c.y),-.2,1+.1*std::cos(c.y)});
        for(std::size_t k=0;k<4;++k)scaled.cells[i][k]=state.cells[i][k]*units[k];
    }
    double maxError=0;
    for(unsigned step=0;step<40;++step) {
        auto controls=activeControls;controls.maximumStep=.002;
        const auto next=advanceEuler2D(mesh,bc,{},state,controls);
        controls.maximumStep/=velocityScale;
        const auto si=advanceEuler2D(mesh,bc,{},scaled,controls);
        for(const auto* result:{&next,&si}) {
            require(result->boundaryFlux[0]==0&&result->boundaryFlux[3]==0,"closed wall mass/energy leak");
            for(const auto& b:bc) {
                const auto& f=result->faceFlux[b.face];const auto n=mesh.faces[b.face].areaVector;
                require(f[0]==0&&f[3]==0,"wall flux symmetry not exact");
                require(std::abs(f[1]*n.y-f[2]*n.x)<=8*std::numeric_limits<double>::epsilon()*std::hypot(n.x,n.y)*std::hypot(f[1],f[2]),"wall has tangential traction");
            }
        }
        for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)
            maxError=std::max(maxError,std::abs(si.state.cells[i][k]/units[k]-next.state.cells[i][k])/(1+std::abs(next.state.cells[i][k])));
        state=next.state;scaled=si.state;
    }
    // Roundoff-only equivalence after 40 accepted steps, normalized back to the
    // original units; this is not a physical solution-accuracy requirement.
    require(maxError<256*std::numeric_limits<double>::epsilon()*state.steps,"Euler depends on chosen physical units");
    std::cout<<"rotated wall SI scaling: max normalized field difference "<<maxError<<'\n';
}
void shockRotationAndCheckpoint() {
    const auto mesh=rectangle(100,4,1);const double angle=.37;const auto rotatedMesh=rotated(mesh,angle);
    auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(auto& b:bc)if(std::abs(mesh.faces[b.face].areaVector.x)>0){b.kind=EulerBoundaryKind2D::Transmissive;b.name="ends";}
    EulerState2D state=constant(mesh,{}),rotatedState;
    for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D(mesh.cells[i].centre.x<.5?EulerPrimitive2D{1,0,0,1}:EulerPrimitive2D{.125,0,0,.1});
    rotatedState=state;const auto initial=state;double localRotationError=0;
    while(state.time<.2) {
        auto controls=activeControls;controls.maximumStep=.2-state.time;
        const auto next=advanceEuler2D(mesh,bc,{},state,controls);
        if(activeControls.order==2&&activeControls.fluxScheme==EulerFluxScheme2D::Hllc) {
            auto exactRotated=state;
            for(auto& q:exactRotated.cells){const auto x=q[1],y=q[2];q[1]=x*std::cos(angle)-y*std::sin(angle);q[2]=x*std::sin(angle)+y*std::cos(angle);}
            auto cap=controls;cap.maximumStep=.9*next.step;
            const auto referenceStep=advanceEuler2D(mesh,bc,{},state,cap),rotatedStep=advanceEuler2D(rotatedMesh,bc,{},exactRotated,cap);
            double local=0;
            for(std::size_t i=0;i<state.cells.size();++i) {
                const auto& q=referenceStep.state.cells[i];const auto& r=rotatedStep.state.cells[i];
                const EulerConservative2D expected{q[0],q[1]*std::cos(angle)-q[2]*std::sin(angle),q[1]*std::sin(angle)+q[2]*std::cos(angle),q[3]};
                for(std::size_t k=0;k<4;++k)local=std::max(local,std::abs(expected[k]-r[k]));
            }
            if(local>localRotationError){localRotationError=local;if(local>1e-10)std::cerr<<"local rotation: step "<<state.steps<<" error "<<local<<" dt "<<referenceStep.step<<' '<<rotatedStep.step<<'\n';}
        }
        controls.maximumStep=next.step;
        rotatedState=advanceEuler2D(rotatedMesh,bc,{},rotatedState,controls).state;state=next.state;
        require(next.minimumDensity>0&&next.minimumPressure>0,"Sod positivity failed");
        if(state.steps==50) {
            std::ostringstream stream;writeEulerCheckpoint2D(stream,mesh,bc,{},state,"Sod");
            std::istringstream input(stream.str());const auto restored=readEulerCheckpoint2D(input,mesh,bc,{},"Sod");
            same(state,restored,1e-15);require(state.time==restored.time&&state.steps==restored.steps,"restart counter differs");
            same(advanceEuler2D(mesh,bc,{},state,controls).state,advanceEuler2D(mesh,bc,{},restored,controls).state,1e-15);
            rejects([&]{std::istringstream in(stream.str());(void)readEulerCheckpoint2D(in,mesh,bc,{1.3,287.05},"Sod");});
            rejects([&]{std::istringstream in(stream.str()+"extra");(void)readEulerCheckpoint2D(in,mesh,bc,{},"Sod");});
            rejects([&]{std::istringstream in(stream.str());(void)readEulerCheckpoint2D(in,rotatedMesh,bc,{},"Sod");});
        }
    }
    double rotationError=0;
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto a=eulerPrimitive2D(state.cells[i]),b=eulerPrimitive2D(rotatedState.cells[i]);
        rotationError=std::max(rotationError,std::max({std::abs(a.density-b.density),std::abs(a.pressure-b.pressure),
            std::abs(a.u*std::cos(angle)-a.v*std::sin(angle)-b.u),std::abs(a.u*std::sin(angle)+a.v*std::cos(angle)-b.v)}));
    }
    std::cout<<"rotation diagnostic "<<std::setprecision(17)<<rotationError<<", max single-step "<<localRotationError<<", times "<<state.time<<' '<<rotatedState.time<<std::endl;
    // Keep the original first-order tolerance. The added nonlinear, two-stage
    // path uses a roundoff accumulation budget proportional to accepted steps,
    // in these nondimensional Sod variables; this is not a physical error gate.
    const double rotationBudget=activeControls.order==1?1e-12:256*std::numeric_limits<double>::epsilon()*static_cast<double>(state.steps);
    require(rotationError<rotationBudget,"Sod rotation covariance failed");
    require(localRotationError<1e-12,"single-step rotation covariance failed");
    std::cout<<"Sod 400 cells: "<<state.steps<<" steps, rotated state max difference "<<rotationError<<'\n';
    auto controls=EulerStepControls2D{};controls.minimumStep=.1;controls.maximumStep=1;
    rejects([&]{(void)advanceEuler2D(mesh,bc,{},initial,controls);});
    require(initial.time==0&&initial.steps==0,"failed step changed input");
}
void strongWaves() {
    const auto mesh=rectangle(80,2,1);
    const auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(bool expansion:{false,true}) {
        auto state=constant(mesh,{});
        for(std::size_t i=0;i<state.cells.size();++i) {
            const bool left=mesh.cells[i].centre.x<.5;
            state.cells[i]=eulerConservative2D(expansion?EulerPrimitive2D{1,left?-2.:2.,0,.4}:
                EulerPrimitive2D{1,0,0,left?1000.:.01});
        }
        double energy=0,mass=0;
        for(std::size_t i=0;i<state.cells.size();++i){energy+=mesh.cells[i].area*state.cells[i][3];mass+=mesh.cells[i].area*state.cells[i][0];}
        for(int j=0;j<200;++j) {
            const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);state=next.state;
            require(next.minimumDensity>0&&next.minimumPressure>0,"strong wave positivity failed");
            require(std::abs(next.afterIntegral[3]-energy)<1e-12*energy&&std::abs(next.afterIntegral[0]-mass)<1e-13,"closed blast/expansion lost energy or mass");
        }
    }
}
void faceFluxIdentities() {
    const auto state=eulerConservative2D({1.2,.4,-.3,2});
    for(auto scheme:{EulerFluxScheme2D::Rusanov,EulerFluxScheme2D::Hllc}) {
        const auto f=eulerFaceFlux2D(state,state,{.6,.8},{},scheme);const double vn=.4*.6-.3*.8;
        const EulerConservative2D exact{state[0]*vn,state[1]*vn+1.2,state[2]*vn+1.6,(state[3]+2)*vn};
        for(std::size_t k=0;k<4;++k)require(std::abs(f.integratedFlux[k]-exact[k])<3e-15,"consistent physical flux failed");
        for(int i=0;i<40;++i) {
            const auto l=eulerConservative2D({.1+.07*i,-3+.17*i,1-.09*i,.2+.13*i});
            const auto r=eulerConservative2D({3-.03*i,2-.19*i,-.7+.04*i,4-.09*i});
            const auto a=eulerFaceFlux2D(l,r,{.3,.4},{},scheme),b=eulerFaceFlux2D(r,l,{-.3,-.4},{},scheme);
            for(std::size_t k=0;k<4;++k)require(std::abs(a.integratedFlux[k]+b.integratedFlux[k])<1e-12*(1+std::abs(a.integratedFlux[k])),"flux owner reversal failed");
            require(a.hllcFallback==b.hllcFallback&&std::abs(a.waveSpeed-b.waveSpeed)<=8*std::numeric_limits<double>::epsilon()*a.waveSpeed,"wave/fallback owner reversal failed");
        }
    }
    for(double velocity:{-5.,0.,5.}) {
        const auto l=eulerConservative2D({1,velocity,2,1}),r=eulerConservative2D({3,velocity,-1,1});
        const auto f=eulerFaceFlux2D(l,r,{1,0},{},EulerFluxScheme2D::Hllc);const auto& q=velocity>=0?l:r;
        const EulerConservative2D exact{q[0]*velocity,q[1]*velocity+1,q[2]*velocity,(q[3]+1)*velocity};
        for(std::size_t k=0;k<4;++k)require(std::abs(f.integratedFlux[k]-exact[k])<2e-13,"HLLC contact/shear resolution failed");
        require(!f.hllcFallback,"admissible contact unexpectedly fell back");
    }
    const auto expansion=eulerFaceFlux2D(eulerConservative2D({1,-2,0,.4}),eulerConservative2D({1,2,0,.4}),{1,0},{},EulerFluxScheme2D::Hllc);
    require(expansion.hllcFallback,"non-positive HLLC star pressure did not report fallback");
    rejects([&]{(void)eulerFaceFlux2D(state,state,{0,0},{},EulerFluxScheme2D::Hllc);});
    rejects([&]{(void)eulerFaceFlux2D(state,state,{1,0},{},static_cast<EulerFluxScheme2D>(99));});
}
void stationaryContact() {
    const auto mesh=rectangle(40,6,2);const auto bc=periodic(mesh);auto state=constant(mesh,{});
    for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D({mesh.cells[i].centre.x<1?1.:3.,0,0,1});
    const auto initial=state;auto controls=activeControls;controls.fluxScheme=EulerFluxScheme2D::Hllc;
    for(int j=0;j<50;++j){const auto step=advanceEuler2D(mesh,bc,{},state,controls);require(step.hllcFallbackEvaluations==0,"contact fallback");state=step.state;}
    same(initial,state,2e-14);
}
void perturbedNormalShock() {
    const double gamma=1.4,amplitude=1e-5;
    for(double mach:{3.,10.}) {
        const auto mesh=rectangle(48,12,2);auto bc=periodic(mesh);
        const EulerPrimitive2D left{1,mach*std::sqrt(gamma),0,1};
        const double ratio=(gamma+1)*mach*mach/((gamma-1)*mach*mach+2);
        const EulerPrimitive2D right{ratio,left.u/ratio,0,1+2*gamma/(gamma+1)*(mach*mach-1)};
        for(auto& b:bc)if(b.name=="periodic-x") {
            const bool inlet=mesh.faces[b.face].areaVector.x<0;
            b.kind=EulerBoundaryKind2D::Farfield;b.partner.reset();b.reference=inlet?left:right;b.name=inlet?"inlet":"outlet";
        }
        auto state=constant(mesh,left);
        for(std::size_t i=0;i<state.cells.size();++i) {
            const auto c=mesh.cells[i].centre;auto q=c.x<1?left:right;
            if(std::abs(c.x-1)<.1)q.density*=1+amplitude*((i/48)%2==0?1:-1);
            state.cells[i]=eulerConservative2D(q);
        }
        EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=2;
        double maximumVariation=0;
        const double end=.3/left.u; // fixed convected distance at either Mach number
        while(state.time<end) {
            controls.maximumStep=end-state.time;const auto next=advanceEuler2D(mesh,bc,{},state,controls);state=next.state;
            require(next.minimumDensity>0&&next.minimumPressure>0&&next.minimumContactRestoration<.9,"shock guard/positivity not active");
            for(int i=0;i<48;++i) {
                double low=INFINITY,high=-INFINITY;
                for(int j=0;j<12;++j){const double rho=state.cells[static_cast<std::size_t>(j*48+i)][0];low=std::min(low,rho);high=std::max(high,rho);}
                maximumVariation=std::max(maximumVariation,(high-low)/ratio);
            }
        }
        std::cout<<"perturbed normal shock: Mach="<<mach<<", maximum transverse density span/postshock density="<<maximumVariation<<'\n';
        // Input density sawtooth spans 2e-5 of the postshock density. Bound its
        // growth by one decade over the stated window; this is not a universal
        // shock-stability claim and does not tune the product shock sensor.
        require(maximumVariation<20*amplitude,"shock transverse perturbation grew by more than one decade");
    }
}
void smoothAccuracy() {
    const double pi=std::acos(-1.),end=.3;std::array<double,2> last{};
    for(unsigned order:{1u,2u}) {
        double prev=0;
        for(int nx:{20,40,80}) {
            const int ny=nx/2;const auto mesh=rectangle(nx,ny,2);const auto bc=periodic(mesh);
            auto state=constant(mesh,{});const double sinc=std::sin(pi/nx)/(pi/nx)*std::sin(pi/ny)/(pi/ny);
            const auto density=[&](Point2D c,double t){return 1+.2*sinc*std::sin(pi*(c.x-.3*t))*std::cos(2*pi*(c.y+.2*t));};
            for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D({density(mesh.cells[i].centre,0),.3,-.2,1});
            EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=order;
            while(state.time<end){controls.maximumStep=end-state.time;state=advanceEuler2D(mesh,bc,{},state,controls).state;}
            double error=0;
            for(std::size_t i=0;i<state.cells.size();++i)error+=mesh.cells[i].area*std::abs(state.cells[i][0]-density(mesh.cells[i].centre,end));
            error/=2;
            std::cout<<"smooth entropy wave: order="<<order<<", nx="<<nx<<", area-weighted density L1="<<error;
            if(prev>0)std::cout<<", observed order="<<std::log2(prev/error);
            std::cout<<'\n';
            if(nx==80) {
                // A smooth periodic transport case tests the implemented order.
                // A BJ limiter reduces accuracy near extrema; 1.4 is a bounded
                // development regression floor, not a general CFD accuracy gate.
                require(prev/error>(order==2?std::pow(2.,1.4):1.5),"smooth grid refinement lost expected order");
                last[order-1]=error;
            }
            prev=error;
        }
    }
    require(last[1]<.3*last[0],"second-order smooth transport did not improve resolution");
}
void heatOperatorChecks() {
    const auto mesh=rotated(rectangle(12,8,2,true),.37);
    constexpr double conductivity=2.7;
    const Vector2D gradient{20,-7};
    const auto temperature=[&](Point2D p){return 300+gradient.x*p.x+gradient.y*p.y;};
    for(bool mixed:{false,true}) {
        std::vector<HeatBoundary2D> bc;
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];if(f.neighbour)continue;
            if(mixed&&id%2)bc.push_back({id,HeatBoundaryKind2D::OutwardFlux,-conductivity*dot(gradient,f.areaVector)/std::hypot(f.areaVector.x,f.areaVector.y),{}});
            else bc.push_back({id,HeatBoundaryKind2D::Temperature,temperature(f.centre),{}});
        }
        HeatConductionOperator2D op(mesh,bc,conductivity);
        std::vector<double> t,capacity(mesh.cells.size(),1000);
        for(const auto& c:mesh.cells)t.push_back(temperature(c.centre));
        const auto heat=op.evaluate(t,capacity);double maximum=0;
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];const double exact=-conductivity*dot(gradient,f.areaVector);
            maximum=std::max(maximum,std::abs(heat.faceHeatFlux[id]-exact)/(conductivity*std::hypot(gradient.x,gradient.y)*std::hypot(f.areaVector.x,f.areaVector.y)));
        }
        require(maximum<1024*std::numeric_limits<double>::epsilon(),"nonorthogonal Fourier operator lost linear exactness");
        // Independently differentiate evaluated residuals, column by column.
        // The operator is affine: a centered 1 K perturbation measures its exact
        // Jacobian up to rounding. Its absolute row sum must include the full
        // nonorthogonal correction, not only a two-point diagonal estimate.
        std::vector<double> measuredNorm(t.size());
        for(std::size_t j=0;j<t.size();++j) {
            auto plus=t,minus=t;plus[j]+=1;minus[j]-=1;
            const auto a=op.evaluate(plus,capacity),b=op.evaluate(minus,capacity);
            for(std::size_t i=0;i<t.size();++i)measuredNorm[i]+=.25*std::abs(a.cellResidual[i]-b.cellResidual[i]);
        }
        double rateError=0;
        for(std::size_t i=0;i<t.size();++i) {
            const double expected=measuredNorm[i]/(mesh.cells[i].area*capacity[i]);
            rateError=std::max(rateError,std::abs(heat.rate[i]-expected)/expected);
        }
        require(rateError<2048*std::numeric_limits<double>::epsilon(),"full heat Jacobian row bound mismatch");
        std::cout<<"heat row norm finite perturbation: "<<rateError<<'\n';
        std::cout<<"heat linear rotated/skew mesh: mixed="<<mixed<<", normalized flux error="<<maximum<<'\n';
    }
    const double pi=std::acos(-1.),end=.03,alpha=.2;
    double previous=0;
    for(int n:{16,32,64}) {
        const auto grid=rectangle(n,4,1);const auto periodicBc=periodic(grid);std::vector<HeatBoundary2D> bc;
        for(const auto& b:periodicBc)bc.push_back({b.face,HeatBoundaryKind2D::Periodic,0,b.partner});
        HeatConductionOperator2D op(grid,bc,.5);
        std::vector<double> t,capacity(grid.cells.size(),2.5);
        const double average=std::sin(pi/n)/(pi/n);
        for(const auto& c:grid.cells)t.push_back(2+.2*average*std::cos(2*pi*c.centre.x));
        double time=0;
        while(time<end) {
            auto first=op.evaluate(t,capacity);double rate=0;for(double r:first.rate)rate=std::max(rate,r);
            const double dt=std::min(end-time,.4/rate);auto stage=t;
            for(std::size_t i=0;i<t.size();++i)stage[i]-=dt*first.cellResidual[i]/(capacity[i]*grid.cells[i].area);
            const auto second=op.evaluate(stage,capacity);
            for(std::size_t i=0;i<t.size();++i)t[i]-=.5*dt*(first.cellResidual[i]+second.cellResidual[i])/(capacity[i]*grid.cells[i].area);
            time+=dt;
        }
        double error=0,total=0;
        for(std::size_t i=0;i<t.size();++i) {
            const double exact=2+.2*average*std::cos(2*pi*grid.cells[i].centre.x)*std::exp(-4*pi*pi*alpha*end);
            error+=grid.cells[i].area*std::abs(t[i]-exact);total+=grid.cells[i].area*t[i];
        }
        require(std::abs(total-2)<256*std::numeric_limits<double>::epsilon(),"periodic diffusion loses energy");
        std::cout<<"heat Fourier decay: nx="<<n<<", temperature L1="<<error;
        if(previous)std::cout<<", observed order="<<std::log2(previous/error);std::cout<<'\n';
        if(n==64)require(previous/error>std::pow(2.,1.8),"Fourier diffusion lost second-order convergence");
        previous=error;
    }
}
void coupledHeatChecks() {
    const auto mesh=rectangle(8,4,1,true);auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    const IdealGas2D gas{1.4,1};const auto initial=constant(mesh,{1,0,0,1});
    EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=2;controls.maximumStep=.001;
    const auto legacy=advanceEuler2D(mesh,bc,gas,initial,controls);
    require(legacy.state.cells==EulerStepper2D(mesh,bc,gas).advance(initial,controls).state.cells,"cached zero-conductivity path differs");
    const auto uniform=EulerStepper2D(mesh,bc,gas,{10}).advance(initial,controls);
    same(uniform.state,initial,128*std::numeric_limits<double>::epsilon());
    for(auto q:uniform.faceHeatFlux)require(q==0,"uniform temperature has false heat flux");
    require(uniform.step<legacy.step&&uniform.combinedCourant<=.4*(1+1e-12),"diffusion not included in time step");
    for(auto& b:bc){b.thermalKind=HeatBoundaryKind2D::OutwardFlux;b.thermalValue=-2;}
    EulerStepper2D solver(mesh,bc,gas,{.5});auto state=initial;
    for(int n=0;n<20;++n) {
        const auto next=solver.advance(state,controls);state=next.state;
        const double balance=next.afterIntegral[3]-next.beforeIntegral[3]+next.step*next.boundaryHeat;
        require(std::abs(balance)<256*std::numeric_limits<double>::epsilon()*next.beforeIntegral[3],"heated closed gas violates total energy");
        for(const auto& b:bc){require(next.faceFlux[b.face][0]==0,"heated wall leaks mass");require(next.faceFlux[b.face][3]==next.faceHeatFlux[b.face],"wall energy differs from prescribed heat");}
    }
    require(state.cells[0][3]>initial.cells[0][3],"negative outward wall flux did not heat fluid");
    std::stringstream saved;writeEulerCheckpoint2D(saved,mesh,bc,gas,state,"sealed",{.5});const auto data=saved.str();
    require(data.starts_with("CM2D_EULER_CHECKPOINT 2"),"conductive checkpoint lacks transport version");
    same(state,readEulerCheckpoint2D(saved,mesh,bc,gas,"sealed",{.5}));
    rejects([&]{std::istringstream in(data);(void)readEulerCheckpoint2D(in,mesh,bc,gas,"sealed",{.6});});
    auto changed=bc;changed[0].thermalValue=-1;
    rejects([&]{std::istringstream in(data);(void)readEulerCheckpoint2D(in,mesh,changed,gas,"sealed",{.5});});
    rejects([&]{(void)EulerStepper2D(mesh,bc,gas,{0});});
    rejects([&]{(void)EulerStepper2D(mesh,bc,gas,{-1});});
    for(auto& b:bc)b.thermalValue=1e8;
    controls.maximumRetries=0;controls.maximumStep=.01;
    rejects([&]{(void)EulerStepper2D(mesh,bc,gas,{.5}).advance(initial,controls);});
    controls.maximumRetries=30;controls.minimumStep=1e-16;
    const auto cooled=EulerStepper2D(mesh,bc,gas,{.5}).advance(initial,controls);
    require(cooled.rejectedCandidates>0&&cooled.minimumPressure>0&&initial.time==0,"cooling positivity retry/initial state preservation failed");
    controls.timeStepControl=EulerTimeStepControl2D::StageGuarded;
    const auto guardedCooling=EulerStepper2D(mesh,bc,gas,{.5}).advance(initial,controls);
    require(guardedCooling.rejectedCandidates>0&&guardedCooling.minimumPressure>0&&initial.time==0,"stage guard bypassed cooling positivity retries");
    controls.maximumRetries=0;
    rejects([&]{(void)EulerStepper2D(mesh,bc,gas,{.5}).advance(initial,controls);});
    std::cout<<"coupled heat: conservation, combined CFL, zero-k identity, wall signs, cooling retries and physical restart binding passed\n";
}
void stageGuardedChecks() {
    const auto mesh=rectangle(12,4,2,true);const auto bc=periodic(mesh);
    const IdealGas2D gas{1.4,1};const EulerTransport2D transport{.003,.002};
    auto initial=constant(mesh,{1,.4,0,1});
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const double wave=.05*std::sin(std::acos(-1.)*mesh.cells[i].centre.x);
        initial.cells[i]=eulerConservative2D({1+wave,.4+wave,.02*std::cos(std::acos(-1.)*mesh.cells[i].centre.x),std::pow(1+wave,gas.gamma)},gas);
    }
    const EulerStepper2D solver(mesh,bc,gas,transport);
    EulerStepControls2D controls;controls.order=2;controls.fluxScheme=EulerFluxScheme2D::Hllc;
    controls.timeStepControl=static_cast<EulerTimeStepControl2D>(91);
    rejects([&]{(void)solver.advance(initial,controls);});
    controls.order=1;controls.timeStepControl=EulerTimeStepControl2D::Legacy;
    const auto forward=solver.advance(initial,controls);
    controls.timeStepControl=EulerTimeStepControl2D::StageGuarded;
    const auto guardedForward=solver.advance(initial,controls);
    require(forward.state.cells==guardedForward.state.cells&&forward.step==guardedForward.step,"stage guard changed first-order stepping");
    controls.order=2;controls.endTime=.8;
    const auto evolve=[&](double cfl) {
        auto state=initial;controls.acousticCourant=cfl;
        while(state.time<*controls.endTime) {
            const auto step=solver.advance(state,controls);
            require(step.step>=controls.minimumStep&&step.step<=controls.maximumStep,"guarded time limits violated");
            require(step.combinedCourant<=cfl*(1+8*std::numeric_limits<double>::epsilon()),"guarded stage exceeded original combined CFL limit");
            require(step.spatialEvaluations>=2&&step.cflRejectedCandidates<=step.rejectedCandidates,"invalid stage work counters");
            state=step.state;
            if(state.steps==8) {
                std::stringstream checkpoint;writeEulerCheckpoint2D(checkpoint,mesh,bc,gas,state,"guarded",transport);
                const auto restored=readEulerCheckpoint2D(checkpoint,mesh,bc,gas,"guarded",transport);
                const auto a=solver.advance(state,controls),b=EulerStepper2D(mesh,bc,gas,transport).advance(restored,controls);
                require(a.state.cells==b.state.cells&&a.state.time==b.state.time&&a.rejectedCandidates==b.rejectedCandidates,"guarded restart depends on unrecorded controller history");
            }
        }
        require(state.time==*controls.endTime,"guarded integration missed the exact horizon");return state;
    };
    const auto coarse=evolve(.4),medium=evolve(.2),fine=evolve(.1),reference=evolve(.025),checkedReference=evolve(.0125);
    const auto error=[&](const EulerState2D& a,const EulerState2D& b){double sum=0;for(std::size_t i=0;i<a.cells.size();++i)for(std::size_t k=0;k<4;++k){const double d=a.cells[i][k]-b.cells[i][k];sum+=mesh.cells[i].area*d*d;}return std::sqrt(sum);};
    const double e0=error(coarse,checkedReference),e1=error(medium,checkedReference),e2=error(fine,checkedReference);
    const double order=std::log2(e1/e2),referenceDifference=error(reference,checkedReference);
    // Same semi-discrete equations/mesh: retain the existing temporal target
    // >1.8, with reference uncertainty below 10% of the finest measured error.
    std::cout<<std::setprecision(17)<<"stage-guarded temporal errors="<<e0<<','<<e1<<','<<e2<<" order="<<order<<" reference difference="<<referenceDifference<<'\n';
    require(e0>e1&&order>1.8&&referenceDifference<.1*e2,"stage-guarded temporal refinement failed");
    controls.acousticCourant=.4;controls.minimumStep=.0001;controls.maximumStep=.001;controls.endTime=.00105;
    const auto penultimate=solver.advance(initial,controls),last=solver.advance(penultimate.state,controls);
    require(penultimate.step>=controls.minimumStep&&last.step>=controls.minimumStep&&last.state.time==*controls.endTime&&last.state.steps==2,"stage guard mishandled a sub-minimum final tail");
    auto hot=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(auto& b:hot){b.thermalKind=HeatBoundaryKind2D::OutwardFlux;b.thermalValue=-1000;}
    controls.minimumStep=1e-14;controls.maximumStep=1;controls.endTime.reset();
    const auto heated=EulerStepper2D(mesh,hot,gas,{.5}).advance(constant(mesh,{1,0,0,1}),controls);
    require(heated.cflRejectedCandidates>0&&heated.minimumPressure>0&&heated.combinedCourant<=controls.acousticCourant*(1+8*std::numeric_limits<double>::epsilon()),"guarded retry did not handle a rapidly increasing stage rate");
    controls.minimumStep=.1;controls.maximumStep=.1;
    const auto before=initial;rejects([&]{(void)solver.advance(initial,controls);});
    require(initial.cells==before.cells&&initial.time==before.time,"minimum-step failure changed the input state");
}
void totalInletAndImplicitChecks() {
    const IdealGas2D gas{1.4,1};const EulerPrimitive2D q{1,.3,.07,1};
    const auto inlet=eulerTotalInletState2D(q,q,{-1,0},gas);
    require(std::abs(inlet.u-q.u)<1e-14&&std::abs(inlet.pressure-q.pressure)<1e-14,"total inlet free stream mismatch");
    const auto altered=eulerTotalInletState2D({1,.32,.07,1},q,{-1,0},gas);
    const auto h=[&](const EulerPrimitive2D& a){return gas.gamma/(gas.gamma-1)*a.pressure/a.density+.5*(a.u*a.u+a.v*a.v);};
    require(std::abs(h(altered)-h(q))<1e-14&&altered.v==q.v,"total inlet prescribed incorrect incoming invariants");
    require(std::abs(altered.pressure/std::pow(altered.density,gas.gamma)-1)<1e-14,"total inlet entropy changed");
    rejects([&]{(void)eulerTotalInletState2D({1,-10,0,1},q,{-1,0},gas);});
    const auto mesh=rectangle(12,4,2,true);const auto bc=periodic(mesh);
    auto initial=constant(mesh,{1,0,0,1});
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const double wave=std::sin(2*std::acos(-1.)*mesh.cells[i].centre.x/2);
        initial.cells[i]=eulerConservative2D({1+.001*wave,.02*wave,.01*wave,1+.0014*wave},gas);
    }
    const EulerStepper2D solver(mesh,bc,gas,{.02,.02});EulerStepControls2D ctl;
    ctl.integrator=EulerTimeIntegrator2D::Sdirk2;ctl.fluxScheme=EulerFluxScheme2D::Hllc;ctl.order=2;ctl.endTime=.08;
    const auto evolve=[&](double dt) {
        ctl.maximumStep=dt;auto state=initial;
        while(state.time<*ctl.endTime) {
            const auto next=solver.advance(state,ctl);
            require(next.maximumCellBalanceError<1e-12&&next.minimumPressure>0&&next.minimumDensity>0,"implicit stage conservation/positivity failed");
            std::vector<double> spectral(mesh.cells.size());
            for(std::size_t f=0;f<mesh.faces.size();++f) {
                const auto& face=mesh.faces[f];
                const double value=next.faceWaveSpeed[f]*std::hypot(face.areaVector.x,face.areaVector.y);
                spectral[face.owner]+=value;if(face.neighbour)spectral[*face.neighbour]+=value;
            }
            double acoustic=0,combined=0;
            for(std::size_t i=0;i<spectral.size();++i) {
                acoustic=std::max(acoustic,next.step*spectral[i]/mesh.cells[i].area);
                combined=std::max(combined,next.step*(spectral[i]/mesh.cells[i].area+next.cellHeatRate[i]+next.cellViscousRate[i]));
            }
            require(acoustic==next.acousticCourant&&combined==next.combinedCourant,"implicit Courant diagnostic differs from exported face maxima");
            if(state.steps==1) {
                std::stringstream file;writeEulerCheckpoint2D(file,mesh,bc,gas,state,"implicit",{.02,.02});
                const auto restored=readEulerCheckpoint2D(file,mesh,bc,gas,"implicit",{.02,.02});
                require(solver.advance(restored,ctl).state.cells==next.state.cells,"implicit restart has hidden history");
            }
            state=next.state;
        }
        return state;
    };
    const auto coarse=evolve(.005),medium=evolve(.0025),fine=evolve(.00125),reference=evolve(.00015625),checkedReference=evolve(.0003125);
    const auto err=[&](const EulerState2D& a){double v=0;for(std::size_t i=0;i<a.cells.size();++i)for(std::size_t k=0;k<4;++k)v+=mesh.cells[i].area*std::pow(a.cells[i][k]-reference.cells[i][k],2);return std::sqrt(v);};
    const double e0=err(coarse),e1=err(medium),e2=err(fine),order=std::log2(e1/e2);
    std::cout<<"SDIRK2 errors="<<e0<<","<<e1<<","<<e2<<" order="<<order<<std::endl;
    require(e0>e1&&order>1.8&&err(checkedReference)<.1*e2,"SDIRK2 temporal refinement/reference failed");
    ctl.endTime=.020000000001;const auto tail=evolve(.01);
    require(tail.time==*ctl.endTime&&tail.steps==3,"implicit tiny physical tail failed");
    ctl.interrupted=[] {return true;};rejects([&]{(void)solver.advance(initial,ctl);});ctl.interrupted={};
    // No field repair when Newton cannot complete: the input is immutable.
    ctl.maximumStep=1;ctl.endTime.reset();ctl.maximumNewtonIterations=1;ctl.maximumRetries=0;
    const auto before=initial;rejects([&]{(void)solver.advance(initial,ctl);});
    require(initial.cells==before.cells&&initial.time==before.time,"failed Newton candidate replaced accepted input");
}
int main() {
    try {
        const EulerPrimitive2D p{1.225,320,-12,101325};const auto round=eulerPrimitive2D(eulerConservative2D(p));
        require(std::abs(round.pressure-p.pressure)<1e-10&&round.density==p.density,"ideal gas round trip failed");
        rejects([]{(void)eulerConservative2D({-1,0,0,1});});
        rejects([]{(void)eulerPrimitive2D({1,10,0,1});});
        rejects([]{validateIdealGas2D({1,287});});
        totalInletAndImplicitChecks();faceFluxIdentities();
        for(auto scheme:{EulerFluxScheme2D::Rusanov,EulerFluxScheme2D::Hllc})for(unsigned order:{1u,2u}) {
            activeControls.fluxScheme=scheme;activeControls.order=order;
            std::cout<<"scheme="<<(scheme==EulerFluxScheme2D::Hllc?"hllc":"rusanov")<<", order="<<order<<std::endl;
            pressureOutletChecks();restAndFreeStream();periodicConservation();wallAndUnitScaling();shockRotationAndCheckpoint();strongWaves();stationaryContact();
        }
        perturbedNormalShock();smoothAccuracy();heatOperatorChecks();coupledHeatChecks();stageGuardedChecks();
        std::cout<<"Euler conservation, acoustic CFL, boundaries, rotation and checkpoint passed\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
